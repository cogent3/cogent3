#!/usr/bin/env python
"""This file controls the central function of EVOLVE, the calculation of the log-
likelihood of an alignment given a phylogenetic tree and substitution model.
The likelihood calculation is done according to Felsenstein's 1981 pruning
algorithm. This file contains a Python implementation of that
algorithm and an interface to a more computationally efficient Pyrex
implementation. The two versions are maintained for the purpose of cross-
validating accuracy. The calculations can be performed for tree's that have polytomies
in addition to binary trees.
"""
import numpy

from cogent3.evolve.likelihood_tree import LikelihoodTreeEdge
from cogent3.evolve.simulate import argpick
from cogent3.maths.markov import SiteClassTransitionMatrix
from cogent3.recalculation.definition import (
    CalcDefn,
    CalculationDefn,
    CallDefn,
    NonParamDefn,
    ProbabilityParamDefn,
    SumDefn,
)


Float = numpy.core.numerictypes.sctype2char(float)


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


class _PartialLikelihoodDefn(CalculationDefn):
    def setup(self, edge_name):
        self.edge_name = edge_name


class LeafPartialLikelihoodDefn(_PartialLikelihoodDefn):
    name = "sequence"

    def calc(self, lh_tree):
        lh_leaf = lh_tree.get_edge(self.edge_name)
        return lh_leaf.input_likelihoods


class PartialLikelihoodProductDefn(_PartialLikelihoodDefn):
    name = "plh"
    recycling = True

    def calc(self, recycled_result, lh_edge, *child_likelihoods):
        if recycled_result is None:
            recycled_result = lh_edge.make_partial_likelihoods_array()
        return lh_edge.sum_input_likelihoodsR(recycled_result, *child_likelihoods)


class PartialLikelihoodProductDefnFixedMotif(PartialLikelihoodProductDefn):
    def calc(self, recycled_result, fixed_motif, lh_edge, *child_likelihoods):
        if recycled_result is None:
            recycled_result = lh_edge.make_partial_likelihoods_array()
        result = lh_edge.sum_input_likelihoodsR(recycled_result, *child_likelihoods)
        if fixed_motif not in [None, -1]:
            for motif in range(result.shape[-1]):
                if motif != fixed_motif:
                    result[:, motif] = 0.0
        return result


class LhtEdgeLookupDefn(CalculationDefn):
    name = "col_index"

    def setup(self, edge_name):
        self.edge_name = edge_name
        # so that it can be found by reconstruct_ancestral_seqs etc:
        if edge_name == "root":
            self.name = "root"

    def calc(self, lht):
        return lht.get_edge(self.edge_name)


def make_partial_likelihood_defns(edge, lht, psubs, fixed_motifs):
    kw = {"edge_name": edge.name}

    if edge.istip():
        plh = LeafPartialLikelihoodDefn(lht, **kw)
    else:
        lht_edge = LhtEdgeLookupDefn(lht, **kw)
        children = []
        for child in edge.children:
            child_plh = make_partial_likelihood_defns(child, lht, psubs, fixed_motifs)
            psub = psubs.select_from_dimension("edge", child.name)
            child_plh = CalcDefn(numpy.inner)(child_plh, psub)
            children.append(child_plh)

        if fixed_motifs:
            fixed_motif = fixed_motifs.select_from_dimension("edge", edge.name)
            plh = PartialLikelihoodProductDefnFixedMotif(
                fixed_motif, lht_edge, *children, **kw
            )
        else:
            plh = PartialLikelihoodProductDefn(lht, *children, **kw)

    return plh


def recursive_lht_build(edge, leaves):
    if edge.istip():
        lhe = leaves[edge.name]
    else:
        lht_children = []
        for child in edge.children:
            lht = recursive_lht_build(child, leaves)
            lht_children.append(lht)
        lhe = LikelihoodTreeEdge(lht_children, edge_name=edge.name)
    return lhe


class LikelihoodTreeDefn(CalculationDefn):
    name = "lht"

    def setup(self, tree):
        self.tree = tree

    def calc(self, leaves):
        return recursive_lht_build(self.tree, leaves)


def make_total_loglikelihood_defn(
    tree, leaves, psubs, mprobs, bprobs, bin_names, locus_names, sites_independent
):

    fixed_motifs = NonParamDefn("fixed_motif", ["edge"])

    lht = LikelihoodTreeDefn(leaves, tree=tree)
    plh = make_partial_likelihood_defns(tree, lht, psubs, fixed_motifs)

    # After the root partial likelihoods have been calculated it remains to
    # sum over the motifs, local sites, other sites (ie: cpus), bins and loci.
    # The motifs are always done first, but after that it gets complicated.
    # If a bin HMM is being used then the sites from the different CPUs must
    # be interleaved first, otherwise summing over the CPUs is done last to
    # minimise inter-CPU communicaton.

    root_mprobs = mprobs.select_from_dimension("edge", "root")
    lh = CalcDefn(numpy.inner, name="lh")(plh, root_mprobs)
    if len(bin_names) > 1:
        if sites_independent:
            site_pattern = CalcDefn(BinnedSiteDistribution, name="bdist")(bprobs)
        else:
            switch = ProbabilityParamDefn("bin_switch", dimensions=["locus"])
            site_pattern = CalcDefn(PatchSiteDistribution, name="bdist")(switch, bprobs)
        blh = CallDefn(site_pattern, lht, name="bindex")
        tll = CallDefn(blh, *lh.across_dimension("bin", bin_names), **dict(name="tll"))
    else:
        lh = lh.select_from_dimension("bin", bin_names[0])
        tll = CalcDefn(log_sum_across_sites, name="logsum")(lht, lh)

    if len(locus_names) > 1:
        # currently has no .make_likelihood_function() method.
        tll = SumDefn(*tll.across_dimension("locus", locus_names))
    else:
        tll = tll.select_from_dimension("locus", locus_names[0])

    return tll


def log_sum_across_sites(root, root_lh):
    return root.get_log_sum_across_sites(root_lh)


class BinnedSiteDistribution(object):
    def __init__(self, bprobs):
        self.bprobs = bprobs

    def get_weighted_sum_lh(self, lhs):
        result = numpy.zeros(lhs[0].shape, lhs[0].dtype.char)
        temp = numpy.empty(result.shape, result.dtype.char)
        for (bprob, lh) in zip(self.bprobs, lhs):
            temp[:] = lh
            temp *= bprob
            result += temp
        return result

    def __call__(self, root):
        return BinnedLikelihood(self, root)

    def emit(self, length, random_series):
        result = numpy.zeros([length], int)
        for i in range(length):
            result[i] = argpick(self.bprobs, random_series)
        return result


class PatchSiteDistribution(object):
    def __init__(self, switch, bprobs):
        half = len(bprobs) // 2
        self.alloc = [0] * half + [1] * (len(bprobs) - half)

        pprobs = numpy.zeros([max(self.alloc) + 1], Float)
        for (b, p) in zip(self.alloc, bprobs):
            pprobs[b] += p

        self.bprobs = [p / pprobs[self.alloc[i]] for (i, p) in enumerate(bprobs)]
        self.transition_matrix = SiteClassTransitionMatrix(switch, pprobs)

    def get_weighted_sum_lhs(self, lhs):
        result = numpy.zeros((2,) + lhs[0].shape, lhs[0].dtype.char)
        temp = numpy.empty(lhs[0].shape, result.dtype.char)
        for (patch, weight, lh) in zip(self.alloc, self.bprobs, lhs):
            temp[:] = lh
            temp *= weight
            result[patch] += temp
        return result

    def __call__(self, root):
        return SiteHmm(self, root)

    def emit(self, length, random_series):
        bprobs = [
            [p for (patch, p) in zip(self.alloc, self.bprobs) if patch == a]
            for a in [0, 1]
        ]
        source = self.transition_matrix.emit(random_series)
        result = numpy.zeros([length], int)
        for i in range(length):
            patch = next(source) - 1
            result[i] = argpick(bprobs[patch], random_series)
        return result


class BinnedLikelihood(object):
    def __init__(self, distrib, root):
        self.distrib = distrib
        self.root = root

    def __call__(self, *lhs):
        result = self.distrib.get_weighted_sum_lh(lhs)
        return self.root.get_log_sum_across_sites(result)

    def get_posterior_probs(self, *lhs):
        # posterior bin probs, not motif probs
        assert len(lhs) == len(self.distrib.bprobs)
        result = numpy.array(
            [
                b * self.root.get_full_length_likelihoods(p)
                for (b, p) in zip(self.distrib.bprobs, lhs)
            ]
        )
        result /= result.sum(axis=0)
        return result


class SiteHmm(object):
    def __init__(self, distrib, root):
        self.root = root
        self.distrib = distrib

    def __call__(self, *lhs):
        plhs = self.distrib.get_weighted_sum_lhs(lhs)
        plhs = numpy.ascontiguousarray(numpy.transpose(plhs))
        matrix = self.distrib.transition_matrix
        return self.root.log_dot_reduce(matrix.StationaryProbs, matrix.Matrix, plhs)

    def get_posterior_probs(self, *lhs):
        plhs = [
            self.root.get_full_length_likelihoods(lh)
            for lh in self.distrib.get_weighted_sum_lhs(lhs)
        ]
        plhs = numpy.transpose(plhs)
        pprobs = self.distrib.transition_matrix.get_posterior_probs(plhs)
        pprobs = numpy.array(numpy.transpose(pprobs))

        lhs = numpy.array(lhs)
        blhs = lhs / numpy.sum(lhs, axis=0)
        blhs = numpy.array(
            [
                b * self.root.get_full_length_likelihoods(p)
                for (b, p) in zip(self.distrib.bprobs, blhs)
            ]
        )

        binsum = numpy.zeros(pprobs.shape, Float)
        for (patch, data) in zip(self.distrib.alloc, blhs):
            binsum[patch] += data

        for (patch, data) in zip(self.distrib.alloc, blhs):
            data *= pprobs[patch] / binsum[patch]

        return blhs
