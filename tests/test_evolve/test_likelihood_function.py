#!/usr/bin/env python
"""
Some tests for the likelihood function class.

tests to do:
    setting of parameters, by coord, by for-all, checking pars sets
    testing the likelihood for specified pars
    getting ancestral probs
    simulating sequence (not possible to verify values as random)
    
    checking that the object resets on tree change, model change, etc
"""
import json
import os
import warnings

from unittest import TestCase, main

import numpy

from numpy import dot, ones
from numpy.testing import assert_allclose

from cogent3 import (
    DNA,
    load_aligned_seqs,
    load_tree,
    make_aligned_seqs,
    make_tree,
)
from cogent3.evolve import ns_substitution_model, predicate, substitution_model
from cogent3.evolve.models import (
    CNFGTR,
    GN,
    GTR,
    HKY85,
    JTT92,
    MG94HKY,
    Y98,
    get_model,
    ssGN,
)
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3.maths.matrix_exponentiation import PadeExponentiator as expm
from cogent3.maths.stats.information_criteria import aic, bic


warnings.filterwarnings("ignore", "Motif probs overspecified")
warnings.filterwarnings("ignore", "Ignoring tree edge lengths")

TimeReversibleNucleotide = substitution_model.TimeReversibleNucleotide
MotifChange = predicate.MotifChange

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Rob Knight",
    "Matthew Wakefield",
    "Brett Easton",
    "Ananias Iliadis",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, "data")

ALIGNMENT = load_aligned_seqs(
    moltype=DNA, filename=os.path.join(data_path, "brca1.fasta")
)

OTU_NAMES = ["Human", "Mouse", "HowlerMon"]

_data = {
    "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
    "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
    "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
}
_aln = make_aligned_seqs(data=_data, moltype=DNA)


########################################################
# some funcs for assembling Q-matrices for 'manual' calc


def isTransition(motif1, motif2):
    position = getposition(motif1, motif2)
    a, b = motif1[position], motif2[position]
    transitions = {("A", "G"): 1, ("C", "T"): 1}
    pair = (min(a, b), max(a, b))

    return pair in transitions


def numdiffs_position(motif1, motif2):
    assert len(motif1) == len(
        motif2
    ), f"motif1[{motif1}] & motif2[{motif2}] have inconsistent length"

    ndiffs, position = 0, -1
    for i in range(len(motif1)):
        if motif1[i] != motif2[i]:
            position = i
            ndiffs += 1

    return ndiffs == 1, position


def isinstantaneous(motif1, motif2):
    if motif1 != motif2 and (
        motif1 == "-" * len(motif1) or motif2 == "-" * len(motif1)
    ):
        return True
    ndiffs, position = numdiffs_position(motif1, motif2)
    return ndiffs


def getposition(motif1, motif2):
    ndiffs, position = numdiffs_position(motif1, motif2)
    return position


##############################################################
# funcs for testing the monomer weighted substitution matrices
def _root_probs(x):
    return dict(
        [(n1 + n2, p1 * p2) for n1, p1 in list(x.items()) for n2, p2 in list(x.items())]
    )


def make_p(length, coord, val):
    """returns a probability matrix with value set at coordinate in
    instantaneous rate matrix"""
    Q = ones((4, 4), float) * 0.25  # assumes equi-frequent mprobs at root
    for i in range(4):
        Q[i, i] = 0.0
    Q[coord] *= val
    row_sum = Q.sum(axis=1)
    scale = 1 / (0.25 * row_sum).sum()
    for i in range(4):
        Q[i, i] -= row_sum[i]
    Q *= scale
    return expm(Q)(length)


def next_pi(pi, P):
    """returns next pi from apply a P matrix to it"""
    return dot(pi, P)


class LikelihoodCalcs(TestCase):
    """tests ability to calculate log-likelihoods for several
    substitution models."""

    def setUp(self):
        self.alignment = ALIGNMENT.take_seqs(OTU_NAMES)[0:42]
        self.tree = make_tree(tip_names=OTU_NAMES)

    def _makeLikelihoodFunction(self, submod, translate=False, **kw):
        alignment = self.alignment
        if translate:
            alignment = alignment.get_translation()
        calc = submod.make_likelihood_function(self.tree, **kw)
        calc.set_alignment(alignment)
        calc.set_param_rule("length", value=1.0, is_constant=True)
        if not translate:
            calc.set_param_rule("kappa", value=3.0, is_constant=True)
        return calc

    def test_no_seq_named_root(self):
        """root is a reserved name"""
        aln = self.alignment.take_seqs(self.alignment.names[:4])
        aln = aln.to_dict()
        one = aln.pop("Mouse")
        aln["root"] = one
        aln = make_aligned_seqs(data=aln)
        submod = get_model("TN93")
        tree = make_tree(f"{str(tuple(aln.names))}")
        lf = submod.make_likelihood_function(tree)
        try:
            lf.set_alignment(aln)
        except AssertionError:
            pass

        collection = aln.degap().named_seqs
        collection.pop("Human")
        tree = make_tree(f"{str(tuple(collection.keys()))}")
        lf = submod.make_likelihood_function(tree, aligned=False)
        try:
            lf.set_sequences(collection)
        except AssertionError:
            pass

    def test_binned_gamma(self):
        """just rate is gamma distributed"""
        submod = substitution_model.TimeReversibleCodon(
            predicates={"kappa": "transition", "omega": "replacement"},
            ordered_param="rate",
            distribution="gamma",
            mprob_model="tuple",
        )
        lf = self._makeLikelihoodFunction(submod, bins=3)
        try:
            values = list(lf.get_param_value_dict(["bin"])["omega_factor"].values())
        except KeyError:
            # there shouldn't be an omega factor
            pass
        values = list(lf.get_param_value_dict(["bin"])["rate"].values())
        obs = round(sum(values) / len(values), 6)
        self.assertEqual(obs, 1.0)
        self.assertEqual(len(values), 3)
        lf.get_param_value("rate_shape")

    def test_binned_gamma_ordered_param(self):
        """rate is gamma distributed omega follows"""
        submod = substitution_model.TimeReversibleCodon(
            predicates={"kappa": "transition", "omega": "replacement"},
            ordered_param="rate",
            partitioned_params="omega",
            distribution="gamma",
            mprob_model="tuple",
        )
        lf = self._makeLikelihoodFunction(submod, bins=3)
        values = list(lf.get_param_value_dict(["bin"])["omega_factor"].values())
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 3)
        lf.get_param_value("rate_shape")

    def test_binned_partition(self):
        submod = substitution_model.TimeReversibleCodon(
            predicates={"kappa": "transition", "omega": "replacement"},
            ordered_param="rate",
            partitioned_params="omega",
            distribution="free",
            mprob_model="tuple",
        )
        lf = self._makeLikelihoodFunction(submod, bins=3)
        values = list(lf.get_param_value_dict(["bin"])["omega_factor"].values())
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 3)

    def test_complex_binned_partition(self):
        submod = substitution_model.TimeReversibleCodon(
            predicates={"kappa": "transition", "omega": "replacement"},
            ordered_param="kappa",
            partitioned_params=["omega"],
            mprob_model="tuple",
        )
        lf = self._makeLikelihoodFunction(submod, bins=["slow", "fast"])
        lf.set_param_rule("kappa", value=1.0, is_constant=True)
        lf.set_param_rule("kappa", edge="Human", init=1.0, is_constant=False)
        values = list(lf.get_param_value_dict(["bin"])["kappa_factor"].values())
        self.assertEqual(round(sum(values) / len(values), 6), 1.0)
        self.assertEqual(len(values), 2)

    def test_codon(self):
        """test a three taxa codon model."""
        submod = substitution_model.TimeReversibleCodon(
            equal_motif_probs=True,
            motif_probs=None,
            predicates={"kappa": "transition", "omega": "replacement"},
            mprob_model="tuple",
        )

        likelihood_function = self._makeLikelihoodFunction(submod)
        likelihood_function.set_param_rule("omega", value=0.5, is_constant=True)
        evolve_lnL = likelihood_function.get_log_likelihood()
        assert_allclose(evolve_lnL, -103.05742415448259)

    def test_nucleotide(self):
        """test a nucleotide model."""
        submod = TimeReversibleNucleotide(
            equal_motif_probs=True, motif_probs=None, predicates={"kappa": "transition"}
        )
        # now do using the evolve
        likelihood_function = self._makeLikelihoodFunction(submod)
        self.assertEqual(likelihood_function.get_num_free_params(), 0)
        evolve_lnL = likelihood_function.get_log_likelihood()
        assert_allclose(evolve_lnL, -148.6455087258624)

    def test_solved_nucleotide(self):
        """test a solved nucleotide model."""
        submod = get_model("TN93", rate_matrix_required=False)
        # now do using the evolve
        lf = submod.make_likelihood_function(self.tree)
        lf.set_alignment(self.alignment)
        self.assertEqual(lf.get_num_free_params(), 5)
        lf.optimise(show_progress=False, max_evaluations=20, limit_action="ignore")
        self.assertTrue(lf.lnL > -152)

    def test_discrete_nucleotide(self):
        """test that partially discrete nucleotide model can be constructed,
        differs from continuous, and has the expected number of free params"""
        submod = TimeReversibleNucleotide(
            equal_motif_probs=True, motif_probs=None, predicates={"kappa": "transition"}
        )
        likelihood_function = self._makeLikelihoodFunction(
            submod, discrete_edges=["Human"]
        )
        self.assertEqual(likelihood_function.get_num_free_params(), 12)
        evolve_lnL = likelihood_function.get_log_likelihood()
        self.assertNotEqual(evolve_lnL, -157.49363874840455)

    def test_dinucleotide(self):
        """test a dinucleotide model."""
        submod = substitution_model.TimeReversibleDinucleotide(
            equal_motif_probs=True,
            predicates={"kappa": "transition"},
            mprob_model="tuple",
        )
        likelihood_function = self._makeLikelihoodFunction(submod)
        evolve_lnL = likelihood_function.get_log_likelihood()
        assert_allclose(evolve_lnL, -118.35045332768402)

    def test_protein(self):
        """test a protein model."""
        submod = substitution_model.TimeReversibleProtein(equal_motif_probs=True)

        likelihood_function = self._makeLikelihoodFunction(submod, translate=True)

        evolve_lnL = likelihood_function.get_log_likelihood()
        assert_allclose(evolve_lnL, -91.35162044257062)


class LikelihoodFunctionTests(TestCase):
    """tests for a tree analysis class. Various tests to create a tree analysis class,
    set parameters, and test various functions.
    """

    def setUp(self):
        self.submodel = TimeReversibleNucleotide(
            equal_motif_probs=True, predicates={"beta": "transition"}
        )

        self.data = load_aligned_seqs(
            filename=os.path.join(data_path, "brca1_5.paml"),
            moltype=self.submodel.moltype,
        )

        self.tree = load_tree(filename=os.path.join(data_path, "brca1_5.tree"))

    def _makeLikelihoodFunction(self, **kw):
        lf = self.submodel.make_likelihood_function(self.tree, **kw)
        lf.set_param_rule("beta", is_independent=True)
        lf.set_alignment(self.data)
        return lf

    def _setLengthsAndBetas(self, likelihood_function):
        for (species, length) in [
            ("DogFaced", 0.1),
            ("NineBande", 0.2),
            ("Human", 0.3),
            ("HowlerMon", 0.4),
            ("Mouse", 0.5),
        ]:
            likelihood_function.set_param_rule(
                "length", value=length, edge=species, is_constant=True
            )
        for (species1, species2, length) in [
            ("Human", "HowlerMon", 0.7),
            ("Human", "Mouse", 0.6),
        ]:
            LCA = self.tree.get_connecting_node(species1, species2).name
            likelihood_function.set_param_rule(
                "length", value=length, edge=LCA, is_constant=True
            )

        likelihood_function.set_param_rule("beta", value=4.0, is_constant=True)

    def test_information_criteria(self):
        """test get information criteria from a model."""
        lf = self._makeLikelihoodFunction()
        nfp = lf.get_num_free_params()
        lnL = lf.get_log_likelihood()
        l = len(self.data)
        assert_allclose(lf.get_aic(), aic(lnL, nfp))
        assert_allclose(lf.get_aic(second_order=True), aic(lnL, nfp, l))

        assert_allclose(lf.get_bic(), bic(lnL, nfp, l))

    def test_result_str(self):
        # actualy more a test of self._setLengthsAndBetas()
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        self.assertEqual(
            str(likelihood_function),
            """Likelihood function statistics
log-likelihood = -250.6867
number of free parameters = 0\n\
======
  beta
------
4.0000
------
=============================
edge         parent    length
-----------------------------
Human        edge.0    0.3000
HowlerMon    edge.0    0.4000
edge.0       edge.1    0.7000
Mouse        edge.1    0.5000
edge.1       root      0.6000
NineBande    root      0.2000
DogFaced     root      0.1000
-----------------------------
====================================
     A         C         G         T
------------------------------------
0.2500    0.2500    0.2500    0.2500
------------------------------------""",
        )

        likelihood_function = self._makeLikelihoodFunction(digits=2, space=2)
        self.assertEqual(
            str(likelihood_function),
            """Likelihood function statistics
log-likelihood = -382.5399
number of free parameters = 14
===============================
edge       parent  length  beta
-------------------------------
Human      edge.0    1.00  1.00
HowlerMon  edge.0    1.00  1.00
edge.0     edge.1    1.00  1.00
Mouse      edge.1    1.00  1.00
edge.1     root      1.00  1.00
NineBande  root      1.00  1.00
DogFaced   root      1.00  1.00
-------------------------------
======================
   A     C     G     T
----------------------
0.25  0.25  0.25  0.25
----------------------""",
        )

    def test_calclikelihood(self):
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        self.assertAlmostEqual(
            -250.686745262, likelihood_function.get_log_likelihood(), places=9
        )

    def test_g_statistic(self):
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        self.assertAlmostEqual(
            230.77670557, likelihood_function.get_G_statistic(), places=6
        )

    def test_ancestralsequences(self):
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        result = likelihood_function.reconstruct_ancestral_seqs()["edge.0"]
        a_column_with_mostly_Ts = -1
        motif_G = 2
        self.assertAlmostEqual(
            2.28460181711e-05, result[a_column_with_mostly_Ts][motif_G], places=8
        )
        lf = self.submodel.make_likelihood_function(self.tree, bins=["low", "high"])
        lf.set_param_rule("beta", bin="low", value=0.1)
        lf.set_param_rule("beta", bin="high", value=10.0)
        lf.set_alignment(self.data)
        result = lf.reconstruct_ancestral_seqs()

    def test_likely_ancestral(self):
        """excercising the most likely ancestral sequences"""
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        likelihood_function.likely_ancestral_seqs()

    def test_simulate_alignment(self):
        "Simulate DNA alignment"
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        simulated_alignment = likelihood_function.simulate_alignment(
            20, exclude_internal=False
        )
        self.assertEqual(len(simulated_alignment), 20)
        self.assertEqual(len(simulated_alignment.names), 8)

    def test_simulateHetergeneousAlignment(self):
        "Simulate substitution-heterogeneous DNA alignment"
        lf = self.submodel.make_likelihood_function(self.tree, bins=["low", "high"])
        lf.set_param_rule("beta", bin="low", value=0.1)
        lf.set_param_rule("beta", bin="high", value=10.0)
        lf.simulate_alignment(100)

    def test_simulatePatchyHetergeneousAlignment(self):
        "Simulate patchy substitution-heterogeneous DNA alignment"
        lf = self.submodel.make_likelihood_function(
            self.tree, bins=["low", "high"], sites_independent=False
        )
        lf.set_param_rule("beta", bin="low", value=0.1)
        lf.set_param_rule("beta", bin="high", value=10.0)
        lf.simulate_alignment(100)

    def test_simulate_alignment1(self):
        "Simulate alignment when no alignment set"
        al = make_aligned_seqs(data={"a": "ggaatt", "c": "cctaat"})
        t = make_tree("(a,c);")
        sm = get_model("F81")
        lf = sm.make_likelihood_function(t)
        # no provided alignment raises an exception
        with self.assertRaises(ValueError):
            lf.simulate_alignment()

        # unless you provide length
        sim_aln = lf.simulate_alignment(sequence_length=10)
        self.assertEqual(len(sim_aln), 10)

    def test_simulate_alignment2(self):
        "Simulate alignment with dinucleotide model"
        al = make_aligned_seqs(data={"a": "ggaatt", "c": "cctaat"})
        t = make_tree("(a,c);")
        sm = substitution_model.TimeReversibleDinucleotide(mprob_model="tuple")
        lf = sm.make_likelihood_function(t)
        lf.set_alignment(al)
        simalign = lf.simulate_alignment()
        self.assertEqual(len(simalign), 6)

    def test_simulate_alignment3(self):
        """Simulated alignment with gap-induced ambiguous positions
        preserved"""
        t = make_tree(treestring="(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1)root;")
        al = make_aligned_seqs(
            data={
                "a": "g--cactat?",
                "b": "---c-ctcct",
                "c": "-a-c-ctat-",
                "d": "-a-c-ctat-",
            }
        )
        sm = TimeReversibleNucleotide(recode_gaps=True)
        lf = sm.make_likelihood_function(t)
        # pc.set_constant_lengths()
        lf.set_alignment(al)
        # print lf.simulate_alignment(sequence_length=10)
        simulated = lf.simulate_alignment()
        self.assertEqual(len(simulated.names), 4)
        import re

        self.assertEqual(re.sub("[ATCG]", "x", simulated.to_dict()["a"]), "x??xxxxxx?")

    def test_simulate_alignment_root_sequence(self):
        """provide a root sequence for simulating an alignment"""

        def use_root_seq(root_sequence):
            al = make_aligned_seqs(data={"a": "ggaatt", "c": "cctaat"})
            t = make_tree(treestring="(a,c);")
            sm = substitution_model.TimeReversibleDinucleotide(mprob_model="tuple")
            lf = sm.make_likelihood_function(t)
            lf.set_alignment(al)
            simalign = lf.simulate_alignment(
                exclude_internal=False, root_sequence=root_sequence
            )
            root = simalign.named_seqs["root"]
            self.assertEqual(str(root), str(root_sequence))

        root_sequence = DNA.make_seq("GTAATT")
        use_root_seq(root_sequence)  # as a sequence instance
        use_root_seq("GTAATC")  # as a string

    def test_pc_initial_parameters(self):
        """Default parameter values from original annotated tree"""
        likelihood_function = self._makeLikelihoodFunction()
        self._setLengthsAndBetas(likelihood_function)
        tree = likelihood_function.get_annotated_tree()
        lf = self.submodel.make_likelihood_function(tree)
        lf.set_alignment(self.data)
        self.assertEqual(lf.get_param_value("length", "Human"), 0.3)
        self.assertEqual(lf.get_param_value("beta", "Human"), 4.0)

    def test_set_par_all(self):
        likelihood_function = self._makeLikelihoodFunction()
        likelihood_function.set_param_rule("length", value=4.0, is_constant=True)
        likelihood_function.set_param_rule("beta", value=6.0, is_constant=True)
        self.assertEqual(
            str(likelihood_function),
            """Likelihood function statistics
log-likelihood = -413.1886
number of free parameters = 0
======
  beta
------
6.0000
------
=============================
edge         parent    length
-----------------------------
Human        edge.0    4.0000
HowlerMon    edge.0    4.0000
edge.0       edge.1    4.0000
Mouse        edge.1    4.0000
edge.1       root      4.0000
NineBande    root      4.0000
DogFaced     root      4.0000
-----------------------------
====================================
     A         C         G         T
------------------------------------
0.2500    0.2500    0.2500    0.2500
------------------------------------""",
        )

    def test_set_param_rule_adjust_bounds(self):
        """check behaviour when modify bound and reset param rule"""
        lf = self._makeLikelihoodFunction()
        lf.set_param_rule(
            "beta", init=4.0, is_independent=True, edges=["DogFaced", "NineBande"]
        )
        lf.set_param_rule("beta", upper=2)
        val = lf.get_param_value("beta", edge="DogFaced")
        self.assertLess(val, 4)  # it will be the average of default and set values

    def test_get_motif_probs(self):
        likelihood_function = self._makeLikelihoodFunction()
        mprobs = likelihood_function.get_motif_probs()
        assert hasattr(mprobs, "keys"), mprobs
        keys = list(mprobs.keys())
        keys.sort()
        obs = self.submodel.get_motifs()
        obs.sort()
        self.assertEqual(obs, keys)

    def test_get_annotated_tree(self):
        lf = self._makeLikelihoodFunction()
        lf.set_param_rule("length", value=4.0, edge="Human", is_constant=True)
        lf.set_param_rule("beta", value=2.0, edge="Human", is_constant=True)

        result = lf.get_annotated_tree()
        human = result.get_node_matching_name("Human")
        self.assertEqual(human.params["length"], 4.0)
        self.assertEqual(human.length, 4.0)
        # specify length as paralinear or ENS does not fail
        ens = lf.get_annotated_tree(length_as="ENS")
        # now check correctly decorate with paralin
        plin = lf.get_annotated_tree(length_as="paralinear")
        plin_metric = lf.get_paralinear_metric()
        human = plin.get_node_matching_name("Human")
        assert_allclose(human.length, plin_metric["Human"])

    def test_getparamsasdict(self):
        likelihood_function = self._makeLikelihoodFunction()
        likelihood_function.set_name("TEST")
        self.assertEqual(
            str(likelihood_function),
            """TEST
log-likelihood = -382.5399
number of free parameters = 14
=======================================
edge         parent    length      beta
---------------------------------------
Human        edge.0    1.0000    1.0000
HowlerMon    edge.0    1.0000    1.0000
edge.0       edge.1    1.0000    1.0000
Mouse        edge.1    1.0000    1.0000
edge.1       root      1.0000    1.0000
NineBande    root      1.0000    1.0000
DogFaced     root      1.0000    1.0000
---------------------------------------
====================================
     A         C         G         T
------------------------------------
0.2500    0.2500    0.2500    0.2500
------------------------------------""",
        )
        self.assertEqual(
            likelihood_function.get_param_value_dict(["edge"]),
            {
                "beta": {
                    "NineBande": 1.0,
                    "edge.1": 1.0,
                    "DogFaced": 1.0,
                    "Human": 1.0,
                    "edge.0": 1.0,
                    "Mouse": 1.0,
                    "HowlerMon": 1.0,
                },
                "length": {
                    "NineBande": 1.0,
                    "edge.1": 1.0,
                    "DogFaced": 1.0,
                    "Human": 1.0,
                    "edge.0": 1.0,
                    "Mouse": 1.0,
                    "HowlerMon": 1.0,
                },
            },
        )

    def test_get_statistics_from_empirical_model(self):
        """should return valid dict from an empirical substitution model"""
        submod = JTT92()
        aln = self.data.get_translation()

        lf = submod.make_likelihood_function(self.tree)
        lf.set_alignment(aln)
        stats = lf.get_param_value_dict(["edge"], params=["length"])

    def test_constant_to_free(self):
        """excercise setting a constant param rule, then freeing it"""
        # checks by just trying to make the calculator
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_alignment(self.data)
        lf.set_param_rule(
            "beta",
            is_constant=True,
            value=2.0,
            edges=["NineBande", "DogFaced"],
            clade=True,
        )
        lf.set_param_rule(
            "beta",
            init=2.0,
            is_constant=False,
            edges=["NineBande", "DogFaced"],
            clade=True,
        )

    def test_get_psub_rate_matrix(self):
        """lf should return consistent rate matrix and psub"""
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_alignment(self.data)
        edge = "NineBande"
        length = 0.5
        lf.set_param_rule("length", edge=edge, init=length)
        Q = lf.get_rate_matrix_for_edge("NineBande", calibrated=False)
        Q2 = lf.get_rate_matrix_for_edge("NineBande", calibrated=True)
        P = lf.get_psub_for_edge("NineBande")
        assert_allclose(expm(Q.array)(1.0), P.array)
        assert_allclose(expm(Q2.array)(length), P.array)

        # should fail for a discrete Markov model
        dm = ns_substitution_model.DiscreteSubstitutionModel(DNA.alphabet)
        lf = dm.make_likelihood_function(self.tree)
        lf.set_alignment(self.data)
        self.assertRaises(Exception, lf.get_rate_matrix_for_edge, "NineBande")

    def test_get_all_psubs_discrete(self):
        """should work for discrete time models"""
        sm = get_model("BH")
        lf = sm.make_likelihood_function(self.tree)
        lf.set_alignment(self.data)
        lf.get_all_psubs()

    def test_get_all_rate_matrices(self):
        """return matrices when just a pair"""
        aln = make_aligned_seqs(
            data={
                "Human": "GGCCTCCTGCGCTCCCTGGCCCGCCACCAG",
                "Opossum": "GGCTCCCTGCGCTCCCTTTCCCGCCGCCGG",
            },
            moltype="dna",
        )
        tree = make_tree(tip_names=aln.names)
        sm = get_model("HKY85")
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        Qs = lf.get_all_rate_matrices(calibrated=False)
        self.assertEqual(len(Qs), 2)

    def test_get_p_q_sitehet_model(self):
        """exercising get psub in phylohmm model"""
        from cogent3.maths.matrix_exponentiation import PadeExponentiator

        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf1 = sm.make_likelihood_function(self.tree, bins=4, sites_independent=False)
        lf1.set_alignment(self.data)
        Qs = lf1.get_all_rate_matrices(calibrated=False)
        Ps = lf1.get_all_psubs()
        self.assertEqual(len(Ps), len(Qs))
        self.assertEqual(set(Ps), set(Qs))
        for key, P in Ps.items():
            Pcomp = PadeExponentiator(Qs[key].to_array())()
            assert_allclose(Pcomp, P.to_array())

    def test_all_rate_matrices_unique(self):
        """exercising this code"""
        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf1 = sm.make_likelihood_function(self.tree, bins=4)
        lf1.set_param_rule("kappa", init=4)
        lf1.set_alignment(self.data)
        lf1.all_rate_matrices_unique()

    def test_make_discrete_markov(self):
        """lf ignores tree lengths if a discrete Markov model"""
        t = make_tree(treestring="(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1)root;")
        dm = ns_substitution_model.DiscreteSubstitutionModel(DNA.alphabet)
        dm.make_likelihood_function(t)

    def test_exercise_set_align(self):
        "lf.set_align should work for different models"
        al = make_aligned_seqs(data={"a": "ggaatt", "c": "cctaat"})
        t = make_tree(treestring="(a,c);")
        for klass in [CNFGTR, Y98, MG94HKY]:
            sm = klass()
            lf = sm.make_likelihood_function(t)
            lf.set_alignment(al)

    def test_get_param_rules(self):
        """correctly return rules that can be used to reconstruct a lf"""
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_alignment(self.data)
        lf.set_param_rule("beta", init=2.0)
        lf.set_param_rule(
            "beta", value=2.0, edges=["Human", "HowlerMon"], is_constant=True
        )
        lf.set_param_rule("length", init=0.5, edges="Human", upper=5)
        lf.set_param_rule("length", value=0.25, edges="HowlerMon", is_constant=True)
        lnL = lf.get_log_likelihood()
        rules = lf.get_param_rules()
        for rule in rules:
            if rule["par_name"] == "length":
                if rule["edge"] == "Human":
                    self.assertEqual(rule["upper"], 5)
                elif rule["edge"] == "HowlerMon":
                    self.assertTrue(rule.get("is_constant", False))
            elif rule["par_name"] == "mprobs":
                self.assertEqual(rule["value"], {b: 0.25 for b in "ACGT"})

        self.assertEqual(len(rules), 10)
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_alignment(self.data)
        with lf.updates_postponed():
            for rule in rules:
                lf.set_param_rule(**rule)
        new_lnL = lf.get_log_likelihood()
        assert_allclose(new_lnL, lnL)

    def test_get_param_rules_multilocus(self):
        """correctly return rules from multilocus lf"""
        data = load_aligned_seqs(
            filename=os.path.join(os.getcwd(), "data", "brca1_5.paml")
        )
        half = len(data) // 2
        aln1 = data[:half]
        aln2 = data[half:]
        loci_names = ["1st-half", "2nd-half"]
        loci = [aln1, aln2]
        tree = make_tree(tip_names=data.names)
        model = get_model("HKY85", optimise_motif_probs=True)
        lf = model.make_likelihood_function(tree, loci=loci_names)
        lf.set_alignment(loci)
        lf.set_param_rule("mprobs", is_independent=False)
        rules = lf.get_param_rules()
        lf2 = model.make_likelihood_function(tree, loci=loci_names)
        lf2.set_alignment(loci)
        lf2.apply_param_rules(rules=rules)
        assert_allclose(lf.lnL, lf2.lnL)

    def test_get_param_rules_discrete(self):
        """discrete time models produce valid rules"""
        sm = get_model("BH")
        aln = self.data.take_seqs(self.data.names[:3])
        tree = self.tree.get_sub_tree(aln.names)
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        lf.optimise(max_evaluations=100, limit_action="ignore")
        rules = lf.get_param_rules()

        new_lf = sm.make_likelihood_function(tree)
        new_lf.set_alignment(aln)
        new_lf.apply_param_rules(rules)
        assert_allclose(new_lf.get_motif_probs().array, lf.get_motif_probs().array)
        for edge in tree.preorder():
            if edge.is_root():
                continue
            orig_p = lf.get_psub_for_edge(edge.name)
            new_p = new_lf.get_psub_for_edge(edge.name)
            assert_allclose(new_p.array, orig_p.array, err_msg=edge.name, atol=1e-5)
        assert_allclose(lf.lnL, new_lf.lnL, atol=1e-4)

    def test_get_param_rules_constrained(self):
        """correctly return rules that reconstruct a lf with constrained length"""
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_alignment(self.data)
        lf.set_param_rule("beta", init=2.0)
        lf.set_param_rule(
            "beta", value=2.0, edges=["Human", "HowlerMon"], is_constant=True
        )
        lf.set_param_rule("length", init=0.5, is_independent=False)
        rules = lf.get_param_rules()
        new = self.submodel.make_likelihood_function(self.tree)
        new.apply_param_rules(rules)
        self.assertEqual(new.nfp, lf.nfp)

    def test_apply_param_rules(self):
        """successfully apply a set of parameter rules"""
        lf = self.submodel.make_likelihood_function(self.tree)
        nfp = lf.get_num_free_params()
        rules = [
            dict(par_name="beta", edges=["Human", "HowlerMon"], init=2),
            dict(par_name="beta", edges=["NineBande", "DogFaced"], init=4),
        ]
        lf.apply_param_rules(rules)
        self.assertEqual(lf.get_num_free_params(), nfp + 2)

        lf = self.submodel.make_likelihood_function(self.tree)
        rules = [
            dict(par_name="beta", edges=["Human", "HowlerMon"], init=2),
            dict(
                par_name="beta",
                edges=["NineBande", "DogFaced"],
                init=4,
                is_independent=True,
            ),
        ]
        lf.apply_param_rules(rules)
        self.assertEqual(lf.get_num_free_params(), nfp + 3)

        lf = self.submodel.make_likelihood_function(self.tree)
        rules = [
            dict(par_name="beta", edges=["Human", "HowlerMon"], init=2),
            dict(
                par_name="beta",
                edges=["NineBande", "DogFaced"],
                value=4,
                is_constant=True,
            ),
        ]
        lf.apply_param_rules(rules)
        self.assertEqual(lf.get_num_free_params(), nfp + 1)

    def test_get_apply_param_rules_site_het_models(self):
        """correctly use and apply param rules from site-het and phyloHMM models"""
        with open("data/site-het-param-rules.json") as infile:
            rules = json.load(infile)
        aln = load_aligned_seqs("data/primates_brca1.fasta", moltype="dna")
        tree = load_tree("data/primates_brca1.tree")
        # gamma distributed length
        rule_lnL = rules.pop("gamma-length")
        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf1 = sm.make_likelihood_function(tree, bins=4)
        lf1.set_alignment(aln)
        lf1.apply_param_rules(rule_lnL["rules"])
        assert_allclose(lf1.lnL, rule_lnL["lnL"])
        lf2 = sm.make_likelihood_function(tree, bins=4)
        lf2.set_alignment(aln)
        lf2.apply_param_rules(lf1.get_param_rules())
        assert_allclose(lf2.lnL, lf1.lnL)
        # gamma distributed kappa
        rule_lnL = rules.pop("gamma-kappa")
        sm = get_model("HKY85", ordered_param="kappa", distribution="gamma")
        lf1 = sm.make_likelihood_function(tree, bins=4)
        lf1.set_alignment(aln)
        lf1.apply_param_rules(rule_lnL["rules"])
        assert_allclose(lf1.lnL, rule_lnL["lnL"])
        lf2 = sm.make_likelihood_function(tree, bins=4)
        lf2.set_alignment(aln)
        lf2.apply_param_rules(lf1.get_param_rules())
        assert_allclose(lf2.lnL, lf1.lnL)
        # free kappa
        rule_lnL = rules.pop("free-kappa")
        sm = get_model("HKY85")
        lf1 = sm.make_likelihood_function(tree, bins=["slow", "fast"])
        lf1.set_alignment(aln)
        lf1.apply_param_rules(rule_lnL["rules"])
        assert_allclose(lf1.lnL, rule_lnL["lnL"])
        lf2 = sm.make_likelihood_function(tree, bins=["slow", "fast"])
        lf2.set_alignment(aln)
        lf2.apply_param_rules(lf1.get_param_rules())
        assert_allclose(lf2.lnL, lf1.lnL)
        # phylo-hmm kappa
        rule_lnL = rules.pop("phylohmm-gamma-kappa")
        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf1 = sm.make_likelihood_function(tree, bins=4, sites_independent=False)
        lf1.set_alignment(aln)
        lf1.apply_param_rules(rule_lnL["rules"])
        assert_allclose(lf1.lnL, rule_lnL["lnL"])
        lf2 = sm.make_likelihood_function(tree, bins=4, sites_independent=False)
        lf2.set_alignment(aln)
        lf2.apply_param_rules(lf1.get_param_rules())
        assert_allclose(lf2.lnL, lf1.lnL)

    def test_set_time_heterogeneity_partial(self):
        """correctly apply partial time heterogeneity of rate terms"""
        lf = self.submodel.make_likelihood_function(self.tree)
        nfp0 = lf.nfp
        edge_sets = [dict(edges=["Human", "HowlerMon"])]
        lf.set_time_heterogeneity(edge_sets=edge_sets, is_independent=False)
        nfp1 = lf.nfp
        self.assertEqual(nfp1 - nfp0, 1)
        lf.set_time_heterogeneity(edge_sets=edge_sets, is_independent=True)
        nfp2 = lf.nfp
        self.assertEqual(nfp2 - nfp1, 1)

    def test_set_time_heterogeneity_multilocus(self):
        """apply time heterogeneity for multilocus function"""
        half = len(self.data) // 2
        aln1 = self.data[:half]
        aln2 = self.data[half:]
        loci_names = ["1st-half", "2nd-half"]
        loci = [aln1, aln2]
        model = get_model("GN", optimise_motif_probs=True)
        # should not fail
        lf = model.make_likelihood_function(self.tree, loci=loci_names)
        assert lf.locus_names == loci_names
        lf.set_alignment(loci)
        edges = ["Human", "HowlerMon"]
        lf = model.make_likelihood_function(
            self.tree,
            loci=loci_names,
            discrete_edges=edges,
        )
        lf.set_time_heterogeneity(upper=100, is_independent=True)
        lf.set_alignment(loci)
        lf.optimise(max_evaluations=10, limit_action="ignore", show_progress=False)
        stats = lf.get_statistics()
        timehet_edge_names = set(
            n for n in self.tree.get_node_names(includeself=False) if n not in edges
        )
        for t in stats:
            if t.title == "edge locus params":
                assert set(t.columns["edge"]) == timehet_edge_names

    def test_getting_pprobs(self):
        """posterior bin probs same length as aln for phylo-HMM model"""
        with open("data/site-het-param-rules.json") as infile:
            rules = json.load(infile)

        aln = load_aligned_seqs("data/primates_brca1.fasta", moltype="dna")
        tree = load_tree("data/primates_brca1.tree")
        rule_lnL = rules.pop("gamma-length")
        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf = sm.make_likelihood_function(tree, bins=4, sites_independent=True)
        lf.set_alignment(aln)
        lf.apply_param_rules(rule_lnL["rules"])
        bprobs = lf.get_bin_probs()
        self.assertEqual(bprobs.shape[1], len(aln))

    def test_bin_probs(self):
        """posterior bin probs same length as aln for rate-het model"""
        aln = load_aligned_seqs("data/primates_brca1.fasta", moltype="dna")
        tree = load_tree("data/primates_brca1.tree")
        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf = sm.make_likelihood_function(tree, bins=4, sites_independent=False)
        lf.set_alignment(aln)
        bprobs = lf.get_bin_probs()
        self.assertEqual(bprobs.shape[1], len(aln))

    def test_time_het_init_from_nested(self):
        """initialise from nested should honour alt model setting"""
        # setting time-het for entire Q
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        gn = GN(optimise_motif_probs=True)
        null = gn.make_likelihood_function(tree)
        null.set_alignment(_aln)
        edge_sets = [dict(edges=["Human", "Mouse"])]
        null.set_time_heterogeneity(edge_sets=edge_sets, is_independent=False)
        nfp_null = null.nfp
        alt = gn.make_likelihood_function(tree)
        alt.set_alignment(_aln)
        alt.set_time_heterogeneity(is_independent=True)
        nfp_alt_0 = alt.nfp
        self.assertEqual(nfp_alt_0 - nfp_null, 11)
        alt.initialise_from_nested(null)
        nfp_alt_1 = alt.nfp
        self.assertEqual(nfp_alt_1 - nfp_null, 11)
        edge_sets = [dict(edges=("Human", "Mouse"))]
        null.set_time_heterogeneity(edge_sets=edge_sets, is_independent=False)

    def test_init_from_nested_genstat(self):
        """initialising a general stationary model from a nested time-reversible model works"""
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        gtr = get_model("GTR")
        gs = GeneralStationary(gtr.alphabet)
        gtr_lf = gtr.make_likelihood_function(tree)
        gtr_lf.set_alignment(_aln)
        mprobs = dict(A=0.1, T=0.2, C=0.3, G=0.4)
        gtr_lf.set_motif_probs(mprobs)
        rate_params = {"A/C": 0.75, "A/G": 3, "A/T": 1.5, "C/G": 0.2, "C/T": 6}
        for par_name, val in rate_params.items():
            gtr_lf.set_param_rule(par_name, init=val)

        gs_lf = gs.make_likelihood_function(tree)
        gs_lf.set_alignment(_aln)
        gs_lf.initialise_from_nested(gtr_lf)
        assert_allclose(gs_lf.lnL, gtr_lf.lnL)

    def test_set_time_heterogeneity(self):
        """correctly apply time heterogeneity of rate terms"""
        lf = self.submodel.make_likelihood_function(self.tree)
        nfp = lf.get_num_free_params()
        # cannot exclude a param that isn't part of the model
        with self.assertRaises(ValueError):
            lf.set_time_heterogeneity(is_independent=True, exclude_params="omega")

        # cannot specify is_constant and init
        with self.assertRaises(ValueError):
            lf.set_time_heterogeneity(is_constant=True, init=2)

        # we should be able to just make the entire lf time-heterogeneous
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_time_heterogeneity(is_independent=True)
        got = lf.get_num_free_params()
        self.assertEqual(got, nfp + len(lf.tree.get_node_names(includeself=False)) - 1)

        # we should be able to specify a set of edges that get treated as a block
        # if not specified, the edges are considered to be not-independent
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_time_heterogeneity(
            edge_sets=dict(edges=["Human", "HowlerMon"]), is_independent=False
        )
        got = lf.get_num_free_params()
        self.assertEqual(got, nfp + 1)

        # making them independent
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_time_heterogeneity(
            edge_sets=dict(edges=["Human", "HowlerMon"]), is_independent=True
        )
        got = lf.get_num_free_params()
        self.assertEqual(got, nfp + 2)

        # making them constant
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_time_heterogeneity(
            edge_sets=dict(edges=["Human", "HowlerMon"]), is_constant=True
        )
        got = lf.get_num_free_params()
        self.assertEqual(got, nfp)

        # providing other settings within the edge_set
        lf = self.submodel.make_likelihood_function(self.tree)
        lf.set_time_heterogeneity(
            edge_sets=[
                dict(edges=["Human", "HowlerMon"], init=3),
                dict(edges=["NineBande", "DogFaced"], value=5, is_constant=True),
            ]
        )
        got = lf.get_num_free_params()
        self.assertEqual(got, nfp + 1)
        for edge, exp in [("Human", 3), ("NineBande", 5)]:
            got = lf.get_param_value("beta", edge=edge)
            self.assertEqual(got, exp)

    def test_initialise_from_nested_diff_scoped(self):
        """non-reversible likelihood initialised from nested, scoped, time-reversible"""
        mprobs = {b: p for b, p in zip(DNA, [0.1, 0.2, 0.3, 0.4])}
        rate_params = {"A/C": 2.0, "A/G": 3.0, "A/T": 4.0, "C/G": 5.0, "C/T": 6.0}

        rate_params1 = {"A/C": 4.0, "A/G": 6.0, "C/T": 3.0}

        simple = GTR()
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree, digits=2)
        slf.set_alignment(_aln)
        slf.set_name("GTR")
        slf.set_motif_probs(mprobs)
        for param, val in rate_params.items():
            slf.set_param_rule(param, init=val)

        for param, val in rate_params1.items():
            slf.set_param_rule(param, init=val, edges=["Human"])

        lengths = {e: v for e, v in zip(tree.get_tip_names(), (0.2, 0.4, 0.1))}
        for e, val in lengths.items():
            slf.set_param_rule("length", edge=e, init=val)

        rich = GN(optimise_motif_probs=True)
        glf = rich.make_likelihood_function(tree, digits=2)
        glf.set_alignment(_aln)
        glf.set_name("GN")
        glf.initialise_from_nested(slf)
        assert_allclose(glf.get_log_likelihood(), slf.get_log_likelihood())

    def test_initialise_from_nested_diff(self):
        """non-reversible likelihood initialised from nested, non-scoped, time-reversible"""
        mprobs = {b: p for b, p in zip(DNA, [0.1, 0.2, 0.3, 0.4])}
        rate_params = {"A/C": 2.0, "A/G": 3.0, "A/T": 4.0, "C/G": 5.0, "C/T": 6.0}

        simple = GTR()
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree, digits=2)
        slf.set_alignment(_aln)
        slf.set_name("GTR")
        slf.set_motif_probs(mprobs)
        for param, val in rate_params.items():
            slf.set_param_rule(param, init=val)
        lengths = {e: v for e, v in zip(tree.get_tip_names(), (0.2, 0.4, 0.1))}
        for e, val in lengths.items():
            slf.set_param_rule("length", edge=e, init=val)

        # set mprobs and then set the rate terms
        rich = GN(optimise_motif_probs=True)
        glf = rich.make_likelihood_function(tree, digits=2)
        glf.set_alignment(_aln)
        glf.set_name("GN")
        glf.initialise_from_nested(slf)
        assert_allclose(glf.get_log_likelihood(), slf.get_log_likelihood())

    def test_initialise_from_nested_diff_stat(self):
        """non-reversible stationary initialised from nested time-reversible"""
        mprobs = {b: p for b, p in zip(DNA, [0.1, 0.2, 0.3, 0.4])}
        rate_params = {"A/C": 2.0, "A/G": 3.0, "A/T": 4.0, "C/G": 5.0, "C/T": 6.0}

        simple = GTR()
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree, digits=2)
        slf.set_alignment(_aln)
        slf.set_name("GTR")
        slf.set_motif_probs(mprobs)
        for param, val in rate_params.items():
            slf.set_param_rule(param, init=val)
        lengths = {e: v for e, v in zip(tree.get_tip_names(), (0.2, 0.4, 0.1))}
        for e, val in lengths.items():
            slf.set_param_rule("length", edge=e, init=val)

        # set mprobs and then set the rate terms
        from cogent3.evolve.ns_substitution_model import GeneralStationary

        rich = GeneralStationary(DNA.alphabet)
        glf = rich.make_likelihood_function(tree, digits=2)
        glf.set_alignment(_aln)
        glf.set_name("GSN")
        glf.initialise_from_nested(slf)
        assert_allclose(glf.get_log_likelihood(), slf.get_log_likelihood())

    def test_initialise_from_nested_same_type_tr(self):
        """time-reversible likelihood initialised from nested, non-scoped, time-reversible"""
        mprobs = {b: p for b, p in zip(DNA, [0.1, 0.2, 0.3, 0.4])}
        rate_params = {"kappa": 6}
        simple = HKY85()
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree, digits=2)
        slf.set_alignment(_aln)
        slf.set_name("HKY85")
        slf.set_motif_probs(mprobs)
        for param, val in rate_params.items():
            slf.set_param_rule(param, init=val)

        lengths = {e: v for e, v in zip(tree.get_tip_names(), (0.2, 0.4, 0.1))}
        for e, val in lengths.items():
            slf.set_param_rule("length", edge=e, init=val)

        rich = GTR()
        glf = rich.make_likelihood_function(tree, digits=2)
        glf.set_alignment(_aln)
        glf.set_name("GTR")
        glf.initialise_from_nested(slf)
        assert_allclose(glf.get_log_likelihood(), slf.get_log_likelihood())

    def test_initialise_from_nested_same_type_tr_scoped(self):
        """time-reversible likelihood initialised from nested, scoped, time-reversible"""
        mprobs = {b: p for b, p in zip(DNA, [0.1, 0.2, 0.3, 0.4])}
        rate_params = {"kappa": 6}
        rate_params1 = {"kappa": 3}
        simple = HKY85()
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree, digits=2)
        slf.set_alignment(_aln)
        slf.set_name("HKY85")
        slf.set_motif_probs(mprobs)
        for param, val in rate_params.items():
            slf.set_param_rule(param, init=val)
        for param, val in rate_params1.items():
            slf.set_param_rule(param, init=val, edges=["Human"])

        lengths = {e: v for e, v in zip(tree.get_tip_names(), (0.2, 0.4, 0.1))}
        for e, val in lengths.items():
            slf.set_param_rule("length", edge=e, init=val)

        rich = GTR()
        glf = rich.make_likelihood_function(tree, digits=2)
        glf.set_alignment(_aln)
        glf.set_name("GTR")
        glf.initialise_from_nested(slf)
        assert_allclose(glf.get_log_likelihood(), slf.get_log_likelihood())

    def test_initialise_from_nested_same_type_nr(self):
        """non-reversible likelihood initialised from nested, non-scoped, non-reversible"""
        mprobs = {b: p for b, p in zip(DNA, [0.1, 0.2, 0.3, 0.4])}
        rate_params = {
            "(A>G | T>C)": 5,
            "(A>T | T>A)": 4,
            "(C>G | G>C)": 3,
            "(C>T | G>A)": 2,
            "(G>T | C>A)": 1,
        }

        simple = ssGN(optimise_motif_probs=True)
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree, digits=2)
        slf.set_alignment(_aln)
        slf.set_name("ssGN")
        slf.set_motif_probs(mprobs)
        for param, val in rate_params.items():
            slf.set_param_rule(param, init=val)

        lengths = {e: v for e, v in zip(tree.get_tip_names(), (0.2, 0.4, 0.1))}
        for e, val in lengths.items():
            slf.set_param_rule("length", edge=e, init=val)

        rich = GN(optimise_motif_probs=True)
        glf = rich.make_likelihood_function(tree, digits=2)
        glf.set_alignment(_aln)
        glf.set_name("GN")
        glf.initialise_from_nested(slf)
        expect = slf.get_log_likelihood()
        got = glf.get_log_likelihood()
        assert_allclose(got, expect)

    def test_initialise_from_nested_same_type_nr_scoped(self):
        """non-reversible likelihood initialised from nested, scoped, non-reversible"""
        mprobs = {b: p for b, p in zip(DNA, [0.1, 0.2, 0.3, 0.4])}
        rate_params = {
            "(A>G | T>C)": 5,
            "(A>T | T>A)": 4,
            "(C>G | G>C)": 3,
            "(C>T | G>A)": 2,
            "(G>T | C>A)": 1,
        }
        rate_params1 = {"(A>G | T>C)": 2, "(A>T | T>A)": 6, "(C>G | G>C)": 1}

        simple = ssGN(optimise_motif_probs=True)
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree, digits=2)
        slf.set_alignment(_aln)
        slf.set_name("ssGN")
        slf.set_motif_probs(mprobs)
        for param, val in rate_params.items():
            slf.set_param_rule(param, init=val)

        for param, val in rate_params1.items():
            slf.set_param_rule(param, init=val, edges=["Human"])

        lengths = {e: v for e, v in zip(tree.get_tip_names(), (0.2, 0.4, 0.1))}
        for e, val in lengths.items():
            slf.set_param_rule("length", edge=e, init=val)

        rich = GN(optimise_motif_probs=True)
        glf = rich.make_likelihood_function(tree, digits=2)
        glf.set_alignment(_aln)
        glf.set_name("GN")
        glf.initialise_from_nested(slf)
        assert_allclose(glf.get_log_likelihood(), slf.get_log_likelihood())

    def test_initialise_from_nested_codon_scoped(self):
        """scoped non-reversible likelihood initialised from nested scoped, non-reversible"""
        simple = get_model("H04GK")
        tree = make_tree(tip_names=["Human", "Mouse", "Opossum"])
        slf = simple.make_likelihood_function(tree)
        slf.set_alignment(_aln)
        slf.set_time_heterogeneity(
            edge_sets=[
                dict(edges=["Opossum"], is_independent=True),
            ],
            exclude_params=["kappa", "omega"],
        )
        slf.optimise(max_evaluations=50, limit_action="ignore", show_progress=False)
        glf = simple.make_likelihood_function(tree)
        glf.set_alignment(_aln)
        glf.set_time_heterogeneity(
            edge_sets=[
                dict(edges=["Opossum"], is_independent=True),
                dict(edges=["Human", "Mouse"], is_independent=True),
            ],
            exclude_params=["kappa", "omega"],
        )
        glf.initialise_from_nested(slf)
        assert_allclose(glf.lnL, slf.lnL)

    def test_get_lengths_as_ens_equal(self):
        """lengths equals ENS for a time-reversible model"""
        moprobs = numpy.array([0.1, 0.2, 0.3, 0.4])
        length = 0.1
        lf = HKY85().make_likelihood_function(make_tree(tip_names=["a", "b", "c"]))
        lf.set_motif_probs(moprobs)
        lf.set_param_rule("kappa", init=1)
        lf.set_param_rule("length", edge="a", init=length)
        len_dict = lf.get_lengths_as_ens()
        assert_allclose(len_dict["a"], length)

    def test_get_lengths_as_ens_not_equal(self):
        """lengths do not equal ENS for a non-reversible model"""
        moprobs = numpy.array([0.1, 0.2, 0.3, 0.4])
        length = 0.1
        lf = GN().make_likelihood_function(make_tree(tip_names=["a", "b", "c"]))
        lf.set_motif_probs(moprobs)
        lf.set_param_rule("length", init=length)
        # setting arbitrary values for GN rate terms
        init = 0.1
        for par_name in lf.model.get_param_list():
            lf.set_param_rule(par_name, init=init)
            init += 0.1

        len_dict = lf.get_lengths_as_ens()
        self.assertNotAlmostEqual(len_dict["a"], length)

    def test_bin_probs(self):
        """bin probs has same length as alignment for monomer alphabet"""
        aln = load_aligned_seqs("data/primates_brca1.fasta", moltype="dna")
        tree = load_tree("data/primates_brca1.tree")
        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf = sm.make_likelihood_function(tree, bins=4, sites_independent=False)
        lf.set_alignment(aln)
        bprobs = lf.get_bin_probs()
        self.assertEqual(bprobs.shape[1], len(aln))

    def test_get_paralinear(self):
        """returns correct paralinear from a lf"""
        from cogent3.maths.measure import paralinear_continuous_time

        aln = load_aligned_seqs("data/primates_brca1.fasta", moltype="dna")
        tree = load_tree("data/primates_brca1.tree")
        names = ["TreeShrew", "Mouse", "Rhesus", "Orangutan", "Human"]
        tree = tree.get_sub_tree(names)
        aln = aln.take_seqs(names)
        sm = get_model("GN")
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        rules = [
            {"par_name": "G>T", "init": 0.6836},
            {"par_name": "G>C", "init": 0.9666},
            {"par_name": "G>A", "init": 4.207},
            {"par_name": "T>C", "init": 3.419},
            {"par_name": "T>A", "init": 1.06},
            {"par_name": "C>G", "init": 1.04},
            {"par_name": "C>T", "init": 2.95},
            {"par_name": "C>A", "init": 0.739},
            {"par_name": "A>G", "init": 3.49},
            {"par_name": "A>T", "init": 0.835},
            {"par_name": "A>C", "init": 1.15},
            {
                "par_name": "mprobs",
                "init": {
                    "T": 0.2416726186347738,
                    "C": 0.1696848907070866,
                    "A": 0.38111966812998566,
                    "G": 0.2075228225281539,
                },
            },
            {"par_name": "length", "edge": "Human", "init": 0.009},
            {"par_name": "length", "edge": "Mouse", "init": 0.292},
            {"par_name": "length", "edge": "Orangutan", "init": 0.008},
            {"par_name": "length", "edge": "Rhesus", "init": 0.025},
            {"par_name": "length", "edge": "TreeShrew", "init": 0.12},
            {"par_name": "length", "edge": "edge.3", "init": 0.011},
            {"par_name": "length", "edge": "edge.4", "init": 0.056},
        ]
        lf.apply_param_rules(rules)
        rhesus = "Rhesus"
        # hand calc a value
        parent = lf.tree.get_node_matching_name(rhesus).parent.name
        all_mprobs = lf.get_motif_probs_by_node()
        pi = all_mprobs[parent]
        P = lf.get_psub_for_edge(rhesus)
        Q = lf.get_rate_matrix_for_edge(rhesus, calibrated=False)
        pl = paralinear_continuous_time(P.array, pi.array, Q.array)
        # check against lf computed one
        plinear = lf.get_paralinear_metric()
        assert_allclose(plinear[rhesus], pl)

    def test_to_rich_dict(self):
        """lf's from different substitution model classes can make rich dict"""
        # mixture of discrete-time, continuous-time models
        names = ["BH", "DT", "CNFGTR", "GN", "WG01"]
        tree = make_tree(tip_names=_aln.names)
        for name in names:
            sm = get_model(name)
            lf = sm.make_likelihood_function(tree)
            lf.set_alignment(_aln)
            _ = lf.to_rich_dict()

        # tests multiple alignments
        half = len(self.data) // 2
        aln1 = self.data[:half]
        aln2 = self.data[half:]
        loci_names = ["1st-half", "2nd-half"]
        loci = [aln1, aln2]
        tree = make_tree(tip_names=self.data.names)
        model = get_model("HKY85")
        lf = model.make_likelihood_function(tree, loci=loci_names)
        lf.set_alignment(loci)
        for i, loci_name in enumerate(loci_names):
            d = lf.to_rich_dict()
            alignment = d["alignment"]
            motif_probs = d["motif_probs"]
            self.assertEqual(alignment[loci_name], loci[i].to_rich_dict())
            self.assertEqual(motif_probs[loci_name], loci[i].get_motif_probs())
        # tests single alignment
        lf = model.make_likelihood_function(tree)
        lf.set_alignment(aln1)
        d = lf.to_rich_dict()
        alignment = d["alignment"]
        motif_probs = d["motif_probs"]
        self.assertEqual(alignment, aln1.to_rich_dict())
        self.assertEqual(motif_probs, aln1.get_motif_probs())

    def test_repr(self):
        """repr should not fail"""
        lf = self._makeLikelihoodFunction()
        got = repr(lf)
        self.assertIn("log-likelihood", got)

    def test_repr_html(self):
        "exercising for jupyter"
        lf = self._makeLikelihoodFunction()
        got = lf._repr_html_()
        self.assertIn("<p>log-likelihood", got)

    def test_get_set_name_properties(self):
        """correctly creates lf name attr"""
        lf = get_model("HKY85").make_likelihood_function(self.tree)
        self.assertEqual(lf.name, lf.model.name)
        lf.name = ""
        self.assertEqual(lf.name, "")


class ComparisonTests(TestCase):
    """comparisons of likelihood calcs with earlier pycogent"""

    def test_simple_codon(self):
        """return same likelihood for simple codon model"""
        # from docs/examples/neutral_test
        aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
        tree = load_tree("data/long_testseqs.tree")
        sm = get_model("MG94GTR")
        lf = sm.make_likelihood_function(tree, digits=3, space=2)
        lf.set_alignment(aln)
        mprobs = {
            "T": 0.23167456556082147,
            "C": 0.18775671406003158,
            "A": 0.36808846761453395,
            "G": 0.21248025276461296,
        }
        lengths = {
            "Human": 0.0930121148197949,
            "HowlerMon": 0.12455050902011401,
            "edge.0": 0.11566361642563996,
            "Mouse": 0.8352888057852214,
            "edge.1": 0.05801595370263309,
            "NineBande": 0.28274573844117873,
            "DogFaced": 0.33986809148384595,
        }
        rates = {
            "A/C": 1.0193255854344692,
            "A/G": 3.360125224439532,
            "A/T": 0.7324800239384959,
            "C/G": 0.9460411801916156,
            "C/T": 3.7077261494484466,
            "omega": 0.8991178277398568,
        }
        with lf.updates_postponed():
            lf.set_motif_probs(mprobs)
            for e, v in lengths.items():
                lf.set_param_rule("length", edge=e, init=v)
            for p, v in rates.items():
                lf.set_param_rule(p, init=v)
        assert_allclose(lf.lnL, -8636.180078337804)
        assert_allclose(lf.nfp, 13)

        # time-het variant
        mprobs = {
            "T": 0.23167456556082147,
            "C": 0.18775671406003158,
            "A": 0.36808846761453395,
            "G": 0.21248025276461296,
        }
        rates = {
            "A/C": 1.025238355386263,
            "A/G": 3.3839111255685093,
            "A/T": 0.7349455378588486,
            "C/G": 0.954931780018286,
            "C/T": 3.7190819957014605,
        }
        omega = {
            "Human": 0.5932371497829559,
            "HowlerMon": 0.9594738352514026,
            "edge.0": 1.1318706411033768,
            "Mouse": 0.9179170934077907,
            "edge.1": 0.39245724788647784,
            "NineBande": 1.2840869353261197,
            "DogFaced": 0.8432882455927212,
        }
        lengths = {
            "Human": 0.09424973246208838,
            "HowlerMon": 0.12400390171673893,
            "edge.0": 0.11480058221117359,
            "Mouse": 0.8347998998652563,
            "edge.1": 0.060388849069219555,
            "NineBande": 0.27863806025656374,
            "DogFaced": 0.3411611571500524,
        }
        with lf.updates_postponed():
            lf.set_motif_probs(mprobs)
            for e in lengths:
                lf.set_param_rule("length", edge=e, init=lengths[e])
                lf.set_param_rule("omega", edge=e, init=omega[e])
            for p, v in rates.items():
                lf.set_param_rule(p, init=v)
        assert_allclose(lf.lnL, -8632.135459058676)
        assert_allclose(lf.nfp, 19)

    def test_codon_rate_het(self):
        """recap rate het likelihoods"""
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype=DNA)
        tree = load_tree("data/primate_brca1.tree")
        cnf = get_model("CNFGTR")
        rate_lf = cnf.make_likelihood_function(
            tree, bins=["neutral", "adaptive"], digits=2, space=3
        )
        rate_lf.set_alignment(aln)

        rules = [
            {
                "par_name": "mprobs",
                "edges": None,
                "value": {
                    "TTT": 0.018732866280840695,
                    "TTC": 0.007767286018885166,
                    "TTA": 0.02116966189460859,
                    "TTG": 0.010813280536095034,
                    "TCT": 0.02512945476698142,
                    "TCC": 0.008224185196466647,
                    "TCA": 0.02208346024977155,
                    "TCG": 0.0015229972586049345,
                    "TAT": 0.010051781906792567,
                    "TAC": 0.002284495887907402,
                    "TGT": 0.020103563813585135,
                    "TGC": 0.0018275967103259215,
                    "TGG": 0.0039597928723728295,
                    "CTT": 0.010508681084374048,
                    "CTC": 0.007767286018885166,
                    "CTA": 0.01370697532744441,
                    "CTG": 0.012488577520560463,
                    "CCT": 0.026347852573865366,
                    "CCC": 0.006244288760280232,
                    "CCA": 0.019494364910143162,
                    "CCG": 0.0006091989034419738,
                    "CAT": 0.02208346024977155,
                    "CAC": 0.005178190679256778,
                    "CAA": 0.019646664636003654,
                    "CAG": 0.02375875723423698,
                    "CGT": 0.0031982942430703624,
                    "CGC": 0.0009137983551629607,
                    "CGA": 0.0010660980810234541,
                    "CGG": 0.002284495887907402,
                    "ATT": 0.019189765458422176,
                    "ATC": 0.007005787389582699,
                    "ATA": 0.0185805665549802,
                    "ATG": 0.01279317697228145,
                    "ACT": 0.028936947913493757,
                    "ACC": 0.004568991775814804,
                    "ACA": 0.022844958879074017,
                    "ACG": 0.0007614986293024672,
                    "AAT": 0.05558939993908011,
                    "AAC": 0.023454157782515993,
                    "AAA": 0.05558939993908011,
                    "AAG": 0.03441973804447152,
                    "AGT": 0.03807493146512336,
                    "AGC": 0.028632348461772768,
                    "AGA": 0.023149558330795003,
                    "AGG": 0.014011574779165398,
                    "GTT": 0.021321961620469083,
                    "GTC": 0.007005787389582699,
                    "GTA": 0.014773073408467865,
                    "GTG": 0.006853487663722205,
                    "GCT": 0.01370697532744441,
                    "GCC": 0.009594882729211088,
                    "GCA": 0.015839171489491318,
                    "GCG": 0.0013706975327444412,
                    "GAT": 0.031526043253122145,
                    "GAC": 0.010508681084374048,
                    "GAA": 0.07554066402680475,
                    "GAG": 0.030307645446238197,
                    "GGT": 0.01325007614986293,
                    "GGC": 0.008985683825769114,
                    "GGA": 0.016143770941212304,
                    "GGG": 0.006701187937861712,
                },
                "is_constant": True,
            },
            {
                "par_name": "omega",
                "edges": [
                    "Chimpanzee",
                    "Galago",
                    "Gorilla",
                    "HowlerMon",
                    "Human",
                    "Orangutan",
                    "Rhesus",
                    "edge.0",
                    "edge.1",
                    "edge.2",
                    "edge.3",
                ],
                "bins": ["adaptive"],
                "init": 1.171804535696805,
                "lower": 1.000001,
                "upper": 100,
            },
            {
                "par_name": "omega",
                "edges": [
                    "Chimpanzee",
                    "Galago",
                    "Gorilla",
                    "HowlerMon",
                    "Human",
                    "Orangutan",
                    "Rhesus",
                    "edge.0",
                    "edge.1",
                    "edge.2",
                    "edge.3",
                ],
                "bins": ["neutral"],
                "init": 0.011244516074953196,
                "lower": 1e-06,
                "upper": 1,
            },
            {
                "par_name": "C/T",
                "edges": None,
                "init": 4.197133081913153,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "C/G",
                "edges": None,
                "init": 1.959573221565973,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "A/T",
                "edges": None,
                "init": 0.7811160220383347,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "A/G",
                "edges": None,
                "init": 3.95572114917098,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "A/C",
                "edges": None,
                "init": 1.0689060203129388,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "length",
                "edges": ["Chimpanzee"],
                "init": 0.008587879544220245,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Galago"],
                "init": 0.5633709139717007,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Gorilla"],
                "init": 0.007511821374501532,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["HowlerMon"],
                "init": 0.141460056857809,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Human"],
                "init": 0.01830472793116926,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Orangutan"],
                "init": 0.023530857462230628,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Rhesus"],
                "init": 0.06742406161034178,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.0"],
                "init": 1.4909085767703391e-12,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.1"],
                "init": 0.010125599780668576,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.2"],
                "init": 0.03492051749140553,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.3"],
                "init": 0.02185436682925549,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "bprobs",
                "init": {"neutral": 0.1362876637533325, "adaptive": 0.8637123362466674},
                "lower": None,
                "upper": None,
                "edges": None,
            },
        ]
        rate_lf.apply_param_rules(rules)
        assert_allclose(rate_lf.lnL, -6755.451985437475)
        assert_allclose(rate_lf.nfp, 19)

        # check posterior bin probs match
        adaptive = [0.8138636520270726, 0.8723917957174725, 0.9922405018465282]
        got = rate_lf.get_bin_probs()
        assert_allclose(got["adaptive"][: len(adaptive)], adaptive, rtol=5e-6)
        self.assertEqual(got.shape, (2, len(aln) // 3))  # /3 because codon model

    def test_time_rate_het(self):
        """recap the zhang model"""
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype=DNA)
        tree = load_tree("data/primate_brca1.tree")
        cnf = get_model("CNFGTR")
        lf = cnf.make_likelihood_function(
            tree, bins=["0", "1", "2a", "2b"], digits=2, space=3
        )
        lf.set_alignment(aln)
        epsilon = 1e-6
        lf.set_param_rule("omega", bins=["0", "2a"], upper=1.0, init=1 - epsilon)
        lf.set_param_rule("omega", bins=["1", "2b"], is_constant=True, value=1.0)
        lf.set_param_rule(
            "omega",
            bins=["2a", "2b"],
            edges=["Chimpanzee", "Human"],
            init=99,
            lower=1.0,
            upper=100.0,
            is_constant=False,
        )

        rules = [
            {
                "par_name": "mprobs",
                "edges": None,
                "value": {
                    "TTT": 0.018732866280840695,
                    "TTC": 0.007767286018885166,
                    "TTA": 0.02116966189460859,
                    "TTG": 0.010813280536095034,
                    "TCT": 0.02512945476698142,
                    "TCC": 0.008224185196466647,
                    "TCA": 0.02208346024977155,
                    "TCG": 0.0015229972586049345,
                    "TAT": 0.010051781906792567,
                    "TAC": 0.002284495887907402,
                    "TGT": 0.020103563813585135,
                    "TGC": 0.0018275967103259215,
                    "TGG": 0.0039597928723728295,
                    "CTT": 0.010508681084374048,
                    "CTC": 0.007767286018885166,
                    "CTA": 0.01370697532744441,
                    "CTG": 0.012488577520560463,
                    "CCT": 0.026347852573865366,
                    "CCC": 0.006244288760280232,
                    "CCA": 0.019494364910143162,
                    "CCG": 0.0006091989034419738,
                    "CAT": 0.02208346024977155,
                    "CAC": 0.005178190679256778,
                    "CAA": 0.019646664636003654,
                    "CAG": 0.02375875723423698,
                    "CGT": 0.0031982942430703624,
                    "CGC": 0.0009137983551629607,
                    "CGA": 0.0010660980810234541,
                    "CGG": 0.002284495887907402,
                    "ATT": 0.019189765458422176,
                    "ATC": 0.007005787389582699,
                    "ATA": 0.0185805665549802,
                    "ATG": 0.01279317697228145,
                    "ACT": 0.028936947913493757,
                    "ACC": 0.004568991775814804,
                    "ACA": 0.022844958879074017,
                    "ACG": 0.0007614986293024672,
                    "AAT": 0.05558939993908011,
                    "AAC": 0.023454157782515993,
                    "AAA": 0.05558939993908011,
                    "AAG": 0.03441973804447152,
                    "AGT": 0.03807493146512336,
                    "AGC": 0.028632348461772768,
                    "AGA": 0.023149558330795003,
                    "AGG": 0.014011574779165398,
                    "GTT": 0.021321961620469083,
                    "GTC": 0.007005787389582699,
                    "GTA": 0.014773073408467865,
                    "GTG": 0.006853487663722205,
                    "GCT": 0.01370697532744441,
                    "GCC": 0.009594882729211088,
                    "GCA": 0.015839171489491318,
                    "GCG": 0.0013706975327444412,
                    "GAT": 0.031526043253122145,
                    "GAC": 0.010508681084374048,
                    "GAA": 0.07554066402680475,
                    "GAG": 0.030307645446238197,
                    "GGT": 0.01325007614986293,
                    "GGC": 0.008985683825769114,
                    "GGA": 0.016143770941212304,
                    "GGG": 0.006701187937861712,
                },
                "is_constant": True,
            },
            {
                "par_name": "omega",
                "edges": [
                    "Chimpanzee",
                    "Galago",
                    "Gorilla",
                    "HowlerMon",
                    "Human",
                    "Orangutan",
                    "Rhesus",
                    "edge.0",
                    "edge.1",
                    "edge.2",
                    "edge.3",
                ],
                "bins": ["0", "2a"],
                "init": 1.00000034130505e-06,
                "lower": 1e-06,
                "upper": 1.0,
            },
            {
                "par_name": "omega",
                "edges": [
                    "Chimpanzee",
                    "Galago",
                    "Gorilla",
                    "HowlerMon",
                    "Human",
                    "Orangutan",
                    "Rhesus",
                    "edge.0",
                    "edge.1",
                    "edge.2",
                    "edge.3",
                ],
                "bins": ["1", "2b"],
                "value": 1.0,
                "is_constant": True,
            },
            {
                "par_name": "omega",
                "edges": ["Chimpanzee", "Human"],
                "bins": ["2a", "2b"],
                "init": 19.00550828862512,
                "lower": 1.0,
                "upper": 100.0,
            },
            {
                "par_name": "C/T",
                "edges": None,
                "init": 4.107017179700691,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "C/G",
                "edges": None,
                "init": 1.96190792766859,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "A/T",
                "edges": None,
                "init": 0.7729320058102641,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "A/G",
                "edges": None,
                "init": 3.9100044579566493,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "A/C",
                "edges": None,
                "init": 1.0654568108482438,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "length",
                "edges": ["Chimpanzee"],
                "init": 0.008556133589133694,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Galago"],
                "init": 0.5601765224662885,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Gorilla"],
                "init": 0.007512028921106184,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["HowlerMon"],
                "init": 0.1409451751006042,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Human"],
                "init": 0.018235814843667496,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Orangutan"],
                "init": 0.023516404248762182,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["Rhesus"],
                "init": 0.06730954128649547,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.0"],
                "init": 1.4909085767703391e-12,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.1"],
                "init": 0.01012135420651205,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.2"],
                "init": 0.0349479495846757,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "length",
                "edges": ["edge.3"],
                "init": 0.02202892799465516,
                "lower": 0.0,
                "upper": 10.0,
            },
            {
                "par_name": "bprobs",
                "init": {
                    "0": 0.034996999380649305,
                    "1": 0.3355065884659966,
                    "2a": 0.0791090296297802,
                    "2b": 0.5503873825235739,
                },
                "lower": None,
                "upper": None,
                "edges": None,
            },
        ]
        lf.apply_param_rules(rules)
        # values from earlier pycogent
        lnL, nfp = -6753.45662012937, 21
        assert_allclose(lf.lnL, lnL, rtol=1e-5)
        assert_allclose(lf.nfp, nfp)

    def test_loci(self):
        """recap multiple-loci"""
        from cogent3.recalculation.scope import ALL, EACH

        aln = load_aligned_seqs("data/long_testseqs.fasta")
        half = len(aln) // 2
        aln1 = aln[:half]
        aln2 = aln[half:]
        loci_names = ["1st-half", "2nd-half"]
        loci = [aln1, aln2]
        tree = make_tree(tip_names=aln.names)
        sm = get_model("HKY85")
        lf = sm.make_likelihood_function(tree, loci=loci_names, digits=2, space=3)
        lf.set_param_rule("length", is_independent=False)
        lf.set_param_rule("kappa", loci=ALL)
        lf.set_alignment(loci)
        # lf.optimise()
        # rules = lf.get_param_rules()
        rules = [
            {
                "par_name": "kappa",
                "init": 3.9795417148915124,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "mprobs",
                "locus": "1st-half",
                "edges": [
                    "DogFaced",
                    "HowlerMon",
                    "Human",
                    "Mouse",
                    "NineBande",
                    "root",
                ],
                "value": {
                    "T": 0.22338072669826226,
                    "C": 0.1843601895734597,
                    "A": 0.3824644549763033,
                    "G": 0.20979462875197472,
                },
                "is_constant": True,
            },
            {
                "par_name": "mprobs",
                "locus": "2nd-half",
                "edges": [
                    "DogFaced",
                    "HowlerMon",
                    "Human",
                    "Mouse",
                    "NineBande",
                    "root",
                ],
                "value": {
                    "T": 0.23996840442338072,
                    "C": 0.1911532385466035,
                    "A": 0.35371248025276464,
                    "G": 0.21516587677725119,
                },
                "is_constant": True,
            },
            {
                "par_name": "length",
                "init": 0.12839234555151016,
                "lower": 0.0,
                "upper": 10.0,
                "is_independent": False,
            },
        ]
        lf.apply_param_rules(rules)
        lnL, nfp = -9168.333116203343, 2
        assert_allclose(lf.lnL, lnL)
        assert_allclose(lf.nfp, nfp)

        lf.set_param_rule("kappa", loci=EACH)
        lf.optimise(local=True, show_progress=False)
        rules = [
            {
                "par_name": "kappa",
                "edges": ["DogFaced", "HowlerMon", "Human", "Mouse", "NineBande"],
                "locus": "1st-half",
                "init": 4.3253540083718285,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "kappa",
                "edges": ["DogFaced", "HowlerMon", "Human", "Mouse", "NineBande"],
                "locus": "2nd-half",
                "init": 3.738155123420754,
                "lower": 1e-06,
                "upper": 1000000.0,
            },
            {
                "par_name": "mprobs",
                "locus": "1st-half",
                "edges": [
                    "DogFaced",
                    "HowlerMon",
                    "Human",
                    "Mouse",
                    "NineBande",
                    "root",
                ],
                "value": {
                    "T": 0.22338072669826226,
                    "C": 0.1843601895734597,
                    "A": 0.3824644549763033,
                    "G": 0.20979462875197472,
                },
                "is_constant": True,
            },
            {
                "par_name": "mprobs",
                "locus": "2nd-half",
                "edges": [
                    "DogFaced",
                    "HowlerMon",
                    "Human",
                    "Mouse",
                    "NineBande",
                    "root",
                ],
                "value": {
                    "T": 0.23996840442338072,
                    "C": 0.1911532385466035,
                    "A": 0.35371248025276464,
                    "G": 0.21516587677725119,
                },
                "is_constant": True,
            },
            {
                "par_name": "length",
                "init": 0.12844512742274064,
                "lower": 0.0,
                "upper": 10.0,
                "is_independent": False,
            },
        ]
        lf.apply_param_rules(rules)
        lnL, nfp = -9167.537271862948, 3
        assert_allclose(lf.lnL, lnL)
        assert_allclose(lf.nfp, nfp)


if __name__ == "__main__":
    main()
