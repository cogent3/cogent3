"""
This file defines a class for controlling the scope and heterogeneity of
parameters involved in a maximum-likelihood based tree analysis.
"""


import pickle
import warnings

import numpy

from cogent3.align import dp_calculation
from cogent3.align.pairwise import AlignableSeq
from cogent3.core.tree import TreeError
from cogent3.evolve import likelihood_calculation
from cogent3.evolve.likelihood_function import LikelihoodFunction as _LF
from cogent3.maths.stats.information_criteria import aic, bic
from cogent3.recalculation.scope import _indexed
from cogent3.util.misc import adjusted_gt_minprob


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Andrew Butterfield", "Peter Maxwell", "Gavin Huttley", "Helen Lindsay"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.ed.au"
__status__ = "Production"


def _category_names(dimension, specified):
    if type(specified) is int:
        cats = [f"{dimension}{i}" for i in range(specified)]
    else:
        cats = tuple(specified)
    assert len(cats) >= 1, cats
    assert len(set(cats)) == len(cats), f"{dimension} names must be unique"
    return list(cats)


def load(filename):
    # first cut at saving pc's
    f = open(filename, "rb")
    (version, info, pc) = pickle.load(f)
    assert version < 2.0, version
    pc.update_intermediate_values()
    return pc


class _LikelihoodParameterController(_LF):
    """A ParameterController works by setting parameter rules. For each
    parameter in the model the edges of the tree are be partitioned into groups
    that share one value.

    For usage see the set_param_rule method.
    """

    # Basically wrapper around the more generic recalculation.ParameterController
    # class, which doesn't know about trees.

    def __init__(
        self,
        model,
        tree,
        bins=1,
        loci=1,
        optimise_motif_probs=False,
        motif_probs_from_align=False,
        **kw,
    ):
        # cache of arguments used to construct
        d = locals()
        exclude = ("self", "__class__", "kw")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}
        self._serialisable.update(kw)

        self.model = self._model = model
        self.tree = self._tree = tree
        self.seq_names = tree.get_tip_names()
        self.locus_names = _category_names("locus", loci)
        self.bin_names = _category_names("bin", bins)
        self.posn_names = [str(i) for i in range(model.word_length)]
        self.motifs = self._motifs = model.get_motifs()
        self._mprob_motifs = list(model.get_mprob_alphabet())
        defn = self.make_likelihood_defn(**kw)
        super(_LikelihoodParameterController, self).__init__(defn)
        self.set_default_param_rules()
        self.set_default_tree_parameter_rules()
        self.mprobs_from_alignment = motif_probs_from_align
        self.optimise_motif_probs = optimise_motif_probs
        self._name = None
        self._format = {}

    def save(self, filename):
        with open(filename, "w") as f:
            temp = {}
            try:
                for d in self.defns:
                    temp[id(d)] = d.values
                    del d.values
                pickle.dump((1.0, None, self), f)
            finally:
                for d in self.defns:
                    if id(d) in temp:
                        d.values = temp[id(d)]

    def set_default_tree_parameter_rules(self):
        """Lengths are set to the values found in the tree (if any), and
        free to be optimised independently.
        Other parameters are scoped based on the unique values found in the
        tree (if any) or default to having one value shared across the whole
        tree"""
        with self.updates_postponed():
            edges = self.tree.get_edge_vector()
            for par_name in self.model.get_param_list():
                try:
                    values = dict(
                        [
                            (edge.name, edge.params[par_name])
                            for edge in edges
                            if not edge.isroot()
                        ]
                    )
                    (uniq, index) = _indexed(values)
                except KeyError:
                    continue  # new parameter
                for (u, value) in enumerate(uniq):
                    group = [edge for (edge, i) in list(index.items()) if i == u]
                    self.set_param_rule(par_name, edges=group, init=value)
            for edge in edges:
                if edge.length is not None:
                    try:
                        self.set_param_rule("length", edge=edge.name, init=edge.length)
                    except KeyError:
                        # hopefully due to being a discrete model
                        warnings.warn("Ignoring tree edge lengths", stacklevel=4)
                        break

    def set_motif_probs_from_data(
        self,
        align,
        locus=None,
        is_constant=None,
        include_ambiguity=False,
        is_independent=None,
        auto=False,
        pseudocount=None,
        **kwargs,
    ):
        counts = self.model.count_motifs(align, include_ambiguity=include_ambiguity)
        if is_constant is None:
            is_constant = not self.optimise_motif_probs

        pseudocount = 0.0 if is_constant or (counts != 0).all() else pseudocount or 0.5
        counts += pseudocount
        mprobs = counts / (1.0 * sum(counts))
        self.set_motif_probs(
            mprobs,
            locus=locus,
            is_constant=is_constant,
            is_independent=is_independent,
            auto=auto,
            **kwargs,
        )

    def set_motif_probs(
        self,
        motif_probs,
        locus=None,
        bin=None,
        is_constant=None,
        is_independent=None,
        auto=False,
        **kwargs,
    ):

        motif_probs = self.model.adapt_motif_probs(motif_probs, auto=auto)
        motif_probs = adjusted_gt_minprob(motif_probs, minprob=1e-6)
        if is_constant is None:
            is_constant = not self.optimise_motif_probs
        self.model.set_param_controller_motif_probs(
            self,
            motif_probs,
            is_constant=is_constant,
            bin=bin,
            locus=locus,
            is_independent=is_independent,
            **kwargs,
        )
        if not auto:
            self.mprobs_from_alignment = False  # should be done per-locus

    def set_expm(self, expm):
        assert expm in ["pade", "either", "eigen", "checked"], expm
        self.set_param_rule("expm", is_constant=True, value=expm)

    def make_calculator(self, **kw):
        return super(_LikelihoodParameterController, self).make_calculator(**kw)

    def _process_scope_info(
        self,
        edge=None,
        tip_names=None,
        edges=None,
        clade=None,
        stem=None,
        outgroup_name=None,
    ):
        """From information specifying the scope of a parameter derive a list of
        edge names"""

        if edges is not None:
            if tip_names or edge:
                raise TreeError("Only ONE of edge, edges or tip_names")
        elif edge is not None:
            if tip_names:
                raise TreeError("Only ONE of edge, edges or tip_names")
            edges = [edge]
        elif tip_names is None:
            edges = None  # meaning all edges
        elif len(tip_names) != 2:
            raise TreeError("tip_names must contain 2 species")
        else:
            (species1, species2) = tip_names
            if stem is None:
                stem = False
            if clade is None:
                clade = not stem
            edges = self.tree.get_edge_names(
                species1, species2, stem=stem, clade=clade, outgroup_name=outgroup_name
            )

        return edges

    def apply_param_rules(self, rules):
        """batch applies a collection of param rules"""
        with self.updates_postponed():
            for rule in rules:
                self.set_param_rule(**rule)

    def set_time_heterogeneity(
        self,
        exclude_params=None,
        edge_sets=None,
        is_independent=None,
        is_constant=False,
        value=None,
        lower=None,
        init=None,
        upper=None,
    ):
        """modifies the scope of substitution model rate params

        Parameters
        ----------
        exclude_params
            name(s) of substitution model predicate(s) to be excluded
        edge_sets
            series of dicts with an 'edges' key. Can also specify
            is_independent, is_constant etc.. If those are not provided, the
            method argument values are applied
        is_independent : bool
            whether edges in all edge sets are to be considered independent.
            default is False
            Overridden by edge_sets values.
        is_constant : bool
            makes constant all rate term parameters for all edge sets.
            Overridden by edge_sets values.
        value
            value for constant parameters, only valid when is_constant.
            Overridden by edge_sets values.
        lower, init, upper
            lower bound, starting value, upper bound for all parameters for
            all edge sets. Only valid if not is_constant.
            Overridden by edge_sets values.
        """
        if is_constant and any([lower, init, upper]):
            raise ValueError("cannot specify bounds or init for a constant param")

        if is_constant:
            kwargs = dict(is_constant=True, value=value)
        else:
            kwargs = dict(
                is_independent=is_independent, init=init, lower=lower, upper=upper
            )

        rate_terms = self._model.get_param_list()
        exclude_params = exclude_params or []
        if exclude_params and type(exclude_params) == str:
            exclude_params = [exclude_params]

        for param in exclude_params:
            if param not in rate_terms:
                raise ValueError(f"'{param}' not a valid rate param")

            rate_terms.remove(param)

        if edge_sets is None:
            # this just makes the following algorithm consistent
            # but we need to exclude any edges assigned to discrete psubs
            dpsubs = self.defn_for.get("dpsubs", None)
            exclude_edges = set()
            if dpsubs:
                dims = dpsubs.valid_dimensions
                index = dims.index("edge")
                for k in dpsubs.assignments:
                    exclude_edges.add(k[index])

            edge_sets = [
                dict(edges=[n])
                for n in self.tree.get_node_names(includeself=False)
                if n not in exclude_edges
            ]
        elif type(edge_sets) == dict:
            edge_sets = [edge_sets]

        # we make param rules
        param_rules = []
        for edge_set in edge_sets:
            edges = edge_set.get("edges", None)
            if type(edges) == str:
                edges = [edges]

            if edges:
                edges = list(edges)
                edge_set["edges"] = edges

            rule_base = kwargs.copy()
            rule_base.update(edge_set)
            for param in rate_terms:
                rule = rule_base.copy()
                rule.update(dict(par_name=param))
                param_rules.append(rule)

        self.apply_param_rules(param_rules)

    def set_param_rule(
        self,
        par_name,
        is_independent=None,
        is_constant=False,
        value=None,
        lower=None,
        init=None,
        upper=None,
        **scope_info,
    ):
        """Define a model constraint for par_name. Parameters can be set
        constant or split according to tree/bin scopes.

        Parameters
        ----------
        par_name
            The model parameter being modified.
        is_constant, value
            if True, the parameter is held constant at
            value, if provided, or the likelihood functions current value.
        is_independent
            whether the partition specified by scope/bin
            arguments are to be considered independent.
        lower, init, upper
            specify the lower bound, initial value and
            upper bound for optimisation. Can be set separately.
        bin, bins
            the name(s) of the bin to apply rule.
        locus, loci
            the name of the locus/loci to apply rule.
        **scope_info
            tree scope arguments

              - edge, edges: The name of the tree edge(s) affected by rule.
              - tip_names: a tuple of two tip names, specifying a tree scope
                to apply rule.
              - outgroup_name: A tip name that, provided along with tip_names,
                ensures a consistently specified tree scope.
              - clade: The rule applies to all edges descending from the most
                recent common ancestor defined by the tip_names+outgroup_name
                arguments.
              - stem: The rule applies to the edge preceding the most recent
                common ancestor defined by the tip_names+outgroup_name
                arguments.
        """

        par_name = str(par_name)

        scopes = {}
        for (single, plural) in [
            ("bin", "bins"),
            ("locus", "loci"),
            ("position", "positions"),
            ("motif", "motifs"),
        ]:
            if single in scope_info:
                v = scope_info.pop(single)
                if v:
                    assert isinstance(v, str), f"{plural}=, maybe?"
                    assert plural not in scope_info
                    scopes[single] = [v]
            elif plural in scope_info:
                v = scope_info.pop(plural)
                if v:
                    scopes[single] = v

        edges = self._process_scope_info(**scope_info)
        if edges:
            scopes["edge"] = edges

        if is_constant:
            assert not (init or lower or upper)
        elif init is not None:
            assert not value
            value = init
        self.assign_all(
            par_name, scopes, value, lower, upper, is_constant, is_independent
        )

    def set_local_clock(self, tip1name, tip2name):
        """Constrain branch lengths for tip1name and tip2name to be equal.
        This is a molecular clock condition. Currently only valid for tips
        connected to the same node.

        Note: This is just a convenient interface to setParameterRule.
        """
        self.set_param_rule(
            "length", tip_names=[tip1name, tip2name], clade=True, is_independent=0
        )

    def set_constant_lengths(self, tree=None, exclude_list=None):
        """Constrains edge lengths to those in the tree.

        Parameters
        ----------
        tree
            must have the same topology as the current model.
            If not provided, the current tree length's are used.
        exclude_list
            a list of edge names whose branch lengths
            will be constrained.

        """

        exclude_list = exclude_list or []

        if tree is None:
            tree = self.tree

        with self.updates_postponed():
            for edge in tree.get_edge_vector():
                if edge.length is None or edge.name in exclude_list:
                    continue
                self.set_param_rule(
                    "length", edge=edge.name, is_constant=1, value=edge.length
                )

    def get_aic(self, second_order=False):
        """returns Aikake Information Criteria

        Parameters
        ----------
        second_order
            if true, the second
            adjusted by the alignment length

        """
        if second_order:
            sequence_length = sum(
                len(self.get_param_value("lht", locus=l).index)
                for l in self.locus_names
            )
        else:
            sequence_length = None

        lnL = self.get_log_likelihood()
        nfp = self.get_num_free_params()
        return aic(lnL, nfp, sequence_length)

    def get_bic(self):
        """returns the Bayesian Information Criteria"""
        sequence_length = sum(
            len(self.get_param_value("lht", locus=l).index) for l in self.locus_names
        )
        lnL = self.get_log_likelihood()
        nfp = self.get_num_free_params()
        return bic(lnL, nfp, sequence_length)


class AlignmentLikelihoodFunction(_LikelihoodParameterController):
    def set_default_param_rules(self):
        try:
            self.assign_all("fixed_motif", None, value=-1, const=True, independent=True)
        except KeyError:
            pass

    def make_likelihood_defn(self, sites_independent=True, discrete_edges=None):
        defns = self.model.make_param_controller_defns(bin_names=self.bin_names)
        if discrete_edges is not None:
            from .discrete_markov import PartialyDiscretePsubsDefn

            defns["psubs"] = PartialyDiscretePsubsDefn(
                self.motifs, defns["psubs"], discrete_edges
            )
        return likelihood_calculation.make_total_loglikelihood_defn(
            self.tree,
            defns["align"],
            defns["psubs"],
            defns["word_probs"],
            defns["bprobs"],
            self.bin_names,
            self.locus_names,
            sites_independent,
        )

    def set_alignment(self, aligns, motif_pseudocount=None):
        """set the alignment to be used for computing the likelihood."""
        if type(aligns) is not list:
            aligns = [aligns]
        assert len(aligns) == len(self.locus_names), len(aligns)
        tip_names = set(self.tree.get_tip_names())
        for index, aln in enumerate(aligns):
            if len(aligns) > 1:
                locus_name = f"for locus '{self.locus_names[index]}'"
            else:
                locus_name = ""
            assert not set(aln.names).symmetric_difference(
                tip_names
            ), "Tree tip names %s and aln seq names %s don't match %s" % (
                self.tree.get_tip_names(),
                aln.names,
                locus_name,
            )
            assert "root" not in aln.names, "'root' is a reserved name."
        with self.updates_postponed():
            for (locus_name, align) in zip(self.locus_names, aligns):
                self.assign_all(
                    "alignment", {"locus": [locus_name]}, value=align, const=True
                )
                if self.mprobs_from_alignment:
                    self.set_motif_probs_from_data(
                        align,
                        locus=locus_name,
                        auto=True,
                        pseudocount=motif_pseudocount,
                    )


class SequenceLikelihoodFunction(_LikelihoodParameterController):
    def set_default_param_rules(self):
        pass

    def make_likelihood_defn(
        self, sites_independent=None, with_indel_params=True, kn=True
    ):
        assert sites_independent is None or not sites_independent
        assert len(self.locus_names) == 1
        return dp_calculation.make_forward_tree_defn(
            self.model,
            self.tree,
            self.bin_names,
            with_indel_params=with_indel_params,
            kn=kn,
        )

    def set_sequences(self, seqs, locus=None):
        from cogent3.core.alignment import SequenceCollection

        leaves = {}

        if isinstance(seqs, SequenceCollection):
            seqs = seqs.named_seqs

        for (name, seq) in list(seqs.items()):
            # if has uniq, probably already a likelihood tree leaf obj already
            if hasattr(seq, "uniq"):
                # XXX more checks - same alphabet as model, name etc ...
                leaf = seq
            else:
                leaf = self.model.convert_sequence(seq, name)
            leaf = AlignableSeq(leaf)
            leaves[name] = leaf
            assert name != "root", "'root' is a reserved name."
        self.set_pogs(leaves, locus=locus)

    def set_pogs(self, leaves, locus=None):
        with self.updates_postponed():
            for (name, pog) in list(leaves.items()):
                self.set_param_rule("leaf", edge=name, value=pog, is_constant=True)
            if self.mprobs_from_alignment:
                counts = numpy.sum(
                    [pog.leaf.get_motif_counts() for pog in list(leaves.values())], 0
                )
                mprobs = counts / (1.0 * sum(counts))
                self.set_motif_probs(mprobs, locus=locus, is_constant=True, auto=True)
