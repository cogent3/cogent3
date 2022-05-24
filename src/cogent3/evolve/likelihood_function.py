import json
import random

from collections import defaultdict
from copy import deepcopy

import numpy

from cogent3.core.alignment import ArrayAlignment
from cogent3.evolve import substitution_model
from cogent3.evolve.simulate import AlignmentEvolver, random_sequence
from cogent3.maths.matrix_exponential_integration import expected_number_subs
from cogent3.maths.matrix_logarithm import is_generator_unique
from cogent3.maths.measure import (
    paralinear_continuous_time,
    paralinear_discrete_time,
)
from cogent3.recalculation.definition import ParameterController
from cogent3.util import table
from cogent3.util.dict_array import DictArrayTemplate
from cogent3.util.misc import adjusted_gt_minprob, get_object_provenance


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Gavin Huttley",
    "Andrew Butterfield",
    "Peter Maxwell",
    "Matthew Wakefield",
    "Rob Knight",
    "Brett Easton",
    "Ben Kaehler",
    "Ananias Iliadis",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


# cogent3.evolve.parameter_controller.LikelihoodParameterController tells the
# recalculation framework to use this subclass rather than the generic
# recalculation Calculator.  It adds methods which are useful for examining
# the parameter, psub, mprob and likelihood values after the optimisation is
# complete.


def _get_keyed_rule_indices(rules):
    """returns {frozesent((par_name, edge1, edge2, ..)): index}"""
    new = {}
    for i, rule in enumerate(rules):
        edges = rule.get("edges", rule.get("edge", None)) or []
        edges = [edges] if type(edges) == str else edges
        par_name = rule["par_name"]
        key = frozenset([par_name] + edges)
        new[key] = i
    return new


def update_rule_value(rich, null):
    """applies value from null rule to rich rule"""
    val_key = "init" if "init" in rich else "value"
    rich[val_key] = null.get("init", null.get("value"))
    return rich


def extend_rule_value(rich, nulls):
    """creates new rich rules from edges in null rules"""
    val_key = "init" if "init" in rich else "value"
    rules = []
    for null in nulls:
        edges = null.get("edges", null.get("edge"))
        edges = [edges] if type(edges) == str else edges
        for edge in edges:
            rule = deepcopy(rich)
            rule["edge"] = edge
            rule[val_key] = null.get("init", null.get("value"))
            rules.append(rule)
    return rules


def update_scoped_rules(rich, null):
    """returns rich rules with values derived from those in null rules"""
    new_rules = []
    rich = deepcopy(rich)
    # we build a dict keyed by frozen set consisting of the param name
    # and affected edges. The dict value is the list index in original.
    richd = _get_keyed_rule_indices(rich)
    nulld = _get_keyed_rule_indices(null)
    common = set(richd) & set(nulld)
    # 1-to-1 mapping, just extract the param value
    for key in common:
        rule = update_rule_value(rich[richd[key]], null[nulld[key]])
        new_rules.append(rule)

    # following rules differing in scope
    rich_remainder = set(richd) - set(nulld)
    null_remainder = set(nulld) - set(richd)
    for rich_key in rich_remainder:
        matches = []
        rich_rule = rich[richd[rich_key]]
        pname = rich_rule["par_name"]
        enames = rich_rule.get("edges", rich_rule.get("edge", None))
        if type(enames) == str:
            enames = [enames]
        enames = None if enames is None else set(enames)
        for null_key in null_remainder:
            null_rule = null[nulld[null_key]]
            if pname != null_rule["par_name"]:
                continue
            # parameter now fully general
            if enames is None:
                matches.append(null_rule)
                continue

            null_enames = null_rule.get("edges", null_rule.get("edge", None))
            null_enames = None if null_enames is None else set(null_enames)
            if None in (enames, null_enames) or null_enames & enames:
                matches.append(null_rule)

        if enames is None:  # rich rule is "free"
            new_rules.extend(extend_rule_value(rich_rule, matches))
            continue

        if len(matches) > 1 and enames is not None:
            raise ValueError(f"{rich_key} has too many mappings {matches}")

        match = matches[0]
        new_rules.append(update_rule_value(rich_rule, match))
    return new_rules


def _get_param_mapping(rich, simple):
    """returns {simple_param_name: rich_param_name, ...}, the mapping of simple
    to rich parameters based on matrix coordinates
    """
    assert len(rich) >= len(simple)
    simple_to_rich = defaultdict(set)
    rich_to_simple = defaultdict(set)
    for simple_param in simple:
        simple_coords = simple[simple_param]
        for rich_param in rich:
            rich_coords = rich[rich_param]
            if rich_coords <= simple_coords:
                simple_to_rich[simple_param].add(rich_param)
                rich_to_simple[rich_param].add(simple_param)

    for rich_param, simple_counterparts in rich_to_simple.items():
        if len(simple_counterparts) == 1:
            continue

        sized_simple = [(len(simple[param]), param) for param in simple_counterparts]
        sized_simple.sort()
        if sized_simple[0][0] == sized_simple[1][0]:
            msg = f"{sized_simple[0][1]} and {sized_simple[1][1]} tied for matrix space"
            raise ValueError(msg)

        _, chosen = sized_simple.pop(0)
        rich_to_simple[rich_param] = [chosen]
        for _, simple_param in sized_simple:
            simple_to_rich[simple_param].remove(rich_param)

    return simple_to_rich


class _ParamProjection:
    """projects parameter names, values between nested models"""

    def __init__(self, simple_model, rich_model, motif_probs, same=True):
        # construct following by calling the functions we wrote
        self._rich_coords = rich_model.get_param_matrix_coords(include_ref_cell=True)
        self._simple_coords = simple_model.get_param_matrix_coords(
            include_ref_cell=True
        )
        self._param_map = _get_param_mapping(self._rich_coords, self._simple_coords)
        self._same = same
        # end of constructing attributes
        self._motif_probs = motif_probs
        self._ref_val = self._set_ref_val(same)
        self.projected_rate = {False: self._rate_not_same}.get(same, self._rate_same)

    def _set_ref_val(self, same):
        """returns the motif prob corresponding to the model reference cell"""
        if same:
            return 1
        else:
            i, j = list(self._rich_coords["ref_cell"])[0]
            return self._motif_probs[j]

    def _rate_not_same(self, simple_param, mle):
        """returns {rich_param: val, ...} from simple_param: val"""
        ref_val = self._ref_val
        new_terms = {}
        for rich_param in self._param_map[simple_param]:
            if rich_param == "ref_cell":
                continue
            for i, j in self._rich_coords[rich_param]:
                new_terms[rich_param] = self._motif_probs[j] * mle / ref_val

        return new_terms

    def _rate_same(self, simple_param, mle):
        new_terms = {}
        for rich_param in self._param_map[simple_param]:
            if rich_param == "ref_cell":
                continue
            for i, j in self._rich_coords[rich_param]:
                new_terms[rich_param] = mle
        return new_terms

    def update_param_rules(self, rules):
        new_rules = []
        if not self._same:
            rules = rules[:] + [dict(par_name="ref_cell", init=1.0, edges=None)]

        for rule in rules:
            # get the param name, mle, call self.projected_rate
            name = rule["par_name"]
            if name in ("mprobs", "length"):
                new_rules.append(rule)
                continue

            if rule.get("is_constant", False):
                par_val_key = "value"
            else:
                par_val_key = "init"

            mle = rule[par_val_key]

            proj_rate = self.projected_rate(name, mle)
            for new_name, new_mle in proj_rate.items():
                rule_dict = rule.copy()
                rule_dict["par_name"] = new_name
                # update it with the new parname and mle and append to new rules
                rule_dict["init"] = new_mle
                new_rules.append(rule_dict)

        return new_rules


def compatible_likelihood_functions(lf1, lf2):
    """returns True if all attributes of the two likelihood functions are compatible
    for mapping parameters, else raises ValueError or an AssertionError"""
    # tree's must have the same topology AND be oriented the same way
    # plus have the same edge names
    if len(lf1.bin_names) != 1 or len(lf1.bin_names) != len(lf2.bin_names):
        raise NotImplementedError("Too many bins")
    if len(lf1.locus_names) != 1 or len(lf1.locus_names) != len(lf2.locus_names):
        raise NotImplementedError("Too many loci")
    if lf1.model.get_motifs() != lf2.model.get_motifs():
        raise AssertionError("Motifs don't match")
    if lf1.tree.get_newick(with_node_names=True) != lf2.tree.get_newick(
        with_node_names=True
    ):
        raise AssertionError("Topology, Orientation or node names don't match")
    return True


class LikelihoodFunction(ParameterController):
    @property
    def lnL(self):
        """log-likelihood"""
        return self.get_log_likelihood()

    def get_log_likelihood(self):
        return self.get_final_result()

    def get_all_psubs(self):
        """returns all psubs as a dict keyed by used dimensions"""
        try:
            defn = self.defn_for["dsubs"]
        except KeyError:
            defn = self.defn_for["psubs"]

        used_dims = defn.used_dimensions()
        vdims = defn.valid_dimensions
        indices = [vdims.index(k) for k in used_dims if k in vdims]
        result = {}
        darr_template = DictArrayTemplate(self._motifs, self._motifs)
        for scope, index in defn.index.items():
            psub = defn.values[index]
            key = tuple(numpy.take(scope, indices))
            result[key] = darr_template.wrap(psub)
        return result

    def get_psub_for_edge(self, name, **kw):
        """returns the substitution probability matrix for the named edge

        Parameters
        ----------
        name : str
            name of the edge

        Returns
        -------
        DictArray
        """
        # todo handle case of multiple loci
        try:
            # For PartialyDiscretePsubsDefn
            array = self.get_param_value("dpsubs", edge=name, **kw)
        except KeyError:
            array = self.get_param_value("psubs", edge=name, **kw)
        return DictArrayTemplate(self._motifs, self._motifs).wrap(array)

    def get_all_rate_matrices(self, calibrated=True):
        """returns all rate matrices (Q) as a dict, keyed by scope

        Parameters
        ----------
        calibrated : bool
            If True, the rate matrix is scaled such that
            ``sum(pi_i * Qii) == 1``. If False, the calibrated matrix is
            multiplied by the length parameter (and the rate parameter for a
            bin if it is a rate heterogeneity model).

        Returns
        -------
        {scope: DictArray, ...}

        Notes
        -----
        If a single rate matrix (e.g. it's a time-homogeneous model), the key
        is an empty tuple.
        """
        defn = self.defn_for["Q"]

        rate_het = self.defn_for.get("rate", False)
        if rate_het:
            bin_index = rate_het.valid_dimensions.index("bin")
            bin_names = [k[bin_index] for k in rate_het.index]
            bin_names = {n: i for i, n in enumerate(bin_names)}
            bin_index = defn.valid_dimensions.index("bin")
        else:
            bin_names = None
            bin_index = None

        used_dims = defn.used_dimensions()
        edge_index = defn.valid_dimensions.index("edge")

        indices = {defn.valid_dimensions.index(k) for k in used_dims}
        if not calibrated:
            indices.add(edge_index)

        if not calibrated and rate_het:
            indices.add(bin_index)

        indices = list(sorted(indices))
        result = {}
        darr_template = DictArrayTemplate(self._motifs, self._motifs)
        for scope, index in defn.index.items():
            q = defn.values[index]  # this gives the appropriate Q
            # from scope we extract only the relevant dimensions
            key = tuple(numpy.take(scope, indices))
            q = q.copy()
            if not calibrated:
                length = self.get_param_value("length", edge=scope[edge_index])
                if rate_het:
                    bdex = bin_names[scope[bin_index]]
                    rate = rate_het.values[bdex]
                    length *= rate
                q *= length
            result[key] = darr_template.wrap(q)
            if not indices and calibrated:
                break  # single rate matrix

        return result

    def get_rate_matrix_for_edge(self, name, calibrated=True, **kw):
        """returns the rate matrix (Q) for the named edge

        Parameters
        ----------
        name : str
            name of the edge
        calibrated : bool
            If True, the rate matrix is scaled such that
            ``sum(pi_i * Qii) == 1``. If False, the calibrated matrix is
            multiplied by the length parameter (and the rate parameter for a
            bin if it is a rate heterogeneity model).

        Notes
        -----
        If ``calibrated=False``, ``expm(Q)`` will give the same result as
        ``self.get_psub_for_edge(name)``
        """
        # todo handle case of multiple loci
        try:
            array = self.get_param_value("Q", edge=name, **kw)
            array = array.copy()
            if not calibrated:
                length = self.get_param_value("length", edge=name, **kw)
                array *= length
        except KeyError as err:
            if err[0] == "Q" and name != "Q":
                raise RuntimeError("rate matrix not known by this model")
            else:
                raise
        return DictArrayTemplate(self._motifs, self._motifs).wrap(array)

    def _getLikelihoodValuesSummedAcrossAnyBins(self, locus=None):
        if self.bin_names and len(self.bin_names) > 1:
            root_lhs = [
                self.get_param_value("lh", locus=locus, bin=bin)
                for bin in self.bin_names
            ]
            bprobs = self.get_param_value("bprobs")
            root_lh = bprobs.dot(root_lhs)
        else:
            root_lh = self.get_param_value("lh", locus=locus)
        return root_lh

    def get_full_length_likelihoods(self, locus=None):
        """Array of [site, motif] likelihoods from the root of the tree"""
        root_lh = self._getLikelihoodValuesSummedAcrossAnyBins(locus=locus)
        root_lht = self.get_param_value("root", locus=locus)
        return root_lht.get_full_length_likelihoods(root_lh)

    def get_G_statistic(self, return_table=False, locus=None):
        """Goodness-of-fit statistic derived from the unambiguous columns"""
        root_lh = self._getLikelihoodValuesSummedAcrossAnyBins(locus=locus)
        root_lht = self.get_param_value("root", locus=locus)
        return root_lht.calc_G_statistic(root_lh, return_table)

    def reconstruct_ancestral_seqs(self, locus=None):
        """computes the conditional probabilities of each state for each node
        in the tree.

        Parameters
        ----------
        locus
            a named locus

        Returns
        -------
        {node_name: DictArray, ...}

        Notes
        -----
        Alignment columns are rows in the DictArray.
        """
        result = {}
        array_template = None
        for restricted_edge in self._tree.get_edge_vector():
            if restricted_edge.istip():
                continue
            try:
                r = []
                for motif in range(len(self._motifs)):
                    self.set_param_rule(
                        "fixed_motif",
                        value=motif,
                        edge=restricted_edge.name,
                        locus=locus,
                        is_constant=True,
                    )
                    likelihoods = self.get_full_length_likelihoods(locus=locus)
                    r.append(likelihoods)
                    if array_template is None:
                        array_template = DictArrayTemplate(
                            likelihoods.shape[0], self._motifs
                        )
            finally:
                self.set_param_rule(
                    "fixed_motif",
                    value=-1,
                    edge=restricted_edge.name,
                    locus=locus,
                    is_constant=True,
                )
            # dict of site x motif arrays
            result[restricted_edge.name] = array_template.wrap(
                numpy.transpose(numpy.asarray(r))
            )
        return result

    def likely_ancestral_seqs(self, locus=None) -> ArrayAlignment:
        """Returns the most likely reconstructed ancestral sequences as an
        alignment.

        Parameters
        ----------
        locus
            a named locus
        """
        prob_array = self.reconstruct_ancestral_seqs(locus=locus)
        seqs = []
        for edge, probs in list(prob_array.items()):
            seq = []
            for row in probs:
                by_p = [(p, state) for state, p in list(row.items())]
                seq.append(max(by_p)[1])
            seqs += [(edge, self.model.moltype.make_seq("".join(seq)))]
        return ArrayAlignment(data=seqs, moltype=self.model.moltype)

    def get_bin_probs(self, locus=None):
        hmm = self.get_param_value("bindex", locus=locus)
        lhs = [
            self.get_param_value("lh", locus=locus, bin=bin) for bin in self.bin_names
        ]
        array = hmm.get_posterior_probs(*lhs)
        return DictArrayTemplate(self.bin_names, array.shape[1]).wrap(array)

    def _valuesForDimension(self, dim):
        # in support of __str__
        if dim == "edge":
            result = [e.name for e in self._tree.get_edge_vector()]
        elif dim == "bin":
            result = self.bin_names[:]
        elif dim == "locus":
            result = self.locus_names[:]
        elif dim.startswith("motif"):
            result = self._mprob_motifs
        elif dim == "position":
            result = self.posn_names[:]
        else:
            raise KeyError(dim)
        return result

    def _valuesForDimensions(self, dims):
        # in support of __str__
        result = [[]]
        for dim in dims:
            new_result = []
            for r in result:
                for cat in self._valuesForDimension(dim):
                    new_result.append(r + [cat])
            result = new_result
        return result

    def _for_display(self):
        """processes statistics tables for display"""
        title = self.name or "Likelihood function statistics"
        result = []
        result += self.get_statistics(with_motif_probs=True, with_titles=True)
        for i, table_ in enumerate(result):
            if (
                "motif" in table_.title
                and table_.shape[1] == 2
                and table_.shape[0] >= 60
            ):  # just sort codon motif probs, then truncate
                table_ = table_.sorted(columns="motif")
                table_.set_repr_policy(head=5, tail=5, show_shape=False)
                result[i] = table_
        return title, result

    def _repr_html_(self):
        """for jupyter notebook display"""
        try:
            lnL = f"<p>log-likelihood = {self.get_log_likelihood():.4f}</p>"
        except ValueError:
            # alignment probably not yet set
            lnL = ""

        nfp = "<p>number of free parameters = %d</p>" % self.get_num_free_params()
        title, results = self._for_display()
        for i, table_ in enumerate(results):
            table_.title = table_.title.capitalize()
            table_.set_repr_policy(show_shape=False)
            results[i] = table_._repr_html_()
        results = [f"<h4>{title}</h4>", lnL, nfp] + results
        return "\n".join(results)

    def __repr__(self):
        return str(self)

    def __str__(self):
        title, results = self._for_display()

        try:
            lnL = f"log-likelihood = {self.get_log_likelihood():.4f}"
        except ValueError:
            # alignment probably not yet set
            lnL = None

        nfp = "number of free parameters = %d" % self.get_num_free_params()
        for table_ in results:
            table_.title = ""

        results = [title, lnL, nfp] + results if lnL else [title, nfp] + results
        return "\n".join(map(str, results))

    def get_annotated_tree(self, length_as=None):
        """returns tree with model attributes on node.params

        length_as : str or None
            replaces 'length' param with either 'ENS' or 'paralinear'.
            'ENS' is the expected number of substitution, (which will be
            different to standard length if the substitution model is
            non-stationary). 'paralinear' is the measure of Lake 1994.

        The other measures are always available in the params dict of each
        node.
        """
        from cogent3.evolve.ns_substitution_model import (
            DiscreteSubstitutionModel,
        )

        is_discrete = isinstance(self.model, DiscreteSubstitutionModel)

        if is_discrete and not length_as == "paralinear":
            raise ValueError(f"{length_as} invalid for discrete time process")

        assert length_as in ("ENS", "paralinear", None)
        d = self.get_param_value_dict(["edge"])
        lengths = d.pop("length", None)
        mprobs = self.get_motif_probs_by_node()
        if not is_discrete:
            ens = self.get_lengths_as_ens(motif_probs=mprobs)

        plin = self.get_paralinear_metric(motif_probs=mprobs)
        if length_as == "ENS":
            lengths = ens
        elif length_as == "paralinear":
            lengths = plin

        tree = self._tree.deepcopy()
        for edge in tree.get_edge_vector():
            if edge.name == "root":
                edge.params["mprobs"] = mprobs[edge.name].to_dict()
                continue

            if not is_discrete:
                edge.params["ENS"] = ens[edge.name]

            edge.params["length"] = lengths[edge.name]
            edge.params["paralinear"] = plin[edge.name]
            edge.params["mprobs"] = mprobs[edge.name].to_dict()
            for par in d:
                val = d[par][edge.name]
                if par == length_as:
                    val = ens[edge.name]
                edge.params[par] = val

        return tree

    def get_motif_probs(self, edge=None, bin=None, locus=None, position=None):
        """
        Parameters
        ----------
        edge : str
            name of edge
        bin : int or str
            name of bin
        locus : str
            name of locus
        position : int or str
            name of position

        Returns
        -------
        If 1D, returns DictArray, else a dict of DictArray
        """
        param_names = self.get_param_names()
        mprob_name = [n for n in param_names if "mprob" in n][0]
        dims = tuple(self.get_used_dimensions(mprob_name))
        mprobs = self.get_param_value_dict(dimensions=dims, params=[mprob_name])
        if len(dims) == 2:
            var = [c for c in dims if c != mprob_name][0]
            key = locals().get(var, None)
            mprobs = mprobs[mprob_name]
            if key is not None:
                mprobs = mprobs.get(str(key), mprobs.get(key))
                mprobs = {mprob_name: mprobs}

        # these can fall below the minimum allowed value due to
        # rounding errors, so I adjust these
        for k, value in mprobs.items():
            value.array = adjusted_gt_minprob(value.array, minprob=1e-6)

        if len(mprobs) == 1:
            mprobs = mprobs[mprob_name]

        return mprobs

    def get_bin_prior_probs(self, locus=None):
        bin_probs_array = self.get_param_value("bprobs", locus=locus)
        return DictArrayTemplate(self.bin_names).wrap(bin_probs_array)

    def get_scaled_lengths(self, predicate, bin=None, locus=None):
        """A dictionary of {scale:{edge:length}}"""
        if not hasattr(self._model, "get_scaled_lengths_from_Q"):
            return {}

        get_value_of = self.get_param_value
        value_of_kw = dict(locus=locus)

        if bin is None:
            bin_names = self.bin_names
        else:
            bin_names = [bin]

        if len(bin_names) == 1:
            bprobs = [1.0]
        else:
            bprobs = get_value_of("bprobs", **value_of_kw)

        mprobs = [get_value_of("mprobs", bin=b, **value_of_kw) for b in bin_names]

        scaled_lengths = {}
        for edge in self._tree.get_edge_vector():
            if edge.isroot():
                continue
            Qs = [
                get_value_of("Qd", bin=b, edge=edge.name, **value_of_kw).Q
                for b in bin_names
            ]
            length = get_value_of("length", edge=edge.name, **value_of_kw)
            scaled_lengths[edge.name] = length * self._model.get_scale_from_Qs(
                Qs, bprobs, mprobs, predicate
            )
        return scaled_lengths

    def get_paralinear_metric(self, motif_probs=None):
        """returns {edge.name: paralinear, ...}
        Parameters
        ----------
        motif_probs : dict or DictArray
            an item for each edge of the tree. Computed if not provided.
        """
        from cogent3.evolve.ns_substitution_model import (
            DiscreteSubstitutionModel,
        )

        is_discrete = isinstance(self.model, DiscreteSubstitutionModel)

        if motif_probs is None:
            motif_probs = self.get_motif_probs_by_node()
        plin = {}
        for edge in self.tree.get_edge_vector(include_root=False):
            parent_name = edge.parent.name
            pi = motif_probs[parent_name]
            P = self.get_psub_for_edge(edge.name)
            if is_discrete:
                para = paralinear_discrete_time(P.array, pi.array)
            else:
                Q = self.get_rate_matrix_for_edge(edge.name, calibrated=False)
                para = paralinear_continuous_time(P.array, pi.array, Q.array)

            plin[edge.name] = para

        return plin

    def get_lengths_as_ens(self, motif_probs=None):
        """returns {edge.name: ens, ...} where ens is the expected number of substitutions

        for a stationary Markov process, this is just branch length
        Parameters
        ----------
        motif_probs : dict or DictArray
            an item for each edge of the tree. Computed if not provided.
        """
        if motif_probs is None:
            motif_probs = self.get_motif_probs_by_node()
        node_names = [n for n in self.tree.get_node_names() if n != "root"]
        lengths = {e: self.get_param_value("length", edge=e) for e in node_names}
        if not isinstance(self.model, substitution_model.Stationary):
            ens = {}
            for e in node_names:
                Q = self.get_rate_matrix_for_edge(e)
                length = expected_number_subs(motif_probs[e], Q, lengths[e])
                ens[e] = length

            lengths = ens

        return lengths

    def get_param_rules(self):
        """returns the [{rule}, ..] that would allow reconstruction"""
        # markov model rate terms
        rules = []
        param_names = self.get_param_names()
        for param_name in param_names:
            defn = self.defn_for[param_name]
            try:
                rules.extend(defn.get_param_rules())
            except AttributeError:
                # aggregate params, like those deriving from gamma shaped rates
                pass

        return rules

    def get_statistics(self, with_motif_probs=True, with_titles=True):
        """returns the parameter values as tables/dict

        Parameters
        ----------
        with_motif_probs
            include the motif probability table
        with_titles
            include a title for each table based on it's
            dimension

        """
        result = []
        group = {}
        param_names = self.get_param_names()

        mprob_name = [n for n in param_names if "mprob" in n]
        mprob_name = mprob_name[0] if mprob_name else ""
        if not with_motif_probs:
            param_names.remove(mprob_name)

        for param in param_names:
            dims = tuple(self.get_used_dimensions(param))
            if dims not in group:
                group[dims] = []
            group[dims].append(param)
        table_order = list(group.keys())
        table_order.sort()
        for table_dims in table_order:
            raw_table = self.get_param_value_dict(
                dimensions=table_dims, params=group[table_dims]
            )
            param_names = group[table_dims]
            param_names.sort()
            if table_dims == ("edge",):
                if "length" in param_names:
                    param_names.remove("length")
                    param_names.insert(0, "length")
                raw_table["parent"] = dict(
                    [
                        (e.name, e.parent.name)
                        for e in self._tree.get_edge_vector()
                        if not e.isroot()
                    ]
                )
                param_names.insert(0, "parent")
            list_table = []
            heading_names = list(table_dims) + param_names
            row_order = self._valuesForDimensions(table_dims)
            for scope in row_order:
                row = {}
                row_used = False
                for param in param_names:
                    d = raw_table[param]
                    try:
                        for part in scope:
                            d = d[part]
                    except KeyError:
                        d = "NA"
                    else:
                        row_used = True
                    row[param] = d
                if row_used:
                    row.update(dict(list(zip(table_dims, scope))))
                    row = [row[k] for k in heading_names]
                    list_table.append(row)
            if table_dims:
                title = ["", f"{' '.join(table_dims)} params"][with_titles]
            else:
                title = ["", "global params"][with_titles]
            row_ids = None
            stat_table = table.Table(
                heading_names,
                list_table,
                max_width=80,
                index_name=row_ids,
                title=title,
                **self._format,
            )
            if group[table_dims] == [mprob_name]:
                # if stat_table.shape
                # if mprobs, we use the motifs as header
                motifs = list(sorted(set(stat_table.tolist("motif"))))
                if stat_table.shape[1] == 2:
                    motif_prob = dict(stat_table.tolist())
                    heading_names = motifs
                    list_table = [motif_prob[m] for m in motifs]
                    list_table = [list_table]
                elif stat_table.shape[1] == 3:
                    rows = []
                    other_col = [
                        c
                        for c in stat_table.header
                        if "motif" not in c and "mprobs" not in c
                    ][0]
                    for val in stat_table.distinct_values(other_col):
                        subtable = stat_table.filtered(
                            lambda x: x == val, columns=other_col
                        )
                        motif_prob = dict(
                            subtable.tolist(
                                [c for c in stat_table.header if c != other_col]
                            )
                        )
                        rows.append([val] + [motif_prob[m] for m in motifs])
                    heading_names = [other_col] + motifs
                    list_table = rows
                stat_table = table.Table(
                    heading_names, list_table, max_width=80, title=title, **self._format
                )

            result.append(stat_table)
        return result

    def to_rich_dict(self):
        """returns detailed info on object, used by to_json"""
        data = deepcopy(self._serialisable)
        for key in ("model", "tree"):
            del data[key]

        tree = self.tree.to_rich_dict()
        edge_attr = tree["edge_attributes"]
        for edge in edge_attr:
            if edge == "root":
                continue
            try:
                edge_attr[edge]["length"] = self.get_param_value("length", edge=edge)
            except KeyError:
                # probably discrete-time model
                edge_attr[edge]["length"] = None

        model = self._model.to_rich_dict(for_pickle=False)

        aln_defn = self.defn_for["alignment"]
        if len(aln_defn.index) == 1:
            alignment = self.get_param_value("alignment").to_rich_dict()
            mprobs = self.get_motif_probs().to_dict()
        else:
            # this is a multi-locus likelihood function
            alignment = {a["locus"]: a["value"] for a in aln_defn.get_param_rules()}
            for k in alignment:
                alignment[k] = alignment[k].to_rich_dict()

            mprobs = self.get_motif_probs()
            if isinstance(mprobs, dict):
                # separate mprobs per locus
                for k in alignment:
                    mprobs[k] = mprobs[k].to_dict()
            else:
                # motif probs are constrained to be the same between loci
                mprobs = self.get_motif_probs().to_dict()

        DLC = self.all_psubs_DLC()
        try:
            unique_Q = self.all_rate_matrices_unique()
        except Exception:
            # there's a mix of assertions
            # for "storage", make this indeterminate in those cases
            unique_Q = None

        data = dict(
            model=model,
            tree=tree,
            alignment=alignment,
            likelihood_construction=data,
            param_rules=self.get_param_rules(),
            lnL=self.get_log_likelihood(),
            nfp=self.get_num_free_params(),
            motif_probs=mprobs,
            DLC=DLC,
            unique_Q=unique_Q,
            type=get_object_provenance(self),
            name=self.get_name(),
            version=__version__,
        )
        return data

    def to_json(self):
        data = self.to_rich_dict()
        data = json.dumps(data)
        return data

    @property
    def name(self):
        if self._name is None:
            self._name = self.model.name or ""

        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    # For tests.  Compat with old LF interface
    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def set_tables_format(self, space=4, digits=4):
        """sets display properties for statistics tables. This affects results
        of str(lf) too."""
        space = [space, 4][type(space) != int]
        digits = [digits, 4][type(digits) != int]
        self._format = dict(space=space, digits=digits)

    def _get_motif_probs_by_node_tr(self, edges=None, bin=None, locus=None):
        """returns motif probs by node for time-reversible models"""
        mprob_rules = [r for r in self.get_param_rules() if "mprob" in r["par_name"]]
        if len(mprob_rules) > 1 or self.model.mprob_model == "monomers":
            raise NotImplementedError

        mprobs = self.get_motif_probs()
        if len(mprobs) != len(self.motifs):
            # a Muse and Gaut model
            expanded = numpy.zeros(len(self.motifs), dtype=float)
            for i, motif in enumerate(self.motifs):
                val = 1.0
                for b in motif:
                    val *= mprobs[b]
                expanded[i] = val
            mprobs = expanded / expanded.sum()
        else:
            mprobs = [mprobs[m] for m in self.motifs]
        edges = []
        values = []
        for e in self.tree.postorder():
            edges.append(e.name)
            values.append(mprobs)

        return DictArrayTemplate(edges, self.motifs).wrap(values)

    def get_motif_probs_by_node(self, edges=None, bin=None, locus=None):
        from cogent3.evolve.substitution_model import TimeReversible

        if isinstance(self.model, TimeReversible):
            return self._get_motif_probs_by_node_tr(edges=edges, bin=bin, locus=locus)

        kw = dict(bin=bin, locus=locus)
        mprobs = self.get_param_value("mprobs", **kw)
        mprobs = self._model.calc_word_probs(mprobs)
        result = self._nodeMotifProbs(self._tree, mprobs, kw)
        if edges is None:
            edges = [name for (name, m) in result]
        result = dict(result)
        values = [result[name] for name in edges]
        return DictArrayTemplate(edges, self._mprob_motifs).wrap(values)

    def _nodeMotifProbs(self, tree, mprobs, kw):
        result = [(tree.name, mprobs)]
        for child in tree.children:
            psub = self.get_psub_for_edge(child.name, **kw)
            child_mprobs = numpy.dot(mprobs, psub)
            result.extend(self._nodeMotifProbs(child, child_mprobs, kw))
        return result

    def simulate_alignment(
        self,
        sequence_length=None,
        random_series=None,
        exclude_internal=True,
        locus=None,
        seed=None,
        root_sequence=None,
    ):
        """
        Returns an alignment of simulated sequences with key's corresponding to
        names from the current attached alignment.

        Parameters
        ----------
        sequence_length
            the legnth of the alignment to be simulated,
            default is the length of the attached alignment.
        random_series
            a random number generator.
        exclude_internal
            if True, only sequences for tips are returned.
        locus
            if fit to multiple alignments, select the values corresponding to
            locus for generating data
        seed
            seed value for the random number generator
        root_sequence
            a sequence from which all others evolve
        """
        orig_ambig = {}
        if sequence_length is None:
            lht = self.get_param_value("lht", locus=locus)
            try:
                sequence_length = len(lht.index)
            except AttributeError:
                raise ValueError(
                    "Must provide sequence_length since no alignment set on self"
                )

            leaves = self.get_param_value("leaf_likelihoods", locus=locus)
            for (seq_name, leaf) in list(leaves.items()):
                orig_ambig[seq_name] = leaf.get_ambiguous_positions()

        if random_series is None:
            random_series = random.Random()
            random_series.seed(seed)

        def psub_for(edge, bin):
            return self.get_psub_for_edge(edge, bin=bin, locus=locus)

        if len(self.bin_names) > 1:
            hmm = self.get_param_value("bdist", locus=locus)
            site_bins = hmm.emit(sequence_length, random_series)
        else:
            site_bins = numpy.zeros([sequence_length], int)

        evolver = AlignmentEvolver(
            random_series,
            orig_ambig,
            exclude_internal,
            self.bin_names,
            site_bins,
            psub_for,
            self._motifs,
        )

        if root_sequence is not None:  # we convert to a vector of motifs
            if isinstance(root_sequence, str):
                root_sequence = self._model.moltype.make_seq(root_sequence)
            motif_len = self._model.get_alphabet().get_motif_len()
            root_sequence = root_sequence.get_in_motif_size(motif_len)
        else:
            mprobs = self.get_param_value("mprobs", locus=locus, edge="root")
            mprobs = self._model.calc_word_probs(mprobs)
            mprobs = dict((m, p) for (m, p) in zip(self._motifs, mprobs))
            root_sequence = random_sequence(random_series, mprobs, sequence_length)

        simulated_sequences = evolver(self._tree, root_sequence)

        return ArrayAlignment(data=simulated_sequences, moltype=self._model.moltype)

    def all_psubs_DLC(self):
        """Returns True if every Psub matrix is Diagonal Largest in Column"""
        all_psubs = self.get_all_psubs()
        for P in all_psubs.values():
            if (P.to_array().diagonal() < P).any():
                return False
        return True

    def all_rate_matrices_unique(self):
        """Returns True if every rate matrix is unique for its Psub matrix"""
        # get all possible Q, as products of t, and any rate-het terms
        all_Q = self.get_all_rate_matrices(calibrated=False)
        for Q in all_Q.values():
            Q = Q.to_array()
            if not is_generator_unique(Q):
                return False
        return True

    def initialise_from_nested(self, nested_lf):
        from cogent3.evolve.substitution_model import Stationary

        assert (
            self.get_num_free_params() > nested_lf.get_num_free_params()
        ), "wrong order for likelihood functions"
        compatible_likelihood_functions(self, nested_lf)

        same = (
            isinstance(self.model, Stationary)
            and isinstance(nested_lf.model, Stationary)
        ) or (
            not isinstance(self.model, Stationary)
            and not isinstance(nested_lf.model, Stationary)
        )

        mprobs = nested_lf.get_motif_probs()
        edge_names = self.tree.get_node_names()
        edge_names.remove("root")
        param_proj = _ParamProjection(nested_lf.model, self.model, mprobs, same=same)
        param_rules = nested_lf.get_param_rules()
        param_rules = param_proj.update_param_rules(param_rules)
        my_rules = self.get_param_rules()
        my_rules = update_scoped_rules(my_rules, param_rules)
        self.apply_param_rules(my_rules)
        return
