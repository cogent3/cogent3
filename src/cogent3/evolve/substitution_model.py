#!/usr/bin/env python
"""
substitution_model.py

Contains classes for defining Markov models of substitution.
These classes depend on an Alphabet class member for defining the set
of motifs that each represent a state in the Markov chain. Examples of
a 'dna' type alphabet motif is 'a', and of a 'codon' type motif is'atg'.

By default all models include the gap motif ('-' for a 'dna' alphabet or
'---' for a 'codon' alphabet). This differs from software such as PAML,
where gaps are treated as ambiguituous states (specifically, as 'n'). The gap
motif state can be excluded from the substitution model using the method
excludeGapMotif(). It is recommended that to ensure the alignment and the
substitution model are defined with the same alphabet that modifications
are done to the substitution model alphabet and this instance is then given
to the alignment.

The model's substitution rate parameters are represented as a dictionary
with the parameter names as keys, and predicate functions as the values.
These predicate functions compare a pair of motifs, returning True or False.
Many such functions are provided as methods of the class. For instance,
the istransition method is pertinent to dna based models. This method returns
True if an 'a'/'g' or 'c'/'t' pair is passed to it, False otherwise. In this
way the positioning of parameters in the instantaneous rate matrix (commonly
called Q) is determined.

>>> model = Nucleotide(equal_motif_probs=True)
>>> model.setparameterrules({'alpha': model.istransition})
>>> parameter_controller = model.make_likelihood_function(tree)
"""

import json
import warnings

from collections.abc import Callable
from copy import deepcopy

import numpy

from numpy.linalg import svd

from cogent3.core import moltype
from cogent3.evolve import motif_prob_model, parameter_controller, predicate
from cogent3.evolve.likelihood_tree import make_likelihood_tree_leaf
from cogent3.evolve.substitution_calculation import (
    AlignmentAdaptDefn,
    ExpDefn,
    LengthDefn,
)
from cogent3.evolve.substitution_calculation import (
    SubstitutionParameterDefn as ParamDefn,
)
from cogent3.recalculation.definition import (
    CalcDefn,
    CallDefn,
    ConstDefn,
    GammaDefn,
    MonotonicDefn,
    NonParamDefn,
    PartitionDefn,
    ProductDefn,
    SelectForDimension,
    WeightedPartitionDefn,
)
from cogent3.util.misc import extend_docstring_from, get_object_provenance


__author__ = "Peter Maxwell, Gavin Huttley and Andrew Butterfield"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__contributors__ = [
    "Gavin Huttley",
    "Andrew Butterfield",
    "Peter Maxwell",
    "Matthew Wakefield",
    "Brett Easton",
    "Rob Knight",
    "Von Bing Yap",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def predicate2matrix(alphabet, pred, mask=None):
    """From a test like istransition() produce an MxM boolean matrix"""
    M = len(alphabet)
    result = numpy.zeros([M, M], int)
    for i in range(M):
        for j in range(M):
            if mask is None or mask[i, j]:
                result[i, j] = pred(alphabet[i], alphabet[j])
    return result


def redundancy_in_predicate_masks(preds):
    # Calculate the nullity of the predicates.  If non-zero
    # there is some redundancy and the model will be overparameterised.
    if len(preds) <= 1:
        return 0
    eqns = 1.0 * numpy.array([list(mask.flat) for mask in list(preds.values())])
    svs = svd(eqns)[1]
    # count non-duplicate non-zeros singular values
    matrix_rank = len([sv for sv in svs if abs(sv) > 1e-8])
    return len(preds) - matrix_rank


def _maxWidthIfTruncated(pars, delim, each):
    # 'pars' is an array of lists of strings, how long would the longest
    # list representation be if the strings were truncated at 'each'
    # characters and joined together with 'delim'.
    return max(
        [
            sum([min(len(par), each) for par in par_list])
            + len(delim) * (len(par_list) - 1)
            for par_list in pars.flat
        ]
    )


def _isSymmetrical(matrix):
    return numpy.alltrue(numpy.alltrue(matrix == numpy.transpose(matrix)))


class _SubstitutionModel(object):
    # Subclasses must provide
    #  .make_param_controller_defns()

    def __init__(
        self,
        alphabet,
        motif_probs=None,
        optimise_motif_probs=False,
        equal_motif_probs=False,
        motif_probs_from_data=None,
        motif_probs_alignment=None,
        mprob_model=None,
        model_gaps=False,
        recode_gaps=False,
        motif_length=1,
        name="",
        motifs=None,
    ):
        """
        Parameters
        ----------
        alphabet
            An Alphabet object
        motif_probs
            Dictionary of probabilities.
        optimise_motif_probs: bool
            Treat like other free parameters. Any values set by the other
            motif_prob options will be used as initial values.
        equal_motif_probs: bool
            Flag to set alignment motif probs equal.
        motif_probs_from_data: bool
            Get motif probabilities from data provided to likelihood function.
        motif_probs_alignment
            An alignment from which motif probs are set.
        mprob_model: str
            'tuple', 'conditional', 'monomer' or 'monomers' to specify how
           tuple-alphabet (including codon) motif probs are used.
        model_gaps: bool
            Whether the gap motif should be included as a state.
        recode_gaps: bool
            Whether gaps in an alignment should be treated as an ambiguous
            state instead.
        motif_length: int
            Based on 'alphabet', uses a tuple alphabet where individual words
            have motif_length number of characters.
        name: str
            Name of this model
        motifs
            Use a subalphabet that only contains those motifs.
        """
        d = locals()
        exclude = ("self", "__class__")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}
        # MISC
        assert len(alphabet) < 65, (
            "Alphabet too big. Try explicitly " "setting alphabet to PROTEIN or DNA"
        )

        self.name = name
        self._optimise_motif_probs = optimise_motif_probs

        # ALPHABET
        if recode_gaps:
            if model_gaps:
                warnings.warn("Converting gaps to wildcards AND modeling gaps")
            else:
                model_gaps = False

        self.recode_gaps = recode_gaps

        self.moltype = alphabet.moltype
        if model_gaps:
            alphabet = alphabet.with_gap_motif()

        if motif_length > 1:
            alphabet = alphabet.get_word_alphabet(motif_length)

        if motifs is not None:
            alphabet = alphabet.get_subset(motifs)
        self.alphabet = alphabet
        self.gapmotif = alphabet.get_gap_motif()
        self._word_length = alphabet.get_motif_len()

        # MOTIF PROB ALPHABET MAPPING
        if mprob_model is None:
            mprob_model = "tuple" if self._word_length == 1 else "conditional"
        elif mprob_model == "word":
            mprob_model = "tuple"

        if model_gaps and mprob_model != "tuple":
            raise ValueError("mprob_model must be 'tuple' to model gaps")

        isinst = self._is_instantaneous
        self._instantaneous_mask = predicate2matrix(self.alphabet, isinst)
        self._instantaneous_mask_f = self._instantaneous_mask * 1.0
        self._mprob_model = mprob_model
        self.mprob_model = motif_prob_model.make_model(
            mprob_model, alphabet, self._instantaneous_mask_f
        )

        # MOTIF PROBS
        if equal_motif_probs:
            assert not (
                motif_probs or motif_probs_alignment
            ), "Motif probs equal or provided but not both"
            motif_probs = self.mprob_model.make_equal_motif_probs()
        elif motif_probs_alignment is not None:
            assert (
                not motif_probs
            ), "Motif probs from alignment or provided but not both"
            motif_probs = self.count_motifs(motif_probs_alignment)
            motif_probs = motif_probs.astype(float) / sum(motif_probs)
            assert len(alphabet) == len(motif_probs)
            motif_probs = dict(list(zip(alphabet, motif_probs)))
        if motif_probs:
            self.adapt_motif_probs(motif_probs)  # to check
            self.motif_probs = motif_probs
            if motif_probs_from_data is None:
                motif_probs_from_data = False
        else:
            self.motif_probs = None
            if motif_probs_from_data is None:
                motif_probs_from_data = True
        self.motif_probs_from_align = motif_probs_from_data

    def __getnewargs_ex__(self, *args, **kw):
        data = self.to_rich_dict(for_pickle=True)
        return (), data

    def to_rich_dict(self, for_pickle=False):
        data = deepcopy(self._serialisable)
        if not for_pickle:
            for key, value in data.items():
                type_ = get_object_provenance(value)
                if type_.startswith("cogent3"):
                    try:
                        value = value.to_rich_dict(for_pickle=False)
                    except AttributeError:
                        pass
                    finally:
                        data[key] = value
            if "predicates" in data and data["predicates"]:
                data["predicates"] = [str(p) for p in data["predicates"]]
            data["type"] = get_object_provenance(self)
            data["version"] = __version__
        return data

    def to_json(self):
        """returns result of json formatted string"""
        data = self.to_rich_dict(for_pickle=False)
        return json.dumps(data)

    def get_param_list(self):
        return []

    def __repr__(self):
        s = []
        s.append(f"name={getattr(self, 'name', None)!r};")
        if hasattr(self, "predicate_masks"):
            parlist = list(self.predicate_masks.keys())
            s.append(f"params={parlist};")
        motifs = self.get_motifs()
        s.append(f"num_motifs={len(motifs)};")
        s.append(f"motifs={motifs})")
        return f"{self.__class__.__name__}({' '.join(s)})"

    def get_alphabet(self):
        return self.alphabet

    def get_mprob_alphabet(self):
        return self.mprob_model.get_input_alphabet()

    def get_motifs(self):
        return list(self.get_alphabet())

    @property
    def word_length(self):
        return self._word_length

    def get_motif_probs(self):
        """Return the dictionary of motif probabilities."""
        return self.motif_probs.copy()

    def set_param_controller_motif_probs(self, pc, mprobs, **kw):
        return self.mprob_model.set_param_controller_motif_probs(pc, mprobs, **kw)

    def make_likelihood_function(
        self,
        tree,
        motif_probs_from_align=None,
        optimise_motif_probs=None,
        aligned=True,
        expm=None,
        digits=None,
        space=None,
        **kw,
    ):

        if motif_probs_from_align is None:
            motif_probs_from_align = self.motif_probs_from_align

        if optimise_motif_probs is None:
            optimise_motif_probs = self._optimise_motif_probs

        kw["optimise_motif_probs"] = optimise_motif_probs
        kw["motif_probs_from_align"] = motif_probs_from_align

        if aligned:
            klass = parameter_controller.AlignmentLikelihoodFunction
        else:
            alphabet = self.get_alphabet()
            assert alphabet.get_gap_motif() not in alphabet
            klass = parameter_controller.SequenceLikelihoodFunction

        result = klass(self, tree, **kw)

        if self.motif_probs is not None:
            result.set_motif_probs(
                self.motif_probs, is_constant=not optimise_motif_probs, auto=True
            )

        if expm is None:
            expm = self._default_expm_setting
        if expm is not None:
            result.set_expm(expm)

        if digits or space:
            result.set_tables_format(digits=digits, space=space)

        return result

    def convert_alignment(self, alignment):
        # this is to support for everything but HMM
        result = {}
        for seq_name in alignment.names:
            sequence = alignment.get_gapped_seq(seq_name, self.recode_gaps)
            result[seq_name] = self.convert_sequence(sequence, seq_name)
        return result

    def convert_sequence(self, sequence, name):
        # make_likelihood_tree_leaf, sort of an indexed profile where duplicate
        # columns stored once, so likelihoods only calc'd once
        return make_likelihood_tree_leaf(sequence, self.get_alphabet(), name)

    def count_motifs(self, alignment, include_ambiguity=False):
        return self.mprob_model.count_motifs(
            alignment, include_ambiguity, self.recode_gaps
        )

    def make_alignment_defn(self, model):
        align = NonParamDefn("alignment", ("locus",))
        # The name of this matters, it's used in likelihood_function.py
        # to retrieve the correct (adapted) alignment.
        return AlignmentAdaptDefn(model, align)

    def adapt_motif_probs(self, motif_probs, auto=False):
        return self.mprob_model.adapt_motif_probs(motif_probs, auto=auto)

    def calc_monomer_probs(self, word_probs):
        # Not presently used, always go monomer->word instead
        return self.mprob_model.calc_monomer_probs(word_probs)

    def calc_word_probs(self, monomer_probs):
        return self.mprob_model.calc_word_probs(monomer_probs)

    def calc_word_weight_matrix(self, monomer_probs):
        return self.mprob_model.calc_word_weight_matrix(monomer_probs)

    def make_param_controller_defns(self, bin_names, endAtQd=False):
        (
            input_probs,
            word_probs,
            mprobs_matrix,
        ) = self.mprob_model.make_motif_word_prob_defns()

        if len(bin_names) > 1:
            bprobs = PartitionDefn(
                [1.0 / len(bin_names) for bin in bin_names],
                name="bprobs",
                dimensions=["locus"],
                dimension=("bin", bin_names),
            )
        else:
            bprobs = None

        defns = {
            "align": self.make_alignment_defn(ConstDefn(self, "model")),
            "bprobs": bprobs,
            "word_probs": word_probs,
        }

        rate_params = self.make_rate_params(bprobs)
        if endAtQd:
            defns["Qd"] = self.make_Qd_defn(word_probs, mprobs_matrix, rate_params)
        else:
            defns["psubs"] = self.make_psubs_defn(
                bprobs, word_probs, mprobs_matrix, rate_params
            )
        return defns


def non_zero_coords(matrix):
    dim = matrix.shape[0]
    return [(i, j) for i in range(dim) for j in range(dim) if matrix[i, j] != 0]


class _ContinuousSubstitutionModel(_SubstitutionModel):
    # subclass must provide:
    #
    # - parameter_order: a list of parameter names corresponding to the
    #   arguments of:
    #
    # - calc_exchangeability_matrix(*params)
    #   convert len(self.parameter_order) params to a matrix

    """A substitution model for which the rate matrix (P) is derived from an
    instantaneous rate matrix (Q).  The nature of the parameters used to define
    Q is up to the subclasses.
    """

    # At some point this can be made variable, and probably
    # the default changed to False
    long_indels_are_instantaneous = True

    _exponentiator = None
    _default_expm_setting = "either"

    @extend_docstring_from(_SubstitutionModel.__init__)
    def __init__(
        self,
        alphabet,
        with_rate=False,
        ordered_param=None,
        distribution=None,
        partitioned_params=None,
        **kw,
    ):
        """
        with_rate: bool
            Add a 'rate' parameter which varies by bin.
        ordered_param: str
            name of a single parameter which distinguishes any bins.
        distribution: str
            choices of 'free' or 'gamma' or an instance of some distribution
        partitioned_params
            names of params to be partitioned across bins
        kw
        """

        _SubstitutionModel.__init__(self, alphabet, **kw)
        d = locals()
        exclude = ("self", "__class__")
        d = {k: v for k, v in d.items() if k not in exclude}
        self._serialisable.update(d)
        alphabet = self.get_alphabet()  # as may be altered by recode_gaps etc.

        # BINS
        if not ordered_param:
            if ordered_param is not None:
                warnings.warn("ordered_param should be a string or None")
                ordered_param = None
            if distribution:
                if with_rate:
                    ordered_param = "rate"
                else:
                    raise ValueError("distribution provided without ordered_param")
        elif not isinstance(ordered_param, str):
            warnings.warn("ordered_param should be a string or None")
            assert len(ordered_param) == 1, "More than one ordered_param"
            ordered_param = ordered_param[0]
            assert ordered_param, "False value hidden in list"
        self.ordered_param = ordered_param

        if distribution == "gamma":
            distribution = GammaDefn
        elif distribution in [None, "free"]:
            distribution = MonotonicDefn
        elif isinstance(distribution, str):
            raise ValueError(f'Unknown distribution "{distribution}"')
        self.distrib_class = distribution

        if not partitioned_params:
            partitioned_params = ()
        elif isinstance(partitioned_params, str):
            partitioned_params = (partitioned_params,)
        else:
            partitioned_params = tuple(partitioned_params)
        if self.ordered_param:
            if self.ordered_param not in partitioned_params:
                partitioned_params += (self.ordered_param,)
        self.partitioned_params = partitioned_params

        if "rate" in partitioned_params:
            with_rate = True
        self.with_rate = with_rate

        # CACHED SHORTCUTS
        self._exponentiator = None
        # self._ident = numpy.identity(len(self.alphabet), float)

    def check_params_exist(self):
        """Raise an error if the parameters specified to be partitioned or
        ordered don't actually exist."""
        for param in self.partitioned_params:
            if param not in self.parameter_order and param != "rate":
                desc = ["partitioned", "ordered"][param == self.ordered_param]
                raise ValueError(f'{desc} param "{param}" unknown')

    def _is_instantaneous(self, x, y):
        diffs = sum([X != Y for (X, Y) in zip(x, y)])
        return diffs == 1 or (
            diffs > 1
            and self.long_indels_are_instantaneous
            and self._is_any_indel(x, y)
        )

    def _is_any_indel(self, x, y):
        """An indel of any length"""
        # Things get complicated when a contigous indel of any length is OK:
        if x == y:
            return False
        gap_start = gap_end = gap_strand = None
        for (i, (X, Y)) in enumerate(zip(x, y)):
            G = self.gapmotif[i]
            if X != Y:
                if X != G and Y != G:
                    return False  # non-gap differences had their chance above
                elif gap_start is None:
                    gap_start = i
                    gap_strand = [X, Y].index(G)
                elif gap_end is not None or [X, Y].index(G) != gap_strand:
                    return False  # can't start a second gap
                else:
                    pass  # extend open gap
            elif gap_start is not None:
                gap_end = i
        return True

    def calcQ(self, word_probs, mprobs_matrix, *params):
        Q = self.calc_exchangeability_matrix(word_probs, *params)
        row_totals = Q.sum(axis=1)
        Q -= numpy.diag(row_totals)
        Q *= 1.0 / (word_probs * row_totals).sum()
        return Q

    def get_reference_cell(self):
        """returns the reference cell of a given model"""
        dim = len(self.alphabet)
        mats = numpy.zeros((dim, dim), dtype=int)
        for m in self.predicate_masks.values():
            mats += m
        ref_mask = self._instantaneous_mask - mats
        return set(non_zero_coords(ref_mask))

    def get_param_matrix_coords(self, include_ref_cell=False):
        """returncoordinates for every predicate"""
        dim = len(self.alphabet)
        param_coords = {}
        for key, m in self.predicate_masks.items():
            coords = [(i, j) for i in range(dim) for j in range(dim) if m[i, j] != 0]
            coords = set(coords)
            param_coords[key] = coords

        if include_ref_cell:
            param_coords["ref_cell"] = self.get_reference_cell()
        return param_coords

    def make_Qd_defn(self, word_probs, mprobs_matrix, rate_params):
        """Diagonalized Q, ie: rate matrix prepared for exponentiation"""
        Q = CalcDefn(self.calcQ, name="Q")(word_probs, mprobs_matrix, *rate_params)
        expm = NonParamDefn("expm")
        exp = ExpDefn(expm)
        return CallDefn(exp, Q, name="Qd")

    def _make_bin_param_defn(self, edge_par_name, bin_par_name, bprob_defn):
        # if no ordered param defined, behaves as old, everything indexed by
        # and edge
        if edge_par_name not in self.partitioned_params:
            return ParamDefn(dimensions=["bin"], name=bin_par_name)

        if edge_par_name == self.ordered_param:
            whole = self.distrib_class(bprob_defn, bin_par_name)
        else:
            # this forces them to average to one, but no forced order
            # this means you can't force a param value to be shared across bins
            # so 1st above approach has to be used
            whole = WeightedPartitionDefn(bprob_defn, bin_par_name + "_partn")
        whole.bin_names = bprob_defn.bin_names
        return SelectForDimension(whole, "bin", name=bin_par_name)

    def make_rate_params(self, bprobs):
        params = []
        for param_name in self.parameter_order:
            if bprobs is None or param_name not in self.partitioned_params:
                defn = ParamDefn(param_name)
            else:
                e_defn = ParamDefn(param_name, dimensions=["edge", "locus"])
                # should be weighted by bprobs*rates not bprobs
                b_defn = self._make_bin_param_defn(
                    param_name, param_name + "_factor", bprobs
                )
                defn = ProductDefn(b_defn, e_defn, name=param_name + "_BE")
            params.append(defn)
        return params

    def make_fundamental_param_controller_defns(self, bin_names):
        """Everything one step short of the psubs, because cogent3.align code
        needs to handle Q*t itself."""
        defns = self.make_param_controller_defns(bin_names, endAtQd=True)
        assert "length" not in defns
        defns["length"] = LengthDefn()
        return defns

    def make_psubs_defn(self, bprobs, word_probs, mprobs_matrix, rate_params):
        distance = self.make_distance_defn(bprobs)
        return self.make_continuous_psub_defn(
            word_probs, mprobs_matrix, distance, rate_params
        )

    def make_distance_defn(self, bprobs):
        length = LengthDefn()
        if self.with_rate and bprobs is not None:
            b_rate = self._make_bin_param_defn("rate", "rate", bprobs)
            distance = ProductDefn(length, b_rate, name="distance")
        else:
            distance = length
        return distance

    def make_continuous_psub_defn(
        self, word_probs, mprobs_matrix, distance, rate_params
    ):
        Qd = self.make_Qd_defn(word_probs, mprobs_matrix, rate_params)
        return CallDefn(Qd, distance, name="psubs")


class StationaryQ:
    "Contains the Original Definition of calcQ"

    def calcQ(self, word_probs, mprobs_matrix, *params):
        Q = self.calc_exchangeability_matrix(word_probs, *params)
        Q *= mprobs_matrix
        row_totals = Q.sum(axis=1)
        Q -= numpy.diag(row_totals)
        Q *= 1.0 / (word_probs * row_totals).sum()
        return Q


class Empirical(StationaryQ, _ContinuousSubstitutionModel):
    """A continuous substitution model with a predefined instantaneous rate
    matrix."""

    @extend_docstring_from(_ContinuousSubstitutionModel.__init__)
    def __init__(self, alphabet, rate_matrix, **kw):
        """
        - rate_matrix: The instantaneous rate matrix
        """
        _ContinuousSubstitutionModel.__init__(self, alphabet, **kw)
        d = locals()
        exclude = ("self", "__class__")
        d = {k: v for k, v in d.items() if k not in exclude}
        self._serialisable.update(d)

        alphabet = self.get_alphabet()  # as may be altered by recode_gaps etc.
        N = len(alphabet)
        assert rate_matrix.shape == (N, N)
        assert numpy.alltrue(numpy.diagonal(rate_matrix) == 0)
        self._instantaneous_mask_f = rate_matrix * 1.0
        self._instantaneous_mask = self._instantaneous_mask_f != 0.0
        self.symmetric = _isSymmetrical(self._instantaneous_mask_f)
        self.parameter_order = []
        self.check_params_exist()

    def calc_exchangeability_matrix(self, mprobs):
        return self._instantaneous_mask_f.copy()


class Parametric(_ContinuousSubstitutionModel):
    """A continuous substitution model with only user-specified substitution
    parameters. This is a general process -- non-stationary and, if specified
    via predicates, non-reversible"""

    @extend_docstring_from(_ContinuousSubstitutionModel.__init__)
    def __init__(self, alphabet, predicates=None, scales=None, **kw):
        """
        predicates: dict
            a dict of {name:predicate}. See cogent3.evolve.predicate
        scales: dict
            scale rules, dict with predicates
        kw
        """
        self._canned_predicates = None
        _ContinuousSubstitutionModel.__init__(self, alphabet, **kw)

        d = locals()
        exclude = ("self", "__class__")
        d = {k: v for k, v in d.items() if k not in exclude}
        self._serialisable.update(d)

        (predicate_masks, predicate_order) = self._adapt_predicates(predicates or [])

        # Check for redundancy in predicates, ie: 1 or more than combine
        # to be equivalent to 1 or more others, or the distance params.
        # Give a clearer error in simple cases like always false or true.
        for (name, matrix) in list(predicate_masks.items()):
            if numpy.alltrue((matrix == 0).flat):
                raise ValueError(f"Predicate {name} is always false.")
        predicates_plus_scale = predicate_masks.copy()
        predicates_plus_scale[None] = self._instantaneous_mask
        for (name, matrix) in list(predicate_masks.items()):
            if numpy.alltrue((matrix == self._instantaneous_mask).flat):
                raise ValueError(f"Predicate {name} is always true.")
        if redundancy_in_predicate_masks(predicate_masks):
            raise ValueError("Redundancy in predicates.")
        if redundancy_in_predicate_masks(predicates_plus_scale):
            raise ValueError(
                "Some combination of predicates is"
                " equivalent to the overall rate parameter."
            )

        self.predicate_masks = predicate_masks
        self.parameter_order = []
        self.predicate_indices = []
        self.symmetric = _isSymmetrical(self._instantaneous_mask)
        for pred in predicate_order:
            mask = predicate_masks[pred]
            if not _isSymmetrical(mask):
                self.symmetric = False
            indices = numpy.nonzero(mask)
            assert numpy.alltrue(mask[indices] == 1)
            self.parameter_order.append(pred)
            self.predicate_indices.append(indices)
        (self.scale_masks, scale_order) = self._adapt_predicates(scales or [])
        self.check_params_exist()

    def calc_exchangeability_matrix(self, mprobs, *params):
        assert len(params) == len(self.predicate_indices), self.parameter_order
        R = self._instantaneous_mask_f.copy()
        for (indices, par) in zip(self.predicate_indices, params):
            R[indices] *= par
        return R

    def ascii_art(self, delim="", delim2="|", max_width=70, return_table=False):
        """An ASCII-art table representing the model.  'delim' delimits
        parameter names, 'delim2' delimits motifs"""
        from cogent3.util.table import Table

        labels = [m for m in self.alphabet]
        pars = self.get_matrix_params()
        rows = []
        for i, row in enumerate(pars):
            r = [labels[i]] + [delim.join(cell) for cell in row]
            r[i + 1] = "*"  # identity
            rows.append(r)

        labels.insert(0, r"From\To")
        if self.name:
            title = f"{self.name} rate matrix"
        else:
            title = "rate matrix"

        t = Table(
            header=labels,
            data=rows,
            max_width=max_width,
            title=title,
            index_name=r"From\To",
        )
        result = t if return_table else t.to_string(center=True)
        return result

    def get_matrix_params(self):
        """Return the parameter assignment matrix."""
        dim = len(self.alphabet)
        Pars = numpy.zeros([dim, dim], object)
        for x, y in [(x, y) for x in range(dim) for y in range(dim)]:
            Pars[x][y] = []  # a limitation of numpy.  [x,y] = [] fails!
            if not self._instantaneous_mask[x, y]:
                continue
            for par in self.predicate_masks:
                if self.predicate_masks[par][x, y]:
                    Pars[x, y].append(par)
            # sort the matrix entry to facilitate scaling calculations
            Pars[x, y].sort()
        return Pars

    def get_param_list(self):
        """Return a list of parameter names."""
        return list(self.predicate_masks.keys())

    def is_instantaneous(self, x, y):
        return self._is_instantaneous(x, y)

    def get_substitution_rate_value_from_Q(self, Q, motif_probs, pred):
        pred_mask = list(self._adapt_predicates([pred])[0].values())[0]
        pred_row_totals = numpy.sum(pred_mask * Q, axis=1)
        inst_row_totals = numpy.sum(self._instantaneous_mask * Q, axis=1)
        r = sum(pred_row_totals * motif_probs)
        t = sum(inst_row_totals * motif_probs)
        pred_size = numpy.sum(pred_mask.flat)
        inst_size = sum(self._instantaneous_mask.flat)
        return (r / pred_size) / ((t - r) / (inst_size - pred_size))

    def get_scaled_lengths_from_Q(self, Q, motif_probs, length):
        lengths = {}
        for rule in self.scale_masks:
            lengths[rule] = length * self.get_scale_from_Qs(
                [Q], [1.0], motif_probs, rule
            )
        return lengths

    def get_scale_from_Qs(self, Qs, bin_probs, motif_probss, rule):
        rule = self.get_predicate_mask(rule)
        weighted_scale = 0.0
        bin_probs = numpy.asarray(bin_probs)
        for (Q, bin_prob, motif_probs) in zip(Qs, bin_probs, motif_probss):
            row_totals = numpy.sum(rule * Q, axis=1)
            motif_probs = numpy.asarray(motif_probs)
            word_probs = self.calc_word_probs(motif_probs)
            scale = sum(row_totals * word_probs)
            weighted_scale += bin_prob * scale
        return weighted_scale

    def get_predefined_predicates(self):
        # overridden in subclasses
        return {"indel": predicate.parse("-/?")}

    def get_predefined_predicate(self, name):
        # Called by predicate parsing code
        if self._canned_predicates is None:
            self._canned_predicates = self.get_predefined_predicates()
        return self._canned_predicates[name].interpret(self)

    def _adapt_predicates(self, rules):
        # dict or list of callables, predicate objects or predicate strings
        if isinstance(rules, dict):
            rules = list(rules.items())
        else:
            rules = [(None, rule) for rule in rules]
        predicate_masks = {}
        order = []
        for (key, pred) in rules:
            (label, mask) = self.adapt_predicate(pred, key)
            if label in predicate_masks:
                raise KeyError(f'Duplicate predicate name "{label}"')
            predicate_masks[label] = mask
            order.append(label)
        return predicate_masks, order

    def adapt_predicate(self, pred, label=None):
        if isinstance(pred, str):
            pred = predicate.parse(pred)
        elif isinstance(pred, Callable):
            pred = predicate.UserPredicate(pred)
        pred_func = pred.make_model_predicate(self)
        label = label or repr(pred)
        mask = predicate2matrix(
            self.get_alphabet(), pred_func, mask=self._instantaneous_mask
        )
        return (label, mask)

    def get_predicate_mask(self, pred):
        if pred in self.scale_masks:
            mask = self.scale_masks[pred]
        elif pred in self.predicate_masks:
            mask = self.predicate_masks[pred]
        else:
            (label, mask) = self.adapt_predicate(pred)
        return mask


class Stationary(StationaryQ, Parametric):
    @extend_docstring_from(Parametric.__init__)
    def __init__(self, *args, **kw):
        """ """
        Parametric.__init__(self, *args, **kw)


class TimeReversible(Stationary):
    @extend_docstring_from(Stationary.__init__)
    def __init__(self, *args, **kw):
        """ """
        Stationary.__init__(self, *args, **kw)
        if not self.symmetric:
            raise ValueError(
                "TimeReversible exchangeability terms must be fully balanced"
            )


class _TimeReversibleNucleotide(TimeReversible):
    def get_predefined_predicates(self):
        return {
            "transition": predicate.parse("R/R") | predicate.parse("Y/Y"),
            "transversion": predicate.parse("R/Y"),
            "indel": predicate.parse("-/?"),
            "kappa": (predicate.parse("R/R") | predicate.parse("Y/Y")).aliased("kappa"),
        }


class TimeReversibleNucleotide(_TimeReversibleNucleotide):
    """A nucleotide substitution model."""

    def __init__(self, *args, **kw):
        _TimeReversibleNucleotide.__init__(self, moltype.DNA.alphabet, *args, **kw)


class TimeReversibleDinucleotide(_TimeReversibleNucleotide):
    """A dinucleotide substitution model."""

    def __init__(self, *args, **kw):
        _TimeReversibleNucleotide.__init__(
            self, moltype.DNA.alphabet, motif_length=2, *args, **kw
        )


class TimeReversibleTrinucleotide(_TimeReversibleNucleotide):
    """A trinucleotide substitution model."""

    def __init__(self, *args, **kw):
        _TimeReversibleNucleotide.__init__(
            self, moltype.DNA.alphabet, motif_length=3, *args, **kw
        )


class TimeReversibleProtein(TimeReversible):
    """base protein substitution model."""

    def __init__(self, with_selenocysteine=False, *args, **kw):
        alph = moltype.PROTEIN.alphabet
        if not with_selenocysteine:
            alph = alph.get_subset("U", excluded=True)
        TimeReversible.__init__(self, alph, *args, **kw)


def EmpiricalProteinMatrix(
    matrix, motif_probs=None, optimise_motif_probs=False, recode_gaps=True, **kw
):
    alph = moltype.PROTEIN.alphabet.get_subset("U", excluded=True)
    return Empirical(
        alph,
        rate_matrix=matrix,
        motif_probs=motif_probs,
        model_gaps=False,
        recode_gaps=recode_gaps,
        optimise_motif_probs=optimise_motif_probs,
        **kw,
    )


class _CodonPredicates:
    """predicates for silent and replacement substitutions"""

    def __init__(self, gc):
        """
        Parameters
        ----------

        gc
            a genetic code instance
        """
        self.gc = gc

    def silent(self, x, y):
        return x != "---" and y != "---" and self.gc[x] == self.gc[y]

    def replacement(self, x, y):
        return x != "---" and y != "---" and self.gc[x] != self.gc[y]


class _Codon:
    long_indels_are_instantaneous = True

    def _is_instantaneous(self, x, y):
        if x == self.gapmotif or y == self.gapmotif:
            return x != y
        else:
            ndiffs = sum([X != Y for (X, Y) in zip(x, y)])
            return ndiffs == 1

    def get_predefined_predicates(self):
        codon_preds = _CodonPredicates(self.get_alphabet().get_genetic_code())

        preds = _TimeReversibleNucleotide.get_predefined_predicates(self)
        preds.update(
            {
                "indel": predicate.parse("???/---"),
                "silent": predicate.UserPredicate(codon_preds.silent),
                "replacement": predicate.UserPredicate(codon_preds.replacement),
                "omega": predicate.UserPredicate(codon_preds.replacement),
            }
        )
        return preds


class TimeReversibleCodon(_Codon, _TimeReversibleNucleotide):
    """Core substitution model for codons"""

    @extend_docstring_from(_TimeReversibleNucleotide.__init__)
    def __init__(self, alphabet=None, gc=None, **kw):
        if gc is not None:
            alphabet = moltype.CodonAlphabet(gc=gc)
        alphabet = alphabet or moltype.STANDARD_CODON
        _TimeReversibleNucleotide.__init__(self, alphabet, **kw)
