#! /usr/bin/env python
"""A collection of pre-defined models.  These are provided for convenience so that
users do not need to keep reconstructing the standard models.  We encourage users
to think about the assumptions in these models and consider if their problem could
benefit from a user defined model.
Note that models that do not traditionally deal with gaps are implemented with
gap recoding that will convert gaps to Ns, and model gaps set to False."""

# this file using functions etc. to allow each model to serve as an example for users
# wishing to construct their own models
from itertools import permutations

# The models are constructed in a strait forward manner with no attempt to condense
import numpy

from cogent3 import DNA
from cogent3.evolve import ns_substitution_model, substitution_model
from cogent3.evolve.predicate import MotifChange, omega
from cogent3.evolve.solved_models import _solved_nucleotide
from cogent3.evolve.substitution_model import _SubstitutionModel
from cogent3.util.table import Table


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Matthew Wakefield", "Peter Maxwell", "Gavin Huttley", "James Kondilios"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"


nucleotide_models = []
codon_models = []
protein_models = []
models = []


_model_types = {
    "nucleotide": nucleotide_models,
    "codon": codon_models,
    "protein": protein_models,
}
_all_models = {}


class register_model:
    """
    decorator for registering functions that construct substitution models.

    The name of the wrapped function becomes the model abbreviation used
    for selecting the model with ``get_models()``

    Parameters
    ----------
    model_type: str
        valid values are 'codon', 'nucleotide', 'protein'
    """

    def __init__(self, model_type: str):
        assert model_type in _model_types, f"{model_type!r} not in {list(_model_types)}"
        self._model_type = model_type

    def __call__(self, func):
        series = _model_types[self._model_type]
        name = func.__name__
        if name in models:
            raise ValueError(f"{name!r} already in models")

        series.append(name)
        models.append(name)
        _all_models[name] = func
        return func


# Substitution model rate matrix predicates

_gtr_preds = [MotifChange(x, y) for x, y in ["AC", "AG", "AT", "CG", "CT"]]
_omega = omega
_kappa_y = MotifChange("T", "C").aliased("kappa_y")
_kappa_r = MotifChange("A", "G").aliased("kappa_r")
_kappa = (_kappa_y | _kappa_r).aliased("kappa")
_cg = MotifChange("CG").aliased("G")
_cg_k = (_cg & _kappa).aliased("G.K")


def _make_gn_preds():
    return [
        MotifChange(f, t, forward_only=True)
        for f, t in permutations("ACTG", 2)
        if f != "T" or t != "G"
    ]


_general_preds = _make_gn_preds()


@register_model("nucleotide")
def BH(optimise_motif_probs=True, **kw):
    """Barry and Hartigan Discrete Time substitution model

    Barry and Hartigan 1987. Biometrics 43: 261–76.
    """
    return DT(
        optimise_motif_probs=optimise_motif_probs, motif_length=1, name="BH", **kw
    )


@register_model("nucleotide")
def DT(optimise_motif_probs=True, motif_length=1, **kw):
    """
    Discrete Time substitution model (non-stationary, non-reversible).
    motif_length=2 makes this a dinucleotide model, motif_length=3 a
    trinucleotide model.
    """
    alpha = DNA.alphabet.get_word_alphabet(motif_length)
    kw["optimise_motif_probs"] = optimise_motif_probs
    kw["mprob_model"] = "tuple"
    kw["name"] = kw.get("name", f"DT-{motif_length}")
    return ns_substitution_model.DiscreteSubstitutionModel(alpha, **kw)


@register_model("nucleotide")
def GN(optimise_motif_probs=True, **kw):
    """General Markov Nucleotide (non-stationary, non-reversible).

    Kaehler, Yap, Zhang, Huttley, 2015, Sys Biol 64 (2): 281–93"""
    required = dict(
        optimise_motif_probs=optimise_motif_probs, name="GN", predicates=_general_preds
    )
    kwargs = dict(recode_gaps=True, model_gaps=False)
    kwargs.update(kw)
    kwargs.update(required)
    return ns_substitution_model.NonReversibleNucleotide(**kwargs)


@register_model("nucleotide")
def ssGN(optimise_motif_probs=True, **kw):
    """strand-symmetric general Markov nucleotide (non-stationary, non-reversible).

    Kaehler, 2017, Journal of Theoretical Biology 420: 144–51"""
    # note the StrandSymmetric class predefines the predicates and name
    return ns_substitution_model.StrandSymmetric(
        optimise_motif_probs=optimise_motif_probs, name="ssGN", **kw
    )


@register_model("nucleotide")
def K80(**kw):
    """Kimura 1980"""
    required = dict(name="K80", equal_motif_probs=True, optimise_motif_probs=False)
    kw["recode_gaps"] = kw.get("recode_gaps", True)
    kwargs = {}
    kwargs.update(kw)
    kwargs.update(required)
    return HKY85(**kwargs)


@register_model("nucleotide")
def JC69(**kw):
    """Jukes and Cantor's 1969 model"""
    required = dict(name="JC69", equal_motif_probs=True, optimise_motif_probs=False)
    kw["recode_gaps"] = kw.get("recode_gaps", True)
    kwargs = {}
    kwargs.update(kw)
    kwargs.update(required)
    return F81(**kwargs)


@register_model("nucleotide")
def GTR(**kw):
    """General Time Reversible nucleotide substitution model."""
    required = dict(
        name="GTR", predicates=_gtr_preds, mprob_model="conditional", model_gaps=False
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleNucleotide(**kwargs)


@register_model("nucleotide")
def TN93(**kw):
    """Tamura and Nei 1993 model"""
    kw["recode_gaps"] = kw.get("recode_gaps", True)
    kw["name"] = "TN93"
    return _solved_nucleotide([_kappa_y, _kappa_r], **kw)


@register_model("nucleotide")
def HKY85(**kw):
    """Hasegawa, Kishino and Yano 1985 model"""
    kw["recode_gaps"] = kw.get("recode_gaps", True)
    # this function called by others, so we don't overwrite name if it exists
    kw["name"] = kw.get("name", "HKY85")
    return _solved_nucleotide([_kappa], **kw)


@register_model("nucleotide")
def F81(**kw):
    """Felsenstein's 1981 model"""
    kw["recode_gaps"] = kw.get("recode_gaps", True)
    # this function called by others, so we don't overwrite name if it exists
    kw["name"] = kw.get("name", "F81")
    return _solved_nucleotide([], **kw)


# Codon Models
@register_model("codon")
def CNFGTR(**kw):
    """Conditional nucleotide frequency codon substitution model, GTR variant
    (with params analagous to the nucleotide GTR model).

    Yap, Lindsay, Easteal and Huttley, 2010, Mol Biol Evol 27: 726-734"""
    required = dict(
        name="CNFGTR",
        predicates=_gtr_preds + [_omega],
        mprob_model="conditional",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def CNFHKY(**kw):
    """Conditional nucleotide frequency codon substitution model, HKY variant
    (with kappa, the ratio of transitions to transversions)

    Yap, Lindsay, Easteal and Huttley, 2010, Mol Biol Evol 27: 726-734"""
    required = dict(
        name="CNFHKY",
        predicates=[_kappa, _omega],
        mprob_model="conditional",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def MG94HKY(**kw):
    """Muse and Gaut 1994 codon substitution model, HKY variant (with kappa,
    the ratio of transitions to transversions)

    Muse and Gaut, 1994, Mol Biol Evol, 11, 715-24"""
    required = dict(
        name="MG94HKY",
        predicates=[_kappa, _omega],
        mprob_model="monomer",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def MG94GTR(**kw):
    """Muse and Gaut 1994 codon substitution model, GTR variant (with params
    analagous to the nucleotide GTR model)

    Muse and Gaut, 1994, Mol Biol Evol, 11, 715-24"""
    required = dict(
        name="MG94GTR",
        predicates=_gtr_preds + [_omega],
        mprob_model="monomer",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def GY94(**kw):
    """Goldman and Yang 1994 codon substitution model.

    N Goldman and Z Yang, 1994, Mol Biol Evol, 11(5):725-36."""
    required = dict(name="GY94")
    kwargs = {}
    kwargs.update(kw)
    kwargs.update(required)
    return Y98(**kwargs)


@register_model("codon")
def Y98(**kw):
    """Yang's 1998 substitution model, a derivative of the GY94.

    Z Yang, 1998, Mol Biol Evol, 15(5):568-73"""

    required = dict(
        predicates=[_kappa, _omega],
        mprob_model="tuple",
        model_gaps=False,
        name=kw.get("name", "Y98"),
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def H04G(**kw):
    """Huttley 2004 CpG substitution model. Includes a term for substitutions
    to or from CpG's.

    GA Huttley, 2004, Mol Biol Evol, 21(9):1760-8"""
    required = dict(
        name="H04G",
        predicates=[_cg, _kappa, _omega],
        mprob_model="tuple",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def H04GK(**kw):
    """Huttley 2004 CpG substitution model. Includes a term for transition
    substitutions to or from CpG's.

    GA Huttley, 2004, Mol Biol Evol, 21(9):1760-8"""
    required = dict(
        name="H04GK",
        predicates=[_cg_k, _kappa, _omega],
        mprob_model="tuple",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def H04GGK(**kw):
    """Huttley 2004 CpG substitution model. Includes a general term for
    substitutions to or from CpG's and an adjustment for CpG transitions.

    GA Huttley, 2004, Mol Biol Evol, 21(9):1760-8"""
    required = dict(
        name="H04GGK",
        predicates=[_cg, _cg_k, _kappa, _omega],
        mprob_model="tuple",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return substitution_model.TimeReversibleCodon(**kwargs)


@register_model("codon")
def GNC(optimise_motif_probs=True, **kw):
    """General Nucleotide Codon, a non-reversible codon model.

    Kaehler, Yap, Huttley, 2017, Gen Biol Evol 9(1): 134–49"""
    required = dict(
        name="GNC",
        optimise_motif_probs=optimise_motif_probs,
        predicates=_general_preds + [_omega],
        mprob_model="tuple",
        model_gaps=False,
    )
    kwargs = dict(recode_gaps=True, motif_probs=None)
    kwargs.update(kw)
    kwargs.update(required)
    return ns_substitution_model.NonReversibleCodon(**kwargs)


# Protein Models

# Empirical Protein Models

DSO78_matrix = numpy.array(
    [
        [
            0.00000000e00,
            3.60000000e01,
            1.20000000e02,
            1.98000000e02,
            1.80000000e01,
            2.40000000e02,
            2.30000000e01,
            6.50000000e01,
            2.60000000e01,
            4.10000000e01,
            7.20000000e01,
            9.80000000e01,
            2.50000000e02,
            8.90000000e01,
            2.70000000e01,
            4.09000000e02,
            3.71000000e02,
            2.08000000e02,
            0.00000000e00,
            2.40000000e01,
        ],
        [
            3.60000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.10000000e01,
            2.80000000e01,
            4.40000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.90000000e01,
            0.00000000e00,
            2.30000000e01,
            1.61000000e02,
            1.60000000e01,
            4.90000000e01,
            0.00000000e00,
            9.60000000e01,
        ],
        [
            1.20000000e02,
            0.00000000e00,
            0.00000000e00,
            1.15300000e03,
            0.00000000e00,
            1.25000000e02,
            8.60000000e01,
            2.40000000e01,
            7.10000000e01,
            0.00000000e00,
            0.00000000e00,
            9.05000000e02,
            1.30000000e01,
            1.34000000e02,
            0.00000000e00,
            9.50000000e01,
            6.60000000e01,
            1.80000000e01,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            1.98000000e02,
            0.00000000e00,
            1.15300000e03,
            0.00000000e00,
            0.00000000e00,
            8.10000000e01,
            4.30000000e01,
            6.10000000e01,
            8.30000000e01,
            1.10000000e01,
            3.00000000e01,
            1.48000000e02,
            5.10000000e01,
            7.16000000e02,
            1.00000000e00,
            7.90000000e01,
            3.40000000e01,
            3.70000000e01,
            0.00000000e00,
            2.20000000e01,
        ],
        [
            1.80000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.50000000e01,
            4.80000000e01,
            1.96000000e02,
            0.00000000e00,
            1.57000000e02,
            9.20000000e01,
            1.40000000e01,
            1.10000000e01,
            0.00000000e00,
            1.40000000e01,
            4.60000000e01,
            1.30000000e01,
            1.20000000e01,
            7.60000000e01,
            6.98000000e02,
        ],
        [
            2.40000000e02,
            1.10000000e01,
            1.25000000e02,
            8.10000000e01,
            1.50000000e01,
            0.00000000e00,
            1.00000000e01,
            0.00000000e00,
            2.70000000e01,
            7.00000000e00,
            1.70000000e01,
            1.39000000e02,
            3.40000000e01,
            2.80000000e01,
            9.00000000e00,
            2.34000000e02,
            3.00000000e01,
            5.40000000e01,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            2.30000000e01,
            2.80000000e01,
            8.60000000e01,
            4.30000000e01,
            4.80000000e01,
            1.00000000e01,
            0.00000000e00,
            7.00000000e00,
            2.60000000e01,
            4.40000000e01,
            0.00000000e00,
            5.35000000e02,
            9.40000000e01,
            6.06000000e02,
            2.40000000e02,
            3.50000000e01,
            2.20000000e01,
            4.40000000e01,
            2.70000000e01,
            1.27000000e02,
        ],
        [
            6.50000000e01,
            4.40000000e01,
            2.40000000e01,
            6.10000000e01,
            1.96000000e02,
            0.00000000e00,
            7.00000000e00,
            0.00000000e00,
            4.60000000e01,
            2.57000000e02,
            3.36000000e02,
            7.70000000e01,
            1.20000000e01,
            1.80000000e01,
            6.40000000e01,
            2.40000000e01,
            1.92000000e02,
            8.89000000e02,
            0.00000000e00,
            3.70000000e01,
        ],
        [
            2.60000000e01,
            0.00000000e00,
            7.10000000e01,
            8.30000000e01,
            0.00000000e00,
            2.70000000e01,
            2.60000000e01,
            4.60000000e01,
            0.00000000e00,
            1.80000000e01,
            2.43000000e02,
            3.18000000e02,
            3.30000000e01,
            1.53000000e02,
            4.64000000e02,
            9.60000000e01,
            1.36000000e02,
            1.00000000e01,
            0.00000000e00,
            1.30000000e01,
        ],
        [
            4.10000000e01,
            0.00000000e00,
            0.00000000e00,
            1.10000000e01,
            1.57000000e02,
            7.00000000e00,
            4.40000000e01,
            2.57000000e02,
            1.80000000e01,
            0.00000000e00,
            5.27000000e02,
            3.40000000e01,
            3.20000000e01,
            7.30000000e01,
            1.50000000e01,
            1.70000000e01,
            3.30000000e01,
            1.75000000e02,
            4.60000000e01,
            2.80000000e01,
        ],
        [
            7.20000000e01,
            0.00000000e00,
            0.00000000e00,
            3.00000000e01,
            9.20000000e01,
            1.70000000e01,
            0.00000000e00,
            3.36000000e02,
            2.43000000e02,
            5.27000000e02,
            0.00000000e00,
            1.00000000e00,
            1.70000000e01,
            1.14000000e02,
            9.00000000e01,
            6.20000000e01,
            1.04000000e02,
            2.58000000e02,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            9.80000000e01,
            0.00000000e00,
            9.05000000e02,
            1.48000000e02,
            1.40000000e01,
            1.39000000e02,
            5.35000000e02,
            7.70000000e01,
            3.18000000e02,
            3.40000000e01,
            1.00000000e00,
            0.00000000e00,
            4.20000000e01,
            1.03000000e02,
            3.20000000e01,
            4.95000000e02,
            2.29000000e02,
            1.50000000e01,
            2.30000000e01,
            9.50000000e01,
        ],
        [
            2.50000000e02,
            1.90000000e01,
            1.30000000e01,
            5.10000000e01,
            1.10000000e01,
            3.40000000e01,
            9.40000000e01,
            1.20000000e01,
            3.30000000e01,
            3.20000000e01,
            1.70000000e01,
            4.20000000e01,
            0.00000000e00,
            1.53000000e02,
            1.03000000e02,
            2.45000000e02,
            7.80000000e01,
            4.80000000e01,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            8.90000000e01,
            0.00000000e00,
            1.34000000e02,
            7.16000000e02,
            0.00000000e00,
            2.80000000e01,
            6.06000000e02,
            1.80000000e01,
            1.53000000e02,
            7.30000000e01,
            1.14000000e02,
            1.03000000e02,
            1.53000000e02,
            0.00000000e00,
            2.46000000e02,
            5.60000000e01,
            5.30000000e01,
            3.50000000e01,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            2.70000000e01,
            2.30000000e01,
            0.00000000e00,
            1.00000000e00,
            1.40000000e01,
            9.00000000e00,
            2.40000000e02,
            6.40000000e01,
            4.64000000e02,
            1.50000000e01,
            9.00000000e01,
            3.20000000e01,
            1.03000000e02,
            2.46000000e02,
            0.00000000e00,
            1.54000000e02,
            2.60000000e01,
            2.40000000e01,
            2.01000000e02,
            8.00000000e00,
        ],
        [
            4.09000000e02,
            1.61000000e02,
            9.50000000e01,
            7.90000000e01,
            4.60000000e01,
            2.34000000e02,
            3.50000000e01,
            2.40000000e01,
            9.60000000e01,
            1.70000000e01,
            6.20000000e01,
            4.95000000e02,
            2.45000000e02,
            5.60000000e01,
            1.54000000e02,
            0.00000000e00,
            5.50000000e02,
            3.00000000e01,
            7.50000000e01,
            3.40000000e01,
        ],
        [
            3.71000000e02,
            1.60000000e01,
            6.60000000e01,
            3.40000000e01,
            1.30000000e01,
            3.00000000e01,
            2.20000000e01,
            1.92000000e02,
            1.36000000e02,
            3.30000000e01,
            1.04000000e02,
            2.29000000e02,
            7.80000000e01,
            5.30000000e01,
            2.60000000e01,
            5.50000000e02,
            0.00000000e00,
            1.57000000e02,
            0.00000000e00,
            4.20000000e01,
        ],
        [
            2.08000000e02,
            4.90000000e01,
            1.80000000e01,
            3.70000000e01,
            1.20000000e01,
            5.40000000e01,
            4.40000000e01,
            8.89000000e02,
            1.00000000e01,
            1.75000000e02,
            2.58000000e02,
            1.50000000e01,
            4.80000000e01,
            3.50000000e01,
            2.40000000e01,
            3.00000000e01,
            1.57000000e02,
            0.00000000e00,
            0.00000000e00,
            2.80000000e01,
        ],
        [
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            7.60000000e01,
            0.00000000e00,
            2.70000000e01,
            0.00000000e00,
            0.00000000e00,
            4.60000000e01,
            0.00000000e00,
            2.30000000e01,
            0.00000000e00,
            0.00000000e00,
            2.01000000e02,
            7.50000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            6.10000000e01,
        ],
        [
            2.40000000e01,
            9.60000000e01,
            0.00000000e00,
            2.20000000e01,
            6.98000000e02,
            0.00000000e00,
            1.27000000e02,
            3.70000000e01,
            1.30000000e01,
            2.80000000e01,
            0.00000000e00,
            9.50000000e01,
            0.00000000e00,
            0.00000000e00,
            8.00000000e00,
            3.40000000e01,
            4.20000000e01,
            2.80000000e01,
            6.10000000e01,
            0.00000000e00,
        ],
    ]
)

DSO78_freqs = {
    "A": 0.087126912873087131,
    "C": 0.033473966526033475,
    "E": 0.04952995047004953,
    "D": 0.046871953128046873,
    "G": 0.088611911388088618,
    "F": 0.039771960228039777,
    "I": 0.036885963114036892,
    "H": 0.033617966382033626,
    "K": 0.08048191951808048,
    "M": 0.014752985247014754,
    "L": 0.085356914643085369,
    "N": 0.040431959568040438,
    "Q": 0.038254961745038257,
    "P": 0.050679949320050689,
    "S": 0.069576930423069588,
    "R": 0.040903959096040908,
    "T": 0.058541941458058543,
    "W": 0.010493989506010494,
    "V": 0.064717935282064723,
    "Y": 0.029915970084029919,
}

JTT92_matrix = numpy.array(
    [
        [
            0.0,
            56.0,
            81.0,
            105.0,
            15.0,
            179.0,
            27.0,
            36.0,
            35.0,
            30.0,
            54.0,
            54.0,
            194.0,
            57.0,
            58.0,
            378.0,
            475.0,
            298.0,
            9.0,
            11.0,
        ],
        [
            56.0,
            0.0,
            10.0,
            5.0,
            78.0,
            59.0,
            69.0,
            17.0,
            7.0,
            23.0,
            31.0,
            34.0,
            14.0,
            9.0,
            113.0,
            223.0,
            42.0,
            62.0,
            115.0,
            209.0,
        ],
        [
            81.0,
            10.0,
            0.0,
            767.0,
            4.0,
            130.0,
            112.0,
            11.0,
            26.0,
            7.0,
            15.0,
            528.0,
            15.0,
            49.0,
            16.0,
            59.0,
            38.0,
            31.0,
            4.0,
            46.0,
        ],
        [
            105.0,
            5.0,
            767.0,
            0.0,
            5.0,
            119.0,
            26.0,
            12.0,
            181.0,
            9.0,
            18.0,
            58.0,
            18.0,
            323.0,
            29.0,
            30.0,
            32.0,
            45.0,
            10.0,
            7.0,
        ],
        [
            15.0,
            78.0,
            4.0,
            5.0,
            0.0,
            5.0,
            40.0,
            89.0,
            4.0,
            248.0,
            43.0,
            10.0,
            17.0,
            4.0,
            5.0,
            92.0,
            12.0,
            62.0,
            53.0,
            536.0,
        ],
        [
            179.0,
            59.0,
            130.0,
            119.0,
            5.0,
            0.0,
            23.0,
            6.0,
            27.0,
            6.0,
            14.0,
            81.0,
            24.0,
            26.0,
            137.0,
            201.0,
            33.0,
            47.0,
            55.0,
            8.0,
        ],
        [
            27.0,
            69.0,
            112.0,
            26.0,
            40.0,
            23.0,
            0.0,
            16.0,
            45.0,
            56.0,
            33.0,
            391.0,
            115.0,
            597.0,
            328.0,
            73.0,
            46.0,
            11.0,
            8.0,
            573.0,
        ],
        [
            36.0,
            17.0,
            11.0,
            12.0,
            89.0,
            6.0,
            16.0,
            0.0,
            21.0,
            229.0,
            479.0,
            47.0,
            10.0,
            9.0,
            22.0,
            40.0,
            245.0,
            961.0,
            9.0,
            32.0,
        ],
        [
            35.0,
            7.0,
            26.0,
            181.0,
            4.0,
            27.0,
            45.0,
            21.0,
            0.0,
            14.0,
            65.0,
            263.0,
            21.0,
            292.0,
            646.0,
            47.0,
            103.0,
            14.0,
            10.0,
            8.0,
        ],
        [
            30.0,
            23.0,
            7.0,
            9.0,
            248.0,
            6.0,
            56.0,
            229.0,
            14.0,
            0.0,
            388.0,
            12.0,
            102.0,
            72.0,
            38.0,
            59.0,
            25.0,
            180.0,
            52.0,
            24.0,
        ],
        [
            54.0,
            31.0,
            15.0,
            18.0,
            43.0,
            14.0,
            33.0,
            479.0,
            65.0,
            388.0,
            0.0,
            30.0,
            16.0,
            43.0,
            44.0,
            29.0,
            226.0,
            323.0,
            24.0,
            18.0,
        ],
        [
            54.0,
            34.0,
            528.0,
            58.0,
            10.0,
            81.0,
            391.0,
            47.0,
            263.0,
            12.0,
            30.0,
            0.0,
            15.0,
            86.0,
            45.0,
            503.0,
            232.0,
            16.0,
            8.0,
            70.0,
        ],
        [
            194.0,
            14.0,
            15.0,
            18.0,
            17.0,
            24.0,
            115.0,
            10.0,
            21.0,
            102.0,
            16.0,
            15.0,
            0.0,
            164.0,
            74.0,
            285.0,
            118.0,
            23.0,
            6.0,
            10.0,
        ],
        [
            57.0,
            9.0,
            49.0,
            323.0,
            4.0,
            26.0,
            597.0,
            9.0,
            292.0,
            72.0,
            43.0,
            86.0,
            164.0,
            0.0,
            310.0,
            53.0,
            51.0,
            20.0,
            18.0,
            24.0,
        ],
        [
            58.0,
            113.0,
            16.0,
            29.0,
            5.0,
            137.0,
            328.0,
            22.0,
            646.0,
            38.0,
            44.0,
            45.0,
            74.0,
            310.0,
            0.0,
            101.0,
            64.0,
            17.0,
            126.0,
            20.0,
        ],
        [
            378.0,
            223.0,
            59.0,
            30.0,
            92.0,
            201.0,
            73.0,
            40.0,
            47.0,
            59.0,
            29.0,
            503.0,
            285.0,
            53.0,
            101.0,
            0.0,
            477.0,
            38.0,
            35.0,
            63.0,
        ],
        [
            475.0,
            42.0,
            38.0,
            32.0,
            12.0,
            33.0,
            46.0,
            245.0,
            103.0,
            25.0,
            226.0,
            232.0,
            118.0,
            51.0,
            64.0,
            477.0,
            0.0,
            112.0,
            12.0,
            21.0,
        ],
        [
            298.0,
            62.0,
            31.0,
            45.0,
            62.0,
            47.0,
            11.0,
            961.0,
            14.0,
            180.0,
            323.0,
            16.0,
            23.0,
            20.0,
            17.0,
            38.0,
            112.0,
            0.0,
            25.0,
            16.0,
        ],
        [
            9.0,
            115.0,
            4.0,
            10.0,
            53.0,
            55.0,
            8.0,
            9.0,
            10.0,
            52.0,
            24.0,
            8.0,
            6.0,
            18.0,
            126.0,
            35.0,
            12.0,
            25.0,
            0.0,
            71.0,
        ],
        [
            11.0,
            209.0,
            46.0,
            7.0,
            536.0,
            8.0,
            573.0,
            32.0,
            8.0,
            24.0,
            18.0,
            70.0,
            10.0,
            24.0,
            20.0,
            63.0,
            21.0,
            16.0,
            71.0,
            0.0,
        ],
    ]
)

JTT92_freqs = {
    "A": 0.076747923252076758,
    "C": 0.019802980197019805,
    "E": 0.061829938170061841,
    "D": 0.05154394845605155,
    "G": 0.073151926848073159,
    "F": 0.040125959874040135,
    "I": 0.053760946239053767,
    "H": 0.022943977056022944,
    "K": 0.058675941324058678,
    "M": 0.023825976174023829,
    "L": 0.091903908096091905,
    "N": 0.042644957355042652,
    "Q": 0.040751959248040752,
    "P": 0.050900949099050907,
    "S": 0.068764931235068771,
    "R": 0.051690948309051694,
    "T": 0.058564941435058568,
    "W": 0.014260985739014262,
    "V": 0.066004933995066004,
    "Y": 0.032101967898032102,
}

AH96_matrix = numpy.array(
    [
        [
            0.0,
            59.93,
            17.67,
            9.77,
            6.37,
            120.71,
            13.9,
            96.49,
            8.36,
            25.46,
            141.88,
            26.95,
            54.31,
            1.9,
            23.18,
            387.86,
            480.72,
            195.06,
            1.9,
            6.48,
        ],
        [
            59.93,
            0.0,
            1.9,
            1.9,
            70.8,
            30.71,
            141.49,
            62.73,
            1.9,
            25.65,
            6.18,
            58.94,
            31.26,
            75.24,
            103.33,
            277.05,
            179.97,
            1.9,
            33.6,
            254.77,
        ],
        [
            17.67,
            1.9,
            0.0,
            583.55,
            4.98,
            56.77,
            113.99,
            4.34,
            2.31,
            1.9,
            1.9,
            794.38,
            13.43,
            55.28,
            1.9,
            69.02,
            28.01,
            1.9,
            19.86,
            21.21,
        ],
        [
            9.77,
            1.9,
            583.55,
            0.0,
            2.67,
            28.28,
            49.12,
            3.31,
            313.86,
            1.9,
            1.9,
            63.05,
            12.83,
            313.56,
            1.9,
            54.71,
            14.82,
            21.14,
            1.9,
            13.12,
        ],
        [
            6.37,
            70.8,
            4.98,
            2.67,
            0.0,
            1.9,
            48.16,
            84.67,
            6.44,
            216.06,
            90.82,
            15.2,
            17.31,
            19.11,
            4.69,
            64.29,
            33.85,
            6.35,
            7.84,
            465.58,
        ],
        [
            120.71,
            30.71,
            56.77,
            28.28,
            1.9,
            0.0,
            1.9,
            5.98,
            22.73,
            2.41,
            1.9,
            53.3,
            1.9,
            6.75,
            23.03,
            125.93,
            11.17,
            2.53,
            10.92,
            3.21,
        ],
        [
            13.9,
            141.49,
            113.99,
            49.12,
            48.16,
            1.9,
            0.0,
            12.26,
            127.67,
            11.49,
            11.97,
            496.13,
            60.97,
            582.4,
            165.23,
            77.46,
            44.78,
            1.9,
            7.08,
            670.14,
        ],
        [
            96.49,
            62.73,
            4.34,
            3.31,
            84.67,
            5.98,
            12.26,
            0.0,
            19.57,
            329.09,
            517.98,
            27.1,
            20.63,
            8.34,
            1.9,
            47.7,
            368.43,
            1222.94,
            1.9,
            25.01,
        ],
        [
            8.36,
            1.9,
            2.31,
            313.86,
            6.44,
            22.73,
            127.67,
            19.57,
            0.0,
            14.88,
            91.37,
            608.7,
            50.1,
            465.58,
            141.4,
            105.79,
            136.33,
            1.9,
            24.0,
            51.17,
        ],
        [
            25.46,
            25.65,
            1.9,
            1.9,
            216.06,
            2.41,
            11.49,
            329.09,
            14.88,
            0.0,
            537.53,
            15.16,
            40.1,
            39.7,
            15.58,
            73.61,
            126.4,
            91.67,
            32.44,
            44.15,
        ],
        [
            141.88,
            6.18,
            1.9,
            1.9,
            90.82,
            1.9,
            11.97,
            517.98,
            91.37,
            537.53,
            0.0,
            65.41,
            18.84,
            47.37,
            1.9,
            111.16,
            528.17,
            387.54,
            21.71,
            39.96,
        ],
        [
            26.95,
            58.94,
            794.38,
            63.05,
            15.2,
            53.3,
            496.13,
            27.1,
            608.7,
            15.16,
            65.41,
            0.0,
            73.31,
            173.56,
            13.24,
            494.39,
            238.46,
            1.9,
            10.68,
            191.36,
        ],
        [
            54.31,
            31.26,
            13.43,
            12.83,
            17.31,
            1.9,
            60.97,
            20.63,
            50.1,
            40.1,
            18.84,
            73.31,
            0.0,
            137.29,
            23.64,
            169.9,
            128.22,
            8.23,
            4.21,
            16.21,
        ],
        [
            1.9,
            75.24,
            55.28,
            313.56,
            19.11,
            6.75,
            582.4,
            8.34,
            465.58,
            39.7,
            47.37,
            173.56,
            137.29,
            0.0,
            220.99,
            54.11,
            94.93,
            19.0,
            1.9,
            38.82,
        ],
        [
            23.18,
            103.33,
            1.9,
            1.9,
            4.69,
            23.03,
            165.23,
            1.9,
            141.4,
            15.58,
            1.9,
            13.24,
            23.64,
            220.99,
            0.0,
            6.04,
            2.08,
            7.64,
            21.95,
            1.9,
        ],
        [
            387.86,
            277.05,
            69.02,
            54.71,
            64.29,
            125.93,
            77.46,
            47.7,
            105.79,
            73.61,
            111.16,
            494.39,
            169.9,
            54.11,
            6.04,
            0.0,
            597.21,
            1.9,
            38.58,
            64.92,
        ],
        [
            480.72,
            179.97,
            28.01,
            14.82,
            33.85,
            11.17,
            44.78,
            368.43,
            136.33,
            126.4,
            528.17,
            238.46,
            128.22,
            94.93,
            2.08,
            597.21,
            0.0,
            204.54,
            9.99,
            38.73,
        ],
        [
            195.06,
            1.9,
            1.9,
            21.14,
            6.35,
            2.53,
            1.9,
            1222.94,
            1.9,
            91.67,
            387.54,
            1.9,
            8.23,
            19.0,
            7.64,
            1.9,
            204.54,
            0.0,
            5.37,
            1.9,
        ],
        [
            1.9,
            33.6,
            19.86,
            1.9,
            7.84,
            10.92,
            7.08,
            1.9,
            24.0,
            32.44,
            21.71,
            10.68,
            4.21,
            1.9,
            21.95,
            38.58,
            9.99,
            5.37,
            0.0,
            26.25,
        ],
        [
            6.48,
            254.77,
            21.21,
            13.12,
            465.58,
            3.21,
            670.14,
            25.01,
            51.17,
            44.15,
            39.96,
            191.36,
            16.21,
            38.82,
            1.9,
            64.92,
            38.73,
            1.9,
            26.25,
            0.0,
        ],
    ]
)

AH96_freqs = {
    "A": 0.071999999999999995,
    "C": 0.0060000000000000001,
    "E": 0.024,
    "D": 0.019,
    "G": 0.056000000000000001,
    "F": 0.060999999999999999,
    "I": 0.087999999999999995,
    "H": 0.028000000000000001,
    "K": 0.023,
    "M": 0.053999999999999999,
    "L": 0.16900000000000001,
    "N": 0.039,
    "Q": 0.025000000000000001,
    "P": 0.053999999999999999,
    "S": 0.071999999999999995,
    "R": 0.019,
    "T": 0.085999999999999993,
    "W": 0.029000000000000001,
    "V": 0.042999999999999997,
    "Y": 0.033000000000000002,
}

AH96_mtmammals_matrix = numpy.array(
    [
        [
            0.00000000e00,
            0.00000000e00,
            1.10000000e01,
            0.00000000e00,
            0.00000000e00,
            7.80000000e01,
            8.00000000e00,
            7.50000000e01,
            0.00000000e00,
            2.10000000e01,
            7.60000000e01,
            2.00000000e00,
            5.30000000e01,
            0.00000000e00,
            3.20000000e01,
            3.42000000e02,
            6.81000000e02,
            3.98000000e02,
            5.00000000e00,
            0.00000000e00,
        ],
        [
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            7.00000000e00,
            0.00000000e00,
            3.05000000e02,
            4.10000000e01,
            0.00000000e00,
            2.70000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.86000000e02,
            3.47000000e02,
            1.14000000e02,
            0.00000000e00,
            6.50000000e01,
            5.30000000e02,
        ],
        [
            1.10000000e01,
            0.00000000e00,
            0.00000000e00,
            5.69000000e02,
            5.00000000e00,
            7.90000000e01,
            1.10000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            8.64000000e02,
            2.00000000e00,
            4.90000000e01,
            0.00000000e00,
            1.60000000e01,
            0.00000000e00,
            1.00000000e01,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            0.00000000e00,
            0.00000000e00,
            5.69000000e02,
            0.00000000e00,
            0.00000000e00,
            2.20000000e01,
            2.20000000e01,
            0.00000000e00,
            2.15000000e02,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            2.74000000e02,
            0.00000000e00,
            2.10000000e01,
            4.00000000e00,
            2.00000000e01,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            0.00000000e00,
            7.00000000e00,
            5.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            5.70000000e01,
            0.00000000e00,
            2.46000000e02,
            1.10000000e01,
            6.00000000e00,
            1.70000000e01,
            0.00000000e00,
            0.00000000e00,
            9.00000000e01,
            8.00000000e00,
            6.00000000e00,
            0.00000000e00,
            6.82000000e02,
        ],
        [
            7.80000000e01,
            0.00000000e00,
            7.90000000e01,
            2.20000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            4.70000000e01,
            0.00000000e00,
            0.00000000e00,
            1.80000000e01,
            1.12000000e02,
            0.00000000e00,
            5.00000000e00,
            0.00000000e00,
            1.00000000e00,
        ],
        [
            8.00000000e00,
            3.05000000e02,
            1.10000000e01,
            2.20000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            2.60000000e01,
            0.00000000e00,
            4.58000000e02,
            5.30000000e01,
            5.50000000e02,
            2.32000000e02,
            2.00000000e01,
            1.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.52500000e03,
        ],
        [
            7.50000000e01,
            4.10000000e01,
            0.00000000e00,
            0.00000000e00,
            5.70000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            6.00000000e00,
            2.32000000e02,
            3.78000000e02,
            1.90000000e01,
            5.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            3.60000000e02,
            2.22000000e03,
            0.00000000e00,
            1.60000000e01,
        ],
        [
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            2.15000000e02,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            6.00000000e00,
            0.00000000e00,
            4.00000000e00,
            5.90000000e01,
            4.08000000e02,
            1.80000000e01,
            2.42000000e02,
            5.00000000e01,
            6.50000000e01,
            5.00000000e01,
            0.00000000e00,
            0.00000000e00,
            6.70000000e01,
        ],
        [
            2.10000000e01,
            2.70000000e01,
            0.00000000e00,
            0.00000000e00,
            2.46000000e02,
            0.00000000e00,
            2.60000000e01,
            2.32000000e02,
            4.00000000e00,
            0.00000000e00,
            6.09000000e02,
            0.00000000e00,
            4.30000000e01,
            2.00000000e01,
            6.00000000e00,
            7.40000000e01,
            3.40000000e01,
            1.00000000e02,
            1.20000000e01,
            2.50000000e01,
        ],
        [
            7.60000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.10000000e01,
            0.00000000e00,
            0.00000000e00,
            3.78000000e02,
            5.90000000e01,
            6.09000000e02,
            0.00000000e00,
            2.10000000e01,
            0.00000000e00,
            2.20000000e01,
            0.00000000e00,
            4.70000000e01,
            6.91000000e02,
            8.32000000e02,
            1.30000000e01,
            0.00000000e00,
        ],
        [
            2.00000000e00,
            0.00000000e00,
            8.64000000e02,
            0.00000000e00,
            6.00000000e00,
            4.70000000e01,
            4.58000000e02,
            1.90000000e01,
            4.08000000e02,
            0.00000000e00,
            2.10000000e01,
            0.00000000e00,
            3.30000000e01,
            8.00000000e00,
            4.00000000e00,
            4.46000000e02,
            1.10000000e02,
            0.00000000e00,
            6.00000000e00,
            1.56000000e02,
        ],
        [
            5.30000000e01,
            0.00000000e00,
            2.00000000e00,
            0.00000000e00,
            1.70000000e01,
            0.00000000e00,
            5.30000000e01,
            5.00000000e00,
            1.80000000e01,
            4.30000000e01,
            0.00000000e00,
            3.30000000e01,
            0.00000000e00,
            5.10000000e01,
            9.00000000e00,
            2.02000000e02,
            7.80000000e01,
            0.00000000e00,
            7.00000000e00,
            8.00000000e00,
        ],
        [
            0.00000000e00,
            0.00000000e00,
            4.90000000e01,
            2.74000000e02,
            0.00000000e00,
            0.00000000e00,
            5.50000000e02,
            0.00000000e00,
            2.42000000e02,
            2.00000000e01,
            2.20000000e01,
            8.00000000e00,
            5.10000000e01,
            0.00000000e00,
            2.46000000e02,
            3.00000000e01,
            0.00000000e00,
            3.30000000e01,
            0.00000000e00,
            5.40000000e01,
        ],
        [
            3.20000000e01,
            1.86000000e02,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.80000000e01,
            2.32000000e02,
            0.00000000e00,
            5.00000000e01,
            6.00000000e00,
            0.00000000e00,
            4.00000000e00,
            9.00000000e00,
            2.46000000e02,
            0.00000000e00,
            3.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.60000000e01,
            0.00000000e00,
        ],
        [
            3.42000000e02,
            3.47000000e02,
            1.60000000e01,
            2.10000000e01,
            9.00000000e01,
            1.12000000e02,
            2.00000000e01,
            0.00000000e00,
            6.50000000e01,
            7.40000000e01,
            4.70000000e01,
            4.46000000e02,
            2.02000000e02,
            3.00000000e01,
            3.00000000e00,
            0.00000000e00,
            6.14000000e02,
            0.00000000e00,
            1.70000000e01,
            1.07000000e02,
        ],
        [
            6.81000000e02,
            1.14000000e02,
            0.00000000e00,
            4.00000000e00,
            8.00000000e00,
            0.00000000e00,
            1.00000000e00,
            3.60000000e02,
            5.00000000e01,
            3.40000000e01,
            6.91000000e02,
            1.10000000e02,
            7.80000000e01,
            0.00000000e00,
            0.00000000e00,
            6.14000000e02,
            0.00000000e00,
            2.37000000e02,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            3.98000000e02,
            0.00000000e00,
            1.00000000e01,
            2.00000000e01,
            6.00000000e00,
            5.00000000e00,
            0.00000000e00,
            2.22000000e03,
            0.00000000e00,
            1.00000000e02,
            8.32000000e02,
            0.00000000e00,
            0.00000000e00,
            3.30000000e01,
            0.00000000e00,
            0.00000000e00,
            2.37000000e02,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            5.00000000e00,
            6.50000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.20000000e01,
            1.30000000e01,
            6.00000000e00,
            7.00000000e00,
            0.00000000e00,
            1.60000000e01,
            1.70000000e01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            1.40000000e01,
        ],
        [
            0.00000000e00,
            5.30000000e02,
            0.00000000e00,
            0.00000000e00,
            6.82000000e02,
            1.00000000e00,
            1.52500000e03,
            1.60000000e01,
            6.70000000e01,
            2.50000000e01,
            0.00000000e00,
            1.56000000e02,
            8.00000000e00,
            5.40000000e01,
            0.00000000e00,
            1.07000000e02,
            0.00000000e00,
            0.00000000e00,
            1.40000000e01,
            0.00000000e00,
        ],
    ]
)

AH96_mtmammals_freqs = {
    "A": 0.069199999999999998,
    "C": 0.0064999999999999997,
    "E": 0.023599999999999999,
    "D": 0.018599999999999998,
    "G": 0.0557,
    "F": 0.061100000000000002,
    "I": 0.090499999999999997,
    "H": 0.027699999999999999,
    "K": 0.022100000000000002,
    "M": 0.056099999999999997,
    "L": 0.16750000000000001,
    "N": 0.040000000000000001,
    "Q": 0.023800000000000002,
    "P": 0.053600000000000002,
    "S": 0.072499999999999995,
    "R": 0.0184,
    "T": 0.086999999999999994,
    "W": 0.0293,
    "V": 0.042799999999999998,
    "Y": 0.034000000000000002,
}

WG01_matrix = numpy.array(
    [
        [
            0.0,
            1.02704,
            0.738998,
            1.58285,
            0.210494,
            1.41672,
            0.316954,
            0.193335,
            0.906265,
            0.397915,
            0.893496,
            0.509848,
            1.43855,
            0.908598,
            0.551571,
            3.37079,
            2.12111,
            2.00601,
            0.113133,
            0.240735,
        ],
        [
            1.02704,
            0.0,
            0.0302949,
            0.021352,
            0.39802,
            0.306674,
            0.248972,
            0.170135,
            0.0740339,
            0.384287,
            0.390482,
            0.265256,
            0.109404,
            0.0988179,
            0.528191,
            1.40766,
            0.512984,
            1.00214,
            0.71707,
            0.543833,
        ],
        [
            0.738998,
            0.0302949,
            0.0,
            6.17416,
            0.0467304,
            0.865584,
            0.930676,
            0.039437,
            0.479855,
            0.0848047,
            0.103754,
            5.42942,
            0.423984,
            0.616783,
            0.147304,
            1.07176,
            0.374866,
            0.152335,
            0.129767,
            0.325711,
        ],
        [
            1.58285,
            0.021352,
            6.17416,
            0.0,
            0.0811339,
            0.567717,
            0.570025,
            0.127395,
            2.58443,
            0.154263,
            0.315124,
            0.947198,
            0.682355,
            5.46947,
            0.439157,
            0.704939,
            0.822765,
            0.588731,
            0.156557,
            0.196303,
        ],
        [
            0.210494,
            0.39802,
            0.0467304,
            0.0811339,
            0.0,
            0.049931,
            0.679371,
            1.05947,
            0.088836,
            2.11517,
            1.19063,
            0.0961621,
            0.161444,
            0.0999208,
            0.102711,
            0.545931,
            0.171903,
            0.649892,
            1.52964,
            6.45428,
        ],
        [
            1.41672,
            0.306674,
            0.865584,
            0.567717,
            0.049931,
            0.0,
            0.24941,
            0.0304501,
            0.373558,
            0.0613037,
            0.1741,
            1.12556,
            0.24357,
            0.330052,
            0.584665,
            1.34182,
            0.225833,
            0.187247,
            0.336983,
            0.103604,
        ],
        [
            0.316954,
            0.248972,
            0.930676,
            0.570025,
            0.679371,
            0.24941,
            0.0,
            0.13819,
            0.890432,
            0.499462,
            0.404141,
            3.95629,
            0.696198,
            4.29411,
            2.13715,
            0.740169,
            0.473307,
            0.118358,
            0.262569,
            3.87344,
        ],
        [
            0.193335,
            0.170135,
            0.039437,
            0.127395,
            1.05947,
            0.0304501,
            0.13819,
            0.0,
            0.323832,
            3.17097,
            4.25746,
            0.554236,
            0.0999288,
            0.113917,
            0.186979,
            0.31944,
            1.45816,
            7.8213,
            0.212483,
            0.42017,
        ],
        [
            0.906265,
            0.0740339,
            0.479855,
            2.58443,
            0.088836,
            0.373558,
            0.890432,
            0.323832,
            0.0,
            0.257555,
            0.934276,
            3.01201,
            0.556896,
            3.8949,
            5.35142,
            0.96713,
            1.38698,
            0.305434,
            0.137505,
            0.133264,
        ],
        [
            0.397915,
            0.384287,
            0.0848047,
            0.154263,
            2.11517,
            0.0613037,
            0.499462,
            3.17097,
            0.257555,
            0.0,
            4.85402,
            0.131528,
            0.415844,
            0.869489,
            0.497671,
            0.344739,
            0.326622,
            1.80034,
            0.665309,
            0.398618,
        ],
        [
            0.893496,
            0.390482,
            0.103754,
            0.315124,
            1.19063,
            0.1741,
            0.404141,
            4.25746,
            0.934276,
            4.85402,
            0.0,
            0.198221,
            0.171329,
            1.54526,
            0.683162,
            0.493905,
            1.51612,
            2.05845,
            0.515706,
            0.428437,
        ],
        [
            0.509848,
            0.265256,
            5.42942,
            0.947198,
            0.0961621,
            1.12556,
            3.95629,
            0.554236,
            3.01201,
            0.131528,
            0.198221,
            0.0,
            0.195081,
            1.54364,
            0.635346,
            3.97423,
            2.03006,
            0.196246,
            0.0719167,
            1.086,
        ],
        [
            1.43855,
            0.109404,
            0.423984,
            0.682355,
            0.161444,
            0.24357,
            0.696198,
            0.0999288,
            0.556896,
            0.415844,
            0.171329,
            0.195081,
            0.0,
            0.933372,
            0.679489,
            1.61328,
            0.795384,
            0.314887,
            0.139405,
            0.216046,
        ],
        [
            0.908598,
            0.0988179,
            0.616783,
            5.46947,
            0.0999208,
            0.330052,
            4.29411,
            0.113917,
            3.8949,
            0.869489,
            1.54526,
            1.54364,
            0.933372,
            0.0,
            3.0355,
            1.02887,
            0.857928,
            0.301281,
            0.215737,
            0.22771,
        ],
        [
            0.551571,
            0.528191,
            0.147304,
            0.439157,
            0.102711,
            0.584665,
            2.13715,
            0.186979,
            5.35142,
            0.497671,
            0.683162,
            0.635346,
            0.679489,
            3.0355,
            0.0,
            1.22419,
            0.554413,
            0.251849,
            1.16392,
            0.381533,
        ],
        [
            3.37079,
            1.40766,
            1.07176,
            0.704939,
            0.545931,
            1.34182,
            0.740169,
            0.31944,
            0.96713,
            0.344739,
            0.493905,
            3.97423,
            1.61328,
            1.02887,
            1.22419,
            0.0,
            4.37802,
            0.232739,
            0.523742,
            0.786993,
        ],
        [
            2.12111,
            0.512984,
            0.374866,
            0.822765,
            0.171903,
            0.225833,
            0.473307,
            1.45816,
            1.38698,
            0.326622,
            1.51612,
            2.03006,
            0.795384,
            0.857928,
            0.554413,
            4.37802,
            0.0,
            1.38823,
            0.110864,
            0.291148,
        ],
        [
            2.00601,
            1.00214,
            0.152335,
            0.588731,
            0.649892,
            0.187247,
            0.118358,
            7.8213,
            0.305434,
            1.80034,
            2.05845,
            0.196246,
            0.314887,
            0.301281,
            0.251849,
            0.232739,
            1.38823,
            0.0,
            0.365369,
            0.31473,
        ],
        [
            0.113133,
            0.71707,
            0.129767,
            0.156557,
            1.52964,
            0.336983,
            0.262569,
            0.212483,
            0.137505,
            0.665309,
            0.515706,
            0.0719167,
            0.139405,
            0.215737,
            1.16392,
            0.523742,
            0.110864,
            0.365369,
            0.0,
            2.48539,
        ],
        [
            0.240735,
            0.543833,
            0.325711,
            0.196303,
            6.45428,
            0.103604,
            3.87344,
            0.42017,
            0.133264,
            0.398618,
            0.428437,
            1.086,
            0.216046,
            0.22771,
            0.381533,
            0.786993,
            0.291148,
            0.31473,
            2.48539,
            0.0,
        ],
    ]
)

WG01_freqs = {
    "A": 0.086627908662790867,
    "C": 0.019307801930780195,
    "E": 0.058058905805890577,
    "D": 0.057045105704510574,
    "G": 0.083251808325180837,
    "F": 0.038431903843190382,
    "I": 0.048466004846600491,
    "H": 0.024431302443130246,
    "K": 0.062028606202860624,
    "M": 0.019502701950270197,
    "L": 0.086209008620900862,
    "N": 0.039089403908940397,
    "Q": 0.036728103672810368,
    "P": 0.045763104576310464,
    "S": 0.069517906951790692,
    "R": 0.043972004397200441,
    "T": 0.061012706101270617,
    "W": 0.014385901438590145,
    "V": 0.070895607089560719,
    "Y": 0.035274203527420354,
}


@register_model("protein")
def DSO78(**kw):
    """Dayhoff et al 1978 empirical protein model
    Dayhoff, MO, Schwartz RM, and Orcutt, BC. 1978
    A model of evolutionary change in proteins. Pp. 345-352.
    Atlas of protein sequence and structure, Vol 5, Suppl. 3.
    National Biomedical Research Foundation,  Washington D. C
    Matrix imported from PAML dayhoff.dat file"""
    return substitution_model.EmpiricalProteinMatrix(
        DSO78_matrix, DSO78_freqs, name="DSO78", **kw
    )


@register_model("protein")
def JTT92(**kw):
    """Jones, Taylor and Thornton 1992 empirical protein model
    Jones DT, Taylor WR, Thornton JM.
    The rapid generation of mutation data matrices from protein sequences.
    Comput Appl Biosci. 1992 Jun;8(3):275-82.
    Matrix imported from PAML jones.dat file"""
    return substitution_model.EmpiricalProteinMatrix(
        JTT92_matrix, JTT92_freqs, name="JTT92", **kw
    )


@register_model("protein")
def AH96(**kw):
    """Adachi and Hasegawa 1996 empirical model for mitochondrial proteins.
    Adachi J, Hasegawa M.
    Model of amino acid substitution in proteins encoded by mitochondrial DNA.
    J Mol Evol. 1996 Apr;42(4):459-68.
    Matrix imported from PAML mtREV24.dat file"""
    return substitution_model.EmpiricalProteinMatrix(
        AH96_matrix, AH96_freqs, name="AH96_mtREV24", **kw
    )


def get_model(name, **kw):
    """returns an instance of the named model

    name is case sensitive.

    Parameters
    ----------
    optimise_motif_probs: bool
        Treat like other free parameters.
    recode_gaps: bool
        Whether gaps in an alignment should be treated as an ambiguous state
        instead.

    Notes
    -----
    See available_models() for the full list.
    """
    if isinstance(name, _SubstitutionModel):
        # already a substitution model
        return name
    if name not in models:
        msg = f'Unknown model "{name}". Model names are case sensitive!'
        raise ValueError(msg)

    return _all_models[name](**kw)


def mtREV(**kw):
    return AH96(**kw)


@register_model("protein")
def AH96_mtmammals(**kw):
    """Adachi and Hasegawa 1996 empirical model for mammalian mitochondrial
    proteins.
    Adachi J, Hasegawa M.
    Model of amino acid substitution in proteins encoded by mitochondrial DNA.
    J Mol Evol. 1996 Apr;42(4):459-68.
    Matrix imported from PAML mtmam.dat file"""
    return substitution_model.EmpiricalProteinMatrix(
        AH96_mtmammals_matrix, AH96_mtmammals_freqs, name="AH96_mtmammals", **kw
    )


def mtmam(**kw):
    return AH96_mtmammals(**kw)


@register_model("protein")
def WG01(**kw):
    """Whelan and Goldman 2001 empirical model for globular proteins.
    Whelan S, Goldman N.
    A general empirical model of protein evolution derived from multiple protein
    families using a maximum-likelihood approach.
    Mol Biol Evol. 2001 May;18(5):691-9.
    Matrix imported from PAML wag.dat file"""
    return substitution_model.EmpiricalProteinMatrix(
        WG01_matrix, WG01_freqs, name="WG01", **kw
    )


def available_models(model_types=None):
    """returns Table listing the pre-defined substitution models"""
    column_headings = ["Model Type", "Abbreviation", "Description"]
    if model_types is not None:
        model_types = model_types if not isinstance(model_types, str) else [model_types]
    else:
        model_types = _model_types.keys()

    rows = []
    for mod_type in model_types:
        for abbrev in _model_types[mod_type]:
            if _all_models[abbrev].__doc__:
                description = " ".join(_all_models[abbrev].__doc__.split())
            else:
                description = ""
            rows.append([mod_type, abbrev, description])

    return Table(
        header=column_headings,
        data=rows,
        title="Specify a model using 'Abbreviation' (case sensitive).",
    )
