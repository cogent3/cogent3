import numpy

from cogent3.core import moltype
from cogent3.evolve.discrete_markov import PsubMatrixDefn
from cogent3.evolve.predicate import MotifChange
from cogent3.maths.optimisers import ParameterOutOfBoundsError

from .substitution_model import (
    Parametric,
    Stationary,
    TimeReversibleNucleotide,
    _Codon,
    _ContinuousSubstitutionModel,
    _SubstitutionModel,
)


__author__ = "Peter Maxwell, Gavin Huttley and Andrew Butterfield"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__contributors__ = ["Gavin Huttley", "Peter Maxwell", "Ben Kaeheler", "Ananias Iliadis"]
__license__ = "BSD-3"
__version__ = "2019.8.23a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def _gen_sym_preds():
    pair = {"A": "T", "T": "A", "G": "C", "C": "G"}
    sym_preds = []
    for f, t in "AG", "AT", "CG", "CT", "GT":
        sym_preds.append(
            MotifChange(f, t, forward_only=True)
            | MotifChange(pair[f], pair[t], forward_only=True)
        )
    return sym_preds


_sym_preds = _gen_sym_preds()


class General(Parametric):
    """A continuous substitution model with one free parameter for each and
    every possible instantaneous substitution."""

    # k = self.param_pick[i,j], 0<=k<=N+1
    # k==0:   not instantaneous, should be 0.0 in Q
    # k<=N:   apply Kth exchangeability parameter
    # k==N+1: no parameter, should be 1.0 in unscaled Q

    def __init__(self, alphabet, **kw):
        Parametric.__init__(self, alphabet, **kw)

        alphabet = self.get_alphabet()  # as may be altered by recode_gaps etc.
        mask = self._instantaneous_mask
        N = len(alphabet)
        self.param_pick = numpy.zeros([N, N], int)
        self.parameter_order = []
        for (i, x) in enumerate(alphabet):
            for j in numpy.flatnonzero(mask[i]):
                y = alphabet[j]
                self.parameter_order.append("%s/%s" % (x, y))
                self.param_pick[i, j] = len(self.parameter_order)
        const_param = self.parameter_order.pop()
        self.symmetric = False
        self.check_params_exist()

    def calc_exchangeability_matrix(self, mprobs, *params):
        return numpy.array((0.0,) + params + (1.0,)).take(self.param_pick)


class GeneralStationary(Stationary):
    """A continuous substitution model with one free parameter for each and
    every possible instantaneous substitution, except the last in each column.
    As general as can be while still having stationary motif probabilities"""

    # @extend_docstring_from(_ContinuousSubstitutionModel)
    def __init__(self, alphabet, **kw):
        Stationary.__init__(self, alphabet, **kw)

        alphabet = self.get_alphabet()  # as may be altered by recode_gaps etc.
        mask = self._instantaneous_mask
        N = len(alphabet)
        self.param_pick = numpy.zeros([N, N], int)
        self.parameter_order = []
        self.last_in_column = []
        for (d, (row, col)) in enumerate(zip(mask, mask.T)):
            row = list(numpy.flatnonzero(row[d:]) + d)
            col = list(numpy.flatnonzero(col[d:]) + d)
            if col:
                self.last_in_column.append((col.pop(), d))
            else:
                assert not row
            inst = [(d, j) for j in row] + [(i, d) for i in col]

            for (i, j) in inst:
                (x, y) = [alphabet[k] for k in [i, j]]
                self.parameter_order.append("%s/%s" % (x, y))
                self.param_pick[i, j] = len(self.parameter_order)
        const_param = self.parameter_order.pop()
        self.symmetric = False
        self.check_params_exist()

    def calc_exchangeability_matrix(self, mprobs, *params):
        R = numpy.array((0.0,) + params + (1.0,)).take(self.param_pick)
        for (i, j) in self.last_in_column:
            assert i > j
            row_total = numpy.dot(mprobs, R[j])
            col_total = numpy.dot(mprobs, R[:, j])
            required = row_total - col_total
            if required < 0.0:
                raise ParameterOutOfBoundsError
            R[i, j] = required / mprobs[i]
        return R


class DiscreteSubstitutionModel(_SubstitutionModel):
    _default_expm_setting = None

    def _is_instantaneous(self, x, y):
        return True

    def get_param_list(self):
        return []

    def make_rate_params(self, bprobs):
        return []

    def make_psubs_defn(self, bprobs, word_probs, mprobs_matrix, rate_params):
        assert len(rate_params) == 0
        assert word_probs is mprobs_matrix, "Must use simple mprob model"
        motifs = tuple(self.get_alphabet())
        return PsubMatrixDefn(
            name="psubs",
            dimension=("motif", motifs),
            default=None,
            dimensions=("locus", "edge"),
        )


class NonReversibleNucleotide(Parametric):
    """A nucleotide substitution model."""

    def __init__(self, *args, **kw):
        Parametric.__init__(self, moltype.DNA.alphabet, *args, **kw)


class NonReversibleDinucleotide(Parametric):
    """A dinucleotide substitution model."""

    def __init__(self, *args, **kw):
        Parametric.__init__(self, moltype.DNA.alphabet, motif_length=2, *args, **kw)


class NonReversibleTrinucleotide(Parametric):
    """A trinucleotide substitution model."""

    def __init__(self, *args, **kw):
        Parametric.__init__(self, moltype.DNA.alphabet, motif_length=3, *args, **kw)


class NonReversibleCodon(_Codon, Parametric):
    def __init__(self, alphabet=None, gc=None, **kw):
        if gc is not None:
            alphabet = moltype.CodonAlphabet(gc=gc)
        alphabet = alphabet or moltype.STANDARD_CODON
        Parametric.__init__(self, alphabet, **kw)


class StrandSymmetric(NonReversibleNucleotide):
    def __init__(self, **kw):
        for argname in ("predicates", "recode_gaps", "model_gaps"):
            kw.pop(argname, None)
        super(StrandSymmetric, self).__init__(
            predicates=_sym_preds, recode_gaps=True, model_gaps=False, **kw
        )


class NonReversibleProtein(Parametric):
    """base protein substitution model."""

    def __init__(self, with_selenocysteine=False, *args, **kw):
        alph = moltype.PROTEIN.alphabet
        if not with_selenocysteine:
            alph = alph.get_subset("U", excluded=True)
        Parametric.__init__(self, alph, *args, **kw)
