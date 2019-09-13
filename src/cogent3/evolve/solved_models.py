"""P matrices for some DNA models can be calculated without going via the
intermediate rate matrix Q.  A Cython implementation of this calculation can
be used when Q is not required, for example during likelihood tree optimisation.
Equivalent pure python code is NOT provided because it is typically slower
than the rate-matrix based alternative and provides no extra functionality.
"""

import numpy

from numpy.testing import assert_allclose

from cogent3.evolve.predicate import MotifChange
from cogent3.evolve.substitution_model import (
    CalcDefn,
    TimeReversibleNucleotide,
)
from cogent3.maths.matrix_exponentiation import FastExponentiator
from cogent3.util.modules import ExpectedImportError, importVersionedModule


try:
    from . import _solved_models

    # _solved_models = importVersionedModule('_solved_models', globals(),
    # (1, 0), "only matrix exponentiating DNA models")
except ImportError:
    _solved_models = None

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.9.13a"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


class PredefinedNucleotide(TimeReversibleNucleotide):
    _default_expm_setting = None

    # Instead of providing calc_exchangeability_matrix this subclass overrrides
    # make_continuous_psub_defn to bypass the Q / Qd step.

    def make_continuous_psub_defn(
        self, word_probs, mprobs_matrix, distance, rate_params
    ):
        # Only one set of mprobs will be used
        assert word_probs is mprobs_matrix
        # Order of bases is assumed later, so check it really is Y,Y,R,R:
        alphabet = self.get_alphabet()
        assert set(list(alphabet)[:2]) == set(["T", "C"])
        assert set(list(alphabet)[2:]) == set(["G", "A"])
        # Should produce the same P as an ordinary Q based model would:
        self.check_psub_calculations_match()
        return CalcDefn(self.calc_psub_matrix, name="psubs")(
            word_probs, distance, *rate_params
        )

    def calc_psub_matrix(self, pi, time, kappa_y=1.0, kappa_r=None):
        """Is F81, HKY83 or TN93 when passed 0, 1 or 2 parameters"""
        if kappa_r is None:
            kappa_r = kappa_y
        result = numpy.empty([4, 4], float)
        _solved_models.calc_TN93_P(pi, time, kappa_y, kappa_r, result)
        return result

    def check_psub_calculations_match(self):
        pi = numpy.array([0.1, 0.2, 0.3, 0.4])
        params = [4, 6][: len(self.parameter_order)]
        Q = self.calcQ(pi, pi, *params)
        P1 = FastExponentiator(Q)(0.5)
        P2 = self.calc_psub_matrix(pi, 0.5, *params)
        assert_allclose(P1, P2)


def _solvedNucleotide(name, predicates, rate_matrix_required=True, **kw):
    if _solved_models is not None and not rate_matrix_required:
        klass = PredefinedNucleotide
    else:
        klass = TimeReversibleNucleotide
    kw["model_gaps"] = False
    return klass(name=name, predicates=predicates, **kw)


kappa_y = MotifChange("T", "C").aliased("kappa_y")
kappa_r = MotifChange("A", "G").aliased("kappa_r")
kappa = (kappa_y | kappa_r).aliased("kappa")


def TN93(**kw):
    """Tamura and Nei 1993 model"""
    kw["recode_gaps"] = True
    return _solvedNucleotide("TN93", [kappa_y, kappa_r], **kw)


def HKY85(**kw):
    """Hasegawa, Kishino and Yanamo 1985 model"""
    kw["recode_gaps"] = True
    return _solvedNucleotide("HKY85", [kappa], **kw)


def F81(**kw):
    """Felsenstein's 1981 model"""
    kw["recode_gaps"] = True
    return _solvedNucleotide("F81", [], **kw)
