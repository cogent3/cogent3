"""P matrices for some DNA models can be calculated without going via the
intermediate rate matrix Q.
"""

import numpy

from numpy.testing import assert_allclose

from cogent3.evolve.substitution_model import (
    CalcDefn,
    TimeReversibleNucleotide,
)
from cogent3.maths.matrix_exponentiation import FastExponentiator

from . import solved_models_numba as _solved_models


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
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


def _solved_nucleotide(predicates, rate_matrix_required=True, **kw):
    if _solved_models is not None and not rate_matrix_required:
        klass = PredefinedNucleotide
    else:
        klass = TimeReversibleNucleotide
    kw["model_gaps"] = False
    return klass(predicates=predicates, **kw)
