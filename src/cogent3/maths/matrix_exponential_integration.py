from numpy import (
    allclose,
    array,
    asarray,
    diag,
    dot,
    exp,
    identity,
    inner,
    maximum,
    zeros,
)
from numpy.linalg import LinAlgError, eig, inv

import cogent3.maths.matrix_exponentiation as cme


__author__ = "Ben Kaehler"
__copyright__ = "Copyright 2007-2014, The Cogent Project"
__credits__ = ["Ben Kaehler", "Von Bing Yap", "Gavin Huttley", "Ananias Iliadis"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Ben Kaehler"
__email__ = "benjamin.kaehler@anu.edu.au"
__status__ = "Production"


class _Exponentiator(object):
    def __init__(self, Q):
        self.Q = Q

    def __repr__(self):
        return f"{self.__class__.__name__}({repr(self.Q)})"


class VanLoanIntegratingExponentiator(_Exponentiator):

    """An exponentiator that evaluates int_0^t exp(Q*s)ds * R
    using the method of Van Loan [1]. Complexity is that of the Exponentiator.
    [1] Van Loan, C. F. (1978). Computing integrals involving the matrix
    exponential. IEEE Trans. Autmat. Control 23(3), 395-404."""

    def __init__(self, Q, R=None, exponentiator=cme.RobustExponentiator):
        """
        Q -- an n x n matrix.
        R -- an n x m matrix. Defaults to the identity matrix. Can be a rank-1
        array.
        exponentiator -- Exponentiator used in Van Loan method. Defaults to
        RobustEstimator.
        """
        super(VanLoanIntegratingExponentiator, self).__init__(Q)
        Qdim = len(Q)
        if R is None:
            self.R = identity(Qdim)
        else:
            if len(R.shape) == 1:  # Be kind to rank-1 arrays
                self.R = R.reshape((R.shape[0], 1))
            else:
                self.R = R
        Cdim = Qdim + self.R.shape[1]
        C = zeros((Cdim, Cdim))
        C[:Qdim, :Qdim] = Q
        C[:Qdim, Qdim:] = self.R
        self.expm = exponentiator(C)

    def __call__(self, t=1.0):
        return self.expm(t)[: len(self.Q), len(self.Q) :]


class VonBingIntegratingExponentiator(_Exponentiator):

    """An exponentiator that evaluates int_0^t exp(Q*s)ds
    using the method of Von Bing."""

    def __init__(self, Q):
        """
        Q -- a diagonisable matrix.
        """
        super(VonBingIntegratingExponentiator, self).__init__(Q)
        self.roots, self.evT = eig(Q)
        self.evI = inv(self.evT.T)
        # Remove following check if performance is a concern
        reQ = inner(self.evT * self.roots, self.evI).real
        if not allclose(Q, reQ):
            raise ArithmeticError("eigendecomposition failed")

    def __call__(self, t=1.0):
        int_roots = array(
            [t if abs(x.real) < 1e-6 else (exp(x * t) - 1) / x for x in self.roots]
        )
        result = inner(self.evT * int_roots, self.evI)
        if result.dtype.kind == "c":
            result = asarray(result.real)
        result = maximum(result, 0.0)
        return result


def expected_number_subs(p0, Q, t):
    """returns the expected number of substitutions

    p0
        initial state frequencies
    Q
        continuous time rate matrix, calibrated such that, if it were a stationary process,
        ENS is 1
    t
        ens is returned for Q * t
    """
    try:
        iexpm = VonBingIntegratingExponentiator(Q)
        result = -dot(dot(p0, iexpm(t)), diag(Q))
    except (ArithmeticError, LinAlgError):
        iexpm = VanLoanIntegratingExponentiator(Q, R=diag(Q))
        result = -dot(p0, iexpm(t))[0]
    return result
