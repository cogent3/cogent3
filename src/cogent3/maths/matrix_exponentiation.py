# 4 implementations of P = exp(Q*t)
# APIs along the lines of:
#   exponentiator = WhateverExponenentiator(Q or Q derivative(s))
#   P = exponentiator(t)
#
#           Class(Q)     instance(t)       Limitations
# Eigen      slow           fast           not too asymm
# SemiSym    slow           fast           mprobs > 0
# Pade       instant        slow
# Taylor     instant        very slow

import warnings

import numpy

from numpy.linalg import eig, inv, solve


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Zongzhi Liu"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class _Exponentiator:
    def __repr__(self):
        return f"{self.__class__.__name__}({repr(self.Q)})"


class EigenExponentiator(_Exponentiator):
    """A matrix ready for fast exponentiation.  P=exp(Q*t)"""

    __slots__ = ["Q", "ev", "roots", "evI", "evT"]

    def __init__(self, Q, roots, ev, evT, evI):
        self.Q = Q
        self.evI = evI
        self.evT = evT
        self.ev = ev
        self.roots = roots

    def __call__(self, t):
        exp_roots = numpy.exp(t * self.roots)
        result = numpy.inner(self.evT * exp_roots, self.evI)
        if result.dtype.kind == "c":
            result = numpy.asarray(result.real)
        result = numpy.maximum(result, 0.0)
        return result


def SemiSymmetricExponentiator(motif_probs, Q):
    """Like EigenExponentiator, but more numerically stable and
    30% faster when the rate matrix (Q/motif_probs) is symmetrical.
    Only usable when all motif probs > 0.  Unlike the others
    it needs to know the motif probabilities."""

    H = numpy.sqrt(motif_probs)
    H2 = numpy.divide.outer(H, H)
    # A = Q * H2
    # assert numpy.allclose(A, numpy.transpose(A)), A
    (roots, R) = eig(Q * H2)
    ev = R.T / H2
    evI = (R * H2).T
    # self.evT = numpy.transpose(self.ev)
    return EigenExponentiator(Q, roots, ev, ev.T, evI)


# These next two are slow exponentiators, they don't get any speed up
# from reusing Q with a new t, but for compatability with the diagonalising
# approach they look like they do.  They are derived from code in SciPy.


class TaylorExponentiator(_Exponentiator):
    def __init__(self, Q):
        self.Q = Q
        self.q = 21

    def __call__(self, t=1.0):
        """Compute the matrix exponential using a Taylor series of order q."""
        A = self.Q * t
        M = A.shape[0]
        eA = numpy.identity(M, float)
        trm = eA
        for k in range(1, self.q):
            trm = numpy.dot(trm, A / float(k))
            eA += trm
        while not numpy.allclose(eA, eA - trm):
            k += 1
            trm = numpy.dot(trm, A / float(k))
            eA += trm
        if k >= self.q:
            warnings.warn(f"Taylor series lengthened from {self.q} to {k + 1}")
            self.q = k + 1
        return eA


class PadeExponentiator(_Exponentiator):
    def __init__(self, Q):
        self.Q = Q

    def __call__(self, t=1.0):
        """Compute the matrix exponential using Pade approximation of order q."""
        A = self.Q * t
        M = A.shape[0]
        # Scale A so that norm is < 1/2
        norm = numpy.maximum.reduce(numpy.sum(numpy.absolute(A), axis=1))
        j = int(numpy.floor(numpy.log(max(norm, 0.5)) / numpy.log(2.0))) + 1
        A = A / 2.0 ** j

        # How many iterations required
        e = 1.0
        q = 0
        qf = 1.0
        while e > 1e-12:
            q += 1
            q2 = 2.0 * q
            qf *= q ** 2 / (q2 * (q2 - 1) * q2 * (q2 + 1))
            e = 8 * (norm / (2 ** j)) ** (2 * q) * qf

        # Pade Approximation for exp(A)
        X = A
        c = 1.0 / 2
        N = numpy.identity(M) + c * A
        D = numpy.identity(M) - c * A
        for k in range(2, q + 1):
            c = c * (q - k + 1) / (k * (2 * q - k + 1))
            X = numpy.dot(A, X)
            cX = c * X
            N = N + cX
            if not k % 2:
                D = D + cX
            else:
                D = D - cX
        F = solve(D, N)
        for k in range(1, j + 1):
            F = numpy.dot(F, F)
        return F


def FastExponentiator(Q):
    roots, evT = eig(Q)
    ev = evT.T
    return EigenExponentiator(Q, roots, ev, evT, inv(ev))


def CheckedExponentiator(Q):
    roots, evT = eig(Q)
    ev = evT.T
    evI = inv(ev)
    reQ = numpy.inner(ev.T * roots, evI).real
    if not numpy.allclose(Q, reQ):
        raise ArithmeticError("eigen failed precision test")
    return EigenExponentiator(Q, roots, ev, evT, evI)


def RobustExponentiator(Q):
    return PadeExponentiator(Q)
