#!/usr/bin/env python
"""Alternate matrix log algorithms. A simple implementation of matrix log,
following Brett Easton's suggestion, and a taylor series expansion approach.

WARNING: The methods are not robust!
"""
from itertools import combinations

from numpy import allclose, argmin, array, diag, dot, exp, eye
from numpy import inner as innerproduct
from numpy import isclose, log, ones, pi, zeros
from numpy.linalg import eig as eigenvectors
from numpy.linalg import inv as inverse
from numpy.linalg import norm


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2014, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Von Bing Yap", "Ben Kaehler"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def _is_Q_ok(Q):
    """Tests whether a square matrix is a valid transition rate matrix"""
    n = Q.shape[0]
    if not allclose(Q.imag, 0.0):
        return False
    offd = Q * (1.0 - eye(n))
    if not allclose(offd[offd < 0.0], 0.0):
        return False
    one = ones(n)
    if not allclose(Q.dot(one), 0.0):
        return False
    return True


def is_generator_unique(Q):
    """Conservatively tests whether a transition rate matrix uniquely yields
    its transition probability matrix"""
    if not Q.shape[0] in (3, 4):
        raise NotImplementedError("Only Q of 3x3 or 4x4 supported")
    assert _is_Q_ok(Q), "Q must be a valid transition rate matrix"

    e, V = eigenvectors(Q)
    n = len(e)

    # Assert that the matrix is diagonalisable
    if not allclose(V.dot(diag(e)).dot(inverse(V)), Q):
        raise ArithmeticError("matrix not diagonalisable")

    # Find the Perron-Frobenius eigenvalue
    PF_EV = argmin([norm(ones(n) / n - v / v.sum()) for v in V.T])
    # Don't mess with the P-F eigenvalue - it has a special job to do
    ix = list(range(0, PF_EV)) + list(range(PF_EV + 1, n))

    real_close = []
    expe = exp(e)
    for i, j in combinations(ix, 2):
        if isclose(e.real[i], e.real[j]):
            real_close.append((i, j))

        # Can't deal with non-primary roots yet
        if isclose(expe[i], expe[j]):
            raise NotImplementedError("non-primary root detected:\n" + repr(Q))

    # If the real parts of the eigenvalues are distinct, we're ok
    # For each candidate complex conjugate pair, check for equivalent Qs
    for i, j in real_close:
        s = zeros(n)
        s[i] = 1.0
        s[j] = -1.0
        gen = 2.0 * pi * complex(0.0, 1.0) * V.dot(diag(s)).dot(inverse(V))
        Qtest = Q + gen
        if _is_Q_ok(Qtest):
            return False
        Qtest = Q - gen
        if _is_Q_ok(Qtest):
            return False

    return True


def logm(P):
    """Returns logarithm of a matrix.

    This method should work if the matrix is positive definite and
    diagonalizable.
    """
    roots, ev = eigenvectors(P)
    evI = inverse(ev.T)
    evT = ev
    if not allclose(P, innerproduct(evT * roots, evI)):
        raise ArithmeticError("eigendecomposition failed")

    log_roots = log(roots)
    return innerproduct(evT * log_roots, evI)


def logm_taylor(P, tol=1e-30):
    """returns the matrix log computed using the taylor series. If the Frobenius
    norm of P-I is > 1, raises an exception since the series is not gauranteed
    to converge. The series is continued until the Frobenius norm of the current
    element is < tol.

    Note: This exit condition is theoretically crude but seems to work
    reasonably well.

    Parameters
    ----------
        tol - the tolerance

    """
    P = array(P)
    I = eye(P.shape[0])
    X = P - I
    assert norm(X, ord="fro") < 1, "Frobenius norm > 1"

    Y = I
    Q = zeros(P.shape, dtype="double")
    i = 1
    while norm(Y / i, ord="fro") > tol:
        Y = dot(Y, X)
        Q += (-1) ** (i - 1) * Y / i
        i += 1
    return Q
