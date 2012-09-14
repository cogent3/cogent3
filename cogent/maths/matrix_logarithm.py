#!/usr/bin/env python
"""Alternate matrix log algorithms. A simple implementation of matrix log,
following Brett Easton's suggestion, and a taylor series expansion approach.

WARNING: The methods are not robust!
"""
from numpy import array, dot, eye, zeros, transpose, log, inner as innerproduct
from numpy.linalg import inv as inverse, eig as eigenvectors, norm

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Von Bing Yap"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def logm(P):
    """Returns logarithm of a matrix.

    This method should work if the matrix is positive definite and
    diagonalizable.
    """
    roots, ev = eigenvectors(P)
    evI = inverse(ev.T)
    evT = ev
    log_roots = log(roots)
    return innerproduct(evT * log_roots, evI)

def logm_taylor(P, tol=1e-30):
    """returns the matrix log computed using the taylor series. If the Frobenius
    norm of P-I is > 1, raises an exception since the series is not gauranteed
    to converge. The series is continued until the Frobenius norm of the current
    element is < tol.
    
    Note: This exit condition is theoretically crude but seems to work
    reasonably well.
    
    Arguments:
        tol - the tolerance
    """
    P = array(P)
    I = eye(P.shape[0])
    X = P - I
    assert norm(X, ord='fro') < 1, "Frobenius norm > 1"
    
    Y = I
    Q = zeros(P.shape, dtype="double")
    i = 1
    while norm(Y/i, ord='fro') > tol:
        Y = dot(Y,X)
        Q += ((-1)**(i-1)*Y/i)
        i += 1
    return Q
