#!/usr/bin/env python
"""Simple implementation of matrix log, following Brett Easton's suggestion.

WARNING: This method is not robust!
"""
from numpy import transpose, log, inner as innerproduct
from numpy.linalg import inv as inverse, eig as eigenvectors

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
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
