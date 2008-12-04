#!/usr/bin/env python
"""transform.py - transform a rectangular frequency matrix.
Ref:
    Legendre, P. and E. Gallagher. 2001.  Ecologically meaningful
    transformations for ordination of species data.  Oecologia: 129: 271-280.

Note:
    The distance here means row by row distances, each row represents the
    ordinate vector of a point.
"""
from numpy import sum, square, sqrt, asmatrix, zeros

__author__ = "Zongzhi Liu"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Zongzhi Liu"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Zongzhi Liu"
__email__ = "zongzhi.liu@gmail.com"
__status__ = "Production"

def sumsq(x, axis=None):
    return sum(square(x), axis)

def norm(x, axis=None):
    return sqrt(sum(square(x), axis))

def trans_chord(m):
    """transform m so that the euclidean dist of m' equals chord dist of m.
    """
    m = asmatrix(m)
    row_norms = norm(m, axis=1)
    result = m / row_norms
    return result

def trans_chisq(m):
    """transform m so that the euclidean dist of m' equals chisq dist of m.
    """
    m = asmatrix(m)
    grand_sum, row_sums, col_sums = m.sum(), m.sum(1), m.sum(0)
    result = m * sqrt(grand_sum)
    result /= row_sums
    result /= sqrt(col_sums)
    return result

def trans_specprof(m):
    """transform m so that the euclid dist of m' equals species profile of m.
    """
    m = asmatrix(m)
    row_sums = sum(m, axis=1)
    result = m / row_sums
    return result

def trans_hellinger(m):
    """transform m so that the euclid dist of m' equals hellinger dist of m.
    """
    m = asmatrix(m)
    row_sums = sum(m, axis=1)
    result = sqrt(m / row_sums)
    return result

def dist_euclidean(m):
    """return a row-row euclidean dist matrix from m"""
    nrow, ncol = m.shape
    result = zeros((nrow, nrow), float)
    for r in range(nrow):
        for c in range(r):
            curr = norm(m[r] - m[c])
            result[r,c] = result[c,r] = curr 
    return result

def dist_chord(m):
    """return a row-row chord dist matrix from m"""
    return dist_euclidean(trans_chord(m))

def dist_chisq(m):
    """return a row-row chisq dist matrix from m"""
    return dist_euclidean(trans_chisq(m))

def dist_hellinger(m):
    """return a row-row hellinger dist matrix from m"""
    return dist_euclidean(trans_hellinger(m))

def dist_specprof(m):
    """return a row-row species profile from m"""
    return dist_euclidean(trans_specprof(m))
