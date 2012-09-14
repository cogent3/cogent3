#!/usr/bin/env python
"""Performs singular-value decomposition on a set of Q-matrices."""

from __future__ import division
from cogent.maths.stats.test import std # numpy.std is biased
from cogent.maths.matrix_exponentiation import FastExponentiator as expm
from cogent.maths.matrix_logarithm import logm
#note: corrcoef and cov assume rows are observations, cols are variables
from numpy import log, newaxis as NewAxis, array, zeros, product, sqrt, ravel,\
                  sum, sort, reshape, corrcoef, cov, mean
from numpy.random import random
from numpy.linalg import svd, eigvals

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

array_type = type(array([0.1]))

def var(x):
    return std(x)**2

def ratio_two_best(eigenvalues):
    """Returns ratio of best to second best eigenvalue (from vector)."""
    try:
        sorted = sort(eigenvalues)
        return sorted[-1]/sorted[-2]
    except TypeError:   #probably complex-valued
        eigs = abs(eigenvalues)
        sorted = sort(eigs)
        return sorted[-1]/sorted[-2]

def ratio_best_to_sum(eigenvalues):
    """Returns ratio of best singular value to sum. Expects a vector.

    Corresponds to fraction of variance explained by best singular value."""
    try:
        sorted = sort(eigenvalues)
        return sorted[-1]/sum(eigenvalues, axis=0)
    except TypeError:   #probably complex-valued
        eigs = abs(eigenvalues)
        sorted = sort(eigs)
        return sorted[-1]/sum(eigenvalues, axis=0)

def euclidean_distance(q1, q2):
    """Returns Euclidean distance between arrays q1 and q2."""
    diff = ravel(q1 - q2)
    return sqrt(sum(diff*diff, axis=0))
    
def euclidean_norm(m):
    """Returns Euclidean norm of an array or matrix m."""
    flattened = ravel(m)
    return sqrt(sum(flattened*flattened, axis=0))

def _dists_from_mean_slow(qs):
    """Returns distance of each item in qs from the mean.
    
    WARNING: Slow method used only for compatibility testing. Do not use.
    """
    n = len(qs)
    average = mean(qs, axis=0)
    result = zeros(n, 'float64')
    for i in range(n):
        result[i] = euclidean_distance(average, qs[i])
    return result

def dists_from_v(a, v=None):
    """Returns vector of distances between each row in a from v.

    If v is None, returns distance between each row and the mean.
    """
    if v is None:
        v = mean(a, axis=0)
    diff = a - v
    return(sqrt(sum(diff*diff, axis=1)))
        
def weiss(eigens):
    """Returns Weiss(20003) statistic, sum(ln(1+i)) for i in vector of eigens."""
    return sum(log(1+eigens), axis=0)

def three_item_combos(items):
    """Iterates over the 3-item sets from items.

    Doesn't check that items are unique.
    """
    total = len(items)
    for i in range(total-2):
        curr_i = items[i]
        for j in range(i+1, total-1):
            curr_j = items[j]
            for k in range(j+1, total):
                yield curr_i, curr_j, items[k]

def two_item_combos(items):
    """Iterates over the 2-item sets from items.

    Doesn't check that items are unique.
    """
    total = len(items)
    for i in range(total-1):
        curr_i = items[i]
        for j in range(i+1, total):
            yield curr_i, items[j]

def pca_qs(flat_qs):
    """Returns Principal Components vector from correlations in flat_qs."""
    return eigvals(corrcoef(flat_qs))

def pca_cov_qs(flat_qs):
    """Returns Principal Components vector from covariance in flat_qs."""
    return eigvals(cov(flat_qs))

def svd_qs(flat_qs):
    """Returns singular vals from flat_qs directly (returns v, ignores u,w)."""
    return svd(flat_qs)[1]


