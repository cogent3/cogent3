#!/usr/bin/env python
"""microarray_normalize.py: provides functions to normalize array data.

Will use for testing normalization methods against simulated data.

Implementation notes:

All normalization methods should use the interface f(a) -> n, where a is an 
array in which the rows are genes (or probes) and the columns are samples 
(i.e. chips); n is an array of the same shape that contains the normalized 
values.

All probe consolidation methods should be f(probes, groups) -> genes where 
probes is the array of probes (rows = probes, cols = samples), and genes
is the array of genes.

All platform consolidation methods should be f([p_1, p_2, ...]) -> g where
the input is a list of probe arrays for each platform, and the output is a
single gene array.

All missing value imputation methods should be f(a, m) -> n where a is an array 
where each row is a gene (or probe), m is a dict such that m[(x,y)] is True if
a[x,y] is missing (we are assuming that missing values are rare), and n is a
dense array containing imputed values at missing positions.

Revision History
10/28/05 Rob Knight: file created.
11/10/05 Micah Hamady: merged w/my implementations
"""
from cogent.maths.stats.distribution import ndtri
from numpy import ceil, arange, argsort, sort, array, log2, zeros, ravel, \
                  transpose, take, mean, std
from numpy.linalg import svd

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Micah Hamady", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def zscores(a):
    """Converts a to zscores in each col, i.e. subtract mean and div by stdev."""
    return (a - mean(a,axis=0)) / std(a, axis=0)

def logzscores(a):
    """Takes log (base 2) of values then computes zscores"""
    return zscores(log2(a))

def ranks(a):
    """Converts a to absolute ranks in each col, 0 = smallest value.
    
    Doesn't break ties: instead, assigns arbitrary ranks.
    """
    return(argsort(argsort(a,0),0))

def quantiles(a):
    """Converts a to quantiles p(x<a) in each col, 0=smallest, (n-1)/n=largest.

    Doesn't break ties.
    """
    return ranks(a)/float(len(a))

def make_quantile_normalizer(dist):
    """Returns f(a) that converts to the quantile value in each col.

    dist should be an array with bins equally spaced from 0 to 1, giving
    the value in each bin (i.e. cumulative prob of f(x) at f(i/len(dist))
    should be stored in dist[i]) -- can generate from distribution or generate
    empirically.
    """
    def qn(a):
        result = (quantiles(a)*len(dist)).astype('i')
        return take(dist, result)
    return qn

def make_normal_quantile_normalizer(mean, sd, bins=1000):
    """returns f(a) that converts a to the specified normal distribution."""
    dist = array([ndtri(i)*sd for i in arange(1.0/bins,1,1.0/bins)])
    dist = (dist * sd) + mean
    return make_quantile_normalizer(dist)

def make_empirical_quantile_normalizer(data):
    """returns f(a) that converts a to the same distribution as data."""
    return make_quantile_normalizer(sort(data))

def geometric_mean(vals):
    """
    Compute geometric mean

    vals: array of values to compute geometric mean of
    """
    if len(vals) < 2:
        raise ValueError, "Not enought values for geometric mean."
    gmean = 1.0 
    root = 1.0 / len(vals)

    return reduce(lambda x, y: x * y, pow(vals, root))

#FIXME: no formal unit tests for remaining functions.

def svd_cols(a):
    """Returns svd by cols. Note: usually want to discard col[0] of result.
    
    WARNING: will fail on redundant rows! Need to filter out redundant rows
    beforehand, then put back in after normalization. Picks the matrix that
    has the same shape as a.
    """
    u, v, w = svd(a)
    if u.shape == a.shape:
        return u
    elif w.shape == a.shape:
        return w
    else:
        raise TypeError, "Neither u nor w same shape as a."

def rank_rows_by_variance(a):
    """Returns array containing indices 0..num_rows, ranked by variance."""
    return argsort(std(a,1))

def n_least_quantile_var(a, n):
    """Finds the best n rows in a with minimum quantile variance."""
    return rank_row_by_variance(quantiles(a))[:n]

def make_constant_rank_normalizer(constant_rank_indices=None):
    """Returns normalizer that divs by constant_rank_indices.
    
    If constant_rank_indices is an int rather than a list, picks best n.
    """
    def constant_rank_normalizer(a):
        if isinstance(int, constant_rank_indices):
            constant_rank_indices = rank_rows_by_variance(a,\
                constant_rank_indices)
        h = take(a, constant_rank_indices)
        return a / mean(h)
    return constant_rank_normalizer

def omit_rows(a, f):
    """Returns copy of a where f(a[i]) == True, and map of [new]->[old] index.
    
    All omit_rows functions will return both these values.
    """
    mask = array(map(f, a))
    coords = nonzero(mask)
    return take(a, coords), coords

def omit_rows_below_mean_threshold(a, t):
    """Returns copy of a without rows that have mean below threshold."""
    return omit_rows(a, lambda f: mean(f) >= t)

def omit_rows_below_max_threshold(a, t):
    """Returns copy of a without rows that have max val below threshold."""
    return omit_rows(a, lambda f: max(f) >= t)

def group_rows(a, groups):
    """Converts a into list of groups, specified as lists of lists of indices.

    i.e. groups should be a list of groups, where each group contains the 
    indices of the rows that belong to it.
    """
    return [take(a, indices) for indices in groups]

def max_per_group(a, groups):
    """Returns array with max item for each col for each group."""
    return array(map(max, group_rows(a, groups)))

def mean_per_group(a, groups):
    """Returns array with mean for each col for each group."""
    return array(map(mean, group_rows(a, groups)))

def housekeeping_gene_normalize(a, housekeeping_gene_indexes):
    """
    Normalize matrix based on mean/std dev of housekeeping genes 

    Need to refactor to use more effiecient array operations (in place)
    
    a: microarray data. expects rows to be genes, columns to be samples 
    housekeeping_gene_indexes: list of indexes of genes (rows) that represent
        the housekeeping set
    """
    norm_values = ravel(take(a, housekeeping_gene_indexes, 1))
    hk_mean = mean(norm_values)
    hk_std_dev = std(norm_values)
 
    if not hk_std_dev:
        raise NormalizationError, "Cannot normalize, std dev is zero." 
    
    return (a - hk_mean) / hk_std_dev


