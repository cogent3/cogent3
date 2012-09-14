#!/usr/bin/env python
"""Procrustes analysis.  Main fn: procrustes

See for example: 
Principles of Multivariate analysis, by Krzanowski
"""

from numpy.linalg import svd
from numpy import array, sqrt, sum, zeros, trace, dot, transpose,\
    divide, square, subtract, shape, any, abs, mean
from numpy import append as numpy_append

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Production"

def procrustes(data1, data2):
    """Procrustes analysis, a similarity test for two data sets.
    
    Each input matrix is a set of points or vectors (the rows of the matrix)
    The dimension of the space is the number of columns of each matrix.
    Given two identially sized matrices, procrustes standardizes both
    such that:
    - trace(AA') = 1  (A' is the transpose, and the product is
    a standard matrix product).
    - Both sets of points are centered around the origin
    
    Procrustes then applies the optimal transform to the second matrix
    (including scaling/dilation, rotations, and reflections) to minimize
    M^2 = sum(square(mtx1 - mtx2)), or the sum of the squares of the pointwise
    differences between the two input datasets
    
    If two data sets have different dimensionality (different number of
    columns), simply add columns of zeros the the smaller of the two.
    
    This function was not designed to handle datasets with different numbers of
    datapoints (rows)
    
    Arguments:
        - data1: matrix, n rows represent points in k (columns) space
        data1 is the reference data, after it is standardised, the data from
        data2 will be transformed to fit the pattern in data1
        - data2: n rows of data in k space to be fit to data1.  Must be the 
        same shape (numrows, numcols) as data1
        - both must have >1 unique points
        
    Returns:
        - mtx1: a standardized version of data1
        - mtx2: the orientation of data2 that best fits data1.
        centered, but not necessarily trace(mtx2*mtx2') = 1
        - disparity: a metric for the dissimilarity of the two datasets,
        disparity = M^2 defined above
        
    Notes:
        - The disparity should not depend on the order of the input matrices, 
        but the output matrices will, as only the first output matrix is
        guaranteed to be scaled such that trace(AA') = 1.
        - duplicate datapoints are generally ok, duplicating a data point will
        increase it's effect on the procrustes fit.
        - the disparity scales as the number of points per input matrix
        
    
    """
    SMALL_NUM = 1e-6 # used to check for zero values in added dimension
    
    # make local copies
#     mtx1 = array(data1.copy(),'d')
#     mtx2 = array(data2.copy(),'d')
    num_rows, num_cols = shape(data1)
    if (num_rows, num_cols) != shape(data2):
        raise ValueError("input matrices must be of same shape")
    if (num_rows == 0 or num_cols == 0):
        raise ValueError("input matrices must be >0 rows, >0 cols")
        
    
    # add a dimension to allow reflections (rotations in n + 1 dimensions)
    mtx1 = numpy_append(data1, zeros((num_rows, 1)), 1)
    mtx2 = numpy_append(data2, zeros((num_rows, 1)), 1)
    
    # standardize each matrix
    mtx1 = center(mtx1)
    mtx2 = center(mtx2)
    
    if ((not any(mtx1)) or (not any(mtx2))):
        raise ValueError("input matrices must contain >1 unique points")

    mtx1 = normalize(mtx1)
    mtx2 = normalize(mtx2)
    
       
    # transform mtx2 to minimize disparity (sum( (mtx1[i,j] - mtx2[i,j])^2) )
    mtx2 = match_points(mtx1, mtx2)
    
    # WARNING: I haven't proven that after matching the matrices, no point has
    # a nonzero component in the added dimension.  I believe it is true,
    # though, since the unchanged matrix has no points extending into 
    # that dimension
    
    if any(abs(mtx2[:,-1]) > SMALL_NUM):
        raise StandardError("we have accidentially added a dimension to \
the matrix, and the vectors have nonzero components in that dimension")
    
    # strip extra dimension which was added to allow reflections
    mtx1 = mtx1[:,:-1]
    mtx2 = mtx2[:,:-1]
    
    disparity = get_disparity(mtx1, mtx2)
    
    return mtx1, mtx2, disparity
    
def center(mtx):
    """translate all data (rows of the matrix) to center on the origin
    
    returns a shifted version of the input data.  The new matrix is such that
    the center of mass of the row vectors is centered at the origin.  
    Returns a numpy float ('d') array
    """
    result = array(mtx, 'd')
    result -= mean(result, 0) 
    # subtract each column's mean from each element in that column
    return result

def normalize(mtx):
    """change scaling of data (in rows) such that trace(mtx*mtx') = 1
    
    mtx' denotes the transpose of mtx """
    result = array(mtx, 'd')
    num_pts, num_dims = shape(result)
    mag = trace(dot(result, transpose(result)))
    norm = sqrt(mag)
    result /= norm
    return result

def match_points(mtx1, mtx2):
    """returns a transformed mtx2 that matches mtx1.
    
    returns a new matrix which is a transform of mtx2.  Scales and rotates
    a copy of mtx 2.  See procrustes docs for details.
    """
    u,s,vh = svd(dot(transpose(mtx1), mtx2))
    q = dot(transpose(vh), transpose(u))
    new_mtx2 = dot(mtx2, q)
    new_mtx2 *= sum(s)
    
    return new_mtx2

def get_disparity(mtx1, mtx2):
    """ returns a measure of the dissimilarity between two data sets
    
    returns M^2 = sum(square(mtx1 - mtx2)), the pointwise sum of squared
    differences"""
    return(sum(square(mtx1 - mtx2)))
