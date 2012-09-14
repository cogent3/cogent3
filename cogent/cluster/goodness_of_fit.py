#!usr/bin/env python
"""Goodness of fit of Multidimensional Scaling

Implements several functions that measure the degree of correspondence
between an MDS and its coordinates and the original input distances.

See for example: * Johnson & Wichern (2002): Applied Multivariate
Statistical Analysis

"""


import numpy

__author__ = "Andreas Wilm"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Andreas Wilm"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Andreas Wilm"
__email__ = "andreas.wilm@ucd.ie"
__status__ = "Production"



class Stress(object):
    """Degree of correspondence between input distances and an MDS
    
    Stress measures the goodness of fit or degree of correspondence
    between distances implied by an MDS mapping and the original
    distances. There are many different variants of Stress; Kruskal's
    Stress or Stress-1 probably being the most popular one.
    """

    def __init__(self, orig_distmat, mds_coords, apply_scaling=True):
        """Setup class by storing pointers to orginal matrix and
        calculating distances implied by MDS

        By default, implied distances (disparities) are scaled to
        match the range of the original distances.

        The calculation of the implied distances is the only really
        time intensive step, so __init__ will take some time for big
        matrices.
        
        Arguments:
        * orig_distmat (numpy array or numpy matrix)
          original distance matrix
        * mds_coords (numpy array)
          mds coordinates
        * apply_scaling (boolean)
          scale distances implied by an MDS mapping to match
          those of original distance matrix
        """

        # don't use coercion here (orig_distmat to matrix, mds_coords
        # to array),to avoid conversion/copy overhead as these
        # arrays/matrices can potentially be quite big
        #
        assert isinstance(orig_distmat, numpy.ndarray), \
               "orig_distmat is not a numpy.ndarray instance"
        assert isinstance(mds_coords, numpy.ndarray), \
               "mds_coords is not a numpy.ndarray instance"

        # check array/matrix shapes and sizes
        assert len(orig_distmat.shape) == 2, \
               "orig_distmat is not a 2D array/matrix."
        assert len(mds_coords.shape) == 2, \
               "mds_coords is not a 2D array."
        assert orig_distmat.shape[0] == mds_coords.shape[0], \
               "orig_distmat and mds_coords do not have the same" \
               " number of rows/objects."
        assert orig_distmat.shape[1] > mds_coords.shape[1], \
               "orig_distmat shape bigger than mds_coords shape." \
               " Possible argument mixup"

        self._orig_distmat = orig_distmat
        
        # compute distances implied by MDS and scale if requested (and
        # necessary) by the ratio of maxima of both matrices
        #
        self._reproduced_distmat = self._calc_pwdist(mds_coords)
        if apply_scaling:
            max_orig = self._orig_distmat.max()
            max_derived = self._reproduced_distmat.max()
            scale = max_derived / max_orig
            if scale != 1.0:
                self._reproduced_distmat = self._reproduced_distmat / scale
        
   
                    
    def calcKruskalStress(self):
        """Calculate Kruskal's Stress AKA Stress-1
        
        Kruskal's Stress or Stress-1 is defined as:
        sqrt( SUM_ij (d'(i,j) - d(i,j))^2 / SUM_ij d(i,j)^2 for i<j

        where d(i,j) is the distance between i and j in the original distance
        matrix (here: self._orig_distmat), and d'(i,j) the distance implied by an
        mds mapping (here: self._reproduced_distmat)

        According to Johnson & Wichern (2002): 'Applied Multivariate
        Statistical Analysis'' p701 (citing Kruskal (1964)) the
        informal interpretation of stress1 is:
        
        Stress [%]  Goodness of fit
        20          Poor
        10          Fair
        5           Good
        2.5         Excellent
        0           Perfect

        AFAIK there is no upper bound for stress1.
        
        Remember that the stress of an MDS mapping also depends on the
        dimensionality of the MDS, i.e. you will get better Stress
        values for higher dimensions.

        Arguments:
        * None
        Returns:
        * Kruskal Stress as float
        """

        # as we have a symmetrical matrix we can restrict ourselves to
        # i<j, but the element-wise computation (looping over i<j) takes
        # ages and gives the same result.

        # use power here for element-wise exponentiation
        # **-operator might fail if we use array instead of matrix
        numerator = numpy.sum(
            numpy.power(self._reproduced_distmat - self._orig_distmat, 2))
        denominator = numpy.sum(numpy.power(self._orig_distmat, 2))
        result = numpy.sqrt(numerator / denominator)

        return result

 

    def calcSstress(self):
        """Calculate SStress
        
        SStress (Takane, 1977) is defined as
        sqrt( SUM_ij (d'(i,j)^2 - d(i,j)^2)^2 / SUM_ij d(i,j)^4 for i<j

        where d(i,j) is the distance between i and j in the original
        matrix (self._orig_distmat), and d'(i,j) the distance implied by an
        mds mapping (self._reproduced_distmat)

        Preferred criterion according to Johnson & Wichern (2002)
        'Applied Multivariate Statistical Analysis'' (p701)

        Value of SStress is always between 0 and 1. Values less then
        0.1 mean good representation.

        Squaring causes S-Stress to emphasize larger dissimliarities more
        than smaller ones (See p204 in 'Modern multidimensional scaling' by
        Borg & Groenen, 1997)

        Arguments:
        * None
        Returns:
        * Sstress as float
        """

        # as we have a symmetrical matrix we can restrict ourselves to
        # i<j, but the element-wise computation (looping over i<j) takes
        # ages and gives the same result. 

        # use power here for element-wise exponentiation
        # **-operator might fail if we use array instead of matrix
        numerator = numpy.power(
            self._reproduced_distmat, 2) - numpy.power(self._orig_distmat, 2)
        numerator = numpy.sum(numpy.power(numerator, 2))
        denominator = numpy.sum(numpy.power(self._orig_distmat, 4))
        result = numpy.sqrt(numerator / denominator)

        return result



    @staticmethod
    def _calc_rowdist(row1, row2):
        """Calculate the euclidean distance between two row vectors.

        NOTE: This function will be called quite often. To safe some
        time, no test will be made if arguments are actually vectors!
        This function will happily work away on full matrices and
        arrays as well!
        
        Arguments:
        * row1 (numpy array or numpy matrix)
          row vector, i.e. shape has to be [1],n or n,[1]
        * row2 (numpy array or numpy matrix)
          row vector, i.e. shape has to be [1],n or n,[1]
        Returns:
        * distance as float
        """
        
        row_diffsq = numpy.power(row1 - row2, 2)
        result = numpy.sqrt(row_diffsq.sum())
        return result



    def _calc_pwdist(self, mds_coords):
        """Compute a full symmetric distance matrix from MDS coordinates
        
        Calculate pairwise distances between m original distances
        (rows) in n-dimensional space (columns)

        Note: this is the same as
            import scipy.spatial.distance as ssd
            ssd.squareform(ssd.pdist(x))

        Arguments:
        * mds_coords (numpy.ndarray)
          mds coordinates, i.e. a recangular 2D array with shape m,n m>n)
        Returns
        * distance as float
        """

        assert isinstance(mds_coords, numpy.ndarray), \
               "mds_coords is not a numpy ndarray"
        
        result = numpy.zeros((mds_coords.shape[0], mds_coords.shape[0]))
        for i in range(mds_coords.shape[0]):
            row_i = mds_coords[i, :]
            for j in range(i+1, mds_coords.shape[0]):
                result[i, j] = self._calc_rowdist(row_i, mds_coords[j, :])
                result[j, i] = result[i, j]
        return result


