#!/usr/bin/env python

"""Unit tests of microarray_normalize.py: code for normalizing microarrays.
"""
from cogent.seqsim.microarray_normalize import (zscores, logzscores, ranks, 
    quantiles, 
    make_quantile_normalizer, make_normal_quantile_normalizer,
    make_empirical_quantile_normalizer,
    geometric_mean )
from cogent.util.unit_test import TestCase, main
from numpy import array, arange, reshape, log2

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Micah Hamady", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class microarray_normalize_tests(TestCase):
    """Tests of top-level functions."""

    def test_zscores(self):
        """zscores should convert array to zscores within each column"""
        a = reshape(arange(15),(5,3))
        z = zscores(a)
        self.assertEqual(z[2], array([0,0,0])) #middle should be mean
        self.assertFloatEqual(z[0], [-1.41421356]*3)
        #check that it works when arrays aren't sorted
        a[0] = a[-1]
        a[1] = a[-2]
        a[0, -1] = 50
        z = zscores(a)
        self.assertEqual(z[0,0],z[-1,0])
        self.assertFloatEqual(z[0,-1], 1.9853692256351525)
        self.assertFloatEqual(z[-1,-1], -0.30544141932848506)

    def test_logzscores(self):
        """logzscores should perform zscores on log of a"""
        a = reshape(arange(1,16),(5,3)) #won't work with zero value
        self.assertFloatEqual(logzscores(a), zscores(log2(a)))
    
    def test_ranks(self):
        """ranks should convert array to ranks within each column"""
        a = array([[10,20,30],[20,10,50],[30,5,10]])
        r = ranks(a)
        self.assertEqual(r, array([[0,2,1],[1,1,2],[2,0,0]]))

    def test_quantiles(self):
        """quantiles should convert array to quantiles within each column"""
        a = array([[10,20,30],[20,10,50],[30,5,10],[40,40,40]])
        q = quantiles(a)
        self.assertEqual(q, \
            array([[0,.5,.25],[.25,.25,.75],[.5,0,0],[.75,.75,.5]]))

    def test_make_quantile_normalizer(self):
        """make_quantile_normalizer should sample from right distribution."""
        dist = array([1,2,3,4])
        qn = make_quantile_normalizer(dist)
        a = array([[10,20,30],[20,10,50],[30,5,10],[40,40,40]])
        q = qn(a)
        self.assertEqual(q, \
            array([[1,3,2],[2,2,4],[3,1,1],[4,4,3]]))
        #check that it works when they don't match in size exactly
        dist = array([2,4,6,7,8,8,8,8])
        qn = make_quantile_normalizer(dist)
        a = array([[10,20,30],[20,10,50],[30,5,10],[40,40,40]])
        q = qn(a)
        self.assertEqual(q, \
            array([[2,8,6],[6,6,8],[8,2,2],[8,8,8]]))

    def test_make_normal_quantile_normalizer(self):
        """make_normal_quantile_normalizer should sample from normal dist."""
        nqn = make_normal_quantile_normalizer(20, 10)
        a = array([[10,20,30],[20,10,50],[30,5,10],[40,40,40]])
        q = nqn(a)
        exp = array([[-289.02323062,   20.        ,  -47.44897502],
       [ -47.44897502,  -47.44897502,   87.44897502],
       [  20.        , -289.02323062, -289.02323062],
       [  87.44897502,   87.44897502,   20.        ]])
        self.assertFloatEqual(q, exp)

    def test_make_empirical_quantile_normalizer(self):
        """make_empirical_quantile_normalizer should convert a to dist of data"""
        dist = array([4,2,3,1]) #note: out of order
        qn = make_empirical_quantile_normalizer(dist)
        a = array([[10,20,30],[20,10,50],[30,5,10],[40,40,40]])
        q = qn(a)
        self.assertEqual(q, \
            array([[1,3,2],[2,2,4],[3,1,1],[4,4,3]]))

    def test_geometric_mean(self):
        """geometric_mean should return geometric mean."""
        a = array([1.05, 1.2, .96])
        gmean = geometric_mean(a)
        self.assertFloatEqual(gmean, 1.065484802091121)

        
  
if __name__ == "__main__":
    main()
