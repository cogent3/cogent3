#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from numpy import array, sqrt, size
from cogent.cluster.nmds import NMDS, metaNMDS
from cogent.maths.distance_transform import dist_euclidean

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

class NMDSTests(TestCase):
    """test the nonmetric_scaling module, using floating point numpy arrays
    """

    def setUp(self):
        """creates inputs"""
        self.mtx = array([[0,3,4,8],
                [3,0,1,27],
                [4,1,0,3.5],
                [8,27,3.5,0]],'d')
        self.nm = NMDS(self.mtx, verbosity=0)
                        
    def test_getStress(self):
        """stress should be small
        
        this is preliminary, better to check for convergence to similar states
        with random starting points enabled"""
        stress = self.nm.getStress()
        self.assertLessThan(stress, 1e-1)

    def test_getPoints(self):
        """points should be of the right number and dimensionality
        
        this is preliminary, better to check for convergence to similar states
        with random starting points enabled"""
        pts = self.nm.getPoints()
        self.assertEqual(size(pts, 0), 4)
        self.assertEqual(size(pts, 1), 2)

    def test_2(self):
        """l19 data should give stress below .13"""
        ptmtx = array(
            [[7,1,0,0,0,0,0,0,0],
            [4,2,0,0,0,1,0,0,0],
            [2,4,0,0,0,1,0,0,0],
            [1,7,0,0,0,0,0,0,0],
            [0,8,0,0,0,0,0,0,0],
            [0,7,1,0,0,0,0,0,0],#idx 5
            [0,4,2,0,0,0,2,0,0],
            [0,2,4,0,0,0,1,0,0],
            [0,1,7,0,0,0,0,0,0],
            [0,0,8,0,0,0,0,0,0],
            [0,0,7,1,0,0,0,0,0],#idx 10
            [0,0,4,2,0,0,0,3,0],
            [0,0,2,4,0,0,0,1,0],
            [0,0,1,7,0,0,0,0,0],
            [0,0,0,8,0,0,0,0,0],
            [0,0,0,7,1,0,0,0,0],#idx 15
            [0,0,0,4,2,0,0,0,4],
            [0,0,0,2,4,0,0,0,1],
            [0,0,0,1,7,0,0,0,0]], 'float')
        distmtx = dist_euclidean(ptmtx)
        nm = NMDS(distmtx, verbosity=0)
        self.assertLessThan(nm.getStress(), .13)

    def test_metaNMDS(self):
        """l19 data should give stress below .13"""
        ptmtx = array(
            [[7,1,0,0,0,0,0,0,0],
            [4,2,0,0,0,1,0,0,0],
            [2,4,0,0,0,1,0,0,0],
            [1,7,0,0,0,0,0,0,0],
            [0,8,0,0,0,0,0,0,0],
            [0,7,1,0,0,0,0,0,0],#idx 5
            [0,4,2,0,0,0,2,0,0],
            [0,2,4,0,0,0,1,0,0],
            [0,1,7,0,0,0,0,0,0],
            [0,0,8,0,0,0,0,0,0],
            [0,0,7,1,0,0,0,0,0],#idx 10
            [0,0,4,2,0,0,0,3,0],
            [0,0,2,4,0,0,0,1,0],
            [0,0,1,7,0,0,0,0,0],
            [0,0,0,8,0,0,0,0,0],
            [0,0,0,7,1,0,0,0,0],#idx 15
            [0,0,0,4,2,0,0,0,4],
            [0,0,0,2,4,0,0,0,1],
            [0,0,0,1,7,0,0,0,0]], 'float')
        distmtx = dist_euclidean(ptmtx)
        nm = metaNMDS(1, distmtx, verbosity=0)
        self.assertLessThan(nm.getStress(), .13)

if __name__ == '__main__':
       main()
