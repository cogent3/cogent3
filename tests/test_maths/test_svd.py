#!/usr/bin/env python
"""Unit tests for the svd-supporting functionality."""

from cogent.util.unit_test import TestCase, main
from cogent.maths.svd import ratio_two_best, ratio_best_to_sum, \
     euclidean_distance, euclidean_norm, _dists_from_mean_slow, \
     dists_from_v, weiss, three_item_combos, two_item_combos
from numpy import array, sqrt

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class functionTests(TestCase):
    """Tests of top-level functions."""
    
    def test_ratio_two_best(self):
        """ratio_two_best should return ratio of two biggest items in list"""
        v = array([3, 2, 5, 2, 4, 10, 3])
        self.assertEqual(ratio_two_best(v), 2)
        #should return 1 if items the same
        v = array([2,2,2,2,2])
        self.assertEqual(ratio_two_best(v), 1)
        #check that it works on floating-point
        v = array([3,2,1])
        self.assertEqual(ratio_two_best(v), 1.5)

    def test_ratio_best_to_sum(self):
        """ratio_best_to_sum should return ratio of biggest item to sum"""
        v = [3, 2, 5, 2, 4, 10, 3]
        self.assertFloatEqual(ratio_best_to_sum(v), 10/29.0)
        v = [2,2,2,2,2]
        self.assertEqual(ratio_best_to_sum(v), 2/10.0)
        #check that it works on floating-point
        v = [3,2,1]
        self.assertEqual(ratio_best_to_sum(v), 0.5)

    def test_euclidean_distance(self):
        """euclidean_distance should return distance between two points"""
        first = array([2, 3, 4])
        second = array([4, 8, 10])
        self.assertEqual(euclidean_distance(first, first), 0)
        self.assertEqual(euclidean_distance(second, second), 0)
        self.assertFloatEqual(euclidean_distance(first, second), sqrt(65))
        self.assertFloatEqual(euclidean_distance(second, first), sqrt(65))

    def test_euclidean_norm(self):
        """euclidean_norm should match hand-calculated results"""
        first = array([3,4])
        self.assertEqual(euclidean_norm(first), 5)

    def test_dists_from_mean_slow(self):
        """_dists_from_mean_slow should return distance of each item from mean"""
        m = [[1,2,3,4],[2,3,4,5],[0,1,2,3]]
        self.assertEqual(_dists_from_mean_slow(m), array([0.0,2.0,2.0]))

    def test_dists_from_v(self):
        """dists_from_v should return distance of each item from v, or mean"""
        m = [[1,2,3,4],[2,3,4,5],[0,1,2,3]]
        #should calculate distances from mean by default
        self.assertEqual(dists_from_v(m), array([0.0,2.0,2.0]))
        #should caculate distances from vector if supplied
        v = array([2,2,2,3])
        self.assertEqual(dists_from_v(m, v), sqrt(array([3,9,5])))
        
    def test_weiss(self):
        """weiss should perform weiss calculation correctly"""
        e = array([12.0, 5.0, 0.1, 1e-3, 1e-15])
        self.assertFloatEqual(weiss(e), 4.453018506827001)

    def test_three_item_combos(self):
        """three_item_combos should return items in correct order"""
        items = list(three_item_combos('abcde'))
        self.assertEqual(items, map(tuple, \
            ['abc','abd','abe','acd','ace','ade','bcd','bce','bde','cde']))

    def test_two_item_combos(self):
        """two_item_combos should return items in correct order"""
        items = list(two_item_combos('abcd'))
        self.assertEqual(items, map(tuple, ['ab','ac','ad','bc','bd','cd']))

    def test_pca_qs(self):
        """pca_qs not tested b/c it just wraps eigenvalues(corrcoef(qs))"""
        pass

    def test_pca_cov_qs(self):
        """pca_cov_qs not tested b/c it just wraps eigenvalues(cov(qs))"""
        pass

    def test_svd_qs(self):
        """svd_qs not tested b/c it just wraps singular_value_decompositon(qs)"""
        pass

if __name__ == '__main__':
    main()
