#!/usr/bin/env python

import os, tempfile
import numpy as np
try:
    from cogent.util.unit_test import TestCase, main
except ImportError:
    from zenpdb.cogent.util.unit_test import TestCase, main

__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2009, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


class KDTreeTest(TestCase):
    """Tests KD-Trees"""

    def setUp(self):
        self.arr = np.random.random(3000).reshape((1000, 3))
        self.point = np.random.random(3)
        self.center = np.array([0.5, 0.5, 0.5])

    def test_0import(self):
        # sort by name
        """tests if can import ckd3 cython extension."""
        global ckd3
        from cogent.maths.spatial import ckd3
        assert 'KDTree' in dir(ckd3)

    def test_instance(self):
        """check if KDTree instance behaves correctly."""
        kdt = ckd3.KDTree(self.arr)
        self.assertEquals(kdt.dims, 3)
        def assig():
            kdt.dims = 4
        self.assertRaises(TypeError, assig)
        self.assertEquals(kdt.dims, 3)
        self.assertEquals(kdt.pnts, 1000)

    def test_knn(self):
        """testing k-nearest neighbors.
        """
        sqd = np.sum(np.power((self.arr - self.point), 2), axis=1)
        sorted_idx = sqd.argsort()
        kdt = ckd3.KDTree(self.arr)
        points, dists = kdt.knn(self.point, 5)
        self.assertEqualItems(sorted_idx[:5], points)

    def test_rn(self):
        """testing neighbors within radius.
        """
        sqd = np.sum(np.power((self.arr - self.point), 2), axis=1)
        sqd = sqd[sqd <= 0.05]
        sqd.sort()
        kdt = ckd3.KDTree(self.arr)
        points, dists = kdt.rn(self.point, 0.05)
        dists.sort()
        self.assertFloatEqual(dists, sqd)


if __name__ == '__main__':
    main()
