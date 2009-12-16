#!/usr/bin/env python
#file test_parse.py
from numpy import array
from cogent.util.unit_test import TestCase, main
from cogent.maths.stats.rarefaction import subsample, naive_histogram, wrap_numpy_histogram, \
    rarefaction

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_subsample(self):
        """subsample should return a random subsample of a vector"""
        a = array([0,5,0])
        self.assertEqual(subsample(a,5), array([0,5,0]))
        self.assertEqual(subsample(a,2), array([0,2,0]))
        b = array([2,0,1])
        obs = subsample(b,2)
        assert (obs == array([1,0,1])).all() or (obs ==  array([2,0,0])).all()

    def test_naive_histogram(self):
        """naive_histogram should produce expected result"""
        vals = array([1,0,0,3])
        self.assertEqual(naive_histogram(vals), array([2,1,0,1]))
        self.assertEqual(naive_histogram(vals, 4), array([2,1,0,1,0]))

    def test_wrap_numpy_histogram(self):
        """wrap_numpy_histogram should provide expected result"""
        vals = array([1,0,0,3])
        h_f = wrap_numpy_histogram(3)
        self.assertEqual(h_f(vals), array([2,1,0,1]))
        h_f = wrap_numpy_histogram(4)
        self.assertEqual(h_f(vals, 4), array([2,1,0,1,0]))

    def test_rarefaction(self):
        """rarefaction should produce expected curve"""
        vals = array([5,0,0,3,0,10], dtype=int)
        res = [r.copy() for r in rarefaction(vals, stride=1)]
        self.assertEqual(len(res), 18)
        for i, r in enumerate(res):
            self.assertEqual(r.sum(), i+1)
            #make sure we didn't add any bad counts
            for pos in [1,2,4]:
                self.assertEqual(r[pos], 0)
        #when we get to end should recapture orig vals
        self.assertEqual(r, vals)
        res = [r.copy() for r in rarefaction(vals, stride=3)]
        self.assertEqual(len(res), 6)
        for i, r in enumerate(res):
            self.assertEqual(r.sum(), 3*(i+1))
            #make sure we didn't add any bad counts
            for pos in [1,2,4]:
                self.assertEqual(r[pos], 0)
        #when we get to end should recapture orig vals
        self.assertEqual(r, vals)

        #repeat everything above using alt. input format
        orig_vals = vals.copy()
        vals = array([0,0,0,0,0,3,3,3,5,5,5,5,5,5,5,5,5,5], dtype=int)
        res = [r.copy() for r in rarefaction(vals, stride=1, is_counts=False)]
        self.assertEqual(len(res), 18)
        for i, r in enumerate(res):
            self.assertEqual(r.sum(), i+1)
            #make sure we didn't add any bad counts
            for pos in [1,2,4]:
                self.assertEqual(r[pos], 0)
        #when we get to end should recapture orig vals
        self.assertEqual(r, orig_vals)
        res = [r.copy() for r in rarefaction(vals, stride=3, is_counts=False)]
        self.assertEqual(len(res), 6)
        for i, r in enumerate(res):
            self.assertEqual(r.sum(), 3*(i+1))
            #make sure we didn't add any bad counts
            for pos in [1,2,4]:
                self.assertEqual(r[pos], 0)
        #when we get to end should recapture orig vals
        self.assertEqual(r, orig_vals)


if __name__ =='__main__':
    main()
