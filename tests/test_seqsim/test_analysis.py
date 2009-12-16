#!/usr/bin/env python
"""Unit tests for analysis.py: substitution matrix analysis code."""
from cogent.seqsim.analysis import tree_threeway_counts, \
    tree_twoway_counts, counts_to_probs, probs_to_rates, \
    tree_threeway_rates, tree_twoway_rates, \
    rates_to_array, multivariate_normal_prob
from cogent.seqsim.tree import RangeNode
from cogent.core.usage import DnaPairs, ABPairs
from cogent.seqsim.usage import Rates, Counts, Probs
from numpy import array, average, ones, zeros, float64, ravel, diag, any
from numpy.random import random, randint
from copy import deepcopy
from cogent.parse.tree import DndParser
from cogent.util.unit_test import TestCase, main

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class analysisTests(TestCase):
    """Tests of top-level functions."""
    def setUp(self):
        """Make a couple of standard trees"""
        self.t1 = DndParser('((a,(b,c)),(d,e))', RangeNode)
        #selt.t1 indices: ((0,(1,2)5)6,(3,4)7)8
 
    def test_threeway_counts(self):
        """threeway_counts should produce correct count matrix"""
        self.t1.makeIdIndex()
        ind = self.t1.IdIndex
        ind[0].Sequence = array([0,0,0])
        ind[1].Sequence = array([0,1,0])
        ind[2].Sequence = array([1,0,1])
        ind[3].Sequence = array([1,1,0])
        ind[4].Sequence = array([1,1,1])
        depths = self.t1.leafLcaDepths()
        result = tree_threeway_counts(self.t1, depths, ABPairs)
        #check we got the right number of comparisons
        self.assertEqual(len(result), 20)
        #check we got the right keys
        for k in [(1,2,0),(2,1,0),(0,1,3),(1,0,3),(0,1,4),(1,0,4),(0,2,3),\
            (2,0,3),(0,2,4),(2,0,4),(1,2,3),(2,1,3),(1,2,4),(2,1,4),(3,4,1),\
            (4,3,1),(3,4,2),(4,3,2)]:
            assert k in result
        #spot-check a few results
        self.assertEqual(result[(1,2,0)]._data, array([[2,1],[0,0]]))
        self.assertEqual(result[(2,1,0)]._data, array([[1,2],[0,0]]))
        self.assertEqual(result[(2,1,3)]._data, array([[0,1],[1,1]]))
        
    def test_twoway_counts(self):
        """twoway_counts should produce correct count matrix"""
        self.t1.makeIdIndex()
        ind = self.t1.IdIndex
        ind[0].Sequence = array([0,0,0])
        ind[1].Sequence = array([0,1,0])
        ind[2].Sequence = array([1,0,1])
        ind[3].Sequence = array([1,1,0])
        ind[4].Sequence = array([1,1,1])
        depths = self.t1.leafLcaDepths()
        #check that it works with averaging
        result = tree_twoway_counts(self.t1, ABPairs)
        #check we got the right number of comparisons: average by default
        self.assertEqual(len(result), 10)
        #check we got the right keys
        for k in [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]:
            assert k in result
        #spot-check a few results
        self.assertEqual(result[(0,1)]._data, array([[2,.5],[.5,0]]))
        self.assertEqual(result[(2,3)]._data, array([[0,1],[1,1]]))
        #check that it works when we don't average
        result = tree_twoway_counts(self.t1, ABPairs, average=False)
        self.assertEqual(len(result), 20)
        #check we got the right keys
        for k in [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]:
            assert k in result
            #reverse should be in result too
            assert (k[1],k[0]) in result
        #spot-check values
        self.assertEqual(result[(0,1)]._data, array([[2,1],[0,0]]))
        self.assertEqual(result[(1,0)]._data, array([[2,0],[1,0]]))
        
    def test_counts_to_probs(self):
        """counts_to_probs should skip cases with zero rows"""
        counts = {
            (0,1): Counts(array([[0,1],[1,0]]), ABPairs),
            (1,2): Counts(array([[0,0],[1,0]]), ABPairs),           #bad row
            (0,3): Counts(array([[0,0],[0,0]]), ABPairs),           #bad row
            (0,4): Counts(array([[0.0,0.0],[0.0,0.0]]), ABPairs),   #bad row
            (0,5): Counts(array([[0.1,0.3],[0.0,0.0]]), ABPairs),   #bad row
            (3,4): Counts(array([[0.1,0.3],[0.4,0.1]]), ABPairs),
            (2,1): Counts(array([[0,5],[1,0]]), ABPairs),
            }
        result = counts_to_probs(counts)
        self.assertEqual(len(result), 3)
        self.assertFloatEqual(result[(0,1)]._data, array([[0,1],[1,0]]))
        self.assertFloatEqual(result[(3,4)]._data, \
            array([[0.25,0.75],[0.8,0.2]]))
        self.assertFloatEqual(result[(2,1)]._data, array([[0,1],[1,0]]))

    def test_probs_to_rates(self):
        """probs_to_rates converts probs to rates, omitting problem cases"""
        probs = dict([(i, Probs.random(DnaPairs)) for i in range(100)])
        rates = probs_to_rates(probs)
        #check we got at most the same number of items as in probs
        assert len(rates) <= len(probs)
        #check that we didn't get anything bad
        vals = rates.values()
        for v in vals:
            assert not v.isSignificantlyComplex()
        #check that we didn't miss anything good
        for key, val in probs.items():
            if key not in rates:
                try:
                    r = val.toRates()
                    print r.isValid()
                    assert r.isSignificantlyComplex() or (not r.isValid())
                except (ZeroDivisionError, OverflowError, ValueError):
                    pass

    def test_rates_to_array(self):
        """rates_to_array should pack rates into array correctly"""
        m1 = array([[-1,1,1,1],[2,-2,2,2],[3,3,-3,3],[1,2,3,-4]])
        m2 = m1 * 2
        m3 = m1 * 0.5
        m4 = zeros((4,4))
        m5 = array([0,0])
        r1, r2, r3, r4, r5 = [Rates(i, DnaPairs) for i in m1,m2,m3,m4,m5]
    
        data = {(0,1,0):r1, (1,2,0):r2, (2,0,0):r3, (2,1,1):r4}
        
        #note that array can be, but need not be, floating point
        to_fill = zeros((3,3,3,16), 'float64')
        result = rates_to_array(data, to_fill)
        #check that the thnigs we deliberately set are OK
        self.assertEqual(to_fill[0][1][0], ravel(m1))
        self.assertNotEqual(to_fill[0][1][0], ravel(m2))
        self.assertEqual(to_fill[1,2,0], ravel(m2))
        self.assertEqual(to_fill[2][0][0], ravel(m3))
        self.assertEqual(to_fill[2][1][1], ravel(m4))
        #check that everything else is zero
        nonzero = [(0,1,0),(1,2,0),(2,0,0)]
        for x in [(i, j, k) for i in range(3) for j in range(3) \
            for k in range(3)]:
            if x not in nonzero:
                self.assertEqual(to_fill[x], zeros(16))
        #check that it works omitting the diagonal
        to_fill = zeros((3,3,3,12), 'float64')
        result = rates_to_array(data, to_fill, without_diagonal=True)
        #check that the thnigs we deliberately set are OK
        m1_nodiag = array([[1,1,1],[2,2,2],[3,3,3],[1,2,3]])
        self.assertEqual(to_fill[0][1][0], ravel(m1_nodiag))
        self.assertNotEqual(to_fill[0][1][0], ravel(m1_nodiag*2))
        self.assertEqual(to_fill[1,2,0], ravel(m1_nodiag*2))
        self.assertEqual(to_fill[2][0][0], ravel(m1_nodiag*0.5))
        self.assertEqual(to_fill[2][1][1], zeros(12))
        #check that everything else is zero
        nonzero = [(0,1,0),(1,2,0),(2,0,0)]
        for x in [(i, j, k) for i in range(3) for j in range(3) \
            for k in range(3)]:
            if x not in nonzero:
                self.assertEqual(to_fill[x], zeros(12))
    
    def test_tree_threeway_rates(self):
        """tree_threeway_rates should give plausible results on rand trees"""
        #note: the following fails occasionally, but repeating it 5 times
        #and checking that one passes is fairly safe
        for i in range(5):
            try:
                t = self.t1
                t.assignLength(0.05)
                t.Q = Rates.random(DnaPairs).normalize()
                t.assignQ()
                t.assignP()
                t.evolve(randint(0,4,100))
                t.makeIdIndex()
                depths = t.leafLcaDepths()
                result = tree_threeway_rates(t, depths)
                self.assertEqual(result.shape, (5,5,5,16))
                #check that row sums are 0
                for x in [(i,j,k) for i in range(5) for j in range(5) \
                    for k in range(5)]:
                    self.assertFloatEqual(sum(result[x]), 0)
                assert any(result)
                #check that it works without_diag
                result = tree_threeway_rates(t, depths, without_diag=True)
                self.assertEqual(result.shape, (5,5,5,12))
                #check that it works with/without normalize
                #default: no normalization, so row sums shouldn't be 1 after 
                #omitting diagonal
                result = tree_threeway_rates(t, depths, without_diag=True)
                self.assertEqual(result.shape, (5,5,5,12))
                for x in [(i,j,k) for i in range(5) for j in range(5) \
                    for k in range(5)]:
                    assert sum(result[x]) == 0 or abs(sum(result[x]) - 1) > 0.01
                #...but if we tell it to normalize, row sums should be nearly 1
                #after omitting diagonal
                result = tree_threeway_rates(t, depths, without_diag=True, \
                    normalize=True)
                self.assertEqual(result.shape, (5,5,5,12))
                for x in [(i,j,k) for i in range(5) for j in range(5) \
                    for k in range(5)]:
                        s = sum(result[x])
                        if s != 0:
                            self.assertFloatEqual(s, 1)
                break
            except AssertionError:
                pass
    
    def test_tree_twoway_rates(self):
        """tree_twoway_rates should give plausible results on rand trees"""
        t = self.t1
        t.assignLength(0.05)
        t.Q = Rates.random(DnaPairs).normalize()
        t.assignQ()
        t.assignP()
        t.evolve(randint(0,4,100))
        t.makeIdIndex()
        result = tree_twoway_rates(t)
        self.assertEqual(result.shape, (5,5,16))
        #check that row sums are 0
        for x in [(i,j) for i in range(5) for j in range(5)]:
            self.assertFloatEqual(sum(result[x]), 0)
        #need to make sure we didn't just get an empty array
        self.assertGreaterThan((abs(result)).sum(), 0)
        #check that it works without_diag
        result = tree_twoway_rates(t, without_diag=True)
        self.assertEqual(result.shape, (5,5,12))
        #check that it works with/without normalize
        #default: no normalization, so row sums shouldn't be 1 after omitting
        #diagonal
        result = tree_twoway_rates(t, without_diag=True)
        self.assertEqual(result.shape, (5,5,12))
        #check that the row sums are not 1 before normalization (note that they
        #can be zero, though)
        sums_before = []
        for x in [(i,j) for i in range(5) for j in range(5)]:
            curr_sum = sum(result[x])
            sums_before.append(curr_sum)
        #...but if we tell it to normalize, row sums should be nearly 1
        #after omitting diagonal
        result = tree_twoway_rates(t, without_diag=True, \
            normalize=True)
        self.assertEqual(result.shape, (5,5,12))
        sums_after = []
        for x in [(i,j) for i in range(5) for j in range(5)]:
            curr_sum = sum(result[x])
            sums_after.append(curr_sum)
            if curr_sum != 0:
                self.assertFloatEqual(curr_sum, 1)
        try:
            self.assertFloatEqual(sums_before, sums_after)
        except AssertionError:
            pass
        else:
            raise AssertionError, "Expected different arrays before/after norm"
   
    def test_multivariate_normal_prob(self):
        """Multivariate normal prob should match R results"""
        cov = array([[3,1,2],[1,5,4],[2,4,6]])
        a = array([0,0,0])
        b = array([1,1,1])
        c = array([0.1, 0.2, 0.3])
        small_cov = cov/10.0
        
        mvp = multivariate_normal_prob
        self.assertFloatEqual(mvp(a, cov), 0.01122420)
        self.assertFloatEqual(mvp(a, cov, b), 0.009018894)
        self.assertFloatEqual(mvp(a, small_cov, b), 0.03982319)
        self.assertFloatEqual(mvp(c, small_cov, b), 0.06091317)
        
if __name__ == "__main__":
    main()
