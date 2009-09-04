#! /usr/bin/env python
import unittest, os

from cogent.phylo.distance import *
from cogent.phylo.nj import nj
from cogent.phylo.least_squares import wls
from cogent import LoadSeqs, LoadTree
from cogent.evolve.models import JC69, HKY85, F81
from cogent.phylo.consensus import majorityRule, weightedMajorityRule

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield",\
        "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def Tree(t):
    return LoadTree(treestring=t)

class ConsensusTests(unittest.TestCase):
    def setUp(self):
        self.trees = [
                (1, Tree("((a,b),(c,d));")),
                (1, Tree("((a,b),(c,d));")),
                (1, Tree("((a,c),(b,d));")),
                (1, Tree("((a,b),c,d);")),]
    
    def test_majorityRule(self):
        """Tests for majority rule consensus trees"""
        trees = [t for (w,t) in self.trees]
        outtrees = majorityRule(trees, strict=False)
        self.assertEqual(len(outtrees), 1)
        self.assert_(outtrees[0].sameTopology(Tree("((c,d),(a,b));")))
        outtrees = majorityRule(trees, strict=True)
        self.assertEqual(len(outtrees), 1)
        self.assert_(outtrees[0].sameTopology(Tree("(c,d,(a,b));")))

class TreeReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.tree = LoadTree(treestring='((a:3,b:4):2,(c:6,d:7):30,e:5)')
        self.dists = self.tree.getDistances()
        
    def assertTreeDistancesEqual(self, t1, t2):
        d1 = t1.getDistances()
        d2 = t2.getDistances()
        for key in d2:
            self.assertAlmostEqual(d1[key], d2[key])

    def test_nj(self):
        """testing nj"""
        reconstructed = nj(self.dists)
        self.assertTreeDistancesEqual(self.tree, reconstructed)
        
    def test_wls(self):
        """testing wls"""
        reconstructed = wls(self.dists)
        self.assertTreeDistancesEqual(self.tree, reconstructed)

    def test_truncated_wls(self):
        """testing wls with order option"""
        order = ['e', 'b', 'c', 'd']
        reconstructed = wls(self.dists, order=order)
        self.assertEqual(set(reconstructed.getTipNames()), set(order))

    def test_limited_wls(self):
        """testing (well, exercising at least), wls with constrained start"""
        init = LoadTree(treestring='((a,c),b,d)')
        reconstructed = wls(self.dists, start=init)
        self.assertEqual(len(reconstructed.getTipNames()), 5)
        init2 = LoadTree(treestring='((a,d),b,c)')
        reconstructed = wls(self.dists, start=[init, init2])
        self.assertEqual(len(reconstructed.getTipNames()), 5)
        init3 = LoadTree(treestring='((a,d),b,e)')
        self.assertRaises(Exception, wls, self.dists, start=[init, init3])

class DistancesTests(unittest.TestCase):
    def setUp(self):
        self.al = LoadSeqs(data = {'a':'GTACGTACGATC',
                            'b':'GTACGTACGTAC',
                            'c':'GTACGTACGTTC',
                            'e':'GTACGTACTGGT'})
        self.collection = LoadSeqs(data = {'a':'GTACGTACGATC',
                            'b':'GTACGTACGTAC',
                            'c':'GTACGTACGTTC',
                            'e':'GTACGTACTGGT'}, aligned=False)
    
    def test_EstimateDistances(self):
        """testing (well, exercising at least), EstimateDistances"""
        d = EstimateDistances(self.al, JC69())
        d.run()
        canned_result = {('b', 'e'): 0.440840,
                        ('c', 'e'): 0.440840,
                        ('a', 'c'): 0.088337,
                        ('a', 'b'): 0.188486,
                        ('a', 'e'): 0.440840,
                        ('c', 'b'): 0.0883373}
        result = d.getPairwiseDistances()
        for key in canned_result:
            self.assertAlmostEqual(canned_result[key],result[key],4)
        
        # excercise writing to file
        d.writeToFile('junk.txt')
        try:
            os.remove('junk.txt')
        except OSError:
            pass # probably parallel
    
    def test_EstimateDistancesWithMotifProbs(self):
        """EstimateDistances with supplied motif probs"""
        motif_probs= {'A':0.1,'C':0.2,'G':0.2,'T':0.5}
        d = EstimateDistances(self.al, HKY85(), motif_probs=motif_probs)
        d.run()
        canned_result = {('a', 'c'): 0.07537,
                        ('c', 'b'): 0.07537,
                        ('a', 'e'): 0.39921,
                        ('a', 'b'): 0.15096,
                        ('b', 'e'): 0.39921,
                        ('c', 'e'): 0.37243}
        result = d.getPairwiseDistances()
        for key in canned_result:
            self.assertAlmostEqual(canned_result[key],result[key],4)
    
    def test_EstimateDistances_fromThreeway(self):
        """testing (well, exercising at least), EsimateDistances fromThreeway"""
        d = EstimateDistances(self.al, JC69(), threeway=True)
        d.run()
        canned_result = {('b', 'e'): 0.495312,
                        ('c', 'e'): 0.479380,
                        ('a', 'c'): 0.089934,
                        ('a', 'b'): 0.190021,
                        ('a', 'e'): 0.495305,
                        ('c', 'b'): 0.0899339}
        result = d.getPairwiseDistances(summary_function="mean")
        for key in canned_result:
            self.assertAlmostEqual(canned_result[key],result[key],4)
    
    def test_EstimateDistances_fromUnaligned(self):
        """Excercising estimate distances from unaligned sequences"""
        d = EstimateDistances(self.collection, JC69(), do_pair_align=True,
                                rigorous_align=True)
        d.run()
        canned_result = {('b', 'e'): 0.440840,
                        ('c', 'e'): 0.440840,
                        ('a', 'c'): 0.088337,
                        ('a', 'b'): 0.188486,
                        ('a', 'e'): 0.440840,
                        ('c', 'b'): 0.0883373}
        result = d.getPairwiseDistances()
        for key in canned_result:
            self.assertAlmostEqual(canned_result[key],result[key],4)
        
        d = EstimateDistances(self.collection, JC69(), do_pair_align=True,
                                rigorous_align=False)
        d.run()
        canned_result = {('b', 'e'): 0.440840,
                        ('c', 'e'): 0.440840,
                        ('a', 'c'): 0.088337,
                        ('a', 'b'): 0.188486,
                        ('a', 'e'): 0.440840,
                        ('c', 'b'): 0.0883373}
        result = d.getPairwiseDistances()
        for key in canned_result:
            self.assertAlmostEqual(canned_result[key],result[key],4)
    
    def test_EstimateDistances_other_model_params(self):
        """test getting other model params from EstimateDistances"""
        d = EstimateDistances(self.al, HKY85(), est_params=['kappa'])
        d.run()
        # this will be a Number object with Mean, Median etc ..
        kappa = d.getParamValues('kappa')
        self.assertAlmostEqual(kappa.Mean, 0.8939, 4)
        # this will be a dict with pairwise instances, it's called by the above
        # method, so the correctness of it's values is already checked
        kappa = d.getPairwiseParam('kappa')
    
    def test_EstimateDistances_modify_lf(self):
        """tests modifying the lf"""
        def constrain_fit(lf):
            lf.setParamRule('kappa', is_const=True)
            lf.optimise(local=True, show_progress=False)
            return lf
        
        d = EstimateDistances(self.al, HKY85(), modify_lf=constrain_fit)
        d.run()
        result = d.getPairwiseDistances()
        d = EstimateDistances(self.al, F81())
        d.run()
        expect = d.getPairwiseDistances()
        for key in result:
            self.assertAlmostEqual(expect[key], result[key],4)
        
    

if __name__ == '__main__':
    unittest.main()
