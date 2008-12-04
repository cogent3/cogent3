#! /usr/bin/env python
import unittest, os

from cogent.phylo.distance import *
from cogent.phylo.nj import nj
from cogent.phylo.least_squares import wls
from cogent import LoadSeqs, LoadTree
from cogent.evolve.models import JC69, HKY85

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class TreeReconstructionTests(unittest.TestCase):
    def setUp(self):
        pass
        
    def _test_tree(self, method, treestring):
        t = LoadTree(treestring=treestring)
        t_distances = t.getDistances()
        reconstructed = method(t_distances)
        distances = reconstructed.getDistances()
        for key in t_distances:
            self.assertAlmostEqual(t_distances[key], distances[key])

    def _test_phylo_method(self, method):
        """testing (well, exercising at least), nj"""
        self._test_tree(method, '((a:3,b:4):20,(c:6,d:7):30,e:5)')
        self._test_tree(method, '((a:3,b:4):0,(c:6,d:7):30,e:5)')
        self._test_tree(method, '((a:3,b:4,c:6,d:7):30,e:5)')

    def test_nj(self):
        """testing (well, exercising at least), nj"""
        self._test_phylo_method(nj)
        
    def test_wls(self):
        """testing (well, exercising at least), wls"""
        self._test_phylo_method(wls)

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
        os.remove('junk.txt')
        
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

if __name__ == '__main__':
    unittest.main()