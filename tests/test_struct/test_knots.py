#!/usr/bin/env python
# test_knots.py

"""Provides tests for classes and functions in the file knots.py
"""

from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.util.dict2d import Dict2D
from cogent.struct.rna2d import Pairs
from cogent.struct.knots import PairedRegion, PairedRegionFromPairs,\
    PairedRegions, PairedRegionsFromPairs, ConflictMatrix,\
    opt_all, contains_true, empty_matrix,\
    pick_multi_best, dp_matrix_multi, matrix_solutions,\
    opt_single_random, opt_single_property,\
    inc_order, inc_length, inc_range,\
    find_max_conflicts, find_min_gain,\
    conflict_elimination, add_back_non_conflicting,\
    num_bps, hydrogen_bonds,\
    nussinov_fill, nussinov_traceback, nussinov_restricted

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Sandra Smit, Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"


class PairedRegionTests(TestCase):
    """Tests for PairedRegion class"""

    def test_init_valid(self):
        """PairedRegion __init__: should work as expected on valid input
        """
        pr = PairedRegion(3,10,2)
        self.assertEqual(pr.Start, 3)
        self.assertEqual(pr.End, 10)
        self.assertEqual(pr.Length, 2)
        self.assertEqual(pr.Pairs, [(3,10),(4,9)])

        pr = PairedRegion(3,10,2, Id=0)
        self.assertEqual(pr.Id, 0)
        pr = PairedRegion(3,10,2, Id='A')
        self.assertEqual(pr.Id, 'A')

        self.assertRaises(ValueError, PairedRegion, 4, 10, 0)

    def test_init_weird(self):
        """PairedRegion __init__: no error checking
        """
        pr = PairedRegion(3,6,4)
        self.assertEqual(pr.Start, 3)
        self.assertEqual(pr.End, 6)
        self.assertEqual(pr.Length, 4)
        self.assertEqual(pr.Pairs, [(3,6),(4,5),(5,4),(6,3)])

    def test_str(self):
        """PairedRegion __str__: should print pairs"""
        pr = PairedRegion(3,10,2)
        p = Pairs([(3,10),(4,9)])
        self.assertEqual(str(pr), str(p))
    
    def test_len(self):
        """PairedRegion __len__: should return number of pairs"""
        pr = PairedRegion(3,10,2)
        self.assertEqual(len(pr), 2)

    def test_eq(self):
        """PairedRegion __eq__: should use pairs and IDs"""
        pr1 = PairedRegion(3,10,2)
        pr2 = PairedRegion(3,10,2)
        pr3 = PairedRegion(3,10,2, Id='A')
        pr4 = PairedRegion(3,10,2, Id='A')
        pr5 = PairedRegion(3,20,4, Id='A')
        self.assertEqual(pr1==pr2, True) # same pairs, no IDs
        self.assertEqual(pr3==pr4, True) # same pairs, same IDs
        self.assertEqual(pr1==pr3, False) # same pairs, diff ID
        self.assertEqual(pr3==pr5, False) # diff pairs, same IDs

    def test_upstream(self):
        """PairedRegions upstream: single and multiple pair(s)"""
        pr = PairedRegion(3,10,2)
        self.assertEqual(pr.upstream(), [3,4])
        pr = PairedRegion(3,10,1)
        self.assertEqual(pr.upstream(), [3])

    def test_downstream(self):
        """PairedRegions downstream: single and multiple pair(s)"""
        pr = PairedRegion(3,10,2)
        self.assertEqual(pr.downstream(), [9,10])
        pr = PairedRegion(3,10,1)
        self.assertEqual(pr.downstream(), [10])

    def test_paired(self):
        """PairedRegion paired: single and multiple pair(s)"""
        pr = PairedRegion(3,10,2)
        self.assertEqual(pr.paired(), [3,4,9,10])
        pr = PairedRegion(3,10,1)
        self.assertEqual(pr.paired(), [3,10])

    def test_regionRange(self):
        """PairedRegion regionRange: single and multiple pair(s)"""
        pr = PairedRegion(3,10,2)
        self.assertEqual(pr.range(), 4) 
        pr = PairedRegion(1,10,4)
        self.assertEqual(pr.range(), 2)
        pr = PairedRegion(1,5,1)
        self.assertEqual(pr.range(), 3)

        # no error checking
        pr = PairedRegion(5,8,3) # 5,6,7,-- 6,7,8
        self.assertEqual(pr.range(), -2)

    def test_overlapping(self):
        """PairedRegion overlapping: identical and different regions"""
        pr1 = PairedRegion(1,10,2)
        pr2 = PairedRegion(3,15,2)
        self.assertEqual(pr1.overlapping(pr2), False)
        self.assertEqual(pr2.overlapping(pr1), False)

        pr1 = PairedRegion(1,10,2)
        pr2 = PairedRegion(2,15,2)
        pr3 = PairedRegion(9,20,4)
        self.assertEqual(pr1.overlapping(pr2), True)
        self.assertEqual(pr2.overlapping(pr1), True)
        self.assertEqual(pr1.overlapping(pr3), True)

        bl1 = PairedRegion(2,10,1)
        bl2 = PairedRegion(12,20,3)
        self.assertEqual(bl1.overlapping(bl2), False)

        pr1 = PairedRegion(1,10,2, Id='A')
        pr2 = PairedRegion(1,10,2, Id='A')
        pr3 = PairedRegion(1,10,2, Id='B')
        self.assertEqual(pr1.overlapping(pr2), True)
        self.assertEqual(pr1.overlapping(pr3), True) # ignore ID
    
    def test_conflicting(self):
        """PairedRegion conflicting: identical, nested and pseudoknot"""
        bl1 = PairedRegion(2,10,3)
        bl2 = PairedRegion(12,20,3)
        # identical blocks are NOT conflicting...
        self.assertEqual(bl1.conflicting(bl1), False) # identical blocks
        
        self.assertEqual(bl1.conflicting(bl2), False) # one after the other
        self.assertEqual(bl2.conflicting(bl1), False) # one after the other

        bl1 = PairedRegion(1,30,2) #[(1,30),(2,29)]
        bl2 = PairedRegion(14,20,2) #[(14,20),(15,19)]
        self.assertEqual(bl1.conflicting(bl2), False) # one inside the other
        self.assertEqual(bl2.conflicting(bl1), False) # one inside the other
        
        bl1 = PairedRegion(1,10,2) #[(1,10),(2,9)]
        bl2 = PairedRegion(4,15,3) #[(4,15),(5,14),(6,13)]
        self.assertEqual(bl1.conflicting(bl2), True) # pseudoknot
        self.assertEqual(bl2.conflicting(bl1), True) # pseudoknot

    def test_score(self):
        """PairedRegion score: should take arbitrary scoring function"""
        f = lambda x: x.Length # scoring function
        bl1 = PairedRegion(2,10,3)
        bl2 = PairedRegion(12,30,4)
        bl1.score(f) # set Score attribute for bl1
        bl2.score(f) # set Score attribute for bl2
        self.assertEqual(bl1.Score, 3)
        self.assertEqual(bl2.Score, 4)

    def test_PairedRegionFromPairs(self):
        """PairedRegionFromPairs: should handle valid input"""
        p = Pairs([(3,10),(4,9),(5,8)])
        pr = PairedRegionFromPairs(p, Id='A')
        self.assertEqual(pr.Start, 3)
        self.assertEqual(pr.End, 10)
        self.assertEqual(pr.Length, 3)
        self.assertEqual(pr.Id, 'A')
        self.assertEqual(pr.Pairs, [(3,10),(4,9),(5,8)])
    
    def test_PairedRegionFromPairs_invalid(self):
        """PairedRegionFromPairs: conflicts and error checking"""
        p = Pairs([(3,10),(4,9),(4,None)])
        self.assertRaises(ValueError, PairedRegionFromPairs, p)
    
        # no error checking on input pairs...
        p = Pairs([(3,10),(4,9),(6,8)]) # not a real paired region
        pr = PairedRegionFromPairs(p, Id='A')
        self.assertEqual(pr.Start, 3)
        self.assertEqual(pr.End, 10)
        self.assertEqual(pr.Length, 3)
        self.assertEqual(pr.Id, 'A')
        # NOTE: Pairs will be different than input because assumption does
        # not hold
        self.assertEqual(pr.Pairs, [(3,10),(4,9),(5,8)])

        self.assertRaises(ValueError, PairedRegionFromPairs, [])

class PairedRegionsTests(TestCase):
    """Tests for PairedRegions class"""

    def test_init(self):
        """PairedRegions __init__: should accept list of PairedRegion objects
        """
        pr1 = PairedRegion(3,10,2)
        pr2 = PairedRegion(12,20,3)
        obs = PairedRegions([pr1, pr2])
        self.assertEqual(obs[0].Start, 3)
        self.assertEqual(obs[1].Id, None)
        self.assertEqual(obs[1].End, 20)

    def test_init_no_validation(self):
        """PairedRegions __init__: does not perform validation
        """
        # can give any list of arbitrary object as input
        obs = PairedRegions([1,2,3])
        self.assertEqual(obs[0], 1)

    def test_str(self):
        """PairedRegions __str__: full and empty list
        """
        pr1 = PairedRegion(3,10,2)
        pr2 = PairedRegion(12,20,3, Id='A')
        prs = PairedRegions([pr1, pr2])
        self.assertEqual(str(prs), "(None:3,10,2; A:12,20,3;)")

    def test_eq(self):
        """PairedRegions __eq__: with/without IDs, in or out of order
        """
        pr1 = PairedRegion(3,10,2)
        pr2 = PairedRegion(12,20,3, Id='A')
        prs1 = PairedRegions([pr1, pr2])
        pr3 = PairedRegion(3,10,3)
        pr4 = PairedRegion(20,30,3, Id='A')
        prs2 = PairedRegions([pr3, pr4])
        pr5 = PairedRegion(3,10,2)
        pr6 = PairedRegion(12,20,3, Id='A')
        prs3 = PairedRegions([pr5, pr6])
        prs4 = PairedRegions([pr6, pr5])
        self.assertEqual(prs1==prs2, False)
        self.assertEqual(prs1==prs1, True)
        self.assertEqual(prs1==prs3, True)
        self.assertEqual(prs1==prs4, True)

    def test_ne(self):
        """PairedRegions __ne__: with/without IDs, in or out of order
        """
        pr1 = PairedRegion(3,10,2)
        pr2 = PairedRegion(12,20,3, Id='A')
        prs1 = PairedRegions([pr1, pr2])
        pr3 = PairedRegion(3,10,3)
        pr4 = PairedRegion(20,30,3, Id='A')
        prs2 = PairedRegions([pr3, pr4])
        pr5 = PairedRegion(3,10,2)
        pr6 = PairedRegion(12,20,3, Id='A')
        prs3 = PairedRegions([pr5, pr6])
        prs4 = PairedRegions([pr6, pr5])
        self.assertEqual(prs1!=prs2, True)
        self.assertEqual(prs1!=prs1, False)
        self.assertEqual(prs1!=prs3, False)
        self.assertEqual(prs1!=prs4, False)

    def test_byId(self):
        """PairedRegions byId: unique IDs and duplicates
        """
        pr1 = PairedRegion(3,10,2, Id='A')
        pr2 = PairedRegion(12,20,3, Id='B')
        prs1 = PairedRegions([pr1, pr2])
        obs = prs1.byId()
        self.assertEqual(obs['A'], pr1)
        self.assertEqual(obs['B'], pr2)
        self.assertEqual(len(obs), 2)

        pr3 = PairedRegion(3,10,2, Id='A')
        pr4 = PairedRegion(12,20,3, Id='A')
        prs2 = PairedRegions([pr3, pr4])
        self.assertRaises(ValueError, prs2.byId)

        pr3 = PairedRegion(3,10,2)
        pr4 = PairedRegion(12,20,3)
        prs2 = PairedRegions([pr3, pr4])
        self.assertRaises(ValueError, prs2.byId)

        self.assertEqual(PairedRegions().byId(), {})

    def test_numberOfRegions(self):
        """PairedRegions numberOfRegions: full and empty"""
        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(11,20,4)
        prs1 = PairedRegions([pr1,pr2])
        prs2 = PairedRegions([pr1,pr1,pr2,pr2])
        self.assertEqual(prs1.numberOfRegions(), 2)
        self.assertEqual(prs2.numberOfRegions(), 4)
        self.assertEqual(PairedRegions().numberOfRegions(), 0)

    def test_totalLength(self):
        """PairedRegions totalLength: full and empty"""
        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(11,20,4)
        prs1 = PairedRegions([pr1,pr2])
        self.assertEqual(prs1.totalLength(), 6)
        
        self.assertEqual(PairedRegions().totalLength(), 0)

    def test_totalScore(self):
        """PairedRegions totalScore: full, empty, None"""
        pr1 = PairedRegion(2,10,2)
        pr1.Score = 3
        pr2 = PairedRegion(11,20,4)
        pr2.Score = 2
        pr3 = PairedRegion(11,20,4)
        pr3.Score = None
        pr4 = PairedRegion(11,20,4)
        pr4.Score = "abc"
        pr5 = PairedRegion(11,20,4)
        
        prs1 = PairedRegions([pr1,pr2])
        prs2 = PairedRegions([pr1,pr3])
        prs3 = PairedRegions([pr1,pr4])
        prs4 = PairedRegions([pr1,pr5])
        
        self.assertEqual(prs1.totalScore(), 5)
        self.assertRaises(ValueError, prs2.totalScore)
        self.assertRaises(ValueError, prs3.totalScore)
        self.assertRaises(ValueError, prs4.totalScore)

    def test_toPairs(self):
        """PairedRegions toPairs: good data"""
        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(11,20,4)
        prs1 = PairedRegions([pr2,pr1])
        exp = [(2,10),(3,9),(11,20),(12,19),(13,18),(14,17)]
        self.assertEqual(prs1.toPairs(), exp)

        prs2 = PairedRegions([pr1,pr1])
        exp = [(2,10),(2,10),(3,9),(3,9)]
        self.assertEqual(prs2.toPairs(), exp)

        self.assertEqual(PairedRegions().toPairs(), Pairs())

    def test_byStartEnd(self):
        """PairedRegions byStartEnd: unique and duplicate keys"""
        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(11,20,4)
        prs1 = PairedRegions([pr2,pr1])
        exp = {(2,10): pr1, (11,20): pr2}
        self.assertEqual(prs1.byStartEnd(), exp)

        pr3 = PairedRegion(2,10,2, Id='A')
        pr4 = PairedRegion(2,10,3, Id='B')
        prs2 = PairedRegions([pr3,pr4])
        self.assertRaises(ValueError, prs2.byStartEnd)

    def test_lowestStart(self):
        """PairedRegions lowestStart: full and empty object"""
        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(11,20,4)
        pr3 = PairedRegion(2,30,5,Id='A')
        prs1 = PairedRegions([pr2,pr1,pr3])
        self.assertEqual(prs1.lowestStart(), 2)

        self.assertEqual(PairedRegions().lowestStart(), None)

    def test_highestEnd(self):
        """PairedRegions highestEnd: full and empty object"""
        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(11,20,4)
        pr3 = PairedRegion(2,30,5,Id='A')
        prs1 = PairedRegions([pr2,pr1,pr3])
        self.assertEqual(prs1.highestEnd(), 30)

        self.assertEqual(PairedRegions().highestEnd(), None)

    def test_sortedIds(self):
        """PairedRegions sortedIds: full and empty list"""
        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(2,30,5,Id='B')
        prs1 = PairedRegions([pr1,pr2,pr3])
        self.assertEqual(prs1.sortedIds(), ['A','B','C'])

        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(2,30,5,Id='C')
        prs1 = PairedRegions([pr1,pr2,pr3])
        self.assertEqual(prs1.sortedIds(), ['A','C','C'])

        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(2,30,5,Id=2)
        prs1 = PairedRegions([pr1,pr2,pr3])
        self.assertEqual(prs1.sortedIds(), [None, 2, 'A'])

    def test_upstream(self):
        """PairedRegions upstream: full and empty"""
        self.assertEqual(PairedRegions().upstream(), [])
        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(4,30,1,Id='B')
        prs1 = PairedRegions([pr1,pr2,pr3])
        exp = [2,3,11,12,13,14,4]
        exp.sort()
        self.assertEqual(prs1.upstream(), exp)

    def test_downstream(self):
        """PairedRegions upstream: full and empty"""
        self.assertEqual(PairedRegions().downstream(), [])
        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(4,30,1,Id='B')
        prs1 = PairedRegions([pr1,pr2,pr3])
        exp = [10,9,20,19,18,17,30]
        exp.sort()
        self.assertEqual(prs1.downstream(), exp)

    def test_pairedPos(self):
        """PairedRegions pairedPos: full and empty"""
        self.assertEqual(PairedRegions().pairedPos(), [])
        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(4,30,1,Id='B')
        prs1 = PairedRegions([pr1,pr2,pr3])
        exp = [2,3,11,12,13,14,4,10,9,20,19,18,17,30]
        exp.sort()
        self.assertEqual(prs1.pairedPos(), exp)

    def test_boundaries(self):
        """PairedRegions boundaries: full and empty"""
        self.assertEqual(PairedRegions().boundaries(), [])
        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(4,30,1,Id='B')
        prs1 = PairedRegions([pr1,pr2,pr3])
        exp = [2,10,11,20,4,30]
        exp.sort()
        self.assertEqual(prs1.boundaries(), exp)

    def test_enumeratedBoundaries(self):
        """PairedRegions enumeratedBoundaries: full and empty"""
        self.assertEqual(PairedRegions().enumeratedBoundaries(), {})
        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(4,30,1,Id='B')
        prs1 = PairedRegions([pr1,pr2,pr3])
        exp = {0:2,2:10,3:11,4:20,1:4,5:30}
        self.assertEqual(prs1.enumeratedBoundaries(), exp)

    def test_invertedEnumeratedBoundaries(self):
        """PairedRegions invertedEnumeratedBoundaries: full and empty"""
        self.assertEqual(PairedRegions().invertedEnumeratedBoundaries(), {})
        pr1 = PairedRegion(2,10,2, Id='C')
        pr2 = PairedRegion(11,20,4, Id='A')
        pr3 = PairedRegion(4,30,1,Id='B')
        prs1 = PairedRegions([pr1,pr2,pr3])
        exp = {2:0,10:2,11:3,20:4,4:1,30:5}
        self.assertEqual(prs1.invertedEnumeratedBoundaries(), exp)

        pr1 = PairedRegion(3,10,2)
        pr2 = PairedRegion(5,10,3)
        prs = PairedRegions([pr1, pr2])
        self.assertRaises(ValueError, prs.invertedEnumeratedBoundaries)

    def test_merge(self):
        """PairedRegions merge: different, duplicates, empty"""
        pr1 = PairedRegion(3,10,2, Id='A')
        pr2 = PairedRegion(11,20,3, Id='B')
        pr3 = PairedRegion(15,25,1, Id='C')
        prs1 = PairedRegions([pr1, pr2])
        prs2 = PairedRegions([pr1, pr3])
        prs3 = PairedRegions()

        exp = PairedRegions([pr1,pr2,pr3])
        self.assertEqual(prs1.merge(prs2), exp)
        self.assertEqual(prs2.merge(prs1), exp)
        self.assertEqual(prs1.merge(prs3), prs1)
        self.assertEqual(prs2.merge(prs3), prs2)

    def test_conflicting_no_ids(self):
        """PairedRegions conflicting: raises error on duplicate IDs
        """
        pr1 = PairedRegion(1,10,2)
        pr2 = PairedRegion(11,20,2)
        prs = PairedRegions([pr1, pr2])
        self.assertRaises(ValueError, prs.conflicting) #conflicting IDs

    def test_conflicting(self):
        """PairedRegions conflicting: works when IDs are set and unique
        """
        pr1 = PairedRegion(3,10,2, Id='A')
        pr2 = PairedRegion(11,20,3, Id='B')
        pr3 = PairedRegion(15,25,1, Id='C')
        prs = PairedRegions([pr1,pr2,pr3])
        self.assertEqual(prs.conflicting(), PairedRegions([pr2,pr3]))

        prs = PairedRegions()
        self.assertEqual(prs.conflicting(), PairedRegions())

    def test_non_conflicting_no_ids(self):
        """PairedRegions nonConflicting: raises error on duplicate IDs
        """
        pr1 = PairedRegion(1,10,2)
        pr2 = PairedRegion(11,20,2)
        prs = PairedRegions([pr1, pr2])
        self.assertRaises(ValueError, prs.nonConflicting) #conflicting IDs

    def test_non_conflicting(self):
        """PairedRegions nonConflicting: works when IDs are set and unique
        """
        pr1 = PairedRegion(3,10,2, Id='A')
        pr2 = PairedRegion(11,20,3, Id='B')
        pr3 = PairedRegion(15,25,1, Id='C')
        prs = PairedRegions([pr1,pr2,pr3])
        self.assertEqual(prs.nonConflicting(), PairedRegions([pr1]))

        prs = PairedRegions()
        self.assertEqual(prs.conflicting(), PairedRegions())

    def test_conflictCliques(self):
        """PairedRegions conflictCliques: should work when IDs are unique"""
        pr1 = PairedRegion(3,10,2, Id='A')
        pr2 = PairedRegion(11,20,3, Id='B')
        pr3 = PairedRegion(15,25,1, Id='C')
        pr4 = PairedRegion(30,40,2, Id='D')
        pr5 = PairedRegion(28,35,1, Id='E')
        prs = PairedRegions([pr1,pr2,pr3,pr4, pr5])
        obs = prs.conflictCliques()
        exp = [PairedRegions([pr2,pr3]),PairedRegions([pr5,pr4])]
        for i in obs:
            self.failUnless(i in exp)
        self.assertEqual(len(obs), len(exp))

        prs = PairedRegions()
        self.assertEqual(prs.conflictCliques(), [])

    def test_PairedRegionsFromPairs(self):
        """PairedRegionsFromPairs: should work on valid input"""
        p = Pairs([(1,10),(2,9),(12,20),(13,19),(14,18)])
        prs = PairedRegionsFromPairs(p)
        self.assertEqual(len(prs), 2)
        self.assertEqual(prs[0].Id, 0)
        self.assertEqual(prs[0].Pairs, [(1,10),(2,9)])
        self.assertEqual(prs[0].Start, 1)
        self.assertEqual(prs[0].End, 10)

        self.assertEqual(PairedRegionsFromPairs(Pairs()), PairedRegions())

    def test_PairedRegionsFromPairs_conflict(self):
        """PairedRegionsFromPairs: should raise error on overlapping pairs"""
        p = Pairs([(2,20),(5,10),(10,15)])
        self.assertRaises(ValueError, PairedRegionsFromPairs, p)

class ConflictMatrixTests(TestCase):
    """Tests for ConflictMatrix class"""

    def test_conflict_matrix_from_pairs(self):
        """ConflixtMatrix __init__: Pairs as input, w/wo conflict """
        
        f = ConflictMatrix
        # conflict free
        d = [(1,10),(2,9),(12,20),(13,19),(14,18)]
        exp = Dict2D({0:{0:False,1:False},1:{0:False,1:False}})
        self.assertEqual(f(d).Matrix, exp)
        self.failIf(not isinstance(f(d).Matrix, Dict2D))
        
        # 1 conflict
        d = [(1,10),(2,9),(12,20),(13,19),(14,18),(15,30),(16,29)]
        exp = Dict2D({0:{0:False,1:False,2:False},\
            1:{0:False,1:False,2:True},\
            2:{0:False,1:True,2:False}})
        self.assertEqual(f(d).Matrix, exp)

        # 1 conflict
        d = Pairs([(1,10),(2,9),(12,20),(13,19),(14,18),(15,30),(16,29)])
        exp = Dict2D({0:{0:False,1:False,2:False},\
            1:{0:False,1:False,2:True},\
            2:{0:False,1:True,2:False}})
        m = f(d).Matrix
        self.assertEqual(m, exp)
        self.assertEqual(m.RowOrder, [0,1,2])
        self.assertEqual(m.ColOrder, [0,1,2])
        
        d = [] # empty input
        exp = Dict2D()
        self.assertEqual(f(d).Matrix, exp)

    def test_ConflictMatrix_Pairs_overlap(self):
        """ConflictMatrix __init__: raises error on overlapping pairs"""
        p = Pairs([(1,10),(2,9),(3,9),(12,20)])
        self.assertRaises(ValueError, ConflictMatrix, p)

    def test_conflict_matrix_from_PairedRegions(self):
        """ConflictMatrix __init__: PairedRegions as input, w/wo conflict
        """
        f = ConflictMatrix
        # conflict free
        pr1 = PairedRegion(1,10,2, Id=0)
        pr2 = PairedRegion(12,20,3, Id=1)
        prs = PairedRegions([pr1,pr2])
        exp = Dict2D({0:{0:False,1:False},1:{0:False,1:False}})
        self.assertEqual(f(prs).Matrix, exp)
        self.failIf(not isinstance(f(prs).Matrix, Dict2D))
        
        pr1 = PairedRegion(1,10,2, Id=0)
        pr2 = PairedRegion(12,20,3, Id=1)
        pr3 = PairedRegion(15,30,2, Id=2)
        prs = PairedRegions([pr1,pr2, pr3])
        # 1 conflict
        exp = Dict2D({0:{0:False,1:False,2:False},\
            1:{0:False,1:False,2:True},\
            2:{0:False,1:True,2:False}})
        self.assertEqual(f(prs).Matrix, exp)
        
        # 1 conflict
        pr1 = PairedRegion(1,10,2, Id=4)
        pr2 = PairedRegion(12,20,3, Id=1)
        pr3 = PairedRegion(15,30,2, Id=9)
        prs = PairedRegions([pr1,pr2, pr3])
        exp = Dict2D({1:{4:False,1:False,9:True},\
            4:{1:False,4:False,9:False},\
            9:{1:True,4:False,9:False}})
        m = f(prs).Matrix
        self.assertEqual(m, exp)
        self.assertEqual(m.RowOrder, [1,4,9])
        self.assertEqual(m.ColOrder, [1,4,9])
        
        prs = PairedRegions()
        exp = Dict2D()
        self.assertEqual(f(prs).Matrix, exp)
        # input some weird data. Other errors might occur.
        self.assertRaises(ValueError, f, 'ABC')
        self.assertRaises(ValueError, f, [('a','b'),('c','d')])

    def test_ConflictMatrix_PairedRegions_overlap(self):
        """ConflictMatrix __init__: raises error on overlapping PairedRegions
        """
        pr1 = PairedRegion(1,10,2, Id='A')
        pr2 = PairedRegion(8,20,2, Id='B')
        prs = PairedRegions([pr1, pr2])
        self.assertRaises(ValueError, ConflictMatrix, prs)
        
    def test_conflictsOf(self):
        """ConflictMatrix conflictsOf: with/without conflicts"""
        p = Pairs([(1,10),(5,15),(20,30),(25,35),(24,32),(0,80)])
        cm = ConflictMatrix(p)
        self.assertEqual(cm.conflictsOf(0), [])
        self.assertEqual(cm.conflictsOf(1), [2])
        self.assertEqual(cm.conflictsOf(2), [1])
        self.assertEqual(cm.conflictsOf(3), [4,5])

        p = Pairs([(1,10),(11,20)])
        cm = ConflictMatrix(p)
        self.assertEqual(cm.conflictsOf(0), [])
        self.assertEqual(cm.conflictsOf(1), [])
        self.assertRaises(KeyError, cm.conflictsOf, 2)
        
    def test_conflicting(self):
        """ConflictMatrix conflicting: full and empty Pairs"""
        p = Pairs([(1,10),(5,15),(20,30),(25,35),(24,32),(0,80)])
        cm = ConflictMatrix(p)
        obs = cm.conflicting()
        exp = [1,2,3,4,5]
        self.assertEqual(obs, exp)

        self.assertEqual(ConflictMatrix(Pairs()).conflicting(), [])

    def test_nonConflicting(self):
        """ConflictMatrix nonConflicting: full and empty Pairs"""
        p = Pairs([(1,10),(5,15),(20,30),(25,35),(24,32),(0,80)])
        cm = ConflictMatrix(p)
        obs = cm.nonConflicting()
        exp = [0]
        self.assertEqual(obs, exp)

        self.assertEqual(ConflictMatrix(Pairs()).nonConflicting(), [])

    def test_conflictCliques(self):
        """ConflictMatrix conflictCliques: full and empty Pairs"""
        p = Pairs([(1,10),(5,15),(20,30),(25,35),(24,32),(0,80)])
        cm = ConflictMatrix(p)
        obs = cm.conflictCliques()
        exp = [[1,2],[3,4,5]]
        self.assertEqual(obs, exp)

        self.assertEqual(ConflictMatrix(Pairs()).conflictCliques(), [])
    
    
class DPTests(TestCase):
    """Tests for opt_all and related functions"""

    def test_num_bps(self):
        """num_bps: should return length of paired region"""
        f = num_bps
        pr1 = PairedRegion(0,10,3)
        self.assertEqual(f(pr1), 3)

    def test_hydrogen_bonds(self):
        """hydrogen_bonds: score GC, AU, and GU base pairs"""
        f = hydrogen_bonds('UACGAAAUGCGUG')
        pr1 = PairedRegion(0,12,5)
        self.assertEqual(f(pr1),10)

        f = hydrogen_bonds('UACGAAA') # sequence too short
        pr1 = PairedRegion(0,12,5)
        self.assertRaises(IndexError, f, pr1)
    
    def test_contains_true(self):
        """contains_true: should return True if True in input"""
        f = contains_true
        self.assertEqual(f([True]), True)
        self.assertEqual(f([True, False]), True)
        self.assertEqual(f([1, 0]), True)
        self.assertEqual(f([1]), True)
        self.assertEqual(f([False]), False)
        self.assertEqual(f([3]), False)
        self.assertEqual(f(["a","b","c"]), False)
        self.assertEqual(f("abc"), False)

    def test_empty_matrix(self):
        """empty_matrix: valid input and error"""
        f = empty_matrix
        p = PairedRegions()
        exp = [[[p],[p]], [[p],[p]]]
        self.assertEqual(f(2), exp)
        
        self.assertEqual(f(1), [[[p]]])

        self.assertRaises(ValueError, f, 0)

    def test_pick_multi_best_max(self):
        """pick_multi_best: max, full and empty list"""
        pr1 = PairedRegion(2,10,2, Id='A')
        pr2 = PairedRegion(4,15,3, Id='B')
        pr3 = PairedRegion(20,40,5, Id='C')
        pr4 = PairedRegion(22,30,3, Id='D')
        for i in [pr1,pr2,pr3,pr4]:
            i.score(num_bps)
        prs1 = PairedRegions([pr1, pr2])
        prs2 = PairedRegions([pr3])
        prs3 = PairedRegions([pr4])
        self.assertEqualItems(pick_multi_best([prs1, prs2, prs3]), [prs1, prs2])

        self.assertEqual(pick_multi_best([]), [PairedRegions()])

    def test_pick_multi_best_min(self):
        """pick_multi_best: min, full and empty list"""
        f = lambda x: -1
        pr1 = PairedRegion(2,10,2)
        pr2 = PairedRegion(4,15,3)
        pr3 = PairedRegion(20,40,5)
        pr4 = PairedRegion(22,30,3)
        for i in [pr1,pr2,pr3,pr4]:
            i.score(f)
        prs1 = PairedRegions([pr1, pr2])
        prs2 = PairedRegions([pr3])
        prs3 = PairedRegions([pr4])
        self.assertEqual(pick_multi_best([prs1, prs2, prs3], goal='min'),\
            [prs1])

        self.assertEqual(pick_multi_best([], goal='min'), [PairedRegions()])

    def test_dp_matrix_multi_toy(self):
        """dp_matrix_multi: test on initial toy example"""
        pr0 = PairedRegion(0, 70, 2, Id='C')
        pr1 = PairedRegion(10, 30, 4, Id='A')
        pr2 = PairedRegion(20, 50, 3, Id='B')
        pr3 = PairedRegion(40, 90, 2, Id='E')
        pr4 = PairedRegion(60, 80, 3, Id='D')
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4])  
        obs = dp_matrix_multi(prs)
        self.assertEqual(obs[0][0], [PairedRegions()])
        self.assertEqual(obs[0][3], [PairedRegions([pr1])])
        self.assertEqual(obs[2][5], [PairedRegions([pr2])])
        self.assertEqual(obs[4][9], [PairedRegions([pr3,pr4])])
        self.assertEqual(obs[2][9], [PairedRegions([pr2,pr4])])
        self.assertEqual(obs[1][8], [PairedRegions([pr1,pr4])])
        self.assertEqual(obs[1][9], [PairedRegions([pr1,pr3,pr4])])
        self.assertEqual(obs[0][9], [PairedRegions([pr1,pr3,pr4])])

    def test_dp_matrix_multi_lsu(self):
        """dp_matrix_multi: test on LSU rRNA domain I case"""

        pr0 = PairedRegion(56, 69, 3, Id=0)
        pr1 = PairedRegion(60, 92, 1, Id=1)
        pr2 = PairedRegion(62, 89, 3, Id=2)
        pr3 = PairedRegion(75, 109, 6, Id=3)
        pr4 = PairedRegion(84, 96, 3, Id=4)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4])
        obs = dp_matrix_multi(prs)
        self.assertEqual(obs[0][0], [PairedRegions()])
        self.assertEqual(obs[0][5], [PairedRegions([pr0])])
        self.assertEqual(obs[1][6], [PairedRegions([pr2])])
        self.assertEqual(obs[1][7], [PairedRegions([pr1,pr2])])
        self.assertEqualItems(obs[2][8],\
            [PairedRegions([pr2]),PairedRegions([pr4])])
        self.assertEqual(obs[1][9], [PairedRegions([pr3,pr4])])
        self.assertEqual(obs[0][9], [PairedRegions([pr0,pr3,pr4])])

    def test_dp_matrix_multi_artificial(self):
        """dp_matrix_multi: test on artificial structure"""

        pr0 = PairedRegion(0, 77, 2, Id=0)
        pr1 = PairedRegion(7, 75, 5, Id=1)
        pr2 = PairedRegion(13, 83, 3, Id=2)
        pr3 = PairedRegion(18, 41, 5, Id=3)
        pr4 = PairedRegion(23, 53, 10, Id=4)
        pr5 = PairedRegion(33, 70, 3, Id=5)
        pr6 = PairedRegion(59, 93, 9, Id=6)
        pr7 = PairedRegion(78, 96, 3, Id=7)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr6, pr7])
        obs = dp_matrix_multi(prs)
        self.assertEqual(obs[0][0], [PairedRegions()])
        self.assertEqual(obs[0][6], [PairedRegions([pr3])])
        self.assertEqual(obs[0][7], [PairedRegions([pr4])])
        self.assertEqual(obs[9][15], [PairedRegions([pr7])])
        self.assertEqual(obs[1][10], [PairedRegions([pr1,pr4])])
        self.assertEqual(obs[0][11], [PairedRegions([pr0,pr1,pr4])])
        self.assertEqual(obs[3][14], [PairedRegions([pr4, pr6])])
        self.assertEqual(obs[3][14], [PairedRegions([pr4, pr6])])
        self.assertEqual(obs[0][13], [PairedRegions([pr0, pr1, pr4])])
        self.assertEqual(obs[0][14], [PairedRegions([pr4, pr6])])
        self.assertEqual(obs[1][15], [PairedRegions([pr4, pr6])])
        self.assertEqual(obs[0][15], [PairedRegions([pr0, pr1, pr4, pr7])])

    def test_pick_multi_best_saturated(self):
        """pick_multi_best: should only include saturated solutions"""
        pr1 = PairedRegion(2,10,2, Id='A')
        pr1.Score = 2
        pr2 = PairedRegion(15,25,2, Id='B')
        pr2.Score = 2
        pr3 = PairedRegion(4,22,4, Id='C')
        pr3.Score = 0
        prs1 = PairedRegions([pr1])
        prs2 = PairedRegions([pr2])
        prs3 = PairedRegions([pr1, pr3])
        self.assertEqualItems(pick_multi_best([prs1, prs2, prs3]),\
            [prs2, prs3])
        self.assertEqual(pick_multi_best([]), [PairedRegions()])

    def test_matrix_solutions(self):
        """matrix_solutions: should return contents of top-right cell"""
        pr0 = PairedRegion(56, 69, 3, Id=0)
        pr1 = PairedRegion(60, 92, 1, Id=1)
        pr2 = PairedRegion(62, 89, 3, Id=2)
        pr3 = PairedRegion(75, 109, 6, Id=3)
        pr4 = PairedRegion(84, 96, 3, Id=4)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4])
        obs = matrix_solutions(prs)
        self.assertEqual(obs, [PairedRegions([pr0,pr3,pr4])])

        # error, size should be at least 1
        prs = PairedRegions()
        self.assertRaises(ValueError, matrix_solutions, prs)

        pr = PairedRegion(2,20, 5, Id='A')
        prs = PairedRegions([pr])
        obs = matrix_solutions(prs)
        self.assertEqual(obs, [prs])

    def test_opt_all_nested(self):
        """opt_all: should return input when already nested"""
        p = Pairs([(1,10),(2,9),(20,30),(22,29)])
        obs = opt_all(p)
        self.assertEqual(len(obs),1)
        self.assertEqual(obs[0], p)

        p = Pairs()
        self.assertEqual(opt_all(p), [[]])

    def test_opt_all_overlap(self):
        """opt_all: should raise error on overlapping pairs"""
        p = Pairs([(1,10),(2,9),(9,30),(22,29),(1,None)])
        self.assertRaises(ValueError, opt_all, p)

    def test_opt_all_knot(self):
        """opt_all: single/multiple solution(s)"""
        p = Pairs([(1,10),(2,9),(3,15),(4,14),(11,20),(12,19),(25,30)])
        obs = opt_all(p)
        exp = Pairs([(1,10),(2,9),(11,20),(12,19),(25,30)])
        exp_rem = [(3,15),(4,14)]
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0], exp)
        self.assertEqual(opt_all(p, return_removed=True)[0][1],\
            exp_rem)
        
        p = Pairs([(1,10),(2,9),(4,14),(3,15)])
        obs = opt_all(p)
        self.assertEqual(len(obs), 2)
        self.assertEqualItems(obs, [Pairs([(1,10),(2,9)]),\
            Pairs([(3,15),(4,14)])])
        exp_rem = [(Pairs([(1,10),(2,9)]),Pairs([(3,15),(4,14)])),\
            (Pairs([(3,15),(4,14)]),Pairs([(1,10),(2,9)]))]
        self.assertEqualItems(opt_all(p, return_removed=True),\
            exp_rem)

    def test_opt_all_some_non_conflicting(self):
        """opt_all: some conflicting, other not"""
        p = Pairs([(30,40),(10,20),(12,17),(13,None),(17,12),(35,45),(36,44)])
        exp = Pairs([(10,20),(12,17),(35,45),(36,44)])
        exp_rem = [(30,40)]
        self.assertEqual(opt_all(p, return_removed=True),\
            [(exp,exp_rem)])

    def test_opt_all_scoring1(self):
        """opt_all: one optimal in bps, both optimal in energy"""
        p = Pairs([(1,10),(2,9),(4,15),(5,14),(6,13)])
        obs_bps = opt_all(p, goal='max', scoring_function=num_bps)
        obs_energy = opt_all(p, goal='max',\
            scoring_function=hydrogen_bonds('CCCAAAUGGGGUCGUUC'))
        exp_bps = [[(4,15),(5,14),(6,13)]]
        exp_energy = [[(1,10),(2,9)],[(4,15),(5,14),(6,13)]]
        self.assertEqualItems(obs_bps, exp_bps)
        self.assertEqualItems(obs_energy, exp_energy)

    def test_opt_all_scoring2(self):
        """opt_all: both optimal in bps, one optimal in energy"""
        p = Pairs([(0,9),(1,8),(2,7),(3,13),(4,12),(5,11)])
        obs_bps = opt_all(p, goal='max', scoring_function=num_bps)
        obs_energy = opt_all(p, goal='max',\
            scoring_function=hydrogen_bonds('CCCAAAAGGGUUUU'))
        exp_bps = [[(0,9),(1,8),(2,7)],[(3,13),(4,12),(5,11)]]
        exp_energy = [[(0,9),(1,8),(2,7)]]
        self.assertEqualItems(obs_bps, exp_bps)
        self.assertEqualItems(obs_energy, exp_energy)

    def test_opt_all_scoring3(self):
        """opt_all: one optimal in bps, the other optimal in energy"""
        p = Pairs([(0,11),(1,10),(2,9),(4,15),(5,14),(6,13),(7,12)])
        obs_bps = opt_all(p, goal='max', scoring_function=num_bps)
        obs_energy = opt_all(p, goal='max',\
            scoring_function=hydrogen_bonds('CCCCAAAAGGGGUUUU'))
        exp_bps = [[(4,15),(5,14),(6,13),(7,12)]]
        exp_energy = [[(0,11),(1,10),(2,9)]]
        self.assertEqualItems(obs_bps, exp_bps)
        self.assertEqualItems(obs_energy, exp_energy)
    
    def test_opt_single_random(self):
        """opt_single_random: should return single solution"""
        p = Pairs ([(10,20),(11,19),(15,25),(16,24)])
        exp1, exp_rem1 = [(10,20),(11,19)], [(15,25),(16,24)]
        exp2, exp_rem2 = [(15,25),(16,24)], [(10,20),(11,19)]
        obs = opt_single_random(p)
        self.failUnless(obs == exp1 or obs == exp2)
        obs = opt_single_random(p, return_removed=True)
        self.failUnless(obs == (exp1, exp_rem1) or obs == (exp2, exp_rem2))
           
    def test_opt_single_property(self):
        """opt_single_property: three properties"""
        # one solution single region, other solution two regions
        p = Pairs ([(10,20),(25,35),(26,34),(27,33),\
            (12,31),(13,30),(14,29),(15,28)])
        exp = [(12,31),(13,30),(14,29),(15,28)]
        exp_rem = [(10,20),(25,35),(26,34),(27,33)]
        self.assertEqual(opt_single_property(p), exp)
        self.assertEqual(opt_single_property(p, return_removed=True),\
            (exp,exp_rem))

        # both two blocks, one shorter average range
        p = Pairs ([(10,20),(22,40),(23,39),(24,38),\
            (17,26),(18,25),(36,43),(37,42)])
        exp = [(17,26),(18,25),(36,43),(37,42)]
        exp_rem = [(10,20),(22,40),(23,39),(24,38)]
        self.assertEqual(opt_single_property(p), exp)
        self.assertEqual(opt_single_property(p, return_removed=True),\
            (exp,exp_rem))

        # both single block over same range, pick lowest start
        p = Pairs([(10,20),(15,25)])
        exp = [(10,20)]
        exp_rem = [(15,25)]
        self.assertEqual(opt_single_property(p), exp)
        self.assertEqual(opt_single_property(p, return_removed=True),\
            (exp,exp_rem))

class EliminationMethodsTests(TestCase):
    """Tests for conflict_elimination and related functions"""
        
    def test_find_max_conflicts(self):
        """find_max_conflicts: simple case"""
        f = find_max_conflicts
        pr0 = PairedRegion(0, 77, 2, Id=0)
        pr1 = PairedRegion(7, 75, 5, Id=1)
        pr2 = PairedRegion(13, 83, 3, Id=2)
        pr3 = PairedRegion(18, 41, 5, Id=3)
        pr4 = PairedRegion(23, 53, 10, Id=4)
        pr5 = PairedRegion(33, 70, 3, Id=5)
        pr6 = PairedRegion(59, 93, 9, Id=6)
        pr7 = PairedRegion(78, 96, 3, Id=7)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr6, pr7])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 6)

        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr7])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 2)

        prs = PairedRegions([pr0, pr1, pr3, pr4, pr5, pr7])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 5)

    def test_find_max_conflicts_on_start(self):
        """find_max_conflicts: in case of equal conflicts and gain"""
        f = find_max_conflicts
        pr0 = PairedRegion(10, 20, 2, Id=0)
        pr1 = PairedRegion(15, 25, 2, Id=1)
        prs = PairedRegions([pr0, pr1])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 1)

    def test_find_min_gain(self):
        """find_min_gain: differentiate on gain only"""
        f = find_min_gain
        pr0 = PairedRegion(0, 77, 2, Id=0)
        pr1 = PairedRegion(7, 75, 5, Id=1)
        pr2 = PairedRegion(13, 83, 3, Id=2)
        pr3 = PairedRegion(18, 41, 5, Id=3)
        pr4 = PairedRegion(23, 53, 10, Id=4)
        pr5 = PairedRegion(33, 70, 3, Id=5)
        pr6 = PairedRegion(59, 93, 9, Id=6)
        pr7 = PairedRegion(78, 96, 3, Id=7)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr6, pr7])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 5)

        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr6, pr7])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 2)

    def test_find_min_gain_conf(self):
        """find_min_gain: in case of equal gain, differentiate on conflicts"""
        f = find_min_gain
        pr0 = PairedRegion(10,30,3, Id=0)
        pr1 = PairedRegion(1,20,6, Id=1)
        pr2 = PairedRegion(22,40,2, Id=2)
        pr3 = PairedRegion(50,80,3, Id=3)
        pr4 = PairedRegion(60,90,8, Id=4)

        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 0)

    def test_find_min_gain_start(self):
        """find_min_gain: in case of equal gain and number of conflicts"""
        f = find_min_gain
        pr0 = PairedRegion(10,30,3, Id=0)
        pr1 = PairedRegion(1,20,6, Id=1)
        pr2 = PairedRegion(22,40,2, Id=2)
        pr3 = PairedRegion(50,80,3, Id=3)
        pr4 = PairedRegion(60,90,7, Id=4)
        pr5 = PairedRegion(45,55,1, Id=5)

        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5])
        id_to_pr = prs.byId()
        cm = ConflictMatrix(prs)
        conf = cm.conflicting()
        self.assertEqual(f(conf, cm, prs.byId()), 3)

    def test_add_back_non_conflicting(self):
        """add_back_non_conflicting: should add all non-confl regions"""
        f = add_back_non_conflicting
        pr0 = PairedRegion(10,20,3, Id=0)
        pr1 = PairedRegion(30,40,2, Id=1)
        pr2 = PairedRegion(50,60,2, Id=2)
        pr3 = PairedRegion(45,55,3, Id=3) # confl with pr1 and pr2
        pr4 = PairedRegion(0,90,7, Id=4) # not confl with 1,2,3
        pr5 = PairedRegion(32,38,2, Id=5) # not confl with 1,2,3
        prs = PairedRegions([pr0, pr1, pr2])
        removed = {3: pr3, 4: pr4, 5: pr5}
        exp_prs = PairedRegions([pr0, pr1, pr2, pr4, pr5])
        exp_rem = {3: pr3}
        self.assertEqual(f(prs, removed), (exp_prs, exp_rem))

    def test_add_back_non_conflicting_order(self):
        """add_back_non_conflicting: should add 5' side first"""
        f = add_back_non_conflicting
        pr0 = PairedRegion(10,20,3, Id=0)
        pr1 = PairedRegion(30,40,2, Id=1)
        pr2 = PairedRegion(50,60,2, Id=2)
        pr3 = PairedRegion(45,55,3, Id=3) # confl with pr1 and pr2
        pr4 = PairedRegion(0,90,7, Id=4) # not confl with 1,2,3
        pr5 = PairedRegion(80,95,2, Id=5) # not confl with 1,2,3
        prs = PairedRegions([pr0, pr1, pr2])
        removed = {3: pr3, 4: pr4, 5: pr5}
        exp_prs = PairedRegions([pr0, pr1, pr2, pr4 ])
        exp_rem = {3: pr3, 5: pr5}
        self.assertEqual(f(prs, removed), (exp_prs, exp_rem))

    def test_elim_most_conflict(self):
        """conflict_elimination: find_max_conflicts, simple case"""
        f = conflict_elimination
        func = find_max_conflicts
        pr0 = PairedRegion(0, 77, 2, Id=0)
        pr1 = PairedRegion(7, 75, 5, Id=1)
        pr2 = PairedRegion(13, 83, 3, Id=2)
        pr3 = PairedRegion(18, 41, 5, Id=3)
        pr4 = PairedRegion(23, 53, 10, Id=4)
        pr5 = PairedRegion(33, 70, 3, Id=5)
        pr6 = PairedRegion(59, 93, 9, Id=6)
        pr7 = PairedRegion(78, 96, 3, Id=7)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr6, pr7]) 
        pairs = prs.toPairs()
        exp = PairedRegions([pr0, pr1, pr4, pr7]).toPairs()
        exp_rem = PairedRegions([pr2, pr3, pr5, pr6]).toPairs()
        self.assertEqual(f(pairs, func), exp)
        self.assertEqual(f(pairs, func, return_removed=True), (exp, exp_rem))

    def test_elim_mc_circular(self):
        """conflict_elimination: find_max_conflicts, circular removal"""
        # simply remove in order of most conflicts, don't add back
        prfp = PairedRegionFromPairs
        f = conflict_elimination
        func = find_max_conflicts
        pr0 = prfp([(13, 65), (14, 64)], Id=0)
        pr1 = prfp([(15, 102), (16, 101), (17, 100), (18, 99), (19, 98)], Id=1)
        pr2 = prfp([(22, 72), (23, 71), (24, 70), (25, 69),\
            (26, 68), (27, 67), (28, 66)], Id=2)
        pr3 = prfp([(31, 147), (32, 146), (33, 145), (34, 144), (35, 143),\
            (36, 142), (37, 141), (38, 140), (39, 139)], Id=3)
        pr4 = prfp([(42, 129), (43, 128), (44, 127)], Id=4)
        pr5 = prfp([(46, 149), (47, 148)], Id=5)
        pr6 = prfp([(49, 92), (50, 91), (51, 90), (52, 89), (53, 88)], Id=6)
        pr7 = prfp([(75, 138), (76, 137), (77, 136), (78, 135)], Id=7)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr6, pr7])
        exp = PairedRegions([pr3, pr6]).toPairs()
        exp_rem = PairedRegions([pr0, pr1, pr2, pr4, pr5, pr7]).toPairs()
        self.assertEqual(f(prs.toPairs(), func, add_back=False,\
            return_removed=True), (exp, exp_rem))

        # add back circular removals
        exp = PairedRegions([pr3, pr4, pr6]).toPairs()
        exp_rem = PairedRegions([pr0, pr1, pr2, pr5, pr7]).toPairs()
        self.assertEqual(f(prs.toPairs(), func, add_back=True,\
            return_removed=True), (exp, exp_rem))

    def test_elim_min_gain(self):
        """conflict_elimination: find_min_gain, simple case"""
        f = conflict_elimination
        func = find_min_gain
        pr0 = PairedRegion(0, 77, 2, Id=0)
        pr1 = PairedRegion(7, 75, 5, Id=1)
        pr2 = PairedRegion(13, 83, 3, Id=2)
        pr3 = PairedRegion(18, 41, 5, Id=3)
        pr4 = PairedRegion(23, 53, 10, Id=4)
        pr5 = PairedRegion(33, 70, 3, Id=5)
        pr6 = PairedRegion(59, 93, 9, Id=6)
        pr7 = PairedRegion(78, 96, 3, Id=7)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr6, pr7]) 
        pairs = prs.toPairs()
        exp = PairedRegions([pr4, pr6]).toPairs()
        exp_rem = PairedRegions([pr0, pr1, pr2, pr3, pr5, pr7]).toPairs()
        self.assertEqual(f(pairs, func), exp)
        self.assertEqual(f(pairs, func, return_removed=True), (exp, exp_rem))

    def test_elim_min_gain_circular(self):
        """conflict_elimination: find_min_gain, circular removal"""
        # simply remove in order of most conflicts, don't add back
        prfp = PairedRegionFromPairs
        f = conflict_elimination
        func = find_min_gain
        pr0 = prfp([(5, 170), (6, 169), (7, 168), (8, 167), (9, 166),\
            (10, 165)], Id=0)
        pr1 = prfp([(25, 62), (26, 61)], Id=1)
        pr2 = prfp([(29, 46), (30, 45), (31, 44)], Id=2)
        pr3 = prfp([(48, 124), (49, 123)], Id=3)
        pr4 = prfp([(67, 183), (68, 182), (69, 181), (70, 180), (71, 179),\
            (72, 178), (73, 177), (74, 176), (75, 175), (76, 174)], Id=4)
        pr5 = prfp([(82, 172), (83, 171)], Id=5)
        pr6 = prfp([(117, 135), (118, 134), (119, 133)], Id=6)
        pr7 = prfp([(151, 199), (152, 198), (153, 197), (154, 196),\
            (155, 195), (156, 194), (157, 193), (158, 192), (159, 191),\
            (160, 190)], Id=7)
        prs = PairedRegions([pr0, pr1, pr2, pr3, pr4, pr5, pr6, pr7])
        exp = PairedRegions([pr1, pr2, pr4, pr6]).toPairs()
        exp_rem = PairedRegions([pr0, pr3, pr5, pr7]).toPairs()
        self.assertEqual(f(prs.toPairs(), func, add_back=False,\
            return_removed=True), (exp, exp_rem))

        # add back circular removals
        exp = PairedRegions([pr1, pr2, pr4, pr5, pr6]).toPairs()
        exp_rem = PairedRegions([pr0, pr3, pr7]).toPairs()
        self.assertEqual(f(prs.toPairs(), func, add_back=True,\
            return_removed=True), (exp, exp_rem))
 

class IncrementalMethodsTests(TestCase):
    """Tests for incremental pseudoknot-removal methods"""
    
    def test_inc_order_forward(self):
        """nested_in_order: starting at 5' end"""
        f = inc_order
        p = Pairs([(1,10),(2,9),(3,15),(4,14),(11,20),(12,19),(25,30)])
        exp = Pairs([(1,10),(2,9),(11,20),(12,19),(25,30)])
        exp_rem = Pairs([(3,15),(4,14)])
        self.assertEqual(f(p,reversed=False), exp)
        self.assertEqual(f(p,return_removed=True), (exp, exp_rem))

        p = Pairs([(1,20),(2,30),(3,29),(4,28),(5,27),(7,24)])
        exp = Pairs([(1,20)])
        exp_rem = Pairs([(2,30),(3,29),(4,28),(5,27),(7,24)])
        self.assertEqual(f(p,reversed=False), exp)
        self.assertEqual(f(p,return_removed=True), (exp, exp_rem))

        self.assertEqual(f([]), [])

        p = [(1,10),(3,13)] # input as list of tuples
        exp = Pairs([(1,10)])
        self.assertEqual(f(p), exp)

        p = [(1,10),(4,7),(2,9),(5,None)] # pseudoknot-free
        exp = [(1,10),(2,9),(4,7)]
        self.assertEqual(f(p), exp)

        p = [(1,10),(2,10)] #conflict
        self.assertRaises(ValueError, f, p)

    def test_inc_order_reversed(self):
        """nested_in_order: starting at 3' end"""
        f = inc_order
        p = Pairs([(1,10),(2,9),(3,15),(4,14),(24,31),(25,30)])
        exp = Pairs([(3,15),(4,14),(24,31),(25,30)])
        exp_rem = Pairs([(1,10),(2,9)])
        self.assertEqual(f(p,reversed=True), exp)
        self.assertEqual(f(p, reversed=True, return_removed=True),\
            (exp, exp_rem))

        p = Pairs([(1,20),(2,30),(3,29),(4,28),(5,27),(7,24)])
        exp = Pairs([(2,30),(3,29),(4,28),(5,27),(7,24)])
        exp_rem = Pairs([(1,20)])
        self.assertEqual(f(p,reversed=True), exp)
        self.assertEqual(f(p, reversed=True, return_removed=True),\
            (exp, exp_rem))

        self.assertEqual(f([], reversed=True), [])

        p = [(1,10),(3,13)] # input as list of tuples
        exp = Pairs([(3,13)])
        self.assertEqual(f(p, reversed=True), exp)

        p = [(1,10),(4,7),(2,9),(5,None)] # pseudoknot-free
        exp = [(1,10),(2,9),(4,7)]
        self.assertEqual(f(p), exp)

        p = [(1,10),(2,10)] #conflict
        self.assertRaises(ValueError, f, p)

    def test_inc_length(self):
        """inc_length: should handle standard input
        """
        f = inc_length
        # All blocks in conflict, start empty, add first
        p = Pairs([(1,10),(2,9),(3,8),(5,13),(6,12),(7,11)])
        exp = Pairs([(1,10),(2,9),(3,8)])
        self.assertEqual(f(p), exp)

        # Start with length 3 and 2, add 1 block
        p = Pairs([(1,10),(2,9),(3,8),(20,30),(21,29),(25,40),(32,38)])
        exp = Pairs([(1,10),(2,9),(3,8),(20,30),(21,29),(32,38)])
        self.assertEqual(f(p), exp)

        p = Pairs([(1,10),(2,9),(3,8),(12,20),(13,19),(15,23),(16,22)])
        exp_5 = Pairs([(1,10),(2,9),(3,8),(12,20),(13,19)])
        exp_3 = Pairs([(1,10),(2,9),(3,8),(15,23),(16,22)])
        self.assertEqual(f(p), exp_5)
        self.assertEqual(f(p, reversed=True), exp_3)
        self.assertEqual(f(p, return_removed=True),(exp_5,[(15,23),(16,22)]))

        p = [(1,10),(4,7),(2,9),(5,None)] # pseudoknot-free
        exp = [(1,10),(2,9),(4,7)]
        self.assertEqual(f(p), exp)
        
        p = [(1,10),(2,10)] #conflict
        self.assertRaises(ValueError, f, p)

    def test_inc_length_rev(self):
        """inc_length: should prefer 3' side when reversed is True
        """
        f = inc_length
        p = Pairs([(1,10),(2,9),(5,20),(6,19)])
        self.assertEqual(f(p), [(1,10),(2,9)])
        self.assertEqual(f(p, reversed=True), [(5,20),(6,19)])

    def test_inc_range(self):
        """inc_range: should handle normal input
        """
        f = inc_range
        p = [(1,5),(4,20),(15,23),(16,22)]
        exp = [(1,5),(15,23),(16,22)]
        self.assertEqual(f(p), exp)
        self.assertEqual(f(p, return_removed=True), (exp, [(4,20)]))

        p = [(1,11),(5,15)] # same range
        self.assertEqual(f(p), [(1,11)]) # 5' wins
        self.assertEqual(f(p, reversed=True), [(5,15)]) # 3' wins

        p = [(1,10),(2,10)] #conflict
        self.assertRaises(ValueError, f, p)

    def test_inc_range_empty(self):
        """inc_range: should handle empty or pseudoknot-free pairs
        """
        f = inc_range
        p = []
        exp = []
        self.assertEqual(f(p), exp)

        p = [(1,10),(4,7),(2,9),(5,None)]
        exp = [(1,10),(2,9),(4,7)]
        self.assertEqual(f(p), exp)

class NussinovTests(TestCase):
    """Tests for restricted nussinov algorithm and related functions"""

    def test_nussinov_fill(self):
        """nussinov_fill: basic test"""
        p = Pairs([(0,7),(1,6),(2,5),(3,9),(4,8)])
        exp = [[0,0,0,0,0,1,2,3,3,3],
                [0,0,0,0,0,1,2,2,2,2],
                [0,0,0,0,0,1,1,1,1,2],
                [0,0,0,0,0,0,0,0,1,2],
                [0,0,0,0,0,0,0,0,1,1,],
                [0,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0]]
        obs = nussinov_fill(p,size=10)
        self.assertEqual(obs, exp)

    def test_nussinov_traceback(self):
        """nussinov_traceback: basic test"""
        p = Pairs([(0,7),(1,6),(2,5),(3,9),(4,8)])
        m = nussinov_fill(p,size=10)
        exp = set([(0,7),(1,6),(2,5)])
        obs = nussinov_traceback(m, 0, 9, p)
        self.assertEqual(obs, exp)

    def test_nussinov_restricted(self):
        """nussinov_restricted: basic test"""
        p = Pairs([(0,7),(1,6),(2,5),(3,9),(4,8)])
        obs = nussinov_restricted(p)
        obs_rem = nussinov_restricted(p, return_removed=True)
        exp = [(0,7),(1,6),(2,5)]
        exp_rem = ([(0,7),(1,6),(2,5)],[(3,9),(4,8)])
        self.assertEqual(obs, exp)
        self.assertEqual(obs_rem, exp_rem)
        
        p = Pairs([(0,7),(1,6),(2,6)])
        self.assertRaises(ValueError, nussinov_restricted, p)

        p = Pairs([(0,7),(1,6),(2,5)])
        exp = Pairs([(0,7),(1,6),(2,5)])
        self.assertEqual(nussinov_restricted(p), exp)

    def test_nussinov_restricted_bi(self):
        """nussinov_restricted: include bifurcation"""
        p = Pairs([(0,7),(1,6),(2,14),(3,13),(4,12),(5,11),\
            (8,17),(9,16),(10,15)])
        obs = nussinov_restricted(p)
        obs_rem = nussinov_restricted(p, return_removed=True)
        exp = [(0,7),(1,6),(8,17),(9,16),(10,15)]
        exp_rem = ([(0,7),(1,6),(8,17),(9,16),(10,15)],\
            [(2,14),(3,13),(4,12),(5,11)])
        self.assertEqual(obs, exp)
        self.assertEqual(obs_rem, exp_rem)



if __name__ == "__main__":
    main()
