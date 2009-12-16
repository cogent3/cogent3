#!/usr/bin/env python
# test_pairs_util.py
"""Provides tests for gapping/ungapping functions and base pair comparison
"""
from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.core.sequence import RnaSequence, ModelSequence, Sequence
from cogent.core.moltype import RNA
from cogent.core.alphabet import CharAlphabet
from cogent.struct.rna2d import Pairs
from cogent.struct.pairs_util import PairsAdjustmentError,\
    adjust_base, adjust_base_structures, adjust_pairs_from_mapping,\
    delete_gaps_from_pairs, insert_gaps_in_pairs, gapped_to_ungapped,\
    get_gap_symbol, get_gap_list, degap_model_seq, degap_seq,\
    ungapped_to_gapped,\
    pairs_intersection, pairs_union, compare_pairs,\
    compare_pairs_mapping,  compare_random_to_correct,\
    sensitivity, selectivity, get_all_pairs, get_counts, extract_seqs,\
    mcc, approximate_correlation, correlation_coefficient, all_metrics
    
__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Shandy Wikman", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"
 

class GappedUngappedTests(TestCase):
    """Provides tests for gapped_to_ungapped and ungapped_to_gapped functions
    """
   
    def setUp(self):
        """setUp: set up method for all tests"""

        self.rna1 = RnaSequence('UCAG-RYN-N', Name='rna1')
        self.m1 = ModelSequence('UCAG-RYN-N', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        self.s1 = 'UCAG-RYN-N'

    def test_adjust_base(self):
        """adjust_base: should work for pairs object or list of pairs"""
        p = Pairs()
        self.assertEqual(adjust_base(p,10),[])

        pairs = [(1,21),(2,15),(3,13),(4,11),(5,10),(6,9)]
        offset = -1
        expected = [(0,20),(1,14),(2,12),(3,10),(4,9),(5,8)]
        obs_pairs = adjust_base(pairs, offset)
        self.assertEqual(obs_pairs, expected)
        
        pairs = Pairs([(0,10),(1,9)])
        self.assertEqual(adjust_base(pairs, -1), Pairs([(-1,9),(0,8)]))
        self.assertEqual(adjust_base(pairs, 5), Pairs([(5,15),(6,14)]))

        self.assertRaises(PairsAdjustmentError, adjust_base, pairs, 3.5)

    def test_adjust_base_structures(self):
        """adjust_pairs_structures: simple structure"""
        p = Pairs([(3,10),(4,9)])
        p2 = Pairs([(2,7), (30,40)])
        self.assertEqual(adjust_base_structures([p,p2], -1),\
            [[(2,9),(3,8)],[(1,6),(29,39)]])

    def test_adjust_base_None(self):
        """adjust_base: should keep Nones or duplicates, ignore conflicts"""
        pairs = Pairs([(2,8),(3,7),(6,None),(None,None),(2,10)])
        expected = Pairs([(1,7),(2,6),(5,None),(None, None),(1,9)])
        self.assertEqual(adjust_base(pairs,-1), expected)

        p = Pairs([(1,2),(2,1),(1,2),(2,None)])
        self.assertEqual(adjust_base(p, 1), [(2,3),(3,2),(2,3),(3,None)])

    def test_adjust_pairs_from_mapping_confl(self):
        """adjust_pairs_from_mapping: should handle conflicts, pseudo, dupl
        """
        f = adjust_pairs_from_mapping
        p = Pairs([(0,6),(1,5),(2,None),(None,None),(1,4),(3,7),(6,0)])
        m = {0:1,1:3,2:6,3:7,4:8,5:10,6:11,7:12}
        exp = Pairs([(1,11),(3,10),(6,None),(None,None),(3,8),(7,12),(11,1)])
        self.assertEqual(f(p, m), exp)

        p = Pairs([(1,11),(3,10),(7,12),(6,None),(None,None),(5,8)])
        m = {1: 0, 3: 1, 6: 2, 7: 3, 8: 4, 10: 5, 11: 6, 12: 7}
        exp = Pairs([(0,6),(1,5),(3,7),(2,None),(None,None)])
        self.assertEqual(f(p,m), exp)

    def test_delete_gaps_from_pairs(self):
        """delete_gaps_from_pairs: should work on standard input"""
        r = delete_gaps_from_pairs
        # empty list
        p = Pairs([])
        self.assertEqual(r(p,[1,2,3]), [])
        # normal list
        p1 = Pairs([(2,8), (3,6)])
        gap_list = [0,1,4,5,7,9]
        self.assertEqualItems(r(p1, gap_list), [(0,3),(1,2)])
        p2 = Pairs([(2,8),(3,6),(4,9)])
        self.assertEqualItems(r(p2, gap_list), [(0,3),(1,2)])
        p3 = Pairs([(2,8),(3,6),(4,10)])
        self.assertEqualItems(r(p3, gap_list), [(0,3),(1,2)])

    def test_delete_gaps_from_pairs_weird(self):
        """delete_gaps_from_pairs: should ignore conflicts etc"""
        r = delete_gaps_from_pairs
        gap_list = [0,1,4,5,7,9]
        p = Pairs([(2,6),(3,8)])
        self.assertEqualItems(r(p, gap_list), [(0,2),(1,3)])
        p = Pairs([(2,6),(3,8),(3,None),(6,2),(3,8), (None, None)])
        self.assertEqualItems(r(p, gap_list),\
            [(0,2),(1,3),(1,None),(2,0),(1,3),(None, None)])

    def test_insert_gaps_in_pairs(self):
        """insert_gaps_in_pairs: should work with normal and conflicts"""
        p = Pairs([(0,3),(1,2),(1,4),(3,None)])
        gaps = [0,1,4,5,7]
        self.assertEqual(insert_gaps_in_pairs(p, gaps),\
            [(2,8),(3,6),(3,9),(8,None)])
        p = Pairs([(0,6),(1,5),(2,None),(3,7),(0,1),(5,1)])
        gaps = [0,2,6,9]
        self.assertEqual(insert_gaps_in_pairs(p, gaps),\
            [(1,10),(3,8),(4,None),(5,11),(1,3),(8,3)])
        gaps = [2,3,4,9]
        self.assertEqual(insert_gaps_in_pairs(p, gaps),\
            [(0,10),(1,8),(5,None),(6,11),(0,1),(8,1)])
        p = Pairs([(0,6),(1,5),(2,None),(3,7),(0,1),(5,1)])
        gaps = []
        self.assertEqual(insert_gaps_in_pairs(p, gaps),\
            [(0,6),(1,5),(2,None),(3,7),(0,1),(5,1)])

    def test_get_gap_symbol(self):
        """get_gap_symbol: Sequence, ModelSequence, old_cogent, string"""
        self.assertEqual(get_gap_symbol(self.rna1), '-')
        self.assertEqual(get_gap_symbol(self.m1), '-')
        self.assertEqual(get_gap_symbol(self.s1), '-')
        self.assertEqual(get_gap_symbol(''), '-')
   
    def test_get_gap_list(self):
        """get_gap_list: Sequence, ModelSequence, old_cogent, string"""
        gs = '-'
        self.assertEqual(get_gap_list(self.rna1), [4,8])
        self.assertEqual(get_gap_list(self.m1), [4,8])
        self.assertEqual(get_gap_list(self.s1,gs),[4,8])
        self.assertEqual(get_gap_list('',gs), [])

    def test_degap_model_seq(self):
        """degap_model_seq: replacement for broken method"""
        self.assertEqual(str(degap_model_seq(self.m1)),'UCAGRYNN')

    def test_degap_seq(self):
        """degap_seq: Sequence, ModelSequence, old_cogent, string"""
        f = degap_seq
        gs = '-'
        self.assertEqual(f(self.rna1, gs), 'UCAGRYNN')
        self.assertEqual(str(f(self.m1, gs)), 'UCAGRYNN')
        self.assertEqual(f(self.s1, gs), 'UCAGRYNN')

    def test_gapped_to_ungapped(self):
        """gapped_to_ungapped: Sequence, ModelSequence, old_cogent, string
        """
        p = Pairs([(0,6),(1,5),(3,9)])
        exp = Pairs([(0,5),(1,4),(3,7)])
        f = gapped_to_ungapped
        self.assertEqual(f(self.rna1, p)[1], exp)
        self.assertEqual(f(self.m1, p)[1], exp)
        self.assertEqual(f(self.s1, p)[1], exp)

    def test_ungapped_to_gapped(self):
        """ungapped_to_gapped: Sequence, ModelSequence, old_cogent, string
        """
        p = Pairs([(0,6),(1,5),(3,9)])
        exp = Pairs([(0,5),(1,4),(3,7)])
        f = ungapped_to_gapped
        self.assertEqual(f(self.rna1, exp)[1], p)
        self.assertEqual(f(self.m1, exp)[1], p)
        self.assertEqual(f(self.s1, exp)[1], p)


class OldAdjustmentFunctionsTests(TestCase):
    """Provides tests for gapped_to_ungapped and ungapped_to_gapped functions
    """
   
    def setUp(self):
        """setUp: set up method for all tests"""
        self.ungapped = 'AGAUGCUAGCUAC'
        self.gapped = '-AGA--UGC-UAG--CUAC'
        
        self.diff_sym = '*AGA**UGC*UAG**CUAC'

        self.simple = Pairs([(2,7),(3,6),(8,12)])
        self.simple_g = Pairs([(3,11),(6,10),(12,18)])
        
        self.out_order = Pairs([(6,10),(4,1),(9,7),(5,11)])
        self.out_order_g = Pairs([(10,16),(7,2),(15,11),(8,17)])

        self.duplicates = Pairs([(3,9),(3,9),(2,10),(0,12)])
        self.duplicates_g = Pairs([(6,15),(6,15),(3,16),(1,18)])

        self.pseudo = Pairs([(0,7),(2,6),(3,10)])
        self.pseudo_g = Pairs([(1,11),(3,10),(6,16)])

    def test_adjust_base(self):
        """adjust_base: should work for pairs object or list of pairs"""
        p = Pairs()
        self.assertEqual(adjust_base(p,10),[])

        pairs = [(1,21),(2,15),(3,13),(4,11),(5,10),(6,9)]
        offset = -1
        expected = [(0,20),(1,14),(2,12),(3,10),(4,9),(5,8)]
        obs_pairs = adjust_base(pairs, offset)
        self.assertEqual(obs_pairs, expected)
        
        pairs = Pairs([(0,10),(1,9)])
        self.assertEqual(adjust_base(pairs, -1), Pairs([(-1,9),(0,8)]))
        self.assertEqual(adjust_base(pairs, 5), Pairs([(5,15),(6,14)]))

        self.assertRaises(PairsAdjustmentError, adjust_base, pairs, 3.5)

    def test_adjust_base_None(self):
        """adjust_base: should keep Nones or duplicates, ignore conflicts"""
        pairs = Pairs([(2,8),(3,7),(6,None),(None,None),(2,10)])
        expected = Pairs([(1,7),(2,6),(5,None),(None, None),(1,9)])
        self.assertEqual(adjust_base(pairs,-1), expected)

        p = Pairs([(1,2),(2,1),(1,2),(2,None)])
        self.assertEqual(adjust_base(p, 1), [(2,3),(3,2),(2,3),(3,None)])

    def test_delete_gaps_from_pairs(self):
        """delete_gaps_from_pairs: should work on standard input"""
        r = delete_gaps_from_pairs
        # empty list
        p = Pairs([])
        self.assertEqual(r(p,[1,2,3]), [])
        # normal list
        p1 = Pairs([(2,8), (3,6)])
        gap_list = [0,1,4,5,7,9]
        self.assertEqualItems(r(p1, gap_list), [(0,3),(1,2)])
        p2 = Pairs([(2,8),(3,6),(4,9)])
        self.assertEqualItems(r(p2, gap_list), [(0,3),(1,2)])
        p3 = Pairs([(2,8),(3,6),(4,10)])
        self.assertEqualItems(r(p3, gap_list), [(0,3),(1,2)])

    def test_delete_gaps_from_pairs_weird(self):
        """delete_gaps_from_pairs: should ignore conflicts etc"""
        r = delete_gaps_from_pairs
        gap_list = [0,1,4,5,7,9]
        p = Pairs([(2,6),(3,8)])
        self.assertEqualItems(r(p, gap_list), [(0,2),(1,3)])
        p = Pairs([(2,6),(3,8),(3,None),(6,2),(3,8), (None, None)])
        self.assertEqualItems(r(p, gap_list),\
            [(0,2),(1,3),(1,None),(2,0),(1,3),(None, None)])

    def test_insert_gaps_in_pairs(self):
        """insert_gaps_in_pairs: should work with normal and conflicts"""
        p = Pairs([(0,3),(1,2),(1,4),(3,None)])
        gaps = [0,1,4,5,7]
        self.assertEqual(insert_gaps_in_pairs(p, gaps),\
            [(2,8),(3,6),(3,9),(8,None)])
        p = Pairs([(0,6),(1,5),(2,None),(3,7),(0,1),(5,1)])
        gaps = [0,2,6,9]
        self.assertEqual(insert_gaps_in_pairs(p, gaps),\
            [(1,10),(3,8),(4,None),(5,11),(1,3),(8,3)])
        gaps = [2,3,4,9]
        self.assertEqual(insert_gaps_in_pairs(p, gaps),\
            [(0,10),(1,8),(5,None),(6,11),(0,1),(8,1)])

    def test_gapped_to_ungapped_simple(self):
        """gapped_to_ungapped: should work for simple case"""
        s = RnaSequence(self.gapped)
        p = self.simple_g
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.simple)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)

   
    def test_gapped_to_ungapped_out_of_order(self):
        """gapped_to_ungapped: should work when pairs are out of order
        """
        s = RnaSequence(self.gapped)
        p = Pairs(self.out_order_g)
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.out_order)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)

    def test_gapped_to_ungapped_duplicates(self):
        """gapped_to_ungapped: should work when pairs contains duplicates
        """
        s = RnaSequence(self.gapped)
        p = Pairs(self.duplicates_g)
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.duplicates)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)

    def test_gapped_to_ungapped_pseudo(self):
        """gapped_to_ungapped: shouldn't care about pseudoknots
        """
        s = RnaSequence(self.gapped)
        p = Pairs(self.pseudo_g)
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.pseudo)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)

    def test_gapped_to_ungapped_no_gaps(self):
        """gapped_to_ungapped: should return same pairs when no gaps
        """
        s = RnaSequence(self.ungapped)
        p = Pairs(self.simple)
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.simple)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)

    def test_ungapped_to_gapped(self):
        """ungapped_to_gapped: should work for basic case
        """
        s = RnaSequence(self.gapped)
        p = self.simple
        obs_seq, obs_pairs = ungapped_to_gapped(s,p)
        assert obs_seq is s
        self.assertEqualItems(obs_pairs, self.simple_g)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)

    def test_ungapped_to_gapped_out_of_order(self):
        """ungapped_to_gapped: should work when pairs out of order
        """
        s = RnaSequence(self.gapped)
        p = self.out_order
        obs_seq, obs_pairs = ungapped_to_gapped(s,p)
        assert obs_seq is s
        self.assertEqualItems(obs_pairs, self.out_order_g)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)
        
    def test_gapped_to_ungapped_simple(self):
        """gapped_to_ungapped: should work on simple case 
        """
        s = self.gapped
        p = self.simple_g
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.simple)
        assert not isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_seq, str)
        assert isinstance(obs_pairs, Pairs)
 
    def test_gapped_to_ungapped_pseudo(self):
        """gapped_to_ungapped: shouldn't care about pseudoknots
        """
        s = self.gapped
        p = self.pseudo_g
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.pseudo)
        assert not isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_seq, str)
        assert isinstance(obs_pairs, Pairs)


    def test_ungapped_to_gapped_simple(self):
        """ungapped_to_gapped: should work on basic case"""
        s = self.gapped
        p = self.simple
        obs_seq, obs_pairs = ungapped_to_gapped(s,p)
        assert obs_seq is s
        self.assertEqualItems(obs_pairs, self.simple_g)
        assert not isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_seq, str)
        assert isinstance(obs_pairs, Pairs)

    def test_ungapped_to_gapped_duplicates(self):
        """ungapped_to_gapped: should work when pairs are duplicated"""
        s = self.gapped
        p = self.duplicates
        obs_seq, obs_pairs = ungapped_to_gapped(s,p)
        assert obs_seq is s
        self.assertEqualItems(obs_pairs, self.duplicates_g)
        assert not isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_seq, str)
        assert isinstance(obs_pairs, Pairs)

    def test_gapped_to_ungapped_general(self):
        """gapped_to_ungapped: should return object of right type
        """
        s = RnaSequence(self.gapped)
        p = self.simple_g

        #in case of RnaSequence
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.simple)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)
        #in case of str input
        s = self.gapped
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.simple)
        assert not isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_seq, str)
        assert isinstance(obs_pairs, Pairs)

    def test_ungapped_to_gapped_general(self):
        """ungapped_to_gapped: should return object of right type
        """
        s = RnaSequence(self.gapped)
        p = self.simple
        #in case of RnaSequence
        obs_seq, obs_pairs = ungapped_to_gapped(s,p)
        assert obs_seq is s
        self.assertEqualItems(obs_pairs, self.simple_g)
        assert isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_pairs, Pairs)
        #in case of str input
        s = self.gapped
        obs_seq, obs_pairs = ungapped_to_gapped(s,p)
        assert obs_seq is s
        self.assertEqualItems(obs_pairs, self.simple_g)
        assert not isinstance(obs_seq, RnaSequence)
        assert isinstance(obs_seq, str)
        assert isinstance(obs_pairs, Pairs)

    def test_gapped_to_ungapped_general_seq(self):
        """gapped_to_ungapped: when input is Sequence obj, treat as string
        """
        s = Sequence(self.gapped)
        p = self.simple_g
        obs_seq, obs_pairs = gapped_to_ungapped(s,p)
        self.assertEqual(obs_seq, self.ungapped)
        self.assertEqualItems(obs_pairs, self.simple)
        #assert not isinstance(obs_seq, Sequence)
        #assert isinstance(obs_seq, str)
        assert isinstance(obs_seq, Sequence)
        assert isinstance(obs_pairs, Pairs)

    def test_adjust_pairs_from_mapping(self):
        """adjust_pairs_from_mapping: should work both ways
        """
        #ungapped to gapped
        r = RnaSequence('UC-AG-UC-CG-A-')
        u_to_g = r.gapMaps()[0] 
        #{0: 0, 1: 1, 2: 3, 3: 4, 4: 6, 5: 7, 6: 9, 7: 10, 8: 12}
        ungapped_pairs = Pairs([(0,8),(1,6),(2,5)])
        exp_pairs = Pairs([(0,12),(1,9),(3,7)])
        self.assertEqualItems(adjust_pairs_from_mapping(ungapped_pairs,\
            u_to_g), exp_pairs)

        #gapped to ungapped
        r = RnaSequence('UC-AG-UC-CG-A-')
        g_to_u = r.gapMaps()[1]
        #{0: 0, 1: 1, 3: 2, 4: 3, 6: 4, 7: 5, 9: 6, 10: 7, 12: 8}
        gapped_pairs = Pairs([(0,12),(1,9),(3,7)])
        exp_pairs = Pairs([(0,8),(1,6),(2,5)])
        self.assertEqualItems(adjust_pairs_from_mapping(gapped_pairs,\
            g_to_u), exp_pairs)

class PairsComparisonTests(TestCase):
    """Provides tests for comparing different Pairs objects"""

    def test_pairs_intersection(self):
        """pairs_intersection: should work on simple case
        """
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([(1,12),(4,9),(5,8)])
        self.assertEqualItems(pairs_intersection(p1,p2),[(4,9),(5,8)])

        #works when one is empty
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([])
        self.assertEqualItems(pairs_intersection(p1,p2),[])

        #works also on lists (not Pairs)
        p1 = [(3,10),(4,9),(5,8),(20,24)]
        p2 = [(1,12),(4,9),(5,8)]
        self.assertEqualItems(pairs_intersection(p1,p2),[(4,9),(5,8)])

    def test_pairs_intersection_duplicates(self):
        """pairs_intersection: should work on flipped pairs and duplicates
        """
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([(10,3),(4,9),(5,8),(9,4),(4,9),(23,30)])
        self.assertEqualItems(pairs_intersection(p1,p2),[(3,10),(4,9),(5,8)])

        # Conflicts, duplicates, None, pseudoknots 
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24),(22,26),(3,2),(9,4),(6,None)])
        p2 = Pairs([(1,12),(4,9),(5,8)])
        self.assertEqualItems(pairs_intersection(p1,p2),\
            [(4,9),(5,8)])

    def test_pairs_union(self):
        """pairs_union: should work on simple case
        """
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([(1,12),(4,9),(5,8)])
        self.assertEqualItems(pairs_union(p1,p2),\
            [(1,12),(3,10),(4,9),(5,8),(20,24)])

        #works when one is empty
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([])
        self.assertEqualItems(pairs_union(p1,p2),p1)

        #works also on lists (not Pairs)
        p1 = [(3,10),(4,9),(5,8),(20,24)]
        p2 = [(1,12),(4,9),(5,8)]
        self.assertEqualItems(pairs_union(p1,p2),\
            [(1,12),(3,10),(4,9),(5,8),(20,24)])
    
    def test_union_duplicates(self):
        """pairs_union: should work on flipped base pairs and duplicates
        """
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([(10,3),(4,9),(5,8),(9,4),(4,9),(23,30)])
        self.assertEqualItems(pairs_union(p1,p2),\
            [(3,10),(4,9),(5,8),(20,24),(23,30)])

        # Conflicts, duplicates, None, pseudoknots 
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24),(22,26),(3,2),(9,4),(6,None)])
        p2 = Pairs([(1,12),(4,9),(5,8)])
        self.assertEqualItems(pairs_union(p1,p2),\
            [(1,12),(3,10),(4,9),(5,8),(20,24),(22,26),(2,3)])

    def test_compare_pairs(self):
        """compare_pairs: should work on simple case"""
        #all the same
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        self.assertEqual(compare_pairs(p1,p2),1)
        
        #all different
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([(1,2),(3,4),(5,6)])
        self.assertEqual(compare_pairs(p1,p2),0)

        #one empty
        p1 = Pairs([(3,10),(4,9),(5,8),(20,24)])
        p2 = Pairs([])
        self.assertEqual(compare_pairs(p1,p2),0)

        #partially different
        p1 = Pairs([(1,2),(3,4),(5,6),(7,8)])
        p2 = Pairs([(1,2),(3,4),(9,10),(11,12)])
        self.assertFloatEqual(compare_pairs(p1,p2),.33333333333333333)

        #partially different
        p1 = Pairs([(1,2),(3,4),(5,6)])
        p2 = Pairs([(1,2),(3,4),(9,10)])
        self.assertFloatEqual(compare_pairs(p1,p2),.5)
    
    def test_compare_pairs_both_empty(self):
        """compare_pairs: should return 1.0 when both lists are empty
        """
        p1 = Pairs([])
        p2 = Pairs([])
        self.assertEqual(compare_pairs(p1,p2),1)

    def test_compare_pairs_weird(self):
        """compare_pairs: should handle conflicts, duplicates, pseudo, None
        """
        #Should raise error on conflict
        p1 = Pairs([(1,2),(3,4),(5,6),(2,None),(4,3),(None,None)])
        p2 = Pairs([(1,2),(3,4),(9,10)])
        self.assertRaises(ValueError, compare_pairs, p1, p2)

        p1 = Pairs([(1,2),(3,4),(5,6),(4,3),(None,None),(10,None)])
        p2 = Pairs([(1,2),(3,4),(9,10)])
        self.assertFloatEqual(compare_pairs(p1,p2),.5)

        p1 = Pairs([(1,8),(2,10),(7,3)])
        p2 = Pairs([(1,8),(10,2),(3,7),(4,6)])
        self.assertFloatEqual(compare_pairs(p1,p2), 0.75)

    def test_compare_pairs_mapping(self):
        """compare_pairs_mapping: should work with correct mapping
        """
        # pos in first seq, base, pos in second seq
        #1 U 0
        #2 C 1
        #3 G 2
        #4 A 3
        #  A 4
        #5 C 5
        #6 C 6
        #7 U 
        #8 G 7

        #all the same
        p1 = Pairs([(3,6),(1,8)])
        p2 = Pairs([(2,6),(0,7)])
        mapping = {1:0,2:1,3:2,4:3,5:5,6:6,7:None, 8:7}
        self.assertEqual(compare_pairs_mapping(p1,p2, mapping),1)

        #all different
        p1 = Pairs([(3,6),(1,8)])
        p2 = Pairs([(1,5),(4,7)])
        mapping = {1:0,2:1,3:2,4:3,5:5,6:6,7:None, 8:7}
        self.assertEqual(compare_pairs_mapping(p1,p2, mapping),0)
        
        #partially the same
        p1 = Pairs([(5,6),(1,4),(2,7)])
        p2 = Pairs([(5,6),(0,3),(4,7)])
        self.assertEqual(compare_pairs_mapping(p1,p2, mapping),.5)

        p1 = Pairs([(1,8),(2,7),(3,6),(4,5)])
        p2 = Pairs([(0,7),(1,6),(2,5),(3,4)])
        self.assertFloatEqual(compare_pairs_mapping(p1,p2, mapping),1/7)

        #one empty
        p1 = Pairs([(1,8),(2,7),(3,6),(4,5)])
        p2 = []
        self.assertEqual(compare_pairs_mapping(p1,p2, mapping),0)

        #both empty
        p1 = []
        p2 = []
        self.assertEqual(compare_pairs_mapping(p1,p2, mapping),1)

    def test_compare_random_to_correct(self):
        """comapre_random_to_correct: should return correct fraction
        """
        p1 = Pairs([(1,8),(2,7),(3,6),(4,5)])
        p2 = Pairs([(1,8)])
        p3 = Pairs([(1,8), (2,7), (4,5)])
        p4 = Pairs([(1,8),(2,7),(9,10),(11,12)])
        self.assertFloatEqual(compare_random_to_correct(p2,p1),1)
        self.assertFloatEqual(compare_random_to_correct(p3,p1),1)
        self.assertFloatEqual(compare_random_to_correct(p4,p1),0.5)
        self.assertFloatEqual(compare_random_to_correct([],p1),0)
        self.assertFloatEqual(compare_random_to_correct(p2,[]),0)
        self.assertFloatEqual(compare_random_to_correct([],[]),1)

class GardnerMetricsTest(TestCase):
    """Tests for the metrics from Gardner & Giegerich 2004"""
    
    def setUp(self):
        """setUp: setup method for all tests"""
        self.true = Pairs([(0,40),(1,39),(2,38),(3,37),(10,20),\
            (11,19),(12,18),(13,17),(26,33),(27,32)])
        self.predicted = Pairs([(0,40),(1,39),(2,38),(3,37),(4,36),\
            (5,35),(10,22),(11,20),(14,29),(15,28)])
        self.seq = ['>seq1\n','agguugaaggggauccgauccacuccccggcuggucaaccu']

    def test_conflicts(self):
        """all metrics should raise error when conflicts in one of the structs
        """
        ref = Pairs([(1,6),(2,5),(3,10),(7,None),(None,None),(5,2),(1,12)])
        pred = Pairs([(6,1),(10,11),(3,12)])
        
        self.assertRaises(ValueError, sensitivity, ref, pred)
        self.assertRaises(ValueError, sensitivity, pred, ref)
        self.assertRaises(ValueError, selectivity, ref, pred)
        self.assertRaises(ValueError, selectivity, pred, ref)
        self.assertRaises(ValueError, approximate_correlation, ref, pred,\
            self.seq)
        self.assertRaises(ValueError, approximate_correlation, pred, ref,\
            self.seq)
        self.assertRaises(ValueError, correlation_coefficient, ref, pred,\
            self.seq)
        self.assertRaises(ValueError, correlation_coefficient, pred, ref,\
            self.seq)
        self.assertRaises(ValueError, mcc, ref, pred, self.seq)
        self.assertRaises(ValueError, mcc, pred, ref, self.seq)

    def test_get_all_pairs(self):
        """get_all_pairs: should return the number of possible pairs"""
        seq = RnaSequence('UCAG-NACGU')
        seq2 = RnaSequence('UAAG-CACGC')
        self.assertEqual(get_all_pairs([seq], min_dist=4), 6)
        self.assertEqual(get_all_pairs([seq2], min_dist=4), 4)
        # when given multiple sequences, should average over all of them
        self.assertEqual(get_all_pairs([seq,seq2], min_dist=4), 5)
        # different min distance
        self.assertEqual(get_all_pairs([seq], min_dist=2),10)
        # error on invalid minimum distance
        self.assertRaises(ValueError, get_all_pairs, [seq], min_dist=-2)

    def test_get_counts(self):
        """get_counts: should work with all parameters"""
        seq = RnaSequence('UCAG-NAUGU')
        seq2 = RnaSequence('UAAG-CACGC')
        p = Pairs([(1,8),(2,7)])
        p2 = Pairs([(1,8),(2,6),(3,6),(4,9),])
        exp = {'TP':1,'TN':0, 'FN':1,'FP':3,\
            'FP_INCONS':0, 'FP_CONTRA':0, 'FP_COMP':0}
        self.assertEqual(get_counts(p, p2), exp)
        exp = {'TP':1,'TN':0, 'FN':1,'FP':3,\
            'FP_INCONS':1, 'FP_CONTRA':1, 'FP_COMP':1}
        self.assertEqual(get_counts(p, p2, split_fp=True), exp)
        seq = RnaSequence('UCAG-NACGU')
        exp = {'TP':1,'TN':7, 'FN':1,'FP':3,\
            'FP_INCONS':1, 'FP_CONTRA':1, 'FP_COMP':1}
        self.assertEqual(get_counts(p, p2, split_fp=True,\
            sequences=[seq], min_dist=2), exp)
        # check against compare_ct.pm
        exp = {'TP':4,'TN':266, 'FN':6,'FP':6,\
            'FP_INCONS':2, 'FP_CONTRA':2, 'FP_COMP':2}
        seq = 'agguugaaggggauccgauccacuccccggcuggucaaccu'.upper()
        self.assertEqual(get_counts(self.true, self.predicted, split_fp=True,\
            sequences=[seq], min_dist=4), exp)

    def test_extract_seqs(self):
        """extract_seqs: should handle different input formats"""
        s1 = ">seq1\nACGUAGC\n>seq2\nGGUAGCG"
        s2 = [">seq1","ACGUAGC",">seq2","GGUAGCG"]
        s3 = ['ACGUAGC','GGUAGCG']
        s4 = [RnaSequence('ACGUAGC'), RnaSequence('GGUAGCG')]
        m1 = ModelSequence('ACGUAGC', Name='rna1',\
            Alphabet=RNA.Alphabets.DegenGapped)
        m2 = ModelSequence('GGUAGCG', Name='rna2',\
            Alphabet=RNA.Alphabets.DegenGapped)
        s5 = [m1, m2]
        f = extract_seqs
        self.assertEqual(f(s1), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s2), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s3), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s4), ['ACGUAGC', 'GGUAGCG'])
        self.assertEqual(f(s5), ['ACGUAGC', 'GGUAGCG'])

    def test_sensitivity(self):
        """sensitivity: check against compare_ct.pm"""
        sen = sensitivity(self.true,self.predicted)
        self.assertEqual(sen, 0.4)

    def test_sensitivity_general(self):
        """sensitivity: should work in general"""
        ref = Pairs([(1,6),(2,5),(3,10)])
        pred = Pairs([(6,1),(10,11),(3,12)])
        # one good prediction
        self.assertFloatEqual(sensitivity(ref, pred), 1/3)
        # over-prediction not penalized
        pred = Pairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        self.assertFloatEqual(sensitivity(ref, pred), 1/3)

    def test_sensitivity_dupl(self):
        """sensitivity: should handle duplicates, pseudo, None"""
        ref = Pairs([(1,6),(2,5),(3,10),(7,None),(None,None),(5,2),(4,9)])
        pred = Pairs([(6,1),(10,11),(3,12)])
        self.assertFloatEqual(sensitivity(ref, pred), 0.25)

        pred = Pairs([(6,1),(10,11),(3,12),(20,None),(None,None),(1,6)])
        self.assertFloatEqual(sensitivity(ref, pred), 0.25)

    def test_sensitivity_empty(self):
        """sensitivity: should work on emtpy Pairs"""
        # both empty
        self.assertFloatEqual(sensitivity(Pairs(), Pairs()), 1)
        pred = Pairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        # prediction emtpy
        self.assertFloatEqual(sensitivity(Pairs(), pred), 0)
        # reference empty
        self.assertFloatEqual(sensitivity(pred, Pairs()), 0)

    def test_selectivity(self):
        """selectivity: check against compare_ct.pm"""
        sel = selectivity(self.true,self.predicted)
        self.assertEqual(sel, 0.5)

    def test_selectivity_general(self):
        """selectivity: should work in general"""
        ref = Pairs([(1,6),(2,5),(10,13)])
        pred = Pairs([(6,1),(3,4),(10,12)])
        # one good prediction
        self.assertFloatEqual(selectivity(ref, pred), 0.5)
        # over-prediction not penalized
        pred = Pairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        self.assertFloatEqual(selectivity(ref, pred), 0.25)

    def test_selectivity_dupl(self):
        """selectivity: duplicates and Nones shouldn't influence the calc.
        """
        ref = Pairs([(1,6),(2,5),(10,13),(6,1),(7,None),(None,None)])
        pred = Pairs([(6,1),(3,4),(10,12)])
        self.assertFloatEqual(selectivity(ref, pred), 0.5)

    def test_selectivity_empty(self):
        """selectivity: should handle empty reference/predicted structure"""
        # both empty
        self.assertFloatEqual(selectivity(Pairs(), Pairs()), 1)
        pred = Pairs([(6,1),(10,11),(3,12),(13,20),(14,19),(15,18)])
        # prediction emtpy
        self.assertFloatEqual(selectivity(Pairs(), pred), 0)
        # reference empty
        self.assertFloatEqual(selectivity(pred, Pairs()), 0)

    def test_approximate_correlation(self):
        """approximate_correlation: check against compare_ct.pm"""
        self.assertFloatEqual(approximate_correlation(self.true,\
            self.predicted, seqs=self.seq), 0.45)

    def test_correlation_coefficient(self):
        """correlation_coefficient: check against compare_ct.pm"""
        self.assertFloatEqual(correlation_coefficient(self.true,\
            self.predicted, seqs=self.seq, min_dist=4), 0.42906394)

    def test_cc_bad_pred(self):
        """correlation_coefficient: should give 0 when TP=0"""
        ref = Pairs([(1,7),(2,5)])
        pred = Pairs([(0,8)])
        seqs = ['CAUCGAUUG']
        self.assertEqual(correlation_coefficient(ref, pred, seqs=seqs), 0.0)

    def test_mcc(self):
        """mcc: check against compare_ct.pm"""
        res = mcc(self.true,self.predicted,self.seq, min_dist=4)
        self.assertFloatEqual(res, 0.42906394)
        
    def test_all_metrics(self):
        """all_metrics: check against compare_ct.pm"""
        exp = {'SENSITIVITY':0.4, 'SELECTIVITY':0.5, 'AC':0.45,\
            'CC':0.42906394, 'MCC':0.42906394}
        obs = all_metrics(self.true, self.predicted, seqs=self.seq, min_dist=4)
        self.assertEqualItems(obs.keys(), exp.keys())
        for k in exp:
            self.assertFloatEqual(obs[k], exp[k])

    def test_get_counts_pseudo(self):
        """get_counts: should work when pseudo in ref -> classification off"""
        # pairs that would normally be compatible, are now contradicting
        ref = Pairs([(0,8),(1,7),(4,10)])
        pred = Pairs([(0,8),(3,6),(4,10)])
        seq = 'GACUGUGUCAU'
        exp = {'TP':2,'TN':13-2-1, 'FN':1,'FP':1,\
            'FP_INCONS':0, 'FP_CONTRA':1, 'FP_COMP':0}
        self.assertEqual(get_counts(ref, pred, split_fp=True,\
            sequences=[seq], min_dist=4), exp)

    def test_all_metrics_pseudo(self):
        """all_metrics: pseudoknot in ref, check against compare_ct.pm"""
        ref = Pairs([(0,8),(1,7),(4,10)])
        pred = Pairs([(0,8),(3,6),(4,10)])
        seq = 'GACUGUGUCAU'
        exp = {'SENSITIVITY':0.6666667, 'SELECTIVITY':0.6666667,\
            'AC':0.6666667, 'CC':0.57575758, 'MCC':0.57575758}
        obs = all_metrics(ref, pred, seqs=[seq], min_dist=4)
        self.assertEqualItems(obs.keys(), exp.keys())
        for k in exp:
            self.assertFloatEqual(obs[k], exp[k])

    def test_all_metrics_weird_input(self):
        """all_metrics: should work when ref or prediction empty or no seqs"""
        ref = Pairs([(3,10)])
        pred = Pairs()
        seqs = ['UACGUAGCUAGCUAGCUACG']
        obs = all_metrics(ref, pred, seqs=[seqs], min_dist=4)
        exp = {'SENSITIVITY':0, 'SELECTIVITY':0,\
            'AC':0, 'CC':0, 'MCC':0}
        for k in exp:
            self.assertFloatEqual(obs[k], exp[k])

        ref = Pairs()
        pred = Pairs()
        seqs = ['UACGUAGCUAGCUAGCUACG']
        obs = all_metrics(ref, pred, seqs=[seqs], min_dist=4)
        exp = {'SENSITIVITY':1, 'SELECTIVITY':1,\
            'AC':1, 'CC':1, 'MCC':1}
        for k in exp:
            self.assertFloatEqual(obs[k], exp[k])

        ref = Pairs([(3,10)])
        pred = Pairs([(1,12)])
        seqs = ['UACGUAGCUAGCUAGCUACG']
        self.assertRaises(ValueError, all_metrics, ref, pred, seqs="",\
            min_dist=4)


if __name__ == "__main__":
    main()
