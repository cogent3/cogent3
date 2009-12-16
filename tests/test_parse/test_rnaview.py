#!/usr/bin/env python
#test_rnaview.py
"""Provides tests for code in the rnaview.py file.
"""

from cogent.util.unit_test import TestCase, main
from cogent.parse.rnaview import is_roman_numeral, is_edge, is_orientation,\
    parse_annotation, parse_uncommon_residues, parse_base_pairs,\
    parse_base_multiplets, parse_pair_counts, MinimalRnaviewParser,\
    RnaviewParser, Base, BasePair, BasePairs, BaseMultiplet, BaseMultiplets,\
    PairCounts, RnaViewObjectError, RnaViewParseError, MinimalRnaviewParser,\
    parse_filename, parse_number_of_pairs, verify_bp_counts,\
    in_chain, is_canonical, is_not_canonical, is_stacked, is_not_stacked,\
    is_tertiary, is_not_stacked_or_tertiary, is_tertiary_base_base

__author__ = "Greg Caporaso and Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Greg Caporaso", "Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

#==============================================================================
# RNAVIEW OBJECTS TESTS
#==============================================================================

class BaseTests(TestCase):
    """Tests for Base class"""
    
    def test_init(self):
        """Base __init__: should initialize on standard data"""
        b = Base('A','30','G')
        self.assertEqual(b.ChainId, 'A')
        self.assertEqual(b.ResId, '30')
        self.assertEqual(b.ResName, 'G')
        #ResId or ResName can't be None or empty string
        self.assertRaises(RnaViewObjectError, Base, None, None, 'G')
        self.assertRaises(RnaViewObjectError, Base, '1', 'A', '')
        self.assertRaises(RnaViewObjectError, Base, None, '', 'C')
        #Can pass RnaViewSeqPos (str)
        b = Base('C','12','A','10')
        self.assertEqual(b.RnaViewSeqPos, '10')

    def test_str(self):
        """Base __str__: should return correct string"""
        b = Base('A','30','G')
        self.assertEqual(str(b), 'A 30 G')

    def test_eq(self):
        """Base ==: functions as expected """
        # Define a standard to compare others
        b = Base('A','30','G')
        # Identical to b
        b_a = Base('A','30','G')
        # Differ in Chain from b
        b_b = Base('B','30','G')
        # Differ in ResId from b
        b_c = Base('A','25','G')
        # Differ in ResName from b
        b_d = Base('A','30','C')
        # Differ in RnaViewSeqPos
        b_e = Base('A','30','G','2')
        # Differ in everything from b
        b_e = Base('C','12','U','1')
        
        self.assertEqual(b == b, True)
        self.assertEqual(b_a == b, True)
        self.assertEqual(b_b == b, False)
        self.assertEqual(b_c == b, False)
        self.assertEqual(b_d == b, False)
        self.assertEqual(b_e == b, False)

    def test_ne(self):
        """Base !=: functions as expected"""
        # Define a standard to compare others
        b = Base('A','30','G')
        # Identical to b
        b_a = Base('A','30','G')
        # Differ in Chain from b
        b_b = Base('B','30','G')
        # Differ in ResId from b
        b_c = Base('A','25','G')
        # Differ in ResName from b
        b_d = Base('A','30','C')
        # Differ in everything from b
        b_e = Base('C','12','U')
        
        self.assertEqual(b != b, False)
        self.assertEqual(b_a != b, False)
        self.assertEqual(b_b != b, True)
        self.assertEqual(b_c != b, True)
        self.assertEqual(b_d != b, True)
        self.assertEqual(b_e != b, True)


class BasePairTests(TestCase):
    """Tests for BasePair object"""

    def setUp(self):
        """setUp method for all tests in BasePairTests"""
        self.b1 = Base('A','30','G')
        self.b2 = Base('A','36','C')
        self.bp = BasePair(self.b1, self.b2, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation=None)

    def test_init(self):
        """BasePair __init__: should initialize on standard data"""
        self.failUnless(self.bp.Up is self.b1)
        self.failUnless(self.bp.Down is self.b2)
        self.failUnless(self.bp.Conformation is None)
        self.assertEqual(self.bp.Edges, 'H/W')
        self.assertEqual(self.bp.Orientation, 'cis')
        self.assertEqual(self.bp.Saenger, 'XI')

    def test_str(self):
        """BasePair __str__: should return correct string"""
        self.assertEqual(str(self.bp),
            "Bases: A 30 G -- A 36 C; Annotation: H/W -- cis -- "+\
            "None -- XI;")

    def test_eq(self):
        "BasePair ==: should function as expected"
        # identical
        up = Base('A','30','G')
        down = Base('A','36','C')
        bp = BasePair(up, down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation=None) 
        self.assertEqual(bp == self.bp, True)
        # diff up base
        diff_up = Base('C','12','A')
        bp = BasePair(diff_up, down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp == self.bp, False)
        # diff down base
        diff_down = Base('D','13','U')
        bp = BasePair(up, diff_down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp == self.bp, False)
        # diff edges
        bp = BasePair(up, down, Edges='W/W', Saenger='XI',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp == self.bp, False)
        # diff orientation
        bp = BasePair(up, down, Edges='H/W', Saenger='XI',\
            Orientation='tran', Conformation=None)
        self.assertEqual(bp == self.bp, False)
        # diff conformation
        bp = BasePair(up, down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation='syn')
        self.assertEqual(bp == self.bp, False)
        # diff saenger
        bp = BasePair(up, down, Edges='H/W', Saenger='XIX',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp == self.bp, False)

    def test_ne(self):
        "BasePair !=: should function as expected"
        # identical
        up = Base('A','30','G')
        down = Base('A','36','C')
        bp = BasePair(up, down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation=None) 
        self.assertEqual(bp != self.bp, False)
        # diff up base
        diff_up = Base('C','12','A')
        bp = BasePair(diff_up, down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp != self.bp, True)
        # diff down base
        diff_down = Base('D','13','U')
        bp = BasePair(up, diff_down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp != self.bp, True)
        # diff edges
        bp = BasePair(up, down, Edges='W/W', Saenger='XI',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp != self.bp, True)
        # diff orientation
        bp = BasePair(up, down, Edges='H/W', Saenger='XI',\
            Orientation='tran', Conformation=None)
        self.assertEqual(bp != self.bp, True)
        # diff conformation
        bp = BasePair(up, down, Edges='H/W', Saenger='XI',\
            Orientation='cis', Conformation='syn')
        self.assertEqual(bp != self.bp, True)
        # diff saenger
        bp = BasePair(up, down, Edges='H/W', Saenger='XIX',\
            Orientation='cis', Conformation=None)
        self.assertEqual(bp != self.bp, True)

    def test_isWC(self):
        """BasePair isWC: should return True for GC or AU pair"""
        bp = BasePair(Base('A','30','G'), Base('A','36','C'))
        self.assertEqual(bp.isWC(), True)
        bp = BasePair(Base('A','30','g'), Base('A','36','C'))
        self.assertEqual(bp.isWC(), True)
        bp = BasePair(Base('A','30','C'), Base('A','36','G'))
        self.assertEqual(bp.isWC(), True)
        bp = BasePair(Base('A','30','U'), Base('A','36','a'))
        self.assertEqual(bp.isWC(), True)
        bp = BasePair(Base('A','30','G'), Base('A','36','U'))
        self.assertEqual(bp.isWC(), False)

    def test_isWobble(self):
        """BasePair isWobble: should return True for GU pair"""
        bp = BasePair(Base('A','30','G'), Base('A','36','U'))
        self.assertEqual(bp.isWobble(), True)
        bp = BasePair(Base('A','30','g'), Base('A','36','U'))
        self.assertEqual(bp.isWobble(), True)
        bp = BasePair(Base('A','30','u'), Base('A','36','g'))
        self.assertEqual(bp.isWobble(), True)
        bp = BasePair(Base('A','30','A'), Base('A','36','U'))
        self.assertEqual(bp.isWobble(), False)

class BasePairsTests(TestCase):
    """Tests for BasePairs object"""

    def setUp(self):
        """setUp method for all BasePairs tests"""
        self.a1 = BasePair(Base('A','30','G'), Base('A','36','C'), Saenger='XX')
        self.a2 = BasePair(Base('A','31','A'), Base('A','35','U'), Saenger='XX')
        self.a3 = BasePair(Base('A','40','G'), Base('A','60','U'), Saenger='V')
        self.a4 = BasePair(Base('A','41','A'), Base('A','58','U'),\
            Saenger=None)
        self.ab1 = BasePair(Base('A','41','A'), Base('B','58','U'))
        self.ac1 = BasePair(Base('A','10','C'), Base('C','3','G'))
        self.bc1 = BasePair(Base('B','41','A'), Base('C','1','U'))
        self.bn1 = BasePair(Base('B','41','A'), Base(None,'1','U'))
        self.cd1 = BasePair(Base('C','41','A'), Base('D','1','U'))

        self.bp1 = BasePair(Base('A','34','U'), Base('A','40','A'))
        self.bp2 = BasePair(Base('A','35','C'), Base('A','39','G'))
        self.bp3 = BasePair(Base('B','32','G'), Base('B','38','U'))
        self.bp4 = BasePair(Base('B','33','G'), Base('B','37','C'))
        self.bp5 = BasePair(Base('A','31','C'), Base('B','41','G'))
        self.bp6 = BasePair(Base('A','32','U'), Base('B','40','A'))
        self.bp7 = BasePair(Base('A','37','U'), Base('B','35','A'))
        
        self.pairs = BasePairs([self.bp1, self.bp2, self.bp3, self.bp4,\
            self.bp5, self.bp6, self.bp7])
        
    def test_init(self):
        """BasePairs __init__: should work with or without Model"""
        # init from list
        bps = BasePairs([self.a1, self.a2])
        self.failUnless(bps[0] is self.a1)
        self.failUnless(bps[1] is self.a2)

        # init from tuple
        bps = BasePairs((self.a1, self.a2))
        self.failUnless(bps[0] is self.a1)
        self.failUnless(bps[1] is self.a2)

    def test_str(self):
        """BasePairs __str__: should produce expected string"""
        b1 = BasePair(Base('A','30','G'), Base('A','36','C'), Saenger='XX')
        b2 = BasePair(Base('A','31','A'), Base('A','35','U'),\
            Orientation='cis',Edges='W/W')
        bps = BasePairs([b1, b2])
        exp_lines = [
        "===================================================================",\
        "Bases: Up -- Down; Annotation: Edges -- Orient. -- Conf. -- Saenger",\
        "===================================================================",\
        "Bases: A 30 G -- A 36 C; Annotation: None -- None -- None -- XX;",\
        "Bases: A 31 A -- A 35 U; Annotation: W/W -- cis -- None -- None;"]
        self.assertEqual(str(bps), '\n'.join(exp_lines))

    def test_select(self):
        """BasePairs select: should work with any good function"""

        def xx(bp):
            if bp.Saenger == 'XX':
                return True
            return False

        bps = BasePairs([self.a1, self.a2, self.a3, self.a4]) 

        obs = bps.select(xx)
        self.assertEqual(len(obs), 2)
        self.failUnless(obs[0] is self.a1)
        self.failUnless(obs[1] is self.a2)
        for i in obs:
            self.assertEqual(i.Saenger, 'XX')

    def test_PresentChains(self):
        """BasePairs PresentChains: should work on single/multiple chain(s)"""
        bps = BasePairs([self.a1, self.a2, self.a3, self.a4])
        self.assertEqual(bps.PresentChains, ['A'])
        bps = BasePairs([self.a1, self.ab1])
        self.assertEqualItems(bps.PresentChains, ['A','B'])
        bps = BasePairs([self.a1, self.ab1, self.ac1, self.bc1])
        self.assertEqualItems(bps.PresentChains, ['A','B', 'C'])
        bps = BasePairs([self.a1, self.ab1, self.bn1])
        self.assertEqualItems(bps.PresentChains, [None, 'A','B'])

    def test_cliques(self):
        """BasePairs cliques: single/multiple chains and cliques"""
        #one chain, one clique
        bps = BasePairs([self.a1, self.a2, self.a3, self.a4])
        obs_cl = list(bps.cliques())
        self.assertEqual(len(obs_cl), 1)
        
        #3 chains, 2 cliques
        bps = BasePairs([self.a1, self.a2, self.cd1])
        obs_cl = list(bps.cliques())
        self.assertEqual(len(obs_cl), 2)
        self.assertEqual(len(obs_cl[0]), 2)
        self.assertEqual(len(obs_cl[1]), 1)
        self.failUnless(obs_cl[1][0] is self.cd1)
        self.assertEqual(obs_cl[1].PresentChains, ['C','D'])

        #5 chains, 1 clique
        bps = BasePairs([self.a1, self.ab1, self.bc1, self.bn1, self.cd1])
        obs_cl = list(bps.cliques())
        self.assertEqual(len(obs_cl), 1)
        self.assertEqual(len(obs_cl[0]), 5)
        self.failUnless(obs_cl[0][0] is self.a1)
        self.assertEqualItems(obs_cl[0].PresentChains, ['A','B','C','D', None])
    
    def test_hasConflicts(self):
        """BasePairs hadConflicts: handle chains and residue IDs"""
        # no conflict
        b1 = BasePair(Base('A','30','G'), Base('A','36','C'), Saenger='XX')
        b2 = BasePair(Base('A','31','A'), Base('A','35','U'),\
            Orientation='cis',Edges='W/W')
        b3 = BasePair(Base('A','15','G'), Base('A','42','C'))
        bps = BasePairs([b1, b2, b3])
        self.assertEqual(bps.hasConflicts(), False)
        self.assertEqual(bps.hasConflicts(return_conflict=True), (False, None))

        # conflict within chain
        b1 = BasePair(Base('A','30','G'), Base('A','36','C'), Saenger='XX')
        b2 = BasePair(Base('A','31','A'), Base('A','35','U'),\
            Orientation='cis',Edges='W/W')
        b3 = BasePair(Base('A','30','G'), Base('A','42','C'))
        bps = BasePairs([b1, b2, b3])
        self.assertEqual(bps.hasConflicts(), True)

        # conflict within chain -- return conflict
        b1 = BasePair(Base('A','30','G'), Base('A','36','C'), Saenger='XX')
        b2 = BasePair(Base('A','31','A'), Base('A','35','U'),\
            Orientation='cis',Edges='W/W')
        b3 = BasePair(Base('A','30','G'), Base('A','42','C'))
        bps = BasePairs([b1, b2, b3])
        self.assertEqual(bps.hasConflicts(return_conflict=True),\
            (True, "A 30 G"))

        # no conflict, same residue ID, different chain
        b1 = BasePair(Base('A','30','G'), Base('A','36','C'), Saenger='XX')
        b2 = BasePair(Base('A','31','A'), Base('A','35','U'),\
            Orientation='cis',Edges='W/W')
        b3 = BasePair(Base('C','30','G'), Base('A','42','C'))
        bps = BasePairs([b1, b2, b3])
        self.assertEqual(bps.hasConflicts(), False)


class BaseMultipletTests(TestCase):
    """Tests for BaseMultiplet object"""

    def test_init(self):
        """BaseMultiplet __init__: should work as expected"""
        b1 = Base('A','30','A')
        b2 = Base('A','35','G')
        b3 = Base('A','360','U')
        bm = BaseMultiplet([b1, b2, b3])
        self.failUnless(bm[0] is b1)
        self.failUnless(bm[2] is b3)
        #should work from tuple also
        bm = BaseMultiplet((b1, b2, b3))
        self.failUnless(bm[0] is b1)
        self.failUnless(bm[2] is b3)

    def test_str(self):
        """BaseMultiplet __str__: should give expected string"""
        b1 = Base('A','30','A')
        b2 = Base('A','35','G')
        b3 = Base('A','360','U')
        bm = BaseMultiplet([b1, b2, b3])
        exp = "A 30 A -- A 35 G -- A 360 U;"
        self.assertEqual(str(bm), exp)

class BaseMultipletsTests(TestCase):
    """Tests for BaseMultiplets object"""
    
    def test_init(self):
        """BaseMultiplets __init__: from list and tuple"""
        b1 = Base('A','30','A')
        b2 = Base('A','35','G')
        b3 = Base('A','360','U')
        bm1 = BaseMultiplet([b1, b2, b3])
        b4 = Base('B','12','C')
        b5 = Base('B','42','U')
        b6 = Base('C','2','A')
        bm2 = BaseMultiplet([b4, b5, b6])
        bms = BaseMultiplets([bm1, bm2])
        self.failUnless(bms[0] is bm1)
        self.failUnless(bms[1] is bm2)
        self.assertEqual(bms[1][2].ResId, '2')
        
        #should work from tuple also
        bms = BaseMultiplets((bm1, bm2))
        self.failUnless(bms[0] is bm1)
        self.failUnless(bms[1] is bm2)
        self.assertEqual(bms[1][2].ResId, '2')

    def test_str(self):
        """BaseMultiplets __str__: should give expected string"""
        b1 = Base('A','30','A')
        b2 = Base('A','35','G')
        b3 = Base('A','360','U')
        bm1 = BaseMultiplet([b1, b2, b3])
        b4 = Base('B','12','C')
        b5 = Base('B','42','U')
        b6 = Base('C','2','A')
        bm2 = BaseMultiplet([b4, b5, b6])
        bms = BaseMultiplets([bm1, bm2])
        exp_lines = [\
        "A 30 A -- A 35 G -- A 360 U;",\
        "B 12 C -- B 42 U -- C 2 A;"]
        self.assertEqual(str(bms), '\n'.join(exp_lines))

class TestPairCounts(TestCase):
    """Tests for PairCounts object.
    
    Contains only test for __init__. Everything should fucntion as a dict.
    """

    def test_init(self):
        """PairCounts __init__: should work as dict"""
        res = PairCounts(\
            {'Standard':1, 'WS--cis':300, 'Bifurcated': 2, 'HS-tran':0})
        self.assertEqual(res['Standard'], 1)
        self.assertEqual(res['WS--cis'], 300)
        self.assertEqual(res['Bifurcated'], 2)
        self.assertEqual(res['HS-tran'], 0)

#==============================================================================
# SELECTION FUNCTIONS TESTS
#==============================================================================

class SelectionFunctionTests(TestCase):
   
    def test_in_chain(self):
        b1, b2 = Base('A','30','C'), Base('A','40','G')
        bp1 = BasePair(b1, b2)
        self.assertEqual(in_chain("A")(bp1), True)
        self.assertEqual(in_chain(["A",'B'])(bp1), True)
        self.assertEqual(in_chain("B")(bp1), False)

        b3, b4 = Base('A','30','C'), Base('B','40','G')
        bp2 = BasePair(b3,b4)
        self.assertEqual(in_chain("A")(bp2), False)
        self.assertEqual(in_chain(["A",'B'])(bp2), True)
        self.assertEqual(in_chain("AB")(bp2), True)
        self.assertEqual(in_chain("AC")(bp2), False)

        b5, b6 = Base('A','30','C'), Base('C','40','G')
        bp3 = BasePair(b5,b6)
        self.assertEqual(in_chain("A")(bp3), False)
        self.assertEqual(in_chain("C")(bp3), False)
        self.assertEqual(in_chain(["A",'B'])(bp3), False)
        self.assertEqual(in_chain("AC")(bp3), True)


    def test_is_canocical(self):
        """is_canonical: work on annotation, not base identity"""
        b1, b2 = Base('A','30','C'), Base('A','40','G')
        bp = BasePair(b1, b2, Edges='+/+')
        self.assertEqual(is_canonical(bp), True)
        bp = BasePair(b1, b2, Edges='-/-')
        self.assertEqual(is_canonical(bp), True)
        bp = BasePair(b1, b2, Edges='W/W')
        self.assertEqual(is_canonical(bp), False)
        bp = BasePair(b1, b2, Edges='W/W', Orientation='cis',Saenger='XXVIII')
        self.assertEqual(is_canonical(bp), True)

    def test_is_not_canocical(self):
        """is_not_canonical: opposite of is_canonical"""
        b1, b2 = Base('A','30','C'), Base('A','40','G')
        bp = BasePair(b1, b2, Edges='+/+')
        self.assertEqual(is_not_canonical(bp), False)
        bp = BasePair(b1, b2, Edges='-/-')
        self.assertEqual(is_not_canonical(bp), False)
        bp = BasePair(b1, b2, Edges='W/W')
        self.assertEqual(is_not_canonical(bp), True)
        bp = BasePair(b1, b2, Edges='W/W', Orientation='cis',Saenger='XXVIII')
        self.assertEqual(is_not_canonical(bp), False)

    def test_is_stacked(self):
        """is_stacked: checks annotation, not base identity"""
        b1, b2 = Base('A','30','C'), Base('A','40','A')
        bp = BasePair(b1, b2, Edges='stacked')
        self.assertEqual(is_stacked(bp), True)
        bp = BasePair(b1, b2, Edges='H/?')
        self.assertEqual(is_stacked(bp), False)

    def test_is_not_stacked(self):
        """is_not_stacked: opposite of is_stacked"""
        b1, b2 = Base('A','30','C'), Base('A','40','A')
        bp = BasePair(b1, b2, Edges='stacked')
        self.assertEqual(is_not_stacked(bp), False)
        bp = BasePair(b1, b2, Edges='H/?')
        self.assertEqual(is_not_stacked(bp), True)

    def test_is_tertiary(self):
        """is_tertiary: checks annotation, not base identity"""
        b1, b2 = Base('A','30','C'), Base('A','40','U')
        bp = BasePair(b1, b2, Saenger='!1H(b_b)')
        self.assertEqual(is_tertiary(bp), True)
        bp = BasePair(b1, b2, Edges='H/?', Saenger='XX')
        self.assertEqual(is_tertiary(bp), False)
        bp = BasePair(b1,b2, Edges='stacked')
        self.assertEqual(is_tertiary(bp), False)

    def test_is_not_stacked_or_tertiary(self):
        """is_not_stacked_or_tertiary: checks annotation, not base identity"""
        b1, b2 = Base('A','30','C'), Base('A','40','U')
        bp = BasePair(b1, b2, Saenger='!1H(b_b)')
        self.assertEqual(is_not_stacked_or_tertiary(bp), False)
        bp = BasePair(b1, b2, Edges='stacked')
        self.assertEqual(is_not_stacked_or_tertiary(bp), False)
        bp = BasePair(b1, b2, Edges='W/W', Saenger='XX')
        self.assertEqual(is_not_stacked_or_tertiary(bp), True)

    def test_is_tertiary_base_base(self):
        """is_tertiary_base_base: checks annotation, not base identity"""
        b1, b2 = Base('A','30','C'), Base('A','40','U')
        bp = BasePair(b1, b2, Saenger='!1H(b_b)')
        self.assertEqual(is_tertiary_base_base(bp), True)
        bp = BasePair(b1, b2, Edges='H/?', Saenger='!(s_s)')
        self.assertEqual(is_tertiary_base_base(bp), False)

#==============================================================================
# RNAVIEW PARSER TESTS
#==============================================================================

class RnaviewParserTests(TestCase):
    """Tests for RnaviewParser and related code"""

    def test_is_roman_numeral(self):
        """is_roman_numeral: should work for all, including comma"""
        self.assertEqual(is_roman_numeral('XIII'),True)
        self.assertEqual(is_roman_numeral('Xiii'),False)
        self.assertEqual(is_roman_numeral('MMCDXXVIII'),True)
        self.assertEqual(is_roman_numeral('XII,XIII'),True)
        self.assertEqual(is_roman_numeral('n/a'),False)

    def test_is_edge(self):
        """is_edge: should identify valid edges correctly"""
        self.assertEqual(is_edge('H/W'),True)
        self.assertEqual(is_edge('./W'),True)
        self.assertEqual(is_edge('+/+'),True)
        self.assertEqual(is_edge('   '),False)
        self.assertEqual(is_edge('P/W'),False)
        self.assertEqual(is_edge('X/W'),True)
        self.assertEqual(is_edge('X/X'),True)

    def test_is_orientation(self):
        """is_orientation: should fail on anything but 'cis' or 'tran'"""
        self.assertEqual(is_orientation('cis'),True)
        self.assertEqual(is_orientation('tran'),True)
        self.assertEqual(is_orientation('tranxxx'),False)

    def test_parse_annotation(self):
        """parse_annotation: should return correct tuple of 4 or raise error
        """
        self.assertEqual(parse_annotation(['W/S', 'tran', 'syn', 'syn',\
            'n/a']), ('W/S','tran','syn syn','n/a'))
        self.assertEqual(parse_annotation(['syn','stacked']),\
            ('stacked', None,'syn',None))
        self.assertEqual(parse_annotation(['W/W','tran','syn','XII,XIII']),\
            ('W/W', 'tran','syn','XII,XIII'))
        self.assertEqual(parse_annotation(['./W','cis','!1H(b_b)']),\
            ('./W', 'cis',None,'!1H(b_b)'))
        self.assertEqual(parse_annotation([]),\
            (None, None, None, None))
        self.assertRaises(RnaViewParseError, parse_annotation, ['X--X'])

    def test_parse_filename(self):
        """parse_filename: should return name of file"""
        lines = ["PDB data file name: pdb1t4l.ent_nmr.pdb"]
        self.assertEqual(parse_filename(lines), 'pdb1t4l.ent_nmr.pdb')
        lines = ["PDB data file name: pdb1t4l.ent_nmr.pdb","other line"]
        self.assertRaises(RnaViewParseError, parse_filename, lines)

    def test_parse_uncommon_residues(self):
        """parse_uncommon_residues: should fail on some missing residue info
        """
        lines = UC_LINES.split('\n')
        self.assertEqual(parse_uncommon_residues(lines),\
            {('D','16','TLN'):'u',('D','17','LCG'):'g',\
            ('0','2588','OMG'):'g',(' ','2621','PSU'):'P'})

        for l in UC_LINES_WRONG.split('\n'):
            self.assertRaises(RnaViewParseError, parse_uncommon_residues, [l])
    
    def test_parse_base_pairs_basic(self):
        """parse_base_pairs: basic input"""
        basic_lines =\
            ['25_437, 0:    34 C-G   448 0: +/+ cis         XIX',\
            '26_436, 0:    35 U-A   447 0: -/- cis         XX']
        bp1 = BasePair(Up=Base('0','34','C','25'),\
            Down=Base('0','448','G','437'),\
            Edges='+/+', Orientation='cis',Conformation=None,Saenger='XIX')
        bp2 = BasePair(Up=Base('0','35','U','26'),\
            Down=Base('0','447','A','436'),\
            Edges='-/-', Orientation='cis',Conformation=None,Saenger='XX')
        bps = BasePairs([bp1,bp2])
        obs = parse_base_pairs(basic_lines)
        for o,e in zip(obs,[bp1,bp2]):
            self.assertEqual(o,e)
        self.assertEqual(len(obs), 2)
    
        basic_lines =\
            ['25_437, 0:    34 c-P   448 0: +/+ cis         XIX',\
            '26_436, 0:    35 U-X   447 0: -/- cis         XX']
        self.assertRaises(RnaViewParseError, parse_base_pairs, basic_lines)
        
        basic_lines =\
            ['25_437, 0:    34 c-P   448 0: +/+ cis         XIX',\
            '26_436, 0:    35 I-A   447 0: -/- cis         XX']
        bp1 = BasePair(Up=Base('0','34','c','25'),\
            Down=Base('0','448','P','437'),\
            Edges='+/+', Orientation='cis',Conformation=None,Saenger='XIX')
        bp2 = BasePair(Up=Base('0','35','I','26'),\
            Down=Base('0','447','A','436'),\
            Edges='-/-', Orientation='cis',Conformation=None,Saenger='XX')
        bps = BasePairs([bp1,bp2])
        obs = parse_base_pairs(basic_lines)
        for o,e in zip(obs,[bp1,bp2]):
            self.assertEqual(o,e)
        self.assertEqual(len(obs), 2)

        
        lines = ['1_2,  :     6 G-G     7  :      stacked',\
            '1_16,  :     6 G-C    35  : +/+ cis         XIX']
        bp1 = BasePair(Up=Base(' ','6','G','1'),\
            Down=Base(' ','7','G','2'), Edges='stacked')
        bp2 = BasePair(Up=Base(' ','6','G','1'),\
            Down=Base(' ','35','C','16'),\
            Edges='+/+', Orientation='cis',Conformation=None,Saenger='XIX')
        obs = parse_base_pairs(lines)
        for o,e in zip(obs,[bp1,bp2]):
            self.assertEqual(o,e)
        self.assertEqual(len(obs), 2)
        
    def test_parse_base_multiplets_basic(self):
        """parse_base_multiplets: basic input"""
        basic_lines =\
            ['235_237_254_| [20 3]  0: 246 G  +  0: 248 A  +  0: 265 U',\
            '273_274_356_| [21 3]  0: 284 C  +  0: 285 A  +  0: 367 G']
        bm1 = BaseMultiplet([Base('0','246','G','235'),\
            Base('0','248','A','237'), Base('0','265','U','254')])
        bm2 = BaseMultiplet([Base('0','284','C','273'),\
            Base('0','285','A','274'), Base('0','367','G','356')])
        bms = BaseMultiplets([bm1,bm2])
        obs = parse_base_multiplets(basic_lines)
        for o,e in zip(obs,bms):
            for base_x, base_y in zip(o,e):
                self.assertEqual(base_x,base_y)
        self.assertEqual(len(obs), 2)
        self.assertEqual(len(obs[0]), 3)
    
        basic_lines =\
            ['235_237_254_| [20 3]  0: 246 G  +  0: 248 A  +  0: 265 I',\
            '273_274_356_| [21 3]  0: 284 P  +  0: 285 a  +  0: 367 G']
        bm1 = BaseMultiplet([Base('0','246','G','235'),\
            Base('0','248','A','237'), Base('0','265','I','254')])
        bm2 = BaseMultiplet([Base('0','284','P','273'),\
            Base('0','285','a','274'), Base('0','367','G','356')])
        bms = BaseMultiplets([bm1,bm2])
        obs = parse_base_multiplets(basic_lines)
        for o,e in zip(obs,bms):
            for base_x, base_y in zip(o,e):
                self.assertEqual(base_x,base_y)
        self.assertEqual(len(obs), 2)
        self.assertEqual(len(obs[0]), 3)

    def test_parse_base_multiplets_errors(self):
        """parse_base_multiplets: error checking"""
        # Unknown base
        basic_lines =\
            ['235_237_254_| [20 3]  0: 246 X  +  0: 248 A  +  0: 265 U',\
            '273_274_356_| [21 3]  0: 284 C  +  0: 285 A  +  0: 367 G']
        self.assertRaises(RnaViewParseError, parse_base_multiplets,\
            basic_lines)
        
        # number of rnaview_seqpos doesn't match number of bases
        basic_lines =\
            ['235_237_| [20 3]  0: 246 X  +  0: 248 A  +  0: 265 U',\
            '273_274_356_| [21 3]  0: 284 C  +  0: 285 A  +  0: 367 G']
        self.assertRaises(RnaViewParseError, parse_base_multiplets,\
            basic_lines)

        # Number of reported bases incorrect
        basic_lines =\
            ['235_237_254_| [20 3]  0: 246 X  +  0: 248 A  +  0: 265 U',\
            '273_274_356_| [21 5]  0: 284 C  +  0: 285 A  +  0: 367 G']
        self.assertRaises(RnaViewParseError, parse_base_multiplets,\
            basic_lines)

    def test_parse_number_of_pairs(self):
        """parse_number_of_pairs: good/bad input"""
        lines = ["The total base pairs =  31 (from   65 bases)"]
        exp = {'NUM_PAIRS':31, 'NUM_BASES':65}
        self.assertEqual(parse_number_of_pairs(lines), exp)
        
        lines = ["The total base pairs =  31 (from   65 bases)","XXX"]
        self.assertRaises(RnaViewParseError, parse_number_of_pairs, lines)

        lines = ["The total base pairs =  31(from   65 bases)"]
        self.assertRaises(RnaViewParseError, parse_number_of_pairs, lines)

    def test_parse_pair_counts(self):
        """parse_pair_counts: should work for even number of lines"""
        lines = PC_COUNTS1.split('\n')
        res = parse_pair_counts(lines)
        self.assertEqual(res['Standard'], 1)
        self.assertEqual(res['WS--cis'], 300)
        self.assertEqual(res['Bifurcated'], 2)
        self.assertEqual(res['HS-tran'], 0)

        lines = PC_COUNTS2.split('\n')
        res = parse_pair_counts(lines)
        self.assertEqual(res['Standard'], 19)
        self.assertEqual(res['WW-tran'], 1)
        self.assertEqual(res['HS-tran'], 0)
        self.failIf('Bifurcated' in res)

        self.assertRaises(RnaViewParseError, parse_pair_counts,\
            PC_COUNTS2.split('\n')[:-1])

        self.assertEqual(parse_pair_counts([]),{})

    def test_verify_bp_counts(self):
        """verify_bp_count: should raise an error if bp counts are wrong"""
        lines = RNAVIEW_PDB_REAL.split('\n')
        obs = RnaviewParser(lines)
        # this shouldn't raise an error
        verify_bp_counts(obs['BP'],11,obs['PC'])
        # reported number isn't right
        self.assertRaises(RnaViewParseError,\
            verify_bp_counts, obs['BP'], 12, obs['PC'])
        
        # No longer checks for the base pair counts reported in the 
        # dictionary, b/c this number doens't match the total when
        # modified bases are present.
        ## PREVIOUS TEST:
        # pair_counts isn't right
        #obs['PC']['Standard'] = 14
        #self.assertRaises(RnaViewParseError,\
        #    verify_bp_counts, obs['BP'], 11, obs['PC'])
        ## NEW TEST
        obs['PC']['Standard'] = 14
        verify_bp_counts(obs['BP'],11,obs['PC'])
        

    def test_MinimalRnaviewParser(self):
        """MinimalRnaviewParser: should divide lines into right classes"""
        exp = {'FN': ['PDB data file name: 1EHZ.pdb'], 'UC':\
        ['uncommon residue   I    1  on chain A [#1] assigned to: I',
        'uncommon residue 2MG   10  on chain A [#10] assigned to: g'],
        'BP':['1_72, A:     1 I-C    72 A: X/X cis         n/a',
        '58_60, A:    58 a-C    60 A: S/S tran   syn    !(s_s)'],
        'BM':['9_12_23_| [1 3]  A: 9 A  +  A: 12 U  +  A: 23 A',
        '13_22_46_| [2 3]  A: 13 C  +  A: 22 G  +  A: 46 g'],
        'PC':['Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran',
        '19        3        1        0        1        0        0',
        'WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran',
        '0        3        0        2        0        0'],
        'NP':['The total base pairs =  30 (from   76 bases)']}
        obs = MinimalRnaviewParser(RNAVIEW_LINES.split('\n'))
        self.assertEqual(len(obs), len(exp))
        self.assertEqual(obs, exp)

    def test_MinimalRnaviewParser_short(self):
        """MinimalRnaviewparser: should leave lists empty if no lines found"""
        lines = RNAVIEW_LINES_SHORT.split('\n')
        res = MinimalRnaviewParser(lines)
        self.assertEqual(len(res['FN']), 1)
        self.assertEqual(res['UC'], [])
        self.assertEqual(res['BM'], [])
        self.assertEqual(len(res['PC']), 4)
        self.assertEqual(len(res['BP']), 11)
        self.assertEqual(len(res['NP']), 1)

    def test_RnaviewParser(self):
        """RnaviewParser: should work with/without model and/or verification
        """
        rnaview_lines = RNAVIEW_PDB_REAL.split('\n')
        
        obs = RnaviewParser(rnaview_lines)
        self.assertEqual(obs['FN'], 'pdb430d.ent')
        self.assertEqual(len(obs['UC']),1)
        self.assertEqual(len(obs['BP']),19)
        self.assertEqual(len(obs['BM']),0)
        self.assertEqual(obs['BM'],BaseMultiplets())
        self.assertEqual(obs['PC']['Standard'],7)
        self.assertEqual(obs['BP'][2].Down.ResName,'c')
        self.assertEqual(obs['BP'][6].Edges,'stacked')
        self.assertEqual(obs['NP']['NUM_PAIRS'], 11)
        self.assertEqual(obs['NP']['NUM_BASES'], 29)

    def test_RnaviewParser_error(self):
        """RnaviewParser: strict or not"""
        lines = RNAVIEW_LINES_ERROR.split('\n')
        self.assertRaises(RnaViewParseError, RnaviewParser, lines, strict=True)

        obs = RnaviewParser(lines, strict=False)
        self.assertEqual(obs['NP'], None)
        self.assertEqual(obs['BP'][1].Up.ResId, '2')
        self.assertEqual(obs['PC']['Standard'], 6)


RNAVIEW_LINES=\
"""PDB data file name: 1EHZ.pdb
uncommon residue   I    1  on chain A [#1] assigned to: I
uncommon residue 2MG   10  on chain A [#10] assigned to: g
BEGIN_base-pair
     1_72, A:     1 I-C    72 A: X/X cis         n/a
    58_60, A:    58 a-C    60 A: S/S tran   syn    !(s_s)
END_base-pair
Summary of triplets and higher multiplets
BEGIN_multiplets
9_12_23_| [1 3]  A: 9 A  +  A: 12 U  +  A: 23 A
13_22_46_| [2 3]  A: 13 C  +  A: 22 G  +  A: 46 g
END_multiplets
  The total base pairs =  30 (from   76 bases)
------------------------------------------------
 Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran
       19        3        1        0        1        0        0
  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran
        0        3        0        2        0        0
------------------------------------------------
"""

RNAVIEW_LINES_SHORT=\
"""PDB data file name: pdb17ra.ent_nmr.pdb
BEGIN_base-pair
     1_21,  :     1 G-C    21  : +/+ cis         XIX
     2_20,  :     2 G-C    20  : +/+ cis         XIX
     3_19,  :     3 C-G    19  : +/+ cis         XIX
     4_18,  :     4 G-U    18  : W/W cis         XXVIII
      5_6,  :     5 U-A     6  :      stacked
     5_17,  :     5 U-A    17  : -/- cis         XX
     7_16,  :     7 A-U    16  : W/W cis         n/a
     8_15,  :     8 G-C    15  : +/+ cis         XIX
     9_14,  :     9 G-C    14  : +/+ cis         XIX
    10_14,  :    10 A-C    14  :      stacked
    11_13,  :    11 U-A    13  : S/W tran        n/a
END_base-pair
  The total base pairs =   9 (from   21 bases)
------------------------------------------------
 Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran
        6        2        0        0        0        0        0
  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran
        0        0        0        1        0        0
------------------------------------------------"""

PC_COUNTS1=\
""" Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran
        1        0        0        0        0        0        0
  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran
       12        0      300        0        0        0
Single-bond  Bifurcated 
        0        2"""
PC_COUNTS2=\
""" Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran
       19        3        1        0        1        0        0
  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran
        0        3        0        2        0        0"""

UC_LINES=\
"""uncommon residue TLN   16  on chain D [#16] assigned to: u
uncommon residue LCG   17  on chain D [#17] assigned to: g
uncommon residue OMG 2588  on chain 0 [#2430] assigned to: g
uncommon residue PSU 2621  on chain   [#2463] assigned to: P"""

UC_LINES_WRONG=\
"""uncommon residue       16  on chain D [#16] assigned to: u
uncommon residue LCG       on chain D [#17] assigned to: g
uncommon residue OMG 2588  on chain 0 [#2430] assigned to: """

RNAVIEW_LINES_TOTAL=\
"""PDB data file name: 1EHZ.pdb
uncommon residue PSU    1  on chain A [#1] assigned to: P
uncommon residue 2MG   10  on chain A [#10] assigned to: g
BEGIN_base-pair
     1_72, A:     1 P-C    20 A: +/+ cis         n/a
     1_72, A:     2 A-U    19 A: H/W cis         n/a
    58_60, A:    10 a-U    60 A: S/S tran   syn    !(s_s)
END_base-pair
Summary of triplets and higher multiplets
BEGIN_multiplets
9_12_23_| [1 3]  A: 1 P  +  A:  2 U  +  A: 60 U
13_22_46_| [2 3]  A: 20 C  +  A: 19 U  +  A: 60 U
END_multiplets
  The total base pairs =  30 (from   76 bases)
------------------------------------------------
 Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran
       19        3        1        0        1        0        0
  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran
        0        3        0        2        0        0
------------------------------------------------
"""

RNAVIEW_PDB_REAL=\
"""PDB data file name: pdb430d.ent
uncommon residue  +C   27  on chain A [#27] assigned to: c
BEGIN_base-pair
     1_29, A:     1 G-C    29 A: +/+ cis         XIX
     2_28, A:     2 G-C    28 A: +/+ cis         XIX
     3_27, A:     3 G-c    27 A: +/+ cis         XIX
     4_26, A:     4 U-A    26 A: -/- cis         XX
     5_25, A:     5 G-C    25 A: +/+ cis         XIX
     6_24, A:     6 C-G    24 A: +/+ cis         XIX
      7_8, A:     7 U-C     8 A:      stacked
      8_9, A:     8 C-A     9 A:      stacked
     9_21, A:     9 A-A    21 A: H/H tran        II
    11_20, A:    11 U-A    20 A: W/H tran        XXIV
    12_19, A:    12 A-G    19 A: H/S tran        XI
    13_18, A:    13 C-G    18 A: +/+ cis         XIX
    14_17, A:    14 G-A    17 A: S/H tran        XI
    25_26, A:    25 C-A    26 A:      stacked
     7_23, A:     7 U-C    23 A: W/W tran        !1H(b_b)
     8_23, A:     8 C-C    23 A: S/H cis         !1H(b_b).
    10_11, A:    10 G-U    11 A: S/H cis         !1H(b_b)
    11_21, A:    11 U-A    21 A: W/H tran        !1H(b_b)
    10_20, A:    10 G-A    20 A: W/. tran        !(s_s)
END_base-pair
  The total base pairs =  11 (from   29 bases)
------------------------------------------------
 Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran
        7        0        0        0        1        0        0
  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran
        0        1        0        0        0        2
------------------------------------------------"""

RNAVIEW_LINES_ERROR=\
"""PDB data file name: pdb17ra.ent_nmr.pdb
BEGIN_base-pair
     1_21,  :     1 G-C    21  : +/+ cis         XIX
     2_20,  :     2 G-C    20  : +/+ cis         XIX
     3_19,  :     3 C-G    19  : +/+ cis         XIX
     4_18,  :     4 G-U    18  : W/W cis         XXVIII
      5_6,  :     5 U-A     6  :      stacked
     5_17,  :     5 U-A    17  : -/- cis         XX
     7_16,  :     7 A-U    16  : W/W cis         n/a
     8_15,  :     8 G-C    15  : +/+ cis         XIX
     9_14,  :     9 G-C    14  : +/+ cis         XIX
    10_14,  :    10 A-C    14  :      stacked
    11_13,  :    11 U-A    13  : S/W tran        n/a
END_base-pair
  The total base pairs =   9(from   21 bases)
------------------------------------------------
 Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran
        6        2        0        0        0        0        0
  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran
        0        0        0        1        0        0
------------------------------------------------"""


#run if called from command-line
if __name__ == "__main__":
    main()
