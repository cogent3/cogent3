#!/usr/bin/env python
"""Tests for ViennaStructure and related classes.
"""
from cogent.util.unit_test import TestCase, main
from cogent.struct.rna2d import ViennaStructure, Vienna, Pairs,\
    Partners, EmptyPartners, WussStructure, wuss_to_vienna, StructureNode, \
    Stem, classify, PairError

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class RnaAlphabet(object):
    Pairs = {
    ('A','U'): True,    #True vs False for 'always' vs 'sometimes' pairing
    ('C','G'): True,
    ('G','C'): True,
    ('U','A'): True,
    ('G','U'): False,
    ('U','G'): False,
}

class Rna(str):
    Alphabet = RnaAlphabet

class StemTests(TestCase):
    """Tests for the Stem object"""

    def test_init_empty(self):
        """Stem should init ok with no parameters."""
        s = Stem()
        self.assertEqual((s.Start, s.End, s.Length), (None, None, 0))

    def test_init(self):
        """Stem should init with Start, End, and Length"""
        s = Stem(Length=3)
        self.assertEqual((s.Start, s.End, s.Length), (None, None, 3))
        #should set Length to 0 if not supplied and unpaired
        s = Stem(Start=3)
        self.assertEqual((s.Start, s.End, s.Length), (3, None, 0))
        s = Stem(End=3)
        self.assertEqual((s.Start, s.End, s.Length), (None, 3, 0))
        #should set Length to 1 if not supplied and paired
        s = Stem(Start=3, End=5)
        self.assertEqual((s.Start, s.End, s.Length), (3, 5, 1))
        #parameters should be in order Start, End, Length
        #note that you're allowed to initialize an invalid stem, like the
        #following one (can't have 7 pairs between 3 and 5); this is often
        #useful when building up a tree that you plan to renumber().
        s = Stem(3, 5, 7)
        self.assertEqual((s.Start, s.End, s.Length), (3, 5, 7))
        #not allowed more than 3 parameters
        self.assertRaises(TypeError, Stem, 1, 2, 3, 4)
    
    def test_len(self):
        """Stem len() should return self.Length"""
        s = Stem()
        self.assertEqual(len(s), 0)
        s.Length = 5
        self.assertEqual(len(s), 5)
        s.Length = None
        self.assertRaises(TypeError, len, s)

    def test_getitem(self):
        """Stem getitem should return a Stem object for the ith pair in the stem"""
        s = Stem()
        self.assertRaises(IndexError, s.__getitem__, 0)
        s.Start = 5
        s.End = 8
        s.Length = 1
        pairs = list(s)
        self.assertEqual(pairs, [Stem(5, 8, 1)])
        s.Length = 2
        pairs = list(s)
        self.assertEqual(pairs, [Stem(5,8,1),Stem(6,7,1)])
        #WARNING: Stem will not complain when iterating over an invalid helix,
        #as per the one below
        s.Length = 5
        pairs = list(s)
        self.assertEqual(pairs, [Stem(5,8,1),Stem(6,7,1),Stem(7,6,1),\
            Stem(8,5,1), Stem(9,4,1)])

    def test_cmp(self):
        """Stems should compare equal when the data is the same"""
        self.assertEqual(Stem(1,2,3), Stem(1,2,3))
        self.assertNotEqual(Stem(1,2,5), Stem(1,2,3))

        a = Stem(1, 10, 2)
        b = Stem(2, 8, 1)
        c = Stem(15, 20, 2)
        l = [c, b, a]
        l.sort()
        self.assertEqual(l, [a,b,c])

    def test_extract(self):
        """Stems extract should return a list of 'pairs' from a sequence"""
        seq = 'UGAGAUUUUCU'
        s = Stem(1, 10, 3)
        self.assertEqual(s.extract(seq), [('G','U'),('A','C'),('G','U')])
        s = Stem(0, 1)
        self.assertEqual(s.extract(seq), [('U','G')])
        #should put in None if either position hasn't been specified
        s = Stem(5)
        self.assertEqual(s.extract(seq), [('U', None)])
        s = Stem()
        self.assertEqual(s.extract(seq), [(None, None)])
        #should raise IndexError if the stem contains bases outside the seq
        s = Stem(50, 60, 5)
        self.assertRaises(IndexError, s.extract, seq)

    def test_hash(self):
        """Stems hash should allow use as dict keys if unchanged"""
        #WARNING: if you change the Stem after putting it in a dict, all bets
        #are off as to behavior. Don't do it!
        s = Stem(1, 5, 2)
        t = Stem(1, 5, 2)
        u = Stem(2, 4, 6)
        v = Stem(2, 4, 6)
        w = Stem(2, 4, 4)
        d = {}

        assert s is not t
        
        for i in (s, t, u, v, w):
            if i in d:
                d[i] += 1
            else:
                d[i] = 1

        self.assertEqual(len(d), 3)
        self.assertEqual(d[Stem(1, 5, 2)], 2)
        self.assertEqual(d[Stem(2, 4, 6)], 2)
        self.assertEqual(d[Stem(2, 4, 4)], 1)
        assert Stem(1,5) not in d

    def test_str(self):
        """Stem str should print Start, End and Length as a tuple"""
        self.assertEqual(str(Stem()), '(None,None,0)')
        self.assertEqual(str(Stem(3)), '(3,None,0)')
        self.assertEqual(str(Stem(3,4)), '(3,4,1)')
        self.assertEqual(str(Stem(3,4,5)), '(3,4,5)')

    def test_nonzero(self):
        """Stem nonzero should return True if paired (length > 0)"""
        assert not Stem()
        assert not Stem(1)
        assert Stem(7, 10)
        assert Stem(1, 7, 1)
        assert Stem(2, 8, 3)
        #go strictly by length; don't check if data is invalid
        assert Stem(0, 0)
        assert Stem(5, None, 10)
        assert not Stem(5, 7, -1)
        
class PartnersTests(TestCase):
    """Tests for Partners object"""

    def test_init(self):
        """Partners should init with empty list and stay free of conflicts"""
        self.assertEqual(Partners([]),[])
        empty = Partners([None]*6)
        self.assertEqual(empty,[None,None,None,None,None,None])
        self.assertRaises(ValueError,empty.__setitem__,2,2) 
        empty[2] = 3
        self.assertEqual(empty,[None,None,3,2,None,None])
        empty[3] = 1
        self.assertEqual(empty,[None,3,None,1,None,None])
        empty[3] = 5
        self.assertEqual(empty,[None,None,None,5,None,3])
        empty[1] = None
        self.assertEqual(empty,[None,None,None,5,None,3])
    
    def test_toPairs(self):
        """Partners toPairs() should return a Pairs object"""
        p = Partners([None,3,None,1,5,4])
        self.assertEqualItems(p.toPairs(),[(1,3),(4,5)])
        assert isinstance(p.toPairs(),Pairs)
        self.assertEqual(Partners([None]*10).toPairs(),[])

    def test_not_implemented(self):
        """Partners not_implemented should raise error for 'naughty' methods"""
        p = Partners([None,3,1,5,4])
        self.assertRaises(NotImplementedError,p.pop)
        self.assertRaises(NotImplementedError,p.sort)
        self.assertRaises(NotImplementedError,p.__delitem__,3)
        
class PairsTests(TestCase):
    """Tests for Pairs object"""

    def setUp(self):
        """Pairs SetUp method for all tests"""
        self.Empty = Pairs([])
        self.OneList = Pairs([[1,2]])
        self.OneTuple = Pairs([(1,2)])
        self.MoreLists = Pairs([[2,4],[3,9],[6,36],[7,49]])
        self.MoreTuples = Pairs([(2,4),(3,9),(6,36),(7,49)])
        self.MulNoOverlap = Pairs([(1,10),(2,9),(3,7),(4,12)])
        self.MulOverlap = Pairs([(1,2),(2,3)])
        self.Doubles = Pairs([[1,2],[1,2],[2,3],[1,3]])
        self.Undirected = Pairs([(2,1),(6,4),(1,7),(8,3)])
        self.UndirectedNone = Pairs([(5,None),(None,3)])
        self.UndirectedDouble = Pairs([(2,1),(1,2)])
    
        self.NoPseudo = Pairs([(1,20),(2,19),(3,7),(4,6),(10,15),(11,14)])
        self.NoPseudo2 = Pairs([(1,3),(4,6)])
        #((.(.)).)
        self.p0 = Pairs([(0,6),(1,5),(3,8)])
        #(.((..(.).).))
        self.p1 = Pairs([(0,9),(2,12),(3,10),(5,7)])
        #((.(.(.).)).)
        self.p2 = Pairs([(0,10),(1,9),(3,12),(5,7)])
        #((.((.(.)).).))
        self.p3 = Pairs([(0,9),(1,8),(3,14),(4,13),(6,11)])
        #(.(((.((.))).)).(((.((((..))).)))).)
        self.p4 = Pairs([(0,35),(2,11),(3,10),(4,9),(6,14),(7,13),(16,28),\
            (17,27),(18,26),(20,33),(21,32),(22,31),(23,30)])
        #(.((.).))
        self.p5 = Pairs([(0,5),(2,8),(3,7)])
        self.p6 = Pairs([(0,19),(2,6),(3,5),(8,14),(9,13),(10,12),\
            (16,22),(17,21)])
        self.p7 = Pairs([(0,20),(2,6),(3,5),(8,14),(9,10),(11,16),(12,15),\
            (17,23),(18,22)])

         
    def test_init(self):
        """Pairs should initalize with both lists and tuples"""
        self.assertEqual(self.Empty,[])
        self.assertEqual(self.OneList,[[1,2]])
        self.assertEqual(self.OneTuple,[(1,2)])
        self.assertEqual(self.MulNoOverlap,[(1,10),(2,9),(3,7),(4,12)])
        self.assertEqual(self.MulOverlap,[(1,2),(2,3)])

    def test_toPartners(self):
        """Pairs toPartners() should return a Partners object"""
        a = Pairs([(1,5),(3,4),(6,9),(7,8)]) #normal
        b = Pairs([(0,4),(2,6)]) #pseudoknot
        c = Pairs([(1,6),(3,6),(4,5)]) #conflict

        self.assertEqual(a.toPartners(10),[None,5,None,4,3,1,9,8,7,6])
        self.assertEqual(a.toPartners(13,3),\
        [None,None,None,None,8,None,7,6,4,12,11,10,9])
        assert isinstance(a.toPartners(10),Partners)
        self.assertEqual(b.toPartners(7),[4,None,6,None,0,None,2])
        self.assertRaises(ValueError,c.toPartners,7)
        self.assertEqual(c.toPartners(7,strict=False),[None,None,None,6,5,4,3])

        #raises an error when try to insert something at non-existing indices
        self.assertRaises(IndexError,c.toPartners,0)

    def test_toVienna(self):
        """Pairs toVienna() should return a ViennaStructure if possible"""
        a = Pairs([(1,5),(3,4),(6,9),(7,8)]) #normal
        b = Pairs([(0,4),(2,6)]) #pseudoknot
        c = Pairs([(1,6),(3,6),(4,5)]) #conflict
        d = Pairs([(1,6),(3,None)])
        e = Pairs([(1,9),(8,2),(7,3)]) #not directed
        f = Pairs([(1,6),(2,5),(10,15),(14,11)]) # not directed

        self.assertEqual(a.toVienna(10),'.(.())(())')
        self.assertEqual(a.toVienna(13,offset=3),'....(.())(())')
        
        self.assertRaises(PairError,b.toVienna,7) #pseudoknot NOT accepted
        self.assertRaises(Exception,b.toVienna,7) #old test for exception
        self.assertRaises(ValueError,c.toVienna,7)
        
        #pairs containging None are being skipped
        self.assertEquals(d.toVienna(7),'.(....)')
        
        #raises error when trying to insert at non-existing indices
        self.assertRaises(IndexError,a.toVienna,3)

        self.assertEqual(Pairs().toVienna(3),'...')
        
        #test when parsing in the sequence
        self.assertEqual(a.toVienna('ACGUAGCUAG'),'.(.())(())')
        self.assertEqual(a.toVienna(Rna('AACCGGUUAGCUA'), offset=3),\
            '....(.())(())')
       
        self.assertEqual(e.toVienna(10),'.(((...)))')
        self.assertEqual(f.toVienna(20),'.((..))...((..))....')

    def test_tuples(self):
        """Pairs tuples() should transform the elements of list to tuples"""
        x = Pairs([])
        x.tuples()
        assert x == []
        
        x = Pairs([[1,2],[3,4]])
        x.tuples()
        assert x == [(1,2),(3,4)]
        
        x = Pairs([(1,2),(3,4)])
        x.tuples()
        assert x == [(1,2),(3,4)]
        assert x != [[1,2],[3,4]]

    def test_unique(self):
        """Pairs unique() should remove double occurences of certain tuples"""
        self.assertEqual(self.Empty.unique(),[])
        self.assertEqual(self.MoreTuples.unique(),self.MoreTuples)
        self.assertEqual(self.Doubles.unique(),Pairs([(1,2),(2,3),(1,3)]))

    def test_directed(self):
        """Pairs directed() should change all pairs so that a<b in (a,b)"""
        self.assertEqual(self.Empty.directed(),[])
        res = self.Undirected.directed()
        res.sort()
        self.assertEqual(res,Pairs([(1,2),(1,7),(3,8),(4,6)]))
        res = self.UndirectedNone.directed()
        self.assertEqual(res,Pairs([]))
        res = self.UndirectedDouble.directed()
        self.assertEqual(res,Pairs([(1,2)]))

    def test_symmetric(self):
        """Pairs symmetric() should add (down,up) for each (up,down)"""
        self.assertEqual(self.Empty.symmetric(),[])
        self.assertEqualItems(self.OneTuple.symmetric(),[(2,1),(1,2)])
        self.assertEqualItems(Pairs([(1,2),(1,2)]).symmetric(),[(1,2),(2,1)])
        self.assertEqualItems(Pairs([(1,2),(3,4)]).symmetric(),\
        [(1,2),(2,1),(3,4),(4,3)])
        self.assertEqualItems(Pairs([(1,None)]).symmetric(),[])

    def test_paired(self):
        """Pairs paired() should omit all pairs containing None"""
        self.assertEqual(self.Empty.paired(),[])
        self.assertEqual(Pairs([(1,2),(2,None),(None,3),(None,None)]).paired()\
        ,[(1,2)])

    def test_hasConflicts(self):
        """Pairs hasConflicts() should return True if there are conflicts"""
        assert not self.Empty.hasConflicts()
        assert not Pairs([(1,2),(3,4)]).hasConflicts()
        assert Pairs([(1,2),(2,3)]).hasConflicts()
        assert Pairs([(1,2),(2,None)]).hasConflicts()

    def test_mismatches(self):
        """Pairs mismatches() should return #pairs that can't be formed"""
        # with plain string
        self.assertEqual(Pairs([(0,1)]).mismatches('AC',{}),1)
        self.assertEqual(Pairs([(0,1)]).mismatches('AC',{('A','C'):None}),0)
        self.assertEqual(Pairs([(0,1)]).mismatches('AC',{('A','G'):None}),1)
        self.assertEqual(Pairs([(0,1),(2,3),(3,1)]).\
        mismatches('ACGU',{('A','U'):None}),3)

        # using sequence with alphabet
        sequence = Rna('ACGUA')
        self.assertEqual(Pairs([(0,1),(0,4),(0,3)]).mismatches(sequence),2)

    def test_hasPseudoknots(self):
        """Pairs hasPseudoknots() should return True if there's a pseudoknot"""
                
        assert not self.NoPseudo.hasPseudoknots()
        assert not self.NoPseudo2.hasPseudoknots()
        #add tests for ((.))() etc
        assert self.p0.hasPseudoknots()
        assert self.p1.hasPseudoknots() 
        assert self.p2.hasPseudoknots()
        assert self.p3.hasPseudoknots()
        assert self.p4.hasPseudoknots()
        assert self.p5.hasPseudoknots()
        assert self.p6.hasPseudoknots()
        assert self.p7.hasPseudoknots()

        
class StructureStringInitTests(TestCase):
    """Tests for initializing structures"""
    
    def setUp(self):
        self.Struct = ViennaStructure
        self.Empty = ''
        self.NoPairs = '.....'
        self.OneHelix = '((((()))))'
        self.ManyHelices = '(..(((...)).((.(((((..))).)))..((((..))))))...)'
        self.Ends = '..(.)..'   #has trailing bases at ends
        self.Internal = '(((...)))..(()).' #has internal non-nested region
        self.TooManyOpen = '(((...))..'
        self.TooManyClosed = '((...)))..'
        self.Invalid = 'fdjklgfj'
    
    def test_init_empty(self):
        """StructureString initialization should accept an empty string"""
        self.assertEqual(self.Empty, str(self.Struct(self.Empty)))
        
    def test_init_no_pairs(self):
        """StructureString should allow a structure with no pairs"""
        self.assertEqual(self.NoPairs, str(self.Struct(self.NoPairs)))

    def test_init_one_helix(self):
        """StructureString should allow a structure with one helix"""
        self.assertEqual(self.OneHelix, str(self.Struct(self.OneHelix)))

    def test_init_long_nested(self):
        """StructureString should accept long structure w/ many nested pairs"""
        self.assertEqual(self.ManyHelices,str(self.Struct(self.ManyHelices)))

    def test_init_too_many_pairs_opened(self):
        """StructureString should raise IndexError if too many open pairs"""
        self.assertRaises(IndexError, self.Struct, self.TooManyOpen)

    def test_init_too_few_pairs_opened(self):
        """StructureString should raise IndexError if too few open pairs"""
        self.assertRaises(IndexError, self.Struct, self.TooManyClosed)

    def test_init_invalid_string(self):
        """StructureString should raise ValueError if invalid chars"""
        self.assertRaises(ValueError, self.Struct, self.Invalid)


class WussStructureInitTests(StructureStringInitTests):
    """Test for initializing WussStructures"""
    def setUp(self):
        self.Struct = WussStructure
        self.Empty = ''
        self.NoPairs = '__--_'
        self.OneHelix = '[{(<()>)}]'
        self.ManyHelices = '<,,(({___})~((:((({{__}}),))),,<<[(__)]>>)):::>'
        self.TooManyOpen = '{<[___]>::'
        self.TooManyClosed = '<<___>>>::'
        self.Invalid = '__!<><>'


class StructureStringTests(TestCase):
    """Test that StructureString methods and properties work.
    
    By default tested on a ViennaStructure. These tests should
    work for every StructureString.
    """
    Struct = ViennaStructure

    def setUp(self):
        """SetUp for all StructureString tests"""
        self.Energy1 = 0.0
        self.Energy2 = -1e-02
        self.NoPairsStr = '.....'
        self.NoPairs = self.Struct('.....',self.Energy1)
        self.OneHelixStr = '((((()))))'
        self.OneHelix = self.Struct('((((()))))',self.Energy2)
        self.TwoHelixStr = '((.))(()).'
        self.TwoHelix = self.Struct('((.))(()).',self.Energy1)
        self.ThreeHelixStr = '(((((..))..(((..)))..)))'
        self.ThreeHelix = self.Struct('(((((..))..(((..)))..)))')
        self.EmptyStr = ''
        self.ManyHelicesStr = '(..(((...)).((.(((((..))).)))..((((..))))))...)'
        self.EndsStr = '..(.)..'   #has trailing bases at ends
        self.InternalStr = '(((...)))..(()).' #has internal non-nested region
    
    def test_len(self):
        """StructureString len() should match structure length"""
        self.assertEqual(len(self.TwoHelix), 10)
        #Do you want the init possible with Vienna() for empty str???
        self.assertEqual(len(Vienna('')), 0)

    def test_getitem(self):
        """StructureString struct[index] should return char at index in struct"""
        self.assertEqual(self.NoPairs[0], self.NoPairsStr[0])  #dot
        self.assertEqual(self.OneHelix[0], self.OneHelixStr[0]) #close pair
        #negative index, end
        self.assertEqual(self.OneHelix[-1], self.OneHelixStr[-1])
        #middle of sequence
        self.assertEqual(self.TwoHelix[4], self.TwoHelixStr[4]) 
                
    def test_getslice(self):
        """StructureString struct[a:b] should return slice from a to b"""
        self.assertEqual(self.OneHelix[0:3], self.OneHelixStr[0:3])
        self.assertEqual(self.OneHelix[0:0], self.OneHelixStr[0:0])
        self.assertEqual(self.OneHelix[:], self.OneHelixStr[:])
        self.assertEqual(self.OneHelix[4:7], self.OneHelixStr[4:7])
        self.assertEqual(self.OneHelix[5:], self.OneHelixStr[5:])
        
    def test_str(self):
        """StructureString str() should print structure and energy, if known"""
        self.assertEqual(str(self.NoPairs),  '..... (0.0)')
        self.assertEqual(str(self.OneHelix), '((((())))) (-0.01)')
        self.assertEqual(str(self.TwoHelix), '((.))(()). (0.0)')
        self.assertEqual(str(self.ThreeHelix), '(((((..))..(((..)))..)))')
        
    def test_toPartners(self):
        """StructureString toPartners() should return Partners object"""
        self.assertEqual(self.NoPairs.toPartners(), [None]*5)
        self.assertEqual(self.OneHelix.toPartners(),[9,8,7,6,5,4,3,2,1,0])
        self.assertEqual(self.TwoHelix.toPartners(),\
        [4,3,None,1,0,8,7,6,5,None])
        self.assertEqual(self.ThreeHelix.toPartners(),[23,22,21,8,7,None,\
        None,4,3,None,None,18,17,16,None,None,13,12,11,None,None,2,1,0])

    def test_toPairs(self):
        """StructureString toPairs() should return Pairs object"""
        self.assertEqual(self.NoPairs.toPairs(), [])
        self.assertEqualItems(self.OneHelix.toPairs(), \
        [(0,9),(1,8),(2,7),(3,6),(4,5)])
        self.assertEqualItems(self.TwoHelix.toPairs(), \
        [(0,4),(1,3),(5,8),(6,7)])
        self.assertEqualItems(self.ThreeHelix.toPairs(), \
        [(0,23),(1,22),(2,21),(3,8),(4,7),(11,18),(12,17),(13,16)])
    
    def test_toTree(self):
        """StructureString toTree should produce correct tree"""
        for struct in [self.EmptyStr, self.NoPairsStr, self.OneHelixStr, \
            self.ManyHelicesStr, self.EndsStr, self.InternalStr]:
            self.assertEqual(str(self.Struct(struct).toTree()), struct)

        test = self.Struct(self.ManyHelicesStr).toTree()
        n = test[0]
        self.assertEqual(n.Start, 0)
        self.assertEqual(n.End, 46)
        n = test[0][2][0][0][1]
        self.assertEqual(n.Start, 7)
        self.assertEqual(n.End, None)
        n = test[0][2][1]
        self.assertEqual(n.Start, 11)
        self.assertEqual(n.End, None)
        n = test[0][2][2][0][1][0][0][0][0][1]
        self.assertEqual(n.Start, 21)
        self.assertEqual(n.End, None)

        test = self.Struct(self.InternalStr).toTree()
        self.assertEqual([i.Start for i in test], [0, 9, 10, 11, 15])
        self.assertEqual([i.End for i in test], [8, None, None, 14, None])

    
class WussStructureTests(StructureStringTests):
    """Test that WussStructure methods and properties work"""
    def setUp(self):
        """Define three standard structures: note differences in whitespace
        and number formatting"""
        self.Struct = WussStructure
        super(WussStructureTests,self).setUp()
        self.WussNoPairs = WussStructure('.-_,~:')
        self.WussOneHelix = WussStructure('[-{<(__)-->}]',-0.01)
        self.WussTwoHelix = WussStructure('{[.]}(<>).',1.11)
        self.WussThreeHelix = WussStructure('::(<<({__}),,([(__)])-->>)')
        self.WussPseudo = WussStructure('<<__AA>>_aa::')
        
    def test_wuss_toPairs(self):
        """WussStructure toPairs() should return a valid Pairs object"""
        self.assertEqual(self.WussNoPairs.toPairs(),[])
        self.assertEqualItems(self.WussOneHelix.toPairs(),\
        [(0,12),(2,11),(3,10),(4,7)])
        self.assertEqualItems(self.WussTwoHelix.toPairs(),\
        [(0,4),(1,3),(5,8),(6,7)])
        self.assertEqualItems(self.WussThreeHelix.toPairs(),\
        [(2,25),(3,24),(4,23),(5,10),(6,9),(13,20),(14,19),(15,18)])
        self.assertEqualItems(self.WussPseudo.toPairs(),\
        [(0,7),(1,6)])
    
    def test_wuss_toPartners(self):
        """WussStructure toPartners() should return valid Partners object"""
        self.assertEqual(self.WussNoPairs.toPartners(),[None]*6)
        self.assertEqualItems(self.WussThreeHelix.toPartners(),\
        [None,None,25,24,23,10,9,None,None,6,5,None,None,20,19,\
        18,None,None,15,14,13,None,None,4,3,2])
        self.assertEqualItems(self.WussPseudo.toPartners(),\
        [7,6,None,None,None,None,1,0,None,None,None,None,None])


class Rna2dTests(TestCase):

    def test_Vienna(self):
        """Vienna should initalize from several formats"""
        self.NoPairs = Vienna('.......... (0.0)')
        self.OneHelix = Vienna('((((()))))    (-1e-02)')
        self.TwoHelix = Vienna('((.))(()). \t(1.11)')
        self.ThreeHelix = Vienna('(((((..))..(((..)))..)))')
        self.GivenEnergy = Vienna('((.))',0.1)
        self.TwoEnergies = Vienna('((.)) (4.6)',2.1)
        
        self.assertEqual(self.NoPairs, '..........')
        self.assertEqual(self.NoPairs.Energy, 0.0)
        self.assertEqual(self.OneHelix, '((((()))))')
        self.assertEqual(self.OneHelix.Energy, -1e-2)
        self.assertEqual(self.TwoHelix, '((.))(()).')
        self.assertEqual(self.TwoHelix.Energy, 1.11)
        self.assertEqual(self.ThreeHelix, '(((((..))..(((..)))..)))')
        self.assertEqual(self.ThreeHelix.Energy, None)
        self.assertEqual(self.GivenEnergy.Energy,0.1)
        self.assertEqual(self.TwoEnergies.Energy,2.1)

    
    def test_EmptyPartners(self):
        """EmptyPartners should return list of 'None's of given length"""
        self.assertEqual(EmptyPartners(0),[])
        self.assertEqual(EmptyPartners(1),[None])
        self.assertEqual(EmptyPartners(10),[None]*10)

    
    def test_wuss_to_vienna(self):
        """wuss_to_vienna() should transform Wuss into Vienna"""
        empty = WussStructure('.....')
        normal = WussStructure('[.{[<...}}}}')
        pseudo = WussStructure('[..AA..]..aa')
        self.assertEqual(wuss_to_vienna(normal),'(.(((...))))')
        self.assertEqual(wuss_to_vienna(empty),'.....')
        self.assertEqual(wuss_to_vienna(pseudo),'(......)....')

    def test_classify(self):
        """classify() should classify valid structures correctly"""
        Empty = '' 
        NoPairs = '.....'
        OneHelix = '((((()))))'
        ManyHelices = '(..(((...)).((.(((((..))).)))..((((..))))))...)'
        Ends = '..(.)..'
        FirstEnd = '..((()))'
        LastEnd = '((..((.))))...'
        Internal = '(((...)))..((.)).'
        #following structure is from p 25 of Eddy's WUSS description manual
        Eddy = '..((((.(((...)))...((.((....))..)).)).))'

        structs = [Empty, NoPairs, OneHelix, ManyHelices, Ends, \
            FirstEnd, LastEnd, Internal, Eddy]

        EmptyResult = '' 
        NoPairsResult = 'EEEEE'
        OneHelixResult = 'SSSSSSSSSS'
        ManyHelicesResult = 'SBBSSSLLLSSJSSBSSSSSLLSSSBSSSJJSSSSLLSSSSSSBBBS'
        EndsResult = 'EESLSEE'
        FirstEndResult = 'EESSSSSS'
        LastEndResult = 'SSBBSSLSSSSEEE'
        InternalResult = 'SSSLLLSSSFFSSLSSE'
        #following structure is from p 25 of Eddy's WUSS description manual
        Eddy = 'EESSSSJSSSLLLSSSJJJSSBSSLLLLSSBBSSJSSBSS'

        results = [EmptyResult, NoPairsResult, OneHelixResult, 
            ManyHelicesResult, EndsResult, FirstEndResult, LastEndResult, 
            InternalResult, Eddy]
      
        for s, r in zip(structs, results):
            c = classify(s)
            self.assertEqual(classify(s), r)

        long_struct = ".((((((((((((((.((((((..((((.....)))))))))).))..))))))))))))....(((.((((.((((((((......((((((.((..(((((((....)))).)))..))))))))...))))))))...........(((((.(..(((((((((......((((((((((((.........))))))))))))....))))).))))..)..)))))..(((((((((((((((((((......(((((((((((((((((((((((((((((((...(((.......((((((((........)))))))).......)))...))))))))))))))))))))))))))))))).((((........(((((((((((((((((((...))))))))))))))))))).......)))).....((((((((((((((((((((((((((((((.(((...))).)))))))))))))))))))))))...........))))))).))))))))))))))))))).....)))).)))......"
        
        #compare standalone method with classification from tree
        c = classify(long_struct)
        d = ViennaStructure(long_struct).toTree().classify()
        self.assertEqual(c,d)

        #Error is raised when trying to classify invalid structures
        invalid_structure = '(((..)).))))(...)(...'
        self.assertRaises(IndexError, classify, invalid_structure)
                    

class ViennaNodeTests(TestCase):
    """Tests of the ViennaNode class."""
    def setUp(self):
        """Instantiate some standard ViennaNodes."""
        self.EmptyStr = '' 
        self.NoPairsStr = '.....'
        self.OneHelixStr = '((((()))))'
        self.ManyHelicesStr = '(..(((...)).((.(((((..))).)))..((((..))))))...)'
        self.EndsStr = '..(.)..'
        self.FirstEndStr = '..((()))'
        self.LastEndStr = '((..((.))))...'
        self.InternalStr = '(((...)))..((.)).'
        #following structure is from p 25 of Eddy's WUSS description manual
        self.EddyStr = '..((((.(((...)))...((.((....))..)).)).))'
        
        #add in the tree versions by deleting trailing 'Str'
        for s in self.__dict__.keys():
            if s.endswith('Str'):
                self.__dict__[s[:-3]] = \
                    ViennaStructure(self.__dict__[s]).toTree()

    def test_str(self):
        """ViennaNode str should return Vienna-format string"""
        for s in [self.EmptyStr, self.NoPairsStr, self.OneHelixStr,
            self.ManyHelicesStr, self.EndsStr, self.InternalStr]:
            self.assertEqual(str(ViennaStructure(s).toTree()), s)

        #test with multiple-base helix in a node
        r = StructureNode()
        r.append(StructureNode())
        r.append(StructureNode(Data=Stem(1,7,5)))
        r[1].append(StructureNode())
        r.append(StructureNode())
        r.append(StructureNode())
        r.renumber()
        self.assertEqual(str(r), '.(((((.)))))..')

    def test_classify(self):
        """ViennaNode classify should return correct classification string"""
        self.assertEqual(self.Empty.classify(), '')
        self.assertEqual(self.NoPairs.classify(), 'EEEEE')
        self.assertEqual(self.OneHelix.classify(), 'SSSSSSSSSS')
        self.assertEqual(self.ManyHelices.classify(), \
            'SBBSSSLLLSSJSSBSSSSSLLSSSBSSSJJSSSSLLSSSSSSBBBS')
        self.assertEqual(self.Ends.classify(), 'EESLSEE')
        self.assertEqual(self.FirstEnd.classify(), 'EESSSSSS')
        self.assertEqual(self.LastEnd.classify(), 'SSBBSSLSSSSEEE')
        self.assertEqual(self.Internal.classify(), 'SSSLLLSSSFFSSLSSE')
        self.assertEqual(self.Eddy.classify(), \
            'EESSSSJSSSLLLSSSJJJSSBSSLLLLSSBBSSJSSBSS')

    def test_renumber(self):
        """ViennaNode renumber should assign correct numbers to nodes"""
        #should have no effect on empty structure
        se = self.Empty
        self.assertEqual(se.renumber(5), 5)
        self.assertEqual((se.Start, se.End, se.Length), (None, None, 0))
        #with no pairs, should number consecutively
        sn = self.NoPairs
        self.assertEqual(sn.renumber(5), 10)
        self.assertEqual([i.Start for i in sn], [5, 6, 7, 8, 9])
        self.assertEqual([i.End for i in sn], [None]*5)
        self.assertEqual([i.Length for i in sn], [0]*5)
        #spot checks on a complex structure
        sm = self.ManyHelices
        self.assertEqual(sm.renumber(5), 52)
        s0 = sm[0]
        self.assertEqual((s0.Start, s0.End, s0.Length), (5, 51, 1))
        s5 = sm[0][2][2][0]
        self.assertEqual(len(s5), 2)
        self.assertEqual((s5.Start, s5.End, s5.Length), (18, 33, 1))
        s6 = s5[0]
        self.assertEqual((s6.Start, s6.End, s6.Length), (19,None,0))
        #test with some helices of different lengths
        root = StructureNode()
        root.extend([StructureNode() for i in range(3)])
        root.insert(1, StructureNode(Data=Stem(3, 7, 5)))
        root.insert(3, StructureNode(Data=Stem(6,2,2)))
        root.append(StructureNode())
        self.assertEqual(root.renumber(0), 18)
        self.assertEqual(len(root), 6)
        curr = root[0]
        self.assertEqual((curr.Start,curr.End,curr.Length), (0, None, 0))
        curr = root[1]
        self.assertEqual((curr.Start, curr.End, curr.Length), (1, 10, 5))
        curr = root[2]
        self.assertEqual((curr.Start, curr.End, curr.Length), (11, None, 0))
        curr = root[3]
        self.assertEqual((curr.Start, curr.End, curr.Length), (12, 15, 2))
        curr = root[4]
        self.assertEqual((curr.Start, curr.End, curr.Length), (16, None, 0))
        curr = root[5]
        self.assertEqual((curr.Start, curr.End, curr.Length), (17, None, 0))

    def test_unpair(self):
        """StructureNode unpair should break a base pair and add correct nodes"""
        i = self.Internal
        self.assertEqual(i[0].unpair(), True)
        self.assertEqual(str(i), '.((...))...((.)).')
        e = self.Ends
        self.assertEqual(e[0].unpair(), False)
        self.assertEqual(str(e), self.EndsStr)
        o = self.OneHelix
        self.assertEqual(o[0].unpair(), True)
        self.assertEqual(str(o), '.(((()))).')
        self.assertEqual(o[1][0][0].unpair(), True)
        self.assertEqual(str(o), '.((.().)).')
        self.assertEqual(o[1].unpair(), True)
        self.assertEqual(str(o), '..(.().)..')
        self.assertEqual(o[2][1].unpair(), True)
        self.assertEqual(str(o), '..(....)..')
        self.assertEqual(o[2].unpair(), True)
        self.assertEqual(str(o), '..........')
        #test with multiple bases in helix
        r = StructureNode()
        r.append(StructureNode(Data=Stem(0,0, 5)))
        r.renumber()
        self.assertEqual(str(r), '((((()))))')
        self.assertEqual(r[0].unpair(), True)
        self.assertEqual(str(r), '.(((()))).')

    def test_pairBefore(self):
        """StructureNode pairBefore should make a pair before the current node"""
        #shouldn't be able to make any pairs if everything is paired already
        o = self.OneHelix
        for i in o:
            self.assertEqual(i.pairBefore(), False)
            
        n = self.NoPairs
        #shouldn't be able to pair at the start...
        self.assertEqual(n[0].pairBefore(), False)
        #...or at the end...
        self.assertEqual(n[-1].pairBefore(), False)
        #...but should work OK in the middle
        self.assertEqual(n[1].pairBefore(), True)
        self.assertEqual(str(n), '(.)..')

        e = self.Ends
        self.assertEqual(e[2].pairBefore(), True)
        self.assertEqual(e[1].pairBefore(), True)
        self.assertEqual(str(e), '(((.)))')
        self.assertEqual((e[0].Start, e[0].End, e[0].Length), (0,6,1))

    def test_pairAfter(self):
        """StructureNode pairAfter should create pairs after a node"""
        n = self.NoPairs
        self.assertEqual(n.pairAfter(), True)
        self.assertEqual(str(n), '(...)')
        self.assertEqual(n[0].pairAfter(), True)
        self.assertEqual(str(n), '((.))')
        self.assertEqual(n[0][0].pairAfter(), False)
        self.assertEqual(str(n), '((.))')
        curr = n[0][0]
        #check that child is correct
        self.assertEqual(len(curr), 1)
        self.assertEqual((curr[0].Start, curr[0].End, curr[0].Length), \
            (2,None,0))
        #check that pair is correct
        self.assertEqual((curr.Start, curr.End, curr.Length), (1,3,1))

        m = self.ManyHelices
        n = m[0][2][0][0]
        self.assertEqual(n.pairAfter(), True)
        self.assertEqual(str(m), \
            '(..((((.))).((.(((((..))).)))..((((..))))))...)')
        self.assertEqual(n[0].pairAfter(), False)
    
    def test_pairChildren(self):
        """StructureNode PairChildren should make the correct pairs"""
        n = ViennaStructure('.....').toTree()   #same as self.NoPairs
        self.assertEqual(n.pairChildren(0, 4), True)
        self.assertEqual(str(n), '(...)')
        n = ViennaStructure('.....').toTree()   #same as self.NoPairs
        self.assertEqual(n.pairChildren(1, 4), True)
        self.assertEqual(str(n), '.(..)')
        n = ViennaStructure('.....').toTree()   #same as self.NoPairs
        #can't pair same object
        self.assertEqual(n.pairChildren(1, 1), False)
        self.assertEqual(str(n), '.....')
        self.assertEqual(n.pairChildren(1, -1), True)
        self.assertEqual(str(n), '.(..)')
        #can't pair something already paired
        self.assertEqual(n.pairChildren(0,1), False)
        #IndexError if out of range
        self.assertRaises(IndexError, n.pairChildren, 0, 5)
        n.append(StructureNode())
        n.append(StructureNode())
        n.renumber()
        self.assertEqual(str(n), '.(..)..')
        self.assertEqual(n.pairChildren(0, -2), True)
        self.assertEqual(str(n), '((..)).')

    def test_expand(self):
        """StructureNode expand should extend helices."""
        s = StructureNode(Data=(Stem(1, 10, 3)))
        s.append(StructureNode())
        #need to make a root node for consistency
        r = StructureNode()
        r.append(s)
        self.assertEqual(str(s), '(((.)))')
        s.expand()
        self.assertEqual(str(s), '(((.)))')
        self.assertEqual((s.Start, s.End, s.Length), (1, 10, 1))
        n = s[0]
        self.assertEqual((n.Start, n.End, n.Length), (2, 9, 1))
        n = s[0][0]
        self.assertEqual((n.Start, n.End, n.Length), (3, 8, 1))
        n = s[0][0][0]
        self.assertEqual((n.Start, n.End, n.Length), (None, None, 0))
        s.renumber()
        self.assertEqual((s.Start, s.End, s.Length), (0, 6, 1))
        n = s[0]
        self.assertEqual((n.Start, n.End, n.Length), (1, 5, 1))
        n = s[0][0]
        self.assertEqual((n.Start, n.End, n.Length), (2, 4, 1))
        n = s[0][0][0]
        self.assertEqual((n.Start, n.End, n.Length), (3, None, 0))
        #check that it's not recursive
        s[0][0].append(StructureNode(Data=Stem(20, 24, 2)))
        s.expand()
        n = s[0][0][-1]
        self.assertEqual((n.Start, n.End, n.Length), (20, 24, 2))
        n.expand()
        self.assertEqual((n.Start, n.End, n.Length), (20, 24, 1))
        n = n[0]
        self.assertEqual((n.Start, n.End, n.Length), (21, 23, 1))
        
    def test_expandAll(self):
        """StructureNode expandAll should act recursively"""
        r = StructureNode()
        r.append(StructureNode(Data=Stem(0, 6, 4)))
        r.append(StructureNode(Data=Stem(0, 6, 3)))
        r.append(StructureNode())
        r[0].append(StructureNode())
        r[0].append(StructureNode(Data=Stem(0,6,2)))
        r[0][-1].append(StructureNode())
        r.renumber()
        self.assertEqual(str(r), '((((.((.))))))((())).')
        r.expandAll()
        self.assertEqual(str(r), '((((.((.))))))((())).')
        expected_nodes = [
            (None, None, 0),
            (0, 13, 1),
            (1, 12, 1),
            (2, 11, 1),
            (3, 10, 1),
            (4, None, 0),
            (5, 9, 1),
            (6, 8, 1),
            (7, None, 0),
            (14, 19, 1),
            (15, 18, 1),
            (16, 17, 1),
            (20, None, 0),
        ]
        for obs, exp in zip(r.traverse(), expected_nodes):
            self.assertEqual((obs.Start, obs.End, obs.Length), exp)

    def test_collapse(self):
        """StructureNode collapse should collapse consecutive pairs from self"""
        one = ViennaStructure('(.)').toTree()
        self.assertEqual(one.collapse(), False)
        self.assertEqual(str(one), '(.)')
        two = ViennaStructure('((.))').toTree()
        #can't collapse root node
        self.assertEqual(two.collapse(), False)
        #should be able to collapse next node
        self.assertEqual(two[0].collapse(), True)
        self.assertEqual((two[0].Start, two[0].End, two[0].Length), (0,4,2))
        self.assertEqual(str(two), '((.))')
        three = ViennaStructure('(((...)))..').toTree()
        self.assertEqual(three[0].collapse(), True)
        self.assertEqual((three[0].Start, three[0].End, three[0].Length), \
            (0,8,3))
        self.assertEqual(str(three), '(((...)))..')
        self.assertEqual(three[0].collapse(), False)
        self.assertEqual(three[-1].collapse(), False)

        oh = self.OneHelix
        self.assertEqual(oh[0].collapse(), True)
        self.assertEqual(str(oh), '((((()))))')
        
    def test_collapseAll(self):
        """StructureNode collapseAll should collapse consecutive pairs"""
        for s in [self.Empty, self.NoPairs, self.OneHelix, self.ManyHelices,\
            self.Ends, self.FirstEnd, self.LastEnd, self.Internal, self.Eddy]:
            before = str(s)
            s.collapseAll()
            after = str(s)
            self.assertEqual(after, before)

        oh = self.OneHelix[0]
        self.assertEqual((oh.Start, oh.End, oh.Length), (0,9,5))
        m_obs = self.ManyHelices.traverse()
        m_exp = [
            (None, None, 0),
            (0, 46, 1),
            (1, None, 0),
            (2, None, 0),
            (3, 42, 1),
            (4, 10, 2),
            (6, None, 0),
            (7, None, 0),
            (8, None, 0),
            (11, None, 0),
            (12, 41, 1),
            (13, 28, 1),
            (14, None, 0),
            (15, 27, 2),
            (17, 24, 3),
            (20, None, 0),
            (21, None, 0),
            (25, None, 0),
            (29, None, 0),
            (30, None, 0),
            (31, 40, 4),
            (35, None, 0),
            (36, None, 0),
            (43, None, 0),
            (44, None, 0),
            (45, None, 0),
            (46, None, 0),
        ]
        for obs, exp in zip([(i.Start, i.End, i.Length) for i in m_obs], m_exp):
            self.assertEqual(obs, exp)

    def test_breakBadPairs(self):
        """StructureNode breakBadPairs should eliminate mispaired bases."""
        oh_str = ViennaStructure(self.OneHelixStr)
        #no change if all pairs valid
        oh = oh_str.toTree()
        oh.breakBadPairs(Rna('CCCCCGGGGG'))
        self.assertEqual(str(oh), str(oh_str))
        #break everything if all pairs invalid
        oh.breakBadPairs(Rna('CCCCCAAAAA'))
        self.assertEqual(str(oh), '..........')
        #break a single pair
        oh = oh_str.toTree()
        oh.breakBadPairs(Rna('GCCCCGGGGG'))
        self.assertEqual(str(oh), '.(((()))).')
        #break two pairs
        oh = oh_str.toTree()
        oh.breakBadPairs(Rna('GCCCCCGGGG'))
        self.assertEqual(str(oh), '.(((..))).')
        #break internal pairs
        oh = oh_str.toTree()
        oh.breakBadPairs(Rna('GCCGCGGGGG'))
        self.assertEqual(str(oh), '.((.().)).')
        #repeat with multiple independent helices
        th_str = ViennaStructure('((.)).((.))')
        th = th_str.toTree()
        th.breakBadPairs(Rna('CCUGGCUUCGG'))
        self.assertEqual(str(th), th_str)
        th.breakBadPairs(Rna('CGUAGCAGUUU'))
        self.assertEqual(str(th), '(...).((.))')
        th = th_str.toTree()
        th.breakBadPairs(Rna('UUUUUUUUUUU'))
        self.assertEqual(str(th), '...........')

    def test_extendHelix(self):
        """StructureNode extendHelix should extend the helix as far as possible
        """
        #single paired node is root[4]
        op_str = ViennaStructure('....(......)...')
        op = op_str.toTree()
        #can't extend if base pairs not allowed
        op[4].extendHelix(Rna('AAAAAAAAAAAAAAA'))
        self.assertEqual(str(op), op_str)
        #should extend a pair 5'
        op[4].extendHelix(Rna('AAACCAAAAAAGGAA'))
        self.assertEqual(str(op), '...((......))..')
        #should extend multiple pairs 5'
        op = op_str.toTree()
        op[4].extendHelix(Rna('CCCCCUUUUUUGGGG'))
        self.assertEqual(str(op), '.((((......))))')
        #should extend a pair 3', but must leave > 2-base loop
        op = op_str.toTree()
        op[4].extendHelix(Rna('AAAACCCCGGGGAAA'))
        self.assertEqual(str(op), '....((....))...')
        op[4][0].insert(1, StructureNode(Data=Stem(Start=1,End=1,Length=1)))
        op.renumber()
        self.assertEqual(str(op), '....((.()...))...')
        op[4][0].extendHelix(Rna( 'AAAACCCUACGGGGAAA'))
        self.assertEqual(str(op), '....(((()..)))...')
        #should extend a pair in both directions if possible
        op = op_str.toTree()
        op[4].extendHelix(Rna('AAACCCAAAAGGGAA'))
        self.assertEqual(str(op), '...(((....)))..')
        
    def test_extendHelices(self):
        """StructureNode extendHelices should extend all helices"""
        e = ViennaStructure('........')
        t = e.toTree()
        t.extendHelices(Rna('CCCCCCCCCC'))
        self.assertEqual(str(t), e)
        #no pairs if sequence can't form them
        s = ViennaStructure('(.....(...)..)...((.....))...')
        r =             Rna('AAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
        t = s.toTree()
        t.extendHelices(r)
        self.assertEqual(str(t), s)
        #should be able to extend a single helix
        s = ViennaStructure('(.....(...)..)...((.....))...')
        r =             Rna('CAAAAACAAAGAAGCCCCCCCAGGGGGGG')
        t = s.toTree()
        t.extendHelices(r)
        self.assertEqual(str(t), '(.....(...)..)((((((...))))))')
        #should be able to extend multiple helices
        s = ViennaStructure('(.....(...)..)...((.....))...')
        r =             Rna('AAAAACCCAGGGUUCCCCCAUAAAGGGAA')
        t = s.toTree()
        t.extendHelices(r)
        self.assertEqual(str(t), '((...((...))))..(((.....)))..')
        
    def test_fitSeq(self):
        """StructureNode fitSeq should adjust structure to match sequence"""
        #this is just a minimal test, since we know that both breakBadPairs()
        #and extendHelices() work fine with more extensive tests.
        s = ViennaStructure('..(((.....)))......(((.....)))...')
        r = Rna(            'UCCCCACUGAGGGGUUUGGGGGGUUUUCGCCCU')
        t = s.toTree()
        t.fitSeq(r)
        self.assertEqual(str(t), '.((((.....))))...(((.((...)).))).')

#run the test suites if invoked as a script from the command line
if __name__ == "__main__":
    main()

