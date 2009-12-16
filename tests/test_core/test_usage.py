#!/usr/bin/env python
"""Tests of classes dealing with base, codon, amino acid usage.
"""
from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import FunctionWrapper
from cogent.core.usage import InfoFreqs, AminoAcidUsage, BaseUsage, CodonUsage,\
    PositionalBaseUsage, UnsafeBaseUsage, EqualBases, DinucUsage
from cogent.core.genetic_code import GeneticCodes, GeneticCode

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class InfoFreqsTests(TestCase):
    """Tests of the InfoFreqs base class"""
    class kp_empty(InfoFreqs):
        pass
        
    class kp(InfoFreqs):
        Mask = FunctionWrapper(int)
        RequiredKeys = dict.fromkeys([1,2,3])

    class has_info(list):
        def __new__(cls, data, *args):
            return list.__new__(cls, data)
        def __init__(self, data, info):
            list.__init__(self, data)
            self.Info = info

    def test_init_empty(self):
        """InfoFreqs empty init should include only required items"""
        self.assertEqual(self.kp_empty(), {})
        self.assertEqual(self.kp(), {1:0, 2:0, 3:0})

    def test_init_data(self):
        """InfoFreqs init with data should count freqs"""
        self.assertEqual(self.kp_empty('1qazaz'), {'1':1,'q':1,'a':2,'z':2})
        self.assertEqual(self.kp([1,2,3,1]), {1:2, 2:1, 3:1})
        #type of exception on invalid conversion currently not specified
        self.assertRaises(Exception, self.kp, '123q')

    def test_init_info(self):
        """InfoFreqs init from object with Info should get ref to Info"""
        s = self.has_info([1,1,1,2,2,3], {'test':'data'})
        obj = self.kp(s)
        self.assertEqual(obj, {1:3, 2:2, 3:1})
        self.assertEqual(obj.Info, {'test':'data'})

        #should use its own info if supplied, though
        obj = self.kp(s, {'new':'info'})
        self.assertEqual(obj, {1:3, 2:2, 3:1})
        self.assertNotEqual(obj.Info, {'test':'data'})
        self.assertEqual(obj.Info, {'new':'info'})

    def test_info(self):
        """InfoFreqs info changes should work as expected."""
        class has_attr(object):
            def __init__(self, attr):
                self.Attr = attr
        s = self.has_info([1,1,1,2,2,3,3], {'test':'data'})
        obj = self.kp(s)
        self.assertEqual(obj, {1:3, 2:2, 3:2})
        self.assertEqual(obj.Info, {'test':'data'})
        self.assertRaises(AttributeError, obj.__getattr__, 'Attr')
        obj.Info = has_attr('Test')
        self.assertEqual(obj.Attr, 'Test')

    def test_normalize(self):
        """InfoFreqs normalize should delete non-required keys, if any required"""
        e = self.kp_empty('1qaa')
        self.assertEqual(e, {'1':1,'q':1,'a':2,})
        e.normalize()
        #should not delete any keys, since kp_empty has no required keys
        self.assertEqual(e, {'1':0.25, 'q':0.25, 'a':0.5})
        k = self.kp([1,2,3,1,4,4,4])
        self.assertEqual(k, {1:2, 2:1, 3:1, 4:3})
        k.normalize()
        #should delete 4, since it's not required
        self.assertEqual(k, {1:0.5, 2:0.25, 3:0.25})

class BaseUsageTests(TestCase):
    """Tests of the BaseUsage class."""
    def test_init_empty(self):
        """BaseUsage should init with empty freqs"""
        b = BaseUsage()
        self.assertEqual(len(b), 4)
        for nt in 'UTACGutacg':
            assert nt in b
            self.assertEqual(b[nt], 0)
        for nt in 'UCAG':
            assert nt in b.keys()

        items = list(iter(b))
        items.sort()
        self.assertEqual(items, ['A', 'C', 'G', 'U'])

    def test_init_data(self):
        """BaseUsage should init with arbitrary data"""
        b = BaseUsage('UUUUUGGGCA')
        self.assertEqual(b, {'U':5, 'G':3, 'C':1, 'A':1})
        b.normalize()
        self.assertEqual(b, {'U':0.5, 'G':0.3, 'C':0.1, 'A':0.1})

    def test_setitem(self):
        """BaseUsage should map keys on setitem"""
        b = BaseUsage()
        b['t'] = 3
        b['G'] = 3
        b.normalize()
        i = b.items()
        i.sort()
        self.assertEqual(i, [('A',0.0),('C',0.0),('G',0.5),('U',0.5)])

    def test_getitem(self):
        """BaseUsage should map key on getitem"""
        b = BaseUsage({'a':3, 'T':2, 'X':1})
        self.assertEqual(b['X'], 1)
        self.assertEqual(b['A'], 3)
        self.assertEqual(b['U'], 2)
        self.assertEqual(b['a'], 3)
        self.assertEqual(b['t'], 2)
        b.normalize()
        assert 'X' not in b
        self.assertFloatEqual(b['A'], 0.6)
        self.assertFloatEqual(b['u'], 0.4)

    def test_getitem_multi(self):
        """BaseUsage should get total frequency on getitem with 2-letter key"""
        b = BaseUsage({'a':3, 'T':2, 'C':5, 'X':1})
        self.assertEqual(b['AT'], 0.5)
        self.assertEqual(b['AC'], 0.8)
        self.assertEqual(b['TC'], 0.7)
        self.assertEqual(b['AX'], 0.3)
    
    def test_copy(self):
        """BaseUsage copy should work correctly"""
        b = BaseUsage({'a':3, 'T':2, 'C':5, 'X':1})
        c = b.copy()
        self.assertEqual(c['AT'], 0.5)

    def test_bases(self):
        """BaseUsage bases() should return same object"""
        b = BaseUsage({'a':3, 'T':2, 'C':5, 'X':1})
        c = b.bases()
        assert b is c

    def test_codons(self):
        """BaseUsage codons should return most likely codon freqs"""
        b = BaseUsage({'a':3, 'T':2, 'X':1})
        c = b.codons()
        known = {
            'AAA' : .6 * .6 * .6,
            'AAU' : .6 * .6 * .4,
            'AUA' : .6 * .4 * .6,
            'AUU' : .6 * .4 * .4,
            'UAA' : .4 * .6 * .6,
            'UAU' : .4 * .6 * .4,
            'UUA' : .4 * .4 * .6,
            'UUU' : .4 * .4 * .4,
        }
        for codon in c:
            if codon in known:
                self.assertFloatEqual(c[codon], known[codon])
            else:
                self.assertEqual(c[codon], 0)

    def test_positionalBases(self):
        """BaseUsage positionalBases should have copy of self at each position"""
        b = BaseUsage('A')
        p = b.positionalBases()
        for i in p:
            assert i is not b
            self.assertEqual(i, {'A':1,'U':0,'C':0,'G':0})

    def test_aminoAcids(self):
        """BaseUsage aminoAcids should give the same results as the codons"""
        known_data = {
            'AAA' : .6 * .6 * .6,
            'AAU' : .6 * .6 * .4,
            'AUA' : .6 * .4 * .6,
            'AUU' : .6 * .4 * .4,
            'UAA' : .4 * .6 * .6,
            'UAU' : .4 * .6 * .4,
            'UUA' : .4 * .4 * .6,
            'UUU' : .4 * .4 * .4,
        }
        known = CodonUsage(known_data)
        b = BaseUsage({'a':3, 'T':2, 'X':1})
        self.assertEqual(b.aminoAcids(), known.aminoAcids())
        #check that the genetic code is passed through correctly
        all_g = GeneticCode('G'*64)
        self.assertEqual(b.aminoAcids(all_g), AminoAcidUsage({'G':1}))
        
    def test_normalize(self):
        """BaseUsage normalize should work when empty"""
        b = BaseUsage()
        b.normalize()
        self.assertEqual(b, {'U':0,'C':0,'A':0,'G':0})
        b = BaseUsage('AACG')
        b.normalize()
        self.assertEqual(b, {'U':0, 'C':0.25, 'A':0.5, 'G':0.25})

    def test_distance(self):
        """BaseUsage distance() should return dist between two BUs"""
        #absolute numbers, will normalize to calculate distance
        self.assertFloatEqual(BaseUsage('GC').distance(BaseUsage('AU')),1)
        self.assertFloatEqual(BaseUsage('AU').distance(BaseUsage('GC')),1)
        self.assertFloatEqual(BaseUsage('GCAU').distance(BaseUsage('AUAU')),.5)
        #should work even against dict with 'T's
        self.assertFloatEqual(BaseUsage('GCAU').distance(\
            BaseUsage({'A':2,'T':2,'C':0,'G':0,})),.5)
        #rounding error
        self.assertEqual(BaseUsage('ACG').distance(BaseUsage('CCGGAA')),0)
        #normalized - as in unit simplex
        ag = BaseUsage('AG')
        ag.normalize()
        uc = BaseUsage('UC')
        uc.normalize()
        self.assertFloatEqual(ag.distance(uc),1)
        self.assertFloatEqual(BaseUsage({'A':0.4,'G':0.1,'C':0.4,'U':0.1})\
        .distance(BaseUsage({'A':0.1,'G':0.4,'U':0.4,'C':0.1})),0.6)
        self.assertFloatEqual(BaseUsage({'A':0.25,'G':0.25,'C':0.25,'U':0.25})\
            .distance(BaseUsage({'A':0.25,'G':0.25,'U':0.25,'C':0.25})),0.0)
        self.assertFloatEqual(BaseUsage({'A':0.245,'G':0.255,'C':0.245,\
            'U':0.255}).distance(BaseUsage({'A':0.255,'G':0.245,'U':0.245,\
            'C':0.255})),0.02)
        self.assertFloatEqual(BaseUsage({'A':0.245,'G':0.255,'C':0.245,\
            'U':0.255}).distance(BaseUsage({'A':0.25,'G':0.25,'U':0.25,\
            'C':0.25})),0.01)
        self.assertFloatEqual(BaseUsage({'A':0.248,'G':0.252,'C':0.248,\
            'U':0.252}).distance(BaseUsage({'A':0.25,'G':0.25,'U':0.25,\
            'C':0.25})),0.004)

    def test_content(self):
        """BaseUsage content should return sum of specified bases."""
        b = BaseUsage('UUUUUCCCAG')
        #should work for combinations
        self.assertEqual(b.content('UCAG'), 10)
        self.assertEqual(b.content('GC'), 4)
        self.assertEqual(b.content('AU'), 6)
        #should map T to U
        self.assertEqual(b.content('AT'), 6)
        #should work for single bases
        self.assertEqual(b.content('U'), 5)
        #shouldn't complain about invalid bases
        self.assertEqual(b.content('X'), 0)

    def test_add(self):
        """BaseUsage add should sum two base usages"""
        b = BaseUsage('U')
        b2 = BaseUsage('C')
        self.assertEqual(b+b2, BaseUsage('UC'))
        b += b2
        self.assertEqual(b, BaseUsage('UC'))

    def test_toCartesian(self):
        """BaseUsage toCartesian should return x, y, z from instance"""
        b = BaseUsage('ACGU')
        self.assertEqual(b.toCartesian(), (0.5,0.5,0.5))
        b = BaseUsage('A')
        self.assertEqual(b.toCartesian(), (0,0,1))
        b = BaseUsage('CGA')
        self.assertEqual(b.toCartesian(), (1/3.0,1/3.0,1/3.0))

    def test_fromCartesian(self):
        """BaseUsage fromCartesian should init instance from x,y,z"""
        b = BaseUsage.fromCartesian(0.5,.5,.5)
        self.assertFloatEqual(b['A'], 0.25)
        self.assertFloatEqual(b['C'], 0.25)
        self.assertFloatEqual(b['G'], 0.25)
        self.assertFloatEqual(b['U'], 0.25)
        b = BaseUsage.fromCartesian(1/3.0, 1/3.0, 1/3.0)
        self.assertEqual(b['U'], 0)
        self.assertEqual(b['A'], 1/3.0)

class UnsafeBaseUsageTests(TestCase):
    """Tests of the UnsafeBaseUsage class."""

    def test_init(self):
        """UnsafeFreqs uses dict init, so must use += after creation for others"""
        self.assertRaises(ValueError, UnsafeBaseUsage, 'acguc')
        u = UnsafeBaseUsage()
        u += 'acgtc'
        #note lack of conversion to uppercase RNA!
        self.assertEqual(u, {'a':1,'c':2,'g':1,'t':1})

    def test_normalize(self):
        """UnsafeFreqs should normalize based on the uppercase RNA alphabet"""
        u = UnsafeBaseUsage()
        #note that the lower-case u will be discarded, not converted
        u += 'AAACCGGGGGu'
        u.normalize()
        self.assertFloatEqual(u, {'A':0.3, 'C':0.2, 'G':0.5, 'U':0.0})

class CodonUsageTests(TestCase):
    """Tests of the CodonUsage class."""
    
    def test_init_empty(self):
        """Empty CodonUsage init should have 64 codons, all 0"""
        u = CodonUsage()
        self.assertEqual(len(u), 64)
        for i in u:
            self.assertEqual(u[i], 0)
        #check that the genetic code is the default
        assert u.GeneticCode is GeneticCodes[1]

    def test_init_string(self):
        """CodonUsage should count codons in string"""
        u = CodonUsage('UUUCCCUUUUUUGA')
        self.assertEqual(u, CodonUsage({'UUU':3, 'CCC':1, 'GA':1}))
        u.normalize()
        self.assertEqual(u, CodonUsage({'UUU':0.75, 'CCC':0.25}))

    def test_getitem(self):
        """CodonUsage should allow lookup as RNA or DNA, case-insensitive"""
        u = CodonUsage()
        rna, dna, lc = 'UCAG', 'TCAG', 'ucag'
        for a in [rna, dna, lc]:
            codons = [i+j+k for i in a for j in a for k in a]
            for c in codons:
                self.assertEqual(u[c], 0)

    def test_bases(self):
        """CodonUsage bases should count bases correctly"""
        u = CodonUsage('UUUCCCUAGCCCGGGAA')
        b = u.bases()
        self.assertEqual(b, BaseUsage('UUUCCCUAGCCCGGGAA'))
        #purge_unwanted should get rid of bad codons
        b = u.bases(purge_unwanted=True)
        self.assertEqual(b, BaseUsage('UUUCCCCCCGGG'))

    def test_codons(self):
        """CodonUsage codons should return same object"""
        u = CodonUsage('abc')
        c = u.codons()
        assert u is c

    def test_aminoAcids(self):
        """CodonUsage aminoAcids should correctly count amino acids"""
        freqs = {'UUC':5, 'AUA':10, 'AUG':10, 'CGC':3, 'AGG':2, 'XYZ':8,
            'UAA':2, 'UGA':1}
        u = CodonUsage(freqs, "test")
        self.assertEqual(u.Info, 'test')
        for key, val in u.items():
            if key in freqs:
                self.assertEqual(val, freqs[key])
            else:
                self.assertEqual(val, 0)
        aa = u.aminoAcids()
        self.assertEqual(aa,
            AminoAcidUsage({'F':5,'I':10,'M':10,'R':5,'*':3,'X':8}))
        #check that it works with a different genetic code
        u.GeneticCode = GeneticCodes['2']
        aa = u.aminoAcids()
        self.assertEqual(aa,
            AminoAcidUsage({'F':5,'I':0,'M':20,'R':3,'*':4,'W':1,'X':8}))
        #check that it works if a genetic code is supplied explicitly
        u.GeneticCode = GeneticCodes[1]
        aa = u.aminoAcids()
        self.assertEqual(aa,
            AminoAcidUsage({'F':5,'I':10,'M':10,'R':5,'*':3,'X':8}))
        aa_2 = u.aminoAcids(2)
        self.assertEqual(aa_2,
            AminoAcidUsage({'F':5,'I':0,'M':20,'R':3,'*':4,'W':1,'X':8}))
        #check that we held onto the info object through the above
        self.assertEqual(aa_2.Info, 'test')

    def test_positionalBases(self):
        """CodonUsage bases should count bases at each position correctly"""
        freqs = {'UUC':5, 'AUA':10, 'AUG':10, 'CGC':3, 'AGG':2, 'XYZ':8,
            'UAA':2, 'UGA':1}
        u = CodonUsage(freqs)
        b = u.positionalBases()
        assert isinstance(b, PositionalBaseUsage)
        first, second, third = b
        self.assertEqual(first, BaseUsage({'U':8,'C':3,'A':22,'X':8}))
        self.assertEqual(second, BaseUsage({'U':25,'C':0,'A':2,'G':6,'Y':8}))
        self.assertEqual(third, BaseUsage({'C':8,'A':13,'G':12,'Z':8}))
        #check that it also works when we purge
        p = u.positionalBases(purge_unwanted=True)
        first, second, third = p
        self.assertEqual(first, BaseUsage({'U':5,'C':3,'A':2}))
        self.assertEqual(second, BaseUsage({'U':5,'G':5}))
        self.assertEqual(third, BaseUsage({'C':8,'G':2}))
        #check that it also works with a different genetic code, and, 
        #incidentally, that the purging didn't affect the original object
        u.GeneticCode = GeneticCodes[2] #mt code: different stop codons
        p = u.positionalBases(purge_unwanted=True)
        first, second, third = p
        self.assertEqual(first, BaseUsage({'U':6,'C':3,'A':20}))
        self.assertEqual(second, BaseUsage({'U':25,'G':4}))
        self.assertEqual(third, BaseUsage({'C':8,'A':11,'G':10}))

    def test_positionalGC(self):
        """CodonUsage positionalGC should give correct GC contents."""
        c = EqualBases.codons()
        self.assertEqual(c.positionalGC(False), [0.5,0.5,0.5,0.5])
        c = EqualBases.codons()
        self.assertNotEqual(c.positionalGC(True), [0.5,0.5,0.5,0.5])

    def test_fingerprint(self):
        """CodonUsage fingerprint should give correct ratios."""
        c = EqualBases.codons()
        f = c.fingerprint()
        self.assertEqual(len(f), 9)
        self.assertEqual(f, \
            [[.5,.5,.125] for i in range(8)] + [[.5,.5,1]])
        #should be able to omit mean...
        f = c.fingerprint(include_mean=False)
        self.assertEqual(f, [[.5,.5,.125] for i in range(8)])
        #...or use all doublets
        f = c.fingerprint(include_mean=False, which_blocks='all')
        self.assertEqual(len(f), 16)
        #...or do just the non-quartet ones
        f = c.fingerprint(include_mean=False, which_blocks='split')
        self.assertEqual(len(f), 6)
        #check that it doesn't fail on an empty codon usage
        c = CodonUsage('')
        f = c.fingerprint()
        self.assertEqual(f[0], [0.5, 0.5, 0])

    def test_pr2bias(self):
        """CodonUsage pr2bias should give correct ratios."""
        c = EqualBases.codons()
        b = c.pr2bias('UU')
        self.assertEqual(len(b), 6)
        self.assertEqual(b, tuple([.5]*6))
        c = CodonUsage()
        c['ACU'] = 10
        c['ACC'] = 5
        c['ACA'] = 15
        c['ACG'] = 20
        self.assertEqual(c.pr2bias('AC'), (20/25,15/25,20/35,5/15,20/30,5/20))

    def test_add(self):
        """CodonUsage add should sum two base usages"""
        c = CodonUsage('UUU')
        c2 = CodonUsage('CCC')
        self.assertEqual(c+c2, CodonUsage('UUUCCC'))
        c += c2
        self.assertEqual(c, CodonUsage('UUUCCC'))

    def test_rscu(self):
        """CodonUsage rscu should calculate synonymous usage correctly"""
        c = CodonUsage({'UUU':3,'UUC':1,'ACA':1})
        c.rscu()
        self.assertEqual(c['UUU'], 0.75)
        self.assertEqual(c['UUC'], 0.25)
        self.assertEqual(c['ACA'], 1)
        self.assertEqual(c['GGG'], 0)


class PositionalBaseUsageTests(TestCase):
    """Tests of the PositionalBaseUsage class."""
    def test_init_empty(self):
        """PositionalBaseUsage init when empty should set all freqs to 0"""
        p = PositionalBaseUsage()
        assert p[0] is not p[1]
        assert p[1] is not p[2]
        assert p[0] is not p[2]
        for i in p:
            self.assertEqual(i, BaseUsage())

    def test_info(self):
        """PositionalBaseUsage info should work as expected"""
        #test a cycle of setting the Info and setting it back again
        p = PositionalBaseUsage()
        self.assertRaises(AttributeError, getattr, p, 'upper')
        p.Info = 'xyz'
        self.assertEqual(p.upper(), 'XYZ')
        p.Info = None
        self.assertRaises(AttributeError, getattr, p, 'upper')

    def test_getitem(self):
        """PositionalBaseUsage getitem should return 1st, 2nd, 3rd in order"""
        a, c, g = BaseUsage('A'), BaseUsage('C'), BaseUsage('G')
        p = PositionalBaseUsage(a, c, g)
        assert p.First is a
        assert p.Second is c
        assert p.Third is g
        #make sure they're not all the same object
        assert p.First is not g
        #test positive indices
        assert p[0] is p.First
        assert p[1] is p.Second
        assert p[2] is p.Third
        #test negative indices
        assert p[-1] is p.Third
        assert p[-2] is p.Second
        assert p[-3] is p.First
        try:
            x = p[3]
        except IndexError:
            pass
        else:
            self.fail("Failed to raise IndexError on bad index")
        #test iteration
        for o, e in zip(p, [a, c, g]):
            assert o is e

    def test_normalize(self):
        """PositionalBaseUsage normalize should normalize each position"""
        a, c, g = BaseUsage('AAGC'), BaseUsage('CCGA'), BaseUsage('GGCA')
        p = PositionalBaseUsage(a, c, g)
        self.assertEqual(p[0], {'A':2, 'C':1, 'G':1, 'U':0})
        p.normalize()
        self.assertEqual(p[0], {'A':0.5, 'C':0.25, 'G':0.25, 'U':0})
        self.assertEqual(p[1], {'A':0.25, 'C':0.5, 'G':0.25, 'U':0})
        self.assertEqual(p[2], {'A':0.25, 'C':0.25, 'G':0.5, 'U':0})
    
    def test_bases(self):
        """PositionalBaseUsage bases should sum bases at each position"""
        a, c, g = BaseUsage('AAGC'), BaseUsage('CCGA'), BaseUsage('GGCA')
        p = PositionalBaseUsage(a, c, g)
        b = p.bases()
        self.assertEqual(b, BaseUsage('AAGCCCGAGGCA'))

    def test_codons(self):
        """PositionalBaseUsage codons should give expected codon freqs"""
        #one of each base should give freqs if 1/64 for everything
        orig = CodonUsage('UUUCCCAAAGGG')
        b = orig.positionalBases()
        final = b.codons()
        self.assertEqual(len(final), 64)
        for i in final:
            self.assertFloatEqual(final[i], 1.0/64)

        #two bases at each position should give correct freqs
        orig = CodonUsage('UCGAGUUCGUCG')
        final = orig.positionalBases().codons()
        exp = {
            'UCG':  0.75 * 0.75 * 0.75,
            'UCU':  0.75 * 0.75 * 0.25,
            'UGG':  0.75 * 0.25 * 0.75,
            'UGU':  0.75 * 0.25 * 0.25,
            'ACG':  0.25 * 0.75 * 0.75,
            'ACU':  0.25 * 0.75 * 0.25,
            'AGG':  0.25 * 0.25 * 0.75,
            'AGU':  0.25 * 0.25 * 0.25,
            }
        for f in final:
            if f in exp:
                self.assertFloatEqual(final[f], exp[f])
            else:
                self.assertEqual(final[f], 0)

    def test_positionalBases(self):
        """PositionalBaseUsage positionalBases should return same object"""
        p = PositionalBaseUsage()
        x = p.positionalBases()
        assert p is x

    def test_aminoAcids(self):
        """PositionalBaseUsage aminoAcids should return correct amino acids"""
        #check hand-calculated values on a particular sequence
        orig = CodonUsage('UCGAGUUCGUCG')
        final = orig.positionalBases().aminoAcids()
        exp = {
            'S':  0.75 * 0.75 * 0.75 + 0.75 * 0.75 * 0.25 + 0.25*0.25*0.25,
            'W':  0.75 * 0.25 * 0.75,
            'C':  0.75 * 0.25 * 0.25,
            'T':  0.25 * 0.75 * 0.75 + 0.25 * 0.75 * 0.25,
            'R':  0.25 * 0.25 * 0.75,
            }
        for f in final:
            if f in exp:
                self.assertFloatEqual(final[f], exp[f])
            else:
                self.assertEqual(final[f], 0)

        #test for unbiased freqs on a couple of different genetic codes
        orig = CodonUsage('UUUCCCAAAGGG')
        final = orig.positionalBases().aminoAcids()
        SGC = GeneticCodes[1]
        for aa in final:
            self.assertEqual(final[aa], len(SGC[aa])/64.0)
        mt = GeneticCodes[2]
        final_mt = orig.positionalBases().aminoAcids(mt)
        self.assertNotEqual(final, final_mt)
        for aa in final_mt:
            self.assertEqual(final_mt[aa], len(mt[aa])/64.0)
        

class AminoAcidUsageTests(TestCase):
    """Tests of the AminoAcidUsage class."""
    def test_init_empty(self):
        """AminoAcidUsage should init with empty freqs"""
        a = AminoAcidUsage()
        for key, val in a.items():
            self.assertEqual(val, 0)
        self.assertEqual(len(a), 21)
        assert 'A' in a
        assert 'a' in a

    def test_init_data(self):
        """AminoAcidUsage should init with data"""
        a = AminoAcidUsage('aadddx')
        self.assertEqual(a['A'], 2)
        self.assertEqual(a['d'], 3)
        self.assertEqual(a['X'], 1)
        a.normalize()
        self.assertEqual(a['a'], 0.4)
        self.assertEqual(a['d'], 0.6)
        assert 'x' not in a

    def test_bases(self):
        """AminoAcidUsage bases should return most likely base freqs"""
        a = AminoAcidUsage('GGG')
        self.assertEqual(a.bases(),
            {'G':9.0/12,'U':1.0/12,'C':1.0/12,'A':1.0/12})

        a = AminoAcidUsage('CAGTWERQWE')
        exp = a.codons().bases()
        exp.normalize()
        self.assertFloatEqual(a.bases(), exp)
        
    def test_codons(self):
        """AminoAcidUsage codons should return most likely codon freqs"""
        a = AminoAcidUsage('GGG')
        c = CodonUsage('GGUGGCGGAGGG')
        c.normalize()
        self.assertEqual(a.codons(), c)
        a = AminoAcidUsage('D')
        c = CodonUsage('GAUGAC')
        c.normalize()
        self.assertEqual(a.codons(), c)

        a = AminoAcidUsage('GDDFMM')
        c = CodonUsage('GGUGGCGGAGGG'+'GAUGAC'*4+'UUUUUC'*2+'AUG'*8)
        c.normalize()
        self.assertEqual(a.codons(), c)

        a = AminoAcidUsage('II*')
        c = CodonUsage('AUUAUCAUA'*2+'UAAUAGUGA')
        c.normalize()
        self.assertEqual(a.codons(), c)

        #check that it works with a nonstandard code
        code = GeneticCode('A'*4+'C'*28+'G'*32)
        a = AminoAcidUsage('AAA')
        c = CodonUsage('UUUUUCUUAUUG')
        c.normalize()
        self.assertEqual(a.codons(code), c)

        #check that it works with unequal codon frequencies
        unequal = CodonUsage({'GGU':5,'GGC':2,'GGA':2,'GGG':1,'UUU':3,'UUC':1})
        a = AminoAcidUsage('GFFF')
        exp = {
            'GGU':0.5*0.25,
            'GGC':0.2*0.25,
            'GGA':0.2*0.25,
            'GGG':0.1*0.25,
            'UUU':0.75*0.75,
            'UUC':0.25*0.75
        }
        obs = a.codons(codon_usage=unequal)
        for codon, freq in obs.items():
            self.assertFloatEqual(freq, exp.get(codon, 0))
    
    def test_positionalBases(self):
        """AminoAcidUsage positionalBases should return best positional bases"""
        a = AminoAcidUsage('WQRSFADDQW')
        exp = a.codons().positionalBases()
        obs = a.positionalBases()
        for o, e in zip(obs, exp):
            self.assertFloatEqual(o, e)

    def test_aminoAcids(self):
        """AminoAcidUsage aminoAcids should return same object"""
        a = AminoAcidUsage('REWQDFTDSF')
        b = a.aminoAcids()
        assert a is b

class DinucUsageTests(TestCase):
    """Tests of the DinucUsage class."""
    
    def test_init_from_seq(self):
        """DinucUsage should init correctly from string."""
        s1 = 'AAAAA'
        s2 = 'ACTACG'
        fd = filter_dict
        self.assertEqual(fd(DinucUsage(s1)), {'AA':4})
        #NOTE: will map DNA seq tp RNA.
        self.assertEqual(fd(DinucUsage(s2)), {'AC':2,'CU':1,'UA':1,'CG':1})
        #check that it works for non-overlapping
        self.assertEqual(fd(DinucUsage(s1, Overlapping=False)), {'AA':2})
        self.assertEqual(fd(DinucUsage(s2, Overlapping=False)), \
            {'AC':1,'UA':1,'CG':1})
        #check that it works for the 3-1 case
        self.assertEqual(fd(DinucUsage(s1, Overlapping='3-1')), {'AA':1})
        self.assertEqual(fd(DinucUsage(s2, Overlapping='3-1')), \
            {'UA':1})
        s3 = 'ACG'*5
        self.assertEqual(fd(DinucUsage(s3, Overlapping='3-1')), \
            {'GA':4})
        s4 = s3 + 'GAA'
        self.assertEqual(fd(DinucUsage(s4, Overlapping='3-1')), \
            {'GA':4,'GG':1})

    def test_distance(self):
        """Dinuc distance should calculate Euclidean dist. correctly"""
        s1 ='AA'+'GG'*10
        s2 = 'AA'*5 + 'GG'*7
        d1 = DinucUsage(s1, Overlapping=False)
        d2 = DinucUsage(s2, Overlapping=False)
        self.assertEqual(d1.distance(d1), 0)
        self.assertEqual(d1.distance(d2), 5)
        self.assertEqual(d2.distance(d1), 5)

def filter_dict(d):
    """Removes zero keys from dict-like object."""
    result = dict(d)
    for k, v in d.items():
        if not v:
            del result[k]
    return result
        

#run if called from command-line
if __name__ == '__main__':
    main()
