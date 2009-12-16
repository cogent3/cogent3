#!/usr/bin/env python

from cogent.core import moltype, sequence
from cogent.core.moltype import AlphabetError, \
    CoreObjectGroup, AlphabetGroup, make_matches, make_pairs, \
    array, MolType, RNA, DNA, PROTEIN, STANDARD_CODON,\
    IUPAC_RNA_chars, \
    IUPAC_RNA_ambiguities, IUPAC_RNA_ambiguities_complements, \
    IUPAC_DNA_chars, IUPAC_DNA_ambiguities, IUPAC_DNA_ambiguities_complements, \
    RnaStandardPairs, DnaStandardPairs

from cogent.util.unit_test import TestCase, main
from cogent.data.molecular_weight import DnaMW, RnaMW, ProteinMW

__author__ = "Gavin Huttley, Peter Maxwell, and Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

#ind some of the standard alphabets to reduce typing
RnaBases = RNA.Alphabets.Base
DnaBases = DNA.Alphabets.Base
AminoAcids = PROTEIN.Alphabets.Base

#the following classes are to preserve compatibility for older test code
#that assumes mixed-case is OK.
RnaMolType = MolType(
    Sequence = sequence.RnaSequence,
    motifset = IUPAC_RNA_chars,
    Ambiguities = IUPAC_RNA_ambiguities,
    label = "rna_with_lowercase",
    MWCalculator = RnaMW,
    Complements = IUPAC_RNA_ambiguities_complements,
    Pairs = RnaStandardPairs,
    add_lower=True,
    preserve_existing_moltypes=True,
    make_alphabet_group=True,
    )
DnaMolType = MolType(
    Sequence = sequence.DnaSequence,
    motifset = IUPAC_DNA_chars,
    Ambiguities = IUPAC_DNA_ambiguities,
    label = "dna_with_lowercase",
    MWCalculator = DnaMW,
    Complements = IUPAC_DNA_ambiguities_complements,
    Pairs = DnaStandardPairs,
    add_lower=True,
    preserve_existing_moltypes=True,
    make_alphabet_group=True,
    )
ProteinMolType = PROTEIN

class make_matches_tests(TestCase):
    """Tests of the make_matches top-level function"""
    def test_init_empty(self):
        """make_matches should init ok with no parameters"""
        self.assertEqual(make_matches(), {})

    def test_init_monomers(self):
        """make_matches with only monomers should produce {(i,i):True}"""
        m = make_matches('')
        self.assertEqual(m, {})
        m = make_matches('qaz')
        self.assertEqual(m, {('q','q'):True,('a','a'):True,('z','z'):True})
    
    def test_init_gaps(self):
        """make_matches with only gaps should match all gaps to each other"""
        m = make_matches('', '~!')
        self.assertEqual(m, {('~','~'):True,('!','!'):True,('!','~'):True,
            ('~','!'):True})

    def test_init_degen(self):
        """make_matches with only degen should work as expected"""
        m = make_matches(None, None, {'x':'ab','y':'bc','z':'cd', 'n':'bcd'})
        self.assertEqual(m, {('x','x'):False, ('x','y'):False, ('x','n'):False,
            ('y','x'):False, ('y','y'):False, ('y','z'):False, ('y','n'):False,
            ('z','y'):False, ('z','z'):False, ('z','n'):False, ('n','x'):False,
            ('n','y'):False, ('n','z'):False, ('n','n'):False})
        self.assertNotContains(m, ('x','z'))

    def test_init_all(self):
        """make_matches with everything should produce correct dict"""
        m = make_matches('ABC',('-','~'),{'X':'AB','Y':('B','C'),'N':list('ABC')})
        exp = {
            ('-','-'):True,
            ('~','~'):True,
            ('-','~'):True,
            ('~','-'):True,
            ('A','A'):True,
            ('B','B'):True,
            ('C','C'):True,
            ('A','X'):False,
            ('X','A'):False,
            ('B','X'):False,
            ('X','B'):False,
            ('B','Y'):False,
            ('Y','B'):False,
            ('C','Y'):False,
            ('Y','C'):False,
            ('A','N'):False,
            ('N','A'):False,
            ('B','N'):False,
            ('N','B'):False,
            ('C','N'):False,
            ('N','C'):False,
            ('X','X'):False,
            ('Y','Y'):False,
            ('N','N'):False,
            ('X','Y'):False,
            ('Y','X'):False,
            ('X','N'):False,
            ('N','X'):False,
            ('Y','N'):False,
            ('N','Y'):False,
            }
        self.assertEqual(m, exp)

class make_pairs_tests(TestCase):
    """Tests of the top-level make_pairs factory function."""
    def setUp(self):
        """Define some standard pairs and other data"""
        self.pairs = {('U','A'):True, ('A','U'):True, ('G','U'):False}
    
    def test_init_empty(self):
        """make_pairs should init ok with no parameters"""
        self.assertEqual(make_pairs(), {})

    def test_init_pairs(self):
        """make_pairs with just pairs should equal the original"""
        self.assertEqual(make_pairs(self.pairs), self.pairs)
        self.assertNotSameObj(make_pairs(self.pairs), self.pairs)
    
    def test_init_monomers(self):
        """make_pairs with pairs and monomers should equal just the pairs"""
        self.assertEqual(make_pairs(self.pairs, 'ABCDEFG'), self.pairs)
        self.assertNotSameObj(make_pairs(self.pairs, 'ABCDEFG'), self.pairs)

    def test_init_gaps(self):
        """make_pairs should add all combinations of gaps as weak pairs"""
        p = make_pairs(self.pairs, None, '-~')
        self.assertNotEqual(p, self.pairs)
        self.pairs.update({('~','~'):False,('-','~'):False,('-','-'):False,
            ('~','-'):False})
        self.assertEqual(p, self.pairs)

    def test_init_degen(self):
        """make_pairs should add in degenerate combinations as weak pairs"""
        p = make_pairs(self.pairs, 'AUG','-', {'R':'AG','Y':'CU','W':'AU'})
        self.assertNotEqual(p, self.pairs)
        self.pairs.update({
        ('-','-'):False,
        ('A','Y'):False,
        ('Y','A'):False,
        ('A','W'):False,
        ('W','A'):False,
        ('U','R'):False,
        ('R','U'):False,
        ('U','W'):False,
        ('W','U'):False,
        ('G','Y'):False,
        ('G','W'):False,
        ('R','Y'):False,
        ('R','W'):False,
        ('Y','R'):False,
        ('Y','W'):False,
        ('W','R'):False,
        ('W','Y'):False,
        ('W','W'):False,
        })
        self.assertEqual(p, self.pairs)

class CoreObjectGroupTests(TestCase):
    """Tests of the CoreObjectGroup class."""
    
    def test_init(self):
        """CoreObjectGroup should init with basic list of objects."""
        class o(object):
            def __init__(self, s):
                self.s = s
        base = o('base')
        c = CoreObjectGroup(base)
        self.assertSameObj(c.Base, base)
        self.assertSameObj(c.Degen, None)
        self.assertSameObj(c.Base.Degen, None)

        base, degen, gap, degengap = map(o, ['base','degen','gap','degengap'])
        c = CoreObjectGroup(base, degen, gap, degengap)
        self.assertSameObj(c.Base, base)
        self.assertSameObj(c.Base.Degen, degen)
        self.assertSameObj(c.Degen.Gapped, degengap)

class AlphabetGroupTests(TestCase):
    """Tests of the AlphabetGroup class."""

    def test_init(self):
        """AlphabetGroup should initialize successfully"""
        chars = 'AB'
        degens = {'C':'AB'}
        g = AlphabetGroup(chars, degens)
        self.assertEqual(''.join(g.Base), 'AB')
        self.assertEqual(''.join(g.Degen), 'ABC')
        self.assertEqual(''.join(g.Gapped), 'AB-')
        self.assertEqual(''.join(g.DegenGapped), 'AB-C?')

class MolTypeTests(TestCase):
    """Tests of the MolType class. Should support same API as old Alphabet."""
    def test_init_minimal(self):
        """MolType should init OK with just monomers"""
        a = MolType('Abc')
        self.assertContains(a.Alphabet, 'A')
        self.assertNotContains(a.Alphabet, 'a')    # case-sensitive
        self.assertContains(a.Alphabet, 'b')
        self.assertNotContains(a.Alphabet, 'B')
        self.assertNotContains(a.Alphabet, 'x')

    def test_init_everything(self):
        """MolType should init OK with all parameters set"""
        k = dict.fromkeys
        a = MolType(k('Abc'), Ambiguities={'d':'bc'}, Gaps=k('~'), \
            Complements={'b':'c','c':'b'}, Pairs={}, add_lower=False)
        for i in 'Abcd~':
            self.assertContains(a, i)
        self.assertEqual(a.complement('b'), 'c')
        self.assertEqual(a.complement('AbcAA'), 'AcbAA')
        self.assertEqual(a.firstDegenerate('AbcdA'), 3)
        self.assertEqual(a.firstGap('a~c'), 1)
        self.assertEqual(a.firstInvalid('Abcx'), 3)

    def test_stripDegenerate(self):
        """MolType stripDegenerate should remove any degenerate bases"""
        s = RnaMolType.stripDegenerate
        self.assertEqual(s('UCAG-'), 'UCAG-')
        self.assertEqual(s('NRYSW'), '')
        self.assertEqual(s('USNG'), 'UG')

    def test_stripBad(self):
        """MolType stripBad should remove any non-base, non-gap chars"""
        s = RnaMolType.stripBad
        self.assertEqual(s('UCAGwsnyrHBND-D'), 'UCAGwsnyrHBND-D')
        self.assertEqual(s('@#^*($@!#&()!@QZX'), '')
        self.assertEqual(s('aaa ggg ---!ccc'), 'aaaggg---ccc')

    def test_stripBadAndGaps(self):
        """MolType stripBadAndGaps should remove gaps and bad chars"""
        s = RnaMolType.stripBadAndGaps
        self.assertEqual(s('UCAGwsnyrHBND-D'), 'UCAGwsnyrHBNDD')
        self.assertEqual(s('@#^*($@!#&()!@QZX'), '')
        self.assertEqual(s('aaa ggg ---!ccc'), 'aaagggccc')

    def test_complement(self):
        """MolType complement should correctly complement sequence"""
        self.assertEqual(RnaMolType.complement('UauCG-NR'), 'AuaGC-NY')
        self.assertEqual(DnaMolType.complement('TatCG-NR'), 'AtaGC-NY')
        self.assertEqual(RnaMolType.complement(''), '')
        self.assertRaises(TypeError, ProteinMolType.complement, 'ACD')
        #if it wasn't a string, result should be a list
        self.assertEqual(RnaMolType.complement(list('UauCG-NR')), 
            list('AuaGC-NY'))
        self.assertEqual(RnaMolType.complement(('a','c')), ('u','g')) 
        #constructor should fail for a dict
        self.assertRaises(ValueError, RnaMolType.complement, {'a':'c'})

    def test_rc(self):
        """MolType rc should correctly reverse-complement sequence"""
        self.assertEqual(RnaMolType.rc('U'), 'A')
        self.assertEqual(RnaMolType.rc(''), '')
        self.assertEqual(RnaMolType.rc('R'), 'Y')
        self.assertEqual(RnaMolType.rc('UauCG-NR'), 'YN-CGauA')
        self.assertEqual(DnaMolType.rc('TatCG-NR'), 'YN-CGatA')
        self.assertRaises(TypeError, ProteinMolType.rc, 'ACD')
        #if it wasn't a string, result should be a list
        self.assertEqual(RnaMolType.rc(list('UauCG-NR')),
            list('YN-CGauA'))
        self.assertEqual(RnaMolType.rc(('a','c')), ('g','u'))
        #constructor should fail for a dict
        self.assertRaises(ValueError, RnaMolType.rc, {'a':'c'})

    def test_contains(self):
        """MolType contains should return correct result"""
        for i in 'UCAGWSMKRYBDHVN-' + 'UCAGWSMKRYBDHVN-'.lower():
            self.assertContains(RnaMolType, i)
        for i in 'x!@#$%^&ZzQq':
            self.assertNotContains(RnaMolType, i)

        a = MolType(dict.fromkeys('ABC'), add_lower=True)
        for i in 'abcABC':
            self.assertContains(a, i)
        self.assertNotContains(a, 'x')
        b = MolType(dict.fromkeys('ABC'), add_lower=False)
        for i in 'ABC':
            self.assertContains(b, i)
        for i in 'abc':
            self.assertNotContains(b, i)

    def test_iter(self):
       """MolType iter should iterate over monomer order"""
       self.assertEqual(list(RnaMolType), ['U','C','A','G', 'u','c','a','g'])
       a = MolType('ZXCV')
       self.assertEqual(list(a), ['Z','X','C','V'])

    def test_isGapped(self):
        """MolType isGapped should return True if gaps in seq"""
        g = RnaMolType.isGapped
        self.assertFalse(g(''))
        self.assertFalse(g('ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN'))
        self.assertTrue(g('-'))
        self.assertTrue(g('--'))
        self.assertTrue(g('CAGUCGUACGUCAGUACGUacucauacgac-caguACUG'))
        self.assertTrue(g('CA--CGUAUGCA-----g'))
        self.assertTrue(g('CAGU-'))

    def test_isGap(self):
        """MolType isGap should return True if char is a gap"""
        g = RnaMolType.isGap
        #True for the empty string
        self.assertFalse(g(''))
        #True for all the standard and degenerate symbols
        s = 'ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN'
        self.assertFalse(g(s))
        for i in s:
            self.assertFalse(g(i))
        #should be true for a single gap
        self.assertTrue(g('-'))
        #note that it _shouldn't_ be true for a run of gaps: use a.isGapped()
        self.assertFalse(g('--'))

    def test_isDegenerate(self):
        """MolType isDegenerate should return True if degen symbol in seq"""
        d = RnaMolType.isDegenerate
        self.assertFalse(d(''))
        self.assertFalse(d('UACGCUACAUGuacgucaguGCUAGCUA---ACGUCAG'))
        self.assertTrue(d('N'))
        self.assertTrue(d('R'))
        self.assertTrue(d('y'))
        self.assertTrue(d('GCAUguagcucgUCAGUCAGUACgUgcasCUAG'))
        self.assertTrue(d('ACGYAUGCUGYEWEWNFMNfuwbybcwuybcjwbeiwfub'))

    def test_isValid(self):
        """MolType isValid should return True if any unknown symbol in seq"""
        v = RnaMolType.isValid
        self.assertFalse(v(3))
        self.assertFalse(v(None))
        self.assertTrue(v('ACGUGCAUGUCAYCAYGUACGcaugacyugc----RYNCYRNC'))
        self.assertTrue(v(''))
        self.assertTrue(v('a'))
        self.assertFalse(v('ACIUBHFWUIXZKLNJUCIHBICNSOWMOINJ'))
        self.assertFalse(v('CAGUCAGUCACA---GACCAUG-_--cgau'))

    def test_isStrict(self):
        """MolType isStrict should return True if all symbols in Monomers"""
        s = RnaMolType.isStrict
        self.assertFalse(s(3))
        self.assertFalse(s(None))
        self.assertTrue(s(''))
        self.assertTrue(s('A'))
        self.assertTrue(s('UAGCACUgcaugcauGCAUGACuacguACAUG'))
        self.assertFalse(s('CAGUCGAUCA-cgaucagUCGAUGAC'))
        self.assertFalse(s('ACGUGCAUXCAGUCAG'))

    def test_firstGap(self):
        """MolType firstGap should return index of first gap symbol, or None"""
        g = RnaMolType.firstGap
        self.assertEqual(g(''), None)
        self.assertEqual(g('a'), None)
        self.assertEqual(g('uhacucHuhacUIUIhacan'), None)
        self.assertEqual(g('-abc'), 0)
        self.assertEqual(g('b-ac'), 1)
        self.assertEqual(g('abcd-'), 4)

    def test_firstDegenerate(self):
        """MolType firstDegenerate should return index of first degen symbol"""
        d = RnaMolType.firstDegenerate
        self.assertEqual(d(''), None)
        self.assertEqual(d('a'), None)
        self.assertEqual(d('UCGACA--CU-gacucaguacgua'), None)
        self.assertEqual(d('nCAGU'), 0)
        self.assertEqual(d('CUGguagvAUG'), 7)
        self.assertEqual(d('ACUGCUAacgud'), 11)

    def test_firstInvalid(self):
        """MolType firstInvalid should return index of first invalid symbol"""
        i = RnaMolType.firstInvalid
        self.assertEqual(i(''), None)
        self.assertEqual(i('A'), None)
        self.assertEqual(i('ACGUNVBuacg-wskmWSMKYRryNnN--'), None)
        self.assertEqual(i('x'), 0)
        self.assertEqual(i('rx'), 1)
        self.assertEqual(i('CAGUNacgunRYWSwx'), 15)

    def test_firstNonStrict(self):
        """MolType firstNonStrict should return index of first non-strict symbol"""
        s = RnaMolType.firstNonStrict
        self.assertEqual(s(''), None)
        self.assertEqual(s('A'), None)
        self.assertEqual(s('ACGUACGUcgaucagu'), None)
        self.assertEqual(s('N'), 0)
        self.assertEqual(s('-'), 0)
        self.assertEqual(s('x'), 0)
        self.assertEqual(s('ACGUcgAUGUGCAUcaguX'), 18)
        self.assertEqual(s('ACGUcgAUGUGCAUcaguX-38243829'), 18)

    def test_disambiguate(self):
        """MolType disambiguate should remove degenerate bases"""
        d = RnaMolType.disambiguate
        self.assertEqual(d(''), '')
        self.assertEqual(d('AGCUGAUGUA--CAGU'),'AGCUGAUGUA--CAGU')
        self.assertEqual(d('AUn-yrs-wkmCGwmrNMWRKY', 'strip'), 'AU--CG')
        self.assertEqual(d(tuple('AUn-yrs-wkmCGwmrNMWRKY'), 'strip'), \
            tuple('AU--CG'))
        s = 'AUn-yrs-wkmCGwmrNMWRKY'
        t = d(s, 'random')
        u = d(s, 'random')
        for i, j in zip(s, t):
            if i in RnaMolType.Degenerates:
                self.assertContains(RnaMolType.Degenerates[i], j)
            else:
                self.assertEquals(i, j)
        self.assertNotEqual(t, u)
        self.assertEqual(d(tuple('UCAG'), 'random'), tuple('UCAG'))
        self.assertEqual(len(s), len(t))
        self.assertSameObj(RnaMolType.firstDegenerate(t), None)
        #should raise exception on unknown disambiguation method
        self.assertRaises(NotImplementedError, d, s, 'xyz')

    def test_degap(self):
        """MolType degap should remove all gaps from sequence"""
        g = RnaMolType.degap
        self.assertEqual(g(''), '')
        self.assertEqual(g('GUCAGUCgcaugcnvuincdks'), 'GUCAGUCgcaugcnvuincdks')
        self.assertEqual(g('----------------'), '')
        self.assertEqual(g('gcuauacg-'), 'gcuauacg')
        self.assertEqual(g('-CUAGUCA'), 'CUAGUCA')
        self.assertEqual(g('---a---c---u----g---'), 'acug')
        self.assertEqual(g(tuple('---a---c---u----g---')), tuple('acug'))
        
    def test_gapList(self):
        """MolType gapList should return correct gap positions"""
        g = RnaMolType.gapList
        self.assertEqual(g(''), [])
        self.assertEqual(g('ACUGUCAGUACGHFSDKJCUICDNINS'), [])
        self.assertEqual(g('GUACGUIACAKJDC-SDFHJDSFK'), [14])
        self.assertEqual(g('-DSHFUHDSF'), [0])
        self.assertEqual(g('UACHASJAIDS-'), [11])
        self.assertEqual(g('---CGAUgCAU---ACGHc---ACGUCAGU---'), \
            [0,1,2,11,12,13,19,20,21,30,31,32])
        a = MolType({'A':1}, Gaps=dict.fromkeys('!@#$%'))
        g = a.gapList
        self.assertEqual(g(''), [])
        self.assertEqual(g('!!!'), [0,1,2])
        self.assertEqual(g('!@#$!@#$!@#$'), range(12))
        self.assertEqual(g('cguua!cgcuagua@cguasguadc#'), [5,14,25])

    def test_gapVector(self):
        """MolType gapVector should return correct gap positions"""
        g = RnaMolType.gapVector
        self.assertEqual(g(''), [])
        self.assertEqual(g('ACUGUCAGUACGHFSDKJCUICDNINS'), [False]*27)
        self.assertEqual(g('GUACGUIACAKJDC-SDFHJDSFK'), 
         map(bool, map(int,'000000000000001000000000')))
        self.assertEqual(g('-DSHFUHDSF'), 
         map(bool, map(int,'1000000000')))
        self.assertEqual(g('UACHASJAIDS-'), 
         map(bool, map(int,'000000000001')))
        self.assertEqual(g('---CGAUgCAU---ACGHc---ACGUCAGU---'), \
         map(bool, map(int,'111000000001110000011100000000111')))
        a = MolType({'A':1}, Gaps=dict.fromkeys('!@#$%'))
        g = a.gapVector
        self.assertEqual(g(''), [])
        self.assertEqual(g('!!!'), map(bool, [1,1,1]))
        self.assertEqual(g('!@#$!@#$!@#$'), [True] * 12)
        self.assertEqual(g('cguua!cgcuagua@cguasguadc#'), 
         map(bool, map(int,'00000100000000100000000001')))

    def test_gapMaps(self):
        """MolType gapMaps should return dicts mapping gapped/ungapped pos"""
        empty = ''
        no_gaps = 'aaa'
        all_gaps = '---'
        start_gaps = '--abc'
        end_gaps = 'ab---'
        mid_gaps = '--a--b-cd---'
        gm = RnaMolType.gapMaps
        self.assertEqual(gm(empty), ({},{}))
        self.assertEqual(gm(no_gaps), ({0:0,1:1,2:2}, {0:0,1:1,2:2}))
        self.assertEqual(gm(all_gaps), ({},{}))
        self.assertEqual(gm(start_gaps), ({0:2,1:3,2:4},{2:0,3:1,4:2}))
        self.assertEqual(gm(end_gaps), ({0:0,1:1},{0:0,1:1}))
        self.assertEqual(gm(mid_gaps), ({0:2,1:5,2:7,3:8},{2:0,5:1,7:2,8:3}))
        
    def test_countGaps(self):
        """MolType countGaps should return correct gap count"""
        c = RnaMolType.countGaps
        self.assertEqual(c(''), 0)
        self.assertEqual(c('ACUGUCAGUACGHFSDKJCUICDNINS'), 0)
        self.assertEqual(c('GUACGUIACAKJDC-SDFHJDSFK'), 1)
        self.assertEqual(c('-DSHFUHDSF'), 1)
        self.assertEqual(c('UACHASJAIDS-'), 1)
        self.assertEqual(c('---CGAUgCAU---ACGHc---ACGUCAGU---'), 12)
        a = MolType({'A':1}, Gaps=dict.fromkeys('!@#$%'))
        c = a.countGaps
        self.assertEqual(c(''), 0)
        self.assertEqual(c('!!!'), 3)
        self.assertEqual(c('!@#$!@#$!@#$'), 12)
        self.assertEqual(c('cguua!cgcuagua@cguasguadc#'), 3)

    def test_countDegenerate(self):
        """MolType countDegenerate should return correct degen base count"""
        d = RnaMolType.countDegenerate
        self.assertEqual(d(''), 0)
        self.assertEqual(d('GACUGCAUGCAUCGUACGUCAGUACCGA'), 0)
        self.assertEqual(d('N'), 1)
        self.assertEqual(d('NRY'), 3)
        self.assertEqual(d('ACGUAVCUAGCAUNUCAGUCAGyUACGUCAGS'), 4)

    def test_possibilites(self):
        """MolType possibilities should return correct # possible sequences"""
        p = RnaMolType.possibilities
        self.assertEqual(p(''), 1)
        self.assertEqual(p('ACGUgcaucagUCGuGCAU'), 1)
        self.assertEqual(p('N'), 4)
        self.assertEqual(p('R'), 2)
        self.assertEqual(p('H'), 3)
        self.assertEqual(p('nRh'), 24)
        self.assertEqual(p('AUGCnGUCAg-aurGauc--gauhcgauacgws'), 96)
       
    def test_MW(self):
        """MolType MW should return correct molecular weight"""
        r = RnaMolType.MW
        p = ProteinMolType.MW
        self.assertEqual(p(''), 0)
        self.assertEqual(r(''), 0)
        self.assertFloatEqual(p('A'), 107.09)
        self.assertFloatEqual(r('A'), 375.17)
        self.assertFloatEqual(p('AAA'), 285.27)
        self.assertFloatEqual(r('AAA'), 1001.59)
        self.assertFloatEqual(r('AAACCCA'), 2182.37)

    def test_canMatch(self):
        """MolType canMatch should return True if all positions can match"""
        m = RnaMolType.canMatch
        self.assertTrue(m('', ''))
        self.assertTrue(m('UCAG', 'UCAG'))
        self.assertFalse(m('UCAG', 'ucag'))
        self.assertTrue(m('UCAG', 'NNNN'))
        self.assertTrue(m('NNNN', 'UCAG'))
        self.assertTrue(m('NNNN', 'NNNN'))
        self.assertFalse(m('N', 'x'))
        self.assertFalse(m('N', '-'))
        self.assertTrue(m('UCAG', 'YYRR'))
        self.assertTrue(m('UCAG', 'KMWS'))

    def test_canMismatch(self):
        """MolType canMismatch should return True on any possible mismatch"""
        m = RnaMolType.canMismatch
        self.assertFalse(m('',''))
        self.assertTrue(m('N', 'N'))
        self.assertTrue(m('R', 'R'))
        self.assertTrue(m('N', 'r'))
        self.assertTrue(m('CGUACGCAN', 'CGUACGCAN'))
        self.assertTrue(m('U', 'C'))
        self.assertTrue(m('UUU', 'UUC'))
        self.assertTrue(m('UUU', 'UUY'))
        self.assertFalse(m('UUU', 'UUU'))
        self.assertFalse(m('UCAG', 'UCAG'))
        self.assertFalse(m('U--', 'U--'))

    def test_mustMatch(self):
        """MolType mustMatch should return True when no possible mismatches"""
        m = RnaMolType.mustMatch
        self.assertTrue(m('',''))
        self.assertFalse(m('N', 'N'))
        self.assertFalse(m('R', 'R'))
        self.assertFalse(m('N', 'r'))
        self.assertFalse(m('CGUACGCAN', 'CGUACGCAN'))
        self.assertFalse(m('U', 'C'))
        self.assertFalse(m('UUU', 'UUC'))
        self.assertFalse(m('UUU', 'UUY'))
        self.assertTrue(m('UU-', 'UU-'))
        self.assertTrue(m('UCAG', 'UCAG'))

    def test_canPair(self):
        """MolType canPair should return True if all positions can pair"""
        p = RnaMolType.canPair
        self.assertTrue(p('', ''))
        self.assertFalse(p('UCAG', 'UCAG'))
        self.assertTrue(p('UCAG', 'CUGA'))
        self.assertFalse(p('UCAG', 'cuga'))
        self.assertTrue(p('UCAG', 'NNNN'))
        self.assertTrue(p('NNNN', 'UCAG'))
        self.assertTrue(p('NNNN', 'NNNN'))
        self.assertFalse(p('N', 'x'))
        self.assertFalse(p('N', '-'))
        self.assertTrue(p('-', '-'))
        self.assertTrue(p('UCAGU', 'KYYRR'))
        self.assertTrue(p('UCAG', 'KKRS'))
        self.assertTrue(p('U', 'G'))

        d = DnaMolType.canPair
        self.assertFalse(d('T', 'G'))

    def test_canMispair(self):
        """MolType canMispair should return True on any possible mispair"""
        m = RnaMolType.canMispair
        self.assertFalse(m('',''))
        self.assertTrue(m('N', 'N'))
        self.assertTrue(m('R', 'Y'))
        self.assertTrue(m('N', 'r'))
        self.assertTrue(m('CGUACGCAN', 'NUHCHUACH'))
        self.assertTrue(m('U', 'C'))
        self.assertTrue(m('U', 'R'))
        self.assertTrue(m('UUU', 'AAR'))
        self.assertTrue(m('UUU', 'GAG'))
        self.assertFalse(m('UUU', 'AAA'))
        self.assertFalse(m('UCAG', 'CUGA'))
        self.assertTrue(m('U--', '--U'))

        d = DnaMolType.canPair
        self.assertTrue(d('TCCAAAGRYY', 'RRYCTTTGGA'))

    def test_mustPair(self):
        """MolType mustPair should return True when no possible mispairs"""
        m = RnaMolType.mustPair
        self.assertTrue(m('',''))
        self.assertFalse(m('N', 'N'))
        self.assertFalse(m('R', 'Y'))
        self.assertFalse(m('A', 'A'))
        self.assertFalse(m('CGUACGCAN', 'NUGCGUACG'))
        self.assertFalse(m('U', 'C'))
        self.assertFalse(m('UUU', 'AAR'))
        self.assertFalse(m('UUU', 'RAA'))
        self.assertFalse(m('UU-', '-AA'))
        self.assertTrue(m('UCAG', 'CUGA'))

        d = DnaMolType.mustPair
        self.assertTrue(d('TCCAGGG', 'CCCTGGA'))
        self.assertTrue(d('tccaggg', 'ccctgga'))
        self.assertFalse(d('TCCAGGG', 'NCCTGGA'))

class RnaMolTypeTests(TestCase):
    """Spot-checks of alphabet functionality applied to RNA alphabet."""
    
    def test_contains(self):
        """RnaMolType should __contain__ the expected symbols."""
        keys = 'ucagrymkwsbhvdn?-'
        for k in keys:
            self.assertContains(RnaMolType, k)
        for k in keys.upper():
            self.assertContains(RnaMolType, k)
        self.assertNotContains(RnaMolType, 'X')

    def test_degenerateFromSequence(self):
        """RnaMolType degenerateFromSequence should give correct results"""
        d = RnaMolType.degenerateFromSequence
        #check monomers
        self.assertEqual(d('a'), 'a')
        self.assertEqual(d('C'), 'C')
        #check seq of monomers
        self.assertEqual(d('aaaaa'), 'a')
        #check some 2- to 4-way cases
        self.assertEqual(d('aaaaag'), 'r')
        self.assertEqual(d('ccgcgcgcggcc'), 's')
        self.assertEqual(d('accgcgcgcggcc'), 'v')
        self.assertEqual(d('aaaaagcuuu'), 'n')
        #check some cases with gaps
        self.assertEqual(d('aa---aaagcuuu'), '?')
        self.assertEqual(d('aaaaaaaaaaaaaaa-'), '?')
        self.assertEqual(d('----------------'), '-')
        #check mixed case example
        self.assertEqual(d('AaAAaa'), 'A')
        #check example with degenerate symbols in set
        self.assertEqual(d('RS'), 'V')
        self.assertEqual(d('RN-'), '?')
        #check that it works for proteins as well
        p = ProteinMolType.degenerateFromSequence
        self.assertEqual(p('A'), 'A')
        self.assertEqual(p('AAA'), 'A')
        self.assertEqual(p('DN'), 'B')
        self.assertEqual(p('---'), '-')
        self.assertEqual(p('ACD'), 'X')
        self.assertEqual(p('ABD'), 'X')
        self.assertEqual(p('ACX'), 'X')
        self.assertEqual(p('AC-'), '?')


class _AlphabetTestCase(TestCase):
    def assertEqualSeqs(self, a, b):
        """For when we don't care about the type, just the elements"""
        self.assertEqual(list(a), list(b))
        
    def assertEqualSets(self, a, b):
        self.assertEqual(set(a), set(b))
        

class DNAAlphabet(_AlphabetTestCase):
    def setUp(self):
        self.alpha = DNA.Alphabet
        
    def test_exclude(self):
        """Nucleotide alphabet testing excluding gap motif"""
        self.assertEqualSeqs(self.alpha, ['T','C','A','G'])
        
    def test_include(self):
        """Nucleotide alphabet testing including gap motif"""
        self.assertEqualSets(self.alpha.withGapMotif(), ['A','C','G','T','-'])
        
    def test_usesubset(self):
        """testing using a subset of motifs."""
        self.assertEqualSets(self.alpha.withGapMotif(), ['A','C','G','T','-'])
        alpha = self.alpha.getSubset(motif_subset = ['A'])
        self.assertEqualSets(alpha, ['A'])
        #self.assertRaises(AlphabetError, self.alpha.getSubset, ['A','C'])
        alpha = DNA.Alphabet
        self.assertEqualSets(alpha, ['T','C','A','G'])
        alpha = alpha.getSubset(motif_subset = ['A','T','G'])
        self.assertEqualSets(alpha, ['A','G','T'])
        
class DinucAlphabet(_AlphabetTestCase):
    def setUp(self):
        self.alpha = DNA.Alphabet.withGapMotif().getWordAlphabet(2)
    
    def test_exclude(self):
        """Dinucleotide alphabet testing excluding gap motif"""
        expected = ['-A', '-C', '-G', '-T',
                    'A-', 'AA', 'AC', 'AG', 'AT',
                    'C-', 'CA', 'CC', 'CG', 'CT',
                    'G-', 'GA', 'GC', 'GG', 'GT',
                    'T-', 'TA', 'TC', 'TG', 'TT']
                    
        self.assertEqualSets(
                self.alpha.getSubset(['--'], excluded=True)
                , expected)
        
    def test_include(self):
        """Dinucleotide alphabet testing including gap motif"""
        expected =  ['--', '-A', '-C', '-G', '-T',
                    'A-', 'AA', 'AC', 'AG', 'AT',
                    'C-', 'CA', 'CC', 'CG', 'CT',
                    'G-', 'GA', 'GC', 'GG', 'GT',
                    'T-', 'TA', 'TC', 'TG', 'TT']
        self.assertEqualSets(self.alpha, expected)
        
    def test_usesubset(self):
        """testing using a subset of motifs."""
        alpha = self.alpha.getSubset(motif_subset = ['AA', 'CA','GT'])
        self.assertEqualSeqs(alpha, ['AA', 'CA','GT'])

        self.assertRaises(AlphabetError, alpha.getSubset, motif_subset = ['AA','CA','GT', 'TT'])
    
    def test_usesubsetbyfreq(self):
        """testing using a subset of motifs by using motif probs."""
        motif_freqs = {'--':0, '-A': 0.0, '-C': 0, '-G': 0, '-T': 0,
                       'A-': 0, 'AA': 1.0, 'AC': 0.0, 'AG': 0, 'AT': 0,
                       'C-': 0, 'CA': 1, 'CC': 0, 'CG': 0, 'CT': 0,
                       'G-': 0, 'GA': 0, 'GC': 0, 'GG': 0, 'GT': 1,
                       'T-': 0, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0}
        
        alpha = self.alpha.getSubset(motif_freqs)
        self.assertEqualSets(alpha, ['AA', 'CA', 'GT'])

class CodonAlphabet(_AlphabetTestCase):
    def setUp(self):
        self.alpha = STANDARD_CODON
    
    def test_ambiguous_gaps(self):
        alpha = self.alpha.withGapMotif()
        self.assertEqual(len(alpha.resolveAmbiguity('AT?')), 4)
        self.assertRaises(Exception, alpha.resolveAmbiguity, 'at-')
        self.assertEqual(len(alpha.resolveAmbiguity('???')), 62)
        self.assertEqual(len(alpha.resolveAmbiguity('---')), 1)
        
        alpha = self.alpha
        self.assertEqual(len(alpha.resolveAmbiguity('AT?')), 4)
        self.assertRaises(Exception, alpha.resolveAmbiguity, 'at-')
        self.assertEqual(len(alpha.resolveAmbiguity('???')), 61)
        self.assertRaises(Exception, alpha.resolveAmbiguity, '---')
        

    def test_exclude(self):
        """testing excluding gap motif"""
        alpha = self.alpha
        expected = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA',
                    'ACC', 'ACG', 'ACT', 'AGA', 'AGC',
                    'AGG', 'AGT', 'ATA', 'ATC', 'ATG',
                    'ATT', 'CAA', 'CAC', 'CAG', 'CAT',
                    'CCA', 'CCC', 'CCG', 'CCT', 'CGA',
                    'CGC', 'CGG', 'CGT', 'CTA', 'CTC',
                    'CTG', 'CTT', 'GAA', 'GAC', 'GAG',
                    'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                    'GGA', 'GGC', 'GGG', 'GGT', 'GTA',
                    'GTC', 'GTG', 'GTT', 'TAC', 'TAT',
                    'TCA', 'TCC', 'TCG', 'TCT', 'TGC',
                    'TGG', 'TGT', 'TTA', 'TTC', 'TTG',
                    'TTT']
        self.assertEqualSets(alpha, expected)
        
    def test_include(self):
        """testing including gap motif"""
        alpha = self.alpha.withGapMotif()
        expected = ['---', 'AAA', 'AAC', 'AAG', 'AAT', 'ACA',
                    'ACC', 'ACG', 'ACT', 'AGA', 'AGC',
                    'AGG', 'AGT', 'ATA', 'ATC', 'ATG',
                    'ATT', 'CAA', 'CAC', 'CAG', 'CAT',
                    'CCA', 'CCC', 'CCG', 'CCT', 'CGA',
                    'CGC', 'CGG', 'CGT', 'CTA', 'CTC',
                    'CTG', 'CTT', 'GAA', 'GAC', 'GAG',
                    'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                    'GGA', 'GGC', 'GGG', 'GGT', 'GTA',
                    'GTC', 'GTG', 'GTT', 'TAC', 'TAT',
                    'TCA', 'TCC', 'TCG', 'TCT', 'TGC',
                    'TGG', 'TGT', 'TTA', 'TTC', 'TTG',
                    'TTT']
        self.assertEqualSets(alpha, expected)
  

if __name__ == '__main__':
    main()
