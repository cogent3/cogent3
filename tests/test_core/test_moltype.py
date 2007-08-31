#!/usr/bin/env python

from cogent.core import moltype, sequence
from cogent.core.moltype import AlphabetError, make_matches, make_pairs, \
    array, MolType, RNA, DNA, PROTEIN, STANDARD_CODON,\
    IUPAC_RNA_chars, \
    IUPAC_RNA_ambiguities, IUPAC_RNA_ambiguities_complements, \
    IUPAC_DNA_chars, IUPAC_DNA_ambiguities, IUPAC_DNA_ambiguities_complements, \
    RnaStandardPairs, DnaStandardPairs

from cogent.util.unit_test import TestCase, main
from cogent.data.molecular_weight import DnaMW, RnaMW, ProteinMW

__author__ = "Gavin Huttley, Peter Maxwell, and Rob Knight"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.0.1"
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
        assert ('x','z') not in m

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
        assert make_pairs(self.pairs) is not self.pairs
    
    def test_init_monomers(self):
        """make_pairs with pairs and monomers should equal just the pairs"""
        self.assertEqual(make_pairs(self.pairs, 'ABCDEFG'), self.pairs)
        assert make_pairs(self.pairs, 'ABCDEFG') is not self.pairs

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

class MolTypeTests(TestCase):
    """Tests of the MolType class. Should support same API as old Alphabet."""
    def test_init_minimal(self):
        """MolType should init OK with just monomers"""
        a = MolType('Abc')
        assert 'A' in a.Alphabet
        assert 'a' not in a.Alphabet    # case-sensitive
        assert 'b' in a.Alphabet
        assert 'B' not in a.Alphabet
        assert 'x' not in a.Alphabet

    def test_init_everything(self):
        """MolType should init OK with all parameters set"""
        k = dict.fromkeys
        a = MolType(k('Abc'), Ambiguities={'d':'bc'}, Gaps=k('~'), \
            Complements={'b':'c','c':'b'}, Pairs={}, add_lower=False)
        for i in 'Abcd~':
            assert i in a
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
            assert i in RnaMolType
        for i in 'x!@#$%^&ZzQq':
            assert i not in RnaMolType

        a = MolType(dict.fromkeys('ABC'), add_lower=True)
        for i in 'abcABC':
            assert i in a
        assert 'x' not in a
        b = MolType(dict.fromkeys('ABC'), add_lower=False)
        for i in 'ABC':
            assert i in b
        for i in 'abc':
            assert i not in b

    def test_iter(self):
       """MolType iter should iterate over monomer order"""
       self.assertEqual(list(RnaMolType), ['U','C','A','G', 'u','c','a','g'])
       a = MolType('ZXCV')
       self.assertEqual(list(a), ['Z','X','C','V'])

    def test_isGapped(self):
        """MolType isGapped should return True if gaps in seq"""
        g = RnaMolType.isGapped
        assert not g('')
        assert not g('ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN')
        assert g('-')
        assert g('--')
        assert g('CAGUCGUACGUCAGUACGUacucauacgac-caguACUG')
        assert g('CA--CGUAUGCA-----g')
        assert g('CAGU-')

    def test_isGap(self):
        """MolType isGap should return True if char is a gap"""
        g = RnaMolType.isGap
        #True for the empty string
        assert not g('')
        #True for all the standard and degenerate symbols
        s = 'ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN'
        assert not g(s)
        for i in s:
            assert not g(i)
        #should be true for a single gap
        assert g('-')
        #note that it _shouldn't_ be true for a run of gaps: use a.isGapped()
        assert not g('--')

    def test_isDegenerate(self):
        """MolType isDegenerate should return True if degen symbol in seq"""
        d = RnaMolType.isDegenerate
        assert not d('')
        assert not d('UACGCUACAUGuacgucaguGCUAGCUA---ACGUCAG')
        assert d('N')
        assert d('R')
        assert d('y')
        assert d('GCAUguagcucgUCAGUCAGUACgUgcasCUAG')
        assert d('ACGYAUGCUGYEWEWNFMNfuwbybcwuybcjwbeiwfub')

    def test_isValid(self):
        """MolType isValid should return True if any unknown symbol in seq"""
        v = RnaMolType.isValid
        assert not v(3)
        assert not v(None)
        assert v('ACGUGCAUGUCAYCAYGUACGcaugacyugc----RYNCYRNC')
        assert v('')
        assert v('a')
        assert not v('ACIUBHFWUIXZKLNJUCIHBICNSOWMOINJ')
        assert not v('CAGUCAGUCACA---GACCAUG-_--cgau')

    def test_isStrict(self):
        """MolType isStrict should return True if all symbols in Monomers"""
        s = RnaMolType.isStrict
        assert not s(3)
        assert not s(None)
        assert s('')
        assert s('A')
        assert s('UAGCACUgcaugcauGCAUGACuacguACAUG')
        assert not s('CAGUCGAUCA-cgaucagUCGAUGAC')
        assert not s('ACGUGCAUXCAGUCAG')

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
                assert j in RnaMolType.Degenerates[i]
            else:
                assert i == j
        self.assertNotEqual(t, u)
        self.assertEqual(d(tuple('UCAG'), 'random'), tuple('UCAG'))
        self.assertEqual(len(s), len(t))
        assert RnaMolType.firstDegenerate(t) is None
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
        assert m('', '')
        assert m('UCAG', 'UCAG')
        assert not m('UCAG', 'ucag')
        assert m('UCAG', 'NNNN')
        assert m('NNNN', 'UCAG')
        assert m('NNNN', 'NNNN')
        assert not m('N', 'x')
        assert not m('N', '-')
        assert m('UCAG', 'YYRR')
        assert m('UCAG', 'KMWS')

    def test_canMismatch(self):
        """MolType canMismatch should return True on any possible mismatch"""
        m = RnaMolType.canMismatch
        assert not m('','')
        assert m('N', 'N')
        assert m('R', 'R')
        assert m('N', 'r')
        assert m('CGUACGCAN', 'CGUACGCAN')
        assert m('U', 'C')
        assert m('UUU', 'UUC')
        assert m('UUU', 'UUY')
        assert not m('UUU', 'UUU')
        assert not m('UCAG', 'UCAG')
        assert not m('U--', 'U--')

    def test_mustMatch(self):
        """MolType mustMatch should return True when no possible mismatches"""
        m = RnaMolType.mustMatch
        assert m('','')
        assert not m('N', 'N')
        assert not m('R', 'R')
        assert not m('N', 'r')
        assert not m('CGUACGCAN', 'CGUACGCAN')
        assert not m('U', 'C')
        assert not m('UUU', 'UUC')
        assert not m('UUU', 'UUY')
        assert m('UU-', 'UU-')
        assert m('UCAG', 'UCAG')

    def test_canPair(self):
        """MolType canPair should return True if all positions can pair"""
        p = RnaMolType.canPair
        assert p('', '')
        assert not p('UCAG', 'UCAG')
        assert p('UCAG', 'CUGA')
        assert not p('UCAG', 'cuga')
        assert p('UCAG', 'NNNN')
        assert p('NNNN', 'UCAG')
        assert p('NNNN', 'NNNN')
        assert not p('N', 'x')
        assert not p('N', '-')
        assert p('-', '-')
        assert p('UCAGU', 'KYYRR')
        assert p('UCAG', 'KKRS')
        assert p('U', 'G')

        d = DnaMolType.canPair
        assert not d('T', 'G')

    def test_canMispair(self):
        """MolType canMispair should return True on any possible mispair"""
        m = RnaMolType.canMispair
        assert not m('','')
        assert m('N', 'N')
        assert m('R', 'Y')
        assert m('N', 'r')
        assert m('CGUACGCAN', 'NUHCHUACH')
        assert m('U', 'C')
        assert m('U', 'R')
        assert m('UUU', 'AAR')
        assert m('UUU', 'GAG')
        assert not m('UUU', 'AAA')
        assert not m('UCAG', 'CUGA')
        assert m('U--', '--U')

        d = DnaMolType.canPair
        assert d('TCCAAAGRYY', 'RRYCTTTGGA')

    def test_mustPair(self):
        """MolType mustPair should return True when no possible mispairs"""
        m = RnaMolType.mustPair
        assert m('','')
        assert not m('N', 'N')
        assert not m('R', 'Y')
        assert not m('A', 'A')
        assert not m('CGUACGCAN', 'NUGCGUACG')
        assert not m('U', 'C')
        assert not m('UUU', 'AAR')
        assert not m('UUU', 'RAA')
        assert not m('UU-', '-AA')
        assert m('UCAG', 'CUGA')

        d = DnaMolType.mustPair
        assert d('TCCAGGG', 'CCCTGGA')
        assert d('tccaggg', 'ccctgga')
        assert not d('TCCAGGG', 'NCCTGGA')

class RnaMolTypeTests(TestCase):
    """Spot-checks of alphabet functionality applied to RNA alphabet."""
    
    def test_contains(self):
        """RnaMolType should __contain__ the expected symbols."""
        keys = 'ucagrymkwsbhvdn?-'
        for k in keys:
            assert k in RnaMolType
        for k in keys.upper():
            assert k in RnaMolType
        assert 'X' not in RnaMolType

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
