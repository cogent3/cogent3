#!/usr/bin/env python
"""Unit tests for Sequence class and its subclasses.
"""

from cogent.core.sequence import Sequence, RnaSequence, DnaSequence, \
    ProteinSequence, ModelSequenceBase, \
    ModelSequence, ModelNucleicAcidSequence, ModelRnaSequence, \
    ModelDnaSequence, ModelProteinSequence, ModelCodonSequence, \
    ModelDnaCodonSequence, ModelRnaCodonSequence
from cogent.core.moltype import RNA, DNA, PROTEIN, ASCII, BYTES, AlphabetError
from cogent.util.unit_test import TestCase, main

import re
from pickle import dumps
from numpy import array

__author__ = "Rob Knight, Gavin Huttley and Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Peter Maxwell",
                    "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class SequenceTests(TestCase):
    """Tests of the Sequence class."""
    SEQ = Sequence
    RNA = RnaSequence
    DNA = DnaSequence
    PROT = ProteinSequence
    def test_init_empty(self):
        """Sequence and subclasses should init correctly."""
        #NOTE: ModelSequences can't be initialized empty because it screws up
        #the dimensions of the array, and not worth special-casing.
        s = self.SEQ()
        self.assertEqual(s, '')
        assert s.MolType in (ASCII, BYTES)

        r = self.RNA()
        assert r.MolType is RNA 

    def test_init_data(self):
        """Sequence init with data should set data in correct location"""
        r = self.RNA('ucagg')
        #no longer preserves case
        self.assertEqual(r, 'UCAGG')

    def test_init_other_seq(self):
        """Sequence init with other seq should preserve name and info."""
        r = self.RNA('UCAGG', Name='x', Info={'z':3})
        s = Sequence(r)
        self.assertEqual(s._seq, 'UCAGG')
        self.assertEqual(s.Name, 'x')
        self.assertEqual(s.Info.z, 3)

    def test_compare_to_string(self):
        """Sequence should compare equal to same string."""
        r = self.RNA('UCC')
        self.assertEqual(r, 'UCC')

    def test_slice(self):
        """Sequence slicing should work as expected"""
        r = self.RNA('UCAGG')
        self.assertEqual(r[0], 'U')
        self.assertEqual(r[-1], 'G')
        self.assertEqual(r[1:3], 'CA')

    def test_conversion(self):
        """Should convert t to u automatically"""
        r = self.RNA('TCAtu')
        self.assertEqual(str(r), 'UCAUU')

        d = self.DNA('UCAtu')
        self.assertEqual(str(d), 'TCATT')

    def test_toDna(self):
        """Returns copy of self as DNA."""
        r = self.RNA('TCA')
        self.assertEqual(str(r), 'UCA')
        self.assertEqual(str(r.toDna()), 'TCA')

    def test_toRna(self):
        """Returns copy of self as RNA."""
        r = self.DNA('UCA')
        self.assertEqual(str(r), 'TCA')
        self.assertEqual(str(r.toRna()), 'UCA')


    def test_toFasta(self):
        """Sequence toFasta() should return Fasta-format string"""
        even = 'TCAGAT'
        odd = even + 'AAA'
        even_dna = self.SEQ(even, Name='even')
        odd_dna = self.SEQ(odd, Name='odd')
        self.assertEqual(even_dna.toFasta(), '>even\nTCAGAT')
        #set line wrap to small number so we can test that it works
        even_dna.LineWrap = 2
        self.assertEqual(even_dna.toFasta(), '>even\nTC\nAG\nAT')
        odd_dna.LineWrap = 2
        self.assertEqual(odd_dna.toFasta(), '>odd\nTC\nAG\nAT\nAA\nA')
        #check that changing the linewrap again works
        even_dna.LineWrap  = 4
        self.assertEqual(even_dna.toFasta(), '>even\nTCAG\nAT')

    def test_serialize(self):
        """Sequence should be serializable"""
        r = self.RNA('ugagg')
        assert dumps(r)

    def test_stripDegenerate(self):
        """Sequence stripDegenerate should remove any degenerate bases"""
        self.assertEqual(self.RNA('UCAG-').stripDegenerate(), 'UCAG-')
        self.assertEqual(self.RNA('NRYSW').stripDegenerate(), '')
        self.assertEqual(self.RNA('USNG').stripDegenerate(), 'UG')

    def test_stripBad(self):
        """Sequence stripBad should remove any non-base, non-gap chars"""
        #have to turn off check to get bad data in; no longer preserves case
        self.assertEqual(self.RNA('UCxxxAGwsnyrHBNzzzD-D', check=False\
            ).stripBad(), 'UCAGWSNYRHBND-D')
        self.assertEqual(self.RNA('@#^*($@!#&()!@QZX', check=False \
            ).stripBad(), '')
        self.assertEqual(self.RNA('aaaxggg---!ccc', check=False).stripBad(), 
            'AAAGGG---CCC')

    def test_stripBadAndGaps(self):
        """Sequence stripBadAndGaps should remove gaps and bad chars"""
        #have to turn off check to get bad data in; no longer preserves case
        self.assertEqual(self.RNA('UxxCAGwsnyrHBNz#!D-D', check=False \
            ).stripBadAndGaps(), 'UCAGWSNYRHBNDD')
        self.assertEqual(self.RNA('@#^*($@!#&()!@QZX', check=False \
            ).stripBadAndGaps(), '')
        self.assertEqual(self.RNA('aaa ggg ---!ccc', check=False \
            ).stripBadAndGaps(), 'AAAGGGCCC')

    def test_shuffle(self):
        """Sequence shuffle should return new random sequence w/ same monomers"""
        r = self.RNA('UUUUCCCCAAAAGGGG')
        s = r.shuffle()
        self.assertNotEqual(r, s)
        self.assertEqualItems(r, s)

    def test_complement(self):
        """Sequence complement should correctly complement sequence"""
        self.assertEqual(self.RNA('UAUCG-NR').complement(), 'AUAGC-NY')
        self.assertEqual(self.DNA('TATCG-NR').complement(), 'ATAGC-NY')
        self.assertEqual(self.DNA('').complement(), '')
        self.assertRaises(TypeError, self.PROT('ACD').complement)

    def test_rc(self):
        """Sequence rc should correctly reverse-complement sequence"""
        #no longer preserves case!
        self.assertEqual(self.RNA('UauCG-NR').rc(), 'YN-CGAUA')
        self.assertEqual(self.DNA('TatCG-NR').rc(), 'YN-CGATA')
        self.assertEqual(self.RNA('').rc(), '')
        self.assertEqual(self.RNA('A').rc(), 'U')
        self.assertRaises(TypeError, self.PROT('ACD').rc)

    def test_contains(self):
        """Sequence contains should return correct result"""
        r = self.RNA('UCA')
        assert 'U' in r
        assert 'CA' in r
        assert 'X' not in r
        assert 'G' not in r

    def test_iter(self):
       """Sequence iter should iterate over sequence"""
       p = self.PROT('QWE')
       self.assertEqual(list(p), ['Q','W','E'])

    def test_isGapped(self):
        """Sequence isGapped should return True if gaps in seq"""
        assert not self.RNA('').isGapped()
        assert not self.RNA('ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN').isGapped()
        assert self.RNA('-').isGapped()
        assert self.PROT('--').isGapped()
        assert self.RNA('CAGUCGUACGUCAGUACGUacucauacgac-caguACUG').isGapped()
        assert self.RNA('CA--CGUAUGCA-----g').isGapped()
        assert self.RNA('CAGU-').isGapped()

    def test_isGap(self):
        """Sequence isGap should return True if char is a valid gap char"""
        r = self.RNA('ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN')
        for char in 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOASDFGHJKLZXCVBNM':
            assert not r.isGap(char)
        assert r.isGap('-')
        #only works on a single literal that's a gap, not on a sequence.
        #possibly, this behavior should change?
        assert not r.isGap('---')
        #check behaviour on self
        assert not self.RNA('CGAUACGUACGACU').isGap()
        assert not self.RNA('---CGAUA----CGUACG---ACU---').isGap()
        assert self.RNA('').isGap()
        assert self.RNA('----------').isGap()


    def test_isDegenerate(self):
        """Sequence isDegenerate should return True if degen symbol in seq"""
        assert not self.RNA('').isDegenerate()
        assert not self.RNA('UACGCUACAUGuacgucaguGCUAGCUA---ACGUCAG').isDegenerate()
        assert self.RNA('N').isDegenerate()
        assert self.RNA('R').isDegenerate()
        assert self.RNA('y').isDegenerate()
        assert self.RNA('GCAUguagcucgUCAGUCAGUACgUgcasCUAG').isDegenerate()
        assert self.RNA('ACGYAUGCUGYWWNMNuwbycwuybcwbwub').isDegenerate()

    def test_isStrict(self):
        """Sequence isStrict should return True if all symbols in Monomers"""
        assert self.RNA('').isStrict()
        assert self.PROT('A').isStrict()
        assert self.RNA('UAGCACUgcaugcauGCAUGACuacguACAUG').isStrict()
        assert not self.RNA('CAGUCGAUCA-cgaucagUCGAUGAC').isStrict()

    def test_firstGap(self):
        """Sequence firstGap should return index of first gap symbol, or None"""
        self.assertEqual(self.RNA('').firstGap(), None)
        self.assertEqual(self.RNA('a').firstGap(), None)
        self.assertEqual(self.RNA('uhacucHuhacUUhacan').firstGap(), None)
        self.assertEqual(self.RNA('-abc').firstGap(), 0)
        self.assertEqual(self.RNA('b-ac').firstGap(), 1)
        self.assertEqual(self.RNA('abcd-').firstGap(), 4)

    def test_firstDegenerate(self):
        """Sequence firstDegenerate should return index of first degen symbol"""
        self.assertEqual(self.RNA('').firstDegenerate(), None)
        self.assertEqual(self.RNA('a').firstDegenerate(), None)
        self.assertEqual(self.RNA('UCGACA--CU-gacucaguacgua'
            ).firstDegenerate(), None)
        self.assertEqual(self.RNA('nCAGU').firstDegenerate(), 0)
        self.assertEqual(self.RNA('CUGguagvAUG').firstDegenerate(), 7)
        self.assertEqual(self.RNA('ACUGCUAacgud').firstDegenerate(), 11)

    def test_firstNonStrict(self):
        """Sequence firstNonStrict should return index of first non-strict symbol"""
        self.assertEqual(self.RNA('').firstNonStrict(), None)
        self.assertEqual(self.RNA('A').firstNonStrict(), None)
        self.assertEqual(self.RNA('ACGUACGUcgaucagu').firstNonStrict(), None)
        self.assertEqual(self.RNA('N').firstNonStrict(), 0)
        self.assertEqual(self.RNA('-').firstNonStrict(), 0)
        self.assertEqual(self.RNA('ACGUcgAUGUGCAUcagu-').firstNonStrict(),18)

    def test_disambiguate(self):
        """Sequence disambiguate should remove degenerate bases"""
        self.assertEqual(self.RNA('').disambiguate(), '')
        self.assertEqual(self.RNA('AGCUGAUGUA--CAGU').disambiguate(),
            'AGCUGAUGUA--CAGU')
        self.assertEqual(self.RNA('AUn-yrs-wkmCGwmrNMWRKY').disambiguate(
            'strip'), 'AU--CG')
        s = self.RNA('AUn-yrs-wkmCGwmrNMWRKY')
        t = s.disambiguate('random')
        u = s.disambiguate('random')
        for i, j in zip(str(s), str(t)):
            if i in s.MolType.Degenerates:
                assert j in s.MolType.Degenerates[i]
            else:
                assert i == j
        self.assertNotEqual(t, u)
        self.assertEqual(len(s), len(t))

    def test_degap(self):
        """Sequence degap should remove all gaps from sequence"""
        #doesn't preserve case
        self.assertEqual(self.RNA('').degap(), '')
        self.assertEqual(self.RNA('GUCAGUCgcaugcnvuncdks').degap(), 
            'GUCAGUCGCAUGCNVUNCDKS')
        self.assertEqual(self.RNA('----------------').degap(), '')
        self.assertEqual(self.RNA('gcuauacg-').degap(), 'GCUAUACG')
        self.assertEqual(self.RNA('-CUAGUCA').degap(), 'CUAGUCA')
        self.assertEqual(self.RNA('---a---c---u----g---').degap(), 'ACUG')
        self.assertEqual(self.RNA('?a-').degap(), 'A')
        
    def test_gapList(self):
        """Sequence gapList should return correct gap positions"""
        self.assertEqual(self.RNA('').gapList(), [])
        self.assertEqual(self.RNA('ACUGUCAGUACGHSDKCUCDNNS').gapList(),[])
        self.assertEqual(self.RNA('GUACGUACAKDC-SDHDSK').gapList(),[12])
        self.assertEqual(self.RNA('-DSHUHDS').gapList(), [0])
        self.assertEqual(self.RNA('UACHASADS-').gapList(), [9])
        self.assertEqual(self.RNA('---CGAUgCAU---ACGHc---ACGUCAGU---'
            ).gapList(), [0,1,2,11,12,13,19,20,21,30,31,32])

    def test_gapVector(self):
        """Sequence gapVector should return correct gap positions"""
        g = lambda x: self.RNA(x).gapVector()
        self.assertEqual(g(''), [])
        self.assertEqual(g('ACUGUCAGUACGHCSDKCCUCCDNCNS'), [False]*27)
        self.assertEqual(g('GUACGUAACAKADC-SDAHADSAK'), 
         map(bool, map(int,'000000000000001000000000')))
        self.assertEqual(g('-DSHSUHDSS'), 
         map(bool, map(int,'1000000000')))
        self.assertEqual(g('UACHASCAGDS-'), 
         map(bool, map(int,'000000000001')))
        self.assertEqual(g('---CGAUgCAU---ACGHc---ACGUCAGU--?'), \
         map(bool, map(int,'111000000001110000011100000000111')))

    def test_gapMaps(self):
        """Sequence gapMaps should return dicts mapping gapped/ungapped pos"""
        empty = ''
        no_gaps = 'aaa'
        all_gaps = '---'
        start_gaps = '--abc'
        end_gaps = 'ab---'
        mid_gaps = '--a--b-cd---'
        gm = lambda x: self.RNA(x).gapMaps()
        self.assertEqual(gm(empty), ({},{}))
        self.assertEqual(gm(no_gaps), ({0:0,1:1,2:2}, {0:0,1:1,2:2}))
        self.assertEqual(gm(all_gaps), ({},{}))
        self.assertEqual(gm(start_gaps), ({0:2,1:3,2:4},{2:0,3:1,4:2}))
        self.assertEqual(gm(end_gaps), ({0:0,1:1},{0:0,1:1}))
        self.assertEqual(gm(mid_gaps), ({0:2,1:5,2:7,3:8},{2:0,5:1,7:2,8:3}))
     
    def test_countGaps(self):
        """Sequence countGaps should return correct gap count"""
        self.assertEqual(self.RNA('').countGaps(), 0)
        self.assertEqual(self.RNA('ACUGUCAGUACGHSDKCUCDNNS').countGaps(),
            0)
        self.assertEqual(self.RNA('GUACGUACAKDC-SDHDSK').countGaps(), 1)
        self.assertEqual(self.RNA('-DSHUHDS').countGaps(), 1)
        self.assertEqual(self.RNA('UACHASADS-').countGaps(), 1)
        self.assertEqual(self.RNA('---CGAUgCAU---ACGHc---ACGUCAGU---'
            ).countGaps(), 12)

    def test_countDegenerate(self):
        """Sequence countDegenerate should return correct degen base count"""
        self.assertEqual(self.RNA('').countDegenerate(), 0)
        self.assertEqual(self.RNA('GACUGCAUGCAUCGUACGUCAGUACCGA'
            ).countDegenerate(), 0)
        self.assertEqual(self.RNA('N').countDegenerate(), 1)
        self.assertEqual(self.PROT('N').countDegenerate(), 0)
        self.assertEqual(self.RNA('NRY').countDegenerate(), 3)
        self.assertEqual(self.RNA('ACGUAVCUAGCAUNUCAGUCAGyUACGUCAGS'
            ).countDegenerate(), 4)

    def test_possibilites(self):
        """Sequence possibilities should return correct # possible sequences"""
        self.assertEqual(self.RNA('').possibilities(), 1)
        self.assertEqual(self.RNA('ACGUgcaucagUCGuGCAU').possibilities(), 1)
        self.assertEqual(self.RNA('N').possibilities(), 4)
        self.assertEqual(self.RNA('R').possibilities(), 2)
        self.assertEqual(self.RNA('H').possibilities(), 3)
        self.assertEqual(self.RNA('nRh').possibilities(), 24)
        self.assertEqual(self.RNA('AUGCnGUCAg-aurGauc--gauhcgauacgws'
            ).possibilities(), 96)
       
    def test_MW(self):
        """Sequence MW should return correct molecular weight"""
        self.assertEqual(self.PROT('').MW(), 0)
        self.assertEqual(self.RNA('').MW(), 0)
        self.assertFloatEqual(self.PROT('A').MW(), 107.09)
        self.assertFloatEqual(self.RNA('A').MW(), 375.17)
        self.assertFloatEqual(self.PROT('AAA').MW(), 285.27)
        self.assertFloatEqual(self.RNA('AAA').MW(), 1001.59)
        self.assertFloatEqual(self.RNA('AAACCCA').MW(), 2182.37)

    def test_canMatch(self):
        """Sequence canMatch should return True if all positions can match"""
        assert self.RNA('').canMatch('')
        assert self.RNA('UCAG').canMatch('UCAG')
        assert not self.RNA('UCAG').canMatch('ucag')
        assert self.RNA('UCAG').canMatch('NNNN')
        assert self.RNA('NNNN').canMatch('UCAG')
        assert self.RNA('NNNN').canMatch('NNNN')
        assert not self.RNA('N').canMatch('x')
        assert not self.RNA('N').canMatch('-')
        assert self.RNA('UCAG').canMatch('YYRR')
        assert self.RNA('UCAG').canMatch('KMWS')

    def test_canMismatch(self):
        """Sequence canMismatch should return True on any possible mismatch"""
        assert not self.RNA('').canMismatch('')
        assert self.RNA('N').canMismatch('N')
        assert self.RNA('R').canMismatch('R')
        assert self.RNA('N').canMismatch('r')
        assert self.RNA('CGUACGCAN').canMismatch('CGUACGCAN')
        assert self.RNA('U').canMismatch('C')
        assert self.RNA('UUU').canMismatch('UUC')
        assert self.RNA('UUU').canMismatch('UUY')
        assert not self.RNA('UUU').canMismatch('UUU')
        assert not self.RNA('UCAG').canMismatch('UCAG')
        assert not self.RNA('U--').canMismatch('U--')

    def test_mustMatch(self):
        """Sequence mustMatch should return True when no possible mismatches"""
        assert self.RNA('').mustMatch('')
        assert not self.RNA('N').mustMatch('N')
        assert not self.RNA('R').mustMatch('R')
        assert not self.RNA('N').mustMatch('r')
        assert not self.RNA('CGUACGCAN').mustMatch('CGUACGCAN')
        assert not self.RNA('U').mustMatch('C')
        assert not self.RNA('UUU').mustMatch('UUC')
        assert not self.RNA('UUU').mustMatch('UUY')
        assert self.RNA('UU-').mustMatch('UU-')
        assert self.RNA('UCAG').mustMatch('UCAG')

    def test_canPair(self):
        """Sequence canPair should return True if all positions can pair"""
        assert self.RNA('').canPair('')
        assert not self.RNA('UCAG').canPair('UCAG')
        assert self.RNA('UCAG').canPair('CUGA')
        assert not self.RNA('UCAG').canPair('cuga')
        assert self.RNA('UCAG').canPair('NNNN')
        assert self.RNA('NNNN').canPair('UCAG')
        assert self.RNA('NNNN').canPair('NNNN')
        assert not self.RNA('N').canPair('x')
        assert not self.RNA('N').canPair('-')
        assert self.RNA('-').canPair('-')
        assert self.RNA('UCAGU').canPair('KYYRR')
        assert self.RNA('UCAG').canPair('KKRS')
        assert self.RNA('U').canPair('G')

        assert not self.DNA('T').canPair('G')

    def test_canMispair(self):
        """Sequence canMispair should return True on any possible mispair"""
        assert not self.RNA('').canMispair('')
        assert self.RNA('N').canMispair('N')
        assert self.RNA('R').canMispair('Y')
        assert self.RNA('N').canMispair('r')
        assert self.RNA('CGUACGCAN').canMispair('NUHCHUACH')
        assert self.RNA('U').canMispair('C')
        assert self.RNA('U').canMispair('R')
        assert self.RNA('UUU').canMispair('AAR')
        assert self.RNA('UUU').canMispair('GAG')
        assert not self.RNA('UUU').canMispair('AAA')
        assert not self.RNA('UCAG').canMispair('CUGA')
        assert self.RNA('U--').canMispair('--U')

        assert self.DNA('TCCAAAGRYY').canMispair('RRYCTTTGGA')

    def test_mustPair(self):
        """Sequence mustPair should return True when no possible mispairs"""
        assert self.RNA('').mustPair('')
        assert not self.RNA('N').mustPair('N')
        assert not self.RNA('R').mustPair('Y')
        assert not self.RNA('A').mustPair('A')
        assert not self.RNA('CGUACGCAN').mustPair('NUGCGUACG')
        assert not self.RNA('U').mustPair('C')
        assert not self.RNA('UUU').mustPair('AAR')
        assert not self.RNA('UUU').mustPair('RAA')
        assert not self.RNA('UU-').mustPair('-AA')
        assert self.RNA('UCAG').mustPair('CUGA')

        assert self.DNA('TCCAGGG').mustPair('CCCTGGA')
        assert self.DNA('tccaggg').mustPair(self.DNA('ccctgga'))
        assert not self.DNA('TCCAGGG').mustPair('NCCTGGA')

    def test_diff(self):
        """Sequence diff should count 1 for each difference between sequences"""
        self.assertEqual(self.RNA('UGCUGCUC').diff(''), 0)
        self.assertEqual(self.RNA('UGCUGCUC').diff('U'), 0)
        self.assertEqual(self.RNA('UGCUGCUC').diff('UCCCCCUC'), 3)
        #case-sensitive!
        self.assertEqual(self.RNA('AAAAA').diff('CCCCC'), 5)
        #raises TypeError if other not iterable
        self.assertRaises(TypeError, self.RNA('AAAAA').diff, 5)
        
    def test_distance(self):
        """Sequence distance should calculate correctly based on function"""
        def f(a, b):
            if a == b:
                return 0
            if (a in 'UC' and b in 'UC') or (a in 'AG' and b in 'AG'):
                return 1
            else:
                return 10
        #uses identity function by default
        self.assertEqual(self.RNA('UGCUGCUC').distance(''), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('U'), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('UCCCCCUC'), 3)
        #case-sensitive!
        self.assertEqual(self.RNA('AAAAA').distance('CCCCC'), 5)
        #should use function if supplied
        self.assertEqual(self.RNA('UGCUGCUC').distance('', f), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('U', f), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('C', f), 1)
        self.assertEqual(self.RNA('UGCUGCUC').distance('G', f), 10)
        self.assertEqual(self.RNA('UGCUGCUC').distance('UCCCCCUC', f), 21)
        #case-sensitive!
        self.assertEqual(self.RNA('AAAAA').distance('CCCCC', f), 50)

    def test_matrixDistance(self):
        """Sequence matrixDistance should look up distances from a matrix"""
        #note that the score matrix must contain 'diagonal' elements m[i][i] 
        #to avoid failure when the sequences match.
        m = {'U':{'U':0, 'C':1, 'A':5}, 'C':{'C':0, 'A':2,'G':4}}
        self.assertEqual(self.RNA('UUUCCC').matrixDistance('UCACGG', m), 14)
        self.assertEqual(self.RNA('UUUCCC').matrixDistance('', m), 0)
        self.assertEqual(self.RNA('UUU').matrixDistance('CAC', m), 7)
        self.assertRaises(KeyError, self.RNA('UUU').matrixDistance, 'CAG', m)

    def test_fracSame(self):
        """Sequence fracSame should return similarity between sequences"""
        s1 = self.RNA('ACGU')
        s2 = self.RNA('AACG')
        s3 = self.RNA('GG')
        s4 = self.RNA('A')
        e = self.RNA('')
        self.assertEqual(s1.fracSame(e), 0)
        self.assertEqual(s1.fracSame(s2), 0.25)
        self.assertEqual(s1.fracSame(s3), 0)
        self.assertEqual(s1.fracSame(s4), 1.0)  #note truncation

    def test_fracDiff(self):
        """Sequence fracDiff should return difference between sequences"""
        s1 = self.RNA('ACGU')
        s2 = self.RNA('AACG')
        s3 = self.RNA('GG')
        s4 = self.RNA('A')
        e = self.RNA('')
        self.assertEqual(s1.fracDiff(e), 0)
        self.assertEqual(s1.fracDiff(s2), 0.75)
        self.assertEqual(s1.fracDiff(s3), 1)
        self.assertEqual(s1.fracDiff(s4), 0)  #note truncation

    def test_fracSameGaps(self):
        """Sequence fracSameGaps should return similarity in gap positions"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('GGGG')
        s3 = self.RNA('----')
        s4 = self.RNA('A-A-')
        s5 = self.RNA('-G-G')
        s6 = self.RNA('UU--')
        s7 = self.RNA('-')
        s8 = self.RNA('GGG')
        e =  self.RNA('')
        self.assertEqual(s1.fracSameGaps(s1), 1)
        self.assertEqual(s1.fracSameGaps(s2), 1)
        self.assertEqual(s1.fracSameGaps(s3), 0)
        self.assertEqual(s1.fracSameGaps(s4), 0.5)
        self.assertEqual(s1.fracSameGaps(s5), 0.5)
        self.assertEqual(s1.fracSameGaps(s6), 0.5)
        self.assertEqual(s1.fracSameGaps(s7), 0)
        self.assertEqual(s1.fracSameGaps(e), 0)
        self.assertEqual(s3.fracSameGaps(s3), 1)
        self.assertEqual(s3.fracSameGaps(s4), 0.5)
        self.assertEqual(s3.fracSameGaps(s7), 1.0)
        self.assertEqual(e.fracSameGaps(e), 0.0)
        self.assertEqual(s4.fracSameGaps(s5), 0.0)
        self.assertEqual(s4.fracSameGaps(s6), 0.5)
        self.assertFloatEqual(s6.fracSameGaps(s8), 2/3.0)
        
    def test_fracDiffGaps(self):
        """Sequence fracDiffGaps should return difference in gap positions"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('GGGG')
        s3 = self.RNA('----')
        s4 = self.RNA('A-A-')
        s5 = self.RNA('-G-G')
        s6 = self.RNA('UU--')
        s7 = self.RNA('-')
        s8 = self.RNA('GGG')
        e =  self.RNA('')
        self.assertEqual(s1.fracDiffGaps(s1), 0)
        self.assertEqual(s1.fracDiffGaps(s2), 0)
        self.assertEqual(s1.fracDiffGaps(s3), 1)
        self.assertEqual(s1.fracDiffGaps(s4), 0.5)
        self.assertEqual(s1.fracDiffGaps(s5), 0.5)
        self.assertEqual(s1.fracDiffGaps(s6), 0.5)
        self.assertEqual(s1.fracDiffGaps(s7), 1)
        self.assertEqual(s1.fracDiffGaps(e), 0)
        self.assertEqual(s3.fracDiffGaps(s3), 0)
        self.assertEqual(s3.fracDiffGaps(s4), 0.5)
        self.assertEqual(s3.fracDiffGaps(s7), 0.0)
        self.assertEqual(e.fracDiffGaps(e), 0.0)
        self.assertEqual(s4.fracDiffGaps(s5), 1.0)
        self.assertEqual(s4.fracDiffGaps(s6), 0.5)
        self.assertFloatEqual(s6.fracDiffGaps(s8), 1/3.0)

    def test_fracSameNonGaps(self):
        """Sequence fracSameNonGaps should return similarities at non-gaps"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('AGGG')
        s3 = self.RNA('GGGG')
        s4 = self.RNA('AG--GA-G')
        s5 = self.RNA('CU--CU-C')
        s6 = self.RNA('AC--GC-G')
        s7 = self.RNA('--------')
        s8 = self.RNA('AAAA----')
        s9 = self.RNA('A-GG-A-C')
        e =  self.RNA('')

        test = lambda x, y, z: self.assertFloatEqual(x.fracSameNonGaps(y), z)
        test(s1, s2, 0.25)
        test(s1, s3, 0)
        test(s2, s3, 0.75)
        test(s1, s4, 0.5)
        test(s4, s5, 0)
        test(s4, s6, 0.6)
        test(s4, s7, 0)
        test(s4, s8, 0.5)
        test(s4, s9, 2/3.0)
        test(e, s4, 0)

    def test_fracDiffNonGaps(self):
        """Sequence fracDiffNonGaps should return differences at non-gaps"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('AGGG')
        s3 = self.RNA('GGGG')
        s4 = self.RNA('AG--GA-G')
        s5 = self.RNA('CU--CU-C')
        s6 = self.RNA('AC--GC-G')
        s7 = self.RNA('--------')
        s8 = self.RNA('AAAA----')
        s9 = self.RNA('A-GG-A-C')
        e =  self.RNA('')

        test = lambda x, y, z: self.assertFloatEqual(x.fracDiffNonGaps(y), z)
        test(s1, s2, 0.75)
        test(s1, s3, 1)
        test(s2, s3, 0.25)
        test(s1, s4, 0.5)
        test(s4, s5, 1)
        test(s4, s6, 0.4)
        test(s4, s7, 0)
        test(s4, s8, 0.5)
        test(s4, s9, 1/3.0)
        test(e, s4, 0)

    def test_fracSimilar(self):
        """Sequence fracSimilar should return the fraction similarity"""
        transitions = dict.fromkeys([ \
            ('A','A'), ('A','G'), ('G','A'), ('G','G'),
            ('U','U'), ('U','C'), ('C','U'), ('C','C')])
        
        s1 = self.RNA('UCAGGCAA')
        s2 = self.RNA('CCAAAUGC')
        s3 = self.RNA('GGGGGGGG')
        e =  self.RNA('')

        test = lambda x, y, z: self.assertFloatEqual( \
            x.fracSimilar(y, transitions), z)

        test(e, e, 0)
        test(s1, e, 0)
        test(s1, s1, 1)
        test(s1, s2, 7.0/8)
        test(s1, s3, 5.0/8)
        test(s2,s3, 4.0/8)

    def test_withTerminiUnknown(self):
        """withTerminiUnknown should reset termini to unknown char"""
        s1 = self.RNA('-?--AC--?-')
        s2 = self.RNA('AC')
        self.assertEqual(s1.withTerminiUnknown(), '????AC????')
        self.assertEqual(s2.withTerminiUnknown(), 'AC')
    
    def test_consistent_gap_degen_handling(self):
        """gap degen character should be treated consistently"""
        # the degen character '?' can be a gap, so when we strip either gaps or
        # degen characters it should be gone too
        raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
        raw_ungapped = re.sub("[-?]", "", raw_seq)
        raw_no_ambigs = re.sub("[N?]+", "", raw_seq)
        dna = self.DNA(raw_seq)
        self.assertEqual(dna.degap(), raw_ungapped)
        self.assertEqual(dna.stripDegenerate(), raw_no_ambigs)
        self.assertEqual(dna.stripBadAndGaps(), raw_ungapped)
    

class SequenceSubclassTests(TestCase):
    """Only one general set of tests, since the subclasses are very thin."""
    def test_DnaSequence(self):
        """DnaSequence should behave as expected"""
        x = DnaSequence('tcag')
        #note: no longer preserves case
        self.assertEqual(x, 'TCAG')
        
        x = DnaSequence('aaa') + DnaSequence('ccc')
        #note: doesn't preserve case
        self.assertEqual(x, 'AAACCC')
        assert x.MolType is DNA
        self.assertRaises(AlphabetError, x.__add__, 'z')
        self.assertEqual(DnaSequence('TTTAc').rc(), 'GTAAA')


class ModelSequenceTests(object):
    """Base class for tests of specific ModelSequence objects."""
    SequenceClass = None   #override in derived classes
    
    def test_toFasta(self):
        """Sequence toFasta() should return Fasta-format string"""
        even = 'TCAGAT'
        odd = even + 'AAA'
        even_dna = self.SequenceClass(even, Name='even')
        odd_dna = self.SequenceClass(odd, Name='odd')
        self.assertEqual(even_dna.toFasta(), '>even\nTCAGAT')
        #set line wrap to small number so we can test that it works
        even_dna.LineWrap = 2
        self.assertEqual(even_dna.toFasta(), '>even\nTC\nAG\nAT')
        odd_dna.LineWrap = 2
        self.assertEqual(odd_dna.toFasta(), '>odd\nTC\nAG\nAT\nAA\nA')
        #check that changing the linewrap again works
        even_dna.LineWrap  = 4
        self.assertEqual(even_dna.toFasta(), '>even\nTCAG\nAT')

    def test_toPhylip(self):
        """Sequence toPhylip() should return one-line phylip string"""
        s = self.SequenceClass('ACG', Name='xyz')
        self.assertEqual(s.toPhylip(), 'xyz'+' '*27+'ACG')
        

class DnaSequenceTests(ModelSequenceTests, TestCase):
    
    class SequenceClass(ModelNucleicAcidSequence):
        Alphabet = DNA.Alphabets.Base

    def test_init(self):
        """Sequence should do round-trip from string"""
        orig = ''
        r = self.SequenceClass(orig)
        self.assertEqual(str(r), orig)
        
        orig = 'TCAGGA'
        r = self.SequenceClass(orig)
        self.assertEqual(r._data, array([0,1,2,3,3,2]))
        self.assertEqual(str(r), orig)

    def test_toKwords(self):
        """Sequence toKwords should give expected counts"""
        orig = 'ATCCCTAGC'
        r = self.SequenceClass(orig)
        #if we use k = 1, should just get the characters
        w = r.toKwords(1)
        self.assertEqual(w, r._data)
        w = r.toKwords(1, overlapping=False)
        self.assertEqual(w, r._data)
        #if we use k = 2, should get overlapping or nonoverlapping k-words
        w = r.toKwords(2)
        self.assertEqual(w, array([8,1,5,5,4,2,11,13]))
        w = r.toKwords(2, overlapping=False)
        self.assertEqual(w, array([8,5,4,11]))
        #check a case with k = 3, i.e. codons
        w = r.toKwords(3, overlapping=False)
        self.assertEqual(w, array([33,20,45]))
    

class CodonSequenceTests(SequenceTests, TestCase):
    class SequenceClass(ModelCodonSequence):
        Alphabet = DNA.Alphabets.Base.Triples
    def test_init(self):
        """Sequence should do round-trip from string"""
        orig = ''
        r = self.SequenceClass(orig)
        self.assertEqual(str(r), orig)
        
        orig = 'TCAGGA'
        r = self.SequenceClass(orig)
        self.assertEqual(r._data, array([6,62]))
        self.assertEqual(str(r), orig)

    def test_toKwords(self):
        """Sequence toKwords should give expected counts"""
        orig = 'ATCCCTAGC'
        r = self.SequenceClass(orig)
        #if we use k = 1, should just get the characters
        w = r.toKwords(1)
        self.assertEqual(w, r._data)
        w = r.toKwords(1, overlapping=False)
        self.assertEqual(w, r._data)
        #if we use k = 2, should get overlapping or nonoverlapping k-words
        w = r.toKwords(2)
        self.assertEqual(w, array([2132,1325]))
        w = r.toKwords(2, overlapping=False)
        self.assertEqual(w, array([2132]))
  
class DnaSequenceGapTests(TestCase):
    """Tests of gapped DNA sequences."""
    class SequenceClass(ModelNucleicAcidSequence):
        Alphabet = DNA.Alphabets.Gapped
        Gap = '-'
    def test_init(self):
        """Gapped sequence should init ok"""
        orig = 'TC---'
        seq = self.SequenceClass(orig)
        self.assertEqual(str(seq), orig)
        
    def test_gaps(self):
        """Gapped sequence gaps() should return correct array"""
        sc = self.SequenceClass
        self.assertEqual(sc('TC').gaps(), array([0,0]))
        self.assertEqual(sc('T-').gaps(), array([0,1]))

    def test_degap(self):
        """Gapped sequence degap() should return correct array"""
        sc = self.SequenceClass
        self.assertEqual(sc('T-').degap(), sc('T'))
  
    def test_nongaps(self):
        """Gapped sequence nongaps() should return correct array"""
        sc = self.SequenceClass
        self.assertEqual(sc('TC').nongaps(), array([1,1]))
        self.assertEqual(sc('T-').nongaps(), array([1,0]))
  
    def test_regap(self):
        """Gapped sequence regap() should return correct sequence"""
        sc = self.SequenceClass
        self.assertEqual(str(sc('TC').regap(sc('A---A-'))), 'T---C-')

class SequenceIntegrationTests(TestCase):
    """Should be able to convert regular to model sequences, and back"""
    
    def test_regular_to_model(self):
        """Regular sequence should convert to model sequence"""
        r = RNA.Sequence('AAA', Name='x')
        s = RNA.ModelSeq(r)
        self.assertEqual(str(s), 'AAA')
        self.assertEqual(s.MolType, RNA)
        self.assertEqual(s.Name, 'x')
    
    def test_model_to_regular(self):
        """Model sequence should convert to regular sequence"""
        r = RNA.ModelSeq('AAA', Name='x')
        s = RNA.Sequence(r)
        self.assertEqual(str(s), 'AAA')
        self.assertEqual(s.MolType, RNA)
        self.assertEqual(s.Name, 'x')
    
    def test_regular_to_regular(self):
        """Regular sequence should convert to regular sequence"""
        r = RNA.Sequence('AAA', Name='x')
        s = RNA.Sequence(r)
        self.assertEqual(str(s), 'AAA')
        self.assertEqual(s.MolType, RNA)
        self.assertEqual(s.Name, 'x')
    
    def test_model_to_model(self):
        """Model sequence should convert to model sequence"""
        r = RNA.ModelSeq('AAA', Name='x')
        s = RNA.ModelSeq(r)
        self.assertEqual(str(s), 'AAA')
        self.assertEqual(s.MolType, RNA)
        self.assertEqual(s.Name, 'x')

    def test_ModelDnaCodonSequence(self):
        """ModelDnaCodonSequence should behave as expected"""
        d = ModelDnaCodonSequence('UUUCGU')
        self.assertEqual(str(d), 'TTTCGT')
        self.assertEqual(d._data, array([0,28]))
        self.assertEqual(str(d.toRna()), 'UUUCGU')
        self.assertEqual(str(d.toDna()), 'TTTCGT')
    
    def test_ModelRnaCodonSequence(self):
        """ModelRnaCodonSequence should behave as expected"""
        r = ModelRnaCodonSequence('UUUCGU')
        self.assertEqual(str(r), 'UUUCGU')
        self.assertEqual(r._data, array([0,28]))
        self.assertEqual(str(r.toRna()), 'UUUCGU')
        self.assertEqual(str(r.toDna()), 'TTTCGT')

class ModelSequenceTests(SequenceTests):
    """Tests of the ModelSequence class's inheritance of SequenceI."""
    SEQ = ModelSequence
    RNA = ModelRnaSequence
    DNA = ModelDnaSequence
    PROT = ModelProteinSequence

    def test_distance_indices(self):
        """ModelSequence distance should work with function of indices"""
        s1 = self.RNA('AUGC')
        s2 = self.RNA('AAGC')
        def f(x,y):
            if x == 2 or y == 2:
                return 10
            return 0
        self.assertEqual(s1.distance(s2, f, use_indices=True), 20)
    
    def test_stripBad(self):
        """Sequence stripBad should remove any non-base, non-gap chars"""
        #have to turn off check to get bad data in; no longer preserves case
        r = self.RNA('UCAGRYU')
        r._data[0] = 31
        r._data[2] = 55
        self.assertEqual(r.stripBad(), 'CGRYU')

    def test_stripBadAndGaps(self):
        """Sequence stripBadAndGaps should remove gaps and bad chars"""
        #have to turn off check to get bad data in; no longer preserves case
        r = self.RNA('ACG--GRN?')
        self.assertEqual(r.stripBadAndGaps(), 'ACGGRN')
        r._data[0] = 99
        self.assertEqual(r.stripBadAndGaps(), 'CGGRN')

    def test_gapArray(self):
        """Sequence gapArray should return array of gaps"""
        r = self.RNA('-?A-?NRY-')
        v = r.gapArray()
        self.assertEqual(v, array([1,1,0,1,1,0,0,0,1]))
        r = self.RNA('AC')
        v = r.gapArray()
        self.assertEqual(v, array([0,0]))
        r = self.RNA('-?')
        v = r.gapArray()
        self.assertEqual(v, array([1,1]))

    def test_gapIndices(self):
        """Sequence gapIndices should return positions of gaps"""
        r = self.RNA('-?A-?NRY-')
        v = r.gapIndices()
        self.assertEqual(v, array([0,1,3,4,8]))
        r = self.RNA('AC')
        v = r.gapIndices()
        self.assertEqual(v, array([])) #note: always returns array
        r = self.RNA('-?')
        v = r.gapIndices()
        self.assertEqual(v, array([0,1]))

 
#run if called from command-line
if __name__ == "__main__":
    main()
