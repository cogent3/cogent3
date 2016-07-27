#!/usr/bin/env python
"""Unit tests for Sequence class and its subclasses.
"""

from cogent3.core.sequence import Sequence, RnaSequence, DnaSequence, \
    ProteinSequence, ModelSequenceBase, \
    ModelSequence, ModelNucleicAcidSequence, ModelRnaSequence, \
    ModelDnaSequence, ModelProteinSequence, ModelCodonSequence, \
    ModelDnaCodonSequence, ModelRnaCodonSequence
from cogent3.core.moltype import RNA, DNA, PROTEIN, ASCII, BYTES, AlphabetError
from cogent3.util.unit_test import TestCase, main

import re
from pickle import dumps
from numpy import array

__author__ = "Rob Knight, Gavin Huttley and Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Peter Maxwell",
               "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
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
        # NOTE: ModelSequences can't be initialized empty because it screws up
        # the dimensions of the array, and not worth special-casing.
        s = self.SEQ()
        self.assertEqual(s, '')
        assert s.MolType in (ASCII, BYTES)

        r = self.RNA()
        assert r.MolType is RNA

    def test_init_data(self):
        """Sequence init with data should set data in correct location"""
        r = self.RNA('ucagg')
        # no longer preserves case
        self.assertEqual(r, 'UCAGG')

    def test_init_other_seq(self):
        """Sequence init with other seq should preserve name and info."""
        r = self.RNA('UCAGG', Name='x', Info={'z': 3})
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

    def test_to_dna(self):
        """Returns copy of self as DNA."""
        r = self.RNA('TCA')
        self.assertEqual(str(r), 'UCA')
        self.assertEqual(str(r.to_dna()), 'TCA')

    def test_to_rna(self):
        """Returns copy of self as RNA."""
        r = self.DNA('UCA')
        self.assertEqual(str(r), 'TCA')
        self.assertEqual(str(r.to_rna()), 'UCA')

    def test_to_fasta(self):
        """Sequence to_fasta() should return Fasta-format string"""
        even = 'TCAGAT'
        odd = even + 'AAA'
        even_dna = self.SEQ(even, Name='even')
        odd_dna = self.SEQ(odd, Name='odd')
        self.assertEqual(even_dna.to_fasta(), '>even\nTCAGAT')
        # set line wrap to small number so we can test that it works
        even_dna.LineWrap = 2
        self.assertEqual(even_dna.to_fasta(), '>even\nTC\nAG\nAT')
        odd_dna.LineWrap = 2
        self.assertEqual(odd_dna.to_fasta(), '>odd\nTC\nAG\nAT\nAA\nA')
        # check that changing the linewrap again works
        even_dna.LineWrap = 4
        self.assertEqual(even_dna.to_fasta(), '>even\nTCAG\nAT')

    def test_serialize(self):
        """Sequence should be serializable"""
        r = self.RNA('ugagg')
        assert dumps(r)

    def test_strip_degenerate(self):
        """Sequence strip_degenerate should remove any degenerate bases"""
        self.assertEqual(self.RNA('UCAG-').strip_degenerate(), 'UCAG-')
        self.assertEqual(self.RNA('NRYSW').strip_degenerate(), '')
        self.assertEqual(self.RNA('USNG').strip_degenerate(), 'UG')

    def test_strip_bad(self):
        """Sequence strip_bad should remove any non-base, non-gap chars"""
        # have to turn off check to get bad data in; no longer preserves case
        self.assertEqual(self.RNA('UCxxxAGwsnyrHBNzzzD-D', check=False
                                  ).strip_bad(), 'UCAGWSNYRHBND-D')
        self.assertEqual(self.RNA('@#^*($@!#&()!@QZX', check=False
                                  ).strip_bad(), '')
        self.assertEqual(self.RNA('aaaxggg---!ccc', check=False).strip_bad(),
                         'AAAGGG---CCC')

    def test_strip_bad_and_gaps(self):
        """Sequence strip_bad_and_gaps should remove gaps and bad chars"""
        # have to turn off check to get bad data in; no longer preserves case
        self.assertEqual(self.RNA('UxxCAGwsnyrHBNz#!D-D', check=False
                                  ).strip_bad_and_gaps(), 'UCAGWSNYRHBNDD')
        self.assertEqual(self.RNA('@#^*($@!#&()!@QZX', check=False
                                  ).strip_bad_and_gaps(), '')
        self.assertEqual(self.RNA('aaa ggg ---!ccc', check=False
                                  ).strip_bad_and_gaps(), 'AAAGGGCCC')

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
        # no longer preserves case!
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
        self.assertEqual(list(p), ['Q', 'W', 'E'])

    def test_is_gapped(self):
        """Sequence is_gapped should return True if gaps in seq"""
        assert not self.RNA('').is_gapped()
        assert not self.RNA('ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN').is_gapped()
        assert self.RNA('-').is_gapped()
        assert self.PROT('--').is_gapped()
        assert self.RNA('CAGUCGUACGUCAGUACGUacucauacgac-caguACUG').is_gapped()
        assert self.RNA('CA--CGUAUGCA-----g').is_gapped()
        assert self.RNA('CAGU-').is_gapped()

    def test_is_gap(self):
        """Sequence is_gap should return True if char is a valid gap char"""
        r = self.RNA('ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN')
        for char in 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOASDFGHJKLZXCVBNM':
            assert not r.is_gap(char)
        assert r.is_gap('-')
        # only works on a single literal that's a gap, not on a sequence.
        # possibly, this behavior should change?
        assert not r.is_gap('---')
        # check behaviour on self
        assert not self.RNA('CGAUACGUACGACU').is_gap()
        assert not self.RNA('---CGAUA----CGUACG---ACU---').is_gap()
        assert self.RNA('').is_gap()
        assert self.RNA('----------').is_gap()

    def test_is_degenerate(self):
        """Sequence is_degenerate should return True if degen symbol in seq"""
        assert not self.RNA('').is_degenerate()
        assert not self.RNA(
            'UACGCUACAUGuacgucaguGCUAGCUA---ACGUCAG').is_degenerate()
        assert self.RNA('N').is_degenerate()
        assert self.RNA('R').is_degenerate()
        assert self.RNA('y').is_degenerate()
        assert self.RNA('GCAUguagcucgUCAGUCAGUACgUgcasCUAG').is_degenerate()
        assert self.RNA('ACGYAUGCUGYWWNMNuwbycwuybcwbwub').is_degenerate()

    def test_is_strict(self):
        """Sequence is_strict should return True if all symbols in Monomers"""
        assert self.RNA('').is_strict()
        assert self.PROT('A').is_strict()
        assert self.RNA('UAGCACUgcaugcauGCAUGACuacguACAUG').is_strict()
        assert not self.RNA('CAGUCGAUCA-cgaucagUCGAUGAC').is_strict()

    def test_first_gap(self):
        """Sequence first_gap should return index of first gap symbol, or None"""
        self.assertEqual(self.RNA('').first_gap(), None)
        self.assertEqual(self.RNA('a').first_gap(), None)
        self.assertEqual(self.RNA('uhacucHuhacUUhacan').first_gap(), None)
        self.assertEqual(self.RNA('-abc').first_gap(), 0)
        self.assertEqual(self.RNA('b-ac').first_gap(), 1)
        self.assertEqual(self.RNA('abcd-').first_gap(), 4)

    def test_first_degenerate(self):
        """Sequence first_degenerate should return index of first degen symbol"""
        self.assertEqual(self.RNA('').first_degenerate(), None)
        self.assertEqual(self.RNA('a').first_degenerate(), None)
        self.assertEqual(self.RNA('UCGACA--CU-gacucaguacgua'
                                  ).first_degenerate(), None)
        self.assertEqual(self.RNA('nCAGU').first_degenerate(), 0)
        self.assertEqual(self.RNA('CUGguagvAUG').first_degenerate(), 7)
        self.assertEqual(self.RNA('ACUGCUAacgud').first_degenerate(), 11)

    def test_first_non_strict(self):
        """Sequence first_non_strict should return index of first non-strict symbol"""
        self.assertEqual(self.RNA('').first_non_strict(), None)
        self.assertEqual(self.RNA('A').first_non_strict(), None)
        self.assertEqual(self.RNA('ACGUACGUcgaucagu').first_non_strict(), None)
        self.assertEqual(self.RNA('N').first_non_strict(), 0)
        self.assertEqual(self.RNA('-').first_non_strict(), 0)
        self.assertEqual(self.RNA('ACGUcgAUGUGCAUcagu-').first_non_strict(), 18)

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
        # doesn't preserve case
        self.assertEqual(self.RNA('').degap(), '')
        self.assertEqual(self.RNA('GUCAGUCgcaugcnvuncdks').degap(),
                         'GUCAGUCGCAUGCNVUNCDKS')
        self.assertEqual(self.RNA('----------------').degap(), '')
        self.assertEqual(self.RNA('gcuauacg-').degap(), 'GCUAUACG')
        self.assertEqual(self.RNA('-CUAGUCA').degap(), 'CUAGUCA')
        self.assertEqual(self.RNA('---a---c---u----g---').degap(), 'ACUG')
        self.assertEqual(self.RNA('?a-').degap(), 'A')

    def test_gap_indices(self):
        """Sequence gap_indices should return correct gap positions"""
        self.assertEqual(self.RNA('').gap_indices(), [])
        self.assertEqual(self.RNA('ACUGUCAGUACGHSDKCUCDNNS').gap_indices(), [])
        self.assertEqual(self.RNA('GUACGUACAKDC-SDHDSK').gap_indices(), [12])
        self.assertEqual(self.RNA('-DSHUHDS').gap_indices(), [0])
        self.assertEqual(self.RNA('UACHASADS-').gap_indices(), [9])
        self.assertEqual(self.RNA('---CGAUgCAU---ACGHc---ACGUCAGU---'
                                  ).gap_indices(), [0, 1, 2, 11, 12, 13, 19, 20, 21, 30, 31, 32])

    def test_gap_vector(self):
        """Sequence gap_vector should return correct gap positions"""
        g = lambda x: self.RNA(x).gap_vector()
        self.assertEqual(g(''), [])
        self.assertEqual(g('ACUGUCAGUACGHCSDKCCUCCDNCNS'), [False] * 27)
        self.assertEqual(g('GUACGUAACAKADC-SDAHADSAK'),
                         list(map(bool, list(map(int, '000000000000001000000000')))))
        self.assertEqual(g('-DSHSUHDSS'),
                         list(map(bool, list(map(int, '1000000000')))))
        self.assertEqual(g('UACHASCAGDS-'),
                         list(map(bool, list(map(int, '000000000001')))))
        self.assertEqual(g('---CGAUgCAU---ACGHc---ACGUCAGU--?'),
                         list(map(bool, list(map(int, '111000000001110000011100000000111')))))

    def test_gap_maps(self):
        """Sequence gap_maps should return dicts mapping gapped/ungapped pos"""
        empty = ''
        no_gaps = 'aaa'
        all_gaps = '---'
        start_gaps = '--abc'
        end_gaps = 'ab---'
        mid_gaps = '--a--b-cd---'
        gm = lambda x: self.RNA(x).gap_maps()
        self.assertEqual(gm(empty), ({}, {}))
        self.assertEqual(gm(no_gaps), ({0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2}))
        self.assertEqual(gm(all_gaps), ({}, {}))
        self.assertEqual(
            gm(start_gaps), ({0: 2, 1: 3, 2: 4}, {2: 0, 3: 1, 4: 2}))
        self.assertEqual(gm(end_gaps), ({0: 0, 1: 1}, {0: 0, 1: 1}))
        self.assertEqual(
            gm(mid_gaps), ({0: 2, 1: 5, 2: 7, 3: 8}, {2: 0, 5: 1, 7: 2, 8: 3}))

    def test_count_gaps(self):
        """Sequence count_gaps should return correct gap count"""
        self.assertEqual(self.RNA('').count_gaps(), 0)
        self.assertEqual(self.RNA('ACUGUCAGUACGHSDKCUCDNNS').count_gaps(),
                         0)
        self.assertEqual(self.RNA('GUACGUACAKDC-SDHDSK').count_gaps(), 1)
        self.assertEqual(self.RNA('-DSHUHDS').count_gaps(), 1)
        self.assertEqual(self.RNA('UACHASADS-').count_gaps(), 1)
        self.assertEqual(self.RNA('---CGAUgCAU---ACGHc---ACGUCAGU---'
                                  ).count_gaps(), 12)

    def test_count_degenerate(self):
        """Sequence count_degenerate should return correct degen base count"""
        self.assertEqual(self.RNA('').count_degenerate(), 0)
        self.assertEqual(self.RNA('GACUGCAUGCAUCGUACGUCAGUACCGA'
                                  ).count_degenerate(), 0)
        self.assertEqual(self.RNA('N').count_degenerate(), 1)
        self.assertEqual(self.PROT('N').count_degenerate(), 0)
        self.assertEqual(self.RNA('NRY').count_degenerate(), 3)
        self.assertEqual(self.RNA('ACGUAVCUAGCAUNUCAGUCAGyUACGUCAGS'
                                  ).count_degenerate(), 4)

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
        self.assertEqual(self.PROT('').mw(), 0)
        self.assertEqual(self.RNA('').mw(), 0)
        self.assertFloatEqual(self.PROT('A').mw(), 89.09)
        self.assertFloatEqual(self.RNA('A').mw(), 375.17)
        self.assertFloatEqual(self.PROT('AAA').mw(), 231.27)
        self.assertFloatEqual(self.RNA('AAA').mw(), 1001.59)
        self.assertFloatEqual(self.RNA('AAACCCA').mw(), 2182.37)

    def test_can_match(self):
        """Sequence can_match should return True if all positions can match"""
        assert self.RNA('').can_match('')
        assert self.RNA('UCAG').can_match('UCAG')
        assert not self.RNA('UCAG').can_match('ucag')
        assert self.RNA('UCAG').can_match('NNNN')
        assert self.RNA('NNNN').can_match('UCAG')
        assert self.RNA('NNNN').can_match('NNNN')
        assert not self.RNA('N').can_match('x')
        assert not self.RNA('N').can_match('-')
        assert self.RNA('UCAG').can_match('YYRR')
        assert self.RNA('UCAG').can_match('KMWS')

    def test_can_mismatch(self):
        """Sequence can_mismatch should return True on any possible mismatch"""
        assert not self.RNA('').can_mismatch('')
        assert self.RNA('N').can_mismatch('N')
        assert self.RNA('R').can_mismatch('R')
        assert self.RNA('N').can_mismatch('r')
        assert self.RNA('CGUACGCAN').can_mismatch('CGUACGCAN')
        assert self.RNA('U').can_mismatch('C')
        assert self.RNA('UUU').can_mismatch('UUC')
        assert self.RNA('UUU').can_mismatch('UUY')
        assert not self.RNA('UUU').can_mismatch('UUU')
        assert not self.RNA('UCAG').can_mismatch('UCAG')
        assert not self.RNA('U--').can_mismatch('U--')

    def test_must_match(self):
        """Sequence must_match should return True when no possible mismatches"""
        assert self.RNA('').must_match('')
        assert not self.RNA('N').must_match('N')
        assert not self.RNA('R').must_match('R')
        assert not self.RNA('N').must_match('r')
        assert not self.RNA('CGUACGCAN').must_match('CGUACGCAN')
        assert not self.RNA('U').must_match('C')
        assert not self.RNA('UUU').must_match('UUC')
        assert not self.RNA('UUU').must_match('UUY')
        assert self.RNA('UU-').must_match('UU-')
        assert self.RNA('UCAG').must_match('UCAG')

    def test_can_pair(self):
        """Sequence can_pair should return True if all positions can pair"""
        assert self.RNA('').can_pair('')
        assert not self.RNA('UCAG').can_pair('UCAG')
        assert self.RNA('UCAG').can_pair('CUGA')
        assert not self.RNA('UCAG').can_pair('cuga')
        assert self.RNA('UCAG').can_pair('NNNN')
        assert self.RNA('NNNN').can_pair('UCAG')
        assert self.RNA('NNNN').can_pair('NNNN')
        assert not self.RNA('N').can_pair('x')
        assert not self.RNA('N').can_pair('-')
        assert self.RNA('-').can_pair('-')
        assert self.RNA('UCAGU').can_pair('KYYRR')
        assert self.RNA('UCAG').can_pair('KKRS')
        assert self.RNA('U').can_pair('G')

        assert not self.DNA('T').can_pair('G')

    def test_can_mispair(self):
        """Sequence can_mispair should return True on any possible mispair"""
        assert not self.RNA('').can_mispair('')
        assert self.RNA('N').can_mispair('N')
        assert self.RNA('R').can_mispair('Y')
        assert self.RNA('N').can_mispair('r')
        assert self.RNA('CGUACGCAN').can_mispair('NUHCHUACH')
        assert self.RNA('U').can_mispair('C')
        assert self.RNA('U').can_mispair('R')
        assert self.RNA('UUU').can_mispair('AAR')
        assert self.RNA('UUU').can_mispair('GAG')
        assert not self.RNA('UUU').can_mispair('AAA')
        assert not self.RNA('UCAG').can_mispair('CUGA')
        assert self.RNA('U--').can_mispair('--U')

        assert self.DNA('TCCAAAGRYY').can_mispair('RRYCTTTGGA')

    def test_must_pair(self):
        """Sequence must_pair should return True when no possible mispairs"""
        assert self.RNA('').must_pair('')
        assert not self.RNA('N').must_pair('N')
        assert not self.RNA('R').must_pair('Y')
        assert not self.RNA('A').must_pair('A')
        assert not self.RNA('CGUACGCAN').must_pair('NUGCGUACG')
        assert not self.RNA('U').must_pair('C')
        assert not self.RNA('UUU').must_pair('AAR')
        assert not self.RNA('UUU').must_pair('RAA')
        assert not self.RNA('UU-').must_pair('-AA')
        assert self.RNA('UCAG').must_pair('CUGA')

        assert self.DNA('TCCAGGG').must_pair('CCCTGGA')
        assert self.DNA('tccaggg').must_pair(self.DNA('ccctgga'))
        assert not self.DNA('TCCAGGG').must_pair('NCCTGGA')

    def test_diff(self):
        """Sequence diff should count 1 for each difference between sequences"""
        self.assertEqual(self.RNA('UGCUGCUC').diff(''), 0)
        self.assertEqual(self.RNA('UGCUGCUC').diff('U'), 0)
        self.assertEqual(self.RNA('UGCUGCUC').diff('UCCCCCUC'), 3)
        # case-sensitive!
        self.assertEqual(self.RNA('AAAAA').diff('CCCCC'), 5)
        # raises TypeError if other not iterable
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
        # uses identity function by default
        self.assertEqual(self.RNA('UGCUGCUC').distance(''), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('U'), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('UCCCCCUC'), 3)
        # case-sensitive!
        self.assertEqual(self.RNA('AAAAA').distance('CCCCC'), 5)
        # should use function if supplied
        self.assertEqual(self.RNA('UGCUGCUC').distance('', f), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('U', f), 0)
        self.assertEqual(self.RNA('UGCUGCUC').distance('C', f), 1)
        self.assertEqual(self.RNA('UGCUGCUC').distance('G', f), 10)
        self.assertEqual(self.RNA('UGCUGCUC').distance('UCCCCCUC', f), 21)
        # case-sensitive!
        self.assertEqual(self.RNA('AAAAA').distance('CCCCC', f), 50)

    def test_matrix_distance(self):
        """Sequence matrix_distance should look up distances from a matrix"""
        # note that the score matrix must contain 'diagonal' elements m[i][i]
        # to avoid failure when the sequences match.
        m = {'U': {'U': 0, 'C': 1, 'A': 5}, 'C': {'C': 0, 'A': 2, 'G': 4}}
        self.assertEqual(self.RNA('UUUCCC').matrix_distance('UCACGG', m), 14)
        self.assertEqual(self.RNA('UUUCCC').matrix_distance('', m), 0)
        self.assertEqual(self.RNA('UUU').matrix_distance('CAC', m), 7)
        self.assertRaises(KeyError, self.RNA('UUU').matrix_distance, 'CAG', m)

    def test_frac_same(self):
        """Sequence frac_same should return similarity between sequences"""
        s1 = self.RNA('ACGU')
        s2 = self.RNA('AACG')
        s3 = self.RNA('GG')
        s4 = self.RNA('A')
        e = self.RNA('')
        self.assertEqual(s1.frac_same(e), 0)
        self.assertEqual(s1.frac_same(s2), 0.25)
        self.assertEqual(s1.frac_same(s3), 0)
        self.assertEqual(s1.frac_same(s4), 1.0)  # note truncation

    def test_frac_diff(self):
        """Sequence frac_diff should return difference between sequences"""
        s1 = self.RNA('ACGU')
        s2 = self.RNA('AACG')
        s3 = self.RNA('GG')
        s4 = self.RNA('A')
        e = self.RNA('')
        self.assertEqual(s1.frac_diff(e), 0)
        self.assertEqual(s1.frac_diff(s2), 0.75)
        self.assertEqual(s1.frac_diff(s3), 1)
        self.assertEqual(s1.frac_diff(s4), 0)  # note truncation

    def test_frac_same_gaps(self):
        """Sequence frac_same_gaps should return similarity in gap positions"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('GGGG')
        s3 = self.RNA('----')
        s4 = self.RNA('A-A-')
        s5 = self.RNA('-G-G')
        s6 = self.RNA('UU--')
        s7 = self.RNA('-')
        s8 = self.RNA('GGG')
        e = self.RNA('')
        self.assertEqual(s1.frac_same_gaps(s1), 1)
        self.assertEqual(s1.frac_same_gaps(s2), 1)
        self.assertEqual(s1.frac_same_gaps(s3), 0)
        self.assertEqual(s1.frac_same_gaps(s4), 0.5)
        self.assertEqual(s1.frac_same_gaps(s5), 0.5)
        self.assertEqual(s1.frac_same_gaps(s6), 0.5)
        self.assertEqual(s1.frac_same_gaps(s7), 0)
        self.assertEqual(s1.frac_same_gaps(e), 0)
        self.assertEqual(s3.frac_same_gaps(s3), 1)
        self.assertEqual(s3.frac_same_gaps(s4), 0.5)
        self.assertEqual(s3.frac_same_gaps(s7), 1.0)
        self.assertEqual(e.frac_same_gaps(e), 0.0)
        self.assertEqual(s4.frac_same_gaps(s5), 0.0)
        self.assertEqual(s4.frac_same_gaps(s6), 0.5)
        self.assertFloatEqual(s6.frac_same_gaps(s8), 2 / 3.0)

    def test_frac_diffGaps(self):
        """Sequence frac_diff_gaps should return difference in gap positions"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('GGGG')
        s3 = self.RNA('----')
        s4 = self.RNA('A-A-')
        s5 = self.RNA('-G-G')
        s6 = self.RNA('UU--')
        s7 = self.RNA('-')
        s8 = self.RNA('GGG')
        e = self.RNA('')
        self.assertEqual(s1.frac_diff_gaps(s1), 0)
        self.assertEqual(s1.frac_diff_gaps(s2), 0)
        self.assertEqual(s1.frac_diff_gaps(s3), 1)
        self.assertEqual(s1.frac_diff_gaps(s4), 0.5)
        self.assertEqual(s1.frac_diff_gaps(s5), 0.5)
        self.assertEqual(s1.frac_diff_gaps(s6), 0.5)
        self.assertEqual(s1.frac_diff_gaps(s7), 1)
        self.assertEqual(s1.frac_diff_gaps(e), 0)
        self.assertEqual(s3.frac_diff_gaps(s3), 0)
        self.assertEqual(s3.frac_diff_gaps(s4), 0.5)
        self.assertEqual(s3.frac_diff_gaps(s7), 0.0)
        self.assertEqual(e.frac_diff_gaps(e), 0.0)
        self.assertEqual(s4.frac_diff_gaps(s5), 1.0)
        self.assertEqual(s4.frac_diff_gaps(s6), 0.5)
        self.assertFloatEqual(s6.frac_diff_gaps(s8), 1 / 3.0)

    def test_frac_same_non_gaps(self):
        """Sequence frac_same_non_gaps should return similarities at non-gaps"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('AGGG')
        s3 = self.RNA('GGGG')
        s4 = self.RNA('AG--GA-G')
        s5 = self.RNA('CU--CU-C')
        s6 = self.RNA('AC--GC-G')
        s7 = self.RNA('--------')
        s8 = self.RNA('AAAA----')
        s9 = self.RNA('A-GG-A-C')
        e = self.RNA('')

        test = lambda x, y, z: self.assertFloatEqual(x.frac_same_non_gaps(y), z)
        test(s1, s2, 0.25)
        test(s1, s3, 0)
        test(s2, s3, 0.75)
        test(s1, s4, 0.5)
        test(s4, s5, 0)
        test(s4, s6, 0.6)
        test(s4, s7, 0)
        test(s4, s8, 0.5)
        test(s4, s9, 2 / 3.0)
        test(e, s4, 0)

    def test_frac_diffNonGaps(self):
        """Sequence frac_diff_non_gaps should return differences at non-gaps"""
        s1 = self.RNA('AAAA')
        s2 = self.RNA('AGGG')
        s3 = self.RNA('GGGG')
        s4 = self.RNA('AG--GA-G')
        s5 = self.RNA('CU--CU-C')
        s6 = self.RNA('AC--GC-G')
        s7 = self.RNA('--------')
        s8 = self.RNA('AAAA----')
        s9 = self.RNA('A-GG-A-C')
        e = self.RNA('')

        test = lambda x, y, z: self.assertFloatEqual(x.frac_diff_non_gaps(y), z)
        test(s1, s2, 0.75)
        test(s1, s3, 1)
        test(s2, s3, 0.25)
        test(s1, s4, 0.5)
        test(s4, s5, 1)
        test(s4, s6, 0.4)
        test(s4, s7, 0)
        test(s4, s8, 0.5)
        test(s4, s9, 1 / 3.0)
        test(e, s4, 0)

    def test_frac_similar(self):
        """Sequence frac_similar should return the fraction similarity"""
        transitions = dict.fromkeys([
            ('A', 'A'), ('A', 'G'), ('G', 'A'), ('G', 'G'),
            ('U', 'U'), ('U', 'C'), ('C', 'U'), ('C', 'C')])

        s1 = self.RNA('UCAGGCAA')
        s2 = self.RNA('CCAAAUGC')
        s3 = self.RNA('GGGGGGGG')
        e = self.RNA('')

        test = lambda x, y, z: self.assertFloatEqual(
            x.frac_similar(y, transitions), z)

        test(e, e, 0)
        test(s1, e, 0)
        test(s1, s1, 1)
        test(s1, s2, 7.0 / 8)
        test(s1, s3, 5.0 / 8)
        test(s2, s3, 4.0 / 8)

    def test_with_termini_unknown(self):
        """with_termini_unknown should reset termini to unknown char"""
        s1 = self.RNA('-?--AC--?-')
        s2 = self.RNA('AC')
        self.assertEqual(s1.with_termini_unknown(), '????AC????')
        self.assertEqual(s2.with_termini_unknown(), 'AC')

    def test_consistent_gap_degen_handling(self):
        """gap degen character should be treated consistently"""
        # the degen character '?' can be a gap, so when we strip either gaps or
        # degen characters it should be gone too
        raw_seq = "---??-??TC-GGCG-GCA-G-GC-?-C-TAN-GCGC-CCTC-AGGA?-???-??--"
        raw_ungapped = re.sub("[-?]", "", raw_seq)
        raw_no_ambigs = re.sub("[N?]+", "", raw_seq)
        dna = self.DNA(raw_seq)
        self.assertEqual(dna.degap(), raw_ungapped)
        self.assertEqual(dna.strip_degenerate(), raw_no_ambigs)
        self.assertEqual(dna.strip_bad_and_gaps(), raw_ungapped)


class SequenceSubclassTests(TestCase):
    """Only one general set of tests, since the subclasses are very thin."""

    def test_DnaSequence(self):
        """DnaSequence should behave as expected"""
        x = DnaSequence('tcag')
        # note: no longer preserves case
        self.assertEqual(x, 'TCAG')

        x = DnaSequence('aaa') + DnaSequence('ccc')
        # note: doesn't preserve case
        self.assertEqual(x, 'AAACCC')
        assert x.MolType is DNA
        self.assertRaises(AlphabetError, x.__add__, 'z')
        self.assertEqual(DnaSequence('TTTAc').rc(), 'GTAAA')


class ModelSequenceTests(object):
    """Base class for tests of specific ModelSequence objects."""
    SequenceClass = None  # override in derived classes

    def test_to_fasta(self):
        """Sequence to_fasta() should return Fasta-format string"""
        even = 'TCAGAT'
        odd = even + 'AAA'
        even_dna = self.SequenceClass(even, Name='even')
        odd_dna = self.SequenceClass(odd, Name='odd')
        self.assertEqual(even_dna.to_fasta(), '>even\nTCAGAT')
        # set line wrap to small number so we can test that it works
        even_dna.LineWrap = 2
        self.assertEqual(even_dna.to_fasta(), '>even\nTC\nAG\nAT')
        odd_dna.LineWrap = 2
        self.assertEqual(odd_dna.to_fasta(), '>odd\nTC\nAG\nAT\nAA\nA')
        # check that changing the linewrap again works
        even_dna.LineWrap = 4
        self.assertEqual(even_dna.to_fasta(), '>even\nTCAG\nAT')

    def test_to_phylip(self):
        """Sequence to_phylip() should return one-line phylip string"""
        s = self.SequenceClass('ACG', Name='xyz')
        self.assertEqual(s.to_phylip(), 'xyz' + ' ' * 27 + 'ACG')


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
        self.assertEqual(r._data, array([0, 1, 2, 3, 3, 2]))
        self.assertEqual(str(r), orig)

    def test_toKwords(self):
        """Sequence toKwords should give expected counts"""
        orig = 'ATCCCTAGC'
        r = self.SequenceClass(orig)
        # if we use k = 1, should just get the characters
        w = r.toKwords(1)
        self.assertEqual(w, r._data)
        w = r.toKwords(1, overlapping=False)
        self.assertEqual(w, r._data)
        # if we use k = 2, should get overlapping or nonoverlapping k-words
        w = r.toKwords(2)
        self.assertEqual(w, array([8, 1, 5, 5, 4, 2, 11, 13]))
        w = r.toKwords(2, overlapping=False)
        self.assertEqual(w, array([8, 5, 4, 11]))
        # check a case with k = 3, i.e. codons
        w = r.toKwords(3, overlapping=False)
        self.assertEqual(w, array([33, 20, 45]))


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
        self.assertEqual(r._data, array([6, 62]))
        self.assertEqual(str(r), orig)

    def test_toKwords(self):
        """Sequence toKwords should give expected counts"""
        orig = 'ATCCCTAGC'
        r = self.SequenceClass(orig)
        # if we use k = 1, should just get the characters
        w = r.toKwords(1)
        self.assertEqual(w, r._data)
        w = r.toKwords(1, overlapping=False)
        self.assertEqual(w, r._data)
        # if we use k = 2, should get overlapping or nonoverlapping k-words
        w = r.toKwords(2)
        self.assertEqual(w, array([2132, 1325]))
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
        self.assertEqual(sc('TC').gaps(), array([0, 0]))
        self.assertEqual(sc('T-').gaps(), array([0, 1]))

    def test_degap(self):
        """Gapped sequence degap() should return correct array"""
        sc = self.SequenceClass
        self.assertEqual(sc('T-').degap(), sc('T'))

    def test_nongaps(self):
        """Gapped sequence nongaps() should return correct array"""
        sc = self.SequenceClass
        self.assertEqual(sc('TC').nongaps(), array([1, 1]))
        self.assertEqual(sc('T-').nongaps(), array([1, 0]))

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
        self.assertEqual(d._data, array([0, 28]))
        self.assertEqual(str(d.to_rna()), 'UUUCGU')
        self.assertEqual(str(d.to_dna()), 'TTTCGT')

    def test_ModelRnaCodonSequence(self):
        """ModelRnaCodonSequence should behave as expected"""
        r = ModelRnaCodonSequence('UUUCGU')
        self.assertEqual(str(r), 'UUUCGU')
        self.assertEqual(r._data, array([0, 28]))
        self.assertEqual(str(r.to_rna()), 'UUUCGU')
        self.assertEqual(str(r.to_dna()), 'TTTCGT')


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

        def f(x, y):
            if x == 2 or y == 2:
                return 10
            return 0
        self.assertEqual(s1.distance(s2, f, use_indices=True), 20)

    def test_strip_bad(self):
        """Sequence strip_bad should remove any non-base, non-gap chars"""
        # have to turn off check to get bad data in; no longer preserves case
        r = self.RNA('UCAGRYU')
        r._data[0] = 31
        r._data[2] = 55
        self.assertEqual(r.strip_bad(), 'CGRYU')

    def test_strip_bad_and_gaps(self):
        """Sequence strip_bad_and_gaps should remove gaps and bad chars"""
        # have to turn off check to get bad data in; no longer preserves case
        r = self.RNA('ACG--GRN?')
        self.assertEqual(r.strip_bad_and_gaps(), 'ACGGRN')
        r._data[0] = 99
        self.assertEqual(r.strip_bad_and_gaps(), 'CGGRN')

    def test_gap_array(self):
        """Sequence gap_array should return array of gaps"""
        r = self.RNA('-?A-?NRY-')
        v = r.gap_array()
        self.assertEqual(v, array([1, 1, 0, 1, 1, 0, 0, 0, 1]))
        r = self.RNA('AC')
        v = r.gap_array()
        self.assertEqual(v, array([0, 0]))
        r = self.RNA('-?')
        v = r.gap_array()
        self.assertEqual(v, array([1, 1]))

    def test_gapIndices(self):
        """Sequence gapIndices should return positions of gaps"""
        r = self.RNA('-?A-?NRY-')
        v = r.gapIndices()
        self.assertEqual(v, array([0, 1, 3, 4, 8]))
        r = self.RNA('AC')
        v = r.gapIndices()
        self.assertEqual(v, array([]))  # note: always returns array
        r = self.RNA('-?')
        v = r.gapIndices()
        self.assertEqual(v, array([0, 1]))


# run if called from command-line
if __name__ == "__main__":
    main()
