"""Unit tests for Sequence class and its subclasses.
"""
import json
import os
import re

from pickle import dumps
from unittest import TestCase

import numpy
import pytest

from numpy import array
from numpy.testing import assert_allclose, assert_equal

import cogent3

from cogent3._version import __version__
from cogent3.core.alignment import Aligned
from cogent3.core.moltype import (
    ASCII,
    BYTES,
    DNA,
    RNA,
    AlphabetError,
    get_moltype,
)
from cogent3.core.sequence import (
    ABSequence,
    ArrayCodonSequence,
    ArrayDnaCodonSequence,
    ArrayDnaSequence,
    ArrayNucleicAcidSequence,
    ArrayProteinSequence,
    ArrayRnaCodonSequence,
    ArrayRnaSequence,
    ArraySequence,
    DnaSequence,
    ProteinSequence,
    RnaSequence,
    Sequence,
    SeqView,
    _coerce_to_seqview,
)
from cogent3.util.misc import get_object_provenance


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
        self.assertEqual(s, "")
        assert s.moltype in (ASCII, BYTES)

        r = self.RNA()
        assert r.moltype is RNA

    def test_init_data(self):
        """Sequence init with data should set data in correct location"""
        r = self.RNA("ucagg")
        # no longer preserves case
        self.assertEqual(r, "UCAGG")

    def test_init_from_bytes(self):
        """correctly convert bytes to str"""
        s = self.SEQ(b"ACGT")
        self.assertEqual(s, "ACGT")

    def test_init_other_seq(self):
        """Sequence init with other seq should preserve name and info."""
        r = self.RNA("UCAGG", name="x", info={"z": 3})
        s = Sequence(r)
        self.assertEqual(str(s), "UCAGG")
        self.assertEqual(s.name, "x")
        self.assertEqual(s.info.z, 3)

    def test_copy(self):
        """correctly returns a copy version of self"""
        s = Sequence("TTTTTTTTTTAAAA", name="test_copy")
        annot1 = s.add_feature(biotype="exon", name="annot1", spans=[(0, 10)])
        annot2 = s.add_feature(biotype="exon", name="annot2", spans=[(10, 14)])
        got = s.copy()
        got_annot1 = list(got.get_features(biotype="exon", name="annot1"))[0]
        got_annot2 = list(got.get_features(biotype="exon", name="annot2"))[0]
        self.assertIsNot(got, s)
        self.assertIsNot(got_annot1, annot1)
        self.assertIsNot(got_annot2, annot2)
        self.assertEqual(got.name, s.name)
        self.assertEqual(got.info, s.info)
        self.assertEqual(str(got), str(s))
        self.assertEqual(got.moltype, s.moltype)
        annot1_slice = str(annot1.get_slice())
        annot2_slice = str(annot2.get_slice())
        got1_slice = str(got_annot1.get_slice())
        got2_slice = str(got_annot2.get_slice())
        self.assertEqual(annot1_slice, got1_slice)
        self.assertEqual(annot2_slice, got2_slice)

    def test_compare_to_string(self):
        """Sequence should compare equal to same string."""
        r = self.RNA("UCC")
        self.assertEqual(r, "UCC")

    def test_slice(self):
        """Sequence slicing should work as expected"""
        r = self.RNA("UCAGG")
        self.assertEqual(r[0], "U")
        self.assertEqual(r[-1], "G")
        self.assertEqual(r[1:3], "CA")

    def test_to_dna(self):
        """Returns copy of self as DNA."""
        r = self.RNA("UCA")
        self.assertEqual(str(r), "UCA")
        self.assertEqual(str(r.to_dna()), "TCA")

    def test_to_rna(self):
        """Returns copy of self as RNA."""
        r = self.DNA("TCA")
        self.assertEqual(str(r), "TCA")
        self.assertEqual(str(r.to_rna()), "UCA")

    def test_to_fasta(self):
        """Sequence to_fasta() should return Fasta-format string"""
        even = "TCAGAT"
        odd = even + "AAA"
        even_dna = self.SEQ(even, name="even")
        odd_dna = self.SEQ(odd, name="odd")
        self.assertEqual(even_dna.to_fasta(), ">even\nTCAGAT\n")
        # set line wrap to small number so we can test that it works
        self.assertEqual(even_dna.to_fasta(block_size=2), ">even\nTC\nAG\nAT\n")
        self.assertEqual(odd_dna.to_fasta(block_size=2), ">odd\nTC\nAG\nAT\nAA\nA\n")
        # check that changing the linewrap again works
        self.assertEqual(even_dna.to_fasta(block_size=4), ">even\nTCAG\nAT\n")

    def test_serialize(self):
        """Sequence should be serializable"""
        r = self.RNA("ugagg")
        assert dumps(r)

    def test_sequence_to_moltype(self):
        """correctly convert to specified moltype"""
        s = Sequence("TTTTTTTTTTAAAA", name="test1")
        s.add_feature(biotype="exon", name="fred", spans=[(0, 10)])
        s.add_feature(biotype="exon", name="trev", spans=[(10, 14)])
        got = s.to_moltype("rna")
        fred = list(got.get_features(name="fred"))[0]
        assert str(got[fred]) == "UUUUUUUUUU"
        trev = list(got.get_features(name="trev"))[0]
        assert str(got[trev]) == "AAAA"

        # calling with a null object should raise an exception
        with self.assertRaises(ValueError):
            s.to_moltype(None)

        with self.assertRaises(ValueError):
            s.to_moltype("")

    def test_strip_degenerate(self):
        """Sequence strip_degenerate should remove any degenerate bases"""
        self.assertEqual(self.RNA("UCAG-").strip_degenerate(), "UCAG-")
        self.assertEqual(self.RNA("NRYSW").strip_degenerate(), "")
        self.assertEqual(self.RNA("USNG").strip_degenerate(), "UG")

    def test_strip_bad(self):
        """Sequence strip_bad should remove any non-base, non-gap chars"""
        # have to turn off check to get bad data in; no longer preserves case
        self.assertEqual(
            self.RNA("UCxxxAGwsnyrHBNzzzD-D", check=False).strip_bad(),
            "UCAGWSNYRHBND-D",
        )
        self.assertEqual(self.RNA("@#^*($@!#&()!@QZX", check=False).strip_bad(), "")
        self.assertEqual(
            self.RNA("aaaxggg---!ccc", check=False).strip_bad(), "AAAGGG---CCC"
        )

    def test_strip_bad_and_gaps(self):
        """Sequence strip_bad_and_gaps should remove gaps and bad chars"""
        # have to turn off check to get bad data in; no longer preserves case
        self.assertEqual(
            self.RNA("UxxCAGwsnyrHBNz#!D-D", check=False).strip_bad_and_gaps(),
            "UCAGWSNYRHBNDD",
        )
        self.assertEqual(
            self.RNA("@#^*($@!#&()!@QZX", check=False).strip_bad_and_gaps(), ""
        )
        self.assertEqual(
            self.RNA("aaa ggg ---!ccc", check=False).strip_bad_and_gaps(), "AAAGGGCCC"
        )

    def test_shuffle(self):
        """Sequence shuffle should return new random sequence w/ same monomers"""
        r = self.RNA("UUUUCCCCAAAAGGGG")
        s = r.shuffle()
        self.assertFalse(r == s)
        self.assertCountEqual(r, s)

    def test_complement(self):
        """Sequence complement should correctly complement sequence"""
        self.assertEqual(self.RNA("UAUCG-NR").complement(), "AUAGC-NY")
        self.assertEqual(self.DNA("TATCG-NR").complement(), "ATAGC-NY")
        self.assertEqual(self.DNA("").complement(), "")
        self.assertRaises(TypeError, self.PROT("ACD").complement)

    def test_rc(self):
        """Sequence rc should correctly reverse-complement sequence"""
        # no longer preserves case!
        self.assertEqual(self.RNA("UauCG-NR").rc(), "YN-CGAUA")
        self.assertEqual(self.DNA("TatCG-NR").rc(), "YN-CGATA")
        self.assertEqual(self.RNA("").rc(), "")
        self.assertEqual(self.RNA("A").rc(), "U")
        self.assertRaises(TypeError, self.PROT("ACD").rc)

    def test_contains(self):
        """Sequence contains should return correct result"""
        r = self.RNA("UCA")
        assert "U" in r
        assert "CA" in r
        assert "X" not in r
        assert "G" not in r

    def test_iter(self):
        """Sequence iter should iterate over sequence"""
        p = self.PROT("QWE")
        self.assertEqual(list(p), ["Q", "W", "E"])

    def test_is_gapped(self):
        """Sequence is_gapped should return True if gaps in seq"""
        assert not self.RNA("").is_gapped()
        assert not self.RNA("ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN").is_gapped()
        assert self.RNA("-").is_gapped()
        assert self.PROT("--").is_gapped()
        assert self.RNA("CAGUCGUACGUCAGUACGUacucauacgac-caguACUG").is_gapped()
        assert self.RNA("CA--CGUAUGCA-----g").is_gapped()
        assert self.RNA("CAGU-").is_gapped()

    def test_is_gap(self):
        """Sequence is_gap should return True if char is a valid gap char"""
        r = self.RNA("ACGUCAGUACGUCAGNRCGAUcaguaguacYRNRYRN")
        for char in "qwertyuiopasdfghjklzxcvbnmQWERTYUIOASDFGHJKLZXCVBNM":
            assert not r.is_gap(char)
        assert r.is_gap("-")
        # only works on a single literal that's a gap, not on a sequence.
        # possibly, this behavior should change?
        assert not r.is_gap("---")
        # check behaviour on self
        assert not self.RNA("CGAUACGUACGACU").is_gap()
        assert not self.RNA("---CGAUA----CGUACG---ACU---").is_gap()
        assert self.RNA("").is_gap()
        assert self.RNA("----------").is_gap()

    def test_is_degenerate(self):
        """Sequence is_degenerate should return True if degen symbol in seq"""
        assert not self.RNA("").is_degenerate()
        assert not self.RNA("UACGCUACAUGuacgucaguGCUAGCUA---ACGUCAG").is_degenerate()
        assert self.RNA("N").is_degenerate()
        assert self.RNA("R").is_degenerate()
        assert self.RNA("y").is_degenerate()
        assert self.RNA("GCAUguagcucgUCAGUCAGUACgUgcasCUAG").is_degenerate()
        assert self.RNA("ACGYAUGCUGYWWNMNuwbycwuybcwbwub").is_degenerate()

    def test_is_strict(self):
        """Sequence is_strict should return True if all symbols in Monomers"""
        assert self.RNA("").is_strict()
        assert self.PROT("A").is_strict()
        assert self.RNA("UAGCACUgcaugcauGCAUGACuacguACAUG").is_strict()
        assert not self.RNA("CAGUCGAUCA-cgaucagUCGAUGAC").is_strict()

    def test_first_gap(self):
        """Sequence first_gap should return index of first gap symbol, or None"""
        self.assertEqual(self.RNA("").first_gap(), None)
        self.assertEqual(self.RNA("a").first_gap(), None)
        self.assertEqual(self.RNA("uhacucHuhacUUhacan").first_gap(), None)
        self.assertEqual(self.RNA("-abc").first_gap(), 0)
        self.assertEqual(self.RNA("b-ac").first_gap(), 1)
        self.assertEqual(self.RNA("abcd-").first_gap(), 4)

    def test_first_degenerate(self):
        """Sequence first_degenerate should return index of first degen symbol"""
        self.assertEqual(self.RNA("").first_degenerate(), None)
        self.assertEqual(self.RNA("a").first_degenerate(), None)
        self.assertEqual(self.RNA("UCGACA--CU-gacucaguacgua").first_degenerate(), None)
        self.assertEqual(self.RNA("nCAGU").first_degenerate(), 0)
        self.assertEqual(self.RNA("CUGguagvAUG").first_degenerate(), 7)
        self.assertEqual(self.RNA("ACUGCUAacgud").first_degenerate(), 11)

    def test_first_non_strict(self):
        """Sequence first_non_strict should return index of first non-strict symbol"""
        self.assertEqual(self.RNA("").first_non_strict(), None)
        self.assertEqual(self.RNA("A").first_non_strict(), None)
        self.assertEqual(self.RNA("ACGUACGUcgaucagu").first_non_strict(), None)
        self.assertEqual(self.RNA("N").first_non_strict(), 0)
        self.assertEqual(self.RNA("-").first_non_strict(), 0)
        self.assertEqual(self.RNA("ACGUcgAUGUGCAUcagu-").first_non_strict(), 18)

    def test_disambiguate(self):
        """Sequence disambiguate should remove degenerate bases"""
        self.assertEqual(self.RNA("").disambiguate(), "")
        self.assertEqual(
            self.RNA("AGCUGAUGUA--CAGU").disambiguate(), "AGCUGAUGUA--CAGU"
        )
        self.assertEqual(
            self.RNA("AUn-yrs-wkmCGwmrNMWRKY").disambiguate("strip"), "AU--CG"
        )
        s = self.RNA("AUn-yrs-wkmCGwmrNMWRKY")
        t = s.disambiguate("random")
        u = s.disambiguate("random")
        for i, j in zip(str(s), str(t)):
            if i in s.moltype.degenerates:
                assert j in s.moltype.degenerates[i]
            else:
                assert i == j
        self.assertFalse(t == u)
        self.assertEqual(len(s), len(t))

    def test_degap(self):
        """Sequence degap should remove all gaps from sequence"""
        # doesn't preserve case
        self.assertEqual(self.RNA("").degap(), "")
        self.assertEqual(
            self.RNA("GUCAGUCgcaugcnvuncdks").degap(), "GUCAGUCGCAUGCNVUNCDKS"
        )
        self.assertEqual(self.RNA("----------------").degap(), "")
        self.assertEqual(self.RNA("gcuauacg-").degap(), "GCUAUACG")
        self.assertEqual(self.RNA("-CUAGUCA").degap(), "CUAGUCA")
        self.assertEqual(self.RNA("---a---c---u----g---").degap(), "ACUG")
        self.assertEqual(self.RNA("?a-").degap(), "A")

    def test_gap_indices(self):
        """Sequence gap_indices should return correct gap positions"""
        self.assertEqual(self.RNA("").gap_indices(), [])
        self.assertEqual(self.RNA("ACUGUCAGUACGHSDKCUCDNNS").gap_indices(), [])
        self.assertEqual(self.RNA("GUACGUACAKDC-SDHDSK").gap_indices(), [12])
        self.assertEqual(self.RNA("-DSHUHDS").gap_indices(), [0])
        self.assertEqual(self.RNA("UACHASADS-").gap_indices(), [9])
        self.assertEqual(
            self.RNA("---CGAUgCAU---ACGHc---ACGUCAGU---").gap_indices(),
            [0, 1, 2, 11, 12, 13, 19, 20, 21, 30, 31, 32],
        )

    def test_gap_vector(self):
        """Sequence gap_vector should return correct gap positions"""

        def g(x):
            return self.RNA(x).gap_vector()

        self.assertEqual(g(""), [])
        self.assertEqual(g("ACUGUCAGUACGHCSDKCCUCCDNCNS"), [False] * 27)
        self.assertEqual(
            g("GUACGUAACAKADC-SDAHADSAK"),
            list(map(bool, list(map(int, "000000000000001000000000")))),
        )
        self.assertEqual(g("-DSHSUHDSS"), list(map(bool, list(map(int, "1000000000")))))
        self.assertEqual(
            g("UACHASCAGDS-"), list(map(bool, list(map(int, "000000000001"))))
        )
        self.assertEqual(
            g("---CGAUgCAU---ACGHc---ACGUCAGU--?"),
            list(map(bool, list(map(int, "111000000001110000011100000000111")))),
        )

    def test_gap_maps(self):
        """Sequence gap_maps should return dicts mapping gapped/ungapped pos"""
        empty = ""
        no_gaps = "aaa"
        all_gaps = "---"
        start_gaps = "--abc"
        end_gaps = "ab---"
        mid_gaps = "--a--b-cd---"

        def gm(x):
            return self.RNA(x).gap_maps()

        self.assertEqual(gm(empty), ({}, {}))
        self.assertEqual(gm(no_gaps), ({0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2}))
        self.assertEqual(gm(all_gaps), ({}, {}))
        self.assertEqual(gm(start_gaps), ({0: 2, 1: 3, 2: 4}, {2: 0, 3: 1, 4: 2}))
        self.assertEqual(gm(end_gaps), ({0: 0, 1: 1}, {0: 0, 1: 1}))
        self.assertEqual(
            gm(mid_gaps), ({0: 2, 1: 5, 2: 7, 3: 8}, {2: 0, 5: 1, 7: 2, 8: 3})
        )

    def test_count_gaps(self):
        """Sequence count_gaps should return correct gap count"""
        self.assertEqual(self.RNA("").count_gaps(), 0)
        self.assertEqual(self.RNA("ACUGUCAGUACGHSDKCUCDNNS").count_gaps(), 0)
        self.assertEqual(self.RNA("GUACGUACAKDC-SDHDSK").count_gaps(), 1)
        self.assertEqual(self.RNA("-DSHUHDS").count_gaps(), 1)
        self.assertEqual(self.RNA("UACHASADS-").count_gaps(), 1)
        self.assertEqual(self.RNA("---CGAUgCAU---ACGHc---ACGUCAGU---").count_gaps(), 12)

    def test_count_degenerate(self):
        """Sequence count_degenerate should return correct degen base count"""
        self.assertEqual(self.RNA("").count_degenerate(), 0)
        self.assertEqual(self.RNA("GACUGCAUGCAUCGUACGUCAGUACCGA").count_degenerate(), 0)
        self.assertEqual(self.RNA("N").count_degenerate(), 1)
        self.assertEqual(self.PROT("N").count_degenerate(), 0)
        self.assertEqual(self.RNA("NRY").count_degenerate(), 3)
        self.assertEqual(
            self.RNA("ACGUAVCUAGCAUNUCAGUCAGyUACGUCAGS").count_degenerate(), 4
        )

    def test_possibilites(self):
        """Sequence possibilities should return correct # possible sequences"""
        self.assertEqual(self.RNA("").possibilities(), 1)
        self.assertEqual(self.RNA("ACGUgcaucagUCGuGCAU").possibilities(), 1)
        self.assertEqual(self.RNA("N").possibilities(), 4)
        self.assertEqual(self.RNA("R").possibilities(), 2)
        self.assertEqual(self.RNA("H").possibilities(), 3)
        self.assertEqual(self.RNA("nRh").possibilities(), 24)
        self.assertEqual(
            self.RNA("AUGCnGUCAg-aurGauc--gauhcgauacgws").possibilities(), 96
        )

    def test_mw(self):
        """Sequence MW should return correct molecular weight"""
        self.assertEqual(self.PROT("").mw(), 0)
        self.assertEqual(self.RNA("").mw(), 0)
        assert_allclose(self.PROT("A").mw(), 89.09)
        assert_allclose(self.RNA("A").mw(), 375.17)
        assert_allclose(self.PROT("AAA").mw(), 231.27)
        assert_allclose(self.RNA("AAA").mw(), 1001.59)
        assert_allclose(self.RNA("AAACCCA").mw(), 2182.37)

    def test_can_match(self):
        """Sequence can_match should return True if all positions can match"""
        assert self.RNA("").can_match("")
        assert self.RNA("UCAG").can_match("UCAG")
        assert not self.RNA("UCAG").can_match("ucag")
        assert self.RNA("UCAG").can_match("NNNN")
        assert self.RNA("NNNN").can_match("UCAG")
        assert self.RNA("NNNN").can_match("NNNN")
        assert not self.RNA("N").can_match("x")
        assert not self.RNA("N").can_match("-")
        assert self.RNA("UCAG").can_match("YYRR")
        assert self.RNA("UCAG").can_match("KMWS")

    def test_can_mismatch(self):
        """Sequence can_mismatch should return True on any possible mismatch"""
        assert not self.RNA("").can_mismatch("")
        assert self.RNA("N").can_mismatch("N")
        assert self.RNA("R").can_mismatch("R")
        assert self.RNA("N").can_mismatch("r")
        assert self.RNA("CGUACGCAN").can_mismatch("CGUACGCAN")
        assert self.RNA("U").can_mismatch("C")
        assert self.RNA("UUU").can_mismatch("UUC")
        assert self.RNA("UUU").can_mismatch("UUY")
        assert not self.RNA("UUU").can_mismatch("UUU")
        assert not self.RNA("UCAG").can_mismatch("UCAG")
        assert not self.RNA("U--").can_mismatch("U--")

    def test_must_match(self):
        """Sequence must_match should return True when no possible mismatches"""
        assert self.RNA("").must_match("")
        assert not self.RNA("N").must_match("N")
        assert not self.RNA("R").must_match("R")
        assert not self.RNA("N").must_match("r")
        assert not self.RNA("CGUACGCAN").must_match("CGUACGCAN")
        assert not self.RNA("U").must_match("C")
        assert not self.RNA("UUU").must_match("UUC")
        assert not self.RNA("UUU").must_match("UUY")
        assert self.RNA("UU-").must_match("UU-")
        assert self.RNA("UCAG").must_match("UCAG")

    def test_can_pair(self):
        """Sequence can_pair should return True if all positions can pair"""
        assert self.RNA("").can_pair("")
        assert not self.RNA("UCAG").can_pair("UCAG")
        assert self.RNA("UCAG").can_pair("CUGA")
        assert not self.RNA("UCAG").can_pair("cuga")
        assert self.RNA("UCAG").can_pair("NNNN")
        assert self.RNA("NNNN").can_pair("UCAG")
        assert self.RNA("NNNN").can_pair("NNNN")
        assert not self.RNA("N").can_pair("x")
        assert not self.RNA("N").can_pair("-")
        assert self.RNA("-").can_pair("-")
        assert self.RNA("UCAGU").can_pair("KYYRR")
        assert self.RNA("UCAG").can_pair("KKRS")
        assert self.RNA("U").can_pair("G")

        assert not self.DNA("T").can_pair("G")

    def test_can_mispair(self):
        """Sequence can_mispair should return True on any possible mispair"""
        assert not self.RNA("").can_mispair("")
        assert self.RNA("N").can_mispair("N")
        assert self.RNA("R").can_mispair("Y")
        assert self.RNA("N").can_mispair("r")
        assert self.RNA("CGUACGCAN").can_mispair("NUHCHUACH")
        assert self.RNA("U").can_mispair("C")
        assert self.RNA("U").can_mispair("R")
        assert self.RNA("UUU").can_mispair("AAR")
        assert self.RNA("UUU").can_mispair("GAG")
        assert not self.RNA("UUU").can_mispair("AAA")
        assert not self.RNA("UCAG").can_mispair("CUGA")
        assert self.RNA("U--").can_mispair("--U")

        assert self.DNA("TCCAAAGRYY").can_mispair("RRYCTTTGGA")

    def test_must_pair(self):
        """Sequence must_pair should return True when no possible mispairs"""
        assert self.RNA("").must_pair("")
        assert not self.RNA("N").must_pair("N")
        assert not self.RNA("R").must_pair("Y")
        assert not self.RNA("A").must_pair("A")
        assert not self.RNA("CGUACGCAN").must_pair("NUGCGUACG")
        assert not self.RNA("U").must_pair("C")
        assert not self.RNA("UUU").must_pair("AAR")
        assert not self.RNA("UUU").must_pair("RAA")
        assert not self.RNA("UU-").must_pair("-AA")
        assert self.RNA("UCAG").must_pair("CUGA")

        assert self.DNA("TCCAGGG").must_pair("CCCTGGA")
        assert self.DNA("tccaggg").must_pair(self.DNA("ccctgga"))
        assert not self.DNA("TCCAGGG").must_pair("NCCTGGA")

    def test_diff(self):
        """Sequence diff should count 1 for each difference between sequences"""
        self.assertEqual(self.RNA("UGCUGCUC").diff(""), 0)
        self.assertEqual(self.RNA("UGCUGCUC").diff("U"), 0)
        self.assertEqual(self.RNA("UGCUGCUC").diff("UCCCCCUC"), 3)
        # case-sensitive!
        self.assertEqual(self.RNA("AAAAA").diff("CCCCC"), 5)
        # raises TypeError if other not iterable
        self.assertRaises(TypeError, self.RNA("AAAAA").diff, 5)

    def test_distance(self):
        """Sequence distance should calculate correctly based on function"""

        def f(a, b):
            if a == b:
                return 0
            if (a in "UC" and b in "UC") or (a in "AG" and b in "AG"):
                return 1
            else:
                return 10

        # uses identity function by default
        self.assertEqual(self.RNA("UGCUGCUC").distance(""), 0)
        self.assertEqual(self.RNA("UGCUGCUC").distance("U"), 0)
        self.assertEqual(self.RNA("UGCUGCUC").distance("UCCCCCUC"), 3)
        # case-sensitive!
        self.assertEqual(self.RNA("AAAAA").distance("CCCCC"), 5)
        # should use function if supplied
        self.assertEqual(self.RNA("UGCUGCUC").distance("", f), 0)
        self.assertEqual(self.RNA("UGCUGCUC").distance("U", f), 0)
        self.assertEqual(self.RNA("UGCUGCUC").distance("C", f), 1)
        self.assertEqual(self.RNA("UGCUGCUC").distance("G", f), 10)
        self.assertEqual(self.RNA("UGCUGCUC").distance("UCCCCCUC", f), 21)
        # case-sensitive!
        self.assertEqual(self.RNA("AAAAA").distance("CCCCC", f), 50)

    def test_matrix_distance(self):
        """Sequence matrix_distance should look up distances from a matrix"""
        # note that the score matrix must contain 'diagonal' elements m[i][i]
        # to avoid failure when the sequences match.
        m = {"U": {"U": 0, "C": 1, "A": 5}, "C": {"C": 0, "A": 2, "G": 4}}
        self.assertEqual(self.RNA("UUUCCC").matrix_distance("UCACGG", m), 14)
        self.assertEqual(self.RNA("UUUCCC").matrix_distance("", m), 0)
        self.assertEqual(self.RNA("UUU").matrix_distance("CAC", m), 7)
        self.assertRaises(KeyError, self.RNA("UUU").matrix_distance, "CAG", m)

    def test_frac_same(self):
        """Sequence frac_same should return similarity between sequences"""
        s1 = self.RNA("ACGU")
        s2 = self.RNA("AACG")
        s3 = self.RNA("GG")
        s4 = self.RNA("A")
        e = self.RNA("")
        self.assertEqual(s1.frac_same(e), 0)
        self.assertEqual(s1.frac_same(s2), 0.25)
        self.assertEqual(s1.frac_same(s3), 0)
        self.assertEqual(s1.frac_same(s4), 1.0)  # note truncation

    def test_frac_diff(self):
        """Sequence frac_diff should return difference between sequences"""
        s1 = self.RNA("ACGU")
        s2 = self.RNA("AACG")
        s3 = self.RNA("GG")
        s4 = self.RNA("A")
        e = self.RNA("")
        self.assertEqual(s1.frac_diff(e), 0)
        self.assertEqual(s1.frac_diff(s2), 0.75)
        self.assertEqual(s1.frac_diff(s3), 1)
        self.assertEqual(s1.frac_diff(s4), 0)  # note truncation

    def test_frac_same_gaps(self):
        """Sequence frac_same_gaps should return similarity in gap positions"""
        s1 = self.RNA("AAAA")
        s2 = self.RNA("GGGG")
        s3 = self.RNA("----")
        s4 = self.RNA("A-A-")
        s5 = self.RNA("-G-G")
        s6 = self.RNA("UU--")
        s7 = self.RNA("-")
        s8 = self.RNA("GGG")
        e = self.RNA("")
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
        assert_allclose(s6.frac_same_gaps(s8), 2 / 3.0)

    def test_frac_diffGaps(self):
        """Sequence frac_diff_gaps should return difference in gap positions"""
        s1 = self.RNA("AAAA")
        s2 = self.RNA("GGGG")
        s3 = self.RNA("----")
        s4 = self.RNA("A-A-")
        s5 = self.RNA("-G-G")
        s6 = self.RNA("UU--")
        s7 = self.RNA("-")
        s8 = self.RNA("GGG")
        e = self.RNA("")
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
        assert_allclose(s6.frac_diff_gaps(s8), 1 / 3.0)

    def test_frac_same_non_gaps(self):
        """Sequence frac_same_non_gaps should return similarities at non-gaps"""
        s1 = self.RNA("AAAA")
        s2 = self.RNA("AGGG")
        s3 = self.RNA("GGGG")
        s4 = self.RNA("AG--GA-G")
        s5 = self.RNA("CU--CU-C")
        s6 = self.RNA("AC--GC-G")
        s7 = self.RNA("--------")
        s8 = self.RNA("AAAA----")
        s9 = self.RNA("A-GG-A-C")
        e = self.RNA("")

        def test(x, y, z):
            return assert_allclose(x.frac_same_non_gaps(y), z)

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
        s1 = self.RNA("AAAA")
        s2 = self.RNA("AGGG")
        s3 = self.RNA("GGGG")
        s4 = self.RNA("AG--GA-G")
        s5 = self.RNA("CU--CU-C")
        s6 = self.RNA("AC--GC-G")
        s7 = self.RNA("--------")
        s8 = self.RNA("AAAA----")
        s9 = self.RNA("A-GG-A-C")
        e = self.RNA("")

        def test(x, y, z):
            return assert_allclose(x.frac_diff_non_gaps(y), z)

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
        transitions = dict.fromkeys(
            [
                ("A", "A"),
                ("A", "G"),
                ("G", "A"),
                ("G", "G"),
                ("U", "U"),
                ("U", "C"),
                ("C", "U"),
                ("C", "C"),
            ]
        )

        s1 = self.RNA("UCAGGCAA")
        s2 = self.RNA("CCAAAUGC")
        s3 = self.RNA("GGGGGGGG")
        e = self.RNA("")

        def test(x, y, z):
            return assert_allclose(x.frac_similar(y, transitions), z)

        test(e, e, 0)
        test(s1, e, 0)
        test(s1, s1, 1)
        test(s1, s2, 7.0 / 8)
        test(s1, s3, 5.0 / 8)
        test(s2, s3, 4.0 / 8)

    def test_with_termini_unknown(self):
        """with_termini_unknown should reset termini to unknown char"""
        s1 = self.RNA("-?--AC--?-")
        s2 = self.RNA("AC")
        self.assertEqual(s1.with_termini_unknown(), "????AC????")
        self.assertEqual(s2.with_termini_unknown(), "AC")

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

    def test_replace(self):
        """replace should convert oldchars to new returning same class"""
        seq = self.SEQ("ACC--GT")
        got = seq.replace("-", "N")
        self.assertEqual(str(got), "ACCNNGT")
        self.assertTrue(isinstance(got, self.SEQ))

    def test_counts(self):
        """count motifs of different sizes, +/- ambiguities"""
        # test DNA seq
        orig = "AACCGGTTAN-T"
        seq = self.DNA(orig)
        # no gaps, no ambiguities
        got = seq.counts()
        expect = dict(A=3, C=2, G=2, T=3)
        self.assertEqual(dict(got), expect)
        # gaps allowed
        got = seq.counts(allow_gap=True)
        expect = dict(A=3, C=2, G=2, T=3)
        expect.update({"-": 1})
        self.assertEqual(dict(got), expect)
        # ambig allowed
        got = seq.counts(include_ambiguity=True)
        expect = dict(A=3, C=2, G=2, T=3, N=1)
        self.assertEqual(dict(got), expect)
        # ambig and gap allowed
        got = seq.counts(include_ambiguity=True, allow_gap=True)
        expect = dict(A=3, C=2, G=2, T=3, N=1)
        expect.update({"-": 1})
        self.assertEqual(dict(got), expect)

        # test DNA seq motif length of 2
        got = seq.counts(motif_length=2)
        expect = dict(AA=1, CC=1, GG=1, TT=1)
        self.assertEqual(dict(got), expect)
        # gap allowed
        got = seq.counts(motif_length=2, allow_gap=True)
        expect = dict(AA=1, CC=1, GG=1, TT=1)
        expect.update({"-T": 1})
        # ambig allowed
        got = seq.counts(motif_length=2, include_ambiguity=True)
        expect = dict(AA=1, CC=1, GG=1, TT=1, AN=1)
        self.assertEqual(dict(got), expect)
        # ambig and gap allowed
        got = seq.counts(motif_length=2, include_ambiguity=True, allow_gap=True)
        expect = dict(AA=1, CC=1, GG=1, TT=1, AN=1)
        expect.update({"-T": 1})
        self.assertEqual(dict(got), expect)

        # test base -- no concept of ambiguity, but understands gap
        orig = "AACCGGTTAN-T"
        seq = self.SEQ(orig)
        got = seq.counts()
        expect = dict(A=3, C=2, G=2, T=3, N=1)
        self.assertEqual(dict(got), expect)

        # handle '?'
        orig = "AACCGGTTAN-T?"
        seq = self.DNA(orig)
        got = seq.counts()
        expect = dict(A=3, C=2, G=2, T=3)
        self.assertEqual(dict(got), expect)
        got = seq.counts(allow_gap=True, include_ambiguity=True)
        expect.update({"-": 1, "N": 1, "?": 1})
        self.assertEqual(dict(got), expect)

    def test_strand_symmetry(self):
        """correctly compute test of strand symmetry"""
        from cogent3 import get_moltype
        from cogent3.core.alignment import Aligned

        seq = DnaSequence("ACGGCTGAAGCGCTCCGGGTTTAAAACG")
        ssym = seq.strand_symmetry(motif_length=1)
        assert_allclose(ssym.observed.array, [[7, 5], [7, 9]])
        assert_allclose(ssym.expected.array, [[6, 6], [8, 8]])

        # RNA too
        seq = seq.to_rna()
        ssym = seq.strand_symmetry(motif_length=1)
        assert_allclose(ssym.observed.array, [[7, 5], [7, 9]])

        # Aligned
        seq = DnaSequence("ACGGCTGAAGCGCTCCGGGTTTAAAACG")
        m, s = seq.parse_out_gaps()
        seq = Aligned(m, s)
        ssym = seq.strand_symmetry(motif_length=1)
        assert_allclose(ssym.observed.array, [[7, 5], [7, 9]])

        with self.assertRaises(TypeError):
            text = get_moltype("text")
            m, s = text.make_seq("ACGGCTGAAGCGCTCCGGGTTTAAAACG").parse_out_gaps()
            s.strand_symmetry(motif_length=1)

        # with motif_length=2
        seq = DnaSequence("AC GG CT GA AG CG CT CC GG GT TT AA AA CG".replace(" ", ""))
        ssym = seq.strand_symmetry(motif_length=2)
        self.assertLessEqual(len(ssym.observed.keys()), 8)
        assert_allclose(ssym.observed["AA"].to_array(), [2, 1])
        assert_allclose(ssym.observed["CC"].to_array(), [1, 2])

    def test_is_annotated(self):
        """is_annotated operates correctly"""
        s = self.SEQ("ACGGCTGAAGCGCTCCGGGTTTAAAACG")
        if hasattr(s, "annotation_db"):
            self.assertFalse(s.is_annotated())
            _ = s.add_feature(biotype="gene", name="blah", spans=[(0, 10)])
            self.assertTrue(s.is_annotated())
        else:
            with self.assertRaises(AttributeError):
                s.is_annotated()

    def test_to_html(self):
        """produce correct html formatted text"""
        seq = DnaSequence("ACGGTGGGGGGGGG")
        got = seq.to_html(wrap=50)
        # ensure balanced tags are in the txt
        for tag in ["<style>", "</style>", "<div", "</div>", "<table>", "</table>"]:
            self.assertTrue(tag in got)

        seq_row = (
            '<tr><td class="label">None</td>'
            '<td><span class="A_dna">A</span>'
            '<span class="C_dna">C</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="T_dna">T</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span>'
            '<span class="G_dna">G</span></td></tr>'
        )

        self.assertTrue(seq_row in got)

    def test_repr_html(self):
        """correctly uses set_repr and the environment variable settings"""
        token = 'class="label"'
        seq = self.SEQ("AAAAA")

        orig = [l for l in seq._repr_html_().splitlines() if token in l][0]
        orig_num = len(re.findall(r"\bA\b", orig))
        self.assertEqual(orig_num, 5)

        # using environment variable
        env_name = "COGENT3_ALIGNMENT_REPR_POLICY"
        os.environ[env_name] = "num_pos=2"
        got = [l for l in seq._repr_html_().splitlines() if token in l][0]
        got_num = len(re.findall(r"\bA\b", got))
        self.assertEqual(got_num, 2)
        os.environ.pop(env_name, None)

    def test_add(self):
        """Test for the add method within sequence"""

        even = "TCAGAT"
        odd = even + "AAA"
        original_sequence = self.SEQ(even, name="even")
        duplicate_sequence = self.SEQ(even, name="even")
        name_only_duplicate = self.SEQ(even, name="odd")
        different_sequence = self.SEQ(odd, name="odd")

        added_duplicates = original_sequence + duplicate_sequence
        added_name_only_duplicate = original_sequence + name_only_duplicate
        different_sequences = original_sequence + different_sequence

        self.assertIsNone(different_sequences.name)
        self.assertIsNotNone(added_duplicates.name)
        self.assertIsNotNone(added_name_only_duplicate)

        self.assertEqual(original_sequence.name, added_duplicates.name)
        self.assertNotEqual(original_sequence.name, added_name_only_duplicate.name)
        self.assertNotEqual(original_sequence.name, different_sequences.name)

    def test_add2(self):
        """name property correctly handled in sequence add"""
        a1 = self.SEQ("AAA", name="1")
        a2 = self.SEQ("CC", name="1")
        a = a1 + a2
        self.assertEqual(a.name, "1")
        self.assertEqual(a, "AAACC")

        b = self.SEQ("GGGG", name="2")
        self._check_mix_add(a1, b)
        c = self.SEQ("TT")
        self._check_mix_add(a1, c)

        e = "AA"
        be = b + e
        self.assertIsNone(be.name)
        self.assertEqual(be, str(b) + e)

    def _check_mix_add(self, s1, s2):
        s1s2 = s1 + s2
        s2s1 = s2 + s1
        self.assertIsNone(s1s2.name)
        self.assertIsNone(s2s1.name)
        self.assertEqual(s1s2, str(s1) + str(s2))
        self.assertEqual(s2s1, str(s2) + str(s1))


class SequenceSubclassTests(TestCase):
    SequenceClass = Sequence

    def test_DnaSequence(self):
        """DnaSequence should behave as expected"""
        x = DnaSequence("tcag")
        # note: no longer preserves case
        self.assertEqual(x, "TCAG")

        x = DnaSequence("aaa") + DnaSequence("ccc")
        # note: doesn't preserve case
        self.assertEqual(x, "AAACCC")
        assert x.moltype is DNA
        self.assertRaises(AlphabetError, x.__add__, "z")
        self.assertEqual(DnaSequence("TTTAc").rc(), "GTAAA")

    def test_get_type(self):
        """returns moltype label"""
        for moltype in ("text", "dna", "bytes"):
            seq = get_moltype(moltype).make_seq("ARCGT")
            self.assertEqual(seq.get_type(), moltype)

    def test_resolved_ambiguities(self):
        seq = get_moltype("dna").make_seq("ARC")
        got = seq.resolved_ambiguities()
        self.assertEqual(got, [("A",), ("A", "G"), ("C",)])

        seq = get_moltype("dna").make_seq("AGC")
        got = seq.resolved_ambiguities()
        self.assertEqual(got, [("A",), ("G",), ("C",)])

    def test_iter_kmers(self):
        """correctly yield all k-mers"""
        from typing import Generator

        orig = "TCAGGA"
        r = self.SequenceClass(orig)
        self.assertIsInstance(r.iter_kmers(k=1), Generator)

        for k in range(1, 7):
            expect = [str(orig[i : i + k]) for i in range(len(orig) - k + 1)]
            got = list(r.iter_kmers(k))
            self.assertEqual(got, expect)

        orig = ""
        r = self.SequenceClass(orig)
        self.assertIsInstance(r.iter_kmers(k=1), Generator)
        got = list(r.iter_kmers(k=1))
        self.assertEqual(got, [])

    def test_iter_kmers_handles_invalid(self):
        """raise exceptions on invalid input to iter_kmers"""
        orig = "TCAGGA"
        r = self.SequenceClass(orig)
        for k in (0, -1, 1.1):
            with self.assertRaises(ValueError):
                _ = list(r.iter_kmers(k))

    def test_get_kmers(self):
        """returns a list of k-mers"""
        orig = "TCAGGA"
        r = self.SequenceClass(orig)

        for k in range(1, 7):
            expect = [str(orig[i : i + k]) for i in range(len(orig) - k + 1)]
            got = r.get_kmers(k)
            self.assertEqual(got, expect)


# TODO move methods of this class onto the single class that inherits from it!
class ModelSequenceTests(object):
    """base class for tests of specific ArraySequence objects."""

    SequenceClass = None  # override in derived classes

    def test_to_fasta(self):
        """Sequence to_fasta() should return Fasta-format string"""
        even = "TCAGAT"
        odd = even + "AAA"
        even_dna = self.SequenceClass(even, name="even")
        odd_dna = self.SequenceClass(odd, name="odd")
        self.assertEqual(even_dna.to_fasta(), ">even\nTCAGAT\n")
        # set line wrap to small number so we can test that it works
        self.assertEqual(even_dna.to_fasta(block_size=2), ">even\nTC\nAG\nAT\n")
        self.assertEqual(odd_dna.to_fasta(block_size=2), ">odd\nTC\nAG\nAT\nAA\nA\n")
        # check that changing the linewrap again works
        self.assertEqual(even_dna.to_fasta(block_size=4), ">even\nTCAG\nAT\n")

    def test_to_phylip(self):
        """Sequence to_phylip() should return one-line phylip string"""
        s = self.SequenceClass("ACG", name="xyz")
        self.assertEqual(s.to_phylip(), "xyz" + " " * 27 + "ACG")


class DnaSequenceTests(ModelSequenceTests, TestCase):
    class SequenceClass(ArrayNucleicAcidSequence):
        alphabet = DNA.alphabets.base

    def test_init(self):
        """Sequence should do round-trip from string"""
        orig = ""
        r = self.SequenceClass(orig)
        self.assertEqual(str(r), orig)

        orig = "TCAGGA"
        r = self.SequenceClass(orig)
        assert_equal(r._data, array([0, 1, 2, 3, 3, 2]))
        self.assertEqual(str(r), orig)


class CodonSequenceTests(SequenceTests, TestCase):
    class SequenceClass(ArrayCodonSequence):
        alphabet = DNA.alphabets.base**3

    def test_init(self):
        """Sequence should do round-trip from string"""
        orig = ""
        r = self.SequenceClass(orig)
        self.assertEqual(str(r), orig)

        orig = "TCAGGA"
        r = self.SequenceClass(orig)
        assert_equal(r._data, array([6, 62]))
        self.assertEqual(str(r), orig)


class DnaSequenceGapTests(TestCase):
    """Tests of gapped DNA sequences."""

    class SequenceClass(ArrayNucleicAcidSequence):
        alphabet = DNA.alphabets.gapped
        gap = "-"

    def test_init(self):
        """gapped sequence should init ok"""
        orig = "TC---"
        seq = self.SequenceClass(orig)
        self.assertEqual(str(seq), orig)

    def test_gaps(self):
        """gapped sequence gaps() should return correct array"""
        sc = self.SequenceClass
        assert_equal(sc("TC").gaps(), array([0, 0]))
        assert_equal(sc("T-").gaps(), array([0, 1]))

    def test_degap(self):
        """gapped sequence degap() should return correct array"""
        sc = self.SequenceClass
        self.assertEqual(sc("T-").degap(), sc("T"))

    def test_nongaps(self):
        """gapped sequence nongaps() should return correct array"""
        sc = self.SequenceClass
        assert_equal(sc("TC").nongaps(), array([1, 1]))
        assert_equal(sc("T-").nongaps(), array([1, 0]))

    def test_regap(self):
        """gapped sequence regap() should return correct sequence"""
        sc = self.SequenceClass
        self.assertEqual(str(sc("TC").regap(sc("A---A-"))), "T---C-")

    def test_degap_name(self):
        """degap preserves name attribute"""
        # todo this should work for any seq class, but is not
        seq = DNA.make_seq("ACG---T", "blah")
        got = seq.degap()
        self.assertEqual(str(got), "ACGT")
        self.assertEqual(got.name, "blah")


class SequenceIntegrationTests(TestCase):
    """Should be able to convert regular to model sequences, and back"""

    def test_regular_to_model(self):
        """Regular sequence should convert to model sequence"""
        r = RNA.make_seq("AAA", name="x")
        s = RNA.make_array_seq(r)
        self.assertEqual(str(s), "AAA")
        self.assertEqual(s.moltype, RNA)
        self.assertEqual(s.name, "x")

    def test_model_to_regular(self):
        """Model sequence should convert to regular sequence"""
        r = RNA.make_array_seq("AAA", name="x")
        s = RNA.make_seq(r)
        self.assertEqual(str(s), "AAA")
        self.assertEqual(s.moltype, RNA)
        self.assertEqual(s.name, "x")

    def test_regular_to_regular(self):
        """Regular sequence should convert to regular sequence"""
        r = RNA.make_seq("AAA", name="x")
        s = RNA.make_seq(r)
        self.assertEqual(str(s), "AAA")
        self.assertEqual(s.moltype, RNA)
        self.assertEqual(s.name, "x")

    def test_model_to_model(self):
        """Model sequence should convert to model sequence"""
        r = RNA.make_array_seq("AAA", name="x")
        s = RNA.make_array_seq(r)
        self.assertEqual(str(s), "AAA")
        self.assertEqual(s.moltype, RNA)
        self.assertEqual(s.name, "x")

    def test_ModelDnaCodonSequence(self):
        """ArrayDnaCodonSequence should behave as expected"""
        d = ArrayDnaCodonSequence("UUUCGU")
        self.assertEqual(str(d), "TTTCGT")
        assert_equal(d._data, array([0, 28]))
        self.assertEqual(str(d.to_rna()), "UUUCGU")
        self.assertEqual(str(d.to_dna()), "TTTCGT")

    def test_ModelRnaCodonSequence(self):
        """ArrayRnaCodonSequence should behave as expected"""
        r = ArrayRnaCodonSequence("UUUCGU")
        self.assertEqual(str(r), "UUUCGU")
        assert_equal(r._data, array([0, 28]))
        self.assertEqual(str(r.to_rna()), "UUUCGU")
        self.assertEqual(str(r.to_dna()), "TTTCGT")


class ModelSequenceTests(SequenceTests):
    """Tests of the ArraySequence class's inheritance of SequenceI."""

    SEQ = ArraySequence
    RNA = ArrayRnaSequence
    DNA = ArrayDnaSequence
    PROT = ArrayProteinSequence
    AB = ABSequence

    def test_distance_indices(self):
        """ArraySequence distance should work with function of indices"""
        s1 = self.RNA("AUGC")
        s2 = self.RNA("AAGC")

        def f(x, y):
            if x == 2 or y == 2:
                return 10
            return 0

        self.assertEqual(s1.distance(s2, f, use_indices=True), 20)

    def test_strip_bad(self):
        """Sequence strip_bad should remove any non-base, non-gap chars"""
        # have to turn off check to get bad data in; no longer preserves case
        r = self.RNA("UCAGRYU")
        r._data[0] = 31
        r._data[2] = 55
        self.assertEqual(r.strip_bad(), "CGRYU")

    def test_strip_bad_and_gaps(self):
        """Sequence strip_bad_and_gaps should remove gaps and bad chars"""
        # have to turn off check to get bad data in; no longer preserves case
        r = self.RNA("ACG--GRN?")
        self.assertEqual(r.strip_bad_and_gaps(), "ACGGRN")
        r._data[0] = 99
        self.assertEqual(r.strip_bad_and_gaps(), "CGGRN")

    def test_gap_array(self):
        """Sequence gap_array should return array of gaps"""
        r = self.RNA("-?A-?NRY-")
        v = r.gap_array()
        assert_equal(v, array([1, 1, 0, 1, 1, 0, 0, 0, 1]))
        r = self.RNA("AC")
        v = r.gap_array()
        assert_equal(v, array([0, 0]))
        r = self.RNA("-?")
        v = r.gap_array()
        assert_equal(v, array([1, 1]))

    def test_gap_indices(self):
        """Sequence gap_indices should return positions of gaps"""
        r = self.RNA("-?A-?NRY-")
        v = r.gap_indices()
        assert_equal(v, array([0, 1, 3, 4, 8]))
        r = self.RNA("AC")
        v = r.gap_indices()
        assert_equal(v, array([]))  # note: always returns array
        r = self.RNA("-?")
        v = r.gap_indices()
        assert_equal(v, array([0, 1]))

    def test_count_ab(self):
        """abseq array seq should count characters"""
        AB = get_moltype("ab")
        seq = AB.make_array_seq("aaba-", alphabet=AB.alphabet.with_gap_motif())
        c = seq.counts()
        self.assertEqual(c.to_dict(), {"a": 3, "b": 1})
        c = seq.counts(allow_gap=True)
        self.assertEqual(c.to_dict(), {"a": 3, "b": 1, "-": 1})


@pytest.mark.parametrize("seq,rc", (("ATGTTT", False), ("AAACAT", True)))
def test_translation(seq, rc):
    seq = DNA.make_seq(seq)
    if rc:
        seq = seq.rc()
    assert str(seq) == "ATGTTT"
    aa = seq.get_translation()
    assert str(aa) == "MF"


def test_get_translation_include_stop():
    s = DNA.make_seq("ATTTAACTT", name="s1")
    aa = s.get_translation(include_stop=True)
    assert str(aa) == "I*L"


def test_get_translation_trim_stop():
    s = DNA.make_seq("ATTTCCTGA", name="s1")
    aa = s.get_translation(trim_stop=True)
    assert str(aa) == "IS"
    # no effect on internal stops
    s = DNA.make_seq("ATTTAACTT", name="s1")
    aa = s.get_translation(include_stop=True, trim_stop=True)
    assert str(aa) == "I*L"


@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
def test_seqview_initialisation(start, stop, step):
    """Initialising a SeqView should work with range of provided values"""
    seq_data = "0123456789"
    got = SeqView(seq_data, start=start, stop=stop, step=step)
    expected = seq_data[start:stop:step]
    assert got.value == expected


def test_seqview_invalid_step():
    "Testing that SeqView raises Value error when initialised with step of 0"
    with pytest.raises(ValueError):
        _ = SeqView(seq="0123456789", step=0)


@pytest.mark.parametrize("index", (-10, -5, 0, 5, 9))  # -10 and 9 are boundary
def test_seqview_index(index):
    """SeqView with default values can be sliced with a single index, when within the length of the sequence"""
    seq_data = "0123456789"
    sv = SeqView(seq_data)
    got = sv[index]
    expected = seq_data[index]
    assert got.value == expected
    assert len(got) == 1


def test_seqview_index_null():
    "Indexing a SeqView of length 0 should return an IndexError"
    sv = SeqView("")
    with pytest.raises(IndexError):
        _ = sv[0]


def test_seqview_step_0():
    "Initialising or slicing a SeqView with a step of 0 should return an IndexError"
    sv = SeqView("0123456789")
    with pytest.raises(ValueError):
        _ = sv[::0]
    with pytest.raises(ValueError):
        _ = SeqView("0123456789", step=0)


@pytest.mark.parametrize("start", (0, 2, 4))
def test_seqview_invalid_index(start):
    "indexing out of bounds with a forward step should raise an IndexError"
    seq = "0123456789"
    length = abs(start - len(seq))
    pos_boundary_index = length
    neg_boundary_index = -length - 1

    sv = SeqView(seq=seq, start=start)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("start", (0, 2, 4))
def test_seqview_invalid_index_positive_step_gt_1(start):
    "boundary condition for indexing out of bounds with a forward step greater than 1"
    seq = "0123456789"
    step = 2
    length = abs((start - len(seq)) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sv = SeqView(seq=seq, start=start, step=step)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("stop", (0, 2, -11))
def test_seqview_invalid_index_reverse_step(stop):
    "boundary condition for indexing out of bounds with a reverse step"
    seq = "0123456789"
    step = -1
    start = len(seq)
    length = abs((start - stop) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sv = SeqView(seq=seq, start=start, stop=stop, step=step)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


@pytest.mark.parametrize("stop", (0, 2, -6))
def test_seqview_invalid_index_reverse_step_gt_1(stop):
    "boundary condition for indexing out of bounds with a reverse step less than -1"
    seq = "0123456789"
    step = -2
    start = len(seq)
    length = abs((start - stop) // step)
    neg_boundary_index = -length - 1
    pos_boundary_index = length

    sv = SeqView(seq=seq, start=start, stop=stop, step=step)
    with pytest.raises(IndexError):
        _ = sv[pos_boundary_index]
    with pytest.raises(IndexError):
        _ = sv[neg_boundary_index]


def test_seqview_slice_null():
    sv = SeqView("")
    assert len(sv) == 0
    got = sv[2:]
    assert len(got) == 0


def test_seqview_start_out_of_bounds():
    "boundary condition for start index out of bounds"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, 1
    boundary = abs((init_start - init_stop) // init_step)
    sv = SeqView(seq=seq, start=init_start, stop=init_stop, step=init_step)
    got = sv[boundary::].value
    assert got == ""


def test_seqview_start_out_of_bounds_step_gt_1():
    "boundary condition for start index out of bounds with step greater than 1"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, 2
    boundary = abs((init_start - init_stop) // init_step)
    sv = SeqView(seq=seq, start=init_start, stop=init_stop, step=init_step)
    got = sv[boundary::].value
    assert got == ""


def test_seqview_start_out_of_bounds_reverse_step():
    "boundary condition for start index out of bounds with reverse step"
    seq = "0123456789"
    init_start, init_stop, init_step = 2, 10, -2
    boundary_pos = abs((init_start - init_stop) // init_step)
    boundary_neg = -abs((init_start - init_stop) // init_step) - 1

    sv = SeqView(seq=seq, start=init_start, stop=init_stop, step=init_step)

    assert sv[boundary_pos::].value == ""
    assert sv[boundary_neg::].value == ""


@pytest.mark.parametrize(
    "simple_slices",
    (
        slice(None, None, 1),
        slice(None, 3, None),
        slice(1, None, None),
        slice(1, 3, None),
        slice(None, None, None),
    ),
)
def test_seqview_defaults(simple_slices):
    """SeqView should accept slices with all combinations of default parameters"""
    seq = "0123456789"
    got = SeqView(seq)[simple_slices]
    expected = seq[simple_slices]
    assert got.value == expected


@pytest.mark.parametrize("index", (-8, -5, 0, 5, 8))
@pytest.mark.parametrize(
    "simple_slices",
    (
        slice(None, None, 1),
        slice(None, 10, None),
        slice(1, None, None),
        slice(1, 10, None),
        slice(1, 10, 1),
        slice(None, None, None),
    ),
)
def test_seqview_sliced_index(index, simple_slices):
    """SeqView that has been sliced with default parameters, can then be indexed"""
    seq = "0123456789"
    sv = SeqView(seq)
    got = sv[simple_slices][index]
    expected = seq[simple_slices][index]
    assert got.value == expected


@pytest.mark.parametrize("first_step", (1, 2, -1, -2))
@pytest.mark.parametrize("second_step", (1, 2, -1, -2))
def test_seqview_reverse_slice(first_step, second_step):
    """subsequent slices may reverse the previous slice"""
    seq = "0123456789"
    sv = SeqView(seq=seq, step=first_step)
    got = sv[::second_step]
    expected = seq[::first_step][::second_step]
    assert got.value == expected


@pytest.mark.parametrize("seq", ("0123456789", "01234567890"))
@pytest.mark.parametrize("index", (-10, -4, 0, 6, 10))
@pytest.mark.parametrize("start", (None, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 9, -10, -11))
@pytest.mark.parametrize("step", (-1, -2))
def test_seqview_rev_sliced_index(index, start, stop, step, seq):
    """SeqView that has been reverse sliced, can then be sliced with a single index"""
    seq_data = seq
    try:  # if python slicing raises an index error, we expect SeqView to also throw error
        expected = seq_data[start:stop:step][index]
    except IndexError:
        with pytest.raises(IndexError):
            _ = SeqView(seq=seq_data, start=start, stop=stop, step=step)[index].value
    else:  # if no index error, SeqView should match python slicing
        got = SeqView(seq=seq_data, start=start, stop=stop, step=step)[index].value
        assert got == expected


@pytest.mark.parametrize("seq", ("0123456789", "012345678"))
@pytest.mark.parametrize("start", (None, 0, 1, 9, -1, -10))
@pytest.mark.parametrize("stop", (None, 0, 10, -7, -11))
@pytest.mark.parametrize("step", (1, 2, -1, -2))
def test_seqview_init_with_negatives(seq, start, stop, step):
    "SeqView initialisation should handle any combination of positive and negative slices"
    got = SeqView(seq, start=start, stop=stop, step=step)
    expected = seq[start:stop:step]
    assert got.value == expected


@pytest.mark.parametrize("seq", ("0123456789", "012345678"))
@pytest.mark.parametrize("start", (None, 0, 1, 9, -1, -10))
@pytest.mark.parametrize("stop", (None, 0, 10, -7, -11))
@pytest.mark.parametrize("step", (1, 2, -1, -2))
def test_seqview_slice_with_negatives(seq, start, stop, step):
    """SeqView should handle any combination of positive and negative slices"""
    sv = SeqView(seq)
    got = sv[start:stop:step]
    expected = seq[start:stop:step]
    assert got.value == expected


@pytest.mark.parametrize("start", (None, 0, 2))
@pytest.mark.parametrize("stop", (None, 5, 7, 10))
@pytest.mark.parametrize("step", (1, 2))
@pytest.mark.parametrize("start_2", (None, 0, 1, 2))
@pytest.mark.parametrize("stop_2", (None, 2, 4, 10))
@pytest.mark.parametrize("step_2", (1, 2))
def test_subsequent_slice_forward(start, stop, step, start_2, stop_2, step_2):
    """SeqView should handle subsequent forward slice"""
    seq = "0123456789"
    sv = SeqView(seq=seq)
    got = sv[start:stop:step][start_2:stop_2:step_2]
    expected = seq[start:stop:step][start_2:stop_2:step_2]
    assert got.value == expected
    assert len(got) == len(expected)


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first stop -ve
        (slice(None, -3, None), slice(None, None, None)),
        # second stop -ve
        (slice(None, None, None), slice(None, -1, None)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -3, None), slice(None, -5, None)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -5, None), slice(None, -3, None)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, None), slice(None, -8, None)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, None), slice(None, -3, None)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, None), slice(None, 7, None)),
        # first stop -ve, second stop +ve, second slice OUTSIDE first
        (slice(None, -6, None), slice(None, 7, None)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 6, None), slice(None, -2, None)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, None), slice(None, -7, None)),
        # WITH FIRST STEP > 1
        # first stop -ve
        (slice(None, -3, 2), slice(None, None, None)),
        # second stop -ve
        (slice(None, None, 2), slice(None, -1, None)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -1, 2), slice(None, -3, None)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -3, 2), slice(None, -2, None)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, 2), slice(None, -8, None)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, 2), slice(None, -3, None)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, 2), slice(None, 3, None)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -6, 2), slice(None, 7, None)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 6, 2), slice(None, -2, None)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, 2), slice(None, -7, None)),
        # WITH SECOND STEP > 1
        # first stop -ve
        (slice(None, -3, None), slice(None, None, 3)),
        # second stop -ve
        (slice(None, None, None), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -2, None), slice(None, -4, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -4, None), slice(None, -3, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, None), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, None), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, None), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -6, None), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 9, None), slice(None, -2, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, None), slice(None, -7, 3)),
        # WITH BOTH STEP > 1
        # first stop -ve
        (slice(None, -3, 2), slice(None, None, 3)),
        # second stop -ve
        (slice(None, None, 2), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(None, -1, 3), slice(None, -2, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(None, -2, 2), slice(None, -1, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(None, -3, 3), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(None, -8, 2), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(None, -2, 3), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(None, -3, 3), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(None, 9, 2), slice(None, -1, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(None, 6, 2), slice(None, -7, 3)),
        # NON-ZERO START
        # first stop -ve
        (slice(1, -3, 2), slice(None, None, 3)),
        # second stop -ve
        (slice(1, None, 2), slice(None, -1, 3)),
        # both stop -ve, (first > second), second slice WITHIN first
        (slice(1, -1, 3), slice(None, -2, 2)),
        # both stop -ve, (first < second), second slice WITHIN first
        (slice(1, -2, 2), slice(None, -1, 2)),
        # both stop -ve, (first > second), second slice OUTSIDE first
        (slice(1, -3, 3), slice(None, -8, 2)),
        # both stop -ve, (first < second), second slice OUTSIDE first
        (slice(1, -8, 2), slice(None, -3, 2)),
        # first stop -ve, second stop +ve, second slice WITHIN first
        (slice(1, -2, 3), slice(None, 7, 2)),
        # first stop -ve, second stop +ve, second slice OVERLAP first
        (slice(1, -3, 3), slice(None, 7, 2)),
        # first stop +ve, second stop -ve, second slice WITHIN first
        (slice(1, 10, 2), slice(None, -1, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
        (slice(1, 6, 2), slice(None, -7, 3)),
        # first stop +ve, second stop -ve, second slice OUTSIDE first
    ),
)
def test_subsequent_slice_neg_stop(slice_1, slice_2):
    """SeqView should handle subsequence slices with >=1 negative stop values,
    subsequent slices may overlap or be within previous slices
    """
    seq_data = "abcdefghijk"
    sv = SeqView(seq_data)
    assert sv[slice_1][slice_2].value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first start -ve
        (slice(-6, None, None), slice(None, None, None)),
        # second start -ve
        (slice(None, None, None), slice(-6, None, None)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-6, None, None), slice(-4, None, None)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-4, None, None), slice(-6, None, None)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-8, None, None), slice(2, None, None)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, None), slice(7, None, None)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, None), slice(-7, None, None)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(5, None, None), slice(-6, None, None)),
        # WITH FIRST STEP > 1
        # first start -ve
        (slice(-6, None, 2), slice(None, None, None)),
        # second start -ve
        (slice(None, None, 2), slice(-6, None, None)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-8, None, 2), slice(-6, None, None)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, 2), slice(-9, None, None)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, 2), slice(2, None, None)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, 2), slice(7, None, None)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, 2), slice(-7, None, None)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, 2), slice(-9, None, None)),
        # WITH SECOND STEP > 1
        # first start -ve
        (slice(-6, None, None), slice(None, None, 2)),
        # second start -ve
        (slice(None, None, None), slice(-6, None, 2)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-8, None, None), slice(-6, None, 2)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, None), slice(-9, None, 2)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, None), slice(2, None, 2)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, None), slice(7, None, 2)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, None), slice(-7, None, 2)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, None), slice(-9, None, 2)),
        # WITH BOTH STEP > 1
        # first start -ve
        (slice(-6, None, 3), slice(None, None, 2)),
        # second start -ve
        (slice(None, None, 3), slice(-6, None, 2)),
        # both start -ve, (first < second), second slice WITHIN first
        (slice(-9, None, 3), slice(-7, None, 2)),
        # both start -ve, (first > second), second slice OUTSIDE first
        (slice(-7, None, 3), slice(-9, None, 2)),
        # first start -ve, second start +ve, second slice WITHIN first
        (slice(-9, None, 3), slice(2, None, 2)),
        # first start -ve, second start +ve, second slice OUTSIDE first
        (slice(-6, None, 2), slice(7, None, 2)),
        # first start +ve, second start -ve, second slice WITHIN first
        (slice(2, None, 3), slice(-7, None, 2)),
        # first start +ve, second start -ve, second slice OUTSIDE first
        (slice(3, None, 3), slice(-9, None, 2)),
        (slice(-9, 7, 3), slice(-2, None, None)),
    ),
)
def test_subsequent_slice_neg_start(slice_1, slice_2):
    """SeqView should handle subsequence slices with >=1 negative start values,
    subsequent slices may or may not overlap or be within previous slices
    """
    seq_data = "abcdefghijk"
    sv = SeqView(seq_data)
    assert sv[slice_1][slice_2].value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "slice_1, slice_2",
    (
        # WITH DEFAULTS
        # first step -ve
        (slice(None, None, -1), slice(None, None, None)),
        # second step -ve
        (slice(None, None, None), slice(None, None, -1)),
        # both step -ve, start/stop -ve, second slice WITHIN first
        (slice(-1, -11, -2), slice(-1, -5, -3)),
        # both step -ve, start/stop -ve, second slice OUTSIDE first
        (slice(-1, -11, -2), slice(-1, -11, -3)),
        # both step -ve, start/stop +ve, second slice WITHIN first
        (slice(10, 0, -2), slice(5, 0, -3)),
        # both step -ve, start/stop +ve, second slice OUTSIDE first
        (slice(10, 0, -2), slice(10, 0, -3)),
        # first step -ve, second step +ve, second slice WITHIN first
        (slice(10, 0, -2), slice(1, 5, 2)),
        # first step -ve, second step +ve, second slice OUTSIDE first
        (slice(10, 0, -2), slice(0, 10, 2)),
        # first step +ve, second step -ve, second slice WITHIN first
        (slice(0, 10, 2), slice(4, 0, -2)),
        # first step +ve, second step -ve, second slice OUTSIDE first
        (slice(0, 10, 3), slice(10, 0, -2)),
        # first step -ve, second step +ve, second start/stop +ve
        (slice(10, 1, -1), slice(-8, 11, 2)),
        # first step -ve, second step +ve, second start/stop +ve
        (slice(10, 1, -1), slice(-19, 0, -2)),
    ),
)
def test_subsequent_slice_neg_step(slice_1, slice_2):
    """SeqView should handle subsequence slices with negative step values,
    subsequent slices may overlap or be within previous slices
    """
    seq_data = "0123456789"
    sv = SeqView(seq_data)
    assert sv[slice_1][slice_2].value == seq_data[slice_1][slice_2]


@pytest.mark.parametrize(
    "sub_slices_triple",
    (
        (slice(None, None, None), slice(None, None, None), slice(None, None, None)),
        (slice(1, 9, 1), slice(2, 8, 1), slice(3, 7, 1)),
        (slice(1, 9, 1), slice(2, 8, 1), slice(3, 9, 1)),
        (slice(1, 9, 1), slice(2, 8, 2), slice(3, 7, -3)),
    ),
)
def test_subslice_3(sub_slices_triple):
    """SeqView should handle three subsequent slices"""
    seq_data = "abcdefghijk"
    sv = SeqView(seq_data)
    slice_1, slice_2, slice_3 = sub_slices_triple
    assert sv[slice_1][slice_2][slice_3].value == seq_data[slice_1][slice_2][slice_3]


@pytest.mark.parametrize("start", (0, 2, -1))
@pytest.mark.parametrize("stop", (7, 10, -11))
@pytest.mark.parametrize("step", (1, -2))
@pytest.mark.parametrize("start_2", (0, 2, -8))
@pytest.mark.parametrize("stop_2", (2, 4))
@pytest.mark.parametrize("step_2", (2, -1))
@pytest.mark.parametrize("start_3", (0, 1, -6))
@pytest.mark.parametrize("stop_3", (4, 10, -10))
@pytest.mark.parametrize("step_3", (2, -2))
def test_triple_slice(
    start, stop, step, start_2, stop_2, step_2, start_3, stop_3, step_3
):
    """SeqView should handle subsequent forward slice"""
    seq = "0123456789"
    sv = SeqView(seq=seq)
    got = sv[start:stop:step][start_2:stop_2:step_2][start_3:stop_3:step_3]
    expected = seq[start:stop:step][start_2:stop_2:step_2][start_3:stop_3:step_3]

    assert got.value == expected
    assert len(got) == len(expected)


def test_seqview_replace():
    """SeqView supports replacements of substrings, however overriding the sequence data"""
    seq_data = "abcdefghijk"
    sv = SeqView(seq_data)
    sv_replaced = sv.replace("a", "u")
    assert sv_replaced.value == seq_data.replace("a", "u")
    assert sv_replaced.replace("u", "a").value == seq_data
    assert sv_replaced.seq == seq_data.replace("a", "u")


def test_seqview_remove_gaps():
    """Replacing strings of different lengths should work, although any previous slices will be lost"""
    seq_data = "abc----def"
    sv = SeqView(seq_data)
    sliced = sv[2:4]
    assert sliced.start == 2
    replaced = sliced.replace("-", "")
    assert replaced.value == seq_data.replace("-", "")
    assert replaced.start == 0  # start should now be zero,
    assert replaced.stop == len(seq_data.replace("-", ""))


def test_seqview_repr():
    # Short sequence, defaults
    seq = "ACGT"
    view = SeqView(seq)
    expected = (
        "SeqView(seq='ACGT', start=0, stop=4, step=1, offset=0, seqid=None, seq_len=4)"
    )
    assert repr(view) == expected

    # Long sequence
    seq = "ACGT" * 10
    view = SeqView(seq)
    expected = "SeqView(seq='ACGTACGTAC...TACGT', start=0, stop=40, step=1, offset=0, seqid=None, seq_len=40)"
    assert repr(view) == expected

    # Non-zero start, stop, and step values
    seq = "ACGT" * 10
    view = SeqView(seq, start=5, stop=35, step=2)
    expected = "SeqView(seq='ACGTACGTAC...TACGT', start=5, stop=35, step=2, offset=0, seqid=None, seq_len=40)"
    assert repr(view) == expected

    # offset
    seq = "ACGT"
    view = SeqView(seq, offset=5)
    expected = (
        "SeqView(seq='ACGT', start=0, stop=4, step=1, offset=5, seqid=None, seq_len=4)"
    )
    assert repr(view) == expected

    # seqid
    seq = "ACGT"
    view = SeqView(seq, seqid="seq1")
    expected = "SeqView(seq='ACGT', start=0, stop=4, step=1, offset=0, seqid='seq1', seq_len=4)"
    assert repr(view) == expected


def test_get_kmers_strict():
    orig = "TCAGGAN"
    r = DnaSequence(orig)

    assert r.get_kmers(1, strict=True) == ["T", "C", "A", "G", "G", "A"]
    assert r.get_kmers(2, strict=True) == ["TC", "CA", "AG", "GG", "GA"]
    assert r.get_kmers(3, strict=True) == ["TCA", "CAG", "AGG", "GGA"]
    assert r.get_kmers(4, strict=True) == ["TCAG", "CAGG", "AGGA"]
    assert r.get_kmers(5, strict=True) == ["TCAGG", "CAGGA"]
    assert r.get_kmers(6, strict=True) == ["TCAGGA"]
    assert r.get_kmers(7, strict=True) == []

    assert r.get_kmers(1, strict=False) == ["T", "C", "A", "G", "G", "A", "N"]
    assert r.get_kmers(2, strict=False) == ["TC", "CA", "AG", "GG", "GA", "AN"]
    assert r.get_kmers(3, strict=False) == ["TCA", "CAG", "AGG", "GGA", "GAN"]
    assert r.get_kmers(4, strict=False) == ["TCAG", "CAGG", "AGGA", "GGAN"]
    assert r.get_kmers(5, strict=False) == ["TCAGG", "CAGGA", "AGGAN"]
    assert r.get_kmers(6, strict=False) == ["TCAGGA", "CAGGAN"]
    assert r.get_kmers(7, strict=False) == ["TCAGGAN"]
    assert r.get_kmers(8, strict=False) == []


def test_get_kmers_strict_RNA():
    orig = "UCAGGAN"
    r = RnaSequence(orig)

    assert r.get_kmers(1, strict=True) == ["U", "C", "A", "G", "G", "A"]
    assert r.get_kmers(2, strict=True) == ["UC", "CA", "AG", "GG", "GA"]
    assert r.get_kmers(3, strict=True) == ["UCA", "CAG", "AGG", "GGA"]
    assert r.get_kmers(4, strict=True) == ["UCAG", "CAGG", "AGGA"]
    assert r.get_kmers(5, strict=True) == ["UCAGG", "CAGGA"]
    assert r.get_kmers(6, strict=True) == ["UCAGGA"]
    assert r.get_kmers(7, strict=True) == []

    assert r.get_kmers(1, strict=False) == ["U", "C", "A", "G", "G", "A", "N"]
    assert r.get_kmers(2, strict=False) == ["UC", "CA", "AG", "GG", "GA", "AN"]
    assert r.get_kmers(3, strict=False) == ["UCA", "CAG", "AGG", "GGA", "GAN"]
    assert r.get_kmers(4, strict=False) == ["UCAG", "CAGG", "AGGA", "GGAN"]
    assert r.get_kmers(5, strict=False) == ["UCAGG", "CAGGA", "AGGAN"]
    assert r.get_kmers(6, strict=False) == ["UCAGGA", "CAGGAN"]
    assert r.get_kmers(7, strict=False) == ["UCAGGAN"]
    assert r.get_kmers(8, strict=False) == []


def test_get_kmers_strict_protein():
    orig = "CEFGMNX"
    r = ProteinSequence(orig)

    assert r.get_kmers(1, strict=True) == ["C", "E", "F", "G", "M", "N"]
    assert r.get_kmers(2, strict=True) == ["CE", "EF", "FG", "GM", "MN"]
    assert r.get_kmers(3, strict=True) == ["CEF", "EFG", "FGM", "GMN"]
    assert r.get_kmers(4, strict=True) == ["CEFG", "EFGM", "FGMN"]
    assert r.get_kmers(5, strict=True) == ["CEFGM", "EFGMN"]
    assert r.get_kmers(6, strict=True) == ["CEFGMN"]

    assert r.get_kmers(1, strict=False) == ["C", "E", "F", "G", "M", "N", "X"]
    assert r.get_kmers(2, strict=False) == ["CE", "EF", "FG", "GM", "MN", "NX"]
    assert r.get_kmers(3, strict=False) == ["CEF", "EFG", "FGM", "GMN", "MNX"]
    assert r.get_kmers(4, strict=False) == ["CEFG", "EFGM", "FGMN", "GMNX"]
    assert r.get_kmers(5, strict=False) == ["CEFGM", "EFGMN", "FGMNX"]
    assert r.get_kmers(6, strict=False) == ["CEFGMN", "EFGMNX"]
    assert r.get_kmers(7, strict=False) == ["CEFGMNX"]


def test_get_kmers_strict_DNA_gaps():
    orig = "TCA-GAT"
    r = DnaSequence(orig)

    assert r.get_kmers(1, strict=True) == ["T", "C", "A", "G", "A", "T"]
    assert r.get_kmers(2, strict=True) == ["TC", "CA", "GA", "AT"]
    assert r.get_kmers(3, strict=True) == ["TCA", "GAT"]
    assert r.get_kmers(4, strict=True) == []

    assert r.get_kmers(1, strict=False) == ["T", "C", "A", "-", "G", "A", "T"]
    assert r.get_kmers(2, strict=False) == ["TC", "CA", "A-", "-G", "GA", "AT"]
    assert r.get_kmers(3, strict=False) == ["TCA", "CA-", "A-G", "-GA", "GAT"]
    assert r.get_kmers(4, strict=False) == ["TCA-", "CA-G", "A-GA", "-GAT"]
    assert r.get_kmers(5, strict=False) == ["TCA-G", "CA-GA", "A-GAT"]
    assert r.get_kmers(6, strict=False) == ["TCA-GA", "CA-GAT"]
    assert r.get_kmers(7, strict=False) == ["TCA-GAT"]
    assert r.get_kmers(8, strict=False) == []


def test_get_kmers_strict_RNA_gaps():
    orig = "UCA-GAU"
    r = RnaSequence(orig)

    assert r.get_kmers(1, strict=True) == ["U", "C", "A", "G", "A", "U"]
    assert r.get_kmers(2, strict=True) == ["UC", "CA", "GA", "AU"]
    assert r.get_kmers(3, strict=True) == ["UCA", "GAU"]
    assert r.get_kmers(4, strict=True) == []

    assert r.get_kmers(1, strict=False) == ["U", "C", "A", "-", "G", "A", "U"]
    assert r.get_kmers(2, strict=False) == ["UC", "CA", "A-", "-G", "GA", "AU"]
    assert r.get_kmers(3, strict=False) == ["UCA", "CA-", "A-G", "-GA", "GAU"]
    assert r.get_kmers(4, strict=False) == ["UCA-", "CA-G", "A-GA", "-GAU"]
    assert r.get_kmers(5, strict=False) == ["UCA-G", "CA-GA", "A-GAU"]
    assert r.get_kmers(6, strict=False) == ["UCA-GA", "CA-GAU"]
    assert r.get_kmers(7, strict=False) == ["UCA-GAU"]
    assert r.get_kmers(8, strict=False) == []


def test_get_kmers_strict_protein_gaps():
    orig = "CEF-GMN"
    r = ProteinSequence(orig)

    assert r.get_kmers(1, strict=True) == ["C", "E", "F", "G", "M", "N"]
    assert r.get_kmers(2, strict=True) == ["CE", "EF", "GM", "MN"]
    assert r.get_kmers(3, strict=True) == ["CEF", "GMN"]
    assert r.get_kmers(4, strict=True) == []

    assert r.get_kmers(1, strict=False) == ["C", "E", "F", "-", "G", "M", "N"]
    assert r.get_kmers(2, strict=False) == ["CE", "EF", "F-", "-G", "GM", "MN"]
    assert r.get_kmers(3, strict=False) == ["CEF", "EF-", "F-G", "-GM", "GMN"]
    assert r.get_kmers(4, strict=False) == ["CEF-", "EF-G", "F-GM", "-GMN"]
    assert r.get_kmers(5, strict=False) == ["CEF-G", "EF-GM", "F-GMN"]
    assert r.get_kmers(6, strict=False) == ["CEF-GM", "EF-GMN"]
    assert r.get_kmers(7, strict=False) == ["CEF-GMN"]
    assert r.get_kmers(8, strict=False) == []


def test_get_kmers_strict_DNA_RNA_Protein_allgap():
    orig = "-------"

    r = DnaSequence(orig)
    assert r.get_kmers(1, strict=True) == []
    assert r.get_kmers(1, strict=False) == ["-", "-", "-", "-", "-", "-", "-"]

    r = RnaSequence(orig)
    assert r.get_kmers(1, strict=True) == []
    assert r.get_kmers(1, strict=False) == ["-", "-", "-", "-", "-", "-", "-"]

    r = ProteinSequence(orig)
    assert r.get_kmers(1, strict=True) == []
    assert r.get_kmers(1, strict=False) == ["-", "-", "-", "-", "-", "-", "-"]


def test_get_kmers_strict_DNA_RNA_Protein_mixed_ambiguities():
    r = DnaSequence("NGASTAH")
    assert r.get_kmers(1, strict=True) == ["G", "A", "T", "A"]
    assert r.get_kmers(2, strict=True) == ["GA", "TA"]
    assert r.get_kmers(3, strict=True) == []

    assert r.get_kmers(1, strict=False) == ["N", "G", "A", "S", "T", "A", "H"]
    assert r.get_kmers(2, strict=False) == ["NG", "GA", "AS", "ST", "TA", "AH"]
    assert r.get_kmers(3, strict=False) == ["NGA", "GAS", "AST", "STA", "TAH"]
    assert r.get_kmers(4, strict=False) == ["NGAS", "GAST", "ASTA", "STAH"]
    assert r.get_kmers(5, strict=False) == ["NGAST", "GASTA", "ASTAH"]
    assert r.get_kmers(6, strict=False) == ["NGASTA", "GASTAH"]
    assert r.get_kmers(7, strict=False) == ["NGASTAH"]
    assert r.get_kmers(8, strict=False) == []

    r = RnaSequence("RGAWUAD")
    assert r.get_kmers(1, strict=True) == ["G", "A", "U", "A"]
    assert r.get_kmers(2, strict=True) == ["GA", "UA"]
    assert r.get_kmers(3, strict=True) == []
    assert r.get_kmers(1, strict=False) == ["R", "G", "A", "W", "U", "A", "D"]
    assert r.get_kmers(2, strict=False) == ["RG", "GA", "AW", "WU", "UA", "AD"]
    assert r.get_kmers(3, strict=False) == ["RGA", "GAW", "AWU", "WUA", "UAD"]
    assert r.get_kmers(4, strict=False) == ["RGAW", "GAWU", "AWUA", "WUAD"]
    assert r.get_kmers(5, strict=False) == ["RGAWU", "GAWUA", "AWUAD"]
    assert r.get_kmers(6, strict=False) == ["RGAWUA", "GAWUAD"]
    assert r.get_kmers(7, strict=False) == ["RGAWUAD"]
    assert r.get_kmers(8, strict=False) == []

    r = ProteinSequence("BQMXNRZ")
    assert r.get_kmers(1, strict=True) == ["Q", "M", "N", "R"]
    assert r.get_kmers(2, strict=True) == ["QM", "NR"]
    assert r.get_kmers(3, strict=True) == []
    assert r.get_kmers(1, strict=False) == ["B", "Q", "M", "X", "N", "R", "Z"]
    assert r.get_kmers(2, strict=False) == ["BQ", "QM", "MX", "XN", "NR", "RZ"]
    assert r.get_kmers(3, strict=False) == ["BQM", "QMX", "MXN", "XNR", "NRZ"]
    assert r.get_kmers(4, strict=False) == ["BQMX", "QMXN", "MXNR", "XNRZ"]
    assert r.get_kmers(5, strict=False) == ["BQMXN", "QMXNR", "MXNRZ"]
    assert r.get_kmers(6, strict=False) == ["BQMXNR", "QMXNRZ"]
    assert r.get_kmers(7, strict=False) == ["BQMXNRZ"]
    assert r.get_kmers(8, strict=False) == []


@pytest.fixture(scope="function")
def one_seq():
    from cogent3 import make_seq

    return make_seq("AACCTGGAACC", moltype="dna")


@pytest.fixture(scope="session")
def worm_seq_path(DATA_DIR):
    return DATA_DIR / "c_elegans_WS199_dna_shortened.fasta"


@pytest.fixture(scope="session")
def worm_gff_path(DATA_DIR):
    return DATA_DIR / "c_elegans_WS199_shortened_gff.gff3"


def test_annotate_from_gff(worm_seq_path, worm_gff_path):
    """correctly annotates a Sequence from a gff file"""
    # duplicated in test_alignment on seq collection
    seq = cogent3.load_seq(worm_seq_path, moltype="dna")

    seq.annotate_from_gff(worm_gff_path)
    matches = list(seq.get_features())
    assert len(matches) == 11
    matches = list(seq.get_features(biotype="gene"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="mRNA"))
    assert len(matches) == 1
    matches = list(matches[0].get_children(biotype="exon"))
    assert len(matches) == 3


@pytest.mark.parametrize("rc", (False, True))
def test_seq_repr(one_seq, rc):
    pat = re.compile("[ACGT]+")
    dna = one_seq.moltype
    if rc:
        expect = dna.rc(str(one_seq))
        seq = one_seq.rc()
    else:
        expect = str(one_seq)
        seq = one_seq

    got = pat.findall(repr(seq))[0]
    assert expect.startswith(got), (expect, got)


def test_annotation_from_slice_with_stride():
    seq = DNA.make_seq("AAACGCGCGAAAAAAA", name="s1")
    seq.add_feature(biotype="exon", name="ex1", spans=[(3, 9)])
    f = list(seq.get_features(name="ex1"))[0]
    assert str(f.get_slice()) == "CGCGCG"
    s1 = seq[1::2]
    f = list(s1.get_features(name="ex1"))[0]
    assert str(f.get_slice()) == "CCC"


def test_absolute_position_base_cases(one_seq):
    """with no offset or view, the absolute index should remain unchanged"""
    got = one_seq._seq.absolute_position(5)
    assert got == 5

    # an index outside the range of the sequence should raise an IndexError
    with pytest.raises(IndexError):
        one_seq._seq.absolute_position(20)

    with pytest.raises(IndexError):
        one_seq._seq.absolute_position(-20)


def test_absolute_position_positive(one_seq):
    # with an offset, the abs index should be offset + index
    one_seq.annotation_offset = 2
    got = one_seq._seq.absolute_position(2)
    assert got == 2 + 2

    # with an offset and start, the abs index should be offset + start + index
    view = one_seq[2::]
    view.annotation_offset = 2  # todo: do we want the annotation_offset to be preserved when slicing? I think yes
    got = view._seq.absolute_position(2)
    assert got == 2 + 2 + 2

    # with an offset, start and step, the abs index should be offset + start + index * step
    view = one_seq[2::2]
    view.annotation_offset = 2
    got = view._seq.absolute_position(2)
    assert got == 2 + 2 + 2 * 2


def test_relative_position_base_cases(one_seq):
    """with no offset or view, the absolute index should remain unchanged"""
    got = one_seq._seq.relative_position(5)
    assert got == 5

    # a -ve index  should raise an IndexError
    with pytest.raises(IndexError):
        one_seq._seq.relative_position(-5)


@pytest.fixture(scope="function")
def integer_seq():
    return SeqView("0123456789")


def test_relative_position(integer_seq):
    """This test checks if the method returns the correct relative positions when
    the given index precedes or exceeds the range of the SeqView."""

    view = integer_seq[1:9:]
    # view = "12345678"
    got = view.relative_position(0)
    # precedes the view, so should return -1
    assert got == -1
    # exceeds the view, but still returns a value
    got = view.relative_position(10)
    assert got == 9


def test_relative_position_step_GT_one(integer_seq):
    """This test checks if the method returns the correct relative positions when
    the given index precedes or exceeds the range of the SeqView with a step greater than one.
    """

    # precedes the view, with step > 1
    view = integer_seq[2:7:2]
    # view = "246", precedes the view by 1 step
    got = view.relative_position(0)
    assert got == -1
    # precedes the view by 0.5 step, default behaviour is to round up to 0
    got = view.relative_position(1)
    assert got == 0
    # exceeds the view by two steps, len(view) + 2 = 4
    got = view.relative_position(10)
    assert got == 4


@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("rev", (False, True))
def test_seqview_copy(sliced, rev, integer_seq):
    raw_data = integer_seq.seq
    if rev:
        integer_seq = integer_seq[::-1]
        raw_data = raw_data[::-1]

    slice_start = 2
    slice_end = 4
    sv = integer_seq[slice_start:slice_end]
    copied = sv.copy(sliced=sliced)

    assert copied.value == raw_data[slice_start:slice_end]
    assert copied.reverse == integer_seq.reverse
    assert sliced and copied.seq is not sv.seq or copied.seq is integer_seq.seq


def test_relative_position_with_remainder(integer_seq):
    """tests relative_position when the index given is excluded from the view as it falls on
    a position that is 'stepped over'"""
    view = integer_seq[1:9:2]
    # view = "1357"
    got = view.relative_position(2)
    # 2 is stepped over in the view, so we return the index of 3 (which is 1)
    assert got == 1

    # setting the arg stop=True will adjust to the largest number, smaller than the given abs value, that is in the view
    got = view.relative_position(8, stop=True)
    # 8 is excluded from the view, so we return the index of 7 (which is 3)
    assert got == 3


@pytest.mark.parametrize("value", (0, 3))
@pytest.mark.parametrize("offset", (None, 1, 2))
@pytest.mark.parametrize("start", (None, 1, 2))
@pytest.mark.parametrize("stop", (None, 10, 11))
@pytest.mark.parametrize("step", (None, 1, 2))
def test_absolute_relative_roundtrip(one_seq, value, offset, start, stop, step):
    # a round trip from relative to absolute then from absolute to relative, should return the same value we began with
    view = one_seq[start:stop:step]
    view.annotation_offset = offset or 0
    abs_val = view._seq.absolute_position(value)
    rel_val = view._seq.relative_position(abs_val)
    assert rel_val == value


@pytest.mark.parametrize("value", (0, 2))
@pytest.mark.parametrize("offset", (None, 1, 2))
@pytest.mark.parametrize("start", (None, -1, -2))
@pytest.mark.parametrize("stop", (None, -10))
@pytest.mark.parametrize("step", (-1, -2))
def test_absolute_relative_roundtrip_reverse(
    integer_seq, value, offset, start, stop, step
):
    # a round trip from relative to absolute then from absolute to relative, should return the same value we began with
    view = integer_seq[start:stop:step]
    view.offset = offset or 0
    abs_val = view.absolute_position(value)
    rel_val = view.relative_position(abs_val)
    assert view.offset == (offset or 0)
    assert (view[rel_val]).value == view[value].value


def test_annotate_gff_nested_features(DATA_DIR):
    """correctly annotate a sequence with nested features"""
    # the synthetic example
    #          1111111111222222222333333333334
    # 1234567890123456789012345678901234567890
    #  **** biological_region
    #                                     ** biological_region
    #                                       * biological_region
    #      *******************************  gene
    #         *********************   mRNA
    #            *********            exon
    #                       *****     exon
    # ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC...
    seq = DNA.make_seq("ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC", name="22")
    gff3_path = DATA_DIR / "ensembl_sample.gff3"
    seq.annotate_from_gff(gff3_path)
    # we have 8 records in the gff file
    assert seq.annotation_db.num_matches() == 8

    # get the gene and check it has a single annotation and that
    # its slice is correct
    ann = list(seq.get_features(biotype="gene"))
    assert len(ann) == 1
    ann_seq = ann[0].get_slice()
    assert str(ann_seq) == "GGAAAATTTTTTTTTAAGGGGGAAAAAAAAA"
    # the gene has 1 transcript
    gene = ann[0]
    mrna = list(gene.get_children(biotype="mRNA"))
    assert len(mrna) == 1
    mrna = mrna[0]
    ann_seq = mrna.get_slice()
    assert str(ann_seq) == "AAAATTTTTTTTTAAGGGGGAAA"

    # the transcript has 2 exons, from the parent feature
    exons = list(mrna.get_children(biotype="exon"))
    assert len(exons) == 2
    # or the sequence
    ann = list(seq.get_features(biotype="exon"))
    assert len(ann) == 2
    exon_seqs = ("TTTTTTTTT", "GGGGG")
    assert tuple(str(ex.get_slice()) for ex in exons) == exon_seqs


def test_to_moltype_dna():
    """to_moltype("dna") ensures conversion from T to U"""
    seq = DNA.make_seq("AAAAGGGGTTT", name="seq1")
    rna = seq.to_moltype("rna")

    assert "T" not in rna


def test_to_moltype_rna():
    """to_moltype("rna") ensures conversion from U to T"""
    seq = RNA.make_seq("AAAAGGGGUUU", name="seq1")
    rna = seq.to_moltype("dna")

    assert "U" not in rna


@pytest.mark.parametrize("cls,with_offset", ((ArraySequence, False), (Sequence, True)))
def test_to_rich_dict(cls, with_offset):
    """Sequence to_dict works"""
    r = cls("AAGGCC", name="seq1")
    got = r.to_rich_dict()
    seq = "AAGGCC"

    if cls == Sequence:
        seq = SeqView(seq, seqid="seq1").to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": r.moltype.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
    }
    if with_offset:
        expect["annotation_offset"] = 0

    assert got == expect


@pytest.mark.parametrize("cls,with_offset", ((ArraySequence, False), (Sequence, True)))
def test_to_json(cls, with_offset):
    """to_json roundtrip recreates to_dict"""
    r = cls("AAGGCC", name="seq1")
    got = json.loads(r.to_json())

    seq = "AAGGCC"
    if cls == Sequence:
        seq = SeqView(seq, seqid="seq1").to_rich_dict()

    expect = {
        "name": "seq1",
        "seq": seq,
        "moltype": r.moltype.label,
        "info": None,
        "type": get_object_provenance(r),
        "version": __version__,
    }
    if with_offset:
        expect["annotation_offset"] = 0

    assert got == expect


def test_offset_with_multiple_slices(DATA_DIR):
    from cogent3.util.deserialise import deserialise_object

    seq = DNA.make_seq("ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC", name="22")
    gff3_path = DATA_DIR / "ensembl_sample.gff3"
    seq.annotate_from_gff(gff3_path)
    rd = seq[2:].to_rich_dict()
    s1 = deserialise_object(rd)
    assert s1.annotation_offset == 2
    rd = s1[3:].to_rich_dict()
    s2 = deserialise_object(rd)
    assert s2.annotation_offset == 5
    expect = {(f.seqid, f.biotype, f.name) for f in seq.get_features(start=5)}
    got = {(f.seqid, f.biotype, f.name) for f in s2.get_features()}
    assert got == expect


def test_seqview_to_rich_dict():
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sv = SeqView(seq=parent)
    plus = sv.to_rich_dict()
    minus = sv[::-1].to_rich_dict()
    plus = plus.pop("init_args")
    minus = minus.pop("init_args")
    assert plus.pop("seq") == minus.pop("seq")
    assert plus["step"] == -minus["step"]
    for coord in ("start", "stop"):
        assert coord not in plus
        assert coord not in minus


@pytest.mark.parametrize("reverse", (False, True))
def test_seqview_round_trip(reverse):
    from cogent3.util.deserialise import deserialise_object

    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sv = SeqView(seq=parent)
    if reverse:
        sv = sv[::-1]

    rd = sv.to_rich_dict()
    got = deserialise_object(rd)
    assert isinstance(got, SeqView)
    assert got.to_rich_dict() == sv.to_rich_dict()


@pytest.mark.parametrize("reverse", (False, True))
def test_sliced_seqview_rich_dict(reverse):
    parent = "ACCCCGGAAAATTTTTTTTTAAGGGGGAAAAAAAAACCCCCCC"
    sl = slice(2, 13)
    sv = SeqView(seq=parent)[sl]
    if reverse:
        sv = sv[::-1]
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seq"] == parent[sl]
    assert rd["init_args"]["offset"] == 2


@pytest.mark.parametrize(
    "sl",
    (
        slice(2, 5, 1),  # positive indices, positive step
        slice(-8, -5, 1),  # negative indices, positive step
        slice(4, 1, -1),  # positive indices, negative step
        slice(-6, -9, -1),  # negative indices, negative step
    ),
)
@pytest.mark.parametrize("offset", (4, 0))
def test_parent_start_stop(sl, offset):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = "234" if sl.step > 0 else "432"
    sv = SeqView(data)
    sv.offset = offset
    sv = sv[sl]
    assert sv.value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (2 + offset, 5 + offset)


@pytest.mark.parametrize(
    "sl",
    (
        slice(None, None, 1),  # slice whole sequence plus strand
        slice(None, None, -1),  # slice whole sequence minus strand
    ),
)
def test_parent_start_stop_limits(sl):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = SeqView(data)
    sv = sv[sl]
    assert sv.value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (0, 10)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_start_stop_empty(rev):
    data = "0123456789"
    # check our slice matches the expectation for rest of test
    expect = ""
    sv = SeqView(data)
    sv = sv[0 : 0 : -1 if rev else 1]
    assert sv.value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (0, 0)


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("index", range(9))
def test_parent_start_stop_singletons(index, rev):
    data = "0123456789"
    start, stop = (-(10 - index), -(10 - index + 1)) if rev else (index, index + 1)
    sl = slice(start, stop, -1 if rev else 1)
    # check our slice matches the expectation for rest of test
    expect = data[sl]
    sv = SeqView(data)
    sv = sv[sl]
    assert sv.value == expect
    # now check that start / stop are always the same
    # irrespective of step sign
    assert (sv.parent_start, sv.parent_stop) == (index, index + 1)


def test_get_drawable(DATA_DIR):
    seq = cogent3.load_seq(DATA_DIR / "annotated_seq.gb")
    seq = seq[2000:4000]
    biotypes = "CDS", "gene", "mRNA"
    for feat in seq.get_features(biotype=biotypes, allow_partial=True):
        draw = feat.get_drawable()
        assert "(incomplete)" in draw.text

    full = seq.get_drawable(biotype=biotypes)
    # should only include elements that overlap the segment
    assert len(full.traces) == len(biotypes)
    # and their names should indicate they're incomplete
    for trace in full.traces:
        assert "(incomplete)" in trace.text


@pytest.mark.parametrize("gc,seq", ((1, "TCCTGA"), (1, "ACGTAA---"), (2, "TCCAGG")))
def test_has_terminal_stop_true(gc, seq):
    gc = cogent3.get_code(gc)
    seq = cogent3.make_seq(seq, moltype="dna")
    assert seq.has_terminal_stop(gc=gc)


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_has_terminal_stop_false(gc, seq):
    gc = cogent3.get_code(gc)
    seq = cogent3.make_seq(seq, moltype="dna")
    assert not seq.has_terminal_stop(gc=gc)


def test_has_terminal_stop_strict():
    gc = cogent3.get_code(1)
    seq = cogent3.make_seq("TCCAG", moltype="dna")
    with pytest.raises(AlphabetError):
        seq.has_terminal_stop(gc=gc, strict=True)


@pytest.mark.parametrize(
    "gc,seq",
    (
        (2, "TCCAGG"),
        (1, "TAATGA"),
        (1, "ACGTGA---"),
        (1, "--AT-CTGA"),
    ),
)
def test_trim_terminal_stop_true(gc, seq):
    gc = cogent3.get_code(gc)
    expect = re.sub("(TGA|AGG)(?=[-]*$)", "---" if "-" in seq else "", seq)

    seq = cogent3.make_seq(seq, moltype="dna")
    got = str(seq.trim_stop_codon(gc=gc))
    assert got == expect


@pytest.mark.parametrize("gc,seq", ((1, "T?CTGC"), (2, "TCCAAG")))
def test_trim_terminal_stop_nostop(gc, seq):
    gc = cogent3.get_code(gc)
    seq = cogent3.make_seq(seq, moltype="dna")
    got = seq.trim_stop_codon(gc=gc)
    assert str(got) == str(seq)
    # since there's no stop, we just return the same object
    assert got is seq


@pytest.mark.parametrize(
    "gc,seq", ((1, "TCCAGG"), (2, "TCCAAA"), (1, "CCTGA"), (2, "CCAGG"))
)
def test_trim_terminal_stop_false(gc, seq):
    gc = cogent3.get_code(gc)
    seq = cogent3.make_seq(seq, moltype="dna")
    assert str(seq.trim_stop_codon(gc=gc)) == str(seq)


def test_trim_terminal_stop_strict():
    gc = cogent3.get_code(1)
    seq = cogent3.make_seq("TCCAG", moltype="dna")
    with pytest.raises(AlphabetError):
        seq.trim_stop_codon(gc=gc, strict=True)


@pytest.mark.parametrize("cast", (int, numpy.int32, numpy.int64, numpy.uint8))
def test_index_a_seq(cast):
    seq = cogent3.make_seq("TCCAG", moltype="dna")
    got = seq[cast(1)]
    assert isinstance(got, Sequence)


@pytest.mark.parametrize("cast", (float, numpy.float32))
def test_index_a_seq_float_fail(cast):
    seq = cogent3.make_seq("TCCAG", moltype="dna")
    index = cast(1)
    with pytest.raises(TypeError):
        seq[index]


@pytest.mark.parametrize("moltype", ("dna", "protein"))
def test_same_moltype(moltype):
    moltype = get_moltype(moltype)
    seq = moltype.make_seq("TCCAG")
    got = seq.to_moltype(moltype)
    assert got is seq


def test_gapped_by_map_segment_iter():
    moltype = get_moltype("dna")
    m, seq = moltype.make_seq("-TCC--AG").parse_out_gaps()
    g = list(seq.gapped_by_map_segment_iter(m, allow_gaps=True, recode_gaps=False))


@pytest.mark.parametrize("rev", (False, True))
@pytest.mark.parametrize("sliced", (False, True))
@pytest.mark.parametrize("start_stop", ((None, None), (3, 7)))
def test_copied_parent_coordinates(sliced, rev, start_stop):
    orig_name = "orig"
    seq = DNA.make_seq("ACGGTGGGAC", name=orig_name)
    start, stop = start_stop
    start = start or 0
    stop = stop or len(seq)
    sl = slice(start, stop)
    seq = seq[sl]
    sliced_name = "sliced"
    seq.name = sliced_name
    assert seq.name == sliced_name
    if rev:
        seq = seq.rc()
    copied = seq.copy(sliced=sliced)
    assert copied.name == sliced_name
    # matches original
    assert copied.parent_coordinates() == seq.parent_coordinates()
    # and expected -- the coordinate name always reflects the underlying sequence
    assert copied.parent_coordinates() == (orig_name, start, stop, -1 if rev else 1)


@pytest.mark.parametrize("rev", (False, True))
def test_parent_coordinates(rev):
    seq = DNA.make_seq("ACGGTGGGAC")
    seq = seq[1:1]
    if rev:
        seq = seq.rc()
    seq.name = "sliced"  # this assignment does not affect the
    # note that when a sequence has zero length, the parent seqid is None
    assert seq.parent_coordinates() == (None, 0, 0, 1)


def test_seqview_seqid():
    sv = SeqView("ACGGTGGGAC")
    assert sv.seqid is None

    sv = SeqView("ACGGTGGGAC", seqid="seq1")
    assert sv.seqid == "seq1"


def test_seqview_rich_dict_round_trip_seqid():
    sv = SeqView("ACGGTGGGAC", seqid="seq1")
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] == "seq1"

    got = SeqView.from_rich_dict(rd)
    assert got.seqid == "seq1"

    sv = SeqView("ACGGTGGGAC")
    rd = sv.to_rich_dict()
    assert rd["init_args"]["seqid"] == None

    got = SeqView.from_rich_dict(rd)
    assert got.seqid == None


def test_seqview_slice_propagates_seqid():
    sv = SeqView("ACGGTGGGAC", seqid="seq1")
    sliced_sv = sv[1:8:2]
    assert sliced_sv.seqid == "seq1"

    copied_sv = sliced_sv.copy(sliced=False)
    assert copied_sv.seqid == "seq1"

    copied_sliced_sv = sliced_sv.copy(sliced=True)
    assert copied_sliced_sv.seqid == "seq1"


@pytest.mark.parametrize("cls", (Aligned, Sequence, SeqView, ArraySequence, str, bytes))
def test_coerce_to_seqview(cls):
    seq = "AC--GGTGGGAC"
    seqid = "seq1"
    if cls in (str, bytes):
        got = _coerce_to_seqview(seq, seqid, preserve_case=True, checker=(lambda x: x))
    elif cls is Aligned:
        got = _coerce_to_seqview(
            cls(*Sequence(seq).parse_out_gaps()),
            seqid,
            preserve_case=True,
            checker=(lambda x: x),
        )
    else:
        got = _coerce_to_seqview(
            cls(seq), seqid, preserve_case=True, checker=(lambda x: x)
        )
    assert got.value == seq
    assert isinstance(got, SeqView)


def test_sequences_propogates_seqid():
    # creating a name Sequence propagates the seqid to the SeqView.
    seq = Sequence("ACGGTGGGAC", name="seq1")
    assert seq._seq.seqid == "seq1"

    # renaming the Sequence deosnt change the seqid of the SeqView.
    seq.name = "seq2"
    assert seq.name == "seq2"
    assert seq._seq.seqid == "seq1"

    # creating a name Aligned propagates the seqid to the SeqView.
    seq = Aligned(*Sequence("ACG--G--GAC", name="seq1").parse_out_gaps())
    assert seq.data._seq.seqid == "seq1"

    # creating a Sequence with a seqview does not change the seqid of the SeqView.
    seq = Sequence(SeqView("ACGGTGGGAC", seqid="parent_name"), name="seq_name")
    assert seq.name == "seq_name"
    assert seq._seq.seqid == "parent_name"

    # creating a Sequence with an unnamed seqview does not name the SeqView.
    seq = Sequence(SeqView("ACGGTGGGAC"), name="seq_name")
    assert seq.name == "seq_name"
    assert seq._seq.seqid == None


def test_make_seq_assigns_to_seqview():
    seq = cogent3.make_seq("ACGT", name="s1")
    assert seq.name == seq._seq.seqid == "s1"


def test_empty_seqview_translate_position():
    sv = SeqView("")
    assert sv.absolute_position(0) == 0
    assert sv.relative_position(0) == 0


@pytest.mark.parametrize("start", (None, 0, 1, 10, -1, -10))
@pytest.mark.parametrize("stop", (None, 10, 8, 1, 0, -1, -11))
@pytest.mark.parametrize("step", (None, 1, 2, -1, -2))
@pytest.mark.parametrize("length", (1, 8, 999))
def test_seqview_seq_len_init(start, stop, step, length):
    # seq_len is length of seq when None
    seq_data = "A" * length
    sv = SeqView(seq_data, start=start, stop=stop, step=step)
    expect = len(seq_data)
    # Check property and slot
    assert sv.seq_len == expect
    assert sv._seq_len == expect


@pytest.mark.parametrize("seq, seq_len", [("A", 0), ("", 1), ("A", 2)])
def test_seqview_seq_len_mismatch(seq, seq_len):
    # If provided, seq_len must match len(seq)
    with pytest.raises(AssertionError):
        SeqView(seq, seq_len=seq_len)


def test_seqview_copy_propagates_seq_len():
    seq = "ACGGTGGGAC"
    sv = SeqView(seq)
    copied = sv.copy()
    assert copied.seq_len == len(seq)


def test_seqview_seq_len_modified_seq():
    seq = "ACGGTGGGAC"
    sv = SeqView(seq)

    sv.seq = "ATGC"  # this should not modify seq_len
    assert sv.seq_len == len(seq)
