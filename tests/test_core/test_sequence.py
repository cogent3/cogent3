#!/usr/bin/env python
"""Unit tests for Sequence class and its subclasses.
"""

import json
import os
import re

from pickle import dumps
from unittest import TestCase, main

import pytest

from numpy import array
from numpy.testing import assert_allclose, assert_equal

from cogent3.core.annotation import Feature, SimpleVariable, Variable
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
)
from cogent3.util.misc import get_object_provenance


__author__ = "Rob Knight, Gavin Huttley and Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Peter Maxwell", "Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
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
        annot1 = s.add_feature("exon", "annot1", [(0, 10)])
        annot2 = s.add_feature("exon", "annot2", [(10, 14)])
        got = s.copy()
        got_annot1 = got.get_features_matching(feature_type="exon", name="annot1")[0]
        got_annot2 = got.get_features_matching(feature_type="exon", name="annot2")[0]
        self.assertIsNot(got, s)
        self.assertIsNot(got_annot1, annot1)
        self.assertIsNot(got_annot2, annot2)
        self.assertEqual(got.name, s.name)
        self.assertEqual(got.info, s.info)
        self.assertEqual(str(got), str(s))
        self.assertEqual(got.moltype, s.moltype)
        annot1_slice = str(annot1.get_slice())
        annot2_slice = str(annot2.get_slice())
        got1_slice = str(got.annotations[0].get_slice())
        got2_slice = str(got.annotations[1].get_slice())
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

    def test_to_rich_dict(self):
        """Sequence to_dict works"""
        r = self.SEQ("AAGGCC", name="seq1")
        got = r.to_rich_dict()
        expect = {
            "name": "seq1",
            "seq": "AAGGCC",
            "moltype": r.moltype.label,
            "info": None,
            "type": get_object_provenance(r),
            "version": __version__,
        }
        self.assertEqual(got, expect)

    def test_to_json(self):
        """to_json roundtrip recreates to_dict"""
        r = self.SEQ("AAGGCC", name="seq1")
        got = json.loads(r.to_json())
        expect = {
            "name": "seq1",
            "seq": "AAGGCC",
            "moltype": r.moltype.label,
            "info": None,
            "type": get_object_provenance(r),
            "version": __version__,
        }
        self.assertEqual(got, expect)

    def test_sequence_to_moltype(self):
        """correctly convert to specified moltype"""
        s = Sequence("TTTTTTTTTTAAAA", name="test1")
        annot1 = s.add_annotation(Feature, "exon", "fred", [(0, 10)])
        annot2 = s.add_annotation(Feature, "exon", "trev", [(10, 14)])
        got = s.to_moltype("rna")
        annot1_slice = str(annot1.get_slice())
        annot2_slice = str(annot2.get_slice())
        got1_slice = str(got.annotations[0].get_slice())
        got2_slice = str(got.annotations[1].get_slice())
        self.assertNotEqual(annot1_slice, got1_slice)
        self.assertEqual(annot2_slice, got2_slice)
        self.assertEqual(got.moltype.label, "rna")
        self.assertEqual(got.name, "test1")

        s = Sequence("AAGGGGAAAACCCCCAAAAAAAAAATTTTTTTTTTAAA", name="test2")
        xx_y = [[[2, 6], 2.4], [[10, 15], 5.1], [[25, 35], 1.3]]
        y_valued = s.add_annotation(Variable, "SNP", "freq", xx_y)
        got = s.to_moltype("rna")
        y_valued_slice = str(y_valued.get_slice())
        got_slice = str(str(got.annotations[0].get_slice()))
        self.assertNotEqual(y_valued_slice, got_slice)
        self.assertEqual(got.moltype.label, "rna")
        self.assertEqual(got.name, "test2")

        s = Sequence("TTTTTTTTTTAAAAAAAAAA", name="test3")
        data = [i for i in range(20)]
        annot4 = s.add_annotation(SimpleVariable, "SNP", "freq", data)
        got = s.to_moltype(RNA)
        annot4_slice = str(annot4.get_slice())
        got_slice = str(str(got.annotations[0].get_slice()))
        self.assertNotEqual(annot4_slice[:10], got_slice[:10])
        self.assertEqual(annot4_slice[10:20], got_slice[10:20])
        self.assertEqual(got.moltype.label, "rna")
        self.assertEqual(got.name, "test3")

        # calling with a null object should raise an exception
        with self.assertRaises(ValueError):
            s.to_moltype(None)

        with self.assertRaises(ValueError):
            s.to_moltype("")

    def test_annotate_from_gff(self):
        """correctly annotates a Sequence from a gff file"""
        from cogent3.parse.fasta import FastaParser

        fasta_path = os.path.join("data/c_elegans_WS199_dna_shortened.fasta")
        gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
        name, seq = next(FastaParser(fasta_path))

        sequence = Sequence(seq)
        sequence.annotate_from_gff(gff3_path)
        matches = [m for m in sequence.get_features_matching()]
        # 13 features with one having 2 parents, so 14 instances should be found
        self.assertEqual(len(matches), 14)

    @pytest.mark.xfail(
        reason="todo: annotate_from_gff needs to be using bound annotation_db"
    )
    def test_annotate_gff_nested_features(self):
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
        gff3_path = os.path.join("data/ensembl_sample.gff3")
        seq.annotate_from_gff(gff3_path)
        # we have 1 "full chromosome" annotation, 3 generic regions and 1 gene
        self.assertEqual(len(seq.annotations), 5)

        # get the gene and check it has a single annotation and that
        # its slice is correct
        ann = seq.get_features_matching("gene")
        self.assertEqual(len(ann), 1)
        self.assertEqual(len(ann[0].annotations), 1)
        seq = ann[0].get_slice()
        self.assertEqual(str(seq), "GGAAAATTTTTTTTTAAGGGGGAAAAAAAAA")

        # the gene has 1 transcript
        ann = seq.get_features_matching("mRNA", extend_query=True)
        self.assertEqual(len(ann), 1)
        self.assertEqual(len(ann[0].annotations), 2)  # 2 exons
        seq = ann[0].get_slice()
        self.assertEqual(str(seq), "AAAATTTTTTTTTAAGGGGGAAA")

        # the transcript has 2 exons
        ann = seq.get_features_matching("exon", extend_query=True)
        self.assertEqual(len(ann), 2)
        exon_seqs = ("TTTTTTTTT", "GGGGG")
        for x in ann:
            self.assertEqual(len(x.annotations), 0)
            self.assertTrue(str(x.get_slice()) in exon_seqs, msg=x.get_slice())

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
        from cogent3.core.annotation import _Annotatable

        s = self.SEQ("ACGGCTGAAGCGCTCCGGGTTTAAAACG")
        annotatable = isinstance(s, _Annotatable)
        if annotatable:
            self.assertFalse(s.is_annotated())
            _ = s.add_feature("gene", "blah", [(0, 10)])
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
        alphabet = DNA.alphabets.base ** 3

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
    # case 1: Short sequence, defaults
    seq = "ACGT"
    view = SeqView(seq)
    expected = "SeqView(seq='ACGT', start=0, stop=4, step=1)"
    assert repr(view) == expected

    # case 2: Long sequence
    seq = "ACGT" * 10
    view = SeqView(seq)
    expected = "SeqView(seq='ACGTACGTAC...TACGT', start=0, stop=40, step=1)"
    assert repr(view) == expected

    # case 3: Non-zero start, stop, and step values
    seq = "ACGT" * 10
    view = SeqView(seq, start=5, stop=35, step=2)
    expected = "SeqView(seq='ACGTACGTAC...TACGT', start=5, stop=35, step=2)"
    assert repr(view) == expected


def test_distance_indices():
    """ArraySequence distance should work with function of indices"""
    s1 = RNA.make_array_seq("UGCUGCUC", name="s1")
    s2 = RNA.make_seq(s1)


# run if called from command-line
if __name__ == "__main__":
    main()
