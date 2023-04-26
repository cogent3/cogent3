import unittest

import pytest

from cogent3 import DNA, make_aligned_seqs, make_unaligned_seqs
from cogent3.core.annotation import Feature, _Feature
from cogent3.core.location import Map, Span, as_map
from cogent3.core.sequence import DnaSequence, RnaSequence


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def makeSampleSequence(name, with_gaps=False):
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    utr = (12, 15)
    if with_gaps:
        raw_seq = f"{raw_seq[:5]}-----{raw_seq[10:-2]}--"
    # name is required for creating annotations
    seq = DNA.make_seq(raw_seq, name=name)
    seq.add_feature(biotype="CDS", name="CDS", spans=[cds])
    seq.add_feature(biotype="5'UTR", name="5' UTR", spans=[utr])
    return seq


def makeSampleAlignment():
    seq1 = makeSampleSequence("FAKE01")
    seq2 = makeSampleSequence("FAKE02", with_gaps=True)
    seqs = {seq1.name: seq1, seq2.name: seq2}
    aln = make_aligned_seqs(data=seqs, array_align=False)
    aln.add_feature(
        biotype="misc_feature", name="misc", spans=[(12, 25)], on_alignment=True
    )
    aln.add_feature(biotype="CDS", name="blue", spans=[(15, 25)], on_alignment=True)
    aln.add_feature(biotype="5'UTR", name="red", spans=[(2, 4)], on_alignment=True)
    aln.add_feature(biotype="LTR", name="fake", spans=[(2, 15)], on_alignment=True)
    return aln


class TestAnnotations(unittest.TestCase):
    def setUp(self):
        self.seq = makeSampleSequence("seq1")
        self.aln = makeSampleAlignment()

    def test_slice_seq_with_annotations(self):
        newseq = self.seq[:5] + self.seq[10:]
        for annot_type in ["CDS", "5'UTR"]:
            orig = str(list(self.seq.get_by_annotation(annot_type))[0])
            new = str(list(newseq.get_by_annotation(annot_type))[0])
            assert orig == new, (annot_type, orig, new)

    def test_add_annotated_seqs_drops_annotations(self):
        # retain link to annotation db as long as a simple slice
        a = self.seq[:5]
        b = self.seq[10:]
        assert a.annotation_db is not None
        assert b.annotation_db is not None
        # but adding two seqs drops the annotation db since too
        # difficult to track coords
        newseq = a + b
        assert newseq.annotation_db is None

    def test_aln_annotations(self):
        """test that annotations to alignment and its' sequences"""
        aln_expecteds = {
            "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
            "LTR": {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"},
        }
        seq_expecteds = {
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "TTT", "FAKE02": "TTT"},
        }
        for annot_type in ["misc_feature", "CDS", "5'UTR", "LTR"]:
            observed = list(self.aln.get_by_annotation(annot_type))[0].to_dict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, (annot_type, expected, observed)
            if annot_type in ["misc_feature", "LTR"]:
                continue  # because seqs haven't been annotated with it
            for name in self.aln.names:
                observed = list(
                    self.aln.named_seqs[name].data.get_by_annotation(annot_type)
                )[0]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, (annot_type, name, expected, observed)

    def test_slice_aln_with_annotations(self):
        """test that annotations of sequences and alignments survive alignment
        slicing."""
        aln_expecteds = {
            "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
            "LTR": {"FAKE01": "CCCTTTTT", "FAKE02": "CCCTTTTT"},
        }
        newaln = self.aln[:5] + self.aln[10:]
        for annot_type in ["LTR", "misc_feature", "CDS", "5'UTR"]:
            feature_list = newaln.get_features_matching(annot_type)
            new = newaln.get_region_covering_all(feature_list).get_slice().to_dict()
            expected = aln_expecteds[annot_type]
            assert expected == new, (annot_type, expected, new)
            if annot_type in ["misc_feature", "LTR"]:
                continue  # because seqs haven't been annotated with it
            for name in self.aln.names:
                orig = str(
                    list(self.aln.get_annotations_from_seq(name, annot_type))[
                        0
                    ].get_slice()
                )
                new = str(
                    list(newaln.get_annotations_from_seq(name, annot_type))[
                        0
                    ].get_slice()
                )
                assert orig == new, (name, annot_type, orig, new)

    def test_feature_projection(self):
        expecteds = {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"}
        aln_ltr = self.aln.get_features_matching("LTR")[0]
        for seq_name in ["FAKE01", "FAKE02"]:
            expected = expecteds[seq_name]
            seq_ltr = self.aln.project_annotation(seq_name, aln_ltr)
            if "-" in expected:
                self.assertRaises(ValueError, seq_ltr.get_slice)
                seq_ltr = seq_ltr.without_lost_spans()
                expected = expected.replace("-", "")
            self.assertEqual(seq_ltr.get_slice(), expected)

    def test_feature_copy_annotations_to(self):
        """test correct copy of annotations"""
        orig = DnaSequence("TTTTTTTTTTAAAA", name="Orig")
        annot = orig.add_annotation(Feature, "exon", "fred", [(0, 14)])
        seq = RnaSequence("UUUUUUUUUUAAAA", name="Test")
        got = annot.copy_annotations_to(seq)
        self.assertEqual(len(orig.annotations), len(got.annotations))
        for src, dest in zip(orig.annotations, got.annotations):
            self.assertEqual(src.get_coordinates(), dest.get_coordinates())
            self.assertIsInstance(src, dest.__class__)
            self.assertIs(dest.parent, seq)
        with self.assertRaises(AssertionError):
            _ = annot.copy_annotations_to(seq[:-2])

    def test_reverse_complement(self):
        """test correct translation of annotations on reverse complement."""
        aln_expecteds = {
            "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
            "LTR": {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"},
        }

        seq_expecteds = {
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "TTT", "FAKE02": "TTT"},
        }

        rc = self.aln.rc()
        # rc'ing an Alignment or Sequence rc's their annotations too. This means
        # slicing returns the same sequence as the non-rc'd alignment/seq
        for annot_type in ["misc_feature", "CDS", "5'UTR", "LTR"]:
            observed = list(self.aln.get_by_annotation(annot_type))[0].to_dict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, ("+", annot_type, expected, observed)
            observed = list(rc.get_by_annotation(annot_type))[0].to_dict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, ("-", annot_type, expected, observed)

            if annot_type in ["misc_feature", "LTR"]:
                continue  # because seqs haven't been annotated with it
            for name in self.aln.names:
                observed = list(
                    self.aln.named_seqs[name].data.get_by_annotation(annot_type)
                )[0]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, (
                    "+",
                    annot_type,
                    name,
                    expected,
                    observed,
                )
                observed = list(rc.named_seqs[name].data.get_by_annotation(annot_type))[
                    0
                ]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, (
                    "-",
                    annot_type,
                    name,
                    expected,
                    observed,
                )


class TestMapSpans(unittest.TestCase):
    """Test attributes of Map & Spans classes critical to annotation
    manipulation."""

    def test_span(self):
        forward = Span(20, 30)
        reverse = Span(70, 80, reverse=True)
        assert forward.reversed_relative_to(100) == reverse
        assert reverse.reversed_relative_to(100) == forward

    def test_map(self):
        """reversing a map with multiple spans should preserve span relative
        order"""
        forward = [Span(20, 30), Span(40, 50)]
        fmap = Map(spans=forward, parent_length=100)
        fmap_reversed = fmap.nucleic_reversed()
        reverse = [Span(70, 80, reverse=True), Span(50, 60, reverse=True)]
        rmap = Map(spans=reverse, parent_length=100)
        for i in range(2):
            self.assertEqual(fmap_reversed.spans[i], rmap.spans[i])


@pytest.mark.parametrize("alignment", (False, True))
def test_constructing_collections(alignment):
    seq1 = makeSampleSequence("FAKE01")
    seq2 = makeSampleSequence("FAKE02", with_gaps=True)
    seqs = {"FAKE01": seq1, "FAKE02": seq2}
    expect = sum(s.annotation_db.num_matches() for s in seqs.values())
    if alignment:
        coll = make_aligned_seqs(data=seqs, moltype="dna", array_align=False)
    else:
        coll = make_unaligned_seqs(data=seqs, moltype="dna")

    assert coll.annotation_db.num_matches() == expect
