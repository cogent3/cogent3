import unittest

from cogent3 import DNA, make_aligned_seqs
from cogent3.core.location import FeatureMap, Span


class MapTest(unittest.TestCase):
    """Testing annotation of and by maps"""

    def test_spans(self):
        # a simple two part map of length 10
        map = FeatureMap([(0, 5), (5, 10)], parent_length=10)
        # try different spans on the above map
        for (start, end), expected in [
            ((0, 4), "[0:4]"),
            ((0, 5), "[0:5]"),
            ((0, 6), "[0:5, 5:6]"),
            ((5, 10), "[5:10]"),
            ((-1, 10), "[-1-, 0:5, 5:10]"),
            ((5, 11), "[5:10, -1-]"),
            ((0, 10), "[0:5, 5:10]"),
            ((10, 0), "[10:5, 5:0]"),
        ]:
            r = repr(Span(start, end, reverse=start > end).remap_with(map))
            if r != expected:
                self.fail(repr((r, expected)))

    def test_get_by_annotation(self):
        seq = DNA.make_seq("ATCGATCGAT" * 5, name="base")
        seq.add_feature(biotype="test_type", name="test_label", spans=[(5, 10)])
        seq.add_feature(biotype="test_type", name="test_label2", spans=[(15, 18)])

        answer = list(seq.get_features(biotype="test_type"))
        self.assertEqual(len(answer), 2)
        self.assertEqual(str(seq[answer[0]]), "TCGAT")
        self.assertEqual(str(seq[answer[1]]), "TCG")

        answer = list(seq.get_features(biotype="test_type", name="test_label"))
        self.assertEqual(len(answer), 1)
        self.assertEqual(str(seq[answer[0]]), "TCGAT")

        # test ignoring of a partial annotation
        sliced_seq = seq[:17]
        answer = list(sliced_seq.get_features(biotype="test_type", allow_partial=False))
        self.assertEqual(len(answer), 1)
        self.assertEqual(str(sliced_seq[answer[0]]), "TCGAT")


def test_get_by_seq_annotation_allow_gaps():
    aln = make_aligned_seqs(
        data={"a": "ATCGAAATCGAT", "b": "ATCGA--TCGAT"}, array_align=False
    )
    # original version was putting annotation directly on seq
    f = aln.add_feature(
        seqid="b",
        biotype="test_type",
        name="test_label",
        spans=[(4, 6)],
        on_alignment=False,
    )
    # in the original get_by_seq_annotations(), inside the method it
    # used self[feature.map.start : feature.map.end],
    # so what was returned was an alignment slice compared to below,
    # which is only non-gap positions of "b"
    # when we set allow_gaps, we get gaps in an alignment
    f = list(aln.get_features(seqid="b", biotype="test_type"))[0]
    got = f.get_slice(allow_gaps=True).to_dict()
    assert got == {"b": "A--T", "a": "AAAT"}
    # otherwise we only get the exact positions represented
    got = f.get_slice(allow_gaps=False).to_dict()
    assert got == {"b": "AT", "a": "AT"}
