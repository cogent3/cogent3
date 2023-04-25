#!/usr/bin/env python

import unittest

from cogent3 import DNA, make_aligned_seqs
from cogent3.core.location import Map, Span


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight", "Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class MapTest(unittest.TestCase):
    """Testing annotation of and by maps"""

    def test_spans(self):
        # a simple two part map of length 10
        map = Map([(0, 5), (5, 10)], parent_length=10)
        # try different spans on the above map
        for ((start, end), expected) in [
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
        seq.add_annotation(Feature, "test_type", "test_label", [(5, 10)])
        seq.add_annotation(Feature, "test_type", "test_label2", [(15, 18)])

        answer = list(seq.get_by_annotation("test_type"))
        self.assertEqual(len(answer), 2)
        self.assertEqual(str(answer[0]), "TCGAT")
        self.assertEqual(str(answer[1]), "TCG")

        answer = list(seq.get_by_annotation("test_type", "test_label"))
        self.assertEqual(len(answer), 1)
        self.assertEqual(str(answer[0]), "TCGAT")

        # test ignoring of a partial annotation
        sliced_seq = seq[:17]
        answer = list(sliced_seq.get_by_annotation("test_type", ignore_partial=True))
        self.assertEqual(len(answer), 1)
        self.assertEqual(str(answer[0]), "TCGAT")

    def test_get_by_seq_annotation(self):
        aln = make_aligned_seqs(
            data={"a": "ATCGAAATCGAT", "b": "ATCGA--TCGAT"}, array_align=False
        )
        b = aln.get_seq("b")
        b.add_annotation(Feature, "test_type", "test_label", [(4, 6)])

        answer = aln.get_by_seq_annotation("b", "test_type")[0].to_dict()
        self.assertEqual(answer, {"b": "A--T", "a": "AAAT"})
