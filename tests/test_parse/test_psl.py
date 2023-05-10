#!/usr/bin/env python
"""Unit tests for the PSL parser.
   Compatible with blat v.34
"""

from unittest import TestCase, main

from cogent3.parse.psl import MinimalPslParser, PslToTable


fname = "data/test.psl"


class Test(TestCase):
    def test_header(self):
        """should return correct header"""
        expect = [
            "match",
            "mis-match",
            "rep. match",
            "N's",
            "Q gap count",
            "Q gap bases",
            "T gap count",
            "T gap bases",
            "strand",
            "Q name",
            "Q size",
            "Q start",
            "Q end",
            "T name",
            "T size",
            "T start",
            "T end",
            "block count",
            "blockSizes",
            "qStarts",
            "tStarts",
        ]
        infile = open(fname)
        parser = MinimalPslParser(infile)
        next(parser)
        header = next(parser)
        infile.close()
        self.assertEqual(header, expect)

    def test_psl_to_table(self):
        PslToTable(fname)

    def test_getting_seq_coords(self):
        """get correct sequence coordinates to produce a trimmed sequence"""
        table = PslToTable(fname)
        for row in table:
            query_name = row["Q name"]
            query_strand = row["strand"]
            q_start = row["Q start"]


if __name__ == "__main__":
    main()
