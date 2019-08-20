#!/usr/bin/env python
"""Unit tests for the PSL parser.
   Compatible with blat v.34
"""

from unittest import TestCase, main

from cogent3 import LoadTable
from cogent3.parse.psl import MinimalPslParser, PslToTable, make_header


__author__ = "Gavin Huttley, Anuj Pahwa"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Gavin Huttley", "Anuj Pahwa"]
__license__ = "BSD-3"
__version__ = "2019.8.20a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Development"

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
        version = next(parser)
        header = next(parser)
        infile.close()
        self.assertEqual(header, expect)

    def test_psl_to_table(self):
        table = PslToTable(fname)

    def test_getting_seq_coords(self):
        """get correct sequence coordinates to produce a trimmed sequence"""
        table = PslToTable(fname)
        for row in table:
            query_name = row["Q name"]
            query_strand = row["strand"]
            q_start = row["Q start"]


if __name__ == "__main__":
    main()
