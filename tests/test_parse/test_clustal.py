#!/usr/bin/env python
"""Unit tests for the clustal parsers.
"""
from unittest import TestCase, main

from cogent3.parse.clustal import (
    MinimalClustalParser,
    delete_trailing_number,
    is_clustal_seq_line,
    last_space,
)
from cogent3.parse.record import RecordError


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

# Note: the data are all strings and hence immutable, so it's OK to define
# them here instead of in setUp and then subclassing everything from that
# base class. If the data were mutable, we'd need to take more precautions
# to avoid crossover between tests.

minimal = "abc\tucag"
two = "abc\tuuu\ndef\tccc\n\n    ***\n\ndef ggg\nabc\taaa\n".split("\n")

real = """CLUSTAL W (1.82) multiple sequence alignment


abc             GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA 60
def             ------------------------------------------------------------
xyz             ------------------------------------------------------------


abc             GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC 119
def             -----------------------------------------CGCGAUGCAUGCAU-CGAU 18
xyz             -------------------------------------CAUGCAUCGUACGUACGCAUGAC 23
                                                         *    * * * *    **

abc             UGACUAGUCAGCUAGCAUCGAUCAGU 145
def             CGAUCAGUCAGUCGAU---------- 34
xyz             UGCUGCAUCA---------------- 33
                *     ***""".split(
    "\n"
)

bad = ["dshfjsdfhdfsj", "hfsdjksdfhjsdf"]

space_labels = ["abc uca", "def ggg ccc"]


class clustalTests(TestCase):
    """Tests of top-level functions."""

    def test_is_clustal_seq_line(self):
        """is_clustal_seq_line should reject blanks and 'CLUSTAL'"""
        ic = is_clustal_seq_line
        assert ic("abc")
        assert ic("abc  def")
        assert not ic("CLUSTAL")
        assert not ic("CLUSTAL W fsdhicjkjsdk")
        assert not ic("  *   *")
        assert not ic(" abc def")
        assert not ic("MUSCLE (3.41) multiple sequence alignment")

    def test_last_space(self):
        """last_space should split on last whitespace"""
        self.assertEqual(last_space("a\t\t\t  b    c"), ["a b", "c"])
        self.assertEqual(last_space("xyz"), ["xyz"])
        self.assertEqual(last_space("  a b"), ["a", "b"])

    def test_delete_trailing_number(self):
        """delete_trailing_number should delete the trailing number if present"""
        dtn = delete_trailing_number
        self.assertEqual(dtn("abc"), "abc")
        self.assertEqual(dtn("a b c"), "a b c")
        self.assertEqual(dtn("a \t  b  \t  c"), "a \t  b  \t  c")
        self.assertEqual(dtn("a b 3"), "a b")
        self.assertEqual(dtn("a b c \t 345"), "a b c")


class MinimalClustalParserTests(TestCase):
    """Tests of the MinimalClustalParser class"""

    def test_null(self):
        """MinimalClustalParser should return empty dict and list on null input"""
        result = MinimalClustalParser([])
        self.assertEqual(result, ({}, []))

    def test_minimal(self):
        """MinimalClustalParser should handle single-line input correctly"""
        result = MinimalClustalParser([minimal])  # expects seq of lines
        self.assertEqual(result, ({"abc": ["ucag"]}, ["abc"]))

    def test_two(self):
        """MinimalClustalParser should handle two-sequence input correctly"""
        result = MinimalClustalParser(two)
        self.assertEqual(
            result, ({"abc": ["uuu", "aaa"], "def": ["ccc", "ggg"]}, ["abc", "def"])
        )

    def test_real(self):
        """MinimalClustalParser should handle real Clustal output"""
        data, labels = MinimalClustalParser(real)
        self.assertEqual(labels, ["abc", "def", "xyz"])
        self.assertEqual(
            data,
            {
                "abc": [
                    "GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA",
                    "GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC",
                    "UGACUAGUCAGCUAGCAUCGAUCAGU",
                ],
                "def": [
                    "------------------------------------------------------------",
                    "-----------------------------------------CGCGAUGCAUGCAU-CGAU",
                    "CGAUCAGUCAGUCGAU----------",
                ],
                "xyz": [
                    "------------------------------------------------------------",
                    "-------------------------------------CAUGCAUCGUACGUACGCAUGAC",
                    "UGCUGCAUCA----------------",
                ],
            },
        )

    def test_bad(self):
        """MinimalClustalParser should reject bad data if strict"""
        result = MinimalClustalParser(bad, strict=False)
        self.assertEqual(result, ({}, []))
        # should fail unless we turned strict processing off
        self.assertRaises(RecordError, MinimalClustalParser, bad)

    def test_space_labels(self):
        """MinimalClustalParser should tolerate spaces in labels"""
        result = MinimalClustalParser(space_labels)
        self.assertEqual(
            result, ({"abc": ["uca"], "def ggg": ["ccc"]}, ["abc", "def ggg"])
        )


if __name__ == "__main__":
    main()
