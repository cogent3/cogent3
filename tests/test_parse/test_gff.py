#!/usr/bin/env python
"""Unit tests for GFF and related parsers.
"""
import os

from io import StringIO
from pathlib import Path
from unittest import TestCase, main

from cogent3.parse.gff import *


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"

headers = [
    """##gff-version 2 
##source-version <source> <version text> 
##date <date> 
##Type <type> [<seqname>] 
##DNA <seqname>
##acggctcggattggcgctggatgatagatcagacgac
##...
##end-DNA
""",
    """##gff-version 2
""",
    "",
]

#    '<seqname>\t<source>\t<feature>\t<start>\t<end>\t<score>\t<strand>\t<frame>\t[attribute]\n'

data_lines = [
    (
        'seq1\tBLASTX\tsimilarity\t101\t235\t87.1\t+\t0\tTarget "HBA_HUMAN" 11 55 ; E_value 0.0003\n',
        (
            "seq1",
            "BLASTX",
            "similarity",
            100,
            235,
            "87.1",
            "+",
            "0",
            'Target "HBA_HUMAN" 11 55 ; E_value 0.0003',
            None,
        ),
    ),
    (
        'dJ102G20\tGD_mRNA\tcoding_exon\t7105\t7201\t.\t-\t2\tSequence "dJ102G20.C1.1"\n',
        (
            "dJ102G20",
            "GD_mRNA",
            "coding_exon",
            7201,
            7104,
            ".",
            "-",
            "2",
            'Sequence "dJ102G20.C1.1"',
            None,
        ),
    ),
    (
        "dJ102G20\tGD_mRNA\tcoding_exon\t7105\t7201\t.\t-\t2\t\n",
        ("dJ102G20", "GD_mRNA", "coding_exon", 7201, 7104, ".", "-", "2", "", None),
    ),
    (
        '12345\tSource with spaces\tfeature with spaces\t-100\t3600000000\t1e-5\t-\t.\tSequence "BROADO5" ; Note "This is a \\t tab containing \\n multi line comment"\n',
        (
            "12345",
            "Source with spaces",
            "feature with spaces",
            3600000000,
            101,
            "1e-5",
            "-",
            ".",
            'Sequence "BROADO5" ; Note "This is a \\t tab containing \\n multi line comment"',
            None,
        ),
    ),
]


class GffTest(TestCase):
    """Setup data for all the GFF parsers."""

    def testGffParserData(self):
        """Test gff_parser with valid data lines"""
        for (line, canned_result) in data_lines:
            result = next(gff_parser(StringIO(line)))
            canned_result = list(canned_result)
            self.assertEqual(result.pop("Attributes")["Info"], canned_result.pop(8))
            self.assertEqual(set(result.values()), set(canned_result))

    def test_gff_parser_headers(self):
        """Test gff_parser with valid data headers"""
        data = "".join([x[0] for x in data_lines])
        for header in headers:
            result = list(gff_parser(StringIO(header + data)))
            lines = [(x[0], list(x[1])) for x in data_lines]
            self.assertEqual(
                [l.pop("Attributes")["Info"] for l in result],
                [x[1].pop(8) for x in lines],
            )
            self.assertEqual(
                [set(l.values()) for l in result], [set(x[1]) for x in lines]
            )

    def test_parse_attributes_gff2(self):
        """Test the parse_attributes_gff2 method"""
        self.assertEqual(
            [
                parse_attributes_gff2(x[1][8], (x[1][3], x[1][4]))["ID"]
                for x in data_lines
            ],
            ["HBA_HUMAN", "dJ102G20.C1.1", "", "BROADO5"],
        )

    def test_gff2_parser_string(self):
        """Test the gff_parser works with a string filepath"""
        filepath = os.path.join("data/gff2_test.gff")
        for i, result in enumerate(gff_parser(filepath)):
            line = list(data_lines[i][1])
            self.assertEqual(result.pop("Attributes")["Info"], line.pop(8))
            self.assertEqual(set(result.values()), set(line))

    def test_gff2_parser_path(self):
        """Test the gff_parser works with a pathlib.Path filepath"""
        filepath = Path("data/gff2_test.gff")
        for i, result in enumerate(gff_parser(filepath)):
            line = list(data_lines[i][1])
            self.assertEqual(result.pop("Attributes")["Info"], line.pop(8))
            self.assertEqual(set(result.values()), set(line))

    def test_gff3_parser(self):
        """Test the gff_parser works on a gff3 file"""
        gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
        for i, result in enumerate(gff_parser(gff3_path)):
            self.assertEqual(len(result), 10)
        # 15 total lines, but 2 comments
        self.assertEqual(i + 1, 15 - 2)


if __name__ == "__main__":
    main()
