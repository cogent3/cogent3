"""Unit tests for GFF and related parsers.
"""
import os
import re

from io import StringIO
from pathlib import Path
from unittest import TestCase, main

import pytest

from cogent3.parse.gff import gff_parser, parse_attributes_gff2


DATA_DIR = Path(__file__).parent.parent / "data"

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
        for line, canned_result in data_lines:
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

    def test_custom_attr_func(self):
        """user provided attr parser"""
        gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
        for result in gff_parser(gff3_path, attribute_parser=lambda x, y: x):
            self.assertIsInstance(result["Attributes"], str)


@pytest.fixture
def multi_seqid_path(tmp_path):
    gff_path = DATA_DIR / "ensembl_sample.gff3"
    data = gff_path.read_text().splitlines()
    # add two new seqid's using the last 4 lines, twice
    change = data[-4:]
    data.extend(_modify_lines(change, "3"))
    data.extend(_modify_lines(change, "4"))
    outpath = tmp_path / "ensembl-edited.gff3"
    outpath.write_text("\n".join(data))
    return outpath


def _modify_lines(change, seqid: str) -> list:
    ident = re.compile(r"ENS[GT]\d+")
    result = []
    for line in change:
        # change all identifiers so each seqid has unique identifiers
        for match in ident.findall(line):
            line = line.replace(match, f"{match}{seqid}")
        line = line.split("\t")
        line[0] = seqid
        line = "\t".join(line)
        result.append(line)
    return result


@pytest.mark.parametrize("seqids", ("22", ("3",), ("3", "4")))
def test_seq_names(multi_seqid_path, seqids):
    got = {r["SeqID"] for r in gff_parser(multi_seqid_path, seqids=seqids)}
    expect = {seqids} if isinstance(seqids, str) else set(seqids)
    assert got == expect


def test_no_seq_names(multi_seqid_path):
    got = {r["SeqID"] for r in gff_parser(multi_seqid_path, seqids=None)}
    expect = {"22", "3", "4"}
    assert got == expect


def test_parse_field_spaces():
    path = DATA_DIR / "simple.gff"
    got = list(gff_parser(path))
    for record in got:
        for value in record.values():
            if isinstance(value, str):
                assert value.strip() == value, "should not have spaces!"
