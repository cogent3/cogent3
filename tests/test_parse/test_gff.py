"""Unit tests for GFF and related parsers."""

import io
import os
import pathlib
import re
from io import StringIO
from pathlib import Path
from unittest import TestCase

import pytest
from cogent3.parse.gff import (
    GffRecord,
    gff_parser,
    is_gff3,
    parse_attributes_gff2,
)

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

    def test_parse_attributes_gff2(self):
        """Test the parse_attributes_gff2 method"""
        self.assertEqual(
            [parse_attributes_gff2(x[1][8])["ID"] for x in data_lines],
            ["HBA_HUMAN", "dJ102G20.C1.1", "", "BROADO5"],
        )

    def test_custom_attr_func(self):
        """user provided attr parser"""
        gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
        for result in gff_parser(gff3_path, attribute_parser=lambda x: x):
            self.assertIsInstance(result["Attributes"], str)


@pytest.fixture
def multi_seqid_path(DATA_DIR, tmp_path):
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


def test_parse_field_spaces(DATA_DIR):
    path = DATA_DIR / "simple.gff"
    got = list(gff_parser(path))
    for record in got:
        for attr in dir(record):
            if attr.startswith("_"):
                continue
            value = getattr(record, attr)
            if isinstance(value, str):
                assert value.strip() == value, f"{attr} should not have spaces!"


@pytest.mark.parametrize("line,canned_result", data_lines)
def test_gff_parser_data(line, canned_result):
    """Test gff_parser with valid data lines"""
    result = next(gff_parser(StringIO(line))).to_dict()
    canned_result = list(canned_result)
    assert result["attrs"]["Info"] == canned_result.pop(8)
    result.pop("attrs")
    assert set(result.values()) == set(canned_result)


def test_gff2_parser_path():
    """Test the gff_parser works with a pathlib.Path filepath"""
    filepath = Path("data/gff2_test.gff")
    for i, result in enumerate(gff_parser(filepath)):
        result = result.to_dict()
        line = list(data_lines[i][1])
        assert result.pop("attrs")["Info"] == line.pop(8)
        assert set(result.values()) == set(line)


def test_gff2_parser_string():
    """Test the gff_parser works with a string filepath"""
    filepath = os.path.join("data/gff2_test.gff")
    for i, result in enumerate(gff_parser(filepath)):
        result = result.to_dict()
        line = list(data_lines[i][1])
        assert result.pop("attrs")["Info"] == line.pop(8)
        assert set(result.values()) == set(line)


def test_gff3_parser():
    """Test the gff_parser works on a gff3 file"""
    gff3_path = os.path.join("data/c_elegans_WS199_shortened_gff.gff3")
    records = list(gff_parser(gff3_path))
    lengths = {len(result.to_dict()) for result in records}
    # add 3 for name, parent_id and comments
    assert lengths == {10 + 3}
    # 15 total lines, but 2 comments
    assert len(records) == 15 - 2


def test_gff_parser_headers():
    """Test gff_parser with valid data headers"""
    data = "".join([x[0] for x in data_lines])
    for header in headers:
        result = list(gff_parser(StringIO(header + data)))
        lines = [(x[0], list(x[1])) for x in data_lines]
        assert [l.attrs["Info"] for l in result] == [x[1].pop(8) for x in lines]
        expect = [set(x[1]) for x in lines]
        got = []
        for l in result:
            l = l.to_dict()
            l.pop("attrs")
            got.append(set(l.values()))
        assert got == expect


def test_is_gff3_invalid():
    with pytest.raises(TypeError):
        _ = is_gff3(b"blah")


def make_gff_text():
    return headers[0] + "".join([l for l, _ in data_lines])


@pytest.mark.parametrize(
    "source",
    (
        data_lines,
        io.StringIO(make_gff_text()),
        pathlib.Path("data/c_elegans_WS199_shortened_gff.gff3"),
        "data/c_elegans_WS199_shortened_gff.gff3",
    ),
)
def test_is_gfg3(source):
    assert isinstance(is_gff3(source), bool)


def test_repr_gff3_record():
    rec = GffRecord(name="1")
    assert "name='1'" in repr(rec)


def test_gff3_record_get():
    attrs = "some text"
    rec = GffRecord(name="1", attrs=attrs)
    assert rec.get("attributes") == attrs


def test_gff3_record_set():
    attrs = "some text"
    rec = GffRecord(name="1")
    rec.attrs = attrs
    assert rec["attributes"] == attrs


def test_gff3_record_update():
    attrs = "some text"
    rec = GffRecord(name="1")
    rec.update({"attrs": attrs})
    assert rec.get("attributes") == attrs


@pytest.fixture
def worm_path(DATA_DIR):
    return DATA_DIR / "c_elegans_WS199_shortened_gff.gff3"


@pytest.mark.parametrize("path_type", (str, pathlib.Path))
def test_gff_parser_path_types(worm_path, path_type):
    got = list(gff_parser(path_type(worm_path), gff3=True))
    assert len(got) == 13


@pytest.mark.parametrize(
    "obj", (make_gff_text().splitlines(), io.StringIO(make_gff_text()))
)
def test_gff_parser_obj_types(obj):
    got = list(gff_parser(obj))
    assert len(got) == 4


def test_gff_parser_make_record_override_attr_parser(worm_path):
    got = list(
        gff_parser(
            worm_path, attribute_parser=lambda x: x.split(";"), make_record=GffRecord
        )
    )
    assert len(got) == 13
    # no split operation
    assert isinstance(got[0].attrs, str)
