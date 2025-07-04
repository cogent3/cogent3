"""Unit tests for FASTA and related parsers."""

import io
import os
import pathlib
from unittest import TestCase

import numpy
import pytest

from cogent3.core.info import Info
from cogent3.parse.fasta import (
    LabelParser,
    MinimalFastaParser,
    RichLabel,
    iter_fasta_records,
)
from cogent3.parse.record import RecordError

# ruff: noqa: SIM905

base_path = os.path.dirname(os.path.dirname(__file__))
data_path = os.path.join(base_path, "data")


class GenericFastaTest(TestCase):
    """Setup data for all the various FASTA parsers."""

    def setUp(self):
        """standard files"""
        self.labels = ">abc\n>def\n>ghi\n".split("\n")
        self.oneseq = ">abc\nUCAG\n".split("\n")
        self.multiline = ">xyz\nUUUU\nCC\nAAAAA\nG".split("\n")
        self.threeseq = ">123\na\n> \t abc  \t \ncag\ngac\n>456\nc\ng".split("\n")
        self.twogood = ">123\n\n> \t abc  \t \ncag\ngac\n>456\nc\ng".split("\n")
        self.oneX = ">123\nX\n> \t abc  \t \ncag\ngac\n>456\nc\ng".split("\n")
        self.nolabels = "GJ>DSJGSJDF\nSFHKLDFS>jkfs\n".split("\n")
        self.empty = []


class MinimalFastaParserTests(GenericFastaTest):
    """Tests of MinimalFastaParser: returns (label, seq) tuples."""

    def test_single(self):
        """MinimalFastaParser should read single record as (label, seq) tuple"""
        f = list(MinimalFastaParser(self.oneseq))
        assert len(f) == 1
        a = f[0]
        assert a == ("abc", "UCAG")

        f = list(MinimalFastaParser(self.multiline))
        assert len(f) == 1
        a = f[0]
        assert a == ("xyz", "UUUUCCAAAAAG")

    def test_multiple(self):
        """MinimalFastaParser should read multiline records correctly"""
        f = list(MinimalFastaParser(self.threeseq))
        assert len(f) == 3
        a, b, c = f
        assert a == ("123", "a")
        assert b == ("abc", "caggac")
        assert c == ("456", "cg")

    def test_parser_from_file(self):
        """passing path should work"""
        path = os.path.join(data_path, "brca1.fasta")
        seqs = dict(p for p in MinimalFastaParser(path))
        assert "Human" in seqs


class LabelParsingTest(TestCase):
    """Test generic fasta label parsing"""

    def test_rich_label(self):
        """rich label correctly constructs label strings"""
        # labels should be equal based on the result of applying their
        # attributes to their string template
        k = RichLabel(Info(species="rat"), "%(species)s")
        l = RichLabel(Info(species="rat", seq_id="xy5"), "%(species)s")
        assert k == l

        # labels should construct from Info components correctly
        k = RichLabel(Info(species="rat", seq_id="xy5"), "%(seq_id)s:%(species)s")
        assert k == "xy5:rat"
        k = RichLabel(Info(species="rat", seq_id="xy5"), "%(species)s:%(seq_id)s")
        assert k == "rat:xy5"

        # extra components should be ignored
        k = RichLabel(Info(species="rat", seq_id="xy5"), "%(species)s")
        assert k == "rat"

        # the label should have Info object
        assert k.info.species == "rat"
        assert k.info.seq_id == "xy5"

        # label should be constructable just like a normal string
        assert RichLabel("a") == "a"

    def test_label_parser(self):
        """label parser factory function cope with mixed structure labels"""
        # the label parser factory function should correctly handle label lines
        # with mixed separators
        make = LabelParser(
            "%(species)s:%(accession)s",
            [[0, "accession", str], [2, "species", str]],
            split_with=": ",
        )
        for label, expect in [
            (">abcd:human:misc", "misc:abcd"),
            ("abcd:human:misc", "misc:abcd"),
            (">abcd:Human misc", "misc:abcd"),
            (">abcd Human:misc", "misc:abcd"),
            (">abcd:Human misc", "misc:abcd"),
        ]:
            assert make(label) == expect

        # should raise an assertion error if template doesn't match at least
        # one field name
        self.assertRaises(
            AssertionError,
            LabelParser,
            "%s:%s",
            [[0, "accession", str], [2, "species", str]],
            split_with=": ",
        )


def test_empty():
    assert not list(MinimalFastaParser([]))


def test_missing_labels():
    nolabels = "GJ>DSJGSJDF\nSFHKLDFS>jkfs\n".split("\n")
    with pytest.raises(RecordError):
        list(MinimalFastaParser(nolabels, strict=True))


def test_no_labels_strict():
    labels = ">abc\n>def\n>ghi\n".split("\n")
    with pytest.raises(RecordError):
        list(MinimalFastaParser(labels, strict=True))


def test_no_labels():
    """MinimalFastaParser should return empty list from file w/o seqs"""
    labels = ">abc\n>def\n>ghi\n".split("\n")
    # if not strict, should skip the records
    got = {l: str(v) for l, v in MinimalFastaParser(labels, strict=False)}
    assert not got


def test_multiple_bad_strict():
    """MinimalFastaParser should complain or skip bad records"""
    twogood = ">123\n\n> \t abc  \t \ncag\ngac\n>456\nc\ng".split("\n")
    with pytest.raises(RecordError):
        list(MinimalFastaParser(twogood, strict=True))


def test_multiple_bad_not_strict():
    twogood = ">123\n\n> \t abc  \t \ncag\ngac\n>456\nc\ng".split("\n")
    f = list(MinimalFastaParser(twogood, strict=False))
    assert len(f) == 2
    expect = [("abc", "caggac"), ("456", "cg")]
    assert f == expect


def test_seq_startswith_gt_bracket_fail():
    """as seq that starts with > will fail under strict"""
    oneseq_w_gt = ">abc\n>CAG\n".split("\n")
    with pytest.raises(RecordError):
        list(MinimalFastaParser(oneseq_w_gt))


@pytest.mark.parametrize("sep", [" ", "\t"])
@pytest.mark.parametrize("strict", [False, True])
def test_fasta_with_spaces(strict, sep):
    data = [">A", sep.join(("gaaaa", "tgatt")), ">B", sep.join(("tttga", "gcagg"))]
    got = dict(MinimalFastaParser(data, strict=strict))
    assert got == {"A": "gaaaatgatt", "B": "tttgagcagg"}


def test_iter_fasta_records(DATA_DIR):
    path = DATA_DIR / "brca1.fasta"
    got = dict(iter_fasta_records(path))
    assert len(got) == 55
    assert "LesserEle" in got
    assert got["LesserEle"][:10] == "TGTGGCACAG"
    assert got["LesserEle"][-10:] == "ATGTA-----"


@pytest.fixture(params=(pathlib.Path, str, io.FileIO))
def fasta_path(DATA_DIR, request):
    path = DATA_DIR / "c_elegans_WS199_dna_shortened.fasta"
    if request.param is not io.FileIO:
        yield request.param(path)
    else:
        handle = path.open()
        yield handle
        handle.close()


def test_iter_fasta_records_path_types(fasta_path):
    got = dict(iter_fasta_records(fasta_path))
    assert len(got) == 1


@pytest.fixture(params=("gz", "xz"))
def fa_gz(DATA_DIR, tmp_path, request):
    from cogent3 import open_

    path = DATA_DIR / "brca1.fasta"
    data = path.read_bytes()
    outpath = tmp_path / f"brca1.fasta.{request.param}"
    with open_(outpath, "wb") as out:
        out.write(data)
    return outpath.as_uri()


def test_iter_fasta_records_compressed_file_uri(fa_gz):
    got = dict(iter_fasta_records(str(fa_gz)))
    assert len(got) == 55


def test_iter_fasta_records_invalid():
    with pytest.raises(TypeError):
        list(iter_fasta_records(numpy.array([">abcd", "ACGG"])))
