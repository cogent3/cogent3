# test_rdb.py
"""Unit test for RDB Parser"""

from unittest import TestCase

import cogent3
from cogent3.core.info import Info
from cogent3.parse.rdb import (
    InfoMaker,
    MinimalRdbParser,
    RdbParser,
    create_acceptable_sequence,
    is_seq_label,
)
from cogent3.parse.record import RecordError


# ruff: noqa: SIM905
class RdbTests(TestCase):
    """Tests for top-level functions in Rdb.py"""

    def test_is_seq_label(self):
        """is_seq_label should return True if a line starts with 'seq:'"""
        seq = "seq:this is a sequence line"
        not_seq = "this is not a sequence line"
        still_not_seq = "this seq: is still not a sequence line"
        assert is_seq_label(seq) is True
        assert is_seq_label(not_seq) is False
        assert is_seq_label(still_not_seq) is False

    def test_create_acceptable_sequence(self):
        """create_acceptable_sequence: should handle 'o' and sec. struct"""
        f = create_acceptable_sequence
        # should keep any char accepted by RNA.alphabet.degen_gapped
        s = "UCAG---NRYBDHKMNSRWVY?"
        assert f(s) == s
        # should replace 'o' by '?'
        s = "UCAG-oo-ACGU"
        assert f(s) == "UCAG-??-ACGU"
        # should strip out secondary info
        s = "{UC^AG-[oo]-A(CG)U}"
        assert f(s) == "UCAG-??-ACGU"
        # should leave other chars untouched
        s = "XYZ1234"
        assert f(s) == "XYZ1234"


class InfoMakerTests(TestCase):
    """Tests for the Constructor InfoMaker. Should return an Info object"""

    def test_empty(self):
        """InfoMaker: should return empty Info from empty header"""
        empty_header = []
        obs = InfoMaker(empty_header)
        exp = Info()
        assert obs == exp

    def test_full(self):
        """InfoMaker should return Info object with name, value pairs"""
        test_header = [
            "acc: X3402",
            "abc:1",
            "mty: ssu",
            "seq: Mit. X3402",
            "",
            "nonsense",
            ":no_name",
        ]
        obs = InfoMaker(test_header)
        exp = Info()
        exp.rRNA = "X3402"
        exp.abc = "1"
        exp.Species = "Mit. X3402"
        exp.Gene = "ssu"
        assert obs == exp


class GenericRdbTest(TestCase):
    "SetUp data for all Rdb parsers"

    def setUp(self):
        self.empty = []
        self.labels = "mty:ssu\nseq:bac\n//\nttl:joe\nseq:mit\n//".split("\n")
        self.nolabels = "ACGUAGCUAGCUAC\nGCUGCAUCG\nAUCG\n//".split("\n")
        self.oneseq = "seq:H.Sapiens\nAGUCAUCUAGAUHCAUHC\n//".split("\n")
        self.multiline = "seq:H.Sapiens\nAGUCAUUAG\nAUHCAUHC\n//".split("\n")
        self.threeseq = "seq:bac\nAGU\n//\nseq:mit\nACU\n//\nseq:pla\nAAA\n//".split(
            "\n",
        )
        self.twogood = "seq:bac\n//\nseq:mit\nACU\n//\nseq:pla\nAAA\n//".split("\n")
        self.oneX = "seq:bac\nX\n//\nseq:mit\nACT\n//\nseq:pla\nAAA\n//".split("\n")
        self.strange = "seq:bac\nACGUXxAaKkoo---*\n//".split("\n")


class MinimalRdbParserTests(GenericRdbTest):
    """Tests of MinimalRdbParser: returns (headerLines,sequence) tuples"""

    def test_empty(self):
        """MinimalRdbParser should return empty list from file w/o seqs"""
        assert list(MinimalRdbParser(self.empty)) == []
        assert list(MinimalRdbParser(self.nolabels, strict=False)) == []
        self.assertRaises(RecordError, list, MinimalRdbParser(self.nolabels))

    def test_only_labels(self):
        """MinimalRdbParser should return empty list from file w/o seqs"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list, MinimalRdbParser(self.labels, strict=True))
        # if not strict, should skip the records
        assert list(MinimalRdbParser(self.labels, strict=False)) == []

    def test_only_sequences(self):
        """MinimalRdbParser should return empty list form file w/o lables"""
        # should fail if strict (the default)
        self.assertRaises(
            RecordError,
            list,
            MinimalRdbParser(self.nolabels, strict=True),
        )
        # if not strict, should skip the records
        assert list(MinimalRdbParser(self.nolabels, strict=False)) == []

    def test_single(self):
        """MinimalRdbParser should read single record as (header,seq) tuple"""
        res = list(MinimalRdbParser(self.oneseq))
        assert len(res) == 1
        first = res[0]
        assert first == (["seq:H.Sapiens"], "AGUCAUCUAGAUHCAUHC")

        res = list(MinimalRdbParser(self.multiline))
        assert len(res) == 1
        first = res[0]
        assert first == (["seq:H.Sapiens"], "AGUCAUUAGAUHCAUHC")

    def test_multiple(self):
        """MinimalRdbParser should read multiple record correctly"""
        res = list(MinimalRdbParser(self.threeseq))
        assert len(res) == 3
        a, b, c = res
        assert a == (["seq:bac"], "AGU")
        assert b == (["seq:mit"], "ACU")
        assert c == (["seq:pla"], "AAA")

    def test_multiple_bad(self):
        """MinimalRdbParser should complain or skip bad records"""
        self.assertRaises(RecordError, list, MinimalRdbParser(self.twogood))
        f = list(MinimalRdbParser(self.twogood, strict=False))
        assert len(f) == 2
        a, b = f
        assert a == (["seq:mit"], "ACU")
        assert b == (["seq:pla"], "AAA")

    def test_strange(self):
        """MRP: handle strange char. according to constr. and strip off '*'"""
        f = list(MinimalRdbParser(self.strange))
        obs = f[0]
        exp = (["seq:bac"], "ACGUXxAaKkoo---")
        assert obs == exp


class RdbParserTests(GenericRdbTest):
    """Tests for the RdbParser. Should return Sequence objects"""

    def test_empty(self):
        """RdbParser should return empty list from 'file' w/o labels"""
        assert list(RdbParser(self.empty)) == []
        assert list(RdbParser(self.nolabels, strict=False)) == []
        self.assertRaises(RecordError, list, RdbParser(self.nolabels))

    def test_only_labels(self):
        """RdbParser should return empty list from file w/o seqs"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list, RdbParser(self.labels, strict=True))
        # if not strict, should skip the records
        assert list(RdbParser(self.labels, strict=False)) == []

    def test_only_sequences(self):
        """RdbParser should return empty list form file w/o lables"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list, RdbParser(self.nolabels, strict=True))
        # if not strict, should skip the records
        assert list(RdbParser(self.nolabels, strict=False)) == []

    def test_single(self):
        """RdbParser should read single record as (header,seq) tuple"""
        res = list(RdbParser(self.oneseq))
        assert len(res) == 1
        first = res[0]
        assert first == cogent3.make_seq(seq="AGUCAUCUAGAUHCAUHC", moltype="rna")
        assert first.info == Info(
            {"Species": "H.Sapiens", "OriginalSeq": "AGUCAUCUAGAUHCAUHC"},
        )

        res = list(RdbParser(self.multiline))
        assert len(res) == 1
        first = res[0]
        assert first == cogent3.make_seq(seq="AGUCAUUAGAUHCAUHC", moltype="rna")
        assert first.info == Info(
            {"Species": "H.Sapiens", "OriginalSeq": "AGUCAUUAGAUHCAUHC"},
        )

    def test_multiple_constructor_bad(self):
        """RdbParser should complain or skip bad records w/ constructor"""

        self.assertRaises(RecordError, list, RdbParser(self.oneX))
        f = list(RdbParser(self.oneX, strict=False))
        assert len(f) == 2
        a, b = f

        assert a == "ACU"
        assert a.info == Info({"Species": "mit", "OriginalSeq": "ACT"})
        assert b == "AAA"
        assert b.info == Info({"Species": "pla", "OriginalSeq": "AAA"})

    def test_full(self):
        """RdbParser: full data, valid and invalid"""
        # when only good record, should work independent of strict
        r1 = cogent3.make_seq(
            "-??GG-UGAA--CGCU---ACGU-N???---",
            info=Info(
                {
                    "Species": "unidentified Thermus OPB AF027020",
                    "Refs": {"rRNA": ["AF027020"]},
                    "OriginalSeq": "-o[oGG-U{G}AA--C^GC]U---ACGU-Nooo---",
                },
            ),
            moltype="rna",
        )
        r2 = cogent3.make_seq(
            "---CGAUCG--UAUACG-N???-",
            info=Info(
                {
                    "Species": "Thermus silvanus X84211",
                    "Refs": {"rRNA": ["X84211"]},
                    "OriginalSeq": "---CGAU[C(G){--UA}U]ACG-Nooo-",
                },
            ),
            moltype="rna",
        )
        obs = list(RdbParser(RDB_LINES_ONLY_GOOD.split("\n"), strict=True))
        assert len(obs) == 2
        assert obs[0] == r1
        assert str(obs[0]) == str(r1)
        assert obs[0].info == r1.info
        assert obs[1] == r2
        assert str(obs[1]) == str(r2)
        assert obs[1].info == r2.info

        obs = list(RdbParser(RDB_LINES_ONLY_GOOD.split("\n"), strict=False))
        assert len(obs) == 2
        assert obs[0] == r1
        assert str(obs[0]) == str(r1)
        assert obs[0].info == r1.info

        # when strict, should raise error on invalid record
        f = RdbParser(RDB_LINES_GOOD_BAD.split("\n"), strict=True)
        self.assertRaises(RecordError, list, f)
        # when not strict, malicious record is skipped
        obs = list(RdbParser(RDB_LINES_GOOD_BAD.split("\n"), strict=False))
        assert len(obs) == 2
        assert obs[0] == r1
        assert str(obs[0]) == str(r1)
        assert obs[0].info == r1.info
        assert obs[1] == r2
        assert str(obs[1]) == str(r2)
        assert obs[1].info == r2.info


RDB_LINES_ONLY_GOOD = """acc:AF027020
seq: unidentified Thermus OPB AF027020
-o[oGG-U{G}AA--C^GC]U---ACGU-Nooo---
//
acc:X84211
seq: Thermus silvanus X84211
---CGAU[C(G){--UA}U]ACG-Nooo-
//
"""

RDB_LINES_GOOD_BAD = """acc:AF027020
seq: unidentified Thermus OPB AF027020
-o[oGG-U{G}AA--C^GC]U---ACGU-Nooo---
//
acc:ABC123
seq: E. coli
---ACGU-Nooo-RYXQ-
//
acc:X84211
seq: Thermus silvanus X84211
---CGAU[C(G){--UA}U]ACG-Nooo-
//
"""
