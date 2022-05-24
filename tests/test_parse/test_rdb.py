#!/usr/bin/env python
# test_rdb.py
"""Unit test for RDB Parser
"""
from unittest import TestCase, main

from cogent3.core.info import Info
from cogent3.core.sequence import DnaSequence, RnaSequence, Sequence
from cogent3.parse.rdb import (
    InfoMaker,
    MinimalRdbParser,
    RdbParser,
    create_acceptable_sequence,
    is_seq_label,
)
from cogent3.parse.record import RecordError


__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"


class RdbTests(TestCase):
    """Tests for top-level functions in Rdb.py"""

    def test_is_seq_label(self):
        """is_seq_label should return True if a line starts with 'seq:'"""
        seq = "seq:this is a sequence line"
        not_seq = "this is not a sequence line"
        still_not_seq = "this seq: is still not a sequence line"
        self.assertEqual(is_seq_label(seq), True)
        self.assertEqual(is_seq_label(not_seq), False)
        self.assertEqual(is_seq_label(still_not_seq), False)

    def test_create_acceptable_sequence(self):
        """create_acceptable_sequence: should handle 'o' and sec. struct"""
        f = create_acceptable_sequence
        # should keep any char accepted by RNA.alphabet.degen_gapped
        s = "UCAG---NRYBDHKMNSRWVY?"
        self.assertEqual(f(s), s)
        # should replace 'o' by '?'
        s = "UCAG-oo-ACGU"
        self.assertEqual(f(s), "UCAG-??-ACGU")
        # should strip out secondary info
        s = "{UC^AG-[oo]-A(CG)U}"
        self.assertEqual(f(s), "UCAG-??-ACGU")
        # should leave other chars untouched
        s = "XYZ1234"
        self.assertEqual(f(s), "XYZ1234")


class InfoMakerTests(TestCase):
    """Tests for the Constructor InfoMaker. Should return an Info object"""

    def test_empty(self):
        """InfoMaker: should return empty Info from empty header"""
        empty_header = []
        obs = InfoMaker(empty_header)
        exp = Info()
        self.assertEqual(obs, exp)

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
        self.assertEqual(obs, exp)


class GenericRdbTest(TestCase):
    "SetUp data for all Rdb parsers" ""

    def setUp(self):
        self.empty = []
        self.labels = "mty:ssu\nseq:bac\n//\nttl:joe\nseq:mit\n//".split("\n")
        self.nolabels = "ACGUAGCUAGCUAC\nGCUGCAUCG\nAUCG\n//".split("\n")
        self.oneseq = "seq:H.Sapiens\nAGUCAUCUAGAUHCAUHC\n//".split("\n")
        self.multiline = "seq:H.Sapiens\nAGUCAUUAG\nAUHCAUHC\n//".split("\n")
        self.threeseq = "seq:bac\nAGU\n//\nseq:mit\nACU\n//\nseq:pla\nAAA\n//".split(
            "\n"
        )
        self.twogood = "seq:bac\n//\nseq:mit\nACU\n//\nseq:pla\nAAA\n//".split("\n")
        self.oneX = "seq:bac\nX\n//\nseq:mit\nACT\n//\nseq:pla\nAAA\n//".split("\n")
        self.strange = "seq:bac\nACGUXxAaKkoo---*\n//".split("\n")


class MinimalRdbParserTests(GenericRdbTest):
    """Tests of MinimalRdbParser: returns (headerLines,sequence) tuples"""

    def test_empty(self):
        """MinimalRdbParser should return empty list from file w/o seqs"""
        self.assertEqual(list(MinimalRdbParser(self.empty)), [])
        self.assertEqual(list(MinimalRdbParser(self.nolabels, strict=False)), [])
        self.assertRaises(RecordError, list, MinimalRdbParser(self.nolabels))

    def test_only_labels(self):
        """MinimalRdbParser should return empty list from file w/o seqs"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list, MinimalRdbParser(self.labels, strict=True))
        # if not strict, should skip the records
        self.assertEqual(list(MinimalRdbParser(self.labels, strict=False)), [])

    def test_only_sequences(self):
        """MinimalRdbParser should return empty list form file w/o lables"""
        # should fail if strict (the default)
        self.assertRaises(
            RecordError, list, MinimalRdbParser(self.nolabels, strict=True)
        )
        # if not strict, should skip the records
        self.assertEqual(list(MinimalRdbParser(self.nolabels, strict=False)), [])

    def test_single(self):
        """MinimalRdbParser should read single record as (header,seq) tuple"""
        res = list(MinimalRdbParser(self.oneseq))
        self.assertEqual(len(res), 1)
        first = res[0]
        self.assertEqual(first, (["seq:H.Sapiens"], "AGUCAUCUAGAUHCAUHC"))

        res = list(MinimalRdbParser(self.multiline))
        self.assertEqual(len(res), 1)
        first = res[0]
        self.assertEqual(first, (["seq:H.Sapiens"], "AGUCAUUAGAUHCAUHC"))

    def test_multiple(self):
        """MinimalRdbParser should read multiple record correctly"""
        res = list(MinimalRdbParser(self.threeseq))
        self.assertEqual(len(res), 3)
        a, b, c = res
        self.assertEqual(a, (["seq:bac"], "AGU"))
        self.assertEqual(b, (["seq:mit"], "ACU"))
        self.assertEqual(c, (["seq:pla"], "AAA"))

    def test_multiple_bad(self):
        """MinimalRdbParser should complain or skip bad records"""
        self.assertRaises(RecordError, list, MinimalRdbParser(self.twogood))
        f = list(MinimalRdbParser(self.twogood, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f
        self.assertEqual(a, (["seq:mit"], "ACU"))
        self.assertEqual(b, (["seq:pla"], "AAA"))

    def test_strange(self):
        """MRP: handle strange char. according to constr. and strip off '*'"""
        f = list(MinimalRdbParser(self.strange))
        obs = f[0]
        exp = (["seq:bac"], "ACGUXxAaKkoo---")
        self.assertEqual(obs, exp)


class RdbParserTests(GenericRdbTest):
    """Tests for the RdbParser. Should return Sequence objects"""

    def test_empty(self):
        """RdbParser should return empty list from 'file' w/o labels"""
        self.assertEqual(list(RdbParser(self.empty)), [])
        self.assertEqual(list(RdbParser(self.nolabels, strict=False)), [])
        self.assertRaises(RecordError, list, RdbParser(self.nolabels))

    def test_only_labels(self):
        """RdbParser should return empty list from file w/o seqs"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list, RdbParser(self.labels, strict=True))
        # if not strict, should skip the records
        self.assertEqual(list(RdbParser(self.labels, strict=False)), [])

    def test_only_sequences(self):
        """RdbParser should return empty list form file w/o lables"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list, RdbParser(self.nolabels, strict=True))
        # if not strict, should skip the records
        self.assertEqual(list(RdbParser(self.nolabels, strict=False)), [])

    def test_single(self):
        """RdbParser should read single record as (header,seq) tuple"""
        res = list(RdbParser(self.oneseq))
        self.assertEqual(len(res), 1)
        first = res[0]
        self.assertEqual(first, Sequence("AGUCAUCUAGAUHCAUHC"))
        self.assertEqual(
            first.info,
            Info({"Species": "H.Sapiens", "OriginalSeq": "AGUCAUCUAGAUHCAUHC"}),
        )

        res = list(RdbParser(self.multiline))
        self.assertEqual(len(res), 1)
        first = res[0]
        self.assertEqual(first, Sequence("AGUCAUUAGAUHCAUHC"))
        self.assertEqual(
            first.info,
            Info({"Species": "H.Sapiens", "OriginalSeq": "AGUCAUUAGAUHCAUHC"}),
        )

    def test_single_constructor(self):
        """RdbParser should use constructors if supplied"""
        to_dna = lambda x, info: DnaSequence(str(x).replace("U", "T"), info=info)
        f = list(RdbParser(self.oneseq, to_dna))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, "AGTCATCTAGATHCATHC")
        self.assertEqual(
            a.info, Info({"Species": "H.Sapiens", "OriginalSeq": "AGUCAUCUAGAUHCAUHC"})
        )

        def alternativeConstr(header_lines):
            info = Info()
            for line in header_lines:
                all = line.strip().split(":", 1)
                # strip out empty lines, lines without name, lines without
                # colon
                if not all[0] or len(all) != 2:
                    continue
                name = all[0].upper()
                value = all[1].strip().upper()
                info[name] = value
            return info

        f = list(RdbParser(self.oneseq, to_dna, alternativeConstr))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, "AGTCATCTAGATHCATHC")
        exp_info = Info(
            {"OriginalSeq": "AGUCAUCUAGAUHCAUHC", "Refs": {}, "SEQ": "H.SAPIENS"}
        )
        self.assertEqual(
            a.info,
            Info({"OriginalSeq": "AGUCAUCUAGAUHCAUHC", "Refs": {}, "SEQ": "H.SAPIENS"}),
        )

    def test_multiple_constructor_bad(self):
        """RdbParser should complain or skip bad records w/ constructor"""

        def dnastrict(x, **kwargs):
            try:
                return DnaSequence(x, **kwargs)
            except Exception:
                raise RecordError("Could not convert sequence")

        self.assertRaises(RecordError, list, RdbParser(self.oneX, dnastrict))
        f = list(RdbParser(self.oneX, dnastrict, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f

        self.assertEqual(a, "ACT")
        self.assertEqual(a.info, Info({"Species": "mit", "OriginalSeq": "ACT"}))
        self.assertEqual(b, "AAA")
        self.assertEqual(b.info, Info({"Species": "pla", "OriginalSeq": "AAA"}))

    def test_full(self):
        """RdbParser: full data, valid and invalid"""
        # when only good record, should work independent of strict
        r1 = RnaSequence(
            "-??GG-UGAA--CGCU---ACGU-N???---",
            info=Info(
                {
                    "Species": "unidentified Thermus OPB AF027020",
                    "Refs": {"rRNA": ["AF027020"]},
                    "OriginalSeq": "-o[oGG-U{G}AA--C^GC]U---ACGU-Nooo---",
                }
            ),
        )
        r2 = RnaSequence(
            "---CGAUCG--UAUACG-N???-",
            info=Info(
                {
                    "Species": "Thermus silvanus X84211",
                    "Refs": {"rRNA": ["X84211"]},
                    "OriginalSeq": "---CGAU[C(G){--UA}U]ACG-Nooo-",
                }
            ),
        )
        obs = list(RdbParser(RDB_LINES_ONLY_GOOD.split("\n"), strict=True))
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0], r1)
        self.assertEqual(str(obs[0]), str(r1))
        self.assertEqual(obs[0].info, r1.info)
        self.assertEqual(obs[1], r2)
        self.assertEqual(str(obs[1]), str(r2))
        self.assertEqual(obs[1].info, r2.info)

        obs = list(RdbParser(RDB_LINES_ONLY_GOOD.split("\n"), strict=False))
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0], r1)
        self.assertEqual(str(obs[0]), str(r1))
        self.assertEqual(obs[0].info, r1.info)

        # when strict, should raise error on invalid record
        f = RdbParser(RDB_LINES_GOOD_BAD.split("\n"), strict=True)
        self.assertRaises(RecordError, list, f)
        # when not strict, malicious record is skipped
        obs = list(RdbParser(RDB_LINES_GOOD_BAD.split("\n"), strict=False))
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0], r1)
        self.assertEqual(str(obs[0]), str(r1))
        self.assertEqual(obs[0].info, r1.info)
        self.assertEqual(obs[1], r2)
        self.assertEqual(str(obs[1]), str(r2))
        self.assertEqual(obs[1].info, r2.info)


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

if __name__ == "__main__":
    main()
