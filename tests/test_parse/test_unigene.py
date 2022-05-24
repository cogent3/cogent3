#!/usr/bin/env python
"""Unit tests for unigene-specific classes
"""
from unittest import TestCase, main

from cogent3.parse.record_finder import GbFinder
from cogent3.parse.unigene import (
    LinesToUniGene,
    UniGeneProtSimRecord,
    UniGeneSeqRecord,
    _read_expression,
    _read_seq,
    _read_sts,
)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class unigeneTests(TestCase):
    """Tests toplevel functions."""

    def test_read_sts(self):
        """_read_sts should perform correct conversions"""
        self.assertEqual(
            _read_sts("ACC=RH128467 UNISTS=211775\n"),
            {"ACC": "RH128467", "UNISTS": "211775"},
        )

    def test_read_expression(self):
        """_read_expression should perform correct conversions"""
        self.assertEqual(
            _read_expression("embryo ; whole body ; mammary gland ; brain\n"),
            ["embryo", "whole body", "mammary gland", "brain"],
        )

    def test_read_seq(self):
        """_read_seq should perform correct conversions"""
        # reset the found fields, since we can't guarantee order of test
        # execution and it's persistent class data
        UniGeneSeqRecord.found_fields = {}
        self.assertEqual(
            _read_seq("ACC=BC025044.1\n"), UniGeneSeqRecord({"ACC": "BC025044.1"})
        )
        self.assertEqual(
            _read_seq(
                "ACC=AI842963.1; NID=g5477176; CLONE=UI-M-AO1-aem-f-10-0-UI; END=3'; LID=1944; SEQTYPE=EST; TRACE=158501677\n"
            ),
            UniGeneSeqRecord(
                {
                    "ACC": "AI842963.1",
                    "NID": "g5477176",
                    "CLONE": "UI-M-AO1-aem-f-10-0-UI",
                    "END": "3'",
                    "LID": "1944",
                    "SEQTYPE": "EST",
                    "TRACE": "158501677",
                }
            ),
        )

    def test_LinesToUniGene(self):
        """LinesToUniGene should give expected results on sample data"""
        fake_file = """ID          Mm.1
TITLE       S100 calcium binder
GENE        S100a10
CYTOBAND    3 41.7 cM
LOCUSLINK   20194
EXPRESS     embryo ; whole body ; mammary gland ; brain
CHROMOSOME  3
STS         ACC=RH128467 UNISTS=211775
STS         ACC=M16465 UNISTS= 178878
PROTSIM     ORG=Homo sapiens; PROTGI=107251; PROTID=pir:JC1139; PCT=91; ALN=97
PROTSIM     ORG=Mus musculus; PROTGI=116487; PROTID=sp:P08207; PCT=100; ALN=97
PROTSIM     ORG=Rattus norvegicus; PROTGI=116489; PROTID=sp:P05943; PCT=94; ALN=94
SCOUNT      5
SEQUENCE    ACC=BC025044.1; NID=g19263549; PID=g19263550; SEQTYPE=mRNA
SEQUENCE    ACC=AA471893.1; NID=g2199884; CLONE=IMAGE:872193; END=5'; LID=539; SEQTYPE=EST
SEQUENCE    ACC=AI842963.1; NID=g5477176; CLONE=UI-M-AO1-aem-f-10-0-UI; END=3'; LID=1944; SEQTYPE=EST; TRACE=158501677
SEQUENCE    ACC=CB595147.1; NID=g29513003; CLONE=IMAGE:30300703; END=5'; LID=12885; MGC=6677832; SEQTYPE=EST
SEQUENCE    ACC=BY144053.1; NID=g26280109; CLONE=L930184D22; END=5'; LID=12267; SEQTYPE=EST
//
ID          Mm.5
TITLE       homeo box A10
GENE        Hoxa10
CYTOBAND    6 26.33 cM
LOCUSLINK   15395
EXPRESS     kidney ; colon ; mammary gland
CHROMOSOME  6
PROTSIM     ORG=Caenorhabditis elegans; PROTGI=7510074; PROTID=pir:T31611; PCT=30; ALN=326
SCOUNT      1
SEQUENCE    ACC=AW990320.1; NID=g8185938; CLONE=IMAGE:1513482; END=5'; LID=1043; SEQTYPE=EST; TRACE=94472873
//
"""
        records = list(GbFinder(fake_file.split("\n")))
        self.assertEqual(len(records), 2)
        first, second = list(map(LinesToUniGene, records))
        self.assertEqual(first.ID, "Mm.1")
        self.assertEqual(first.TITLE, "S100 calcium binder")
        self.assertEqual(first.GENE, "S100a10")
        self.assertEqual(first.CYTOBAND, "3 41.7 cM")
        self.assertEqual(first.CHROMOSOME, "3")
        self.assertEqual(first.LOCUSLINK, 20194)
        self.assertEqual(
            first.EXPRESS, ["embryo", "whole body", "mammary gland", "brain"]
        )
        self.assertEqual(
            first.STS,
            [
                {"ACC": "RH128467", "UNISTS": "211775"},
                {"ACC": "M16465", "UNISTS": "178878"},
            ],
        )
        exp_prot_sim = list(
            map(
                UniGeneProtSimRecord,
                [
                    {
                        "ORG": "Homo sapiens",
                        "PROTGI": "107251",
                        "PROTID": "pir:JC1139",
                        "PCT": "91",
                        "ALN": "97",
                    },
                    {
                        "ORG": "Mus musculus",
                        "PROTGI": "116487",
                        "PROTID": "sp:P08207",
                        "PCT": "100",
                        "ALN": "97",
                    },
                    {
                        "ORG": "Rattus norvegicus",
                        "PROTGI": "116489",
                        "PROTID": "sp:P05943",
                        "PCT": "94",
                        "ALN": "94",
                    },
                ],
            )
        )
        for obs, exp in zip(first.PROTSIM, exp_prot_sim):
            self.assertEqual(obs, exp)
        self.assertEqual(first.SCOUNT, 5)
        exp_seqs = list(
            map(
                UniGeneSeqRecord,
                [
                    {
                        "ACC": "BC025044.1",
                        "NID": "g19263549",
                        "PID": "g19263550",
                        "SEQTYPE": "mRNA",
                    },
                    {
                        "ACC": "AA471893.1",
                        "NID": "g2199884",
                        "END": "5'",
                        "CLONE": "IMAGE:872193",
                        "LID": "539",
                        "SEQTYPE": "EST",
                    },
                    {
                        "ACC": "AI842963.1",
                        "NID": "g5477176",
                        "CLONE": "UI-M-AO1-aem-f-10-0-UI",
                        "END": "3'",
                        "LID": "1944",
                        "SEQTYPE": "EST",
                        "TRACE": "158501677",
                    },
                    {
                        "ACC": "CB595147.1",
                        "NID": "g29513003",
                        "CLONE": "IMAGE:30300703",
                        "END": "5'",
                        "LID": "12885",
                        "MGC": "6677832",
                        "SEQTYPE": "EST",
                    },
                    {
                        "ACC": "BY144053.1",
                        "NID": "g26280109",
                        "CLONE": "L930184D22",
                        "END": "5'",
                        "LID": "12267",
                        "SEQTYPE": "EST",
                    },
                ],
            )
        )
        for obs, exp in zip(first.SEQUENCE, exp_seqs):
            self.assertEqual(obs, exp)
        self.assertEqual(second.ID, "Mm.5")
        self.assertEqual(second.TITLE, "homeo box A10")
        self.assertEqual(second.GENE, "Hoxa10")
        self.assertEqual(second.CYTOBAND, "6 26.33 cM")
        self.assertEqual(second.LOCUSLINK, 15395)
        self.assertEqual(second.EXPRESS, ["kidney", "colon", "mammary gland"])
        self.assertEqual(second.CHROMOSOME, "6")
        self.assertEqual(
            second.PROTSIM,
            list(
                map(
                    UniGeneProtSimRecord,
                    [
                        {
                            "ORG": "Caenorhabditis elegans",
                            "PROTGI": "7510074",
                            "PROTID": "pir:T31611",
                            "PCT": "30",
                            "ALN": "326",
                        }
                    ],
                )
            ),
        )
        self.assertEqual(second.SCOUNT, 1)
        self.assertEqual(second.STS, [])
        self.assertEqual(
            second.SEQUENCE,
            list(
                map(
                    UniGeneSeqRecord,
                    [
                        {
                            "ACC": "AW990320.1",
                            "NID": "g8185938",
                            "CLONE": "IMAGE:1513482",
                            "END": "5'",
                            "LID": "1043",
                            "SEQTYPE": "EST",
                            "TRACE": "94472873",
                        }
                    ],
                )
            ),
        )

        # test that the synonym mapping works OK
        self.assertEqual(second.SequenceIds[0].NucleotideId, "g8185938")


if __name__ == "__main__":
    main()
