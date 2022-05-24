from unittest import TestCase, main

from cogent3.parse.blast import (
    FastacmdTaxonomyParser,
    GenericBlastParser9,
    LastProteinIds9,
    PsiBlastFinder,
    PsiBlastParser9,
    PsiBlastQueryFinder,
    PsiBlastTableParser,
    QMEBlast9,
    QMEPsiBlast9,
    TableToValues,
    fastacmd_taxonomy_splitter,
    is_blast_junk,
    is_blat_junk,
    iter_finder,
    iteration_set_finder,
    make_label,
    query_finder,
)


__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Micah Hamady", "Rob Knight"]
__license__ = "GPL"
__version__ = "2022.5.25a1"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Production"

from numpy.testing import assert_allclose, assert_equal


class BlastTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Define some standard data"""
        self.rec = """# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-06	52.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 2
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-08	59.0
ece:Z4181	sfl:CP0138	33.98	103	57	2	8	110	6	97	6e-06	50.5
ece:Z4181	spt:SPA2730	37.50	72	45	0	39	110	30	101	1e-05	49.8
ece:Z4181	sec:SC2804	37.50	72	45	0	39	110	30	101	1e-05	49.8
ece:Z4181	stm:STM2872	37.50	72	45	0	39	110	30	101	1e-05	49.8""".split(
            "\n"
        )

        self.rec2 = """# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-06	52.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 2
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-08	59.0
ece:Z4181	sfl:CP0138	33.98	103	57	2	8	110	6	97	6e-06	50.5
ece:Z4181	spt:SPA2730	37.50	72	45	0	39	110	30	101	1e-05	49.8
ece:Z4181	sec:SC2804	37.50	72	45	0	39	110	30	101	1e-05	49.8
ece:Z4181	stm:STM2872	37.50	72	45	0	39	110	30	101	1e-05	49.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4182
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4182	ece:Z4182	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4182	ecs:ECs3718	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4182	cvi:CV2422	41.67	72	42	0	39	110	29	100	2e-06	52.8""".split(
            "\n"
        )

        self.rec3 = """# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ece:Z4181	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4181	spt:SPA2730	37.50	72	45	0	39	110	30	101	1e-05	49.8
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 2
# Query: ece:Z4181
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4181	ecs:ECs3717	100.00	110	0	0	1	110	1	110	3e-54	 211
ece:Z4181	cvi:CV2421	41.67	72	42	0	39	110	29	100	2e-08	59.0
# BLASTP 2.2.10 [Oct-19-2004]
# Iteration: 1
# Query: ece:Z4182
# Database: db/everything.faa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ece:Z4182	ece:Z4182	100.00	110	0	0	1	110	1	110	3e-47	 187
ece:Z4182	cvi:CV2422	41.67	72	42	0	39	110	29	100	2e-06	52.8""".split(
            "\n"
        )

    def test_iter_finder(self):
        """iter_finder should split on lines starting with '# Iteration:'"""
        lines = "abc\n# Iteration: 3\ndef".splitlines()
        self.assertEqual(list(map(iter_finder, lines)), [False, True, False])

    def test_query_finder(self):
        """query_finder should split on lines starting with '# Query:'"""
        lines = "abc\n# Query: dfdsffsd\ndef".split("\n")
        self.assertEqual(list(map(query_finder, lines)), [False, True, False])

    def test_iteration_set_finder(self):
        """iter_finder should split on lines starting with '# Iteration:'"""
        lines = "abc\n# Iteration: 3\ndef\n# Iteration: 1".split("\n")
        self.assertEqual(
            list(map(iteration_set_finder, lines)), [False, False, False, True]
        )

    def test_is_junk(self):
        """is_junk should reject an assortment of invalid lines"""
        # Note: testing two functions that call it instead of function itself
        lines = "abc\n# BLAST blah blah\n   \n# BLAT blah\n123".split("\n")
        self.assertEqual(
            list(map(is_blast_junk, lines)), [False, True, True, False, False]
        )
        self.assertEqual(
            list(map(is_blat_junk, lines)), [False, False, True, True, False]
        )

    def test_make_label(self):
        """make_label should turn comment lines into (key, val) pairs"""
        a = "this test will fail: no # at start"
        b = "#this test will fail because no colon"
        c = "# Iteration: 1"
        d = "# Query: ece:Z4147  ygdP; putative invasion protein [EC:3.6.1.-]"
        e = "#Iteration: 1"  # no space after the hash
        self.assertRaises(ValueError, make_label, a)
        self.assertRaises(ValueError, make_label, b)
        # Note that we _do_ map the data type of known values value, so the
        # value of the iteration will be 1, not '1'
        self.assertEqual(make_label(c), ("ITERATION", 1))
        self.assertEqual(
            make_label(d),
            ("QUERY", "ece:Z4147  ygdP; putative invasion protein [EC:3.6.1.-]"),
        )
        self.assertEqual(make_label(e), ("ITERATION", 1))

    def test_TableToValues(self):
        """TableToValues should convert itself into the correct type."""
        constructors = {"a": int, "b": float, "c": str}
        table = [["c", "b", "a", "d"], ["1.5", "3.5", "2", "2.5"], ["1", "2", "3", "4"]]
        self.assertEqual(
            TableToValues(table, constructors),
            ([["1.5", 3.5, 2, "2.5"], ["1", 2.0, 3, "4"]], ["c", "b", "a", "d"]),
        )
        # check that it works with supplied header
        self.assertEqual(
            TableToValues(table[1:], constructors, list("cbad")),
            ([["1.5", 3.5, 2, "2.5"], ["1", 2.0, 3, "4"]], ["c", "b", "a", "d"]),
        )

    def test_PsiBlastTableParser(self):
        """PsiBlastTableParser should wrap values in table."""
        fields = [
            v.strip()
            for v in "Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score".split(
                ","
            )
        ]

        table = [
            v.split()
            for v in """ece:Z4147       ece:Z4147       100.00  176     0       0       1       176     1       176     2e-89    328
        ece:Z4147       ecs:ECs3687     100.00  176     0       0       1       176     1       176     2e-89    328
        ece:Z4147       ecc:c3425       100.00  176     0       0       1       176     1       176     2e-89    328
        ece:Z4147       sfl:SF2840      100.00  176     0       0       1       176     1       176     2e-89    328""".splitlines()
        ]
        headed_table = [fields] + table
        new_table, new_fields = PsiBlastTableParser(headed_table)
        self.assertEqual(new_fields, fields)
        self.assertEqual(len(new_table), 4)
        self.assertEqual(
            new_table[1],
            ["ece:Z4147", "ecs:ECs3687", 100.0, 176, 0, 0, 1, 176, 1, 176, 2e-89, 328],
        )

    def test_GenericBlastParser9(self):
        """GenericBlastParser9 should read blast's tabular format (#9)."""
        rec = self.rec
        p = GenericBlastParser9(rec, PsiBlastFinder)
        result = list(p)
        self.assertEqual(len(result), 2)
        first, second = result
        self.assertEqual(
            first[0],
            {
                "ITERATION": 1,
                "QUERY": "ece:Z4181",
                "DATABASE": "db/everything.faa",
                "FIELDS": "Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score",
            },
        )
        self.assertEqual(len(first[1]), 3)
        self.assertEqual(second[0]["ITERATION"], 2)
        self.assertEqual(len(second[1]), 7)
        self.assertEqual(
            second[1][-1],
            "ece:Z4181   stm:STM2872 37.50   72  45  0   39  110 30  101 1e-05   49.8".split(),
        )

    def test_PsiBlastParser9(self):
        """PsiBlastParser9 should provide convenient results for format #9."""
        result = PsiBlastParser9(self.rec2)
        self.assertEqual(len(result), 2)
        assert "ece:Z4181" in result
        assert "ece:Z4182" in result
        first = result["ece:Z4181"]
        second = result["ece:Z4182"]
        self.assertEqual(len(first), 2)
        self.assertEqual(len(second), 1)
        iter_1 = first[0]
        iter_2 = first[1]
        self.assertEqual(len(iter_1), 3)
        self.assertEqual(len(iter_2), 7)
        iter_1_2 = second[0]
        self.assertEqual(len(iter_1_2), 3)
        self.assertEqual(len(result["ece:Z4181"][1][3]), 12)
        self.assertEqual(result["ece:Z4181"][1][3]["ALIGNMENT LENGTH"], 103)

    def test_LastProteinIds9(self):
        """LastProteinIds9 should give last protein ids in iter"""
        result = LastProteinIds9(self.rec)
        self.assertEqual(
            result,
            [
                "ece:Z4181",
                "ecs:ECs3717",
                "cvi:CV2421",
                "sfl:CP0138",
                "spt:SPA2730",
                "sec:SC2804",
                "stm:STM2872",
            ],
        )
        # should also work if threshold set
        result = LastProteinIds9(self.rec, False, threshold=8e-6)
        self.assertEqual(
            result, ["ece:Z4181", "ecs:ECs3717", "cvi:CV2421", "sfl:CP0138"]
        )
        # should work on multiple records
        result = list(map(LastProteinIds9, PsiBlastQueryFinder(self.rec2)))
        self.assertEqual(len(result), 2)
        self.assertEqual(
            result[0],
            [
                "ece:Z4181",
                "ecs:ECs3717",
                "cvi:CV2421",
                "sfl:CP0138",
                "spt:SPA2730",
                "sec:SC2804",
                "stm:STM2872",
            ],
        )
        self.assertEqual(result[1], ["ece:Z4182", "ecs:ECs3718", "cvi:CV2422"])

    def test_QMEBlast9(self):
        """QMEBlast9 should return expected lines from all iterations"""
        expect = list(
            zip(
                *[
                    ("ece:Z4181", "ece:Z4181", 3e-47),
                    ("ece:Z4181", "ecs:ECs3717", 3e-47),
                    ("ece:Z4181", "spt:SPA2730", 1e-5),
                    ("ece:Z4181", "ecs:ECs3717", 3e-54),  # WARNING: allows duplicates
                    ("ece:Z4181", "cvi:CV2421", 2e-8),
                    ("ece:Z4182", "ece:Z4182", 3e-47),
                    ("ece:Z4182", "cvi:CV2422", 2e-6),
                ],
            )
        )
        got = list(zip(*QMEBlast9(self.rec3)))
        assert_equal(got[:-1], expect[:-1])
        assert_allclose(got[-1], expect[-1])

    def test_QMEPsiBlast9(self):
        """QMEPsiBlast9 should only return items from last iterations"""
        expect = list(
            zip(
                *[
                    ("ece:Z4181", "ecs:ECs3717", 3e-54),
                    ("ece:Z4181", "cvi:CV2421", 2e-8),
                    ("ece:Z4182", "ece:Z4182", 3e-47),
                    ("ece:Z4182", "cvi:CV2422", 2e-6),
                ]
            )
        )
        got = list(zip(*QMEPsiBlast9(self.rec3)))
        assert_equal(got[:-1], expect[:-1])
        assert_allclose(got[-1], expect[-1])

    def test_fastacmd_taxonomy_splitter(self):
        """fastacmd_taxonomy_splitter should split records into groups"""
        text = """NCBI sequence id: gi|3021565|emb|AJ223314.1|PSAJ3314
NCBI taxonomy id: 3349
Common name: Scots pine
Scientific name: Pinus sylvestris

NCBI sequence id: gi|37777029|dbj|AB108787.1|
NCBI taxonomy id: 228610
Common name: cf. Acremonium sp. KR21-2
Scientific name: cf. Acremonium sp. KR21-2

""".splitlines()
        recs = list(fastacmd_taxonomy_splitter(text))
        self.assertEqual(len(recs), 2)
        self.assertEqual(recs[0], text[:5])  # includes trailing blank

    def test_FastaCmdTaxonomyParser(self):
        """FastaCmdTaxonomyParser should parse taxonomy record to dict"""
        text = """NCBI sequence id: gi|3021565|emb|AJ223314.1|PSAJ3314
NCBI taxonomy id: 3349
Common name: Scots pine
Scientific name: Pinus sylvestris

NCBI sequence id: gi|37777029|dbj|AB108787.1|
NCBI taxonomy id: 228610
Common name: cf. Acremonium sp. KR21-2
Scientific name: cf. Acremonium sp. KR21-2

""".splitlines()
        recs = list(FastacmdTaxonomyParser(text))
        self.assertEqual(len(recs), 2)
        for r in recs:
            self.assertEqual(
                sorted(r.keys()), ["common_name", "scientific_name", "seq_id", "tax_id"]
            )
        r0, r1 = recs
        self.assertEqual(r0["tax_id"], "3349")
        self.assertEqual(r0["common_name"], "Scots pine")
        self.assertEqual(r0["scientific_name"], "Pinus sylvestris")
        self.assertEqual(r0["seq_id"], "gi|3021565|emb|AJ223314.1|PSAJ3314")
        self.assertEqual(r1["tax_id"], "228610")


if __name__ == "__main__":
    main()
