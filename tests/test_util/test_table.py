#!/usr/bin/env python

"""Unit tests for table.
"""
import os
import pathlib
import pickle

from tempfile import TemporaryDirectory
from unittest import TestCase, main, skipIf

import numpy

from numpy.testing import assert_equal

from cogent3 import load_table, make_table
from cogent3.util.table import (
    Table,
    cast_str_to_array,
    cast_to_array,
    formatted_array,
)


try:
    from pandas import DataFrame
except ImportError:
    DataFrame = None


__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La", "Christopher Bradley"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class TrapOutput:
    def __call__(self, data, *args, **kwargs):
        self.data, _ = data._get_repr_()
        self.output = repr(data)


class TableTests(TestCase):
    """Tests of individual functions in table"""

    # test table 1
    t1_header = ["chrom", "stableid", "length"]
    t1_rows = [
        ["X", "ENSG00000005893", 1353],
        ["A", "ENSG00000019485", 1827],
        ["A", "ENSG00000019102", 999],
        ["X", "ENSG00000012174", 1599],
        ["X", "ENSG00000010671", 1977],
        ["A", "ENSG00000019186", 1554],
        ["A", "ENSG00000019144", 4185],
        ["X", "ENSG00000008056", 2307],
        ["A", "ENSG00000018408", 1383],
        ["A", "ENSG00000019169", 1698],
    ]

    # test table 2
    t2_header = ["id", "foo", "bar"]
    t2_rows = [
        [1, "abc", 11],
        [2, "bca", 22],
        [3, "cab", 33],
        [4, "abc", 44],
        [5, "bca", 55],
    ]

    # test table 3
    t3_header = ["id", "foo", "bar"]
    t3_rows = [[6, "abc", 66], [7, "bca", 77]]

    # test table 4
    t4_header = ["id", "foo", "bar"]
    t4_rows = [[8, "abc", 88], [9, "bca", 99]]

    # test table 5
    t5_header = ["a", "b", "c", "d"]
    t5_rows = [[1, 1, 1, 1], [2, 0, 1, 1], [1, 3, 2, 2]]

    # test table 6
    t6_header = ["id", "foo", "bar"]
    t6_rows = [["60", " | ", "666"], ["70", "bca", "777"]]

    # test table 7
    t7_header = ["chrom", "gene", "stat"]
    t7_rows = [
        ["X", "ENSG00000005893", 0.1353],
        ["A", "ENSG00000019485", 1827.558],
        ["A", "ENSG00000019102", 999.125],
    ]

    t8_header = ["edge.name", "edge.parent", "length", "x", "y", "z"]
    t8_rows = [
        ["Human", "edge.0", 4.0, 1.0, 3.0, 6.0],
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
    ]

    def test_empty(self):
        """create table with no data"""
        for data in [None, [], {}, ()]:
            t = Table(header=["col 1", "col 2"], data=data)
            self.assertEqual(len(t), 0)
            self.assertEqual(t.shape, (0, 2), f"failed with {data}")

    def test_keys_are_str(self):
        """all column headers converted to str"""
        t = Table(header=["col 1", 2], data=[[0, 1]])
        self.assertEqual(t.header, ("col 1", "2"))

    def test_no_index_name(self):
        """assigning None has no effect"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        t.index_name = None
        self.assertIs(t.index_name, None)

    def test_index_name(self):
        """correctly assigns"""
        t = Table(header=self.t3_header, data=self.t3_rows, row_ids="foo")
        self.assertEqual(t.index_name, "foo")
        # fails if not an existing column
        with self.assertRaises(ValueError):
            t.index_name = "missing"

        data = t.columns.to_dict()
        # correctly handled when provided on construction
        with self.assertRaises(ValueError):
            t = Table(data=data, row_ids="missing")

        t = Table(data=data, row_ids="foo")
        self.assertEqual(t.index_name, "foo")

        # ... prior to providing columns
        t = Table(row_ids="foo")
        for c, v in data.items():
            t.columns[c] = v
        self.assertEqual(t.index_name, "foo")
        t = Table(row_ids="missing")
        for c, v in data.items():
            t.columns[c] = v

        with self.assertRaises(ValueError):
            t.index_name

    def test_column_repr_str(self):
        """repr and str of Columns"""
        t = Table(header=list("abcdefg"), data=[[0, 1.1, 2, 3, 4, 5.5, 6.0]])
        r = repr(t.columns)
        s = str(t.columns)
        self.assertIsInstance(r, str)
        self.assertEqual(r, s)

    def test_slicing_columns(self):
        """works using names, ints, bool array"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        n = t.columns["b":"d"]
        assert_equal(n, [[1, 0, 3], [1, 1, 2]])
        n = t.columns[1:3]
        assert_equal(n, [[1, 0, 3], [1, 1, 2]])
        indices = numpy.array([False, True, True, False])
        n = t.columns[indices]
        assert_equal(n, [[1, 0, 3], [1, 1, 2]])

    def test_indexing_columns(self):
        """works using names or ints"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        n = t.columns["b"]
        assert_equal(n, [1, 0, 3])
        n = t.columns[1]
        assert_equal(n, [1, 0, 3])

    def test_indexing_rows(self):
        """works using names or ints"""
        t = Table(header=self.t7_header, data=self.t7_rows, row_ids="gene")
        self.assertEqual(t["ENSG00000019485", "chrom"], "A")

    def test_immutability_cells(self):
        """table cells are immutable"""
        t = Table(header=self.t7_header, data=self.t7_rows, row_ids="gene")
        with self.assertRaises(TypeError):
            t["ENSG00000019485", "chrom"] = "D"

        # even via column instance
        with self.assertRaises(ValueError):
            t.columns["chrom"]["ENSG00000019485"] = "D"

    def test_slicing_table(self):
        """works using column names, ints, bool array"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        n = t[:, "b":"d"]
        self.assertEqual(n.columns.order, ("b", "c"))
        assert_equal(n.array, [[1, 1], [0, 1], [3, 2]])
        # using numpy bool array
        rows = numpy.array([True, False, True])
        n = t[rows]
        self.assertEqual(n.shape, (2, 4))
        assert_equal(n.columns["a"], [1, 1])
        rows = numpy.array([True, False, True])
        columns = numpy.array([True, False, False, True])
        n = t[rows, columns]
        assert_equal(n.header, numpy.array(t.header)[columns])
        self.assertEqual(n.shape, (2, 2))

    def test_specifying_space(self):
        """controls spacing in simple format"""
        space = "        "
        t4 = Table(header=self.t1_header, data=self.t1_rows)
        orig = len(str(t4).splitlines()[1])
        t4 = Table(header=self.t1_header, data=self.t1_rows, space=space)
        got1 = len(str(t4).splitlines()[1])
        self.assertTrue(got1 > orig)
        # repr is same
        got2 = len(repr(t4).splitlines()[1])
        self.assertEqual(got1, got2)

    def test_construct_from_dict2d(self):
        """construction from a 2D dict"""
        data = {
            "edge.parent": {
                "NineBande": "root",
                "edge.1": "root",
                "DogFaced": "root",
                "Human": "edge.0",
                "edge.0": "edge.1",
                "Mouse": "edge.1",
                "HowlerMon": "edge.0",
            },
            "x": {
                "NineBande": 1.0,
                "edge.1": 1.0,
                "DogFaced": 1.0,
                "Human": 1.0,
                "edge.0": 1.0,
                "Mouse": 1.0,
                "HowlerMon": 1.0,
            },
        }
        t = Table(data=data)
        self.assertEqual(t.shape, (len(data["x"]), len(data)))

    def test_wrapping_tables_row_ids_1row(self):
        """correctly wraps table to <= maximum width"""
        row_ids = "A/C"
        h = ["A/C", "A/G", "A/T", "C/A"]
        rows = [[0.0425, 0.1424, 0.0226, 0.0391]]
        t = Table(header=h, data=rows, max_width=30, row_ids=row_ids)
        wrapped = str(t)
        # index column occurs twice for these conditions
        for c in h:
            expect = 2 if c == row_ids else 1
            self.assertEqual(wrapped.count(c), expect)

    def test_wrapping_tables_row_ids_multirow(self):
        """correctly wraps table to <= maximum width"""
        # multi-row table
        d2D = {
            "edge.parent": {"NineBande": "root", "Human": "edge.0",},
            "x": {"NineBande": 1.0, "Human": 1.0,},
            "length": {"NineBande": 4.0, "Human": 4.0,},
            "y": {"NineBande": 3.0, "Human": 3.0,},
            "z": {"NineBande": 6.0, "Human": 6.0,},
            "edge.name": {"Human": "Human", "NineBande": "NineBande",},
        }
        row_order = list(d2D["edge.name"])
        t = Table(
            ["edge.name", "edge.parent", "length", "x", "y", "z"],
            d2D,
            row_order=row_order,
            space=8,
            max_width=50,
            row_ids="edge.name",
            title="My title",
            legend="legend: this is a nonsense example.",
        )
        wrapped = str(t)
        for line in wrapped.splitlines():
            len(line)

    def test_wrapping_tables(self):
        """correctly wraps table to <= maximum width"""
        h = ["A/C", "A/G", "A/T", "C/A"]
        rows = [[0.0425, 0.1424, 0.0226, 0.0391]]
        t = Table(header=h, data=rows, max_width=30)
        wrapped = str(t)
        # index column occurs twice for these conditions
        for c in h:
            self.assertEqual(wrapped.count(c), 1)

        # multi-row table
        data = {
            "edge.parent": {"NineBande": "root", "edge.1": "root",},
            "x": {"NineBande": 1.0, "edge.1": 1.0,},
            "length": {"NineBande": 4.0, "edge.1": 4.0,},
            "y": {"NineBande": 3.0, "edge.1": 3.0,},
            "z": {"NineBande": 6.0, "edge.1": 6.0,},
        }
        t = Table(data=data, max_width=30)
        wrapped = str(t)
        for c in data:
            self.assertEqual(wrapped.count(c), 1)

    def test_format_array(self):
        """correctly format array data"""
        f = (2.53, 12.426, 9.9, 7.382e-08)
        # with default format_spec
        g, l = formatted_array(numpy.array(f), "LR", precision=2)
        self.assertTrue(l.endswith("LR"))
        for e in g:
            v = e.split(".")
            self.assertEqual(len(v[-1]), 2, v)
        # handles bool
        g, l = formatted_array(numpy.array([True, False, True]), "LR", precision=2)
        self.assertEqual(g[0].strip(), "True")
        # title is always right aligned
        _, l = formatted_array(numpy.array(f), "LR", format_spec=">.1f")
        self.assertTrue(l.endswith("LR"))
        _, l = formatted_array(numpy.array(f), "LR", format_spec="<.1f")
        self.assertTrue(l.endswith("LR"))

        # using format_spec with right alignment character
        g, l = formatted_array(numpy.array(f), "   blah", format_spec=">.1f")
        for e in g:
            # padded with spaces
            self.assertTrue(e.startswith(" "), e)
            self.assertFalse(e.endswith(" "), e)

        # using format_spec with left alignment character
        g, l = formatted_array(numpy.array(f), "    blah", format_spec="<.1f")
        for e in g:
            # padded with spaces
            self.assertTrue(e.endswith(" "), e)
            self.assertFalse(e.startswith(" "), e)

        # using format_spec with center alignment character
        g, l = formatted_array(numpy.array(f), "    blah", format_spec="^.1f")
        for e in g:
            # padded with spaces
            self.assertTrue(e.endswith(" "), e)
            self.assertTrue(e.startswith(" "), e)

        g, _ = formatted_array(numpy.array(f), "blah", format_spec=".4f")
        for e in g:
            v = e.split(".")
            self.assertEqual(len(v[-1]), 4, v)

        # cope with C-style format strings
        g, _ = formatted_array(numpy.array(f), "blah", format_spec="%.4f")
        for e in g:
            v = e.split(".")
            self.assertEqual(len(v[-1]), 4, v)

        # handle a formatter function
        def formatcol(value):
            if isinstance(value, float):
                val = "%.2f" % value
            else:
                val = str(value)
            return val

        o = [3, "abc", 3.456789]
        g, _ = formatted_array(numpy.array(o, dtype="O"), "blah", format_spec=formatcol)
        self.assertEqual(g[0], "   3", g[0])
        self.assertEqual(g[1], " abc", g[1])
        self.assertEqual(g[2], "3.46", g)

    def test_cast_to_array(self):
        """correctly cast to numpy array"""
        b = (True, False, True)
        a = cast_to_array(b)
        self.assertTrue("bool" in a.dtype.name)
        self.assertEqual(a.tolist(), list(b))
        s = (
            "NP_003077_hs_mm_rn_dna",
            "NP_004893_hs_mm_rn_dna",
            "NP_005079_hs_mm_rn_dna",
            "NP_005500_hs_mm_rn_dna",
            "NP_055852_hs_mm_rn_dna",
        )
        a = cast_to_array(s)
        self.assertTrue("str" in a.dtype.name)
        self.assertEqual(a.tolist(), list(s))

        f = (2.53, 12.426, 9.9, 7.382e-08)
        a = cast_to_array(f)
        self.assertTrue("float" in a.dtype.name, a.dtype.name)
        self.assertEqual(a.tolist(), list(f))

        d = [3, 4, 5]
        a = cast_to_array(d)
        self.assertTrue("int" in a.dtype.name, a.dtype.name)
        self.assertEqual(a.tolist(), list(d))

        o = [3, "abc", 3.4]
        a = cast_to_array(o)
        self.assertTrue("object" in a.dtype.name, a.dtype.name)
        self.assertEqual(a.tolist(), list(o))

    def test_make_table(self):
        """makes a table"""
        data = {
            "edge.parent": {
                "NineBande": "root",
                "edge.1": "root",
                "DogFaced": "root",
                "Human": "edge.0",
            },
            "x": {"NineBande": 1.0, "edge.1": 1.0, "DogFaced": 1.0, "Human": 1.0,},
            "length": {"NineBande": 4.0, "edge.1": 4.0, "DogFaced": 4.0, "Human": 4.0,},
            "y": {"NineBande": 3.0, "edge.1": 3.0, "DogFaced": 3.0, "Human": 3.0,},
            "z": {"NineBande": 6.0, "edge.1": 6.0, "DogFaced": 6.0, "Human": 6.0,},
            "edge.names": {
                "NineBande": "NineBande",
                "edge.1": "edge.1",
                "DogFaced": "DogFaced",
                "Human": "Human",
            },
        }
        t = make_table(data=data)
        self.assertEqual(t.shape, (4, 6))
        # if index column not specified
        with self.assertRaises(IndexError):
            _ = t["Human", "edge.parent"]

        # applies row_ids as an index
        t = make_table(data=data, row_ids="edge.names")
        # index col is the first one, and the data can be indexed
        self.assertEqual(t.columns.order[0], "edge.names")
        self.assertEqual(t["Human", "edge.parent"], "edge.0")

    def test_modify_title_legend(self):
        """reflected in persistent attrs"""
        rows = (
            ("NP_003077_hs_mm_rn_dna", "Con", 2.5386013224378985),
            ("NP_004893_hs_mm_rn_dna", "Con", 0.12135142635634111e06),
        )
        t = Table(["Gene", "Type", "LR"], rows)
        t.title = "a new one"
        self.assertEqual(t._get_persistent_attrs()["title"], "a new one")
        t.legend = "a new 2"
        self.assertEqual(t._get_persistent_attrs()["legend"], "a new 2")

    def test_dunder_repr_eq_str(self):
        """dunder str and repr methods should produce same"""
        rows = (
            ("NP_003077_hs_mm_rn_dna", "Con", 2.5386013224378985),
            ("NP_004893_hs_mm_rn_dna", "Con", 0.12135142635634111e06),
        )
        t = Table(["Gene", "Type", "LR"], rows)
        t.format_column("LR", ".4e")
        s = str(t)
        r = repr(t)
        self.assertTrue(r.startswith(s))

    @skipIf(DataFrame is None, "pandas not installed")
    def test_make_table_from_dataframe(self):
        """makes a table from a pandas data frame"""
        df = DataFrame(data=[[0, 1], [3, 7]], columns=["a", "b"])
        t = make_table(data_frame=df)
        assert_equal(t.columns["a"], [0, 3])
        assert_equal(t.columns["b"], [1, 7])
        with self.assertRaises(TypeError):
            make_table(data_frame="abcde")

    def test_appended(self):
        """test the table appended method"""
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        t4 = Table(header=self.t4_header, data=self.t4_rows)

        append_0 = t2.appended("foo2", [], title="self")
        self.assertEqual(append_0.shape[0], t2.shape[0])
        # test the title feature
        self.assertEqual(append_0.title, "self")

        append_1 = t2.appended("foo2", [t3])
        self.assertEqual(append_1.shape[0], t2.shape[0] + t3.shape[0])
        # test the new_column feature
        self.assertEqual(append_1.shape[1], 4)
        self.assertEqual(append_1.header[0], "foo2")

        append_2 = t2.appended("foo2", [t3, t4])
        self.assertEqual(append_2.shape[0], t2.shape[0] + t3.shape[0] + t4.shape[0])

    def test_count(self):
        """test the table count method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        self.assertEqual(t1.count('chrom == "X"'), 4)
        self.assertEqual(t1.count('stableid.endswith("6")'), 2)
        self.assertEqual(t1.count("length % 2 == 0"), 2)
        self.assertEqual(t1.count('chrom == "Y"'), 0)
        self.assertEqual(t1.count('length % 2 == 0 and chrom == "A"'), 2)
        self.assertEqual(t1.count('length % 2 == 0 or chrom == "X"'), 6)

        t2 = Table(header=self.t2_header, data=self.t2_rows)
        self.assertEqual(t2.count('foo == "abc"'), 2)
        self.assertEqual(t2.count('foo == "cab"'), 1)
        self.assertEqual(t2.count("bar % 2 == 0"), 2)
        self.assertEqual(t2.count("id == 0"), 0)

    def test_count_unique(self):
        """correctly computes unique values"""
        data = {
            "Project_Code": [
                "Ovary-AdenoCA",
                "Liver-HCC",
                "Panc-AdenoCA",
                "Panc-AdenoCA",
            ],
            "Donor_ID": ["DO46416", "DO45049", "DO51493", "DO32860"],
            "Variant_Classification": ["IGR", "Intron", "Intron", "Intron"],
        }
        table = make_table(data=data)
        co = table.count_unique(["Project_Code", "Variant_Classification"])
        self.assertEqual(co[("Panc-AdenoCA", "Intron")], 2)
        self.assertEqual(co[("Liver-HCC", "IGR")], 0)
        co = table.count_unique("Variant_Classification")
        self.assertEqual(co["Intron"], 3)
        self.assertEqual(co["IGR"], 1)

    def test_distinct_values(self):
        """test the table distinct_values method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        self.assertEqual(len(t1.distinct_values("chrom")), 2)
        self.assertEqual(len(t1.distinct_values("stableid")), 10)
        self.assertEqual(len(t1.distinct_values("length")), 10)

        t2 = Table(header=self.t2_header, data=self.t2_rows)
        self.assertEqual(len(t2.distinct_values("id")), 5)
        self.assertEqual(len(t2.distinct_values("foo")), 3)
        self.assertEqual(len(t2.distinct_values("bar")), 5)
        d = t2.distinct_values("foo")
        self.assertEqual(d, {"cab", "bca", "abc"})

    def test_filtered(self):
        """test the table filtered method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        self.assertEqual(t1.filtered('chrom == "X"').shape[0], 4)
        self.assertEqual(t1.filtered('stableid.endswith("6")').shape[0], 2)
        self.assertEqual(t1.filtered("length % 2 == 0").shape[0], 2)
        self.assertEqual(t1.filtered('chrom == "Y"').shape[0], 0)
        self.assertEqual(t1.filtered('length % 2 == 0 and chrom == "A"').shape[0], 2)
        self.assertEqual(t1.filtered('length % 2 == 0 or chrom == "X"').shape[0], 6)

        t2 = Table(header=self.t2_header, data=self.t2_rows)
        self.assertEqual(t2.filtered('foo == "abc"').shape[0], 2)
        self.assertEqual(t2.filtered('foo == "cab"').shape[0], 1)
        self.assertEqual(t2.filtered("bar % 2 == 0").shape[0], 2)
        self.assertEqual(t2.filtered("id == 0").shape[0], 0)

    def test_filtered_by_column(self):
        """test the table filtered_by_column method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        t2 = Table(header=self.t2_header, data=self.t2_rows)

        def is_numeric(values):
            try:
                sum(values)
            except TypeError:
                return False
            return True

        self.assertEqual(t1.filtered_by_column(is_numeric).shape[1], 1)
        self.assertEqual(t2.filtered_by_column(is_numeric).shape[1], 2)

    def test_get_columns(self):
        """test the table get_columns method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        self.assertEqual(t1.get_columns("chrom").shape[0], t1.shape[0])
        self.assertEqual(t1.get_columns("chrom").shape[1], 1)

        self.assertEqual(t1.get_columns(["chrom", "length"]).shape[0], t1.shape[0])
        self.assertEqual(t1.get_columns(["chrom", "length"]).shape[1], 2)
        # if name_index, includes that in return
        t1 = Table(header=self.t1_header, data=self.t1_rows, row_ids="stableid")
        r = t1.get_columns(["length"])
        self.assertEqual(r.header, ("stableid", "length"))

    def test_joined(self):
        """test the table joined method"""
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        # inner join with defaults
        got = t2.joined(t3)
        self.assertEqual(got.shape[0], 0)

        # inner join test
        self.assertEqual(
            t2.joined(t3, columns_self="foo", columns_other="foo").shape[0], 4
        )
        # merged 'foo' column, so (6-1) columns in join
        self.assertEqual(
            t2.joined(t3, columns_self="foo", columns_other="foo").shape[1], 5
        )
        # non-inner join test (cartesian product of rows)
        got = t2.joined(t3, inner_join=False)
        self.assertEqual(got.shape[0], t2.shape[0] * t3.shape[0])
        self.assertEqual(
            t2.joined(t3, inner_join=False).shape[1], t2.shape[1] + t3.shape[1]
        )

    def test_joined_diff_indexing(self):
        """join handles different indexing"""
        a = Table(
            header=["index", "col2", "col3"],
            data=[[1, 2, 3], [2, 3, 1], [2, 6, 5]],
            title="A",
        )
        b = Table(
            header=["index", "col2", "col3"],
            data=[[1, 2, 3], [2, 2, 1], [3, 6, 3]],
            title="B",
        )
        # index by int
        j1 = a.joined(b, [0, 2])
        # index by column names
        j2 = a.joined(b, ["index", "col3"])
        self.assertEqual(j1.header, j2.header)
        self.assertEqual(str(j1), str(j2))

    def test_normalized(self):
        """test the table normalized method"""
        t5 = Table(header=self.t5_header, data=self.t5_rows)
        self.assertEqual(
            t5.normalized().tolist(t5.header),
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.5, 0.0, 0.25, 0.25],
                [0.125, 0.375, 0.25, 0.25],
            ],
        )
        self.assertEqual(
            t5.normalized(by_row=False).tolist(t5.header),
            [[0.25, 0.25, 0.25, 0.25], [0.5, 0.0, 0.25, 0.25], [0.25, 0.75, 0.5, 0.5]],
        )

    def test_sorted(self):
        """test the table sorted method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        got = t1.sorted("length")
        self.assertEqual(
            got.tolist("length"),
            [999, 1353, 1383, 1554, 1599, 1698, 1827, 1977, 2307, 4185],
        )

        t5 = Table(header=self.t5_header, data=self.t5_rows)
        self.assertEqual(t5.sorted("b").tolist("b"), [0, 1, 3])
        self.assertEqual(t5.sorted().tolist("a"), [1, 1, 2])
        self.assertEqual(t5.sorted(reverse="a").tolist("a"), [2, 1, 1])

        table = Table(
            data={
                "chrom": ("X", "A", "A", "X", "X", "A", "A", "X", "A", "A"),
                "stableid": (
                    "ENSG00000005893",
                    "ENSG00000019485",
                    "ENSG00000019102",
                    "ENSG00000012174",
                    "ENSG00000010671",
                    "ENSG00000019186",
                    "ENSG00000019144",
                    "ENSG00000008056",
                    "ENSG00000018408",
                    "ENSG00000019169",
                ),
                "length": [1353, 1827, 999, 1599, 1977, 1554, 4185, 2307, 1383, 1698],
            }
        )

        table = table.sorted(columns=["chrom", "stableid"])
        last_index = len(table) - 1
        self.assertEqual(table[0, "stableid"], "ENSG00000018408")
        self.assertEqual(table[last_index, "stableid"], "ENSG00000012174")

        table = table.sorted(reverse="stableid")
        self.assertEqual(table[0, "stableid"], "ENSG00000019485")
        self.assertEqual(table[last_index, "stableid"], "ENSG00000005893")

        table = table.sorted(reverse="chrom", columns="length")
        self.assertEqual(table[0, "stableid"], "ENSG00000019102")
        self.assertEqual(table[last_index, "stableid"], "ENSG00000019144")

    def test_summed(self):
        """test the table summed method"""
        t5 = Table(header=self.t5_header, data=self.t5_rows)
        self.assertEqual(t5.summed(), [4, 4, 4, 4])
        self.assertEqual(t5.summed(col_sum=False), [4, 4, 8])
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        self.assertEqual(t2.summed(indices=2), 165)

        mix = Table(header=["A", "B"], data=[[0, ""], [1, 2], [3, 4]])

        self.assertEqual(mix.summed("B", strict=False), 6)
        self.assertEqual(mix.summed(0, col_sum=False, strict=False), 0)
        self.assertEqual(mix.summed(1, col_sum=False), 3)
        self.assertEqual(mix.summed(strict=False), [4, 6])
        self.assertEqual(mix.summed(col_sum=False, strict=False), [0, 3, 7])
        with self.assertRaises(TypeError):
            _ = mix.summed(strict=True)

    def test_tolist(self):
        """test the table tolist method"""
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        self.assertEqual(t3.tolist("id"), [6, 7])
        self.assertEqual(t3.tolist("foo"), ["abc", "bca"])

    def test_tolist_column_order(self):
        """column order of input reflected in result"""
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        rev_order = ["id", "foo", "bar"]
        rev_order.reverse()
        result = t3.tolist(rev_order)
        self.assertEqual(result[0], list(reversed(self.t3_rows[0][:])))

    def test_to_dict(self):
        """cast to 2D dict"""
        t = Table(header=self.t7_header, data=self.t7_rows, digits=1)

    def test_transposed(self):
        """test the table transposed method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        # note the transpose has one less row
        got = t1.transposed("", select_as_header="stableid")
        self.assertEqual(got.shape[0], t1.shape[1] - 1)
        # note the transpose has an extra column
        self.assertEqual(got.shape[1], t1.shape[0] + 1)
        # specifying a column without unique values not supported
        with self.assertRaises(ValueError):
            _ = t1.transposed("", select_as_header="chrom")

    def test_transposed_numeric(self):
        """table transposed with numeric header converted to str's"""
        t = Table(header=self.t2_header, data=self.t2_rows)
        # note the transpose has one less row
        got = t.transposed("", select_as_header="bar")
        self.assertEqual(got.shape[0], t.shape[1] - 1)
        # note the transpose has an extra column
        self.assertEqual(got.shape[1], t.shape[0] + 1)
        self.assertEqual(got.header, ("", "11", "22", "33", "44", "55"))
        r = str(got)  # this should not fail!

    def test_del_column(self):
        """correctly removes the column"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        columns = list(t.columns)
        expect = tuple(columns[1:])
        del t.columns[columns[0]]
        self.assertEqual(t.columns.order, expect)

    def test_take_columns(self):
        """correctly takes columns"""
        t = Table(header=self.t4_header, data=self.t4_rows)
        columns = list(t.columns)
        expect = tuple(columns[1:])
        n = t.columns.take_columns(expect)
        self.assertEqual(n.order, expect)
        n = t.columns.take_columns(columns[0])
        self.assertEqual(n.order, (columns[0],))
        n = t.columns.take_columns(1)
        self.assertEqual(n.order, (columns[1],))

    def test_with_new_column(self):
        """test the table with_new_column method"""
        t5 = Table(header=self.t5_header, data=self.t5_rows)
        t5_row_sum = t5.with_new_column("sum", sum, t5.header)
        self.assertEqual(t5_row_sum.get_columns("sum").tolist(), [4, 4, 8])
        # now using a string expression
        t8 = Table(header=self.t8_header, data=self.t8_rows, row_ids="edge.name")
        n = t8.with_new_column("YZ", callback="y+z")
        assert_equal(n.columns["YZ"], [9.0, 9.0])
        # if the new column alreayb exists, the new table has the newest column
        n2 = t8.with_new_column("YZ", callback="y*z")
        assert_equal(n2.columns["YZ"], [18.0, 18.0])
        self.assertNotEqual(id(n), id(n2))
        # bu the column arrays that have not changed should be equal
        for c in n.columns:
            if c == "YZ":
                self.assertNotEqual(id(n.columns[c]), id(n2.columns[c]))
            else:
                self.assertEqual(id(n.columns[c]), id(n2.columns[c]))

    def test_with_new_header(self):
        """test the table with_new_header method"""
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        t2 = t2.with_new_header("id", "no")
        self.assertEqual(t2.header[0], "no")
        t2 = t2.with_new_header("foo", "moo")
        self.assertEqual(t2.header[1], "moo")
        t2 = t2.with_new_header("moo", "foo")
        self.assertEqual(t2.header[1], "foo")

    def test_formatted_mutable(self):
        """returns a mutable object"""
        # required by formatting functions
        t = Table(self.t1_header, self.t1_rows)
        fmt = t._formatted()
        fmt[0][0] = "24"
        self.assertEqual(fmt[0][0], "24")

    def test_formatted_precision(self):
        """applies numerical precision"""
        # required by formatting functions
        t = Table(self.t7_header, self.t7_rows, digits=1)
        fmt = t._formatted()
        # last column should have single place after decimal
        for l in fmt[1:]:
            decimal = l[-1].strip().split(".")[-1]
            self.assertEqual(len(decimal), 1, l[-1])

    def test_str_md_format(self):
        """str() produces markdown table"""
        md_table = make_table(
            header=["a", "b"], data=[["val1", "val2"], ["has | symbol", "val4"]]
        )
        md = md_table.to_string(format="md")
        self.assertTrue(r"has \| symbol" in md)

    def test_str_tex_format(self):
        """str() produces latex tabular table"""
        tex_table = make_table(
            header=["a", "b"], data=[["val1", "val2"], ["val3", "val4"]]
        )
        tex = tex_table.to_string(format="tex")
        self.assertFalse("caption" in tex)
        # with a title
        tex_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["val3", "val4"]],
            title="a title",
        )
        tex = tex_table.to_string(format="tex")
        tex = tex.splitlines()
        self.assertEqual(tex[-2], r"\caption{a title}")

        tex = tex_table.to_string(format="tex", label="tab:first")
        tex = tex.splitlines()
        self.assertEqual(tex[-3], r"\caption{a title}")
        self.assertEqual(tex[-2], r"\label{tab:first}")

        # with a legend, no title
        tex_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["val3", "val4"]],
            legend="a legend",
        )
        tex = tex_table.to_string(format="tex")
        tex = tex.splitlines()
        # because it's treated as a title by default
        self.assertEqual(tex[-2], r"\caption{a legend}")
        # unless you say not to
        tex = tex_table.to_string(format="tex", concat_title_legend=False)
        tex = tex.splitlines()
        self.assertEqual(tex[-2], r"\caption*{a legend}")
        tex_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["val3", "val4"]],
            title="a title.",
            legend="a legend",
        )
        tex = tex_table.to_string(format="tex")
        tex = tex.splitlines()
        self.assertEqual(tex[-2], r"\caption{a title. a legend}")
        tex = tex_table.to_string(format="tex", concat_title_legend=False)
        tex = tex.splitlines()
        self.assertEqual(tex[2], r"\caption{a title.}")
        self.assertEqual(tex[-2], r"\caption*{a legend}")
        tex = tex_table.to_string(
            format="tex", concat_title_legend=False, label="table"
        )
        tex = tex.splitlines()
        self.assertEqual(tex[2], r"\caption{a title.}")
        self.assertEqual(tex[3], r"\label{table}")

    def test_phylip(self):
        """generates phylip format"""
        rows = [
            ["a", "", 0.088337278874079342, 0.18848582712597683, 0.44084000179091454],
            ["c", 0.088337278874079342, "", 0.088337278874079342, 0.44083999937417828],
            ["b", 0.18848582712597683, 0.088337278874079342, "", 0.44084000179090932],
            ["e", 0.44084000179091454, 0.44083999937417828, 0.44084000179090932, ""],
        ]
        header = ["seq1/2", "a", "c", "b", "e"]
        dist = Table(header=header, data=rows, row_ids="seq1/2")
        r = dist.to_string(format="phylip")
        r = r.splitlines()
        self.assertEqual(r[0].strip(), "4")
        for line in r[1:]:
            line = line.split()
            self.assertTrue(line[0] in dist.header)
            self.assertTrue(line[-1][-1].isdigit())

        line = r[1].split()
        self.assertEqual(line[1], "0.0000", line)

    def test_pickle_unpickle(self):
        """roundtrip via pickling"""
        data = {
            "edge.parent": {"NineBande": "root", "edge.1": "root",},
            "x": {"NineBande": 1.0, "edge.1": 1.0,},
            "length": {"NineBande": 4.0, "edge.1": 4.0,},
            "y": {"NineBande": 3.0, "edge.1": 3.0,},
            "z": {"NineBande": 6.0, "edge.1": 6.0,},
            "edge.name": {"NineBande": "NineBande", "edge.1": "edge.1",},
        }
        t = Table(
            data=data,
            max_width=50,
            row_ids="edge.name",
            title="My title",
            legend="blah",
        )
        # via string
        s = pickle.dumps(t)
        r = pickle.loads(s)
        self.assertEqual(str(t), str(r))
        # via file
        with TemporaryDirectory(".") as dirname:
            path = pathlib.Path(dirname) / "table.pickle"
            t.write(str(path))
            r = load_table(path)
            self.assertEqual(str(t), str(r))

    def test_load_mixed(self):
        """load data with mixed data type columns"""
        t = Table(
            header=["abcd", "data", "float"],
            data=[[str([1, 2, 3, 4, 5]), "0", 1.1], ["x", 5.0, 2.1], ["y", "", 3.1]],
        )
        with TemporaryDirectory(".") as dirname:
            path = pathlib.Path(dirname) / "table.tsv"
            t.write(str(path))
            r = load_table(path)
            self.assertEqual(str(t), str(r))
            self.assertTrue("float", r.columns["float"].dtype.name)

    def test_load_mixed_static(self):
        """load data, mixed data type columns remain as string"""
        t = make_table(header=["A", "B"], data=[[1, 1], ["a", 2]])
        with TemporaryDirectory(".") as dirname:
            path = pathlib.Path(dirname) / "table.txt"
            t.write(str(path), sep="\t")
            # if static types, then mixed columns become strings
            r = load_table(path, sep="\t", static_column_types=True)
            self.assertTrue("str" in r.columns["A"].dtype.name)

    def test_load_mixed_row_lengths(self):
        """skip_inconsistent skips rows that have different length to header"""
        h = list("ABCDE")
        r = [list("12345"), list("000"), list("12345")]
        text = "\n".join(["\t".join(l) for l in [h] + r])
        with TemporaryDirectory(".") as dirname:
            path = pathlib.Path(dirname) / "table.tsv"
            with open(path, "w") as out:
                out.write(text)
            r = load_table(path, skip_inconsistent=True)
            self.assertEqual(r.shape, (2, 5))
            self.assertEqual(r.header, tuple(h))
            self.assertEqual(r.array.tolist(), [list(range(1, 6))] * 2)
            # loading without skip_inconsistent raise ValueError
            with self.assertRaises(ValueError):
                r = load_table(path, skip_inconsistent=False)

    def test_load_table_returns_static_columns(self):
        """for static data, load_table gives same dtypes for static_columns_type=True/False"""
        t = load_table("data/sample.tsv", sep="\t", static_column_types=False)
        is_false = {t.columns[c].dtype.name for c in t.columns}
        t = load_table("data/sample.tsv", sep="\t", static_column_types=True)
        is_true = {t.columns[c].dtype.name for c in t.columns}
        self.assertEqual(is_true, is_false)

    def test_formats(self):
        """exercising the different supported formats"""
        last = ""
        for format, startwith in (
            ("md", "|"),
            ("rst", "+"),
            ("latex", r"\begin"),
            ("markdown", "|"),
            ("csv", "id,"),
            ("tsv", "id\t"),
            ("simple", "="),
        ):
            t3 = Table(
                header=self.t3_header,
                data=self.t3_rows,
                format=format,
                title="A title",
                legend="A legend",
            )
            got = str(t3).splitlines()
            query_line = 1 if format == "simple" else 0
            got = got[query_line]
            self.assertIsInstance(got, str)
            self.assertNotEqual(got, last)
            self.assertTrue(got.startswith(startwith), f"{format}: {got[:10]}")
            last = got

    def test_grid_table_format(self):
        """test the table grid_table_format method"""
        from cogent3.format.table import grid_table_format

        formatted_grid = grid_table_format(
            self.t6_header, self.t6_rows, title="Test", legend="Units"
        )
        self.assertEqual(len(formatted_grid.split("\n")), len(self.t6_rows) * 2 + 7)

        formatted_grid = grid_table_format(
            self.t6_header,
            self.t6_rows,
            title="Really Long Title",
            legend="Extra Long Legend",
        )
        self.assertEqual(len(formatted_grid.split("\n")), len(self.t6_rows) * 2 + 7 + 2)

    def test_markdown(self):
        """Exercising the table markdown method"""
        from cogent3.format.table import markdown

        markdown_table = markdown(self.t6_header, self.t6_rows, justify="crl")
        markdown_list = markdown_table.split("\n")
        self.assertEqual(markdown_list[2].count(r"|"), 5)
        # the pipe symbol should have been escaped
        self.assertEqual(markdown_list[2].count(r"\|"), 1)

        with self.assertRaises(ValueError):
            _ = markdown(self.t6_header, self.t6_rows, justify="cr1")

    def test_repr_html_(self):
        """should produce html"""
        # without an index
        t = Table(header=self.t8_header, data=self.t8_rows)
        html = t._repr_html_()
        # without an index
        t = Table(header=self.t8_header, data=self.t8_rows, row_ids="edge.name")
        html = t._repr_html_()

    def test_separator_format(self):
        """testing separator_format with title and legend, and contents that match the separator"""
        from cogent3.format.table import separator_format

        with self.assertRaises(RuntimeError):
            _ = separator_format(self.t6_header, self.t6_rows)
        separated_table = separator_format(
            self.t6_header, self.t6_rows, sep=" | ", title="Test", legend="Units"
        )
        self.assertEqual(len(separated_table.split("\n")), len(self.t6_rows) + 3)

    def test_separator_format_writer(self):
        """exercising separator_format_writer"""
        from cogent3.format.table import separator_formatter

        t3 = Table(header=self.t3_header, data=self.t3_rows)
        comma_sep = t3.to_string(sep=",").splitlines()
        writer = separator_formatter(sep=" | ")
        formatted = [
            f for f in writer([l.split(",") for l in comma_sep], has_header=True)
        ]
        expected_format = ["id | foo | bar", " 6 | abc |  66", " 7 | bca |  77"]
        self.assertEqual(formatted, expected_format)

    def test_set_repr_policy(self):
        """exercising setting repr policy"""
        t = Table(header=self.t2_header, data=self.t2_rows)
        t.set_repr_policy(random=2)
        r = repr(t)
        self.assertIsInstance(r, str)
        r, _ = t._get_repr_()
        self.assertEqual(r.shape[0], 2)
        t.set_repr_policy(head=1)
        r, _ = t._get_repr_()
        self.assertEqual(r.shape[0], 1)
        t.set_repr_policy(tail=3)
        r, _ = t._get_repr_()
        self.assertEqual(r.shape[0], 3)

    def test_head(self):
        """returns the head of the table!"""
        from cogent3.util import table

        display = table.display
        head = TrapOutput()
        table.display = head
        t = Table(header=self.t1_header, data=self.t1_rows)
        t.head(nrows=3)
        self.assertEqual(head.data.shape[0], 3)
        self.assertEqual(len(head.output.splitlines()), 9)
        self.assertEqual(head.data.tolist(), self.t1_rows[:3])
        table.display = display

    def test_tail(self):
        """returns the tail of the table!"""
        from cogent3.util import table

        display = table.display
        tail = TrapOutput()
        table.display = tail
        t = Table(header=self.t1_header, data=self.t1_rows)
        t.tail(nrows=3)
        self.assertEqual(tail.data.shape[0], 3)
        self.assertEqual(len(tail.output.splitlines()), 9)
        self.assertEqual(tail.data.tolist(), self.t1_rows[-3:])
        table.display = display

    @skipIf(DataFrame is None, "pandas not installed")
    def test_to_dataframe(self):
        """produces a dataframe"""
        t = Table(header=self.t1_header, data=self.t1_rows)
        df = t.to_dataframe()
        self.assertIsInstance(df, DataFrame)
        data = df.to_numpy()
        self.assertEqual(data.tolist(), self.t1_rows)

    def test_load_table(self):
        """exercising load table"""
        path = os.path.dirname(os.path.dirname(__file__))
        path = os.path.join(path, "data/sample.tsv")
        table = load_table(path)
        self.assertEqual(table.shape, (10, 3))

    def test_cast_str_to_array(self):
        """handle processing string series"""
        d = [".123|.345", "123"]
        r = cast_str_to_array(d, static_type=True)
        self.assertTrue("str" in r.dtype.name)
        r = cast_str_to_array(d, static_type=False)
        self.assertEqual(r.dtype.name, "object")
        d = [".123|.345", "123", "()"]
        r = cast_str_to_array(d, static_type=False)
        self.assertEqual(r[-1], ())


if __name__ == "__main__":
    main()
