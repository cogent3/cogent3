"""Unit tests for table."""

import json
import os
import pathlib
import pickle
from collections import defaultdict
from unittest import TestCase, skipIf

import numpy
import pytest
from numpy import arange
from numpy.exceptions import ComplexWarning
from numpy.testing import assert_equal

from cogent3 import load_table, make_table, open_
from cogent3.core.table import (
    Table,
    cast_str_to_array,
    cast_str_to_numeric,
    cast_to_array,
)
from cogent3.format.table import (
    formatted_array,
    get_continuation_tables_headers,
    is_html_markup,
)
from cogent3.parse.table import FilteringParser
from cogent3.util.misc import get_object_provenance

try:
    from pandas import DataFrame
except ImportError:
    DataFrame = None

TEST_ROOT = pathlib.Path(__file__).parent.parent


class TrapOutput:
    def __call__(self, data, *args, **kwargs):
        self.data, _, _ = data._get_repr_()
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
            assert len(t) == 0
            assert t.shape == (0, 2), f"failed with {data}"

    def test_input_containers(self):
        """should not fail on defaultdict"""
        raw = {"a": [1, 2, 3], "b": [3, 4, 6]}
        data = defaultdict(list)
        data.update(raw)
        t = Table(data=data)
        assert t.shape == (3, 2)

    def test_keys_are_str(self):
        """all column headers converted to str"""
        t = Table(header=["col 1", 2], data=[[0, 1]])
        assert t.header == ("col 1", "2")

    def test_no_index_name(self):
        """assigning None has no effect"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        t.index_name = None
        assert t.index_name is None

    def test_index_name(self):
        """correctly assigns"""
        t = Table(header=self.t3_header, data=self.t3_rows, index_name="foo")
        assert t.index_name == "foo"
        # fails if not an existing column
        with pytest.raises(ValueError):
            t.index_name = "missing"

        data = t.columns.to_dict()
        # correctly handled when provided on construction
        with pytest.raises(ValueError):
            t = Table(data=data, index_name="missing")

        t = Table(data=data, index_name="foo")
        assert t.index_name == "foo"

        # correctly reset when assigned None
        t.index_name = None
        assert t.index_name is None
        assert t.columns.index_name is None
        assert t._template is None

        # ... prior to providing columns
        t = Table(index_name="foo")
        for c, v in data.items():
            t.columns[c] = v
        assert t.index_name == "foo"
        t = Table(index_name="missing")
        for c, v in data.items():
            t.columns[c] = v

        with pytest.raises(ValueError):
            t.index_name

    def test_table_data_int_keys(self):
        """correctly construct table from dict with int's as keys"""
        head = ["", 0, 1]
        data = {0: [2, 2], 1: [2, 2], "": [0, 1]}
        t = Table(head, data=data)
        assert_equal(t.array.tolist(), [[0, 2, 2], [1, 2, 2]])

    def test_table_with_empty_string_index(self):
        """handle an index of empty string"""
        d = {
            "": ["Chimpanzee", "Galago", "Gorilla"],
            "Chimpanzee": [0.0, 0.19, 0.005],
            "Galago": [0.19, 0.0, 0.19],
        }
        table = make_table(data=d, index_name="")
        val = table["Galago", "Chimpanzee"]
        assert val == 0.19

    def test_column_repr_str(self):
        """repr and str of Columns"""
        t = Table(header=list("abcdefg"), data=[[0, 1.1, 2, 3, 4, 5.5, 6.0]])
        r = repr(t.columns)
        s = str(t.columns)
        assert isinstance(r, str)
        assert r == s

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
        t = Table(header=self.t7_header, data=self.t7_rows, index_name="gene")
        got = t["ENSG00000019485", "chrom"]
        assert got == "A"

    def test_immutability_cells(self):
        """table cells are immutable"""
        t = Table(header=self.t7_header, data=self.t7_rows, index_name="gene")
        with pytest.raises(TypeError):
            t["ENSG00000019485", "chrom"] = "D"

        # even via column instance
        with pytest.raises(ValueError):
            t.columns["chrom"]["ENSG00000019485"] = "D"

    def test_slicing_table(self):
        """works using column names, ints, bool array"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        n = t[:, "b":"d"]
        assert n.columns.order == ("b", "c")
        assert_equal(n.array, [[1, 1], [0, 1], [3, 2]])
        # using numpy bool array
        rows = numpy.array([True, False, True])
        n = t[rows]
        assert n.shape == (2, 4)
        assert_equal(n.columns["a"], [1, 1])
        rows = numpy.array([True, False, True])
        columns = numpy.array([True, False, False, True])
        n = t[rows, columns]
        assert_equal(n.header, numpy.array(t.header)[columns])
        assert n.shape == (2, 2)

        # column formatting copied on slice
        t = Table(header=self.t5_header, data=self.t5_rows)
        t.format_column("c", "%.2e")
        n = t[:, 1:]
        assert n._column_templates == t._column_templates

    def test_slicing_using_numpy_indexing(self):
        """support numpy advanced indexing"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        indices = t.columns["b"] != 0
        got = t[indices]
        expect = t.array[[0, 2], :]
        assert_equal(got.array, expect)
        got = t[indices, [True, False, True, True]]
        expect = expect[:, [0, 2, 3]]
        assert_equal(got.array, expect)

        # using numpy arrays for rows and columns
        got_np = t[indices, numpy.array([True, False, True, True])]
        assert_equal(got_np.array, got.array)

    def test_slicing_with_index(self):
        """different slice types work when index_name defined"""
        # slicing by int works with index_name too
        t = Table(header=self.t8_header, data=self.t8_rows, index_name="edge.name")
        got = t[[1]]
        assert got.columns["edge.name"] == "NineBande"
        assert got.shape == (1, t.shape[1])
        for v, dtype in [(1, None), (1, object), ("NineBande", "U")]:
            got = t[numpy.array([v], dtype=dtype)]
            assert got.columns["edge.name"] == "NineBande"
            assert got.shape == (1, t.shape[1])

        # works if, for some reason, the index_name column has floats
        t = Table(header=self.t7_header, data=self.t7_rows, index_name="stat")
        got = t[[1827.5580]]
        assert got.shape == (1, t.shape[1])
        got = t[numpy.array([1827.5580])]
        assert got.shape == (1, t.shape[1])

    def test_specifying_space(self):
        """controls spacing in simple format"""
        space = "        "
        t4 = Table(header=self.t1_header, data=self.t1_rows)
        orig = len(str(t4).splitlines()[1])
        t4 = Table(header=self.t1_header, data=self.t1_rows, space=space)
        got1 = len(str(t4).splitlines()[1])
        assert got1 > orig
        # repr is same
        got2 = len(repr(t4).splitlines()[1])
        assert got1 == got2

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
        assert t.shape == (len(data["x"]), len(data))

    def test_wrapping_tables_index_1row(self):
        """correctly wraps table to <= maximum width"""
        index = "A/C"
        h = ["A/C", "A/G", "A/T", "C/A"]
        rows = [[0.0425, 0.1424, 0.0226, 0.0391]]
        t = Table(header=h, data=rows, max_width=30, index_name=index)
        wrapped = str(t)
        # index_name column occurs twice for these conditions
        for c in h:
            expect = 2 if c == index else 1
            assert wrapped.count(c) == expect

    def test_wrapping_tables_index_multirow(self):
        """correctly wraps table to <= maximum width"""
        # multi-row table
        d2D = {
            "edge.parent": {
                "NineBande": "root",
                "Human": "edge.0",
            },
            "x": {
                "NineBande": 1.0,
                "Human": 1.0,
            },
            "length": {
                "NineBande": 4.0,
                "Human": 4.0,
            },
            "y": {
                "NineBande": 3.0,
                "Human": 3.0,
            },
            "z": {
                "NineBande": 6.0,
                "Human": 6.0,
            },
            "edge.name": {
                "Human": "Human",
                "NineBande": "NineBande",
            },
        }
        row_order = list(d2D["edge.name"])
        t = Table(
            ["edge.name", "edge.parent", "length", "x", "y", "z"],
            d2D,
            row_order=row_order,
            space=8,
            max_width=50,
            index_name="edge.name",
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
        # index_name column occurs twice for these conditions
        for c in h:
            assert wrapped.count(c) == 1

        # multi-row table
        data = {
            "edge.parent": {
                "NineBande": "root",
                "edge.1": "root",
            },
            "x": {
                "NineBande": 1.0,
                "edge.1": 1.0,
            },
            "length": {
                "NineBande": 4.0,
                "edge.1": 4.0,
            },
            "y": {
                "NineBande": 3.0,
                "edge.1": 3.0,
            },
            "z": {
                "NineBande": 6.0,
                "edge.1": 6.0,
            },
        }
        t = Table(data=data, max_width=30)
        wrapped = str(t)
        for c in data:
            assert wrapped.count(c) == 1

    def test_format_array(self):
        """correctly format array data"""
        f = (2.53, 12.426, 9.9, 7.382e-08)
        # with default format_spec
        g, l, w = formatted_array(numpy.array(f), "LR", precision=2)
        assert l.endswith("LR")
        for e in g:
            v = e.split(".")
            assert len(v[-1]) == 2, v
        # handles bool
        g, l, w = formatted_array(numpy.array([True, False, True]), "LR", precision=2)
        assert g[0].strip() == "True"
        # title is always right aligned
        _, l, _ = formatted_array(numpy.array(f), "LR", format_spec=">.1f")
        assert l.endswith("LR")
        _, l, _ = formatted_array(numpy.array(f), "LR", format_spec="<.1f")
        assert l.startswith("LR")

        # using format_spec with right alignment character
        g, l, w = formatted_array(numpy.array(f), "   blah", format_spec=">.1f")
        for e in g:
            # padded with spaces
            assert e.startswith(" "), e
            assert not e.endswith(" "), e

        # using format_spec with left alignment character
        g, l, w = formatted_array(numpy.array(f), "    blah", format_spec="<.1f")
        for e in g:
            # padded with spaces
            assert e.endswith(" "), e
            assert not e.startswith(" "), e

        # using format_spec with center alignment character
        g, l, w = formatted_array(numpy.array(f), "    blah", format_spec="^.1f")
        for e in g:
            # padded with spaces
            assert e.endswith(" "), e
            assert e.startswith(" "), e

        g, _, _ = formatted_array(numpy.array(f), "blah", format_spec=".4f")
        for e in g:
            v = e.split(".")
            assert len(v[-1]) == 4, v

        # cope with C-style format strings
        g, _, _ = formatted_array(numpy.array(f), "blah", format_spec="%.4f")
        for e in g:
            v = e.split(".")
            assert len(v[-1]) == 4, v

        # handle a formatter function
        def formatcol(value):
            return f"{value:.2f}" if isinstance(value, float) else str(value)

        o = [3, "abc", 3.456789]
        g, _, _ = formatted_array(
            numpy.array(o, dtype="O"),
            "blah",
            format_spec=formatcol,
        )
        assert g[0] == "   3", g[0]
        assert g[1] == " abc", g[1]
        assert g[2] == "3.46", g

        # don't pad
        g, l, w = formatted_array(numpy.array(f), "    blah", format_spec="<.1f")
        g, l, w = formatted_array(
            numpy.array(f),
            "    blah",
            format_spec="<.1f",
            pad=False,
        )
        assert l == "blah"
        for v in g:
            assert " " not in v

        # use the align argument, 'c'
        g, l, w = formatted_array(
            numpy.array(f),
            "  blah  ",
            precision=1,
            pad=True,
            align="c",
        )
        for v in g:
            assert v.startswith(" ")
            assert v.endswith(" ")

        # use the align argument, 'l'
        g, l, w = formatted_array(
            numpy.array(f),
            "  blah  ",
            precision=1,
            pad=True,
            align="l",
        )
        for v in g:
            assert not v.startswith(" ")
            assert v.endswith(" ")

        # use the align argument, 'r'
        col_title = "  blah  "
        g, l, w = formatted_array(
            numpy.array(f),
            col_title,
            precision=1,
            pad=True,
            align="r",
        )
        for v in g:
            assert v.startswith(" ")
            assert not v.endswith(" ")

        assert w == len(col_title)

        # raises error if align invalid value
        with pytest.raises(ValueError):
            formatted_array(
                numpy.array(f),
                "  blah  ",
                precision=1,
                pad=True,
                align="blah",
            )

    def test_get_continuation_tables_headers(self):
        """correctly identify the columns for subtables"""
        cols_widths = [("", 10), ("b", 5), ("c", 3), ("d", 14), ("e", 15)]
        got = get_continuation_tables_headers(cols_widths)
        # no subtables, returns list of lists
        expect = [[c for c, _ in cols_widths]]
        assert got == expect
        # fails if any column has a width < max_width
        with pytest.raises(ValueError):
            get_continuation_tables_headers(cols_widths, max_width=5)

        # or if the sum of the index_name width and column is > max_width
        with pytest.raises(ValueError):
            get_continuation_tables_headers(cols_widths, index_name="", max_width=24)

        got = get_continuation_tables_headers(cols_widths, max_width=25)
        expect = [["", "b", "c"], ["d"], ["e"]]
        assert got == expect

        # with an index_name column
        got = get_continuation_tables_headers(cols_widths, index_name="", max_width=27)
        expect = [["", "b", "c"], ["", "d"], ["", "e"]]
        assert got == expect

        cols_widths = [("a", 10), ("b", 5), ("c", 3), ("d", 14), ("e", 15)]
        got = get_continuation_tables_headers(cols_widths, index_name="a", max_width=27)
        expect = [["a", "b", "c"], ["a", "d"], ["a", "e"]]
        assert got == expect

        # space has an affect
        got = get_continuation_tables_headers(cols_widths, max_width=25, space=4)
        expect = [["a", "b"], ["c", "d"], ["e"]]
        assert got == expect

    def test_cast_to_array(self):
        """correctly cast to numpy array"""
        b = (True, False, True)
        a = cast_to_array(b)
        assert "bool" in a.dtype.name
        assert a.tolist() == list(b)
        s = (
            "NP_003077_hs_mm_rn_dna",
            "NP_004893_hs_mm_rn_dna",
            "NP_005079_hs_mm_rn_dna",
            "NP_005500_hs_mm_rn_dna",
            "NP_055852_hs_mm_rn_dna",
        )
        a = cast_to_array(s)
        assert "str" in a.dtype.name
        assert a.tolist() == list(s)

        f = (2.53, 12.426, 9.9, 7.382e-08)
        a = cast_to_array(f)
        assert "float" in a.dtype.name, a.dtype.name
        assert a.tolist() == list(f)

        d = [3, 4, 5]
        a = cast_to_array(d)
        assert "int" in a.dtype.name, a.dtype.name
        assert a.tolist() == list(d)

        o = [3, "abc", 3.4]
        a = cast_to_array(o)
        assert "object" in a.dtype.name, a.dtype.name
        assert a.tolist() == list(o)

    def test_make_table(self):
        """makes a table"""
        data = {
            "edge.parent": {
                "NineBande": "root",
                "edge.1": "root",
                "DogFaced": "root",
                "Human": "edge.0",
            },
            "x": {
                "NineBande": 1.0,
                "edge.1": 1.0,
                "DogFaced": 1.0,
                "Human": 1.0,
            },
            "length": {
                "NineBande": 4.0,
                "edge.1": 4.0,
                "DogFaced": 4.0,
                "Human": 4.0,
            },
            "y": {
                "NineBande": 3.0,
                "edge.1": 3.0,
                "DogFaced": 3.0,
                "Human": 3.0,
            },
            "z": {
                "NineBande": 6.0,
                "edge.1": 6.0,
                "DogFaced": 6.0,
                "Human": 6.0,
            },
            "edge.names": {
                "NineBande": "NineBande",
                "edge.1": "edge.1",
                "DogFaced": "DogFaced",
                "Human": "Human",
            },
        }
        t = make_table(data=data)
        assert t.shape == (4, 6)
        # if index_name column not specified
        with pytest.raises(IndexError):
            _ = t["Human", "edge.parent"]

        # use an index_name
        t = make_table(data=data, index_name="edge.names")
        # index_name col is the first one, and the data can be indexed
        assert t.columns.order[0] == "edge.names"
        assert t["Human", "edge.parent"] == "edge.0"

        # providing path raises TypeError
        with pytest.raises(TypeError):
            make_table("some_path.tsv")

        with pytest.raises(TypeError):
            make_table(header="some_path.tsv")

        with pytest.raises(TypeError):
            make_table(data="some_path.tsv")

    def test_modify_title_legend(self):
        """reflected in persistent attrs"""
        rows = (
            ("NP_003077_hs_mm_rn_dna", "Con", 2.5386013224378985),
            ("NP_004893_hs_mm_rn_dna", "Con", 0.12135142635634111e06),
        )
        t = Table(["Gene", "Type", "LR"], rows)
        t.title = "a new one"
        assert t._get_persistent_attrs()["title"] == "a new one"
        t.legend = "a new 2"
        assert t._get_persistent_attrs()["legend"] == "a new 2"

    def test_dunder_repr_eq_str(self):
        """dunder str and repr methods should produce same"""
        rows = (
            ("NP_003077_hs_mm_rn_dna", "Con", 2.5386013224378985),
            ("NP_004893_hs_mm_rn_dna", "Con", 0.12135142635634111e06),
        )
        t = Table(["Gene", "Type", "LR"], rows)
        t.format_column("LR", "%.4e")
        s = str(t)
        r = repr(t)
        assert r.startswith(s)

    @skipIf(DataFrame is None, "pandas not installed")
    def test_make_table_from_dataframe(self):
        """makes a table from a pandas data frame"""
        df = DataFrame(data=[[0, 1], [3, 7]], columns=["a", "b"])
        t = make_table(data_frame=df)
        assert_equal(t.columns["a"], [0, 3])
        assert_equal(t.columns["b"], [1, 7])
        with pytest.raises(TypeError):
            make_table(data_frame="abcde")

    def test_appended(self):
        """test the table appended method"""
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        t4 = Table(header=self.t4_header, data=self.t4_rows)

        append_0 = t2.appended("foo2", [], title="self")
        assert append_0.shape[0] == t2.shape[0]
        # test the title feature
        assert append_0.title == "self"

        append_1 = t2.appended("foo2", [t3])
        assert append_1.shape[0] == t2.shape[0] + t3.shape[0]
        # test the new_column feature
        assert append_1.shape[1] == 4
        assert append_1.header[0] == "foo2"

        append_2 = t2.appended("foo2", [t3, t4])
        assert append_2.shape[0] == t2.shape[0] + t3.shape[0] + t4.shape[0]

        append_3 = t2.appended("", [t3, t4])
        assert append_3.shape[0] == t2.shape[0] + t3.shape[0] + t4.shape[0]
        assert append_3.shape[1] == t2.shape[1] + 1

    def test_appended_mixed_dtypes(self):
        """handles table columns with different dtypes"""
        t1 = Table(header=["a", "b"], data={"a": [1], "b": ["s"]})
        t2 = Table(header=["a", "b"], data={"a": [1.2], "b": [4]})
        appended = t1.appended(None, t2)
        assert "float" in appended.columns["a"].dtype.name
        assert "object" in appended.columns["b"].dtype.name

    def test_count(self):
        """test the table count method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        assert t1.count('chrom == "X"') == 4
        assert t1.count('stableid.endswith("6")') == 2
        assert t1.count("length % 2 == 0") == 2
        assert t1.count('chrom == "Y"') == 0
        assert t1.count('length % 2 == 0 and chrom == "A"') == 2
        assert t1.count('length % 2 == 0 or chrom == "X"') == 6

        t2 = Table(header=self.t2_header, data=self.t2_rows)
        assert t2.count('foo == "abc"') == 2
        assert t2.count('foo == "cab"') == 1
        assert t2.count("bar % 2 == 0") == 2
        assert t2.count("id == 0") == 0

    def test_count_empty(self):
        """empty table count method returns 0"""
        t1 = Table(header=self.t1_header)
        assert t1.count('chrom == "X"') == 0
        assert t1.count(lambda x: x == "X", columns="chrom") == 0

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
        assert co["Panc-AdenoCA", "Intron"] == 2
        assert co["Liver-HCC", "IGR"] == 0
        co = table.count_unique("Variant_Classification")
        assert co["Intron"] == 3
        assert co["IGR"] == 1

    def test_distinct_values(self):
        """test the table distinct_values method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        assert len(t1.distinct_values("chrom")) == 2
        assert len(t1.distinct_values("stableid")) == 10
        assert len(t1.distinct_values("length")) == 10

        t2 = Table(header=self.t2_header, data=self.t2_rows)
        assert len(t2.distinct_values("id")) == 5
        assert len(t2.distinct_values("foo")) == 3
        assert len(t2.distinct_values("bar")) == 5
        d = t2.distinct_values("foo")
        assert d == {"cab", "bca", "abc"}

    def test_filtered(self):
        """test the table filtered method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        assert t1.filtered('chrom == "X"').shape[0] == 4
        assert t1.filtered('stableid.endswith("6")').shape[0] == 2
        assert t1.filtered("length % 2 == 0").shape[0] == 2
        assert t1.filtered('chrom == "Y"').shape[0] == 0
        assert t1.filtered('length % 2 == 0 and chrom == "A"').shape[0] == 2
        assert t1.filtered('length % 2 == 0 or chrom == "X"').shape[0] == 6

        t2 = Table(header=self.t2_header, data=self.t2_rows)
        assert t2.filtered('foo == "abc"').shape[0] == 2
        assert t2.filtered('foo == "cab"').shape[0] == 1
        assert t2.filtered("bar % 2 == 0").shape[0] == 2
        assert t2.filtered("id == 0").shape[0] == 0

    def test_filtered_empty(self):
        """test the table filtered method"""
        t1 = Table(header=self.t1_header)
        assert t1.shape[0] == 0
        got = t1.filtered('chrom == "X"')
        assert got.shape[0] == 0
        got = t1.filtered(lambda x: x == "X", columns="chrom")
        assert got.shape[0] == 0

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

        assert t1.filtered_by_column(is_numeric).shape[1] == 1
        assert t2.filtered_by_column(is_numeric).shape[1] == 2

    def test_get_columns(self):
        """test the table get_columns method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        assert t1.get_columns("chrom").shape[0] == t1.shape[0]
        assert t1.get_columns("chrom").shape[1] == 1

        assert t1.get_columns(["chrom", "length"]).shape[0] == t1.shape[0]
        assert t1.get_columns(["chrom", "length"]).shape[1] == 2
        # if index_name, includes that in return
        t1 = Table(header=self.t1_header, data=self.t1_rows, index_name="stableid")
        r = t1.get_columns(["length"])
        assert r.header == ("stableid", "length")
        # if index_name, unless excluded
        r = t1.get_columns(["length"], with_index=False)
        assert r.index_name is None

    def test_joined(self):
        """test the table joined method"""
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        # inner join with defaults
        got = t2.joined(t3)
        assert got.shape[0] == 0

        # inner join test
        assert t2.joined(t3, columns_self="foo", columns_other="foo").shape[0] == 4
        # merged 'foo' column, so (6-1) columns in join
        assert t2.joined(t3, columns_self="foo", columns_other="foo").shape[1] == 5
        # non-inner join test (cartesian product of rows)
        got = t2.joined(t3, inner_join=False)
        assert got.shape[0] == t2.shape[0] * t3.shape[0]
        assert t2.joined(t3, inner_join=False).shape[1] == t2.shape[1] + t3.shape[1]

        got = t2.cross_join(t3)
        assert got.shape[0] == t2.shape[0] * t3.shape[0]
        assert t2.joined(t3, inner_join=False).shape[1] == t2.shape[1] + t3.shape[1]

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
        assert j1.header == j2.header
        assert str(j1) == str(j2)

    def test_normalized(self):
        """test the table normalized method"""
        t5 = Table(header=self.t5_header, data=self.t5_rows)
        assert t5.normalized().to_list(t5.header) == [
            [0.25, 0.25, 0.25, 0.25],
            [0.5, 0.0, 0.25, 0.25],
            [0.125, 0.375, 0.25, 0.25],
        ]
        assert t5.normalized(by_row=False).to_list(t5.header) == [
            [0.25, 0.25, 0.25, 0.25],
            [0.5, 0.0, 0.25, 0.25],
            [0.25, 0.75, 0.5, 0.5],
        ]

    def test_sorted(self):
        """test the table sorted method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        got = t1.sorted("length")
        assert got.to_list("length") == [
            999,
            1353,
            1383,
            1554,
            1599,
            1698,
            1827,
            1977,
            2307,
            4185,
        ]

        t5 = Table(header=self.t5_header, data=self.t5_rows)
        assert t5.sorted("b").to_list("b") == [0, 1, 3]
        assert t5.sorted().to_list("a") == [1, 1, 2]
        assert t5.sorted(reverse="a").to_list("a") == [2, 1, 1]

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
            },
        )

        table = table.sorted(columns=["chrom", "stableid"])
        last_index = len(table) - 1
        assert table[0, "stableid"] == "ENSG00000018408"
        assert table[last_index, "stableid"] == "ENSG00000012174"

        table = table.sorted(reverse="stableid")
        assert table[0, "stableid"] == "ENSG00000019485"
        assert table[last_index, "stableid"] == "ENSG00000005893"

        table = table.sorted(reverse="chrom", columns="length")
        assert table[0, "stableid"] == "ENSG00000019102"
        assert table[last_index, "stableid"] == "ENSG00000019144"

        # providing reversed argument name raises TypeError
        with pytest.raises(TypeError):
            table.sorted(reversed="chrom")

    def test_summed(self):
        """test the table summed method"""
        t5 = Table(header=self.t5_header, data=self.t5_rows)
        assert t5.summed() == [4, 4, 4, 4]
        assert t5.summed(col_sum=False) == [4, 4, 8]
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        assert t2.summed(indices=2) == 165

        mix = Table(header=["A", "B"], data=[[0, ""], [1, 2], [3, 4]])

        assert mix.summed("B", strict=False) == 6
        assert mix.summed(0, col_sum=False, strict=False) == 0
        assert mix.summed(1, col_sum=False) == 3
        assert mix.summed(strict=False) == [4, 6]
        assert mix.summed(col_sum=False, strict=False) == [0, 3, 7]
        with pytest.raises(TypeError):
            _ = mix.summed(strict=True)

    def test_to_list(self):
        """test the table to_list method"""
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        assert t3.to_list("id") == [6, 7]
        assert t3.to_list("foo") == ["abc", "bca"]

    def test_to_list_column_order(self):
        """column order of input reflected in result"""
        t3 = Table(header=self.t3_header, data=self.t3_rows)
        rev_order = ["id", "foo", "bar"]
        rev_order.reverse()
        result = t3.to_list(rev_order)
        assert result[0] == list(reversed(self.t3_rows[0][:]))

    def test_to_dict(self):
        """cast to 2D dict"""
        Table(header=self.t7_header, data=self.t7_rows, digits=1)

    def test_transposed(self):
        """test the table transposed method"""
        t1 = Table(header=self.t1_header, data=self.t1_rows)
        # note the transpose has one less row
        got = t1.transposed("", select_as_header="stableid")
        assert got.shape[0] == t1.shape[1] - 1
        # note the transpose has an extra column
        assert got.shape[1] == t1.shape[0] + 1
        # specifying a column without unique values not supported
        with pytest.raises(ValueError):
            _ = t1.transposed("", select_as_header="chrom")

    def test_transposed_numeric(self):
        """table transposed with numeric header converted to str's"""
        t = Table(header=self.t2_header, data=self.t2_rows)
        # note the transpose has one less row
        got = t.transposed("", select_as_header="bar")
        assert got.shape[0] == t.shape[1] - 1
        # note the transpose has an extra column
        assert got.shape[1] == t.shape[0] + 1
        assert got.header == ("", "11", "22", "33", "44", "55")
        str(got)  # this should not fail!

    def test_transposed_forgets_index(self):
        """transposed table defaults to no row index_name"""
        data = {
            "": [0, 1, 2, 3, 4, 5, 6],
            "T": [2, 10, 1, 6, 1, 5, 0],
            "C": [0, 0, 0, 0, 0, 0, 1],
            "A": [8, 0, 9, 4, 9, 4, 4],
            "G": [0, 0, 0, 0, 0, 1, 5],
        }
        t = Table(header=["", "T", "C", "A", "G"], data=data, index_name="")
        tr = t.transposed("Base", select_as_header="")
        assert tr.index_name is None

        # but you can set a new one
        tr = t.transposed("Base", select_as_header="", index_name="Base")
        assert tr.index_name == "Base"
        assert tr["G", "5"] == 1

    def test_del_column(self):
        """correctly removes the column"""
        t = Table(header=self.t5_header, data=self.t5_rows)
        columns = list(t.columns)
        expect = tuple(columns[1:])
        del t.columns[columns[0]]
        assert t.columns.order == expect

    def test_take_columns(self):
        """correctly takes columns"""
        t = Table(header=self.t4_header, data=self.t4_rows)
        columns = list(t.columns)
        expect = tuple(columns[1:])
        n = t.columns.take_columns(expect)
        assert n.order == expect
        n = t.columns.take_columns(columns[0])
        assert n.order == (columns[0],)
        n = t.columns.take_columns(1)
        assert n.order == (columns[1],)

    def test_with_new_column(self):
        """test the table with_new_column method"""
        t5 = Table(header=self.t5_header, data=self.t5_rows)
        t5_row_sum = t5.with_new_column("sum", sum, t5.header)
        assert t5_row_sum.get_columns("sum").to_list() == [4, 4, 8]
        # now using a string expression
        t8 = Table(header=self.t8_header, data=self.t8_rows, index_name="edge.name")
        n = t8.with_new_column("YZ", callback="y+z")
        assert_equal(n.columns["YZ"], [9.0, 9.0])
        # if the new column alreayb exists, the new table has the newest column
        n2 = t8.with_new_column("YZ", callback="y*z")
        assert_equal(n2.columns["YZ"], [18.0, 18.0])
        assert id(n) != id(n2)
        # bu the column arrays that have not changed should be equal
        for c in n.columns:
            if c == "YZ":
                assert id(n.columns[c]) != id(n2.columns[c])
            else:
                assert id(n.columns[c]) == id(n2.columns[c])

    def test_with_new_header(self):
        """test the table with_new_header method"""
        t2 = Table(header=self.t2_header, data=self.t2_rows)
        t2 = t2.with_new_header("id", "no")
        assert t2.header[0] == "no"
        t2 = t2.with_new_header("foo", "moo")
        assert t2.header[1] == "moo"
        t2 = t2.with_new_header("moo", "foo")
        assert t2.header[1] == "foo"

    def test_formatted_mutable(self):
        """returns a mutable object"""
        # required by formatting functions
        t = Table(self.t1_header, self.t1_rows)
        fmt = t._formatted()
        fmt[0][0] = "24"
        assert fmt[0][0] == "24"

    def test_formatted_precision(self):
        """applies numerical precision"""
        # required by formatting functions
        t = Table(self.t7_header, self.t7_rows, digits=1)
        fmt = t._formatted()
        # last column should have single place after decimal
        for l in fmt[1:]:
            decimal = l[-1].strip().split(".")[-1]
            assert len(decimal) == 1, l[-1]

    def test_str_empty(self):
        """empty table returns empty str"""
        table = make_table()
        assert str(table) == ""

    def test_repr_empty(self):
        """empty table returns empty str"""
        table = make_table()
        got = repr(table)
        assert got == "0 rows x 0 columns"

    def test_str_zero_rows(self):
        """table with no rows returns column heads"""
        table = make_table(header=["a"])
        assert str(table) == "=\na\n-\n-"

    def test_str_object_col(self):
        """str works when a column has complex object"""
        # data has tuples in an array
        data = {
            "key": numpy.array([("a", "c"), ("b", "c"), ("a", "d")], dtype="O"),
            "count": [1, 3, 2],
        }
        t = Table(data=data)
        got = str(t)
        assert len(got.splitlines()) == 7

    def test_str_md_format(self):
        """str() produces markdown table"""
        md_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["has | symbol", "val4"]],
        )
        md = md_table.to_string(format="md")
        assert "has \\| symbol" in md

    def test_str_tex_format(self):
        """str() produces latex tabular table"""
        tex_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["val3", "val4"]],
        )
        tex = tex_table.to_string(format="tex", justify="cr")
        assert tex_table.to_string(format="tex", justify="cr") == tex_table.to_latex(
            justify="cr",
        )
        assert tex.splitlines()[2] == "\\begin{tabular}{ c r }"
        assert "caption" not in tex
        # with a title
        tex_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["val3", "val4"]],
            title="a title",
        )
        tex = tex_table.to_string(format="tex")
        tex = tex.splitlines()
        assert tex[-2] == "\\caption{a title}"

        tex = tex_table.to_string(format="tex", label="tab:first")
        tex = tex.splitlines()
        assert tex[-3] == "\\caption{a title}"
        assert tex[-2] == "\\label{tab:first}"

        # with a legend, no title
        tex_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["val3", "val4"]],
            legend="a legend",
        )
        tex = tex_table.to_string(format="tex")
        tex = tex.splitlines()
        # because it's treated as a title by default
        assert tex[-2] == "\\caption{a legend}"
        # unless you say not to
        tex = tex_table.to_string(format="tex", concat_title_legend=False)
        tex = tex.splitlines()
        assert tex[-2] == "\\caption*{a legend}"
        tex_table = make_table(
            header=["a", "b"],
            data=[["val1", "val2"], ["val3", "val4"]],
            title="a title.",
            legend="a legend",
        )
        tex = tex_table.to_string(format="tex")
        tex = tex.splitlines()
        assert tex[-2] == "\\caption{a title. a legend}"
        tex = tex_table.to_string(format="tex", concat_title_legend=False)
        tex = tex.splitlines()
        assert tex[2] == "\\caption{a title.}"
        assert tex[-2] == "\\caption*{a legend}"
        tex = tex_table.to_string(
            format="tex",
            concat_title_legend=False,
            label="table",
        )
        tex = tex.splitlines()
        assert tex[2] == "\\caption{a title.}"
        assert tex[3] == "\\label{table}"

    def test_to_html(self):
        """generates html table within c3table div"""
        # with no index_name, or title, or legend
        import re

        t = Table(header=self.t8_header, data=self.t8_rows)
        got = t.to_html()
        # make sure tags are matched
        for tag in ("div", "style", "table", "thead"):
            assert len(re.findall(f"<[/]*{tag}.*>", got)) == 2

        assert len(re.findall("<[/]*tr>", got)) == 4
        # 2 columns should be left aligned, 4 right aligned
        # adding 1 for the CSS style definition
        assert got.count("c3col_left") == 4 + 1
        assert got.count("c3col_right") == 8 + 1
        assert got.count("cell_title") == 1  # CSS defn only
        num_spans = got.count("span")
        num_caption = got.count("caption")

        t = Table(header=self.t8_header, data=self.t8_rows, title="a title")
        got = t.to_html()
        assert got.count("cell_title") == 2
        # number of spans increases by 2 to enclose the title
        assert got.count("span") == num_spans + 2
        assert got.count("caption") == num_caption + 2
        # no <br> element
        assert "<br>" not in got

        t = Table(header=self.t8_header, data=self.t8_rows, legend="a legend")
        got = t.to_html()
        assert got.count("cell_title") == 1
        # cell_legend not actually defined in CSS yet
        assert got.count("cell_legend") == 1
        # number of spans increases by 2 to enclose the title
        assert got.count("span") == num_spans + 2
        assert got.count("caption") == num_caption + 2
        # no <br> element
        assert "<br>" not in got

        t = Table(
            header=self.t8_header,
            data=self.t8_rows,
            title="a title",
            legend="a legend",
        )
        got = t.to_html()
        assert got.count("cell_title") == 2
        # cell_legend not actually defined in CSS yet
        assert got.count("cell_legend") == 1
        assert got.count("caption") == num_caption + 2
        # has <br> element
        assert "<br>" in got

    def test_invalid_format(self):
        """should raise value error"""
        t = make_table(self.t2_header, data=self.t2_rows)
        with pytest.raises(ValueError):
            t.format = "blah"

    def test_phylip(self):
        """generates phylip format"""
        rows = [
            ["a", "", 0.088337278874079342, 0.18848582712597683, 0.44084000179091454],
            ["c", 0.088337278874079342, "", 0.088337278874079342, 0.44083999937417828],
            ["b", 0.18848582712597683, 0.088337278874079342, "", 0.44084000179090932],
            ["e", 0.44084000179091454, 0.44083999937417828, 0.44084000179090932, ""],
        ]
        header = ["seq1/2", "a", "c", "b", "e"]
        dist = Table(header=header, data=rows, index_name="seq1/2")
        r = dist.to_string(format="phylip")
        r = r.splitlines()
        assert r[0].strip() == "4"
        for line in r[1:]:
            line = line.split()
            assert line[0] in dist.header
            assert line[-1][-1].isdigit()

        line = r[1].split()
        assert line[1] == "0.0000", line

    def test_load_table_invalid_type(self):
        """raises TypeError if filename invalid type"""
        with pytest.raises(TypeError):
            load_table({"a": [0, 1]})

    def test_make_table_white_space_in_column(self):
        """strips white space from column headers"""
        # matching header and data keys
        t = make_table(header=[" a"], data={" a": [0, 2]}, sep="\t")
        assert t.columns["a"].tolist() == [0, 2]
        assert isinstance(t.to_string(), str)

        # data key has a space
        t = make_table(data={" a": [0, 2]}, sep="\t")
        assert t.columns["a"].tolist() == [0, 2]
        assert isinstance(t.to_string(), str)

    def test_load_table_limit(self):
        """limit argument to function works"""
        t = load_table("data/sample.tsv", limit=2)
        assert t.shape[0] == 2

    def test_load_table_returns_static_columns(self):
        """for static data, load_table gives same dtypes for static_columns_type=True/False"""
        t = load_table("data/sample.tsv", sep="\t", static_column_types=False)
        is_false = {t.columns[c].dtype.name for c in t.columns}
        t = load_table("data/sample.tsv", sep="\t", static_column_types=True)
        is_true = {t.columns[c].dtype.name for c in t.columns}
        assert is_true == is_false

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
            assert isinstance(got, str)
            assert got != last
            assert got.startswith(startwith), f"{format}: {got[:10]}"
            last = got

    def test_grid_table_format(self):
        """test the table grid_table_format method"""
        from cogent3.format.table import grid_table_format

        formatted_grid = grid_table_format(
            self.t6_header,
            self.t6_rows,
            title="Test",
            legend="Units",
        )
        assert len(formatted_grid.split("\n")) == len(self.t6_rows) * 2 + 7

        formatted_grid = grid_table_format(
            self.t6_header,
            self.t6_rows,
            title="Really Long Title",
            legend="Extra Long Legend",
        )
        assert len(formatted_grid.split("\n")) == len(self.t6_rows) * 2 + 7 + 2

    def test_to_markdown(self):
        """Exercising the table markdown method"""
        table = make_table(self.t6_header, self.t6_rows, format="md")
        markdown_table = table.to_markdown(justify="crl")
        markdown_list = markdown_table.split("\n")
        assert markdown_list[2].count("|") == 5
        # the pipe symbol should have been escaped
        assert markdown_list[2].count("\\|") == 1

        with pytest.raises(ValueError):
            _ = table.to_markdown(justify="cr1")

    def test_to_csv(self):
        """successfully create csv formatted string"""
        table = Table(
            header=self.t3_header,
            data=self.t3_rows,
            title="A title",
            legend="A legend",
        )
        sv = table.to_csv()
        expect = ["id,foo,bar", "6,abc,66", "7,bca,77"]
        assert sv.splitlines() == expect
        sv = table.to_csv(with_title=True)
        assert sv.splitlines() == ["A title", *expect]
        sv = table.to_csv(with_legend=True)
        assert sv.splitlines() == [*expect, "A legend"]
        sv = table.to_csv(with_title=True, with_legend=True)
        assert sv.splitlines() == ["A title", *expect, "A legend"]

    def test_to_tsv(self):
        """successfully create csv formatted string"""
        table = Table(
            header=self.t3_header,
            data=self.t3_rows,
            title="A title",
            legend="A legend",
        )

        sv = table.to_tsv()
        expect = ["id\tfoo\tbar", "6\tabc\t66", "7\tbca\t77"]
        assert sv.splitlines() == expect
        sv = table.to_tsv(with_title=True)
        assert sv.splitlines() == ["A title", *expect]
        sv = table.to_tsv(with_legend=True)
        assert sv.splitlines() == [*expect, "A legend"]
        sv = table.to_tsv(with_title=True, with_legend=True)
        assert sv.splitlines() == ["A title", *expect, "A legend"]

    def test_to_delim(self):
        """successfully create separated format with arbitrary character"""
        table = Table(
            header=self.t3_header,
            data=self.t3_rows,
        )
        sv = table.to_string(sep=";")
        expect = ["id;foo;bar", "6;abc;66", "7;bca;77"]
        assert sv.splitlines() == expect

    def test_to_rst_grid(self):
        """generates a rst grid table"""
        table = Table(header=["a", "b"], data=[[1, 2]], title="A title")
        got = table.to_rst(csv_table=False).splitlines()
        assert table.title in got[1]
        assert set(got[0]) == {"-", "+"}
        assert set(got[4]) == {"=", "+"}

    def test_to_rst_csv(self):
        """generates a rst csv-table"""
        table = Table(
            header=["a", "b"],
            data=[[1, 2]],
            title="A title",
            legend="A legend",
        )
        got = table.to_rst(csv_table=True)
        assert got.splitlines() == [
            ".. csv-table:: A title A legend",
            '    :header: "a", "b"',
            "",
            "    1, 2",
        ]
        # try without a title/legend
        table = Table(header=["a", "b"], data=[[1, 2]])
        got = table.to_rst(csv_table=True)
        assert got.splitlines() == [
            ".. csv-table::",
            '    :header: "a", "b"',
            "",
            "    1, 2",
        ]

    def test_get_repr_(self):
        """handles single column case"""
        t = make_table(self.t2_header, data=self.t2_rows)
        t = t[:, 0]
        # the next line was previously failing
        t._get_repr_()

        table = Table(header=["a", "b"], data=[[1, 2]])
        table, _, unset_columns = table._get_repr_()
        assert table.shape == (1, 2)
        assert unset_columns is None

        table = make_table(header=["a", "b"])
        table.columns["a"] = ["a"]
        table, _, unset_columns = table._get_repr_()
        assert table.shape == (1, 1)
        assert "b" in unset_columns

    def test_repr_html_(self):
        """should produce html"""
        # no index_name
        t = Table(header=self.t8_header, data=self.t8_rows)
        _ = t._repr_html_()

        # with an index_name
        t = Table(header=self.t8_header, data=self.t8_rows, index_name="edge.name")
        got = t._repr_html_()
        # and the index_name column should contain "index_name" css class
        assert got.count("index") == t.shape[0] + 1  # add 1 for CSS style sheet

        # data has tuples in an array
        data = {
            "key": numpy.array([("a", "c"), ("b", "c"), ("a", "d")], dtype="O"),
            "count": [1, 3, 2],
        }
        t = Table(data=data)
        _ = t._repr_html_()

        # some columns without data
        table = make_table(header=["a", "b"])
        table.columns["a"] = ["a"]
        _ = t._repr_html_()

        # single column with a single value should not fail
        table = make_table(data={"kappa": [3.2]}, title="a title")
        _ = table._repr_html_()

        # set head and tail, introduces ellipsis row class
        table = make_table(data={"A": list("abcdefghijk"), "B": list(range(11))})
        table.set_repr_policy(head=8, tail=1)
        got = table._repr_html_().splitlines()
        num_rows = 0
        for l in got:
            if "<tr>" in l:
                num_rows += 1
                if "ellipsis" in l:
                    break

        assert num_rows == 9

    def test_array(self):
        """should produce array"""
        # data has tuples in an array
        data = {
            "key": numpy.array([("a", "c"), ("b", "c"), ("a", "d")], dtype="O"),
            "count": [1, 3, 2],
        }
        expect = [list(v) for v in zip(data["key"][:], data["count"], strict=False)]
        t = Table(data=data)
        arr = t.array
        assert_equal(arr.tolist(), expect)

    def test_separator_format(self):
        """testing separator_format with title and legend, and contents that match the separator"""
        from cogent3.format.table import separator_format

        with pytest.raises(RuntimeError):
            _ = separator_format(self.t6_header, self.t6_rows)
        separated_table = separator_format(
            self.t6_header,
            self.t6_rows,
            sep=" | ",
            title="Test",
            legend="Units",
        )
        assert len(separated_table.split("\n")) == len(self.t6_rows) + 3

    def test_separator_format_writer(self):
        """exercising separator_format_writer"""
        from cogent3.format.table import separator_formatter

        t3 = Table(header=self.t3_header, data=self.t3_rows)
        comma_sep = t3.to_string(sep=",").splitlines()
        writer = separator_formatter(sep=" | ")
        formatted = list(writer([l.split(",") for l in comma_sep], has_header=True))
        expected_format = ["id | foo | bar", "6 | abc | 66", "7 | bca | 77"]
        assert formatted == expected_format

    def test_set_repr_policy(self):
        """exercising setting repr policy"""
        t = Table(header=self.t2_header, data=self.t2_rows)
        t.set_repr_policy(random=2)
        r = repr(t)
        assert isinstance(r, str)
        r, _, _ = t._get_repr_()
        assert r.shape[0] == 2
        t.set_repr_policy(head=1)
        r, _, _ = t._get_repr_()
        assert r.shape[0] == 1
        t.set_repr_policy(tail=3)
        r, _, _ = t._get_repr_()
        assert r.shape[0] == 3
        t.set_repr_policy(show_shape=False)
        r = repr(t)
        assert f"\n{t.shape[0]:,} rows x {t.shape[1]:,} columns" not in r
        r = t._repr_html_()
        assert f"\n{t.shape[0]:,} rows x {t.shape[1]:,} columns" not in r

    def test_head(self):
        """returns the head of the table!"""
        from cogent3.core import table

        display = table.display
        head = TrapOutput()
        table.display = head
        t = Table(header=self.t1_header, data=self.t1_rows)
        t.head(nrows=3)
        assert head.data.shape[0] == 3
        assert len(head.output.splitlines()) == 9
        assert head.data.to_list() == self.t1_rows[:3]
        # tests when number of rows < default
        t = make_table(data={"a": ["a"], "b": ["b"]})
        t.head()
        assert head.data.shape[0] == 1
        assert len(head.output.splitlines()) == 7
        assert head.data.to_list() == [["a", "b"]]
        table.display = display

    def test_tail(self):
        """returns the tail of the table!"""
        from cogent3.core import table

        display = table.display
        tail = TrapOutput()
        table.display = tail
        t = Table(header=self.t1_header, data=self.t1_rows)
        t.tail(nrows=3)
        assert tail.data.shape[0] == 3
        assert len(tail.output.splitlines()) == 9
        assert [int(v) for v in tail.data[:, -1].to_list()] == [
            r[-1] for r in self.t1_rows[-3:]
        ]
        # tests when number of rows < default
        t = make_table(data={"a": ["a"], "b": ["b"]})
        t.tail()
        assert tail.data.shape[0] == 1
        assert len(tail.output.splitlines()) == 7
        assert tail.data.to_list() == [["a", "b"]]
        table.display = display

    @skipIf(DataFrame is None, "pandas not installed")
    def test_to_pandas(self):
        """produces a dataframe"""
        t = Table(header=self.t1_header, data=self.t1_rows)
        df = t.to_pandas()
        assert isinstance(df, DataFrame)
        data = df.to_numpy()
        assert data.tolist() == self.t1_rows

    def test_load_table(self):
        """exercising load table"""
        path = os.path.dirname(os.path.dirname(__file__))
        path = os.path.join(path, "data/sample.tsv")
        table = load_table(path)
        assert table.shape == (10, 3)

    def test_cast_str_to_numerical(self):
        """correctly converts a series of strings to numeric values"""
        d = arange(4, step=0.1)
        r = cast_str_to_numeric(d)
        assert_equal(d, r)

        with numpy.testing.suppress_warnings() as sup:
            # we know that converting to real loses imaginary
            sup.filter(ComplexWarning)
            for d_type in [numpy.int64, numpy.complex128, numpy.float64]:
                d = d.astype(d_type)
                r = cast_str_to_numeric(d)
                assert isinstance(r[0], type(d[0]))

        d = d.astype(str)
        r = cast_str_to_numeric(d)
        assert isinstance(r[0], numpy.float64)
        d = numpy.array(d, dtype="U")
        r = cast_str_to_numeric(d)
        assert isinstance(r[0], numpy.float64)
        d = numpy.array(d, dtype="S")
        r = cast_str_to_numeric(d)
        assert isinstance(r[0], numpy.float64)

    def test_cast_str_to_array(self):
        """handle processing string series"""
        d = [".123|.345", "123"]
        r = cast_str_to_array(d, static_type=True)
        assert "str" in r.dtype.name
        r = cast_str_to_array(d, static_type=False)
        assert r.dtype.name == "object"
        d = [".123|.345", "123", "()"]
        r = cast_str_to_array(d, static_type=False)
        assert r[-1] == ()

    def test_filtering_parser(self):
        """filters rows"""
        expect = []
        for r in self.t1_rows:
            row = [str(e) for e in r]
            expect.append(row)

        t = make_table(self.t1_header, data=self.t1_rows)
        lines = t.to_csv().splitlines()
        # no limit set
        reader = FilteringParser(
            row_condition=lambda x: x[0] == "A",
            with_header=True,
            sep=",",
        )
        got = list(reader(lines))
        assert got[0] == self.t1_header
        assert got[1:] == [r for r in expect if r[0] == "A"]

        # limit set
        reader = FilteringParser(
            lambda x: x[0] == "A",
            with_header=True,
            sep=",",
            limit=2,
        )
        got = list(reader(lines))
        assert got[0] == self.t1_header
        assert got[1:] == [r for r in expect if r[0] == "A"][:2]

        # negate
        reader = FilteringParser(
            lambda x: x[0] == "A",
            negate=True,
            with_header=True,
            sep=",",
        )
        got = list(reader(lines))
        assert got[0] == self.t1_header
        assert got[1:] == [r for r in expect if r[0] == "X"]

        # parser works with load_table
        path = TEST_ROOT / "data" / "sample.tsv"
        reader = FilteringParser(lambda x: x[0] == "A", with_header=True, sep="\t")
        table = load_table(path, reader=reader)
        assert list(table.header) == self.t1_header
        assert table.array.tolist() == [r for r in self.t1_rows if r[0] == "A"]

        # parser works if called on path as Path
        reader = FilteringParser(lambda x: x[0] == "A", with_header=True, sep="\t")
        got = list(reader(path))
        assert len(got) == 7

        # parser works if called on path as str
        got = list(reader(str(path)))
        assert len(got) == 7

        # parser works with no conditions
        reader = FilteringParser(with_header=True, sep="\t")
        t = load_table(path, reader=reader)
        assert t.shape == (10, 3)

    def test_filtering_parser_filter_columns(self):
        """filters columns"""
        path = TEST_ROOT / "data" / "sample.tsv"
        # specified by int
        reader = FilteringParser(columns=0, with_header=True, sep="\t")
        got = load_table(path, reader=reader)
        assert got.shape == (10, 1)

        # specified by str
        reader = FilteringParser(columns="length", with_header=True, sep="\t")
        got = load_table(path, reader=reader)
        assert got.shape == (10, 1)

        # specified by index_name
        reader = FilteringParser(columns=[0, 2], with_header=True, sep="\t")
        got = load_table(path, reader=reader)
        assert got.shape == (10, 2)

        # specified by name
        reader = FilteringParser(
            columns=["chrom", "length"],
            with_header=True,
            sep="\t",
        )
        got2 = load_table(path, reader=reader)
        assert got2.shape == (10, 2)
        assert_equal(got.array, got2.array)

        # raises value error if column name doesn't exist
        reader = FilteringParser(columns=["blah", "length"], with_header=True, sep="\t")
        with pytest.raises(ValueError):
            _ = load_table(path, reader=reader)

        # raises IndexError if column index_name doesn't exist
        reader = FilteringParser(columns=[0, 10], with_header=True, sep="\t")
        with pytest.raises(IndexError):
            _ = load_table(path, reader=reader)

        # raises ValueError if names given and with_header is False
        with pytest.raises(ValueError):
            _ = FilteringParser(columns=["blah"], with_header=False)

    def test_set_column_format(self):
        """fails if invalid format spec provided"""
        data = {
            "Gene": [
                "NP_003077_hs_mm_rn_dna",
                "NP_004893_hs_mm_rn_dna",
                "NP_005079_hs_mm_rn_dna",
                "NP_005500_hs_mm_rn_dna",
                "NP_055852_hs_mm_rn_dna",
            ],
            "Type": ["Con", "Con", "Con", "Con", "Con"],
            "LR": [
                2.5386013224378985,
                121351.42635634111,
                9516594.978886133,
                7.382703020266491e-08,
                10933217.708952725,
            ],
        }
        t = make_table(data=data)
        with pytest.raises(ValueError):
            t.format_column("LR", ".4e")

    def test_table_format(self):
        """table formating function doesn't fail with different input values"""
        from cogent3.format.table import formatted_cells

        data = [[230, "acdef", None], [6, "cc", 1.9876]]
        head = ["one", "two", "three"]
        _ = formatted_cells(data, header=head)
        data = [[230, "acdef", 1.3], [6, "cc", numpy.array([1.9876, 2.34])]]
        _ = formatted_cells(data, header=head)

    def test_to_categorical(self):
        """correctly construct contingency table"""
        data = {"Ts": [31, 58], "Tv": [36, 138], "": ["syn", "nsyn"]}
        table = make_table(header=["", "Ts", "Tv"], data=data)
        with pytest.raises(ValueError):
            # did not set an index_name
            table.to_categorical(columns=["Ts", "Tv"])

        got = table.to_categorical(columns=["Ts", "Tv"], index_name="")
        assert_equal(got.observed, table[:, 1:].array)

        got = table.to_categorical(["Ts"])
        mean = got.observed.array.mean()
        expected = numpy.array([[mean], [mean]])
        assert_equal(got.expected, expected)

        # works if index_name included
        got = table.to_categorical(columns=["Ts", "Tv", ""])
        assert_equal(got.observed, table[:, 1:].array)

        # works if no columns specified
        got = table.to_categorical()
        assert_equal(got.observed, table[:, 1:].array)

        data = {
            "": numpy.array(["syn", "nsyn"], dtype=object),
            "Ts": numpy.array([31, 58], dtype=object),
            "Tv": numpy.array([36, 138], dtype=object),
        }

        table = make_table(header=["", "Ts", "Tv"], data=data, index_name="")
        with pytest.raises(TypeError):
            table.to_categorical(columns=["Ts", "Tv"])

    def test_is_html_markup(self):
        """format function confirms correctly specified html"""
        assert is_html_markup("<table>blah</table>")
        assert is_html_markup("<table>blah<table>blah</table></table>")
        assert is_html_markup("<i>blah</i><sub>blah</sub>")
        assert is_html_markup("<i>blah</i>\n<sub>blah</sub>")
        assert not is_html_markup("<table>blah</tabl>")
        assert not is_html_markup("<table>")
        assert not is_html_markup("blah < blah")
        assert not is_html_markup("blah > blah")


@pytest.mark.parametrize("data", [[[0, 1, 2], [0, 1]], [[0, 1, 2], [0, 1, 2, 3]]])
def test_mixed_row_lengths(data):
    """should fail"""
    with pytest.raises(ValueError):
        _ = Table(header=["a", "b", "c"], data=data)


@pytest.fixture
def t2():
    header = ["id", "foo", "bar"]
    rows = [
        [1, "abc", 11],
        [2, "bca", 22],
        [3, "cab", 33],
        [4, "abc", 44],
        [5, "bca", 55],
    ]
    return Table(header=header, data=rows, index_name="id")


@pytest.fixture
def t3():
    header = ["id", "foo2", "bar2"]
    rows = [[6, "abc", 66], [7, "bca", 77]]
    return Table(header=header, data=rows, index_name="id")


@pytest.fixture
def t4():
    header = ["id", "foo2", "bar2"]
    rows = [[6, "abc", 66], [7, "bca", 77]]
    return Table(header=header, data=rows, index_name="id")


def test_inner_join_col_naming(t2, t3):
    """test the table joined method"""
    # inner join with defaults
    got = t2.inner_join(t3)
    expect = list(t2.header) + [f"right_{c}" for c in t3.header[1:]]
    assert list(got.header) == expect

    # inner join with no col_prefix
    got = t2.inner_join(t3, col_prefix="")
    expect = list(t2.header) + list(t3.header[1:])
    assert list(got.header) == expect


def test_outer_join_col_naming(t2, t3):
    """test the table joined method"""
    # cross join with a new col prefix
    col_prefix = "cp_"
    got = t2.cross_join(t3, col_prefix=col_prefix)
    expect = list(t2.header) + [f"{col_prefix}{c}" for c in t3.header]
    assert list(got.header) == expect
    assert got.shape[0] == t2.shape[0] * t3.shape[0]


def test_repr_html_continuation():
    # should have an ellipsis row for every c ontinued table
    t8_header = ["edge.name", "edge.parent", "length", "x", "y", "z"]
    t8_rows = [
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
        ["NineBande", "root", 4.0, 1.0, 3.0, 6.0],
    ]
    table = make_table(header=t8_header, data=t8_rows, max_width=30, title="a title")
    # with 30 characters wide there are 3 subtables, so we expect 3 ellipsis rows
    table.set_repr_policy(head=2, tail=2)
    assert table._repr_html_().count('<tr class="ellipsis">') == 3


def test_pickle_unpickle(tmp_path):
    """roundtrip via pickling"""
    data = {
        "edge.parent": {
            "NineBande": "root",
            "edge.1": "root",
        },
        "x": {
            "NineBande": 1.0,
            "edge.1": 1.0,
        },
        "length": {
            "NineBande": 4.0,
            "edge.1": 4.0,
        },
        "y": {
            "NineBande": 3.0,
            "edge.1": 3.0,
        },
        "z": {
            "NineBande": 6.0,
            "edge.1": 6.0,
        },
        "edge.name": {
            "NineBande": "NineBande",
            "edge.1": "edge.1",
        },
    }
    t = Table(
        data=data,
        max_width=50,
        index_name="edge.name",
        title="My title",
        legend="blah",
    )
    # via string
    s = pickle.dumps(t)
    r = pickle.loads(s)  # noqa: S301
    assert str(t) == str(r)
    # via file
    path = tmp_path / "table.pickle"
    t.write(str(path))
    r = load_table(path)
    assert str(t) == str(r)


def test_load_mixed(tmp_path):
    """load data with mixed data type columns"""
    t = Table(
        header=["abcd", "data", "float"],
        data=[[str([1, 2, 3, 4, 5]), "0", 1.1], ["x", 5.0, 2.1], ["y", "", 3.1]],
    )
    path = tmp_path / "table.tsv"
    t.write(str(path))
    r = load_table(path)
    assert str(t) == str(r)
    assert "float", r.columns["float"].dtype.name


def test_load_mixed_static(tmp_path):
    """load data, mixed data type columns remain as string"""
    t = make_table(header=["A", "B"], data=[[1, 1], ["a", 2]])
    path = tmp_path / "table.txt"
    t.write(str(path), sep="\t")
    # if static types, then mixed columns become strings
    r = load_table(path, sep="\t", static_column_types=True)
    assert "str" in r.columns["A"].dtype.name


def test_load_mixed_row_lengths(tmp_path):
    """skip_inconsistent skips rows that have different length to header"""
    h = list("ABCDE")
    r = [list("12345"), list("000"), list("12345")]
    text = "\n".join(["\t".join(l) for l in [h, *r]])
    path = tmp_path / "table.tsv"
    with open(path, "w") as out:
        out.write(text)
    r = load_table(path, skip_inconsistent=True)
    assert r.shape == (2, 5)
    assert r.header == tuple(h)
    assert r.array.tolist() == [list(range(1, 6))] * 2
    # loading without skip_inconsistent raise ValueError
    with pytest.raises(ValueError):
        r = load_table(path, skip_inconsistent=False)


def test_write_to_json(tmp_path):
    """tests writing to json file"""
    t = load_table("data/sample.tsv")
    path = tmp_path / "table.json"
    t.write(path)
    with open_(path) as fn:
        got = json.loads(fn.read())
        assert got["type"] == get_object_provenance(Table)
        data = got["data"]
        assert tuple(data["order"]) == t.header
        assert t.shape == (
            len(next(iter(data["columns"].items()))[1]["values"]),
            len(data["columns"]),
        )
        assert t.array.T.tolist() == [v["values"] for v in data["columns"].values()]


def test_write_compressed(tmp_path):
    """tests writing to compressed format"""
    t = load_table("data/sample.tsv")
    with open("data/sample.tsv") as infile:
        expect = infile.read()

    path = tmp_path / "table.txt"
    # using the compressed option
    t.write(path, sep="\t", compress=True)
    with open_(f"{path}.gz") as infile:
        got = infile.read()
    assert got == expect

    # specifying via a suffix
    t.write(f"{path}.gz", sep="\t")
    with open_(f"{path}.gz") as infile:
        got = infile.read()
    assert got == expect


def test_load_table_from_json(tmp_path):
    """tests loading a Table object from json file"""
    json_path = tmp_path / "table.json"
    t = load_table("data/sample.tsv")
    t.write(json_path)

    got = load_table(json_path)
    assert got.shape == t.shape
    assert got.header == t.header
    assert_equal(got.array, t.array)


def test_load_table_filename_case(tmp_path):
    """load_table insensitive to file name case"""
    outpath = tmp_path / "temp.CSV"
    outpath.write_text("a,b,c\n0,2,abc\n1,3,efg")

    table = load_table(outpath)
    data = table.columns.to_dict()
    assert data == {"a": [0, 1], "b": [2, 3], "c": ["abc", "efg"]}
