import json
import os

from tempfile import TemporaryDirectory
from unittest import TestCase, main

import numpy

from numpy.testing import assert_allclose

from cogent3 import DNA
from cogent3.util.dict_array import (
    DictArray,
    DictArrayTemplate,
    convert2DDict,
    convert2Ddistance,
    convert_1D_dict,
    convert_dict,
    convert_for_dictarray,
    convert_series,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class DictArrayTest(TestCase):
    a = numpy.identity(3, int)

    def test_convert_series(self):
        """convert_series produces valid template input"""
        vals, row_keys, col_keys = convert_series([[4], [5]], ["A", "B"], ["a"])
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[4], [5]])
        data = [[245, 599]]
        vals, row_keys, col_keys = convert_series(data)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), data)

        vals, row_keys, col_keys = convert_series(data[0])
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), data[0])

    def test_convert_dict(self):
        """convert_dict produces valid template input"""
        twoDdict = dict(a=dict(b=4, c=5))
        vals, row_keys, col_keys = convert_dict(twoDdict)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[4, 5]])

    def test_convert2DDict(self):
        """convert2DDict produces valid template input"""
        data = dict(a=dict(b=4, c=5))
        vals, row_keys, col_keys = convert2DDict(data)
        self.assertEqual(set(row_keys), set(["a"]))
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[4, 5]])
        # row keys, then column
        self.assertEqual(b.template.names, [["a"], ["b", "c"]])

        data = {
            "a": {"a": 0, "b": 1, "e": 0},
            "b": {"a": 1, "b": 0, "e": 4},
            "e": {"a": 0, "b": 4, "e": 0},
        }
        vals, row_keys, col_keys = convert2DDict(data)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        got = b.to_dict()
        self.assertEqual(got, data)
        self.assertEqual(b.template.names, [["a", "b", "e"], ["a", "b", "e"]])

        data = dict(a=dict(b=4, c=5))
        vals, row_keys, col_keys = convert2DDict(data, make_symmetric=True)
        self.assertEqual(row_keys, col_keys)
        self.assertEqual(vals, [[0, 4, 5], [4, 0, 0], [5, 0, 0]])

    def test_convert2Ddistance(self):
        """convert2Ddistance produces valid template input"""
        data = {("a", "b"): 4, ("a", "c"): 5}
        vals, row_keys, col_keys = convert2Ddistance(data)
        b = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(b.array.tolist(), [[0, 4, 5], [4, 0, 0], [5, 0, 0]])

    def test_convert_1D_dict(self):
        """convert_1D_dict produces valid template input"""
        data = dict(a=0, b=35, c=45)
        vals, keys = convert_1D_dict(data)
        b = DictArrayTemplate(keys)
        b = b.wrap(vals)
        self.assertEqual(b.array.tolist(), [0, 35, 45])

    def test_construct_both_dim_str(self):
        """correctly construct when keys for both dimensions are str"""
        b = DictArrayTemplate("abc", "ABC").wrap(self.a)
        self.assertEqual(b[0].array.tolist(), [1, 0, 0])
        self.assertEqual(b["a"].array.tolist(), [1, 0, 0])
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    def test_key_levels(self):
        """DictArray both levels have keys."""
        b = DictArrayTemplate("abc", "ABC").wrap(self.a)
        self.assertEqual(b.keys(), ["a", "b", "c"])
        self.assertEqual(b["a"].keys(), ["A", "B", "C"])
        self.assertEqual(list(b["a"]), [1, 0, 0])
        self.assertEqual(sum(b["a"]), 1)

    def test_int_labels(self):
        """DictArray with no labels."""
        b = DictArrayTemplate(3, 3).wrap(self.a)
        self.assertEqual(b.keys(), [0, 1, 2])
        self.assertEqual(b[0].keys(), [0, 1, 2])
        self.assertEqual(sum(b[0]), 1)

    def test_str_labels(self):
        """DictArray with str labels or numpy U dtype"""
        b = DictArrayTemplate(["Ab", "Bb", "Cb"], ["Db", "Eb", "Fb"]).wrap(self.a)
        self.assertEqual(b.keys(), ["Ab", "Bb", "Cb"])
        self.assertEqual(b[0].keys(), ["Db", "Eb", "Fb"])
        self.assertEqual(b["Ab", "Eb"], 0)
        self.assertEqual(b["Bb", "Eb"], 1)

        b = DictArrayTemplate(
            numpy.array(["Ab", "Bb", "Cb"], dtype="U"),
            numpy.array(["Db", "Eb", "Fb"], dtype="U"),
        ).wrap(self.a)
        self.assertEqual(b.keys(), ["Ab", "Bb", "Cb"])
        self.assertEqual(b[0].keys(), ["Db", "Eb", "Fb"])
        self.assertEqual(b["Ab", "Eb"], 0)
        self.assertEqual(b["Bb", "Eb"], 1)

    def test_with_mixed_label_types(self):
        """DictArray constructed with mixed label types."""
        b = DictArrayTemplate("ABC", 3).wrap(self.a)
        self.assertEqual(b.keys(), ["A", "B", "C"])
        self.assertEqual(b["A"].keys(), [0, 1, 2])

    def test_numpy_ops(self):
        """DictArray should work properly in numpy operations."""
        darr = DictArrayTemplate(list(DNA), list(DNA)).wrap(
            [
                [0.7, 0.1, 0.1, 0.1],
                [0.1, 0.7, 0.1, 0.1],
                [0.1, 0.1, 0.7, 0.1],
                [0.1, 0.1, 0.1, 0.7],
            ]
        )
        mprobs = numpy.array([0.25, 0.25, 0.25, 0.25])
        assert_allclose(mprobs.dot(darr), [0.25, 0.25, 0.25, 0.25])
        assert_allclose(numpy.dot(mprobs, darr), [0.25, 0.25, 0.25, 0.25])

    def test_to_dict(self):
        """DictArray should convert 1D / 2D arrays with/without named row"""
        # 1D data, only 1D keys provided
        data = [0, 35, 45]
        keys = "a", "b", "c"
        darr = DictArrayTemplate(keys).wrap(data)
        self.assertEqual(darr.to_dict(), dict(zip(keys, data)))
        # 2D data, 2D keys, both string, provided
        data = [[0, 35, 45]]
        darr = DictArrayTemplate(["0"], keys).wrap(data)
        darr.to_dict()
        self.assertEqual(darr.to_dict(), {"0": {"a": 0, "b": 35, "c": 45}})
        # 2D data, 2D keys, one int, one string, provided
        darr = DictArrayTemplate([1], keys).wrap(data)
        self.assertEqual(darr.to_dict(), {1: {"a": 0, "b": 35, "c": 45}})
        darr = DictArrayTemplate([0], keys).wrap(data)
        self.assertEqual(darr.to_dict(), {0: {"a": 0, "b": 35, "c": 45}})

    def test_to_dict_1d(self):
        """should successfully produce a 1D dict"""
        data = {
            "ABAYE2984": {
                "ABAYE2984": 0,
                "Atu3667": None,
                "Avin_42730": 0.6381173875591908,
                "BAA10469": None,
            },
            "Atu3667": {
                "ABAYE2984": None,
                "Atu3667": 0,
                "Avin_42730": 2.3682377869318993,
                "BAA10469": None,
            },
            "Avin_42730": {
                "ABAYE2984": 0.6381173875591908,
                "Atu3667": 2.3682377869318993,
                "Avin_42730": 0,
                "BAA10469": 1.8515731266342546,
            },
            "BAA10469": {
                "ABAYE2984": None,
                "Atu3667": None,
                "Avin_42730": 1.8515731266342546,
                "BAA10469": 0,
            },
        }
        darr = DictArray(data, dtype="O")
        expect = {
            (n1, n2): darr[n1, n2]
            for n1 in darr.template.names[0]
            for n2 in darr.template.names[1]
        }
        self.assertEqual(darr.to_dict(flatten=True), expect)

        darr = DictArrayTemplate(["s1", "s2"], ["s1", "s2"]).wrap(
            [[0.0, 0.25], [0.25, 0.0]]
        )
        self.assertEqual(
            darr.to_dict(flatten=True),
            {
                ("s1", "s2"): 0.25,
                ("s2", "s1"): 0.25,
                ("s1", "s1"): 0.0,
                ("s2", "s2"): 0.0,
            },
        )

    def test_to_dict_nested(self):
        """DictArray.to_dict() should convert nested DictArray instances to
        dict's too."""
        a = numpy.identity(3, int)
        b = DictArrayTemplate("abc", "ABC")
        b = b.wrap(a)
        self.assertEqual(b.array.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        c = DictArrayTemplate("de", "DE").wrap([[b, b], [b, b]])
        self.assertTrue(isinstance(c.to_dict()["d"], dict))

    def test_to_dict_values(self):
        """values from to_dict should be python types"""
        keys = "a", "b", "c", "d"
        for data, _type in [
            ([0, 35, 45, 3], int),
            (["abc", "def", "jkl;", "aa"], str),
            ([0.1, 0.2, 0.3, 0.4], float),
        ]:
            darr = DictArrayTemplate(keys).wrap(data)
            got = {type(v) for v in darr.to_dict().values()}
            self.assertEqual(got, {_type})

        for data, _type in [
            ([0, 35, 45, 3], int),
            (["abc", "def", "jkl;", "aa"], str),
            ([0.1, 0.2, 0.3, 0.4], float),
        ]:
            darr = DictArrayTemplate(keys[:2], keys[2:]).wrap([data[:2], data[2:]])
            got = {type(v) for d in darr.to_dict().values() for v in d.values()}
            self.assertEqual(got, {_type})

    def test_to_dict_json(self):
        """should be able to json.dumps result of to_dict"""
        keys = "a", "b", "c", "d"
        for data in [
            [0, 35, 45, 3],
            ["abc", "def", "jkl;", "aa"],
            [0.1, 0.2, 0.3, 0.4],
        ]:
            darr = DictArrayTemplate(keys).wrap(data)
            got = json.dumps(darr.to_dict())
            self.assertIsInstance(got, str)

        for data in [
            [0, 35, 45, 3],
            ["abc", "def", "jkl;", "aa"],
            [0.1, 0.2, 0.3, 0.4],
        ]:
            darr = DictArrayTemplate(keys[:2], keys[2:]).wrap([data[:2], data[2:]])
            got = json.dumps(darr.to_dict())
            self.assertIsInstance(got, str)

    def test_to_dict_roundtrip(self):
        """roundtrip of DictArray.to_dict() should produce same order."""
        d1 = dict(a=dict(k=1, l=2, m=3), b=dict(k=4, l=5, m=6))
        darr1 = DictArray(d1)
        d2 = darr1.to_dict()
        self.assertEqual(d1, d2)
        d3 = DictArray(d2)
        self.assertEqual(d1, d3)

    def test_convert_for_dictarray(self):
        """successfully delegates when constructed from a DictArray"""
        a = numpy.identity(3, int)
        b = DictArrayTemplate("abc", "ABC").wrap(a)
        vals, row_keys, col_keys = convert_for_dictarray(b)
        got = DictArrayTemplate(row_keys, col_keys).wrap(vals)
        self.assertEqual(got.array.tolist(), b.array.tolist())
        # the wrap method creates a new array
        self.assertIsNot(got.array, b.array)

    def test_convert_for_dictarray2(self):
        """convert_for_dictarray correctly delegates"""
        b = DictArrayTemplate("abc", "ABC").wrap(self.a)
        data_types = (
            [[245, 599]],
            dict(a=dict(b=4, c=5)),
            {("a", "b"): 4, ("a", "c"): 5},
            dict(a=0, b=35, c=45),
            b,
        )
        for data in data_types:
            vals, row_keys, col_keys = convert_for_dictarray(data)
            _ = DictArrayTemplate(row_keys, col_keys).wrap(vals)

    def test_direct_construction(self):
        """directly construct a dict array"""
        b = DictArrayTemplate("abc", "ABC").wrap(self.a)
        data_types = (
            [[245, 599]],
            dict(a=dict(b=4, c=5)),
            {("a", "b"): 4, ("a", "c"): 5},
            dict(a=0, b=35, c=45),
            b,
        )
        # exercising construction
        for data in data_types:
            _ = DictArray(data)

    def test_getitem(self):
        """correctly slices"""
        darr = DictArrayTemplate(list(DNA), list(DNA)).wrap(
            [
                [0.7, 0.1, 0.1, 0.1],
                [0.1, 0.7, 0.1, 0.1],
                [0.1, 0.1, 0.7, 0.1],
                [0.1, 0.1, 0.1, 0.7],
            ]
        )
        r = darr[:, "A":"G"]
        assert_allclose(r.to_array(), [[0.1], [0.1], [0.7], [0.1]])
        r = darr[2:, "A":"G"]
        assert_allclose(r.to_array(), [[0.7], [0.1]])

    def test_to_normalized(self):
        """computes frequencies across correct dimension"""
        data = [[3, 7], [2, 8], [5, 5]]
        darr = DictArrayTemplate(list("ABC"), list("ab")).wrap(data)
        row_normal = darr.to_normalized(by_row=True)
        assert_allclose(row_normal.array, [[0.3, 0.7], [0.2, 0.8], [0.5, 0.5]])
        col_normal = darr.to_normalized(by_column=True)
        assert_allclose(col_normal.array, [[0.3, 7 / 20], [0.2, 8 / 20], [0.5, 5 / 20]])
        # trying to do both raises AssertionError
        with self.assertRaises(AssertionError):
            darr.to_normalized(by_row=True, by_column=True)

    def test_col_sum(self):
        """correctly computes column sum"""
        data = [[3, 7], [2, 8], [5, 5]]
        darr = DictArrayTemplate(list("ABC"), list("ab")).wrap(data)
        col_sum = darr.col_sum()
        assert_allclose(col_sum.array, [10, 20])

    def test_row_sum(self):
        """correctly computes row sum"""
        data = [[3, 7], [2, 8], [5, 5]]
        darr = DictArrayTemplate(list("ABC"), list("ab")).wrap(data)
        row_sum = darr.row_sum()
        assert_allclose(row_sum.array, [10, 10, 10])

    def test_get_repr_html(self):
        """exercising method used by parent classes for nice Jupyter display"""
        data = [[3, 7], [2, 8], [5, 5]]
        darr = DictArrayTemplate(list("ABC"), list("ab")).wrap(data)
        got = darr._repr_html_()
        self.assertIsInstance(got, str)
        self.assertTrue(len(got), 100)

        # case where 1D array
        a = [4, 6, 4, 2]
        darr = DictArrayTemplate(["A", "C", "G", "T"]).wrap(a)
        got = darr._repr_html_()
        self.assertTrue('class="index"' not in got)

        # case of 3D array
        d3 = numpy.arange(8).reshape((2, 2, 2))
        darr = DictArrayTemplate(2, 2, 2).wrap(d3)
        got = darr._repr_html_()
        self.assertIn("3 dimensional", got)

    def test_write(self):
        """exercising write method"""
        data = [[3, 7], [2, 8], [5, 5]]
        darr = DictArrayTemplate(list("ABC"), list("ab")).wrap(data)
        with TemporaryDirectory(dir=".") as dirname:
            outpath = os.path.join(dirname, "delme.tsv")
            darr.write(outpath)
            with open(outpath) as infile:
                contents = [l.strip().split() for l in infile]
            header = contents.pop(0)
            self.assertEqual(header, ["dim-1", "dim-2", "value"])
            got = {(k1, k2): int(v) for k1, k2, v in contents}
            self.assertEqual(got, darr.to_dict(flatten=True))

    def test_to_string(self):
        darr = DictArrayTemplate(2, 2).wrap(
            [[3.123456789, 2 * 3.123456789], [3 * 3.123456789, 4 * 3.123456789]]
        )
        self.assertEqual(
            darr.to_string(sep=","),
            "dim-1,dim-2,value\n0,0,3.123456789\n0,1,6.246913578\n1,0,9.370370367\n1,1,12.493827156",
        )
        self.assertEqual(
            darr.to_string(),
            "dim-1\tdim-2\tvalue\n0\t0\t3.123456789\n0\t1\t6.246913578\n1\t0\t9.370370367\n1\t1\t12.493827156",
        )
        self.assertEqual(
            darr.to_string(sep=" "),
            "dim-1 dim-2 value\n0 0 3.123456789\n0 1 6.246913578\n1 0 9.370370367\n1 1 12.493827156",
        )
        with self.assertRaises(ValueError):
            darr.to_string(format="md"),

    def test_to_table(self):
        """creates Table when ndim <= 2"""
        from cogent3.util.table import Table

        a1D = DictArrayTemplate(["a", "b"]).wrap([0, 1])
        t = a1D.to_table()
        self.assertIsInstance(t, Table)
        # 1D tables don't get an index_name column
        self.assertEqual(t.index_name, None)
        a2D = DictArrayTemplate(["a", "b"], ["c", "d"]).wrap(
            numpy.array([0, 1, 2, 3]).reshape((2, 2))
        )
        t = a2D.to_table()
        self.assertIsInstance(t, Table)
        self.assertEqual(t.shape, (2, 3))  # because index_name column added
        # make sure the 2D variant has an index_name column, name is empty string
        self.assertEqual(t.index_name, "")
        self.assertEqual(t.columns[""].tolist(), a2D.template.names[0])
        # which works
        self.assertEqual(t["b", "d"], 3)

        a3D = DictArrayTemplate(["a", "b"], ["c", "d"], ["e", "f"]).wrap(
            numpy.array([0, 1, 2, 3, 4, 5, 6, 7]).reshape((2, 2, 2))
        )
        with self.assertRaises(ValueError):
            _ = a3D.to_table()

    def test_interpret_index(self):
        """correctly handles just explicitly defined indices"""
        n = ["ab", "dna", "rna"]
        a1D = DictArrayTemplate(n)
        got = a1D.interpret_index(["ab", "rna"])
        self.assertEqual(got[0], ([0, 2],))
        got = a1D.interpret_index([0, 2])
        self.assertEqual(got[0], ([0, 2],))

    def test_slicing_combos(self):
        """different mixtures of slicing should work"""
        darr = DictArrayTemplate(list(DNA), list(DNA)).wrap(
            [
                [0.7, 0.1, 0.2, 0.3],
                [0.1, 0.7, 0.1, 0.3],
                [0.3, 0.2, 0.6, 0.3],
                [0.4, 0.1, 0.1, 0.7],
            ]
        )
        got = darr["C":"G", "C":"G"]
        assert_allclose(got.array, numpy.array([[0.7, 0.1], [0.2, 0.6]]))

        got = darr[[1, 2], [1, 2]]
        assert_allclose(got.array, numpy.array([[0.7, 0.1], [0.2, 0.6]]))

        got = darr[[2, 3], "C"]
        assert_allclose(got.array, numpy.array([0.2, 0.1]))
        got = darr["C", [2, 3]]
        assert_allclose(got.array, numpy.array([0.1, 0.3]))

        got = darr[[1, 2], "T":"A"]
        assert_allclose(got.array, numpy.array([[0.1, 0.7], [0.3, 0.2]]))

        got = darr["T":"A", [1, 2]]
        assert_allclose(got.array, numpy.array([[0.1, 0.2], [0.7, 0.1]]))

        # make sure we cope with keys that are int's
        nums = list(range(1, 5))
        darr = DictArrayTemplate(nums, nums).wrap(
            [
                [0.7, 0.1, 0.2, 0.3],
                [0.1, 0.7, 0.1, 0.3],
                [0.3, 0.2, 0.6, 0.3],
                [0.4, 0.1, 0.1, 0.7],
            ]
        )
        got = darr[[1, 2], [1, 2]]
        assert_allclose(got.array, numpy.array([[0.7, 0.1], [0.2, 0.6]]))

    def test_add(self):
        """can add compatible dict arrays"""
        data = numpy.array([[7, 1], [1, 7]])
        darr1 = DictArrayTemplate(list("AB"), list("CD")).wrap(data)
        darr2 = DictArrayTemplate(list("AB"), list("CD")).wrap(data)
        darr3 = darr1 + darr2
        assert_allclose(darr3.array, 2 * data)
        self.assertEqual(darr3.template.names, darr1.template.names)
        # must be correct type
        with self.assertRaises(TypeError):
            darr1 + data

        # must be equal dimensions
        with self.assertRaises(ValueError):
            darr1 + DictArrayTemplate(list("CD"), list("AB")).wrap(data)


if __name__ == "__main__":
    main()
