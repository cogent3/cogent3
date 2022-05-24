#!/usr/bin/env python

"""Unit tests for union_dict.
"""
from unittest import TestCase, main

from cogent3.util.union_dict import UnionDict


__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class UnionDictTests(TestCase):
    """Tests of individual functions in union_dict"""

    def test_attr(self):
        """test the "." read/write functionality"""
        d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
        self.assertEqual(d.a, 1)
        self.assertEqual(d.b, 2)
        self.assertEqual(d.d.e, 5)
        d.c = 0
        d.d.f = 0
        self.assertEqual(d.c, 0)
        self.assertEqual(d.d.f, 0)

    def test_construction(self):
        """should handle deeply nested dict"""
        data = {"width": 600.0, "xaxis": {"title": {"text": "Alignment Position"}}}
        d = UnionDict(data)
        self.assertEqual(d.xaxis.title.text, "Alignment Position")

    def test_construct_from_empty(self):
        """successfully define from an empty"""
        data = {"width": 600.0, "xaxis": {"title": {"text": "Alignment Position"}}}
        # empty object
        d = UnionDict()
        self.assertTrue(len(d) == 0)
        # using update
        d.update(data)
        self.assertEqual(d.xaxis.title.text, "Alignment Position")

    def test_construct_from_kwargs(self):
        """successfully define from an kwargs"""
        data = {"width": 600.0, "xaxis": {"title": {"text": "Alignment Position"}}}
        # empty object
        d = UnionDict(**data)
        self.assertEqual(d.xaxis.title.text, "Alignment Position")

    def test_union(self):
        """correctly adjust a prob vector so all values > minval"""
        d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
        e = UnionDict({"b": 0, "d": {"f": 0, "g": 7}})
        d |= e
        self.assertEqual(d.a, 1)
        self.assertEqual(d.b, 0)
        self.assertEqual(d.d.e, 5)
        self.assertEqual(d.d.f, 0)
        self.assertEqual(d.d.g, 7)

    def test_or(self):
        """should not modify original"""
        d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
        e = UnionDict({"b": 0, "d": {"f": 0, "g": 7}})
        f = d | e
        self.assertEqual(f.a, 1)
        self.assertEqual(f.b, 0)
        self.assertEqual(f.d.e, 5)
        self.assertEqual(f.d.f, 0)
        self.assertEqual(f.d.g, 7)
        self.assertTrue(f.d is not e.d)

    def test_union_value_dict(self):
        """replacing union or of a value with a dict should be dict"""
        d = UnionDict({"A": {"B": "Blah"}})
        e = UnionDict({"A": "Blah"})
        f = UnionDict(d.copy())
        f |= e
        self.assertNotEqual(d, f)
        e |= d
        self.assertEqual(d, e)

    def test_union_with_empty_sub_dict(self):
        """unioning with a dict that has an empty sub-dict"""
        d = UnionDict({"title": {}})
        e = UnionDict({"title": {"text": "Alignment Position"}})
        f = UnionDict(e.copy())
        e |= d
        self.assertEqual(e, f)

    def test_sub_dicts_are_union(self):
        """checks if UnionDict is propogated to children"""
        d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
        d.e = {"g": 7}
        d.e.g = {"h": 8}
        self.assertTrue(isinstance(d, UnionDict))
        self.assertTrue(isinstance(d.d, UnionDict))
        self.assertTrue(isinstance(d.e, UnionDict))
        self.assertTrue(isinstance(d.e.g, UnionDict))

    def test_get_subattr(self):
        """_getsubattr_ returns nested values via key"""
        d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
        self.assertEqual(d._getsubattr_([], "a"), 1)
        self.assertEqual(d._getsubattr_([], "d"), UnionDict({"e": 5, "f": 6}))
        self.assertEqual(d._getsubattr_(["d"], "e"), 5)

    def test_setitem(self):
        """should work via property or key"""
        d = UnionDict()
        d.a = 23
        d.b = dict(c=42)
        self.assertEqual(d.a, 23)
        self.assertEqual(d["a"], 23)
        self.assertEqual(d.b, dict(c=42))
        self.assertEqual(d.b.c, 42)
        self.assertEqual(d["b"], dict(c=42))
        self.assertIsInstance(d.b, UnionDict)


if __name__ == "__main__":
    main()
