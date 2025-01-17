#!/usr/bin/env python
"""Unit tests for Info class and associated objects (DbRef, DbRefs, etc.)."""

import warnings
from unittest import TestCase

from cogent3.core.info import DbRef, DbRefs, Info, _make_list


class DbRefTests(TestCase):
    """Tests of the DbRef object."""

    def setUp(self):
        """Define a standard DbRef object"""
        self.data = {
            "Accession": "xyz",
            "Db": "abc",
            "name": "qwe",
            "Description": "blah",
            "Data": list(range(20)),
        }
        self.db = DbRef(**self.data)

    def test_init_minimal(self):
        """DbRef minimal init should fill fields as expected"""
        d = DbRef("abc")
        assert d.Accession == "abc"
        assert d.Db == ""
        assert d.name == ""
        assert d.Description == ""
        assert d.Data is None
        # empty init not allowed
        self.assertRaises(TypeError, DbRef)

    def test_init(self):
        """DbRef init should insert correct data"""
        for attr, val in list(self.data.items()):
            assert getattr(self.db, attr) == val

    def test_str(self):
        """DbRef str should be the same as the accession str"""
        assert str(self.db) == "xyz"
        self.db.Accession = 12345
        assert str(self.db) == "12345"

    def test_int(self):
        """DbRef int should be the same as the accession int"""
        self.assertRaises(ValueError, int, self.db)
        self.db.Accession = "12345"
        assert int(self.db) == 12345

    def test_cmp(self):
        """DbRef cmp should first try numeric, then alphabetic, cmp."""
        assert DbRef("abc") < DbRef("xyz")
        assert DbRef("abc") == DbRef("abc")
        assert DbRef("123") > DbRef("14")
        assert DbRef("123") < DbRef("abc")
        # check that it ignores other attributes
        assert DbRef("x", "y", "z", "a", "b") == DbRef("x")


class infoTests(TestCase):
    """Tests of top-level functions."""

    def test_make_list(self):
        """_make_list should always return a list"""
        assert _make_list("abc") == ["abc"]
        assert _make_list([]) == []
        assert _make_list(None) == [None]
        assert _make_list({"x": "y"}) == [{"x": "y"}]
        assert _make_list([1, 2, 3]) == [1, 2, 3]


class DbRefsTests(TestCase):
    """Tests of the DbRefs class."""

    def test_init_empty(self):
        """DbRefs empty init should work as expected"""
        assert DbRefs() == {}

    def test_init_data(self):
        """DbRefs init with data should produce expected results"""
        d = DbRefs({"GenBank": "ab", "GO": (3, 44), "PDB": ["asdf", "ghjk"]})
        assert d == {"GenBank": ["ab"], "GO": [3, 44], "PDB": ["asdf", "ghjk"]}
        d.GenBank = "xyz"
        assert d["GenBank"] == ["xyz"]


class InfoTests(TestCase):
    """Tests of the Info class."""

    def test_init_empty(self):
        """Info empty init should work as expected"""
        d = Info()
        assert len(d) == 1
        assert "Refs" in d
        assert d.Refs == DbRefs()
        assert isinstance(d.Refs, DbRefs)

    def test_init_data(self):
        """Info init with data should put items in correct places"""
        # need to check init, setting, and resetting of attributes that belong
        # in the Info object and attributes that belong in Info.Refs. Also need
        # to check __getitem__, __setitem__, and __contains__.
        d = Info({"x": 3, "GO": 12345})
        assert d.x == 3
        assert d.GO == [12345]
        assert d.Refs.GO == [12345]
        try:
            del d.Refs
        except AttributeError:
            pass
        else:
            msg = "Failed to prevent deletion of required key Refs"
            raise Exception(msg)
        d.GenBank = ("qaz", "wsx")
        assert d.GenBank == ["qaz", "wsx"]
        assert "GenBank" in d.Refs
        assert "GenBank" in d
        d.GenBank = "xyz"
        assert d.GenBank == ["xyz"]
        assert d.GenBank is d.Refs.GenBank
        d.GO = "x"
        assert d.GO == ["x"]
        d.GO.append("y")
        assert d.GO == ["x", "y"]
        d.ZZZ = "zzz"
        assert d.ZZZ == "zzz"
        assert "ZZZ" not in d.Refs
        assert "XXX" not in d
        assert d.XXX is None

    def test_identity(self):
        """Info should get its own new Refs when created"""
        i = Info()
        j = Info()
        assert i is not j
        assert i.Refs is not j.Refs

    def test_update(self):
        """update should warn the user of overlapping keys"""
        with warnings.catch_warnings(record=True) as w:
            d1 = Info({"key1": "value1", "key2": "value2", "key3": "value3"})
            d2 = Info({"key2": "value2", "key3": "value3", "key4": "value4"})
            d1.update(d2)
            assert len(w) == 1


# run the following if invoked from command-line
