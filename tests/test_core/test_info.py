#!/usr/bin/env python
"""Unit tests for Info class and associated objects (DbRef, DbRefs, etc.).
"""
import warnings

from unittest import TestCase, main

from cogent3.core.info import DbRef, DbRefs, Info, _make_list


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class DbRefTests(TestCase):
    """Tests of the DbRef object."""

    def setUp(self):
        """Define a standard DbRef object"""
        self.data = dict(
            Accession="xyz",
            Db="abc",
            name="qwe",
            Description="blah",
            Data=list(range(20)),
        )
        self.db = DbRef(**self.data)

    def test_init_minimal(self):
        """DbRef minimal init should fill fields as expected"""
        d = DbRef("abc")
        self.assertEqual(d.Accession, "abc")
        self.assertEqual(d.Db, "")
        self.assertEqual(d.name, "")
        self.assertEqual(d.Description, "")
        self.assertEqual(d.Data, None)
        # empty init not allowed
        self.assertRaises(TypeError, DbRef)

    def test_init(self):
        """DbRef init should insert correct data"""
        for attr, val in list(self.data.items()):
            self.assertEqual(getattr(self.db, attr), val)

    def test_str(self):
        """DbRef str should be the same as the accession str"""
        self.assertEqual(str(self.db), "xyz")
        self.db.Accession = 12345
        self.assertEqual(str(self.db), "12345")

    def test_int(self):
        """DbRef int should be the same as the accession int"""
        self.assertRaises(ValueError, int, self.db)
        self.db.Accession = "12345"
        self.assertEqual(int(self.db), 12345)

    def test_cmp(self):
        """DbRef cmp should first try numeric, then alphabetic, cmp."""
        self.assertLess(DbRef("abc"), DbRef("xyz"))
        self.assertEqual(DbRef("abc"), DbRef("abc"))
        self.assertGreater(DbRef("123"), DbRef("14"))
        self.assertLess(DbRef("123"), DbRef("abc"))
        # check that it ignores other attributes
        self.assertEqual(DbRef("x", "y", "z", "a", "b"), DbRef("x"))


class infoTests(TestCase):
    """Tests of top-level functions."""

    def test_make_list(self):
        """_make_list should always return a list"""
        self.assertEqual(_make_list("abc"), ["abc"])
        self.assertEqual(_make_list([]), [])
        self.assertEqual(_make_list(None), [None])
        self.assertEqual(_make_list({"x": "y"}), [{"x": "y"}])
        self.assertEqual(_make_list([1, 2, 3]), [1, 2, 3])


class DbRefsTests(TestCase):
    """Tests of the DbRefs class."""

    def test_init_empty(self):
        """DbRefs empty init should work as expected"""
        self.assertEqual(DbRefs(), {})

    def test_init_data(self):
        """DbRefs init with data should produce expected results"""
        d = DbRefs({"GenBank": "ab", "GO": (3, 44), "PDB": ["asdf", "ghjk"]})
        self.assertEqual(d, {"GenBank": ["ab"], "GO": [3, 44], "PDB": ["asdf", "ghjk"]})
        d.GenBank = "xyz"
        self.assertEqual(d["GenBank"], ["xyz"])


class InfoTests(TestCase):
    """Tests of the Info class."""

    def test_init_empty(self):
        """Info empty init should work as expected"""
        d = Info()
        self.assertEqual(len(d), 1)
        self.assertIn("Refs", d)
        self.assertEqual(d.Refs, DbRefs())
        self.assertTrue(isinstance(d.Refs, DbRefs))

    def test_init_data(self):
        """Info init with data should put items in correct places"""
        # need to check init, setting, and resetting of attributes that belong
        # in the Info object and attributes that belong in Info.Refs. Also need
        # to check __getitem__, __setitem__, and __contains__.
        d = Info({"x": 3, "GO": 12345})
        self.assertEqual(d.x, 3)
        self.assertEqual(d.GO, [12345])
        self.assertEqual(d.Refs.GO, [12345])
        try:
            del d.Refs
        except AttributeError:
            pass
        else:
            raise Exception("Failed to prevent deletion of required key Refs" "")
        d.GenBank = ("qaz", "wsx")
        self.assertEqual(d.GenBank, ["qaz", "wsx"])
        self.assertIn("GenBank", d.Refs)
        self.assertIn("GenBank", d)
        d.GenBank = "xyz"
        self.assertEqual(d.GenBank, ["xyz"])
        self.assertIs(d.GenBank, d.Refs.GenBank)
        d.GO = "x"
        self.assertEqual(d.GO, ["x"])
        d.GO.append("y")
        self.assertEqual(d.GO, ["x", "y"])
        d.ZZZ = "zzz"
        self.assertEqual(d.ZZZ, "zzz")
        self.assertNotIn("ZZZ", d.Refs)
        self.assertNotIn("XXX", d)
        self.assertEqual(d.XXX, None)

    def test_identity(self):
        """Info should get its own new Refs when created"""
        i = Info()
        j = Info()
        self.assertIsNot(i, j)
        self.assertIsNot(i.Refs, j.Refs)

    def test_update(self):
        """update should warn the user of overlapping keys"""
        with warnings.catch_warnings(record=True) as w:
            d1 = Info({"key1": "value1", "key2": "value2", "key3": "value3"})
            d2 = Info({"key2": "value2", "key3": "value3", "key4": "value4"})
            d1.update(d2)
            self.assertEqual(len(w), 1)


# run the following if invoked from command-line
if __name__ == "__main__":
    main()
