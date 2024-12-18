#!/usr/bin/env python
"""Unit tests for parser support libraries dealing with records."""

from unittest import TestCase

from cogent3.parse.record import (
    DelimitedSplitter,
    FieldError,
    FieldMorpher,
    FieldWrapper,
    GenericRecord,
    Grouper,
    LineOrientedConstructor,
    MappedRecord,
    StrictFieldWrapper,
    TypeSetter,
    bool_setter,
    dict_adder,
    int_setter,
    list_adder,
    list_extender,
    raise_unknown_field,
    string_and_strip,
)


class recordsTests(TestCase):
    """Tests of top-level functionality in records."""

    def test_string_and_strip(self):
        """string_and_strip should convert all items to strings and strip them"""
        assert string_and_strip() == []
        assert string_and_strip("\t", " ", "\n\t") == ["", "", ""]
        assert string_and_strip("\ta\tb", 3, "   cde   e", None) == [
            "a\tb",
            "3",
            "cde   e",
            "None",
        ]

    def test_raise_unknown_field(self):
        """raise_unknown_field should always raise FieldError"""
        self.assertRaises(FieldError, raise_unknown_field, "xyz", 123)


class GrouperTests(TestCase):
    """Tests of the Grouper class."""

    def test_call(self):
        """Grouper should return lists containing correct number of groups"""
        empty = []
        s3 = "abc"
        s10 = list(range(10))
        g1 = Grouper(1)
        g2 = Grouper(2)
        g5 = Grouper(5)
        assert list(g1(empty)) == []
        assert list(g2(empty)) == []
        assert list(g5(empty)) == []
        assert list(g1(s3)) == [["a"], ["b"], ["c"]]
        assert list(g2(s3)) == [["a", "b"], ["c"]]
        assert list(g5(s3)) == [["a", "b", "c"]]
        assert list(g1(s10)) == [[i] for i in range(10)]
        assert list(g2(s10)) == [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]]
        assert list(g5(s10)) == [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]

    def test_call_bad(self):
        """Grouper call should raise ValueError if NumItems is not an int"""
        g_none = Grouper(None)
        g_neg = Grouper(-1)
        g_zero = Grouper(0)
        g_alpha = Grouper("abc")
        for g in (g_none, g_neg, g_zero, g_alpha):
            iterator = g("abcd")
            self.assertRaises(ValueError, list, iterator)


class DelimitedSplitterTests(TestCase):
    """Tests of the DelimitedSplitter factory function."""

    def test_parsers(self):
        """DelimitedSplitter should return function with correct behavior"""
        empty = DelimitedSplitter()
        space = DelimitedSplitter(None)
        semicolon = DelimitedSplitter(";")
        twosplits = DelimitedSplitter(";", 2)
        allsplits = DelimitedSplitter(";", None)
        lastone = DelimitedSplitter(";", -1)
        lasttwo = DelimitedSplitter(";", -2)

        assert empty("a   b  c") == ["a", "b  c"]
        assert empty("abc") == ["abc"]
        assert empty("   ") == []

        assert empty("a  b  c") == space("a  b  c")
        assert semicolon("  a  ; b   ;  c  d") == ["a", "b   ;  c  d"]
        assert twosplits("  a  ; b   ;  c  d") == ["a", "b", "c  d"]
        assert allsplits(" a ;  b  ; c;;d;e  ;") == ["a", "b", "c", "", "d", "e", ""]
        assert lastone(" a ;  b  ; c;;d;e  ;") == ["a ;  b  ; c;;d;e", ""]
        assert lasttwo(" a ;  b  ; c;;d;e  ;") == ["a ;  b  ; c;;d", "e", ""]
        assert lasttwo("") == []
        assert lasttwo("x") == ["x"]
        assert lasttwo("x;") == ["x", ""]


class GenericRecordTests(TestCase):
    """Tests of the GenericRecord class"""

    class gr(GenericRecord):
        Required = {"a": "x", "b": [], "c": {}}

    def test_init(self):
        """GenericRecord init should work OK empty or with data"""
        assert GenericRecord() == {}
        assert GenericRecord({"a": 1}) == {"a": 1}
        assert isinstance(GenericRecord(), GenericRecord)

    def test_init_subclass(self):
        """GenericRecord subclass init should include required data"""
        assert self.gr() == {"a": "x", "b": [], "c": {}}
        assert self.gr({"a": []}) == {"a": [], "b": [], "c": {}}
        assert isinstance(self.gr(), self.gr)
        assert isinstance(self.gr(), GenericRecord)

    def test_delitem(self):
        """GenericRecord delitem should fail if item required"""
        g = self.gr()
        g["d"] = 3
        assert g == {"a": "x", "b": [], "c": {}, "d": 3}
        del g["d"]
        assert g == {"a": "x", "b": [], "c": {}}
        self.assertRaises(AttributeError, g.__delitem__, "a")
        g["c"][3] = 4
        assert g["c"] == {3: 4}

    def test_copy(self):
        """GenericRecord copy should include attributes and set correct class"""
        g = self.gr()
        g["a"] = "abc"
        g.X = "y"
        h = g.copy()
        assert g == h
        assert isinstance(h, self.gr)
        assert h.X == "y"
        assert h == {"a": "abc", "b": [], "c": {}}


class MappedRecordTests(TestCase):
    """Tests of the MappedRecord class"""

    def setUp(self):
        """Define a few standard MappedRecords"""
        self.empty = MappedRecord()
        self.single = MappedRecord({"a": 3})
        self.several = MappedRecord(a=4, b=5, c="a", d=[1, 2, 3])

    def test_init_empty(self):
        """MappedRecord empty init should work OK"""
        g = MappedRecord()
        assert g == {}

    def test_init_data(self):
        """MappedRecord should work like normal dict init"""
        exp = {"a": 3, "b": 4}
        assert MappedRecord({"a": 3, "b": 4}) == exp
        assert MappedRecord(a=3, b=4) == exp
        assert MappedRecord([["a", 3], ["b", 4]]) == exp

    def test_init_subclass(self):
        """MappedRecord subclasses should behave as expected"""

        class rec(MappedRecord):
            Required = {"a": {}, "b": "xyz", "c": 3}
            Aliases = {"B": "b"}

        r = rec()
        assert r == {"a": {}, "b": "xyz", "c": 3}
        # test that subclassing is correct
        s = r.copy()
        assert isinstance(s, rec)
        # test Aliases
        s.B = 0
        assert s == {"a": {}, "b": 0, "c": 3}
        # test Required
        try:
            del s.B
        except AttributeError:
            pass
        else:
            msg = "Subclass failed to catch requirement"
            raise AssertionError(msg)

    def test_getattr(self):
        """MappedRecord getattr should look in dict after real attrs"""
        s = self.several
        assert s.Aliases == {}
        assert s.a == 4
        assert s.d == [1, 2, 3]
        for key in s:
            assert getattr(s, key) == s[key]
        assert "xyz" not in s
        assert s.xyz is None
        assert s["xyz"] is None
        s.Aliases = {"xyz": "a"}
        assert s["xyz"] == 4

    def test_setattr(self):
        """MappedRecord setattr should add to dict"""
        s = self.single
        # check that we haven't screwed up normal attribute setting
        assert "Aliases" not in s
        s.Aliases = {"x": "y"}
        assert "Aliases" not in s
        assert s.Aliases == {"x": "y"}
        s.x = 5
        assert "x" in s
        assert s["x"] == 5
        assert s.x == 5
        s.Aliases = {"XYZ": "b"}
        s.XYZ = 3
        assert s.b == 3

    def test_delattr(self):
        """MappedRecord delattr should work for 'normal' and other attributes"""
        s = self.single
        s.__dict__["x"] = "y"
        assert "x" not in s
        assert s.x == "y"
        del s.x
        assert s.x is None
        assert s == {"a": 3}
        # try it for an internal attribute: check it doesn't delete anything
        # else
        s.b = 4
        assert s == {"a": 3, "b": 4}
        del s.a
        assert s == {"b": 4}
        del s.abc
        assert s == {"b": 4}
        s.Required = {"b": True}
        try:
            del s.b
        except AttributeError:
            pass
        else:
            msg = "Allowed deletion of required attribute"
            raise AssertionError(msg)
        s.a = 3
        assert s.a == 3
        s.Aliases = {"xyz": "a"}
        del s.xyz
        assert s.a is None

    def test_getitem(self):
        """MappedRecord getitem should work only for keys, not attributes"""
        s = self.single
        assert s["Required"] is None
        assert s["a"] == 3
        assert s["xyz"] is None
        assert s[list("abc")] is None
        s.Aliases = {"xyz": "a"}
        assert s["xyz"] == 3

    def test_setitem(self):
        """MappedRecord setitem should work only for keys, not attributes"""
        s = self.single
        s["Required"] = None
        assert s == {"a": 3, "Required": None}
        assert s.Required == {}
        assert s.Required is not None
        s["c"] = 5
        assert s == {"a": 3, "c": 5, "Required": None}
        # still not allowed unhashable objects as keys
        self.assertRaises(TypeError, s.__setitem__, list(range(3)))
        s.Aliases = {"C": "c"}
        s["C"] = 3
        assert s == {"a": 3, "c": 3, "Required": None}

    def test_delitem(self):
        """MappedRecord delitem should only work for keys, not attributes"""
        s = self.single
        del s["Required"]
        assert s.Required == {}
        s.Required = {"a": True}
        try:
            del s["a"]
        except AttributeError:
            pass
        else:
            msg = "Allowed deletion of required item"
            raise AssertionError(msg)
        s.Aliases = {"B": "b"}
        s.b = 5
        assert s.b == 5
        del s.B
        assert s.b is None

    def test_contains(self):
        """MappedRecord contains should use aliases, but not apply to attrs"""
        s = self.single
        assert "a" in s
        assert "b" not in s
        s.b = 5
        assert "b" in s
        assert "Required" not in s
        assert "A" not in s
        s.Aliases = {"A": "a"}
        assert "A" in s

    def test_get(self):
        """MappedRecord get should be typesafe against unhashables"""
        s = self.single
        assert s.get(1, 6) == 6
        assert s.get("a", "xyz") == 3
        assert s.get("ABC", "xyz") == "xyz"
        s.Aliases = {"ABC": "a"}
        assert s.get("ABC", "xyz") == 3
        assert s.get([1, 2, 3], "x") == "x"

    def test_setdefault(self):
        """MappedRecord setdefault should not be typesafe against unhashables"""
        s = self.single
        x = s.setdefault("X", "xyz")
        assert x == "xyz"
        assert s == {"a": 3, "X": "xyz"}
        self.assertRaises(TypeError, s.setdefault, ["a", "b"], "xyz")

    def test_update(self):
        """MappedRecord update should transparently convert keys"""
        s = self.single
        s.b = 999
        s.Aliases = {"XYZ": "x", "ABC": "a"}
        d = {"ABC": 111, "CVB": 222}
        s.update(d)
        assert s == {"a": 111, "b": 999, "CVB": 222}

    def test_copy(self):
        """MappedRecord copy should return correct class"""
        s = self.single
        t = s.copy()
        assert isinstance(t, MappedRecord)
        s.Aliases = {"XYZ": "x"}
        u = s.copy()
        u.Aliases["ABC"] = "a"
        assert s.Aliases == {"XYZ": "x"}
        assert t.Aliases == {}
        assert u.Aliases == {"XYZ": "x", "ABC": "a"}

    def test_subclass(self):
        """MappedRecord subclassing should work correctly"""

        class ret3(MappedRecord):
            DefaultValue = 3
            ClassData = "xyz"

        x = ret3({"ABC": 777, "DEF": "999"})
        assert x.ZZZ == 3
        assert x.ABC == 777
        assert x.DEF == "999"
        assert x.ClassData == "xyz"
        x.ZZZ = 6
        assert x.ZZZ == 6
        assert x.ZZ == 3
        x.ClassData = "qwe"
        assert x.ClassData == "qwe"
        assert ret3.ClassData == "xyz"

    def test_DefaultValue(self):
        """MappedRecord DefaultValue should give new copy when requested"""

        class m(MappedRecord):
            DefaultValue = []

        a = m()
        b = m()
        assert a["abc"] is not b["abc"]
        assert a["abc"] == b["abc"]


class dummy:
    """Do-nothing class whose attributes can be freely abused."""


class TypeSetterTests(TestCase):
    """Tests of the TypeSetter class"""

    def test_setter_empty(self):
        """TypeSetter should set attrs to vals on empty init"""
        d = dummy()
        ident = TypeSetter()
        ident(d, "x", "abc")
        assert d.x == "abc"
        ident(d, "y", 3)
        assert d.y == 3
        ident(d, "x", 2)
        assert d.x == 2

    def test_setter_typed(self):
        """TypeSetter should set attrs to constructor(val) when specified"""
        d = dummy()
        i = TypeSetter(int)
        i(d, "zz", 3)
        assert d.zz == 3
        i(d, "xx", "456")
        assert d.xx == 456


class TypeSetterLikeTests(TestCase):
    """Tests of the functions that behave similarly to TypeSetter products"""

    def test_list_adder(self):
        """list_adder should add items to list, creating if necessary"""
        d = dummy()
        list_adder(d, "x", 3)
        assert d.x == [3]
        list_adder(d, "x", "abc")
        assert d.x == [3, "abc"]
        list_adder(d, "y", [2, 3])
        assert d.x == [3, "abc"]
        assert d.y == [[2, 3]]

    def test_list_extender(self):
        """list_adder should add items to list, creating if necessary"""
        d = dummy()
        list_extender(d, "x", "345")
        assert d.x == ["3", "4", "5"]
        list_extender(d, "x", "abc")
        assert d.x == ["3", "4", "5", "a", "b", "c"]
        list_extender(d, "y", [2, 3])
        assert d.x == ["3", "4", "5", "a", "b", "c"]
        assert d.y == [2, 3]
        list_extender(d, "y", None)
        assert d.y == [2, 3, None]

    def test_dict_adder(self):
        """dict_adder should add items to dict, creating if necessary"""
        d = dummy()
        dict_adder(d, "x", 3)
        assert d.x == {3: None}
        dict_adder(d, "x", "ab")
        assert d.x == {3: None, "a": "b"}
        dict_adder(d, "x", ["a", 0])
        assert d.x == {3: None, "a": 0}
        dict_adder(d, "y", None)
        assert d.x == {3: None, "a": 0}
        assert d.y == {None: None}


class LineOrientedConstructorTests(TestCase):
    """Tests of the LineOrientedConstructor class"""

    def test_init_empty(self):
        """LOC empty init should succeed with expected defaults"""
        l = LineOrientedConstructor()
        assert l.Lines == []
        assert l.LabelSplitter(" ab  cd  ") == ["ab", "cd"]
        assert l.FieldMap == {}
        assert l.Constructor == MappedRecord
        assert l.Strict is False

    def test_empty_LOC(self):
        """LOC empty should fail if strict, fill fields if not strict"""
        data = ["abc   def", "3  n", "\t  abc   \txyz\n\n", "fgh   "]
        l = LineOrientedConstructor()
        result = l()
        assert result == {}
        result = l([])
        assert result == {}
        result = l(["   ", "\n\t   "])
        assert result == {}
        result = l(data)
        assert result == {"abc": "xyz", "3": "n", "fgh": None}

    def test_full_LOC(self):
        """LOC should behave as expected when initialized with rich data"""
        data = [
            "abc\t def",
            " 3 \t n",
            "  abc   \txyz\n\n",
            "x\t5",
            "fgh   ",
            "x\t3    ",
        ]

        class rec(MappedRecord):
            Required = {"abc": []}

        maps = {"abc": list_adder, "x": int_setter, "fgh": bool_setter}
        label_splitter = DelimitedSplitter("\t")
        constructor = rec
        strict = True
        loc_bad = LineOrientedConstructor(
            data,
            label_splitter,
            maps,
            constructor,
            strict,
        )
        self.assertRaises(FieldError, loc_bad)
        strict = False
        loc_good = LineOrientedConstructor(
            data,
            label_splitter,
            maps,
            constructor,
            strict,
        )
        result = loc_good()
        assert isinstance(result, rec)
        assert result == {"abc": ["def", "xyz"], "3": "n", "fgh": False, "x": 3}


class fake_dict(dict):
    """Test that constructors return the correct subclass"""


class FieldWrapperTests(TestCase):
    """Tests of the FieldWrapper factory function"""

    def test_default(self):
        """Default FieldWrapper should wrap fields and labels"""
        fields = list("abcde")
        f = FieldWrapper(fields)
        assert f("") == {}
        assert f("xy za ") == {"a": "xy", "b": "za"}
        assert f("1   2\t\t 3  \n4 5 6") == {
            "a": "1",
            "b": "2",
            "c": "3",
            "d": "4",
            "e": "5",
        }

    def test_splitter(self):
        """FieldWrapper with splitter should use that splitter"""
        fields = ["label", "count"]
        splitter = DelimitedSplitter(":", -1)
        f = FieldWrapper(fields, splitter)
        assert f("") == {}
        assert f("nknasd:") == {"label": "nknasd", "count": ""}
        assert f("n:k:n:a:sd  ") == {"label": "n:k:n:a", "count": "sd"}

    def test_constructor(self):
        """FieldWrapper with constructor should use that constructor"""
        fields = list("abc")
        f = FieldWrapper(fields, constructor=fake_dict)
        assert f("x y") == {"a": "x", "b": "y"}
        assert isinstance(f("x y"), fake_dict)


class StrictFieldWrapperTests(TestCase):
    """Tests of the StrictFieldWrapper factory function"""

    def test_default(self):
        """Default StrictFieldWrapper should wrap fields if count correct"""
        fields = list("abcde")
        f = StrictFieldWrapper(fields)
        assert f("1   2\t\t 3  \n4 5 ") == {
            "a": "1",
            "b": "2",
            "c": "3",
            "d": "4",
            "e": "5",
        }
        self.assertRaises(FieldError, f, "")
        self.assertRaises(FieldError, f, "xy za ")

    def test_splitter(self):
        """StrictFieldWrapper with splitter should use that splitter"""
        fields = ["label", "count"]
        splitter = DelimitedSplitter(":", -1)
        f = StrictFieldWrapper(fields, splitter)
        assert f("n:k:n:a:sd  ") == {"label": "n:k:n:a", "count": "sd"}
        assert f("nknasd:") == {"label": "nknasd", "count": ""}
        self.assertRaises(FieldError, f, "")

    def test_constructor(self):
        """StrictFieldWrapper with constructor should use that constructor"""
        fields = list("ab")
        f = StrictFieldWrapper(fields, constructor=fake_dict)
        assert f("x y") == {"a": "x", "b": "y"}
        assert isinstance(f("x y"), fake_dict)


class FieldMorpherTests(TestCase):
    """Tests of the FieldMorpher class."""

    def test_default(self):
        """FieldMorpher default should use correct constructors"""
        fm = FieldMorpher({"a": int, "b": str})
        assert fm({"a": "3", "b": 456}) == {"a": 3, "b": "456"}

    def test_default_error(self):
        """FieldMorpher default should raise FieldError on unknown fields"""
        fm = FieldMorpher({"a": int, "b": str})
        self.assertRaises(FieldError, fm, {"a": "3", "b": 456, "c": "4"})

    def test_altered_default(self):
        """FieldMorpher with default set should apply it"""

        def func(x, y):
            return str(x), float(y) - 0.5

        fm = FieldMorpher({"3": str, 4: int}, func)
        # check that recognized values aren't tampered with
        assert fm({3: 3, 4: "4"}) == {"3": "3", 4: 4}
        # check that unrecognized values get the appropriate conversion
        assert fm({3: 3, 5: "5"}) == {"3": "3", "5": 4.5}
