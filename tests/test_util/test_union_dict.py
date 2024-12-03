"""Unit tests for union_dict."""

from cogent3.util.union_dict import UnionDict


def test_attr():
    """test the "." read/write functionality"""
    d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
    assert d.a == 1
    assert d.b == 2
    assert d.d.e == 5
    d.c = 0
    d.d.f = 0
    assert d.c == 0
    assert d.d.f == 0


def test_construction():
    """should handle deeply nested dict"""
    data = {"width": 600.0, "xaxis": {"title": {"text": "Alignment Position"}}}
    d = UnionDict(data)
    assert d.xaxis.title.text == "Alignment Position"


def test_construct_from_empty():
    """successfully define from an empty"""
    data = {"width": 600.0, "xaxis": {"title": {"text": "Alignment Position"}}}
    # empty object
    d = UnionDict()
    assert len(d) == 0
    # using update
    d.update(data)
    assert d.xaxis.title.text == "Alignment Position"


def test_construct_from_kwargs():
    """successfully define from an kwargs"""
    data = {"width": 600.0, "xaxis": {"title": {"text": "Alignment Position"}}}
    # empty object
    d = UnionDict(**data)
    assert d.xaxis.title.text == "Alignment Position"


def test_union():
    """correctly merges two UnionDicts"""
    d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
    e = UnionDict({"b": 0, "d": {"f": 0, "g": 7}})
    d |= e
    assert d.a == 1
    assert d.b == 0
    assert d.d.e == 5
    assert d.d.f == 0
    assert d.d.g == 7


def test_or():
    """should not modify original"""
    d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
    e = UnionDict({"b": 0, "d": {"f": 0, "g": 7}})
    f = d | e
    assert f.a == 1
    assert f.b == 0
    assert f.d.e == 5
    assert f.d.f == 0
    assert f.d.g == 7
    assert f.d is not e.d


def test_union_value_dict():
    """replacing union or of a value with a dict should be dict"""
    d = UnionDict({"A": {"B": "Blah"}})
    e = UnionDict({"A": "Blah"})
    f = UnionDict(d.copy())
    f |= e
    assert d != f
    e |= d
    assert d == e


def test_union_with_empty_sub_dict():
    """unioning with a dict that has an empty sub-dict"""
    d = UnionDict({"title": {}})
    e = UnionDict({"title": {"text": "Alignment Position"}})
    f = UnionDict(e.copy())
    e |= d
    assert e == f


def test_sub_dicts_are_union():
    """checks if UnionDict is propogated to children"""
    d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
    d.e = {"g": 7}
    d.e.g = {"h": 8}
    assert isinstance(d, UnionDict)
    assert isinstance(d.d, UnionDict)
    assert isinstance(d.e, UnionDict)
    assert isinstance(d.e.g, UnionDict)


def test_get_subattr():
    """_getsubattr_ returns nested values via key"""
    d = UnionDict({"a": 1, "b": 2, "c": 3, "d": {"e": 5, "f": 6}})
    assert d._getsubattr_([], "a") == 1
    assert d._getsubattr_([], "d") == UnionDict({"e": 5, "f": 6})
    assert d._getsubattr_(["d"], "e") == 5


def test_setitem():
    """should work via property or key"""
    d = UnionDict()
    d.a = 23
    d.b = {"c": 42}
    assert d.a == 23
    assert d["a"] == 23
    assert d.b == {"c": 42}
    assert d.b.c == 42
    assert d["b"] == {"c": 42}
    assert isinstance(d.b, UnionDict)
