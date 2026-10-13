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


def test_or_with_plain_dict():
    d = UnionDict({"a": 1, "d": {"e": 5}})
    f = d | {"d": {"g": {"h": 8}}}
    assert isinstance(f, UnionDict)
    assert f.d.e == 5
    assert f.d.g.h == 8
    assert "g" not in d.d


def test_two_level_shared_merge():
    d = UnionDict({"d": {"e": {"k": 1}}})
    f = d | {"d": {"e": {"m": 2}}}
    assert f.d.e == {"k": 1, "m": 2}
    assert d.d.e == {"k": 1}


def test_or_leaves_right_operand_unaliased():
    d = UnionDict({"a": 1})
    e = UnionDict({"b": {"c": 1}, "p": {"q": {"r": 1}}})
    f = d | e
    f.b.c = 99
    f.p.q.r = 99
    assert e.b.c == 1
    assert e.p.q.r == 1


def test_union_none_replaces_dict():
    d = UnionDict({"a": {"b": 1}})
    d |= {"a": None}
    assert d.a is None


def test_union_dict_replaces_nested_value():
    d = UnionDict({"a": {"b": 1}})
    d |= {"a": {"b": {"c": 2}}}
    assert isinstance(d.a.b, UnionDict)
    assert d.a.b.c == 2


def test_construction_leaves_input_unchanged():
    data = {"a": {"b": {"c": 1}}}
    d = UnionDict(data)
    d.a.b.c = 2
    assert type(data["a"]) is dict
    assert type(data["a"]["b"]) is dict
    assert data["a"]["b"]["c"] == 1


def test_nested_ior_keeps_held_references():
    layout = UnionDict({"xaxis": {"title": {"text": "a"}}})
    xaxis = layout.xaxis
    title = layout.xaxis.title
    layout.xaxis |= {"range": [0, 1]}
    title.text = "z"
    assert layout.xaxis is xaxis
    assert layout.xaxis.title.text == "z"


def test_assigning_union_dict_stores_it():
    d = UnionDict()
    inner = UnionDict({"c": 1})
    d.b = inner
    assert d.b is inner


def test_construction_preserves_shared_children():
    shared = {"c": 1}
    d = UnionDict({"a": shared, "b": shared})
    assert d.a is d.b


def test_construction_from_self_referencing_dict():
    data: dict[str, object] = {"x": 1}
    data["me"] = data
    d = UnionDict(data)
    assert d.me is d
    assert d.me.x == 1


def test_union_with_self_referencing_dicts():
    data: dict[str, object] = {"x": 1}
    data["me"] = data
    d = UnionDict(data)
    d |= data
    assert d.x == 1
    assert d.me.me.x == 1


def test_update_kwargs_converts_children():
    d = UnionDict()
    d.update(b={"c": {"e": 2}})
    assert d.b.c.e == 2


def test_ior_with_plain_dict_converts_children():
    d = UnionDict({"a": 1})
    d |= {"b": {"c": {"e": 2}}}
    assert isinstance(d.b, UnionDict)
    assert isinstance(d.b.c, UnionDict)
    assert d.b.c.e == 2


def test_union_method():
    d = UnionDict({"d": {"e": 5, "f": 6}})
    d.union({"d": {"f": 0}})
    assert d.d == {"e": 5, "f": 0}


def _make_nested(depth: int, leaf: dict[str, object]) -> dict[str, object]:
    data = leaf
    for _ in range(depth):
        data = {"x": data}
    return data


def test_deeply_nested_construction():
    depth = 2000
    d = UnionDict(_make_nested(depth, {"leaf": 1}))
    node = d
    for _ in range(depth):
        node = node.x
    assert isinstance(node, UnionDict)
    assert node.leaf == 1


def test_deeply_nested_union():
    depth = 2000
    d = UnionDict(_make_nested(depth, {"a": 1}))
    d |= _make_nested(depth, {"b": 2})
    node = d
    for _ in range(depth):
        node = node.x
    assert node == {"a": 1, "b": 2}
