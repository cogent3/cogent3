from types import MappingProxyType

import pytest

from cogent3.core.metadata import MetaData


def make_md(stored=None, local=None):
    stored = MappingProxyType(stored or {})
    return MetaData(stored=stored, local=local)


def test_get_set_contains():
    md = make_md({"a": 1}, {"b": 2})
    md["c"] = 3
    assert md["a"] == 1
    assert md["b"] == 2
    assert md["c"] == 3
    assert "a" in md
    assert "c" in md
    assert "z" not in md
    assert md.get("z") is None
    assert md.get("z", 5) == 5
    assert len(md) == 3


def test_missing_key_raises():
    md = make_md({"a": 1})
    with pytest.raises(KeyError):
        md["z"]
    with pytest.raises(KeyError):
        del md["z"]


def test_local_wins():
    stored = {"a": 1}
    md = make_md(stored)
    md["a"] = 0
    assert md["a"] == 0
    assert md.to_dict() == {"a": 0}
    assert stored == {"a": 1}


def test_local_none_is_a_value():
    md = make_md({"a": 1})
    md["a"] = None
    assert md["a"] is None
    assert md.to_dict() == {"a": None}


def test_delete_stored_hides():
    stored = MappingProxyType({"a": 1})
    md = MetaData(stored=stored)
    del md["a"]
    assert "a" not in md
    assert md.to_dict() == {}
    assert len(md) == 0
    assert list(md) == []
    assert dict(stored) == {"a": 1}
    with pytest.raises(KeyError):
        md["a"]
    with pytest.raises(KeyError):
        del md["a"]


def test_delete_hidden_key_then_set():
    md = make_md({"a": 1})
    del md["a"]
    md["a"] = 2
    assert md.to_dict() == {"a": 2}


def test_delete_local():
    md = make_md(local={"b": 2})
    del md["b"]
    assert "b" not in md
    with pytest.raises(KeyError):
        del md["b"]


def test_delete_locally_overridden_stored_key():
    md = make_md({"a": 1})
    md["a"] = 0
    del md["a"]
    assert "a" not in md
    assert md.get("a", 5) == 5
    assert md.to_dict() == {}


def test_iteration_order():
    md = make_md({"a": 1, "b": 2})
    md["c"] = 3
    md["a"] = 0
    assert list(md) == ["a", "b", "c"]
    assert list(md.to_dict()) == ["a", "b", "c"]
    assert md.to_dict() == {"a": 0, "b": 2, "c": 3}
    assert len(md) == 3


def test_iteration_order_stored_then_local():
    md = make_md({"a": 1, "b": 2}, {"d": 4})
    md["c"] = 3
    del md["a"]
    md["a"] = 5
    assert list(md) == ["a", "b", "d", "c"]
    assert list(md.to_dict()) == ["a", "b", "d", "c"]


def test_write_during_iteration():
    md = make_md({"a": 1}, {"b": 2})
    seen = []
    for key in md:
        seen.append(key)
        md["y"] = 2
    assert seen == ["a", "b"]


def test_default_construction_is_empty():
    md = MetaData()
    assert md.to_dict() == {}
    assert len(md) == 0


def test_local_argument_is_copied():
    local = {"b": 2}
    md = make_md(local=local)
    md["b"] = 3
    assert local == {"b": 2}


def test_clear():
    stored = MappingProxyType({"a": 1, "b": 2})
    md = MetaData(stored=stored, local={"c": 3})
    md.clear()
    assert md.to_dict() == {}
    assert len(md) == 0
    assert dict(stored) == {"a": 1, "b": 2}
    md["b"] = 0
    assert md.to_dict() == {"b": 0}


class _UnreadableValues(dict):
    def __getitem__(self, key):
        raise AssertionError(key)


def test_clear_does_not_read_values():
    md = MetaData(stored=_UnreadableValues(a=1))
    md.clear()
    assert len(md) == 0


def test_pop():
    md = make_md({"a": 1})
    assert md.pop("a") == 1
    assert "a" not in md
    assert md.pop("a", None) is None


def test_repr_shows_contents():
    md = make_md({"a": 1}, {"b": "x"})
    assert repr(md) == repr({"a": 1, "b": "x"})


def test_non_str_key_raises():
    md = make_md()
    with pytest.raises(TypeError):
        md[1] = 2


def test_non_str_key_in_local_argument_raises():
    with pytest.raises(TypeError):
        MetaData(local={1: 2})
