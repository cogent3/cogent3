import copy
import functools
import pickle
import re
import threading
import typing
import warnings
from types import MappingProxyType

import numpy
import pytest

from cogent3.core import metadata
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


def test_nested_stored_edit_stays_local():
    stored = {"p": {"m": "HKY85", "x": 1}}
    md = make_md(stored)
    md["p"]["m"] = "GTR"
    assert md.to_dict() == {"p": {"m": "GTR", "x": 1}}
    assert stored == {"p": {"m": "HKY85", "x": 1}}


def test_nested_stored_dict_is_metadata():
    md = make_md({"p": {"m": "HKY85"}})
    assert isinstance(md["p"], MetaData)
    assert md["p"] is md["p"]


def test_nested_stored_delete_hides():
    stored = {"p": {"m": 1, "x": 2}}
    md = make_md(stored)
    del md["p"]["m"]
    assert md.to_dict() == {"p": {"x": 2}}
    assert stored == {"p": {"m": 1, "x": 2}}


def test_reading_leaves_to_dict_unchanged():
    md = make_md({"p": {"m": 1}, "l": [1, 2]})
    md["p"]
    md["l"]
    assert md.to_dict() == {"p": {"m": 1}, "l": [1, 2]}


def test_set_nested_value_not_shared():
    md1 = make_md({"p": {"m": "HKY85"}})
    md2 = make_md()
    md2["p"] = md1["p"]
    md2["p"]["m"] = "X"
    assert md1.to_dict() == {"p": {"m": "HKY85"}}
    assert md2.to_dict() == {"p": {"m": "X"}}


def test_set_value_is_copied():
    value = {"m": [1]}
    md = make_md()
    md["p"] = value
    value["m"].append(2)
    value["n"] = 3
    assert md.to_dict() == {"p": {"m": [1]}}


def test_local_dict_edit_through_index():
    md = make_md()
    md["p"] = {"m": 1}
    assert isinstance(md["p"], MetaData)
    md["p"]["m"] = 2
    md["p"]["q"] = {"r": 1}
    md["p"]["q"]["r"] = 5
    assert md.to_dict() == {"p": {"m": 2, "q": {"r": 5}}}


def test_local_dict_replaces_stored_dict():
    md = make_md({"p": {"m": 1, "x": 2}})
    md["p"] = {"m": 3}
    assert md.to_dict() == {"p": {"m": 3}}
    assert "x" not in md["p"]


def test_to_dict_returns_copies():
    md = make_md({"p": {"m": 1}, "l": [1]}, {"q": {"n": [2]}})
    result = md.to_dict()
    result["p"]["m"] = 9
    result["l"].append(9)
    result["q"]["n"].append(9)
    assert md.to_dict() == {"p": {"m": 1}, "l": [1], "q": {"n": [2]}}


class _DictSubclass(dict):
    pass


def test_dict_subclass_value_comes_back_plain():
    md = make_md()
    md["p"] = _DictSubclass(m=_DictSubclass(n=1))
    result = md.to_dict()
    assert type(result["p"]) is dict
    assert type(result["p"]["m"]) is dict
    assert result == {"p": {"m": {"n": 1}}}


@pytest.mark.parametrize(
    ("key", "value", "path"),
    [
        ("sm", {("A", "T"): -1}, "sm"),
        ("l", [{1: 2}], "l[0]"),
        ("p", {"q": {1: 2}}, "p.q"),
        ("t", ({"a": [{2: 1}]},), "t[0].a[0]"),
    ],
)
def test_non_str_nested_key_names_path(key, value, path):
    md = make_md({key: 0})
    with pytest.raises(TypeError, match=re.escape(f"in {path}") + "$"):
        md[key] = value
    assert md.to_dict() == {key: 0}


def test_self_referencing_value_raises():
    value = []
    value.append(value)
    md = make_md()
    with pytest.raises(ValueError, match="itself"):
        md["l"] = value


def test_shared_value_is_allowed():
    shared = [1]
    md = make_md()
    md["p"] = {"a": shared, "b": shared}
    assert md.to_dict() == {"p": {"a": [1], "b": [1]}}


def test_tuple_value_kept_as_tuple():
    md = make_md()
    md["t"] = (1, [2], (3,))
    assert md.to_dict() == {"t": (1, [2], (3,))}


def test_other_objects_are_shared():
    obj = object()
    md = make_md()
    md["o"] = obj
    assert md["o"] is obj
    assert md.to_dict()["o"] is obj


def test_stored_list_copied_on_read():
    stored = {"l": [1, 2], "d": [{"a": 1}]}
    md = make_md(stored)
    md["l"].append(3)
    md["d"][0]["a"] = 5
    assert md.to_dict() == {"l": [1, 2, 3], "d": [{"a": 5}]}
    assert stored == {"l": [1, 2], "d": [{"a": 1}]}


def test_to_dict_copies_stored_objects():
    stored = {"s": {1}, "l": [{2}]}
    md = make_md(stored)
    result = md.to_dict()
    result["s"].add(9)
    result["l"][0].add(9)
    assert stored == {"s": {1}, "l": [{2}]}


def test_read_stored_list_holding_read_only_mapping():
    md = make_md({"l": [MappingProxyType({"a": 1})]})
    md["l"][0]["a"] = 2
    assert md.to_dict() == {"l": [{"a": 2}]}


def test_stored_non_str_key_not_checked():
    md = make_md({"sm": {("A", "T"): -1}})
    assert md.to_dict() == {"sm": {("A", "T"): -1}}


class _Point(typing.NamedTuple):
    x: int
    y: int


def test_named_tuple_value_kept():
    md = make_md()
    md["pt"] = _Point(1, 2)
    assert type(md.to_dict()["pt"]) is _Point


def test_setdefault_returns_stored_value():
    md = make_md()
    md.setdefault("l", []).append(1)
    assert md.to_dict() == {"l": [1]}
    md.setdefault("p", {})["m"] = 2
    assert md.to_dict() == {"l": [1], "p": {"m": 2}}
    assert md.setdefault("l", [5]) == [1]


def _paused_read_copy(read_copy, barrier, value):
    result = read_copy(value)
    barrier.wait()
    return result


def _read_and_edit(md, results, index, key):
    nested = md["p"]
    nested[key] = index
    results[index] = nested


def test_concurrent_reads_share_one_copy(monkeypatch):
    # both threads make their copy before either installs it
    md = make_md({"p": {"m": 1, "x": 2}})
    barrier = threading.Barrier(2, timeout=5)
    paused = functools.partial(_paused_read_copy, metadata._read_copy, barrier)
    monkeypatch.setattr(metadata, "_read_copy", paused)
    results = [None, None]
    first = threading.Thread(target=_read_and_edit, args=(md, results, 0, "a"))
    second = threading.Thread(target=_read_and_edit, args=(md, results, 1, "b"))
    first.start()
    second.start()
    first.join()
    second.join()
    assert results[0] is not None
    assert results[0] is results[1]
    assert md.to_dict() == {"p": {"m": 1, "x": 2, "a": 0, "b": 1}}


def test_attribute_read_nested():
    md = make_md({"p": {"m": "HKY85"}})
    assert md.p.m == "HKY85"


def test_attribute_write_nested():
    stored = {"p": {"m": "HKY85", "x": 1}}
    md = make_md(stored)
    md.p.m = "GTR"
    md.n = 2
    assert md.n == 2
    assert md.to_dict() == {"p": {"m": "GTR", "x": 1}, "n": 2}
    assert stored == {"p": {"m": "HKY85", "x": 1}}


def test_attribute_on_local_dict():
    md = make_md()
    md.q = {"r": 1}
    md.q.r = 2
    assert md.q.r == 2
    assert md.to_dict() == {"q": {"r": 2}}


def test_missing_attribute():
    md = make_md({"a": 1})
    with pytest.raises(AttributeError):
        _ = md.missing
    assert md.get("missing") is None
    assert not hasattr(md, "missing")


def test_hidden_key_attribute_missing():
    md = make_md({"a": 1})
    del md.a
    assert not hasattr(md, "a")
    assert md.to_dict() == {}
    with pytest.raises(AttributeError):
        del md.a


def test_underscore_names_not_keys():
    md = make_md({"_p": 1})
    assert md["_p"] == 1
    with pytest.raises(AttributeError):
        _ = md._p
    with pytest.raises(AttributeError):
        md._p = 2
    assert md.to_dict() == {"_p": 1}


@pytest.mark.parametrize("name", ["register", "mro"])
def test_metaclass_names_are_keys(name):
    md = make_md({name: 1})
    assert getattr(md, name) == 1
    setattr(md, name, 2)
    assert md.to_dict() == {name: 2}
    delattr(md, name)
    assert md.to_dict() == {}


@pytest.mark.parametrize("name", ["keys", "to_dict", "get", "update"])
def test_method_names_only_by_index(name):
    md = make_md({name: 1})
    assert md[name] == 1
    assert callable(getattr(md, name))
    with pytest.raises(AttributeError):
        setattr(md, name, 2)
    with pytest.raises(AttributeError):
        delattr(md, name)
    assert md.to_dict() == {name: 1}


def _merge_ior(md, other):
    md |= other


def _merge_update(md, other):
    md.update(other)


@pytest.fixture(params=[_merge_ior, _merge_update])
def merge(request):
    return request.param


@pytest.mark.parametrize("as_local", [False, True])
def test_merge_nested(merge, as_local):
    value = {"p": {"m": 1, "x": 2, "q": {"a": 1, "b": 2}}}
    md = make_md(local=value) if as_local else make_md(value)
    merge(md, {"p": {"m": 3, "q": {"b": 3}}, "n": 4})
    assert md.to_dict() == {"p": {"m": 3, "x": 2, "q": {"a": 1, "b": 3}}, "n": 4}


def test_merge_replaces_non_dicts(merge):
    md = make_md({"p": {"m": 1}, "a": 1})
    merge(md, {"p": 5, "a": {"b": 1}})
    assert md.to_dict() == {"p": 5, "a": {"b": 1}}


def test_merge_into_nested_with_hidden_key(merge):
    md = make_md({"p": {"m": 1, "x": 2}})
    del md["p"]["x"]
    merge(md, {"p": {"y": 1}})
    assert md.to_dict() == {"p": {"m": 1, "y": 1}}


def test_merge_metadata_not_shared(merge):
    md1 = make_md({"p": {"m": 1, "l": [1]}})
    md2 = make_md()
    merge(md2, md1)
    md2["p"]["m"] = "X"
    md2["p"]["l"].append(2)
    assert md1.to_dict() == {"p": {"m": 1, "l": [1]}}
    assert md2.to_dict() == {"p": {"m": "X", "l": [1, 2]}}


def test_merge_bad_key_leaves_unchanged(merge):
    md = make_md({"p": {"m": 1}})
    with pytest.raises(TypeError, match=re.escape("in p.q") + "$"):
        merge(md, {"a": 1, "p": {"m": 2, "q": {1: 2}}})
    assert md.to_dict() == {"p": {"m": 1}}


def test_merge_bad_top_level_key_message(merge):
    md = make_md()
    with pytest.raises(TypeError, match=r"found 1 at the top level$"):
        merge(md, {1: 2})


class _NoDeepCopy:
    def __deepcopy__(self, memo):
        raise AssertionError(memo)


def test_merge_dict_onto_stored_list_does_not_copy_it(merge):
    md = make_md({"a": [_NoDeepCopy()]})
    merge(md, {"a": {"b": 1}})
    assert md.to_dict() == {"a": {"b": 1}}


def test_set_list_edited_in_place_is_checked():
    md = make_md(local={"a": []})
    value = md["a"]
    value.append({1: 2})
    with pytest.raises(TypeError, match=re.escape("in a[0]") + "$"):
        md["a"] = value


def test_merge_does_not_warn(merge):
    md = make_md({"p": {"m": 1}})
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        merge(md, {"p": {"m": 2}})


def test_ior_returns_self():
    md = make_md()
    result = md.__ior__({"a": 1})
    assert result is md


def test_update_argument_forms():
    md = make_md()
    md.update([("a", 1)], b={"c": 2})
    md.update(b={"d": 3})
    assert md.to_dict() == {"a": 1, "b": {"c": 2, "d": 3}}


def test_nested_ior_keeps_nested_object():
    md = make_md({"p": {"m": 1, "x": 2}})
    nested = md["p"]
    md["p"] |= {"m": 3}
    md.p |= {"y": 4}
    assert md["p"] is nested
    nested["z"] = 5
    assert md.to_dict() == {"p": {"m": 3, "x": 2, "y": 4, "z": 5}}


def test_nested_ior_on_local_dict_keeps_it():
    md = make_md(local={"p": {"m": 1}})
    nested = md["p"]
    md["p"] |= {"m": 3}
    nested["z"] = 5
    assert md.to_dict() == {"p": {"m": 3, "z": 5}}


def test_reading_is_not_modified():
    md = make_md({"p": {"m": 1}, "l": [1, 2]})
    md["p"]
    md["l"]
    assert not md.modified


def test_in_place_list_edit_is_modified():
    stored = {"l": [1, 2]}
    md = make_md(stored)
    md["l"].append(3)
    assert md.modified
    assert stored == {"l": [1, 2]}


@pytest.mark.parametrize(
    "edit",
    [
        lambda md: md.__setitem__("b", 2),
        lambda md: md.__delitem__("a"),
        lambda md: md["p"].__setitem__("m", 2),
    ],
)
def test_edits_are_modified(edit):
    md = make_md({"a": 1, "p": {"m": 1}})
    edit(md)
    assert md.modified


def test_equal_value_is_not_modified():
    stored = MappingProxyType({"a": 1})
    md1 = MetaData(stored=stored)
    md2 = MetaData(stored=stored)
    md1["a"] = 1
    assert not md1.modified
    assert md1 == md2


class _RaisesOnEq:
    def __eq__(self, other):
        raise ValueError(other)

    __hash__ = None


def test_comparison_that_raises_is_modified():
    md = make_md({"x": [_RaisesOnEq()]})
    assert not md.modified
    # reading copies the list into local, so the comparison reaches __eq__
    md["x"]
    assert md.modified


def test_untouched_stored_array_is_not_modified():
    md = make_md({"a": numpy.array([1, 2]), "p": {"b": numpy.array([3])}})
    md["p"]
    assert not md.modified
    assert md == copy.copy(md)


def test_connected_to_unbound_raises():
    loaded = _round_trip(make_md({"a": 1}))
    with pytest.raises(ValueError, match="unbound"):
        MetaData.connected(loaded)


def test_unbound_equal_only_to_itself():
    first = _round_trip(make_md({"a": 1}))
    second = _round_trip(make_md({"z": 9}))
    assert first == first  # noqa: PLR0124
    assert first != second


def test_comparison_that_raises_is_not_equal():
    stored = MappingProxyType({})
    md1 = MetaData(stored=stored, local={"x": _RaisesOnEq()})
    md2 = MetaData(stored=stored, local={"x": _RaisesOnEq()})
    assert md1 != md2


def test_not_equal_to_dict():
    md = make_md({"a": 1})
    assert md != md.to_dict()
    assert not md == md.to_dict()  # noqa: SIM201


def test_equality_needs_same_stored():
    md1 = make_md({"a": 1})
    md2 = make_md({"a": 1})
    assert md1.to_dict() == md2.to_dict()
    assert md1 != md2


def test_not_hashable():
    with pytest.raises(TypeError):
        hash(make_md())


@pytest.mark.parametrize("copier", [copy.copy, copy.deepcopy])
def test_copy_shares_stored_and_separates_local(copier):
    md = make_md({"p": {"m": 1}, "l": [1]}, {"q": {"r": 1}})
    md["p"]["m"] = 2
    del md["l"]
    dup = copier(md)
    assert dup == md
    dup["p"]["m"] = 3
    dup["q"]["r"] = 3
    dup["l"] = [5]
    assert md.to_dict() == {"p": {"m": 2}, "q": {"r": 1}}
    assert dup.to_dict() == {"p": {"m": 3}, "q": {"r": 3}, "l": [5]}


def _round_trip(md, protocol=pickle.HIGHEST_PROTOCOL):
    return pickle.loads(pickle.dumps(md, protocol=protocol))  # noqa: S301


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_pickle_unbound_then_bind(protocol):
    stored = MappingProxyType({"p": {"m": 1, "x": 2}, "a": 1, "b": 2})
    md = MetaData(stored=stored)
    md["p"]["m"] = 3
    del md["a"]
    md["c"] = [4]
    loaded = _round_trip(md, protocol)
    # an unbound MetaData reads its stored layer as empty
    assert loaded.to_dict() == {"p": {"m": 3}, "c": [4]}
    loaded.bind(stored)
    assert loaded.to_dict() == {"p": {"m": 3, "x": 2}, "b": 2, "c": [4]}
    assert loaded == md


def test_bind_bound_raises():
    md = make_md({"a": 1})
    with pytest.raises(ValueError, match="bound"):
        md.bind(MappingProxyType({}))


def test_pickle_assigned_nested_value():
    md1 = make_md({"p": {"m": "HKY85"}})
    md2 = make_md()
    md2["p"] = md1["p"]
    md2["p"]["m"] = "X"
    assert _round_trip(md2).to_dict() == {"p": {"m": "X"}}


@pytest.mark.parametrize("copier", [copy.deepcopy, _round_trip])
def test_hidden_key_survives(copier):
    stored = MappingProxyType({"a": 1, "p": {"m": 1, "x": 2}})
    md = MetaData(stored=stored)
    del md["a"]
    del md["p"]["x"]
    dup = copier(md)
    if copier is _round_trip:
        dup.bind(stored)
    assert dup.to_dict() == {"p": {"m": 1}}


def test_connected_reads_parent_stored():
    parent = make_md({"a": 1, "p": {"m": 1, "x": 2}})
    parent["p"]["m"] = 2
    parent["b"] = 3
    child = MetaData.connected(parent)
    assert child.to_dict() == parent.to_dict()
    assert child.modified
    child["p"]["x"] = 5
    child["b"] = 6
    assert parent.to_dict() == {"a": 1, "p": {"m": 2, "x": 2}, "b": 3}


def test_connected_with_inherited_keys_hidden_is_not_modified():
    child = MetaData.connected(make_md({"a": 1}))
    del child["a"]
    assert not child.modified


def test_connected_with_empty_stored_and_no_changes_is_not_modified():
    assert not MetaData.connected(make_md()).modified


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_connected_survives_pickle(protocol):
    parent = make_md({"a": 1, "p": {"m": 1, "x": 2}, "h": 0})
    child = MetaData.connected(parent)
    child["p"]["m"] = 3
    del child["h"]
    loaded = _round_trip(child, protocol)
    assert loaded.to_dict() == {"a": 1, "p": {"m": 3, "x": 2}}
    assert loaded.modified
    loaded["p"]["x"] = 4
    assert loaded.to_dict() == {"a": 1, "p": {"m": 3, "x": 4}}


def test_connected_survives_copy():
    child = MetaData.connected(make_md({"a": 1}))
    dup = copy.deepcopy(child)
    assert dup.modified
    assert dup == child
