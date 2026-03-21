from pathlib import Path
from typing import Union

import pytest

from cogent3.app.data_store import DataMemberABC
from cogent3.app.typing import (
    AlignedSeqsType,
    IdentifierType,
    SerialisableType,
    TabularType,
    check_type_compatibility,
    defined_types,
    get_type_display_names,
    resolve_type_hint,
)

# --- Tests for resolve_type_hint ---


def test_resolve_type_hint_typevar_bound():
    """TypeVar with string bound resolves to actual class"""
    from cogent3.core.alignment import Alignment

    resolved = resolve_type_hint(AlignedSeqsType)
    assert resolved is Alignment


def test_resolve_type_hint_typevar_constraints():
    """TypeVar with string constraints resolves to Union of classes"""
    from cogent3.core.table import Table
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.util.dict_array import DictArray

    resolved = resolve_type_hint(TabularType)
    # Should be a Union of Table, DictArray, DistanceMatrix
    from typing import get_args

    args = set(get_args(resolved))
    assert args == {Table, DictArray, DistanceMatrix}


def test_resolve_type_hint_concrete_class():
    """Concrete classes pass through unchanged"""
    resolved = resolve_type_hint(int)
    assert resolved is int


def test_resolve_type_hint_protocol():
    """Protocol classes pass through unchanged"""
    resolved = resolve_type_hint(SerialisableType)
    assert resolved is SerialisableType


def test_resolve_type_hint_union():
    """Union types are resolved recursively"""
    resolved = resolve_type_hint(IdentifierType)
    from typing import get_args

    args = set(get_args(resolved))
    assert args == {str, Path, DataMemberABC}


def test_resolve_type_hint_container():
    """Container types resolved recursively"""
    from cogent3.core.alignment import Alignment

    resolved = resolve_type_hint(list[AlignedSeqsType])
    assert resolved == list[Alignment]


def test_resolve_type_hint_unresolvable():
    """Unresolvable string raises TypeError"""
    with pytest.raises(TypeError, match="cannot resolve"):
        resolve_type_hint("NoSuchType")


def test_resolve_type_hint_user_module_globals():
    """module_globals are checked first for resolution"""

    class MyCustomType:
        pass

    resolved = resolve_type_hint("MyCustomType", {"MyCustomType": MyCustomType})
    assert resolved is MyCustomType


# --- Tests for SerialisableType Protocol ---


def test_serialisable_type_isinstance():
    """objects with to_rich_dict satisfy SerialisableType"""
    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs({"s1": "ACGT", "s2": "ACGT"}, moltype="dna")
    assert isinstance(aln, SerialisableType)


def test_serialisable_type_not_isinstance():
    """objects without to_rich_dict do not satisfy SerialisableType"""
    assert not isinstance("hello", SerialisableType)
    assert not isinstance(42, SerialisableType)


def test_serialisable_type_custom_class():
    """custom class with to_rich_dict satisfies SerialisableType"""

    class MyObj:
        def to_rich_dict(self) -> dict:
            return {}

    assert isinstance(MyObj(), SerialisableType)


# --- Tests for get_type_display_names ---


def test_get_type_display_names_concrete():
    """concrete class returns its name"""
    names = get_type_display_names(int)
    assert names == frozenset({"int"})


def test_get_type_display_names_union():
    """Union returns names of all constituents"""
    names = get_type_display_names(Union[str, int])  # noqa: UP007
    assert names == frozenset({"str", "int"})


def test_get_type_display_names_protocol():
    """Protocol returns its own name"""
    names = get_type_display_names(SerialisableType)
    assert names == frozenset({"SerialisableType"})


def test_get_type_display_names_resolved_typevar():
    """resolved TypeVar (now a class) returns the class name"""

    resolved = resolve_type_hint(AlignedSeqsType)
    names = get_type_display_names(resolved)
    assert names == frozenset({"Alignment"})


# --- Tests for check_type_compatibility ---


def test_check_type_compatibility_same_type():
    """same type is compatible"""
    from cogent3.core.alignment import Alignment

    assert check_type_compatibility(Alignment, Alignment) is True


def test_check_type_compatibility_subclass():
    """subclass is compatible"""
    from cogent3.core.alignment import Alignment, CollectionBase

    assert check_type_compatibility(Alignment, CollectionBase) is True


def test_check_type_compatibility_protocol_input():
    """Protocol on input side is lenient"""
    assert check_type_compatibility(int, SerialisableType) is True


def test_check_type_compatibility_protocol_return():
    """Protocol on return side is lenient"""
    assert check_type_compatibility(SerialisableType, int) is True


def test_check_type_compatibility_incompatible():
    """incompatible concrete types"""
    assert check_type_compatibility(int, str) is False


# --- Tests for defined_types ---


def test_defined_types():
    types = defined_types()
    # TabularType should resolve to Table, DictArray, DistanceMatrix
    includes = types["TabularType"][0, "includes"]
    assert len(includes.split(",")) == 3
