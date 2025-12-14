from typing import Union

import pytest

from cogent3.app.typing import (
    AlignedSeqsType,
    IdentifierType,
    SeqsCollectionType,
    SerialisableType,
    TabularType,
    UnalignedSeqsType,
    defined_types,
    get_constraint_names,
    type_tree,
)


def test_get_constraint_names():
    """returns the correct names"""
    from cogent3.core.alignment import (
        Alignment,
        SequenceCollection,
    )
    from cogent3.core.table import Table
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.util.dict_array import DictArray

    got = get_constraint_names(AlignedSeqsType)
    assert got == {obj.__name__ for obj in (Alignment,)}
    got = get_constraint_names(UnalignedSeqsType)
    assert got == {SequenceCollection.__name__}
    got = get_constraint_names(SeqsCollectionType)
    assert got == {obj.__name__ for obj in (Alignment, SequenceCollection)}
    got = get_constraint_names(TabularType)
    assert got == {obj.__name__ for obj in (Table, DictArray, DistanceMatrix)}


@pytest.mark.parametrize("constraint", [Union[str, bytes], str | bytes])  # noqa: UP007
def test_get_constraint_names_builtins(constraint):
    """handles built-in types"""
    expected = {"str", "bytes"}

    got = get_constraint_names(constraint)
    assert got == expected


def test_get_constraint_names_serilisable():
    """SerialisableType does not define any compatible types"""

    got = get_constraint_names(SerialisableType)
    assert got == {"SerialisableType"}


def test_get_constraint_names_identifiertype():
    """IdentifierType does not define any compatible types"""

    got = get_constraint_names(IdentifierType)
    assert got == {"IdentifierType"}


def test_get_constraint_names_mixed_serilisable_identifiertype():
    """SerialisableType does not define any compatible types"""
    expected = {"SerialisableType", "IdentifierType", "Alignment"}

    got = get_constraint_names(Union[SerialisableType, IdentifierType, AlignedSeqsType])
    assert got == expected

    if sys.version_info.minor > 9:
        got = get_constraint_names(SerialisableType | IdentifierType | AlignedSeqsType)
        assert got == expected


@pytest.mark.parametrize(
    "constraint",
    [Union[SerialisableType, "DnaSequence"], SerialisableType | "DnaSequence"],
)
def test_hints_resolved_from_str(constraint):
    got = get_constraint_names("DnaSequence")
    assert got == {"DnaSequence"}

    expected = {"SerialisableType", "DnaSequence"}
    got = get_constraint_names(constraint)
    assert expected <= got


@pytest.mark.parametrize("container", [list, tuple, set])
def test_hints_from_container_type(container):
    got = get_constraint_names(container[AlignedSeqsType])
    assert got == {"Alignment"}


@pytest.mark.parametrize("container", [list, set, tuple])
def test_hints_from_container_type_obj(container):
    got = get_constraint_names(container[AlignedSeqsType])
    assert got == {"Alignment"}


def test_hint_inherited_class():
    from collections.abc import MutableSequence

    class dummy(MutableSequence): ...

    got = get_constraint_names(dummy)
    assert got == frozenset(["dummy"])


@pytest.mark.parametrize(
    ("hint", "expect"),
    [(int, 1), (set[int], 2), (list[list[set[float]]], 4)],
)
def test_typing_tree_depth(hint, expect):
    d, _ = type_tree(hint)
    assert d == expect, (d, expect)


@pytest.mark.parametrize(
    ("hint", "expect"),
    [
        (int, int),
        (set[int], (set, (int,))),
        (list[set[int]], (list, (set, (int,)))),
    ],
)
def test_type_tree(hint, expect):
    _, t = type_tree(hint)
    assert t == expect, (t, expect)


def test_defined_types():
    types = defined_types()
    # we are checking a single value which we know has 3 entries
    # also indexing by the type name
    assert len(types["TabularType"][0, "includes"].split(",")) == 3
