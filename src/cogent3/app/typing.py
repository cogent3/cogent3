"""defined type hints for app composability"""

# todo write more extensive docstring explaining limited use of these types
from __future__ import annotations

import inspect
import re
import sys
from typing import ForwardRef, TypeVar, Union

from typing_extensions import get_args, get_origin

if sys.version_info.minor >= 10:
    from types import UnionType

    NESTED_HINTS = (Union, UnionType, list, tuple, set)
else:
    NESTED_HINTS = (Union, list, tuple, set)

AlignedSeqsType = TypeVar("AlignedSeqsType", "Alignment", "ArrayAlignment")
UnalignedSeqsType = TypeVar("UnalignedSeqsType", bound="SequenceCollection")
SeqsCollectionType = Union[AlignedSeqsType, UnalignedSeqsType]
SeqType = TypeVar(
    "SeqType",
    "Sequence",
    "DnaSequence",
    "RnaSequence",
    "ByteSequence",
    "ProteinSequence",
    "ProteinWithStopSequence",
)
PairwiseDistanceType = TypeVar("PairwiseDistanceType", bound="DistanceMatrix")
TabularType = TypeVar("TabularType", "Table", "DictArray", "DistanceMatrix")
TreeType = TypeVar("TreeType", "TreeNode", "PhyloNode")
SerialisableType = TypeVar("SerialisableType")
BootstrapResultType = TypeVar("BootstrapResultType", bound="bootstrap_result")
HypothesisResultType = TypeVar("HypothesisResultType", bound="hypothesis_result")
ModelCollectionResultType = TypeVar(
    "ModelCollectionResultType", bound="model_collection_result"
)
ModelResultType = TypeVar("ModelResultType", bound="model_result")
TabularResultType = TypeVar("TabularResultType", bound="tabular_result")
GenericResultType = TypeVar("GenericResultType", bound="generic_result")
ResultType = Union[
    GenericResultType,
    BootstrapResultType,
    HypothesisResultType,
    ModelResultType,
    TabularResultType,
]

IdentifierType = TypeVar("IdentifierType")


def _is_type(text):
    p = re.compile("[A-Z][a-z]+")
    matches = list(p.finditer(text))
    if len(matches) <= 1 or matches[0].start() != 0:
        return False

    return matches[-1].group() == "Type"


_all_types = {n: t for n, t in locals().items() if _is_type(n)}


def get_constraint_names(*hints) -> set[str | type]:
    """returns the set of named constraints of a type hint"""
    all_hints = set()
    for hint in hints:
        if hint in (SerialisableType, IdentifierType) or (
            inspect.isclass(hint) and get_origin(hint) not in (list, tuple, set)
        ):
            all_hints.add(hint.__name__)
            continue

        if getattr(hint, "__bound__", None):
            all_hints.add(hint.__bound__)
            continue

        if getattr(hint, "__constraints__", None):
            all_hints.update(hint.__constraints__)
            continue

        if get_origin(hint) in NESTED_HINTS:
            all_hints.update(get_constraint_names(*get_args(hint)))

        if type(hint) == type:
            all_hints.add(hint.__name__)
        elif type(hint) == ForwardRef:
            all_hints.add(hint.__forward_arg__)
        elif type(hint) == str:
            all_hints.add(hint)

    return {h.__forward_arg__ if type(h) == ForwardRef else h for h in all_hints}


def type_tree(hint, depth=0) -> tuple:
    """compute the order of types"""
    level_type = get_origin(hint)
    if not level_type:
        return depth + 1, hint

    levels = []
    depths = []
    for arg in get_args(hint):
        d, t = type_tree(arg, depth=depth)
        levels.append(t)
        depths.append(d)
    depth = max(depths) + 1

    if len(levels) == 1:
        levels = levels[0]

    try:
        levels = tuple(levels)
    except TypeError:
        levels = (levels,)

    return depth, (level_type, levels)


def defined_types():
    """returns a table of the type hints and the cogent3 classes they represent

    Notes
    -----
    These (or standard Python) types are required to annotate argument and
    return values from cogent3 apps. They define the compatability of apps.
    """
    from cogent3.util.table import Table

    rows = [[n, ", ".join(get_constraint_names(t))] for n, t in _all_types.items()]
    title = "To use a type hint, from cogent3.app import typing"
    legend = (
        "An app which uses one of these hints is compatible with the indicated types."
    )
    table = Table(
        header=["type hint", "includes"],
        data=rows,
        title=title,
        legend=legend,
        index_name="type hint",
    )
    table.set_repr_policy(show_shape=False)
    return table
