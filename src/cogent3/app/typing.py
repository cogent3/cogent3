"""defined type hints for app composability"""

# TODO write more extensive docstring explaining limited use of these types
from __future__ import annotations

import inspect
import re
from pathlib import Path
from types import UnionType
from typing import (
    TYPE_CHECKING,
    Any,
    ForwardRef,
    Protocol,
    TypeVar,
    Union,
    get_args,
    get_origin,
    runtime_checkable,
)

from cogent3.app.data_store import DataMemberABC
from cogent3.util.warning import deprecated_callable

if TYPE_CHECKING:  # pragma: no cover
    from cogent3.app.result import (
        bootstrap_result,
        generic_result,
        hypothesis_result,
        model_collection_result,
        model_result,
        tabular_result,
    )
    from cogent3.core.alignment import Alignment, SequenceCollection
    from cogent3.core.sequence import (
        ByteSequence,
        DnaSequence,
        ProteinSequence,
        ProteinWithStopSequence,
        RnaSequence,
        Sequence,
    )
    from cogent3.core.table import Table
    from cogent3.core.tree import PhyloNode
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.util.dict_array import DictArray


NESTED_HINTS = (Union, UnionType, list, tuple, set)


@runtime_checkable
class HasSource(Protocol):
    @property
    def source(self) -> Any: ...


@runtime_checkable
class HasInfo(Protocol):
    @property
    def info(self) -> HasSource: ...


AlignedSeqsType = TypeVar("AlignedSeqsType", bound="Alignment")
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
TreeType = TypeVar("TreeType", bound="PhyloNode")
BootstrapResultType = TypeVar("BootstrapResultType", bound="bootstrap_result")
HypothesisResultType = TypeVar("HypothesisResultType", bound="hypothesis_result")
ModelCollectionResultType = TypeVar(
    "ModelCollectionResultType",
    bound="model_collection_result",
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


@runtime_checkable
class SerialisableType(Protocol):
    def to_rich_dict(self) -> dict: ...


IdentifierType = Union[str, Path, DataMemberABC]


def _is_type(text):
    p = re.compile("[A-Z][a-z]+")
    matches = list(p.finditer(text))
    if len(matches) <= 1 or matches[0].start() != 0:
        return False

    return matches[-1].group() == "Type"


_all_types = {n: t for n, t in locals().items() if _is_type(n)}


@deprecated_callable(
    version="2026.6",
    reason="use get_type_display_names(resolve_type_hint(hint)) instead",
    new="get_type_display_names",
)
def get_constraint_names(*hints) -> set[str | type]:  # pragma: no cover
    """returns the set of named constraints of a type hint"""
    return _get_constraint_names(*hints)


def _get_constraint_names(*hints) -> set[str | type]:  # pragma: no cover
    """returns the set of named constraints of a type hint (internal, no warning)"""
    all_hints = set()
    for hint in hints:
        # SerialisableType is now a Protocol class
        if hint is SerialisableType:
            all_hints.add("SerialisableType")
            continue

        # IdentifierType is now a Union
        if hint is IdentifierType:
            all_hints.add("IdentifierType")
            continue

        if inspect.isclass(hint) and get_origin(hint) not in (list, tuple, set):
            all_hints.add(hint.__name__)
            continue

        if getattr(hint, "__bound__", None):
            all_hints.add(hint.__bound__)
            continue

        if getattr(hint, "__constraints__", None):
            all_hints.update(hint.__constraints__)
            continue

        if get_origin(hint) in NESTED_HINTS:
            all_hints.update(_get_constraint_names(*get_args(hint)))

        if type(hint) == type:
            all_hints.add(hint.__name__)
        elif type(hint) == ForwardRef:
            all_hints.add(hint.__forward_arg__)
        elif type(hint) == str:
            all_hints.add(hint)

    return {h.__forward_arg__ if type(h) == ForwardRef else h for h in all_hints}


@deprecated_callable(
    version="2026.6",
    reason="use typeguard.check_type() instead",
    new="typeguard.check_type",
)
def type_tree(hint, depth=0) -> tuple:  # pragma: no cover
    """compute the order of types"""
    return _type_tree(hint, depth=depth)


def _type_tree(hint, depth=0) -> tuple:  # pragma: no cover
    """compute the order of types (internal, no warning)"""
    level_type = get_origin(hint)
    if not level_type:
        return depth + 1, hint

    levels = []
    depths = []
    for arg in get_args(hint):
        d, t = _type_tree(arg, depth=depth)
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


_resolution_ns: dict | None = None


def _get_resolution_namespace() -> dict[str, type]:
    """lazily imports and caches all cogent3 types for forward-ref resolution"""
    global _resolution_ns
    if _resolution_ns is not None:
        return _resolution_ns

    from cogent3.app.result import (
        bootstrap_result,
        generic_result,
        hypothesis_result,
        model_collection_result,
        model_result,
        tabular_result,
    )
    from cogent3.core.alignment import Alignment, SequenceCollection
    from cogent3.core.sequence import (
        ByteSequence,
        DnaSequence,
        ProteinSequence,
        ProteinWithStopSequence,
        RnaSequence,
        Sequence,
    )
    from cogent3.core.table import Table
    from cogent3.core.tree import PhyloNode
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.util.dict_array import DictArray

    _resolution_ns = {
        "Alignment": Alignment,
        "SequenceCollection": SequenceCollection,
        "Sequence": Sequence,
        "DnaSequence": DnaSequence,
        "RnaSequence": RnaSequence,
        "ByteSequence": ByteSequence,
        "ProteinSequence": ProteinSequence,
        "ProteinWithStopSequence": ProteinWithStopSequence,
        "Table": Table,
        "PhyloNode": PhyloNode,
        "DistanceMatrix": DistanceMatrix,
        "DictArray": DictArray,
        "bootstrap_result": bootstrap_result,
        "generic_result": generic_result,
        "hypothesis_result": hypothesis_result,
        "model_collection_result": model_collection_result,
        "model_result": model_result,
        "tabular_result": tabular_result,
    }
    return _resolution_ns


def _resolve_name(name: str, module_globals: dict | None = None) -> type:
    """resolves a string name to a type, checking module_globals first,
    then the cogent3 resolution namespace"""
    if module_globals and name in module_globals:
        result = module_globals[name]
        if isinstance(result, type):
            return result

    ns = _get_resolution_namespace()
    if name in ns:
        return ns[name]

    msg = f"cannot resolve type name {name!r}"
    raise TypeError(msg)


def resolve_type_hint(hint, module_globals=None):
    """walks a type hint tree and resolves all forward references to classes

    Parameters
    ----------
    hint
        a type hint (TypeVar, Union, ForwardRef, str, or concrete class)
    module_globals
        optional dict of the module where the hint was defined, used
        to resolve forward references from user code
    """
    # Protocol classes (like SerialisableType) — return as-is
    if getattr(hint, "_is_protocol", False) and hint is not Protocol:
        return hint

    # TypeVar with __bound__ -> resolve bound class
    if isinstance(hint, TypeVar):
        if hint.__bound__:
            bound = hint.__bound__
            if isinstance(bound, str):
                bound = _resolve_name(bound, module_globals)
            elif isinstance(bound, ForwardRef):
                bound = _resolve_name(bound.__forward_arg__, module_globals)
            return bound
        if hint.__constraints__:
            resolved = tuple(
                resolve_type_hint(c, module_globals) for c in hint.__constraints__
            )
            return Union[resolved]  # type: ignore
        msg = f"unconstrained TypeVar {hint!r} cannot be resolved"
        raise TypeError(msg)

    # Union / UnionType -> recurse
    origin = get_origin(hint)
    if origin is Union or isinstance(hint, UnionType):
        args = tuple(resolve_type_hint(a, module_globals) for a in get_args(hint))
        return Union[args]  # type: ignore

    # Container types (list[X], tuple[X,Y], set[X])
    if origin in (list, tuple, set):
        args = tuple(resolve_type_hint(a, module_globals) for a in get_args(hint))
        return origin[args] if args else hint

    # ForwardRef
    if isinstance(hint, ForwardRef):
        return _resolve_name(hint.__forward_arg__, module_globals)

    # plain str
    return _resolve_name(hint, module_globals) if isinstance(hint, str) else hint


def get_type_display_names(hint) -> frozenset[str]:
    """extracts human-readable class names from a resolved type hint

    Parameters
    ----------
    hint
        a resolved type hint (one that has been through resolve_type_hint)
    """
    names = set()
    origin = get_origin(hint)

    if origin is Union or isinstance(hint, UnionType) or origin in (list, tuple, set):
        for arg in get_args(hint):
            names |= get_type_display_names(arg)
    elif isinstance(hint, type):
        names.add(hint.__name__)
    elif isinstance(hint, TypeVar):
        # fallback for unresolved TypeVars — shouldn't normally happen
        names.add(hint.__name__)

    return frozenset(names)


def _get_concrete_classes(hint) -> set[type]:
    """extracts concrete classes from a resolved type hint, walking Unions"""
    classes = set()
    origin = get_origin(hint)

    if origin is Union or isinstance(hint, UnionType):
        for arg in get_args(hint):
            classes |= _get_concrete_classes(arg)
    elif origin in (list, tuple, set):
        classes.add(origin)
    elif isinstance(hint, type):
        classes.add(hint)

    return classes


def _is_protocol(hint) -> bool:
    """checks if a type hint is or contains a runtime-checkable Protocol"""
    if getattr(hint, "_is_protocol", False) and hint is not Protocol:
        return True

    origin = get_origin(hint)
    if origin is Union or isinstance(hint, UnionType):
        return any(_is_protocol(a) for a in get_args(hint))

    return False


def check_type_compatibility(return_hint, input_hint) -> bool:
    """composition-time check: is the return type compatible with the input type?

    Parameters
    ----------
    return_hint
        resolved return type of the upstream app
    input_hint
        resolved input type of the downstream app

    Returns
    -------
    True if the types are compatible, False otherwise
    """
    # typing.Any is compatible with everything
    if return_hint is Any or input_hint is Any:
        return True

    # If either side is or contains a Protocol, be lenient — runtime check_type
    # provides the real safety net
    if _is_protocol(return_hint) or _is_protocol(input_hint):
        return True

    return_classes = _get_concrete_classes(return_hint)
    input_classes = _get_concrete_classes(input_hint)

    # Check if any return class is a subclass of any input class (or vice versa)
    for ret_cls in return_classes:
        for inp_cls in input_classes:
            try:
                if issubclass(ret_cls, inp_cls) or issubclass(inp_cls, ret_cls):
                    return True
            except TypeError:
                # issubclass can fail for some types
                if ret_cls is inp_cls:
                    return True

    return False


def defined_types():
    """returns a table of the type hints and the cogent3 classes they represent

    Notes
    -----
    These (or standard Python) types are required to annotate argument and
    return values from cogent3 apps. They define the compatability of apps.
    """
    from cogent3.core.table import Table

    rows = []
    for n, t in _all_types.items():
        try:
            resolved = resolve_type_hint(t)
            names = get_type_display_names(resolved)
            rows.append([n, ", ".join(sorted(names))])
        except TypeError:
            rows.append([n, n])

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
