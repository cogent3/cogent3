from __future__ import annotations

from collections.abc import Callable, Iterable, Sized
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    Protocol,
    TypeVar,
    cast,
    overload,
    runtime_checkable,
)
from uuid import UUID, uuid4

if TYPE_CHECKING:
    from pathlib import Path

    from cogent3.app.comp_new import NotCompleted, _AppBaseClass
    from cogent3.app.data_store import DataMemberABC, DataStoreABC

T = TypeVar("T")
U = TypeVar("U")
R = TypeVar("R")


class source_proxy(Generic[T]):
    __slots__ = ("_obj", "_src", "_uuid")

    def __init__(
        self, obj: T, src: object | str | Path | None = None, uuid: UUID | None = None
    ) -> None:
        self._obj: T = obj
        self._src: object | str | Path = obj if src is None else src
        self._uuid: UUID = uuid4() if uuid is None else uuid

    def __hash__(self) -> int:
        return hash(self._uuid)

    @property
    def obj(self) -> T:
        return self._obj

    def set_obj(self, obj: T) -> None:
        self._obj = obj

    def propagate_source(self, new_obj: R) -> "source_proxy[R]":
        return source_proxy(new_obj, self._src, self._uuid)

    @property
    def source(self) -> object | str | Path:
        """origin of this object"""
        return self._src

    @source.setter
    def source(self, src: T | str | Path) -> None:
        # need to check whether src is hashable, how to cope if it isn't?
        # might need to make this instance hashable perhaps using a uuid?
        self._src = src

    @property
    def uuid(self) -> str:
        """unique identifier for this object"""
        return str(self._uuid)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._obj, name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_"):
            super().__setattr__(name, value)
        else:
            setattr(self._obj, name, value)

    def __bool__(self) -> bool:
        return bool(self._obj)

    def __repr__(self) -> str:
        return self.obj.__repr__()

    def __str__(self) -> str:
        return self.obj.__str__()

    def __eq__(self, other: object) -> bool:
        return self.obj.__eq__(other)

    def __len__(self) -> int:
        if not isinstance(self.obj, Sized):
            msg = f"Object of type {type(self.obj).__name__} has no len()"
            raise TypeError(msg)
        return len(self.obj)

    # pickling induces infinite recursion on python 3.10
    # only on Windows, so implementing the following methods explicitly
    def __getstate__(self) -> tuple[T, object | str | Path, UUID]:
        return self._obj, self._src, self._uuid

    def __setstate__(self, state: tuple[T, object | str | Path, UUID]) -> None:
        self._obj, self._src, self._uuid = state


@runtime_checkable
class HasSource(Protocol[T]):
    source: str | Path | T


@overload
def proxy_input(dstore: Iterable[str | Path]) -> list[source_proxy[str | Path]]: ...
@overload
def proxy_input(
    dstore: DataStoreABC | Iterable[DataMemberABC],
) -> list[source_proxy[DataMemberABC]]: ...
@overload
def proxy_input(dstore: Iterable[HasSource[T]]) -> list[HasSource[T]]: ...


def proxy_input(
    dstore: Iterable[str | Path]
    | Iterable[DataMemberABC]
    | DataStoreABC
    | Iterable[HasSource[T]],
) -> (
    list[source_proxy[str | Path]]
    | list[source_proxy[DataMemberABC]]
    | list[HasSource[T]]
):
    inputs: list[source_proxy[str | Path | DataMemberABC] | HasSource[T]] = []
    for e in dstore:
        if not e:
            continue
        if not isinstance(e, source_proxy):
            e = e if isinstance(e, HasSource) else source_proxy(e)
        inputs.append(e)

    return cast(
        "list[source_proxy[str | Path]] | list[source_proxy[DataMemberABC]] | list[HasSource[T]]",
        inputs,
    )


class propagate_source(Generic[T, R]):
    """retains result association with source

    Notes
    -----
    Returns the unwrapped result if it has a .source instance,
    otherwise returns the original source_proxy with the .obj
    updated with result.
    """

    def __init__(
        self,
        app: _AppBaseClass[T, Any, R],
        id_from_source: Callable[[R | HasSource[Any] | NotCompleted], str | None],
    ) -> None:
        self.app = app
        self.id_from_source = id_from_source

    def __call__(
        self, value: source_proxy[T] | T
    ) -> source_proxy[R] | R | NotCompleted:
        from cogent3.app.comp_new import NotCompleted

        if not isinstance(value, source_proxy):
            return self.app(value)

        value = cast("source_proxy[Any]", value)

        result = self.app(value.obj)
        if self.id_from_source(result):
            return result
        if isinstance(result, NotCompleted):
            result.source = value.source
            return result

        return value.propagate_source(result)
