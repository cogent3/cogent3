"""Metadata attached to sequence collections and their storage."""

import copy
import functools
from collections.abc import Iterable, Iterator, Mapping, MutableMapping
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final, Self

if TYPE_CHECKING:
    from _typeshed import SupportsKeysAndGetItem

_NO_METADATA: Mapping[str, Any] = MappingProxyType({})


class _Hidden:
    __slots__ = ()

    def __reduce__(self) -> str:
        return "_HIDDEN"

    def __repr__(self) -> str:
        return "_HIDDEN"


# a local value that hides the stored key of the same name
_HIDDEN: Final = _Hidden()
_MISSING: Final = object()


_IMMUTABLE_SCALARS: Final = (str, bytes, int, float, complex, type(None))


def _check_key(key: object) -> None:
    if not isinstance(key, str):
        msg = f"metadata keys must be str, not {type(key).__name__} {key!r}"
        raise TypeError(msg)


def _copy_leaf(value: Any, stored: bool) -> Any:  # noqa: ANN401
    if isinstance(value, MetaData):
        return value.to_dict()
    if stored and not isinstance(value, _IMMUTABLE_SCALARS):
        return copy.deepcopy(value)
    return value


def _copy_value(value: Any, where: str, stored: bool = False) -> Any:  # noqa: ANN401
    """Return a copy of value with every mapping a plain dict.

    Lists and tuples are copied, a MetaData becomes its to_dict(). Other
    objects are shared, or deep-copied when stored is True. Raises TypeError
    naming the path to a mapping key that is not a str, unless stored is True,
    and ValueError for a value that contains itself.
    """
    root: list[Any] = [None]
    tuples: list[tuple[list[Any], Any, Any]] = []
    stack: list[tuple[Any, Any, Any, str, frozenset[int]]] = [
        (value, root, 0, where, frozenset()),
    ]
    while stack:
        src, parent, slot, path, ancestors = stack.pop()
        if isinstance(src, MetaData) or not (
            isinstance(src, (Mapping, list)) or type(src) is tuple
        ):
            parent[slot] = _copy_leaf(src, stored)
            continue
        if id(src) in ancestors:
            msg = f"metadata value {path} contains itself"
            raise ValueError(msg)
        inner = ancestors | {id(src)}
        if isinstance(src, Mapping):
            mapping: dict[str, Any] = {}
            parent[slot] = mapping
            for key, item in src.items():
                if not stored and not isinstance(key, str):
                    location = f"in {path}" if path else "at the top level"
                    msg = f"metadata keys must be str, found {key!r} {location}"
                    raise TypeError(msg)
                mapping[key] = None
                item_path = f"{path}.{key}" if path else key
                stack.append((item, mapping, key, item_path, inner))
            continue
        sequence: list[Any] = [None] * len(src)
        parent[slot] = sequence
        if isinstance(src, tuple):
            tuples.append((sequence, parent, slot))
        stack.extend(
            (item, sequence, index, f"{path}[{index}]", inner)
            for index, item in enumerate(src)
        )
    # a tuple is built after its items, and items are created after their parent
    for sequence, parent, slot in reversed(tuples):
        parent[slot] = tuple(sequence)
    return root[0]


@functools.cache
def _class_attributes(cls: type) -> frozenset[str]:
    # attribute lookup on an instance searches the MRO, not the metaclass
    return frozenset(name for klass in cls.__mro__ for name in vars(klass))


def _class_attribute_error(name: str) -> AttributeError:
    return AttributeError(f"{name!r} is a MetaData attribute, use md[{name!r}]")


def _no_key_error(name: str) -> AttributeError:
    return AttributeError(f"MetaData has no key {name!r}")


def _read_copy(value: Any) -> Any:  # noqa: ANN401
    if isinstance(value, Mapping):
        return MetaData(stored=value)
    return _copy_value(value, "", stored=True)


class MetaData(MutableMapping[str, Any]):
    """The metadata of one sequence collection or sequence.

    Reads check this object's local changes, then the stored layer, a
    read-only mapping shared with the storage. Writes and deletions go to the
    local layer only, so the stored layer is never changed. Deleting a stored
    key hides it.

    A ``MetaData`` written by several threads needs a lock held by the caller.
    """

    __slots__ = ("_inherited", "_local", "_stored")

    def __init__(
        self,
        stored: Mapping[str, Any] | None = None,
        local: Mapping[str, Any] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        stored
            the read-only mapping held by the storage, never copied
        local
            initial local values, which are copied
        """
        self._stored: Mapping[str, Any] = _NO_METADATA if stored is None else stored
        self._local: dict[str, Any] = {}
        self._inherited = False
        if local:
            self.update(local)

    @classmethod
    def _sharing(cls, local: dict[str, Any]) -> "MetaData":
        # edits through the result change local, which its owner holds
        result = cls()
        result._local = local
        return result

    def __getitem__(self, key: str) -> Any:  # noqa: ANN401
        value = self._local.get(key, _MISSING)
        if value is _HIDDEN:
            raise KeyError(key)
        if value is _MISSING:
            value = self._stored[key]
            if isinstance(value, _IMMUTABLE_SCALARS):
                return value
            # one dict operation, so concurrent readers get the same copy
            value = self._local.setdefault(key, _read_copy(value))
        return MetaData._sharing(value) if type(value) is dict else value

    def __setitem__(self, key: str, value: Any) -> None:  # noqa: ANN401
        _check_key(key)
        current = self._local.get(key, _MISSING)
        # md[key] |= other assigns md[key] back to itself
        if isinstance(value, MetaData) and (
            value is current or value._local is current
        ):
            return
        self._local[key] = _copy_value(value, key)

    def __ior__(self, other: Mapping[str, Any]) -> Self:
        """Merge other into this object.

        Nested dicts are merged, and other's values win everywhere else.
        """
        if not isinstance(other, Mapping):
            return NotImplemented
        self._merge(_copy_value(other, ""))
        return self

    def _merge(self, incoming: dict[str, Any]) -> None:
        # incoming is already a validated copy, so it is installed as it is
        stack: list[tuple[MetaData, dict[str, Any]]] = [(self, incoming)]
        while stack:
            target, values = stack.pop()
            for key, value in values.items():
                if type(value) is dict and target._holds_dict(key):
                    stack.append((target[key], value))
                else:
                    target._local[key] = value

    def _holds_dict(self, key: str) -> bool:
        # checks without reading, so a stored list is not copied only to be replaced
        value = self._local.get(key, _MISSING)
        if value is _MISSING:
            value = self._stored.get(key)
        return isinstance(value, (MetaData, Mapping))

    def update(
        self,
        other: "SupportsKeysAndGetItem[str, Any] | Iterable[tuple[str, Any]]" = (),
        /,
        **kwargs: object,
    ) -> None:
        """Merge into this object as ``|=`` does.

        Takes the same arguments as ``dict.update``.
        """
        incoming = _copy_value(other if isinstance(other, Mapping) else dict(other), "")
        incoming.update(_copy_value(kwargs, ""))
        self._merge(incoming)

    def __getattr__(self, name: str) -> Any:  # noqa: ANN401
        if name.startswith("_"):
            raise _no_key_error(name)
        try:
            return self[name]
        except KeyError:
            raise _no_key_error(name) from None

    def __setattr__(self, name: str, value: Any) -> None:  # noqa: ANN401
        if name.startswith("_"):
            object.__setattr__(self, name, value)
        elif name in _class_attributes(type(self)):
            raise _class_attribute_error(name)
        else:
            self[name] = value

    def __delattr__(self, name: str) -> None:
        if name.startswith("_"):
            object.__delattr__(self, name)
        elif name in _class_attributes(type(self)):
            raise _class_attribute_error(name)
        else:
            try:
                del self[name]
            except KeyError:
                raise _no_key_error(name) from None

    def __delitem__(self, key: str) -> None:
        if key in self._stored:
            if self._local.get(key) is _HIDDEN:
                raise KeyError(key)
            self._local[key] = _HIDDEN
        elif self._local.pop(key, _HIDDEN) is _HIDDEN:
            raise KeyError(key)

    def __contains__(self, key: object) -> bool:
        value = self._local.get(key, _MISSING)  # type: ignore[call-overload]
        if value is _HIDDEN:
            return False
        return value is not _MISSING or key in self._stored

    def __iter__(self) -> Iterator[str]:
        local = self._local.copy()
        for key in self._stored:
            if local.get(key) is not _HIDDEN:
                yield key
        for key, value in local.items():
            if value is not _HIDDEN and key not in self._stored:
                yield key

    def __len__(self) -> int:
        return sum(1 for _ in self)

    def __repr__(self) -> str:
        return repr(self.to_dict())

    def clear(self) -> None:
        """Remove every key."""
        local = self._local
        for key in list(local):
            if key not in self._stored:
                local.pop(key, None)
        for key in self._stored:
            local[key] = _HIDDEN

    def setdefault(self, key: str, default: Any = None) -> Any:  # noqa: ANN401
        """Return md[key], first setting it to a copy of default if absent."""
        if key not in self:
            self[key] = default
        return self[key]

    def to_dict(self) -> dict[str, Any]:
        """Return the visible contents as a plain nested dict, a copy."""
        result: dict[str, Any] = {}
        stack: list[tuple[MetaData, dict[str, Any]]] = [(self, result)]
        while stack:
            md, target = stack.pop()
            local = md._local.copy()
            stored = md._stored
            visible = [
                (key, local.get(key, value), key not in local)
                for key, value in stored.items()
                if local.get(key) is not _HIDDEN
            ]
            visible.extend(
                (key, value, False)
                for key, value in local.items()
                if value is not _HIDDEN and key not in stored
            )
            for key, value, from_stored in visible:
                if isinstance(value, MetaData):
                    target[key] = {}
                    stack.append((value, target[key]))
                else:
                    target[key] = _copy_value(value, key, stored=from_stored)
        return result
