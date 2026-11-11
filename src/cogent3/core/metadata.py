"""Metadata attached to sequence collections and their storage."""

from collections.abc import Iterator, Mapping, MutableMapping
from types import MappingProxyType
from typing import Any, Final

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


def _check_key(key: object) -> None:
    if not isinstance(key, str):
        msg = f"metadata keys must be str, not {type(key).__name__} {key!r}"
        raise TypeError(msg)


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

    def __getitem__(self, key: str) -> Any:  # noqa: ANN401
        value = self._local.get(key, _MISSING)
        if value is _HIDDEN:
            raise KeyError(key)
        if value is not _MISSING:
            return value
        return self._stored[key]

    def __setitem__(self, key: str, value: Any) -> None:  # noqa: ANN401
        _check_key(key)
        self._local[key] = value

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

    def to_dict(self) -> dict[str, Any]:
        """Return the visible contents as a plain dict."""
        local = self._local.copy()
        result = {}
        for key, stored_value in self._stored.items():
            value = local.get(key, _MISSING)
            if value is _MISSING:
                result[key] = stored_value
            elif value is not _HIDDEN:
                result[key] = value
        result.update(
            (key, value)
            for key, value in local.items()
            if value is not _HIDDEN and key not in self._stored
        )
        return result
