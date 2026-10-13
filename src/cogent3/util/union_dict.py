"""UnionDict extension of dict."""

from typing import Any, Self


def _convert_children(root: "UnionDict", source: Any) -> None:
    on_path = {id(source)}
    frames = [(root, iter(root.items()), id(source))]
    while frames:
        node, items, source_id = frames[-1]
        for key, value in items:
            if not isinstance(value, dict):
                continue
            if id(value) in on_path:
                msg = f"cannot copy {key!r}, it refers back to a dict containing it"
                raise ValueError(msg)
            child = UnionDict()
            dict.update(child, value)
            dict.__setitem__(node, key, child)
            on_path.add(id(value))
            frames.append((child, iter(child.items()), id(value)))
            break
        else:
            frames.pop()
            on_path.discard(source_id)


class UnionDict(dict):
    """dictionary class that can be updated as a union with another dict
    entries are also accessible directly as attributes on the object"""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """

        Parameters
        ----------
        data : dict
            nested dicts are copied, never shared with data
        """
        assert not (args and kwargs)
        if args:
            assert len(args) == 1
            kwargs = args[0]
        super().__init__(kwargs)
        _convert_children(self, kwargs)

    def __getattr__(self, item):
        if item in self:
            return self.get(item)

        try:
            return super().__getattr__(item)
        except AttributeError:
            msg = f"'{item}' not a key or attribute"
            raise AttributeError(msg)

    def __setattr__(self, key: str, value: Any) -> None:
        self[key] = value

    def __setitem__(self, key: Any, value: Any) -> None:
        # a |= on an attribute assigns the same object back, which must not copy it
        stored_as_is = value is self.get(key) and isinstance(value, UnionDict)
        if isinstance(value, dict) and not stored_as_is:
            value = UnionDict(value)
        super().__setitem__(key, value)

    def setdefault(self, key: Any, default: Any = None) -> Any:
        if key not in self:
            self[key] = default
        return self[key]

    def __or__(self, other: Any) -> Self:
        result = self.__class__(self)
        result.union(other)
        return result

    def __ior__(self, other: Any) -> Self:
        self.union(other)
        return self

    def union(self, other: Any) -> None:
        """merge other into self in place

        keys unique to other are introduced, keys in self shared with other
        are updated with the value from other, except where both values are
        dicts, which are merged in the same way"""
        stack: list[tuple[dict[Any, Any], dict[Any, Any]]] = [(self, UnionDict(other))]
        while stack:
            target, source = stack.pop()
            for key, value in source.items():
                current = target.get(key)
                if isinstance(value, dict) and isinstance(current, dict):
                    stack.append((current, value))
                else:
                    dict.__setitem__(target, key, value)

    def update(self, *args: Any, **kwargs: Any) -> None:
        """Converts to UnionDict."""
        assert not (args and kwargs)
        if args:
            assert len(args) == 1
            kwargs = args[0]
        super().update(UnionDict(kwargs))
