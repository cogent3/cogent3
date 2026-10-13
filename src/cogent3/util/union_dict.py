"""UnionDict extension of dict."""

from typing import Any, Self


def _convert_children(root: "UnionDict", source: Any) -> None:
    # keyed by id so a dict shared or referenced from within itself is copied once
    copies: dict[int, UnionDict] = {id(source): root}
    stack = [root]
    while stack:
        node = stack.pop()
        for key, value in node.items():
            if not isinstance(value, dict):
                continue
            child = copies.get(id(value))
            if child is None:
                child = UnionDict()
                dict.update(child, value)
                copies[id(value)] = child
                stack.append(child)
            dict.__setitem__(node, key, child)


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
        if isinstance(value, dict) and not isinstance(value, UnionDict):
            value = UnionDict(value)
        super().__setitem__(key, value)

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
        other_copy = UnionDict(other)
        merged = {(id(self), id(other_copy))}
        stack: list[tuple[dict[Any, Any], dict[Any, Any]]] = [(self, other_copy)]
        while stack:
            target, source = stack.pop()
            for key, value in source.items():
                current = target.get(key)
                if not (isinstance(value, dict) and isinstance(current, dict)):
                    dict.__setitem__(target, key, value)
                elif (id(current), id(value)) not in merged:
                    merged.add((id(current), id(value)))
                    stack.append((current, value))

    def update(self, *args: Any, **kwargs: Any) -> None:
        """Converts to UnionDict."""
        assert not (args and kwargs)
        if args:
            assert len(args) == 1
            kwargs = args[0]
        super().update(UnionDict(kwargs))
