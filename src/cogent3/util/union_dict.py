"""UnionDict extension of dict.
"""


__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class UnionDict(dict):
    """dictionary class that can be updated as a union with another dict
    entries are also accessible directly as attributes on the object"""

    def __init__(self, *args, **kwargs):
        """

        Parameters
        ----------
        data : dict or None
        """
        assert not (args and kwargs)
        if args:
            assert len(args) == 1
            kwargs = args[0]
        super().__init__(kwargs)

        # force all children to to UnionDicts
        for key in self:
            value = self._getsubattr_([], key)
            if not isinstance(value, dict):
                continue

            if not isinstance(value, UnionDict):
                self.update({key: UnionDict(value)})

    def __getattr__(self, item):
        if item in self:
            return self.get(item)

        try:
            return super().__getattr__(item)
        except AttributeError:
            raise AttributeError(f"'{item}' not a key or attribute")

    def __setattr__(self, key, value):
        if isinstance(value, dict):
            value = UnionDict(value)

        self.update({key: value})

    def __setitem__(self, key, value):
        if isinstance(value, dict):
            value = UnionDict(value)

        self.update({key: value})

    def __or__(self, other):
        result = self.__class__(self)
        result.union(other)
        return result

    def __ior__(self, other):
        self.union(other)
        return self

    def union(self, other):
        """returns the union of self with other

        keys unique to other are introduced, keys in self shared with other
        are updated with the value from other"""

        if not isinstance(other, UnionDict):
            other = UnionDict(other)

        for a in list(other):

            if self.get(a) is None or not (
                isinstance(other.get(a), dict) and isinstance(self.get(a), dict)
            ):
                self.update({a: other.get(a)})
                continue

            path = [a]
            while len(path) > 0:
                name = path.pop()
                cs = self._getsubattr_(path, name)
                cd = other._getsubattr_(path, name)
                path.append(name)

                for attr in list(cd):
                    if isinstance(cd.get(attr), dict) and isinstance(
                        cs.get(attr), dict
                    ):
                        # go into dict
                        path.append(attr)
                    else:
                        cs.update({attr: cd.get(attr)})
                path.pop()

    def _getsubattr_(self, path, name):
        """returns nested values"""
        d = self
        for node in path:
            d = d.get(node)
        return d.get(name)

    def update(self, *args, **kwargs):
        """Converts to UnionDict."""
        assert not (args and kwargs)
        if args:
            assert len(args) == 1
            kwargs = UnionDict(args[0])
        if not isinstance(kwargs, UnionDict):
            kwargs = UnionDict(kwargs)
        super().update(kwargs)
