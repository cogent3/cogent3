#!/usr/bin/env python
"""UnionDict extension of dict.
"""

__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class UnionDict(dict):

    def __init__(self, d):
        super().__init__(d)

        # force all children to to UnionDicts
        for a in list(self):

            a_value = self.get_sub_attr([], a)
            if not isinstance(a_value, dict):
                continue

            if not isinstance(a_value, UnionDict):
                self.update({a: UnionDict(a_value)})

            path = [a]
            while len(path) > 0:
                name = path.pop()
                cs = self.get_sub_attr(path, name)
                path.append(name)

                for attr in list(cs):
                    if isinstance(cs.get(attr), dict):
                        path.append(attr)
                        if not isinstance(cs.get(name), UnionDict):
                            cs.update({name: UnionDict(cs.get(name))})
                path.pop()

    def __getattr__(self, item):
        if item in list(self):
            return self.get(item)
        else:
            return super().__getattr__(item)

    def __setattr__(self, key, value):
        if isinstance(value, dict):
            value = UnionDict(value)
        if key in list(self):
            self.update({key: value})
        elif key in dir(self):
            super.__setattr__(key, value)
        else:
            self.update({key: value})

    def __or__(self, other):
        self.union(other)
        return self

    def union(self, d):

        if not isinstance(d, UnionDict):
            d = UnionDict(d)

        for a in list(d):

            if self.get(a) is None or not isinstance(d.get_sub_attr([], a), dict):
                self.update({a: d.get(a)})
                continue

            path = [a]
            while len(path) > 0:
                name = path.pop()
                cs = self.get_sub_attr(path, name)
                cd = d.get_sub_attr(path, name)
                path.append(name)

                for attr in list(cd):
                    if isinstance(cd.get(attr), dict) and cs.get(attr) is not None:
                        # go into dict
                        path.append(attr)
                    else:
                        cs.update({attr: cd.get(attr)})
                path.pop()

    def get_sub_attr(self, path, name):
        d = self
        for node in path:
            d = d.get(node)
        return d.get(name)
