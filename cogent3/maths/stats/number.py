from collections.abc import MutableMapping, Mapping
import numpy

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class SummaryStatBase:
    @property
    def mean(self):
        stat = 0
        if len(self) > 0:
            stat = numpy.mean(self.expanded_values())
        return stat

    @property
    def std(self):
        stat = self.var
        if stat > 0:
            stat = numpy.sqrt(stat)
        return stat

    @property
    def var(self):
        stat = 0
        if len(self) > 0:
            stat = numpy.var(self.expanded_values(), ddof=1)
        return stat

    def quantile(self, q):
        stat = 0
        if len(self) > 0:
            stat = numpy.quantile(self.expanded_values(), q=q)
        return stat

    @property
    def median(self):
        stat = 0
        if len(self) > 0:
            stat = numpy.median(self.expanded_values())
        return stat

    @property
    def mode(self):
        _, mode = max((v, k) for k, v in self.items())
        return mode

    @property
    def sum(self):
        stat = 0
        if len(self) > 0:
            stat = numpy.sum(self.expanded_values())
        return stat


class CategoryCounter(MutableMapping, SummaryStatBase):
    """counting class with summary statistic attributes"""

    def __init__(self, data):
        for element in data:
            self[element] += 1

    def __len__(self):
        return sum(self.values())

    def __add__(self, other):
        self[other] += 1
        return self

    def __sub__(self, other):
        self[other] -= 1
        if self[other] == 0:
            del self[other]
        return self

    def tolist(self, keys=None):
        """return values for these keys as a list"""
        if keys is None:
            keys = list(self)
        result = [self[key] for key in keys]
        return result

    def toarray(self, keys=None):
        """return values for these keys as an array"""
        data = self.tolist(keys=keys)
        data = numpy.array(data, dtype=int)
        return data

    def expanded_values(self):
        values = list(self.values())
        return values

    def __setitem__(self, key, val):
        self.__dict__[key] = val

    def __getitem__(self, key):
        val = 0 if key not in self.__dict__ else self.__dict__[key]
        return val

    def __delitem__(self, key):
        del (self.__dict__[key])

    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        return iter(self.__dict__)

    def __repr__(self):
        return repr(self.__dict__)

    def todict(self):
        return dict(self)

    def tolist(self, keys=None):
        """return values for these keys as a list"""
        if keys is None:
            keys = list(self)
        result = [self[key] for key in keys]
        return result

    def toarray(self, keys=None):
        """return just these keys as an array"""
        if keys is None:
            keys = list(self)
        data = self.tolist(keys)
        data = numpy.array(data, dtype=float)
        return data


class NumberCounter(CategoryCounter):
    def __init__(self, data):
        super(NumberCounter, self).__init__(data)

    @property
    def valid(self):
        types = set(map(type, self))
        if types <= {int, float, complex}:
            result = True
        else:
            key = next(iter(self))
            try:  # if a numpy type
                result = key.dtype.kind in 'uifc'
            except AttributeError:
                result = False

        return result

    def expanded_values(self, check=False):
        # todo memory footprint can be improved by directly computing the
        #  summary statistics
        if check and not self.valid:
            raise ValueError('non-numeric keys')
        values = []
        for k, v in self.items():
            values.extend([k] * v)
        return values

    def __len__(self):
        return sum(self.values())

            self[k] += v
