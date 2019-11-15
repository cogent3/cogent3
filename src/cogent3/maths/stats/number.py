from collections.abc import Mapping, MutableMapping

import numpy

from numpy.testing import assert_allclose


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
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

    def __init__(self, data=None):
        if data is not None:
            if isinstance(data, dict):
                self.update_from_counts(data)
            else:
                self.update_from_series(data)

    def update_from_counts(self, data):
        """updates values of self using counts dict"""
        for k, v in data.items():
            self[k] += v

    def update_from_series(self, data):
        """updates values of self from raw series"""
        for element in data:
            self[element] += 1

    def expand(self):
        """returns list of [[k] * val, ..]"""
        result = []
        for k in self:
            result.extend([k] * self[k])
        return result

    def expanded_values(self):
        values = list(self.values())
        return values

    def copy(self):
        data = self.to_dict().copy()
        new = self.__class__(data)
        return new

    def __setitem__(self, key, val):
        self.__dict__[key] = val

    def __getitem__(self, key):
        val = 0 if key not in self.__dict__ else self.__dict__[key]
        return val

    def __delitem__(self, key):
        del self.__dict__[key]

    def __len__(self):
        return sum(self.values())

    def __iter__(self):
        return iter(self.__dict__)

    def __add__(self, other):
        self[other] += 1
        return self

    def __sub__(self, other):
        self[other] -= 1
        if self[other] == 0:
            del self[other]
        return self

    def __repr__(self):
        return repr(self.__dict__)

    def to_dict(self):
        return dict(self)

    def tolist(self, keys=None):
        """return values for these keys as a list"""
        if keys is None:
            keys = list(self)
        result = [self[key] for key in keys]
        return result

    def to_array(self, keys=None):
        """return values for these keys as an array"""
        data = self.tolist(keys=keys)
        data = numpy.array(data, dtype=int)
        return data

    @property
    def entropy(self):
        data = self.to_array()
        data = data / self.sum
        return -(data * numpy.log2(data)).sum()

    def to_freqs(self):
        """returns dict of {key: val/total, ..}"""
        result = CategoryFreqs(self, total=self.sum)
        return result


class CategoryFreqs(MutableMapping, SummaryStatBase):
    """category frequencies with summary statistic attributes"""

    def __init__(self, data=None, total=None, assert_unity=False):
        """
        Parameters
        ----------
        data
            data series or dict
        total
            if provided, and data is not None, elements divided by this
        assert_unity : bool
            checks sum of values (post construction) equals 1
        """
        data = data or None
        if total:
            assert data is not None
            for k, v in data.items():
                self[k] = v / total
        elif data is not None:
            for k, v in data.items():
                self[k] = v

        if assert_unity and data is not None:
            assert_allclose(self.sum, 1)

    def expanded_values(self):
        values = list(self.values())
        return values

    def copy(self):
        data = self.to_dict().copy()
        new = self.__class__(data=data)
        return new

    def __setitem__(self, key, val):
        self.__dict__[key] = val

    def __getitem__(self, key):
        val = 0 if key not in self.__dict__ else self.__dict__[key]
        return val

    def __delitem__(self, key):
        del self.__dict__[key]

    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        return iter(self.__dict__)

    def __repr__(self):
        return repr(self.__dict__)

    def to_dict(self):
        return dict(self)

    def tolist(self, keys=None):
        """return values for these keys as a list"""
        if keys is None:
            keys = list(self)
        result = [self[key] for key in keys]
        return result

    def to_array(self, keys=None):
        """return just these keys as an array"""
        data = self.tolist(keys=keys)
        data = numpy.array(data, dtype=float)
        return data

    @property
    def entropy(self):
        data = self.to_array()
        return -(data * numpy.log2(data)).sum()

    def to_normalized(self):
        """returns rescaled self so sum is 1"""
        result = CategoryFreqs(self, total=self.sum, assert_unity=True)
        return result


class NumberCounter(CategoryCounter):
    def __init__(self, data=None):
        super(NumberCounter, self).__init__(data=data)

    @property
    def valid(self):
        types = set(map(type, self))
        if types <= {int, float, complex}:
            result = True
        else:
            key = next(iter(self))
            try:  # if a numpy type
                result = key.dtype.kind in "uifc"
            except AttributeError:
                result = False

        return result

    def expanded_values(self, check=False):
        # todo memory footprint can be improved by directly computing the
        #  summary statistics
        if check and not self.valid:
            raise ValueError("non-numeric keys")
        values = []
        for k, v in self.items():
            values.extend([k] * v)
        return values

    def __len__(self):
        return sum(self.values())

    @property
    def mean(self):
        mean = sum(k * self[k] for k in self)
        return mean / len(self)

    @property
    def var(self):
        # we scale the variance contribution of a number by its occurrence
        mean = self.mean
        var = sum(self[k] * (k - mean) ** 2 for k in self)
        return var / (len(self) - 1)

    @property
    def std(self):
        var = self.var
        return numpy.sqrt(var)

    def update_from_counts(self, data):
        """updates values of self using counts dict"""
        for k, v in data.items():
            try:
                k ** 2
            except TypeError:
                raise TypeError(f"key {k} is not numeric")
            self[k] += v
