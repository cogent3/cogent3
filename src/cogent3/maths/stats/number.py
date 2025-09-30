from abc import ABC, abstractmethod
from collections import Counter, defaultdict
from collections.abc import (
    Hashable,
    ItemsView,
    Iterable,
    Iterator,
    KeysView,
    Mapping,
    MutableMapping,
    ValuesView,
)
from collections.abc import Sequence as PySeq
from typing import TYPE_CHECKING, Any, Generic, TypeVar, cast

import numpy
import numpy.typing as npt
from numpy.testing import assert_allclose
from typing_extensions import Self

from cogent3.util.dict_array import DictArray

if TYPE_CHECKING:
    from cogent3.core.table import Table
    from cogent3.maths.stats.contingency import CategoryCounts

K = TypeVar("K", bound=Hashable)
V = TypeVar("V", float, int)

NumpyArray = npt.NDArray[Any]


class SummaryStatBase(Generic[K, V], ABC, MutableMapping[K, V]):
    @abstractmethod
    def expanded_values(self) -> list[V]: ...

    @abstractmethod
    def __len__(self) -> int: ...

    @abstractmethod
    def items(self) -> ItemsView[K, V]: ...

    @property
    def mean(self) -> numpy.floating | float | int:
        return numpy.mean(self.expanded_values()) if len(self) > 0 else 0

    @property
    def std(self) -> numpy.floating | float | int:
        stat = self.var
        if stat > 0:
            stat = numpy.sqrt(stat)
        return stat

    @property
    def var(self) -> numpy.floating | float | int:
        return numpy.var(self.expanded_values(), ddof=1) if len(self) > 0 else 0

    def quantile(self, q: float) -> numpy.floating | int:
        return numpy.quantile(self.expanded_values(), q=q) if len(self) > 0 else 0

    @property
    def median(self) -> numpy.floating | int:
        return numpy.median(self.expanded_values()) if len(self) > 0 else 0

    @property
    def mode(self) -> K:
        _, mode = max((v, k) for k, v in self.items())
        return mode

    @property
    def sum(self) -> numpy.floating | int:
        return numpy.sum(self.expanded_values()) if len(self) > 0 else 0


class CategoryCounter(SummaryStatBase[K, int]):
    """counting class with summary statistic attributes"""

    def __init__(self, data: dict[K, int] | Iterable[K] | None = None) -> None:
        self._counts: dict[K, int] = {}
        if data is not None:
            if isinstance(data, dict):
                self.update_from_counts(cast("dict[K, int]", data))
            else:
                self.update_from_series(data)

    def update_from_counts(self, data: dict[K, int]) -> None:
        """updates values of self using counts dict"""
        for k, v in data.items():
            self[k] += v

    def update_from_series(self, data: Iterable[K]) -> None:
        """updates values of self from raw series"""
        for element in data:
            self[element] += 1

    def expand(self) -> list[K]:
        """returns list of [[k] * val, ..]"""
        result: list[K] = []
        for k in self:
            result.extend([k] * self[k])
        return result

    def expanded_values(self) -> list[int]:
        return list(self.values())

    def copy(self) -> Self:
        data = self.to_dict().copy()
        return self.__class__(data)

    def __setitem__(self, key: K, val: int) -> None:
        self._counts[key] = val

    def __getitem__(self, key: K) -> int:
        return self._counts.get(key, 0)

    def __delitem__(self, key: K) -> None:
        del self._counts[key]

    def __len__(self) -> int:
        return sum(self.values())

    def __iter__(self) -> Iterator[K]:
        return iter(self._counts)

    def __add__(self, other: K) -> Self:
        self[other] += 1
        return self

    def __sub__(self, other: K) -> Self:
        self[other] -= 1
        if self[other] == 0:
            del self[other]
        return self

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._counts!r})"

    def keys(self) -> KeysView[K]:
        return self._counts.keys()

    def values(self) -> ValuesView[int]:
        return self._counts.values()

    def items(self) -> ItemsView[K, int]:
        return self._counts.items()

    def to_dict(self) -> dict[K, int]:
        return dict(self)

    def tolist(self, keys: Iterable[K] | None = None) -> list[int]:
        """return values for these keys as a list"""
        if keys is None:
            keys = list(self)
        return [self[key] for key in keys]

    def to_array(self, keys: Iterable[K] | None = None) -> npt.NDArray[numpy.integer]:
        """return values for these keys as an array"""
        data = self.tolist(keys=keys)
        return numpy.array(data, dtype=int)

    def to_dictarray(self) -> DictArray:
        """construct fully enumerated dictarray

        Returns
        -------
        DictArray with dtype of int

        Notes
        -----
        Unobserved combinations have zeros. Result can can be indexed as if it was a numpy array using key values
        """
        from itertools import product

        from cogent3.util.dict_array import DictArrayTemplate

        key = next(iter(self))
        try:
            ndim = 1 if isinstance(key, str) else len(key)  # type: ignore[arg-type]
        except TypeError:
            ndim = 1

        if ndim == 1:
            names: list[K] = sorted(self)  # type: ignore[type-var]
            vals = [self[n] for n in names]
            return DictArrayTemplate(cast("list[int]", names)).wrap(vals, dtype=int)

        categories = [sorted(set(labels)) for labels in zip(*self, strict=False)]
        shape = tuple(len(c) for c in categories)
        darr = DictArrayTemplate(*categories).wrap(numpy.zeros(shape, dtype=int))
        for comb in product(*categories):
            indices = [[categories[i].index(c)] for i, c in enumerate(comb)]
            darr.array[tuple(indices)] = self[cast("K", comb)]

        return darr

    def to_categorical(self) -> "CategoryCounts":
        """create CategoryCount object

        Notes
        -----
        Supports only 1 or dimensional data
        """
        from cogent3.maths.stats.contingency import CategoryCounts

        darr = self.to_dictarray()
        return CategoryCounts(darr)

    def to_table(
        self, column_names: str | PySeq[str] | None = None, **kwargs: Any
    ) -> "Table":
        """converts to Table

        Parameters
        ----------
        column_names
            the column name(s) for the key, defaults to "key". If a series, must
            match dimensions of keys, e.g. for (a, b) keys, column_names=['A', 'B']
            will result in a table with 3 columns ('A', 'B', 'count').
        kwargs
            passed to table constructor

        Returns
        -------
        cogent3 Table instance
        """
        from cogent3.core.table import Table

        if (
            not column_names
            or isinstance(column_names, str)
            or not hasattr(column_names, "__len__")
        ):
            key = column_names if column_names is not None else "key"
            data: dict[Any, PySeq[Any] | NumpyArray] = {
                c[0]: c[1:]
                for c in zip([key, "count"], *list(self.items()), strict=False)
            }
            header: list[str] = cast("list[str]", [key, "count"])
            # if keys are tuples, construct the numpy array manually so the
            # elements remain as tuples. numpy's object type casting converts
            # these to lists otherwise
            if isinstance(next(iter(self)), tuple):
                num = len(data[key])
                arr = numpy.empty(num, dtype=object)
                for i in range(num):
                    arr[i] = data[key][i]
                data[key] = arr
            return Table(header=header, data=data, **kwargs)

        key = cast("str", next(iter(self)))
        assert len(key) == len(column_names), "mismatched dimensions"
        d_data: dict[Any, list[Any]] = defaultdict(list)
        for key, count in self.items():  # type: ignore[assignment]
            for c, e in zip(column_names, cast("str", key), strict=False):
                d_data[c].append(e)
            d_data["count"].append(count)
        header = [*list(column_names), "count"]
        data = dict(d_data)
        return Table(header=header, data=data, **kwargs)

    @property
    def entropy(self) -> numpy.floating:
        data = self.to_array()
        data = data / self.sum
        return -(data * numpy.log2(data)).sum()

    def to_freqs(self) -> "CategoryFreqs[K]":
        """returns dict of {key: val/total, ..}"""
        return CategoryFreqs(self, total=float(self.sum))

    def count(self, indices: int | Iterable[int]) -> Self:
        """
        Parameters
        ----------
        indices
            select element(s) from a multi-element tuple keys, must be int or
            series of ints

        Returns
        -------
        CategoryCounter
        """
        if isinstance(indices, int):
            indices = [indices]

        counts: Counter[Any] = Counter()
        for key in self:
            try:
                sub_key = tuple(cast("PySeq[Any]", key)[i] for i in indices)
                sub_key = sub_key[0] if len(sub_key) == 1 else sub_key
            except IndexError:
                msg = f"indices {indices} too big for key {key}"
                raise IndexError(msg)
            counts[sub_key] += self[key]

        return self.__class__(data=counts)


class CategoryFreqs(SummaryStatBase[K, float]):
    """category frequencies with summary statistic attributes"""

    def __init__(
        self,
        data: Mapping[K, float] | None = None,
        total: float | None = None,
        assert_unity: bool = False,
    ) -> None:
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
        self._freqs: dict[K, float] = {}
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

    def expanded_values(self) -> list[float]:
        return list(self.values())

    def copy(self) -> Self:
        data = self.to_dict().copy()
        return self.__class__(data=data)

    def __setitem__(self, key: K, val: float) -> None:
        self._freqs[key] = val

    def __getitem__(self, key: K) -> float:
        return self._freqs.get(key, 0)

    def __delitem__(self, key: K) -> None:
        del self._freqs[key]

    def __len__(self) -> int:
        return len(self._freqs)

    def __iter__(self) -> Iterator[K]:
        return iter(self._freqs)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._freqs!r})"

    def keys(self) -> KeysView[K]:
        return self._freqs.keys()

    def values(self) -> ValuesView[float]:
        return self._freqs.values()

    def items(self) -> ItemsView[K, float]:
        return self._freqs.items()

    def to_dict(self):
        return dict(self)

    def tolist(self, keys: Iterable[K] | None = None) -> list[float]:
        """return values for these keys as a list"""
        if keys is None:
            keys = list(self)
        return [self[key] for key in keys]

    def to_array(self, keys: Iterable[K] | None = None) -> npt.NDArray[numpy.floating]:
        """return just these keys as an array"""
        data = self.tolist(keys=keys)
        return numpy.array(data, dtype=float)

    @property
    def entropy(self) -> numpy.floating:
        data = self.to_array()
        return -(data * numpy.log2(data)).sum()

    def to_normalized(self) -> "CategoryFreqs[K]":
        """returns rescaled self so sum is 1"""
        return CategoryFreqs[K](self, total=float(self.sum), assert_unity=True)


class NumberCounter(CategoryCounter[int]):
    """counts occurrences of numbers"""

    def __init__(self, data=None) -> None:
        super().__init__(data=data)

    @property
    def valid(self):
        types = set(map(type, self))
        if types <= {int, float, complex}:
            result = True
        else:
            key = next(iter(self))
            try:  # if a numpy type
                result = key.dtype.kind in "uifc"  # type: ignore[attr]
            except AttributeError:
                result = False

        return result

    def expanded_values(self, check: bool = False) -> list[int]:
        # TODO memory footprint can be improved by directly computing the
        #  summary statistics
        if check and not self.valid:
            msg = "non-numeric keys"
            raise ValueError(msg)
        values: list[int] = []
        for k, v in self.items():
            values.extend([k] * v)
        return values

    def __len__(self) -> int:
        return sum(self.values())

    @property
    def mean(self) -> float:
        mean = sum(k * self[k] for k in self)
        return mean / len(self)

    @property
    def var(self) -> float:
        """unbiased estimate of the variance"""
        # we scale the variance contribution of a number by its occurrence
        mean = self.mean
        var = sum(self[k] * (k - mean) ** 2 for k in self)
        return var / (len(self) - 1)

    @property
    def std(self) -> numpy.floating:
        return numpy.sqrt(self.var)

    def update_from_counts(self, data: dict[int, int]) -> None:
        """updates values of self using counts dict"""
        for k, v in data.items():
            try:
                k**2
            except TypeError:
                msg = f"key {k} is not numeric"
                raise TypeError(msg)
            self[k] += v
