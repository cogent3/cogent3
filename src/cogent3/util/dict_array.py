"""Wrapper for numpy arrays so that they can be indexed by name

>>> a = numpy.identity(3, int)
>>> b = DictArrayTemplate("abc", "ABC").wrap(a)
>>> b[0]
===========
A    B    C
-----------
1    0    0
-----------
>>> b["a"]
===========
A    B    C
-----------
1    0    0
-----------
>>> b.keys()
['a', 'b', 'c']
>>> b["a"].keys()
['A', 'B', 'C']
"""

import contextlib
import json
from collections import defaultdict
from collections.abc import Hashable, Iterable, Iterator
from collections.abc import Sequence as PySeq
from itertools import combinations, product
from typing import TYPE_CHECKING, Any, Protocol, SupportsInt, TypeVar, cast

import numpy
import numpy.typing as npt
from typing_extensions import Self

from cogent3._version import __version__
from cogent3.util.deserialise import get_class, register_deserialiser
from cogent3.util.io import PathType, atomic_write
from cogent3.util.misc import get_object_provenance

if TYPE_CHECKING:
    from cogent3.core.table import Table


NumpyArray = npt.NDArray[numpy.number]

PySeqStr = PySeq[str]


class SupportsLTHashable(Hashable, Protocol):
    def __lt__(self, other: object) -> bool: ...


T = TypeVar("T")
K = TypeVar("K", bound=SupportsLTHashable)


def convert_1D_dict(
    data: dict[K, list[float]], row_order: PySeq[K] | None = None
) -> tuple[list[list[float]], PySeq[K]]:
    """returns a 1D list and header as dict keys

    Parameters
    ----------
    data : dict
        a 1D dict
    row_order
        series with column headings. If not provided, the sorted top level dict
        keys are used.
    """
    if row_order is None:
        row_order = sorted(data)

    rows = [data[c] for c in row_order]
    return rows, row_order


def convert2Ddistance(
    dists: dict[tuple[K, K], float],
    header: PySeq[K] | None = None,
    row_order: PySeq[K] | None = None,
) -> tuple[list[list[float]], PySeq[K], PySeq[K]]:
    """returns a 2 dimensional list, header and row order

    Parameters
    ----------
    dists : dict
        a 1Ddict with {(a, b): dist, ..}
    header
        series with column headings. If not provided, the sorted top level dict
        keys are used.
    row_order
        a specified order to generate the rows

    Returns
    -------
    2D list, header and row_order. If a dist not present, it's set to 0, or
    the symmetric value e.g. (a, b) -> (b, a).
    """
    if header is None:
        names: set[K] = set()
        for pair in dists:
            names.update(set(pair))
        header = sorted(names)

    rows: list[list[float]] = []
    for i in range(len(header)):
        n1 = header[i]
        row: list[float] = []
        for j in range(len(header)):
            n2 = header[j]
            dist = dists.get((n1, n2), dists.get((n2, n1), 0))
            row.append(dist)
        rows.append(row)

    row_order = header[:]
    return rows, row_order, header


def convert2DDict(
    twoDdict: dict[K, dict[K, float]],
    header: PySeq[K] | None = None,
    row_order: PySeq[K] | None = None,
    make_symmetric: bool = False,
) -> tuple[list[list[float]], PySeq[K], PySeq[K]]:
    """returns a 2 dimensional list, header and row order

    Parameters
    ----------
    twoDdict : dict
        a 2 dimensional dict with top level keys corresponding to column
        headings, lower level keys correspond to row headings
    header
        series with column headings. If not provided, the sorted top level dict
        keys are used.
    row_order
        a specified order to generate the rows
    make_symmetric : bool
        if True, twoDdict[a][b] == twoDdict[b][a]
    """
    if not row_order:
        row_order = list(twoDdict.keys())
        row_order.sort()

    if not header:  # we assume columns consistent across dict
        header = list(twoDdict[row_order[0]].keys())
        header.sort()

    if make_symmetric:
        combined = sorted(set(header) | set(row_order))
        header = row_order = combined
        data: dict[K, dict[K, float]] = defaultdict(dict)

        for k1, k2 in combinations(combined, 2):
            if k1 in twoDdict:
                val = twoDdict[k1].get(k2, 0)
            elif k2 in twoDdict:
                val = twoDdict[k2].get(k1, 0)
            else:
                val = 0
            data[k1][k2] = data[k2][k1] = val
        for k in data:
            data[k][k] = 0
        twoDdict = data

    # make list of lists
    rows: list[list[float]] = []
    for row in row_order:
        elements: list[float] = []
        for column in header:
            elements.append(twoDdict[row][column])
        rows.append(elements)

    return rows, row_order, header


def convert_dict(
    data: dict[K, list[float]] | dict[K, dict[K, float]] | dict[tuple[K, K], float],
    header: PySeq[K] | None = None,
    row_order: PySeq[K] | None = None,
) -> tuple[list[list[float]], PySeq[K], PySeq[K] | None]:
    """returns a list, DictArrayTemplate args

    Parameters
    ----------
    data : dict
        a 1D or 2D dict
    header
        series with column headings. If not provided, the sorted top level dict
        keys are used.
    row_order
        a specified order to generate the rows
    """
    rows: list[list[float]] | list[float]
    first_key = next(iter(data))
    if isinstance(first_key, tuple) and len(first_key) == 2:
        data = cast("dict[tuple[K, K], float]", data)
        rows, row_order, header = convert2Ddistance(data, header, row_order)
    elif hasattr(data[first_key], "keys"):  # type: ignore[index]
        data = cast("dict[K, dict[K, float]]", data)
        rows, row_order, header = convert2DDict(data, header, row_order)
    else:
        data = cast("dict[K, list[float]]", data)
        rows, row_order = convert_1D_dict(data, header)
    return rows, row_order, header


def convert_series(
    data: PySeq[float] | PySeq[PySeq[float]],
    row_order: PySeq[K] | int | None = None,
    header: PySeq[K] | int | None = None,
) -> tuple[PySeq[float] | PySeq[PySeq[float]], PySeq[K] | int, PySeq[K] | int | None]:
    """returns a list, header and row order

    Parameters
    ----------
    data : dict
        a 1D or 2D dict
    header
        series with column headings. If not provided, the sorted top level dict
        keys are used.
    row_order
        a specified order to generate the rows
    """
    first_element = data[0]
    nrows = len(data)
    try:
        ncols = len(first_element)  # type: ignore[arg-type]
    except TypeError:
        ncols = 1

    if header is not None:
        dim_h = header if isinstance(header, int) else len(header)
    else:
        dim_h = None

    if row_order is not None:
        dim_r = row_order if isinstance(row_order, int) else len(row_order)
    else:
        dim_r = None

    if nrows == 1 and ncols > 1:
        if dim_h is not None and dim_h != ncols:
            msg = (
                f"mismatch between number columns={dim_h} "
                f"and number of elements in data={ncols}"
            )
            raise ValueError(
                msg,
            )
        if dim_r is not None and dim_r != 1:
            msg = (
                f"mismatch between number rows={dim_r} "
                f"and number of rows in data={ncols}"
            )
            raise ValueError(
                msg,
            )

    if not header:
        header = None if ncols == 1 else ncols
    row_order = row_order if row_order else nrows

    return data, row_order, header


def convert_for_dictarray(
    data: "DictArray | PySeq[float] | PySeq[PySeq[float]] | dict[K, list[float]] | dict[K, dict[K, float]] | dict[tuple[K, K], float]",
    header: PySeq[K] | None = None,
    row_order: PySeq[K] | None = None,
) -> tuple[
    NumpyArray | PySeq[float] | PySeq[PySeq[float]], PySeq[K] | None, PySeq[K] | None
]:
    """returns a list, header and row order from data

    Parameters
    ----------
    data : iterable
        data series, dictarray, dict, etc..
    header
        series with column headings. If not provided, the sorted top level dict
        keys are used.
    row_order
        a specified order to generate the rows
    """
    result_data: NumpyArray | PySeq[float] | PySeq[PySeq[float]]
    if isinstance(data, DictArray):
        header = data.template.names[0]
        row_order = data.template.names[1]
        result_data = data.array.copy()
    elif hasattr(data, "keys"):  # dictlike, it could be defaultdict
        data = cast(
            "dict[K, list[float]] | dict[K, dict[K, float]] | dict[tuple[K, K], float]",
            data,
        )
        result_data, row_order, header = convert_dict(data, header, row_order)
    else:
        data = cast("PySeq[float] | PySeq[PySeq[float]]", data)
        result_data, row_order, header = convert_series(data, header, row_order)  # type: ignore[assignment]

    return result_data, row_order, header


class NumericKey(int):
    """a distinct numerical type for use as a DictArray key"""

    def __new__(cls, val: SupportsInt) -> "NumericKey":
        return int.__new__(cls, val)


NestedFloatList = list[float] | list["NestedFloatList"]


class DictArrayTemplate:
    def __init__(
        self, *dimensions: Iterable[str] | Iterable[int] | int | None
    ) -> None:  # Consider replacing Iterable[] with Generic specifier
        self.names: list[list[int] | list[NumericKey] | list[str]] = []
        self.ordinals: list[dict[int | NumericKey | str, int]] = []
        for names in dimensions:
            if names is None:
                continue
            if isinstance(names, int):
                names = list(range(names))
            else:
                names = cast(
                    "list[NumericKey] | list[str]",
                    [NumericKey(v) if isinstance(v, int) else v for v in names],
                )

            self.names.append(names)

            ordinals = cast(
                "dict[int | NumericKey | str, int]",
                {c: i for (i, c) in enumerate(names)},
            )

            self.ordinals.append(ordinals)
        self._shape = tuple(len(keys) for keys in self.names)

    def __eq__(self, other: object) -> bool:
        return self is other or (
            isinstance(other, DictArrayTemplate) and self.names == other.names
        )

    def _dict2list(
        self, value: NestedFloatList | dict[Any, NestedFloatList], depth: int = 0
    ) -> NestedFloatList:
        # Unpack (possibly nested) dictionary into correct order of elements
        if depth < len(self._shape):
            value = cast("dict[Any, NestedFloatList]", value)
            return [self._dict2list(value[key], depth + 1) for key in self.names[depth]]
        return cast("NestedFloatList", value)

    def unwrap(
        self,
        value: "DictArray | dict[Any, NestedFloatList] | Iterable[float] | NestedFloatList",
    ) -> NumpyArray:
        """Convert to a simple numpy array"""
        if isinstance(value, DictArray):
            if value.template == self:
                value = value.array
            else:
                raise ValueError  # used to return None, which can't be right
        elif isinstance(value, dict):
            value = self._dict2list(cast("dict[Any, NestedFloatList]", value))
        value = numpy.asarray(value)
        assert value.shape == self._shape, (value.shape, self._shape)
        return value

    def wrap(
        self,
        array: NumpyArray | PySeq[float] | PySeq[PySeq[float]],
        dtype: numpy.dtype[Any] | type[int] | None = None,
    ) -> "DictArray":
        if hasattr(array, "keys"):
            if len(self._shape) == 2:
                r, h = self.names[:2]
            else:
                r, h = self.names[0], None
            array, _, _ = convert_for_dictarray(
                cast(
                    "dict[Any, list[float]] | dict[Any, dict[Any, float]] | dict[tuple[Any, Any], float]",
                    array,
                ),
                h,
                r,
            )
        array = numpy.asarray(array, dtype=dtype)
        for dim, categories in enumerate(self.names):
            assert len(categories) == numpy.shape(array)[dim], (
                f"cats={categories}; dim={dim}"
            )
        return DictArray(array, self)

    def interpret_index(
        self,
        names: NumpyArray
        | tuple[
            int | slice | list[int | NumericKey | str] | NumpyArray | NumericKey | str,
            ...,
        ]
        | int
        | slice
        | list[int | NumericKey | str]
        | NumericKey
        | str,
    ) -> tuple[
        tuple[int | slice | list[int | NumericKey] | NumericKey, ...],
        Self | None,
    ]:
        if isinstance(names, numpy.ndarray) and "int" in names.dtype.name:
            # the numpy item() method casts to the nearest Python type
            names = tuple(v.item() for v in names)

        if not isinstance(names, tuple):
            names = (names,)

        index: list[int | slice | list[int | NumericKey] | NumericKey] = []
        remaining: list[list[NumericKey] | list[int] | list[str]] = []
        for ordinals, allnames, name in zip(
            self.ordinals,
            self.names,
            names,
            strict=False,
        ):
            if not isinstance(name, (int, slice, list, numpy.ndarray)):
                name = ordinals[name]
            elif isinstance(name, slice):
                start = name.start
                stop = name.stop
                with contextlib.suppress(ValueError):
                    # either None, or it's an int index
                    start = allnames.index(start)
                with contextlib.suppress(ValueError):
                    # as above
                    stop = allnames.index(stop)
                name = slice(start, stop, name.step)
                remaining.append(allnames.__getitem__(name))
            elif isinstance(name, (list, numpy.ndarray)):
                name = [n if isinstance(n, int) else ordinals[n] for n in name]
                remaining.append(
                    cast(
                        "list[NumericKey] | list[int] | list[str]",
                        [allnames[cast("int", i)] for i in name],
                    )
                )

            index.append(
                cast(
                    "int | slice | list[int | NumericKey] | NumericKey",
                    name,
                )
            )
        remaining.extend(self.names[len(index) :])
        klass = type(self)(*remaining) if remaining else None
        return (tuple(index), klass)


class DictArray:
    """Wraps a numpy array so that it can be indexed with strings. Behaves
    like nested dictionaries (only ordered).

    Notes
    -----
    Used for things like substitution matrices and bin probabilities.

    Indexing can be done via conventional integer based operations, using
    keys, lists of int/keys.

    Behaviour differs from numpy array indexing when you provide lists of
    indices. Such indexing is applied sequentially, e.g. darr[[0, 2], [1, 2]]
    will return the intersection of rows [0, 2] with columns [1, 2]. In numpy,
    the result would instead be the elements at [0, 1], [2, 2].
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """allow alternate ways of creating for time being"""
        if len(args) == 1:
            vals, row_keys, col_keys = convert_for_dictarray(args[0])
            dtype = kwargs.get("dtype")
            self.array: npt.NDArray[numpy.number] = numpy.asarray(vals, dtype=dtype)
            self.template = DictArrayTemplate(row_keys, col_keys)
        elif len(args) == 2:
            if not isinstance(args[1], DictArrayTemplate):
                raise NotImplementedError
            self.array = args[0]
            self.template = args[1]
        else:
            if "dtype" in kwargs or "typecode" in kwargs:
                dtype = kwargs["dtype"]
                kwargs.pop("dtype", None)
                kwargs.pop("typecode", None)
            else:
                dtype = None
            create_new = DictArrayTemplate(*args[1:]).wrap(args[0], dtype=dtype)
            self.__dict__ = create_new.__dict__
        self.shape = self.array.shape

    @classmethod
    def from_array_names(
        cls,
        array: npt.NDArray[numpy.number],
        *names: PySeqStr,
    ) -> "DictArray":
        """creates instance directly from a numpy array and series of names

        Parameters
        ----------
        array
            any data type
        names
            must match the array dimensions
        """

        if len(names) != array.ndim or any(
            len(labels) != array.shape[dim] for dim, labels in enumerate(names)
        ):
            msg = "names must match array dimensions"
            raise ValueError(msg)

        template = DictArrayTemplate(*names)
        return DictArray(array, template)

    def to_array(self) -> NumpyArray:
        return self.array

    def __add__(self, other: Self) -> "DictArray":
        if not isinstance(other, type(self)):
            msg = f"Incompatible types: {type(self)} and {type(other)}"
            raise TypeError(msg)

        if other.template.names != self.template.names:
            msg = f"unequal dimension names {self.template.names} != {other.template.names}"
            raise ValueError(
                msg,
            )

        return self.template.wrap(self.array + other.array)

    def __array__(
        self,
        dtype: numpy.dtype[Any] | None = None,
        copy: bool = False,
    ) -> NumpyArray:
        array = self.array
        if dtype is not None:
            array = array.astype(dtype)
        return array

    def to_dict(
        self, flatten: bool = False
    ) -> dict[int | NumericKey | str | tuple[int | NumericKey | str, ...], Any]:
        """returns data as a dict

        Parameters
        ----------
        flatten : bool
            returns a 1D dictionary
        """
        names = self.template.names
        shape = self.shape
        result: dict[
            int | NumericKey | str | tuple[int | NumericKey | str, ...], Any
        ] = {}
        if len(names) == 1:
            result = {
                names[0][i]: v.item() if hasattr(v, "item") else v
                for i, v in enumerate(self.array)
            }
        elif flatten:
            for indices in product(*[range(n) for n in shape]):
                value = self.array[indices]
                value = value.item() if hasattr(value, "item") else value
                coord = tuple(n[i] for n, i in zip(names, indices, strict=False))
                result[coord] = value
        else:
            for indices in product(*[range(n) for n in shape]):
                value = self.array[indices]
                value = value.item() if hasattr(value, "item") else value
                coord = tuple(n[i] for n, i in zip(names, indices, strict=False))
                current = result
                nested = coord[0]
                for nested in coord[:-1]:
                    current[nested] = current.get(nested, {})
                current[nested][coord[-1]] = value

        return result

    def to_rich_dict(self) -> dict[str, Any]:
        data = self.array.tolist()
        return {
            "type": get_object_provenance(self.template),
            "array": data,
            "names": self.template.names,
            "version": __version__,
        }

    def to_json(self) -> str:
        return json.dumps(self.to_rich_dict())

    def __getitem__(
        self,
        names: NumpyArray
        | tuple[
            int | slice | list[int | NumericKey | str] | NumpyArray | NumericKey | str,
            ...,
        ]
        | int
        | slice
        | list[int | NumericKey | str]
        | NumericKey
        | str,
    ) -> "NumpyArray | DictArray":
        index, remaining = self.template.interpret_index(names)
        if list in {type(v) for v in index}:
            result = self.array
            for dim, indices in enumerate(index):
                if isinstance(indices, slice):
                    actual_indices = (
                        (indices,)
                        if dim == 0
                        else (slice(None, None),) * dim + (indices,)
                    )
                    result = result[tuple(actual_indices)]
                    continue

                if isinstance(indices, int):
                    indices = [indices]

                result = result.take(indices, axis=dim)

        else:
            result = self.array[index]

        if remaining is None:
            return result
        return self.__class__(result.reshape(remaining._shape), remaining)

    def __iter__(self) -> Iterator["float | DictArray"]:
        _, remaining = self.template.interpret_index(0)
        for elt in self.array:
            if remaining is None:
                yield elt
            else:
                yield remaining.wrap(elt)

    def __len__(self) -> int:
        return len(self.template.names[0])

    def keys(self) -> list[int] | list[str] | list[NumericKey]:
        return self.template.names[0][:]

    def items(self) -> list[tuple[int | str | NumericKey, Any]]:
        return [(n, self[n]) for n in self.keys()]

    def __repr__(self) -> str:
        if self.array.ndim > 2:
            return f"{self.array.ndim} dimensional {type(self).__name__}"

        t = self.to_table()
        t.set_repr_policy(show_shape=False)
        return str(t)

    def __ne__(self, other: object) -> bool:
        return not self.__eq__(other)

    def __eq__(self, other: object) -> bool:
        if self is other:
            return True
        if isinstance(other, DictArray):
            return bool(
                self.template == other.template
                and numpy.all(
                    self.array == other.array,
                )
            )
        if isinstance(other, type(self.array)):
            return self.array == other
        if isinstance(other, dict):
            return self.to_dict() == other
        return False

    def to_normalized(
        self, by_row: bool = False, by_column: bool = False
    ) -> "DictArray":
        """returns a DictArray as frequencies

        Parameters
        ----------
        by_row
            rows sum to 1
        by_col
            columns sum to 1
        """
        assert not (by_row and by_column)
        # TODO need to check there are two dimension!
        if by_row:
            axis = 1
        elif by_column:
            axis = 0
        else:
            axis = None

        result = self.array / self.array.sum(axis=axis, keepdims=True)
        return self.template.wrap(result)

    def row_sum(self) -> "DictArray":
        """returns DictArray summed across rows"""
        axis = 1 if len(self.shape) == 2 else 0
        result = self.array.sum(axis=axis)
        template = DictArrayTemplate(self.template.names[0])
        return template.wrap(result)

    def col_sum(self) -> "DictArray":
        """returns DictArray summed across columns"""
        result = self.array.sum(axis=0)
        template = DictArrayTemplate(self.template.names[1])
        return template.wrap(result)

    def _repr_html_(self) -> str:
        if self.array.ndim > 2:
            return f"{self.array.ndim} dimensional {type(self).__name__}"

        t = self.to_table()
        t.set_repr_policy(show_shape=False)
        return t._repr_html_()

    def to_string(self, format_name: str = "tsv", sep: str | None = None) -> str:
        """Return the data as a formatted string.

        Parameters
        ----------
        format_name
            possible formats are 'csv', or 'tsv' (default).
        sep
            A string separator for delineating columns, e.g. ',' or
            '\t'. Overrides format.
        """
        if format_name.lower() not in ("tsv", "csv"):
            msg = f"'{format_name}' not supported"
            raise ValueError(msg)

        sep = sep or {"tsv": "\t", "csv": ","}[format_name.lower()]

        data = cast(
            "dict[tuple[int | NumericKey | str, ...], Any]",
            self.to_dict(flatten=True),
        )
        rows = [[f"dim-{i + 1}" for i in range(self.array.ndim)] + ["value"]] + [
            [str(x) for x in row] for row in [[*list(k), v] for k, v in data.items()]
        ]
        return "\n".join([sep.join(row) for row in rows])

    def to_table(self) -> "Table":
        """return Table instance

        Notes
        -----
        Raises ValueError if number of dimensions > 2
        """
        ndim = self.array.ndim
        if ndim > 2:
            msg = f"cannot make 2D table from {ndim}D array"
            raise ValueError(msg)

        from cogent3.core.table import Table

        header = self.template.names[0] if ndim == 1 else self.template.names[1]
        index = "" if ndim == 2 else None
        if ndim == 1:
            data = {c: [v] for c, v in zip(header, self.array, strict=False)}
        else:
            data = {c: self.array[:, i].tolist() for i, c in enumerate(header)}
            data[""] = self.template.names[0]

        return Table(
            header=cast("list[str]", header),
            data=cast("list[list[Any]]", data),  # FIXME: Incorrect cast
            index_name=index,
        )

    def write(self, path: PathType, format_name: str = "tsv", sep: str = "\t") -> None:
        """writes a flattened version to path

        Parameters
        ----------
        path
        format_name
            possible formats are 'rest'/'rst', 'markdown'/'md',
            'latex', 'html', 'phylip', 'bedgraph', 'csv', 'tsv', or 'simple'
            (default).
        sep
            used to split fields, will be inferred from path suffix if not
            provided
        """
        data = self.to_string(format_name=format_name, sep=sep)
        with atomic_write(path, mode="wt") as outfile:
            outfile.write(data)


@register_deserialiser(
    get_object_provenance(DictArrayTemplate),
)
def deserialise_dict_array(data: dict[str, Any]) -> DictArray:
    """deserialising DictArray, Table instances"""
    data.pop("version", None)
    type_ = data.pop("type")
    klass = get_class(type_)
    named_dims = data.pop("names")
    array = data.pop("array")
    template = klass(*named_dims)
    return template.wrap(array)
