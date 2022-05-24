"""Wrapper for numpy arrays so that they can be indexed by name

    >>> a = numpy.identity(3, int)
    >>> b = DictArrayTemplate('abc', 'ABC').wrap(a)
    >>> b[0]
    ===========
    A    B    C
    -----------
    1    0    0
    -----------
    >>> b['a']
    ===========
    A    B    C
    -----------
    1    0    0
    -----------
    >>> b.keys()
    ['a', 'b', 'c']
    >>> b['a'].keys()
    ['A', 'B', 'C']
"""
import json

from collections import defaultdict
from itertools import combinations, product

import numpy

from cogent3.util.io import atomic_write
from cogent3.util.misc import get_object_provenance


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Ben Kaehler"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def convert_1D_dict(data, row_order=None):
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
        row_order = list(sorted(data))

    rows = [data[c] for c in row_order]
    return rows, row_order


def convert2Ddistance(dists, header=None, row_order=None):
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
        names = set()
        for pair in dists:
            names.update(set(pair))
        header = list(sorted(names))

    rows = []
    for i in range(len(header)):
        n1 = header[i]
        row = []
        for j in range(len(header)):
            n2 = header[j]
            dist = dists.get((n1, n2), dists.get((n2, n1), 0))
            row.append(dist)
        rows.append(row)

    row_order = header[:]
    return rows, row_order, header


def convert2DDict(twoDdict, header=None, row_order=None, make_symmetric=False):
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
        combined = list(sorted(set(header) | set(row_order)))
        header = row_order = combined
        data = defaultdict(dict)

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
    rows = []
    for row in row_order:
        elements = []
        for column in header:
            elements.append(twoDdict[row][column])
        rows.append(elements)

    return rows, row_order, header


def convert_dict(data, header=None, row_order=None):
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
    first_key = list(data)[0]
    if type(first_key) == tuple and len(first_key) == 2:
        rows, row_order, header = convert2Ddistance(data, header, row_order)
    elif hasattr(data[first_key], "keys"):
        rows, row_order, header = convert2DDict(data, header, row_order)
    else:
        rows, row_order = convert_1D_dict(data, header)
    return rows, row_order, header


def convert_series(data, row_order=None, header=None):
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
        ncols = len(first_element)
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
            raise ValueError(
                f"mismatch between number columns={dim_h} "
                f"and number of elements in data={ncols}"
            )
        elif dim_r is not None and dim_r != 1:
            raise ValueError(
                f"mismatch between number rows={dim_r} "
                f"and number of rows in data={ncols}"
            )

    if not header:
        header = None if ncols == 1 else ncols
    row_order = row_order if row_order else nrows

    return data, row_order, header


def convert_for_dictarray(data, header=None, row_order=None):
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
    if isinstance(data, DictArray):
        header = data.template.names[0]
        row_order = data.template.names[1]
        data = data.array.copy()
    elif hasattr(data, "keys"):  # dictlike, it could be defaultdict
        data, row_order, header = convert_dict(data, header, row_order)
    else:
        data, row_order, header = convert_series(data, header, row_order)

    return data, row_order, header


class NumericKey(int):
    """a distinct numerical type for use as a DictArray key"""

    def __new__(cls, val):
        return int.__new__(cls, val)


class DictArrayTemplate(object):
    def __init__(self, *dimensions):
        self.names = []
        self.ordinals = []
        for names in dimensions:
            if names is None:
                continue
            elif isinstance(names, int):
                names = list(range(names))
            else:
                names = [NumericKey(v) if type(v) == int else v for v in names]

            self.names.append(names)
            self.ordinals.append(dict((c, i) for (i, c) in enumerate(names)))
        self._shape = tuple(len(keys) for keys in self.names)

    def __eq__(self, other):
        return self is other or (
            isinstance(other, DictArrayTemplate) and self.names == other.names
        )

    def _dict2list(self, value, depth=0):
        # Unpack (possibly nested) dictionary into correct order of elements
        if depth < len(self._shape):
            return [self._dict2list(value[key], depth + 1) for key in self.names[depth]]
        else:
            return value

    def unwrap(self, value):
        """Convert to a simple numpy array"""
        if isinstance(value, DictArray):
            if value.template == self:
                value = value.array
            else:
                raise ValueError  # used to return None, which can't be right
        elif isinstance(value, dict):
            value = self._dict2list(value)
        value = numpy.asarray(value)
        assert value.shape == self._shape, (value.shape, self._shape)
        return value

    def wrap(self, array, dtype=None):
        if hasattr(array, "keys"):
            if len(self._shape) == 2:
                r, h = self.names[:2]
            else:
                r, h = self.names[0], None
            array, _, _ = convert_for_dictarray(array, h, r)
        array = numpy.asarray(array, dtype=dtype)
        for (dim, categories) in enumerate(self.names):
            assert len(categories) == numpy.shape(array)[dim], "cats=%s; dim=%s" % (
                categories,
                dim,
            )
        return DictArray(array, self)

    def interpret_index(self, names):
        if isinstance(names, numpy.ndarray) and "int" in names.dtype.name:
            # the numpy item() method casts to the nearest Python type
            names = tuple(v.item() for v in names)

        if not isinstance(names, tuple):
            names = (names,)

        index = []
        remaining = []
        for (ordinals, allnames, name) in zip(self.ordinals, self.names, names):
            if type(name) not in (int, slice, list, numpy.ndarray):
                name = ordinals[name]
            elif isinstance(name, slice):
                start = name.start
                stop = name.stop
                try:
                    start = allnames.index(start)
                except ValueError:
                    # either None, or it's an int index
                    pass
                try:
                    stop = allnames.index(stop)
                except ValueError:
                    # as above
                    pass
                name = slice(start, stop, name.step)
                remaining.append(allnames.__getitem__(name))
            elif type(name) in (list, numpy.ndarray):
                name = [n if type(n) == int else ordinals[n] for n in name]
                remaining.append([allnames[i] for i in name])

            index.append(name)
        remaining.extend(self.names[len(index) :])
        klass = type(self)(*remaining) if remaining else None
        return (tuple(index), klass)


class DictArray(object):
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

    def __init__(self, *args, **kwargs):
        """allow alternate ways of creating for time being"""
        if len(args) == 1:
            vals, row_keys, col_keys = convert_for_dictarray(args[0])
            dtype = kwargs.get("dtype", None)
            self.array = numpy.asarray(vals, dtype=dtype)
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

    def to_array(self):
        return self.array

    def __add__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError(f"Incompatible types: {type(self)} and {type(other)}")

        if other.template.names != self.template.names:
            raise ValueError(
                f"unequal dimension names {self.template.names} != {other.template.names}"
            )

        return self.template.wrap(self.array + other.array)

    def __array__(self, dtype=None):
        array = self.array
        if dtype is not None:
            array = array.astype(dtype)
        return array

    def to_dict(self, flatten=False):
        """returns data as a dict

        Parameters
        ----------
        flatten : bool
            returns a 1D dictionary
        """
        names = self.template.names
        shape = self.shape
        result = {}
        if len(names) == 1:
            result = {
                names[0][i]: v.item() if hasattr(v, "item") else v
                for i, v in enumerate(self.array)
            }
        elif flatten:
            for indices in product(*[range(n) for n in shape]):
                value = self.array[indices]
                value = value.item() if hasattr(value, "item") else value
                coord = tuple(n[i] for n, i in zip(names, indices))
                result[coord] = value
        else:
            for indices in product(*[range(n) for n in shape]):
                value = self.array[indices]
                value = value.item() if hasattr(value, "item") else value
                coord = tuple(n[i] for n, i in zip(names, indices))
                current = result
                nested = coord[0]
                for nested in coord[:-1]:
                    current[nested] = current.get(nested, {})
                current[nested][coord[-1]] = value

        return result

    def to_rich_dict(self):
        data = self.array.tolist()
        return {
            "type": get_object_provenance(self.template),
            "array": data,
            "names": self.template.names,
            "version": __version__,
        }

    def to_json(self):
        return json.dumps(self.to_rich_dict())

    def __getitem__(self, names):
        (index, remaining) = self.template.interpret_index(names)
        if list in {type(v) for v in index}:
            result = self.array
            for dim, indices in enumerate(index):
                if isinstance(indices, slice):
                    indices = (
                        (indices,)
                        if dim == 0
                        else (slice(None, None),) * dim + (indices,)
                    )
                    result = result[tuple(indices)]
                    continue

                if isinstance(indices, int):
                    indices = [indices]

                result = result.take(indices, axis=dim)

        else:
            result = self.array[index]

        if remaining is not None:
            result = self.__class__(result.reshape(remaining._shape), remaining)
        return result

    def __iter__(self):
        (index, remaining) = self.template.interpret_index(0)
        for elt in self.array:
            if remaining is None:
                yield elt
            else:
                yield remaining.wrap(elt)

    def __len__(self):
        return len(self.template.names[0])

    def keys(self):
        return self.template.names[0][:]

    def items(self):
        return [(n, self[n]) for n in list(self.keys())]

    def __repr__(self):
        if self.array.ndim > 2:
            return f"{self.array.ndim} dimensional {type(self).__name__}"

        t = self.to_table()
        t.set_repr_policy(show_shape=False)
        return str(t)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        if self is other:
            return True
        elif isinstance(other, DictArray):
            return self.template == other.template and numpy.all(
                self.array == other.array
            )
        elif isinstance(other, type(self.array)):
            return self.array == other
        elif isinstance(other, dict):
            return self.to_dict() == other
        else:
            return False

    def to_normalized(self, by_row=False, by_column=False):
        """returns a DictArray as frequencies

        Parameters
        ----------
        by_row
            rows sum to 1
        by_col
            columns sum to 1
        """
        assert not (by_row and by_column)
        # todo need to check there are two dimension!
        if by_row:
            axis = 1
        elif by_column:
            axis = 0
        else:
            axis = None

        result = self.array / self.array.sum(axis=axis, keepdims=True)
        return self.template.wrap(result)

    def row_sum(self):
        """returns DictArray summed across rows"""
        axis = 1 if len(self.shape) == 2 else 0
        result = self.array.sum(axis=axis)
        template = DictArrayTemplate(self.template.names[0])
        return template.wrap(result)

    def col_sum(self):
        """returns DictArray summed across columns"""
        result = self.array.sum(axis=0)
        template = DictArrayTemplate(self.template.names[1])
        return template.wrap(result)

    def _repr_html_(self):
        if self.array.ndim > 2:
            return f"{self.array.ndim} dimensional {type(self).__name__}"

        t = self.to_table()
        t.set_repr_policy(show_shape=False)
        return t._repr_html_()

    def to_string(self, format="tsv", sep=None):
        """Return the data as a formatted string.

        Parameters
        ----------
        format
            possible formats are 'csv', or 'tsv' (default).
        sep
            A string separator for delineating columns, e.g. ',' or
            '\t'. Overrides format.
        """
        if format.lower() not in ("tsv", "csv"):
            msg = f"'{format}' not supported"
            raise ValueError(msg)

        sep = sep or {"tsv": "\t", "csv": ","}[format.lower()]

        data = self.to_dict(flatten=True)
        rows = [[f"dim-{i + 1}" for i in range(self.array.ndim)] + ["value"]] + [
            list(map(lambda x: str(x), row))
            for row in [list(k) + [v] for k, v in data.items()]
        ]
        return "\n".join([sep.join(row) for row in rows])

    def to_table(self):
        """return Table instance

        Notes
        -----
        Raises ValueError if number of dimensions > 2
        """
        ndim = self.array.ndim
        if ndim > 2:
            raise ValueError(f"cannot make 2D table from {ndim}D array")

        from .table import Table

        header = self.template.names[0] if ndim == 1 else self.template.names[1]
        index = "" if ndim == 2 else None
        if ndim == 1:
            data = {c: [v] for c, v in zip(header, self.array)}
        else:
            data = {c: self.array[:, i].tolist() for i, c in enumerate(header)}
            data[""] = self.template.names[0]

        return Table(header=header, data=data, index_name=index)

    def write(self, path, format="tsv", sep="\t"):
        """writes a flattened version to path

        Parameters
        ----------
        path : str
        format
            possible formats are 'rest'/'rst', 'markdown'/'md',
            'latex', 'html', 'phylip', 'bedgraph', 'csv', 'tsv', or 'simple'
            (default).
        sep : str
            used to split fields, will be inferred from path suffix if not
            provided
        """
        data = self.to_string(format=format, sep=sep)
        with atomic_write(path, mode="wt") as outfile:
            outfile.write(data)
