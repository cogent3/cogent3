"""
A light-weight Table class for manipulating 2D data and representing it as
text, or writing to file for import into other packages.

Current output formats include pickle (pythons serialisation format),
restructured text (keyed by 'rest'), latex, html, delimited columns, and a
simple text format.

Table can read pickled and delimited formats.
"""

import csv
import json
import pathlib
import pickle
import re

from collections import defaultdict
from collections.abc import Callable, MutableMapping
from itertools import product
from xml.sax.saxutils import escape

import numpy

from cogent3.format import bedgraph
from cogent3.format import table as table_format
from cogent3.util.dict_array import DictArray, DictArrayTemplate
from cogent3.util.io import atomic_write, get_format_suffixes
from cogent3.util.misc import extend_docstring_from, get_object_provenance
from cogent3.util.union_dict import UnionDict


try:
    from IPython.display import display
except ImportError:
    display = lambda x: print(repr(x))

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Felix Schill", "Sheng Koh"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

# making reversed characters for use in reverse order sorting
_all_chrs = [chr(i) for i in range(256)]
_all_chrs.reverse()
_reversed_chrs = "".join(_all_chrs)


def _reverse_str(x):
    """returns reverse translation of x"""
    return x.translate(_reversed_chrs)


def _reverse_num(x):
    """returns reversed val of x"""
    return x * -1


def _numeric_sum(data):
    """returns sum of all numeric values"""
    try:
        result = numpy.sum(data)
        result / 3
        return result
    except TypeError:
        pass

    total = 0
    valid = False
    for v in data:
        try:
            total += v
            valid = True
        except TypeError:
            pass
    result = total if valid else numpy.nan
    return result


def _callback(callback, row, num_columns=None):
    if isinstance(callback, Callable):
        if num_columns == 1:
            row = row[0]
        return callback(row)
    else:
        return eval(callback, {}, row)


_num_type = re.compile("^(float|int|complex)").search


def array_is_num_type(data):
    """whether data has a dtype for int, float or complex"""
    return _num_type(data.dtype.name) != None


def cast_str_to_numeric(values):
    """converts a series of strings to numeric values"""
    if not (isinstance(values[0], str) or isinstance(values[0], bytes)):
        return numpy.array(values)

    if not isinstance(values, numpy.ndarray):
        values = numpy.array(values, dtype="U")

    for typ in (int, float, complex):
        try:
            values = values.astype(typ)
            break
        except (ValueError, TypeError):
            pass
    return values


def cast_str_to_array(values, static_type=False):
    """converts a series of strings to numeric values"""
    values = numpy.array(values, dtype="U")
    result = cast_str_to_numeric(values)
    if static_type or result is not values:
        return result

    # we handle mixed types by using eval
    result = []
    all_fail = True
    for v in values:
        try:
            v = eval(v)
            all_fail = False
        except (TypeError, NameError, SyntaxError):
            # syntax error from empty strings
            pass
        result.append(v)

    if not all_fail:
        result = numpy.array(result, dtype="O")
    else:
        result = values

    return result


_numeric_types = {int, float, complex}


def cast_to_array(values):
    """converts a series to a general array type"""
    if isinstance(values, numpy.ndarray):
        return values

    types = {type(v) for v in values}
    if len(types) != 1 and not types <= _numeric_types:
        return numpy.array(values, dtype=object)

    # force unicode if str type, otherwise try None
    dtype = "U" if types == {str} else None
    try:
        result = numpy.array(values, dtype=dtype)
    except Exception:
        result = numpy.array(values, dtype=object)

    return result


def cast_2d_to_1d_dict(data, row_order=None):
    """converts a 2D dict to a 1D dict"""
    if not row_order:
        key = list(data.keys())[0]
        row_order = list(data[key])

    return {c: [data[c][r] for r in row_order] for c in data}


def cast_to_1d_dict(data, row_order=None):
    """Returns a 2 dimensional list.

    Parameters
    ----------
    data : dict
        may be 2D
    row_order
        a specified order to generate the rows.
    """
    val_types = {type(v): v for v in data.values()}
    if dict in val_types:
        result = cast_2d_to_1d_dict(data, row_order=row_order)
    else:
        result = data
    return result


class Columns(MutableMapping):
    """Collection of columns. iter operates over columns."""

    def __init__(self):
        self._order = ()
        self._num_rows = 0
        self._template = None
        self._index_name = None

    def _get_key_(self, value):
        """returns string corresponding to column"""

        if isinstance(value, int):
            try:
                value = self._order[value]
            except IndexError:
                raise KeyError(f"no key corresponding to index {value}")

        return value

    def _get_keys_(self, key):
        """returns series of str corresponding to columns"""
        if isinstance(key, str) or isinstance(key, int):
            key = self._get_key_(key)
            return key

        if isinstance(key, slice):
            key, _ = self._template.interpret_index(key)
            key = self._order[key[0]]

        if type(key) in (list, tuple):
            if all(type(e) == bool for e in key) and len(key) == len(self.order):
                key = [k for k, b in zip(self.order, key) if b]
            else:
                key = [self._get_key_(k) for k in key]

            return key

        if isinstance(key, numpy.ndarray):
            # we try slicing by array
            cols = numpy.array(self.order, dtype="U")
            try:
                key = cols[key]
            except Exception:
                msg = f"{key} could not be used to slice columns"
                raise KeyError(msg)

            return key

        raise KeyError(f"{key}")

    def __contains__(self, key):
        return key in self._order

    def __getitem__(self, key):
        if isinstance(key, str) or isinstance(key, int):
            key = self._get_key_(key)
            return self.__dict__[key]

        if isinstance(key, slice):
            key, _ = self._template.interpret_index(key)
            key = self._order[key[0]]

        if isinstance(key, numpy.ndarray):
            key = numpy.array(self._order)[key].tolist()

        if type(key) in (list, tuple):
            result = [self.__dict__[self._get_key_(k)] for k in key]
        else:
            raise KeyError(f"{key}")

        return result

    def __delitem__(self, key):
        key = self._get_key_(key)
        del self.__dict__[key]
        self._order = tuple(k for k in self._order if k != key)
        self._template = DictArrayTemplate(self._order)

    def __iter__(self):
        return iter(k for k in self._order)

    def __len__(self):
        return len(self._order)

    def __setitem__(self, key, val):
        key = str(key).strip()
        if isinstance(val, str):
            val = [val]
        try:
            _ = len(val)
        except TypeError:
            val = [val]

        if self._num_rows == 0:
            self._num_rows = len(val)
        elif len(val) != self._num_rows:
            raise ValueError("number rows incorrect")

        if key not in self._order:
            self._order += (key,)
            self._template = DictArrayTemplate(self._order)

        if not isinstance(val, numpy.ndarray):
            val = cast_to_array(val)

        # make immutable, sort of
        val.flags.writeable = False

        self.__dict__[key] = val

    def __getstate__(self):
        # note that index name is captured by the Table
        result = {"order": list(self.order), "columns": {}}
        for c in self:
            v = self[c]
            dtype = v.dtype.name
            if dtype.startswith("str"):
                dtype = dtype.replace("str", "U")
            result["columns"][c] = dict(values=v.tolist(), dtype=dtype)
        return result

    def __setstate__(self, data):
        new = self.__class__()
        for k in ("type", "version"):
            data.pop(k, None)

        order = data.pop("order")
        columns = data.pop("columns")
        for c in order:
            values, dtype = columns[c]["values"], columns[c]["dtype"]
            values = numpy.array(values, dtype=dtype)
            new[c] = values

        self.__dict__.update(new.__dict__)

    def __repr__(self):
        d = [f"'{c}': {v.dtype}" for c, v in self.items()]
        num = len(d)
        v = d[:5]
        if num > 5:
            v.append(f"... + {num - 5} more")
        return f"{self.__class__.__name__}({', '.join(v)})"

    def __str__(self):
        return repr(self)

    def iter_rows(self):
        columns = [self[c] for c in self]
        for row in zip(*columns):
            yield self._template.wrap(row, dtype=object)

    @property
    def index_name(self):
        """column name whose values can be used to index table rows"""
        return self._index_name

    @index_name.setter
    def index_name(self, name):
        if name is None:
            self._index_name = None
            return

        if name not in self:
            raise ValueError(f"'{name}' unknown, index_name must be an existing column")

        # make sure index_name has unique values
        unique = set(self[name])
        if len(unique) != self._num_rows:
            raise ValueError(
                f"cannot use '{name}' as index_name, not all values unique"
            )

        self._index_name = name
        order = [name] + [c for c in self._order if c != name]
        self._order = tuple(order)

    def add_column_from_str(self, name, values):
        """adds a column from series of str

        Parameters
        ----------
        name : str
            column name
        values : series
            any type, cast to numpy array
        """
        values = cast_str_to_numeric(values)
        self[name] = values

    def take_columns(self, columns):
        """returns new Columns instance with just columns"""
        result = self.__class__()
        if type(columns) in {int, str}:
            columns = [columns]

        columns = self._get_keys_(columns)

        for c in columns:
            result[c] = self[c]

        return result

    @property
    def array(self):
        """object array of all columns"""
        arr = numpy.empty((len(self), self._num_rows), dtype="O")
        for i, c in enumerate(self.order):
            try:
                arr[i] = self[c]
            except ValueError:
                # this can happen of elements of array are tuples, for example
                v = numpy.empty(self._num_rows, dtype="O")
                for j, e in enumerate(self[c]):
                    v[j] = e
                arr[i] = v

        return arr.T

    @property
    def order(self):
        """column order"""
        # if index_name not first, we re-order
        if self._index_name is not None and self._order[0] != self._index_name:
            order = [self._index_name] + [
                c for c in self._order if c != self._index_name
            ]
            self._order = tuple(order)
        return self._order

    def to_dict(self):
        """returns column based dict"""
        return {c: self[c].tolist() for c in self}

    def to_rich_dict(self):
        data = self.__getstate__()
        data["type"] = get_object_provenance(self)
        data["version"] = None  # todo
        return data


class Table:
    """Tabular data. iter operates over rows. Columns are available as an attribute."""

    def __init__(
        self,
        header=None,
        data=None,
        index_name=None,
        title="",
        legend="",
        digits=4,
        space=4,
        max_width=1e100,
        column_templates=None,
        format="simple",
        missing_data="",
        **kwargs,
    ):
        """

        Parameters
        ----------
        header
            column headings
        data
            a 2D dict, list or tuple. If a dict, it must have column
            headings as top level keys, and common row labels as keys in each
            column.
        index_name
            column name with values to be used as row identifiers and keys
            for slicing. All column values must be unique.
        legend
            table legend
        title
            as implied
        digits
            floating point resolution
        space
            number of spaces between columns or a string
        max_width
            maximum column width for printing
        column_templates
            dict of column headings
            or a function that will handle the formatting.
        format
            output format when using str(Table)
        missing_data
            replace missing data with this
        """
        attrs = {
            k: v
            for k, v in locals().items()
            if k not in ("self", "__class__", "data", "header", "kwargs")
        }

        attrs.update(kwargs)

        self._persistent_attrs = attrs

        self.columns = Columns()
        self._template = None
        self._index_name = None

        if isinstance(data, dict):
            # convert containers like a defaultdict to a standard dict
            data = dict(data)
            row_data = False
        else:
            try:
                len(data[0])
                row_data = True
            except (TypeError, IndexError, KeyError):
                row_data = False

        if header and row_data:
            hlen = len(header)
            dcols = len(data[0])
            if hlen != dcols:
                raise ValueError(
                    f"different number of elements in header {hlen} and data row 0 {dcols}"
                )

            data = {c: v for c, v in zip(header, zip(*data))}

        if header is None:
            header = list(data) if isinstance(data, dict) else []
        has_index = index_name is not None
        if has_index and not isinstance(index_name, str):
            raise TypeError(
                f"only str type supported for index_name, not {type(index_name)}"
            )

        if len(data) if hasattr(data, "__len__") else 0:
            row_order = kwargs.get("row_order", None)
            data = cast_to_1d_dict(data, row_order=row_order)
            if has_index:
                try:
                    self.columns[index_name] = data[index_name]
                except KeyError:
                    raise ValueError(f"'{index_name}' not in data")

            for c in header:
                if c == index_name:
                    continue
                self.columns[c] = data[c]

        elif header:
            # empty table
            for c in header:
                self.columns[c] = []

        # this assignment triggers creation of row template if index_name specified
        # but only if we have data
        if len(self.columns) > 0:
            self.index_name = index_name
        elif has_index:
            self._index_name = index_name

        # default title / legend to be empty strings
        self._title = str(title) if title else ""
        self._legend = str(legend) if legend else ""
        try:
            self._space = " " * space
        except TypeError:
            self._space = space
        self._digits = digits
        self._max_width = max_width

        # some attributes are not preserved in any file format, so always based
        # on args
        self._column_templates = column_templates or {}
        # define the repr() display policy
        random = 0
        self._repr_policy = dict(head=None, tail=None, random=random, show_shape=True)
        self.format = format
        self._missing_data = missing_data

    def __iter__(self):
        return iter(self.columns.iter_rows())

    def __len__(self):
        return self.columns._num_rows

    def __getitem__(self, names):
        # this is funky, but a side-effect of construction allowing setting
        # prior to having assigned the index_name column
        self.index_name

        if isinstance(names, tuple):
            rows, columns = names
        else:
            rows = names
            columns = self.columns.order

        if type(columns) in (str, int):
            columns = [columns]
        else:
            columns = self.columns._get_keys_(columns)

        columns = [self.columns._get_key_(c) for c in columns]

        # if a index_name has been specified we need to interpret
        # the provided values using the template
        if self._template:
            rows, _ = self._template.interpret_index(rows)
            rows = rows[0]

        if not hasattr(rows, "__len__") and not isinstance(rows, slice):
            rows = (rows,)

        if isinstance(rows, numpy.ndarray):
            rows = [i for i, v in enumerate(rows) if v]

        # if the length of rows and columns are both 1, return a single value
        if not isinstance(rows, slice) and len(rows) == len(columns) == 1:
            return self.columns[columns[0]][rows]

        attr = self._get_persistent_attrs()
        index_name = attr.pop("index_name")
        result = self.__class__(**attr)
        for c in columns:
            if len(self.columns[c]) == 0:
                continue

            result.columns[c] = self.columns[c][rows]

        if index_name in result.columns:
            result.index_name = index_name

        for c in self._column_templates:
            if c in result.columns:
                result._column_templates[c] = self._column_templates[c]

        return result

    def __getstate__(self):
        attrs = self._get_persistent_attrs()
        data = dict(init_table=attrs)
        cols = self.columns.to_rich_dict()
        data["data"] = cols
        return data

    def __setstate__(self, data):
        # we're not using these right now
        for k in ("type", "version"):
            data.pop(k, None)

        kwargs = data.pop("init_table")
        index = kwargs.pop("index_name")
        table = self.__class__(**kwargs)
        table.columns.__setstate__(data["data"])
        table.index_name = index
        self.__dict__.update(table.__dict__)

    def __repr__(self):
        if self.shape == (0, 0):
            return "0 rows x 0 columns"

        table, shape_info, unset_columns = self._get_repr_()
        if self.shape[0] == 0:
            return "\n".join([shape_info, unset_columns])

        if not self._repr_policy["show_shape"]:
            shape_info = ""
        return (
            "\n".join([str(table), shape_info, unset_columns])
            if unset_columns
            else "\n".join([str(table), shape_info])
        )

    def __str__(self):
        if self.shape == (0, 0):
            return ""

        return self.to_string(self.format)

    def _get_repr_(self):
        """returns a table for __repr__"""
        rn = self._repr_policy["random"]
        head = self._repr_policy["head"]
        tail = self._repr_policy["tail"]
        if not any([head, tail]):
            head, tail = (self.shape[0], None) if self.shape[0] < 50 else (5, 5)
            self._repr_policy["head"] = head
            self._repr_policy["tail"] = tail

        shape_info = ""
        if rn:
            indices = numpy.random.choice(self.shape[0], size=rn, replace=False)
            indices = list(sorted(indices))
            shape_info = f"Random selection of {rn} rows from"
        elif all([head, tail]):
            indices = list(range(head)) + list(
                range(self.shape[0] - tail, self.shape[0])
            )
            if head + tail < self.shape[0]:
                shape_info = f"Top {head} and bottom {tail} rows from"
        elif head:
            indices = list(range(head))
            if head < self.shape[0]:
                shape_info = f"Top {head} rows from"
        elif tail:
            indices = list(range(self.shape[0] - tail, self.shape[0]))
            if tail < self.shape[0]:
                shape_info = f"Bottom {tail} rows from"
        else:
            indices = list(range(self.shape[0]))

        shape_info += f"\n{self.shape[0]:,} rows x {self.shape[1]:,} columns"
        unset_columns = [c for c in self.header if not len(self.columns[c])]
        unset_columns = (
            f"unset columns: {', '.join(map(repr, unset_columns))}"
            if unset_columns
            else None
        )
        table = self[indices] if self.shape[0] > 0 else None

        return table, shape_info, unset_columns

    def _repr_html_(self):
        """returns html, used by Jupyter"""
        table, shape_info, unset_columns = self._get_repr_()
        if isinstance(table, numpy.ndarray):
            # single row / column
            table = self

        shape_info = (
            f"<p>{shape_info}; unset columns={unset_columns}</p>"
            if unset_columns
            else f"<p>{shape_info}</p>"
        )
        if not self._repr_policy["show_shape"]:
            shape_info = ""

        if self.shape[0] == 0:
            return shape_info

        html = table.to_html()
        # add elipsis if head + row < self.shape[0]
        html = html.splitlines()
        head = self._repr_policy.get("head") or self.shape[0]
        tail = self._repr_policy.get("tail") or self.shape[0]
        if head + tail < self.shape[0] and head and tail:
            HE = table_format.HtmlElement
            ellipsis = []
            for c in table.columns:
                if array_is_num_type(table.columns[c]):
                    css_class = "c3col_right"
                else:
                    css_class = "c3col_left"

                ellipsis.append(
                    str(HE(HE("...", "span", css_classes=[css_class]), "td"))
                )

            ellipsis = str(HE("".join(ellipsis), "tr", css_classes="ellipsis"))
            num_rows = 0
            for idx in range(len(html)):
                item = html[idx]
                if "<tr>" in item:
                    num_rows += 1
                if num_rows == head:
                    html.insert(idx + 1, ellipsis)
                    break

        html.insert(-1, shape_info)
        html = "\n".join(html)
        return html

    def _get_persistent_attrs(self):
        return UnionDict(self._persistent_attrs.copy())

    @property
    def title(self):
        return self._title

    @title.setter
    def title(self, value):
        self._title = value
        self._persistent_attrs["title"] = value

    @property
    def legend(self):
        return self._legend

    @legend.setter
    def legend(self, value):
        self._legend = value
        self._persistent_attrs["legend"] = value

    @property
    def space(self):
        return self._space

    @space.setter
    def space(self, value):
        try:
            self._space = " " * value
        except TypeError:
            self._space = value

        self._persistent_attrs["space"] = value

    def set_repr_policy(self, head=None, tail=None, random=0, show_shape=True):
        """specify policy for repr(self)

        Parameters
        ----------

        head : int
            number of top rows to included in represented display
        tail : int
            number of bottom rows to included in represented display
        random : int
            number of rows to sample randomly (supercedes head/tail)
        show_shape : bool
            boolean to determine if table shape info is displayed
        """
        if not any([head, tail, random]):
            self._repr_policy["show_shape"] = show_shape
            return
        if random:
            assert (
                type(random) == int and random > 0
            ), "random must be a positive integer"
            head = tail = None
        self._repr_policy = dict(
            head=head, tail=tail, random=random, show_shape=show_shape
        )

    @property
    def format(self):
        """the str display format"""
        return self._format

    @format.setter
    def format(self, new="simple"):
        """the str display format"""
        new = new.lower()
        if new not in table_format.known_formats:
            msg = (
                f"{new} not a supported format, see cogent3.format.table.known_formats"
            )
            raise ValueError(msg)

        self._format = new

    def format_column(self, column_head, format_template):
        """Provide a formatting template for a named column.

        Parameters
        ----------
        column_head
            the column label.
        format_template
            string formatting template or a function that will handle the formatting.
        """
        test_val = self.columns[column_head].tolist()[0]
        try:
            _ = (
                format_template(test_val)
                if callable(format_template)
                else format_template % test_val
            )
        except Exception as err:
            msg = f"{format_template} invalid for {column_head}: {err.args[0]}"
            raise ValueError(msg)

        self._column_templates[column_head] = format_template

    def head(self, nrows=5):
        """displays top nrows"""
        repr_policy = self._repr_policy
        nrows = min(nrows, self.shape[0])
        show_shape = self._repr_policy["show_shape"]
        self._repr_policy = dict(
            head=nrows, tail=None, random=None, show_shape=show_shape
        )
        display(self)
        self._repr_policy = repr_policy

    def tail(self, nrows=5):
        """displays bottom nrows"""
        repr_policy = self._repr_policy
        nrows = min(nrows, self.shape[0])
        show_shape = self._repr_policy["show_shape"]
        self._repr_policy = dict(
            head=None, tail=nrows, random=None, show_shape=show_shape
        )
        display(self)
        self._repr_policy = repr_policy

    @property
    def index_name(self):
        """column name whose values can be used to index table rows"""
        if self._index_name is not None and not self._template:
            self.columns.index_name = self._index_name
            self.index_name = self._index_name

        return self._index_name

    @index_name.setter
    def index_name(self, name):
        self.columns.index_name = name
        self._index_name = name
        self._persistent_attrs["index_name"] = name
        self._template = None if name is None else DictArrayTemplate(self.columns[name])

    @property
    def header(self):
        return self.columns.order

    @property
    def shape(self):
        return (self.columns._num_rows, len(self.columns))

    @property
    def array(self):
        return self.columns.array

    def cross_join(self, other, **kwargs):
        """cross join, or full outer join, of self with other

        Notes
        -----
        The column headers of the output are made unique by prepending
        other column headers with _, e.g. 'Name' becomes _Name'.
        """
        self_range = range(self.shape[0])
        other_range = range(other.shape[0])
        self_selected, other_selected = list(zip(*product(self_range, other_range)))
        joined_data = {c: self.columns[c].take(self_selected) for c in self.columns}
        col_prefix = "right" if not other.title else other.title
        other_data = {
            f"{col_prefix}_{c}": other.columns[c].take(other_selected)
            for c in other.columns
        }

        joined_data.update(other_data)
        new_header = list(self.columns.order) + [
            f"{col_prefix}_{c}" for c in other.columns
        ]
        attrs = self._get_persistent_attrs()
        attrs.pop("title", None)
        attrs |= kwargs
        joined = self.__class__(**attrs)
        for c in new_header:
            joined.columns[c] = joined_data[c]
        return joined

    def inner_join(
        self,
        other,
        columns_self=None,
        columns_other=None,
        use_index=True,
        **kwargs,
    ):
        """inner join of self with other

        Parameters
        ----------
        other
            A table object which will be joined with this
            table. other must have a title.
        columns_self, columns_other
            indices of key columns that will be compared in the join operation.
            Can be either column index, or a string matching the column header.
            The order matters, and the dimensions of columns_self and
            columns_other have to match. A row will be included in the output iff
            self[row, columns_self]==other[row, columns_other] for all i
        use_index
            if no columns specified and both self and other have a nominated
            index_name, this will be used.

        Notes
        -----
        The column headers of the output are made unique by prepending
        other column headers with _, e.g. 'Name' becomes _Name'.
        """
        col_prefix = "right" if not other.title else other.title

        if columns_self:
            columns_self = self.columns._get_keys_(columns_self)

        if columns_other:
            columns_other = other.columns._get_keys_(columns_other)

        columns_self = [columns_self] if isinstance(columns_self, str) else columns_self
        columns_other = (
            [columns_other] if isinstance(columns_other, str) else columns_other
        )
        columns_self = [columns_self] if isinstance(columns_self, str) else columns_self
        columns_other = (
            [columns_other] if isinstance(columns_other, str) else columns_other
        )
        if columns_self is columns_other is None and not use_index:
            # we do the natural inner join
            shared = set(self.columns) & set(other.columns)
            columns_self = [c for c in self.columns if c in shared]
            columns_other = [c for c in other.columns if c in shared]
        elif columns_self is columns_other is None:
            if not (self.index_name and other.index_name):
                msg = (
                    "indexes not specified, set use_index=False for natural inner join"
                )
                raise ValueError(msg)
            columns_self = [self.index_name]
            columns_other = [other.index_name]
        elif columns_self is None or columns_other is None:
            # the same column labels will be used for both tables
            columns_self = columns_self or columns_other
            columns_other = columns_self or columns_other

        if len(columns_self) != len(columns_other):
            raise RuntimeError(
                "Error during table join: key columns have different dimensions!"
            )

        output_mask = [c for c in other.columns if c not in columns_other]

        other_row_index = defaultdict(list)
        subtable = other[:, columns_other]
        for row_index, row in enumerate(subtable.columns.array):
            # insert new entry for each row
            other_row_index[tuple(row)].append(row_index)

        other_selected = []
        self_selected = []
        subtable = self[:, columns_self]
        for row_index, row in enumerate(subtable.columns.array):
            # assemble key for query of other
            key = tuple(row)
            if key not in other_row_index:
                continue

            self_selected.extend([row_index] * len(other_row_index[key]))
            other_selected.extend(other_row_index[key])

        joined_data = {c: self.columns[c][self_selected] for c in self.columns}
        other_data = {
            f"{col_prefix}_{c}": other.columns[c][other_selected] for c in output_mask
        }

        joined_data.update(other_data)
        new_header = list(self.columns.order) + [
            f"{col_prefix}_{c}" for c in output_mask
        ]
        attr = self._get_persistent_attrs()
        attr.pop("title", None)
        attr |= kwargs
        joined = self.__class__(**attr)
        for c in new_header:
            joined.columns[c] = joined_data[c]
        return joined

    def joined(
        self,
        other,
        columns_self=None,
        columns_other=None,
        inner_join=True,
        **kwargs,
    ):
        """returns a new table containing the join of this table and
        other. See docstring for inner_join, or cross_join
        """
        if not inner_join:
            assert (
                columns_self is columns_other is None
            ), "Cannot specify column indices for a cross join"
            return self.cross_join(other, **kwargs)

        return self.inner_join(
            other=other,
            columns_self=columns_self,
            columns_other=columns_other,
            use_index=False,
            **kwargs,
        )

    # todo check the type info
    # todo implement negate argument
    # todo implement check that callable returns bool
    def get_row_indices(self, callback, columns, negate=False):
        """returns boolean array of callback values given columns"""
        subset = self[:, columns]
        data = subset if not isinstance(callback, Callable) else subset.array
        num_columns = len(columns)
        match = not negate
        return numpy.array(
            [
                True
                if _callback(callback, row=row, num_columns=num_columns) == match
                else False
                for row in data
            ]
        )

    def filtered(self, callback, columns=None, **kwargs):
        """Returns a table with rows satisfying the provided callback function.

        Parameters
        ----------
        columns
            the columns whose values determine whether a row is to be included.
        callback
            Can be a function, which takes rows and returns True/False, or a
            string representing valid python code to be evaluated.

        Notes
        -----
        Row data provided to callback is a 1D list if more than one column,
        single value (row[col]) otherwise.
        """
        # no point filtering if no rows, justv return self
        if self.shape[0] == 0:
            return self

        if isinstance(columns, str):
            columns = (columns,)

        if columns is None:
            columns = self.columns.order

        indices = self.get_row_indices(callback=callback, columns=columns)
        attr = self._get_persistent_attrs()
        attr |= kwargs
        result = self.__class__(**attr)
        for c in self.columns:
            result.columns[c] = self.columns[c][indices]
        return result

    def filtered_by_column(self, callback, **kwargs):
        """Returns a table with columns identified by callback

        Parameters
        ----------
        callback
            A function which takes the columns delimited by columns and returns
            True/False, or a string representing valid python code to be evaluated.
        """
        columns = [c for c in self.columns if callback(self.columns[c])]
        attr = self._get_persistent_attrs()
        attr |= kwargs
        result = self.__class__(**attr)
        for c in columns:
            result.columns[c] = self.columns[c]
        return result

    def count(self, callback, columns=None, **kwargs):
        """Returns number of rows for which the provided callback
        function returns True when passed row data from columns. Row data
        is a 1D list if more than one column, raw row[col] value otherwise.

        Parameters
        ----------
        columns
            the columns whose values determine whether a row is to
            be included.
        callback
            Can be a function, which takes the sub
            by columns and returns True/False, or a string representing valid
            python code to be evaluated.

        """
        # no rows, value must be 0
        if self.shape[0] == 0:
            return 0

        if isinstance(columns, str):
            columns = (columns,)

        if columns is None:
            columns = self.columns.order

        indices = self.get_row_indices(callback=callback, columns=columns)
        return indices.sum()

    def count_unique(self, columns=None):
        """count occurrences of unique combinations of columns

        Parameters
        ----------
        columns
            name of one or more columns. If None, all columns are used

        Returns
        -------
        CategoryCounter instance
        """
        from cogent3.maths.stats.number import CategoryCounter

        if columns is None:
            columns = self.columns.order

        subset = self.columns.take_columns(columns)
        if len(subset) == 1:
            data = subset[0].tolist()
        else:
            data = subset.array
            data = [tuple(e) for e in data]

        return CategoryCounter(data=data)

    def distinct_values(self, columns):
        """returns the set of distinct values for the named column(s)"""
        data = [tuple(r) for r in self[:, columns].array.tolist()]
        result = set(data)
        result = {d[0] if len(d) == 1 else d for d in result}
        return result

    def appended(self, new_column, *tables, **kwargs):
        """Concatenates an arbitrary number of tables together

        Parameters
        ----------
        new_column
            provide a heading for the new column, each tables
            title will be placed in it. If value is false, the result is no
            additional column.
        tables
            series of Table instances

        Notes
        -----
        All tables must have the same columns. If a column dtype differs between tables,
        dtype for that column in result is determined by numpy.
        """
        if new_column is not None:
            assert new_column not in self.columns, f"'{new_column}' already exists"
        # default title is no title
        kwargs["title"] = kwargs.get("title", "")
        attr = self._get_persistent_attrs()
        attr |= kwargs
        result = self.__class__(**attr)
        # convert series of tables
        if isinstance(tables[0], tuple) or isinstance(tables[0], list):
            tables = tuple(tables[0])
        # for each table, determine it's number of rows and create an
        # equivalent length vector of its title
        columns = set(self.columns.order)
        new_col = []
        table_series = (self,) + tables
        raw_data = defaultdict(list)
        dtypes = defaultdict(set)
        for table in table_series:
            assert set(table.columns.order) == columns, "columns don't match"
            if new_column is not None:
                new_col.extend([table.title] * table.shape[0])

            for c in table.columns:
                dtypes[c].add(table.columns[c].dtype)

            data = table.columns.to_dict()
            for c, v in data.items():
                raw_data[c].extend(v)

        if new_column is not None:
            columns = (new_column,) + self.columns.order
            raw_data[new_column] = new_col
            dtypes[new_column] = {"U"}
        else:
            columns = self.columns.order

        for c in columns:
            data = (
                raw_data[c]
                if len(dtypes[c]) != 1
                else numpy.array(raw_data[c], dtype=dtypes[c].pop())
            )
            result.columns[c] = data
        return result

    def get_columns(self, columns, with_index=True):
        """select columns from self with index_name unless excluded

        Parameters
        ----------
        columns : string or sequence of strings
            names of columns
        with_index : bool
            If index_name is set, includes with columns.

        Returns
        -------
        Table
        """
        if self.index_name and with_index:
            columns = [self.index_name] + [c for c in columns if c != self.index_name]
        return self[:, columns]

    def with_new_column(self, new_column, callback, columns=None, dtype=None, **kwargs):
        """Returns new table with an additional column, computed using callback.

        Parameters
        ----------
        new_column
            new column heading
        columns
            the columns whose values determine whether a row is to be included.
        callback
            Can be a function, which takes the subtable by columns and returns
            True/False, or a string representing valid python code to be evaluated.
        dtype
            numpy type of result
        """
        attr = self._get_persistent_attrs()
        index = attr.pop("index_name")
        attr |= kwargs
        result = self.__class__(**attr)
        for c in self.columns:
            if c == new_column:
                continue
            result.columns[c] = self.columns[c]

        if columns is None:
            columns = self.columns.order

        if isinstance(columns, str):
            columns = (columns,)

        subset = self[:, columns]
        data = subset if not isinstance(callback, Callable) else subset.array
        num_columns = len(columns)
        values = numpy.array(
            [_callback(callback, row=row, num_columns=num_columns) for row in data]
        )

        if dtype:
            values = numpy.array(values, dtype=dtype)

        result.columns[new_column] = values

        if index in result.columns:
            result.index_name = index

        return result

    # todo deprecate this method
    def with_new_header(self, old, new, **kwargs):
        """returns a new Table with old header labels replaced by new

        Parameters
        ----------
        old
            the old column header(s). Can be a string or series of them.
        new
            the new column header(s). Can be a string or series of them.
        """
        if isinstance(old, str):
            old = [old]
            new = [new]

        assert len(old) == len(new), "Mismatched number of old/new labels"
        attr = self._get_persistent_attrs()
        attr |= kwargs
        result = self.__class__(**attr)
        for c in self.columns:
            key = c
            if c in old:
                index = old.index(c)
                key = new[index]
            result.columns[key] = self.columns[c]
        return result

    def sum_columns(self, columns=None, strict=True):
        """return sum of indicated columns

        Parameters
        ----------
        columns
            column name(s) or indices
        strict
            if False, ignores cells with non column/row.

        """
        if columns is None:
            columns = self.columns.order

        if isinstance(columns, str) or isinstance(columns, int):
            columns = [columns]

        columns = self.columns[columns]
        func = numpy.sum if strict else _numeric_sum
        result = [func(c) for c in columns]
        if len(result) == 1:
            result = result[0]
        return result

    def sum_rows(self, indices=None, strict=True):
        """return sum of indicated rows

        Parameters
        ----------
        indices
            row indices
        strict
            if False, ignores cells with non numeric values.
        """
        data = self.array if indices is None else self[indices, :].array
        # a multi-rowed result
        if strict:
            result = data.sum(axis=1).tolist()
        else:
            result = [_numeric_sum(row) for row in data]

        if len(result) == 1:
            result = result[0]

        return result

    # todo change indices to columns
    def summed(self, indices=None, col_sum=True, strict=True):
        """returns the sum of numerical values for column(s)/row(s)

        Parameters
        ----------
        indices
            column name(s) or indices or row indices
        col_sum
            sums values in the indicated column, the default. If
            False, returns the row sum.
        strict
            if False, ignores cells with non column/row.
        """
        if col_sum:
            return self.sum_columns(columns=indices, strict=strict)

        return self.sum_rows(indices=indices, strict=strict)

    def normalized(self, by_row=True, denominator_func=None, **kwargs):
        """returns a table with elements expressed as a fraction according
        to the results from func

        Parameters
        ----------
        by_row
            normalisation done by row
        denominator_func
            a callback function that takes an array and
            returns a value to be used as the denominator. Default is sum.

        """
        attr = self._get_persistent_attrs()
        attr |= kwargs
        result = self.__class__(**attr)
        denominator_func = denominator_func if callable(denominator_func) else numpy.sum
        if not by_row:
            for c in self.columns:
                v = self.columns[c]
                result.columns[c] = v / denominator_func(v)

            return result

        totals = numpy.array([denominator_func(r) for r in self.array])
        for i, c in enumerate(self.columns):
            result.columns[c] = self.columns[i] / totals

        return result

    def sorted(self, columns=None, reverse=None, **kwargs):
        """Returns a new table sorted according to columns order.

        Parameters
        ----------
        columns
            column headings, their order determines the sort order.
        reverse
            column headings, these columns will be reverse sorted.

        Notes
        -----
        Either can be provided as just a single string, or a series of
        strings. If only reverse is provided, that order is used.
        """
        if "reversed" in kwargs:
            raise TypeError(
                "got an unexpected keyword argument 'reversed', " "use 'reverse'"
            )

        reverse = reverse if reverse is not None else []
        if reverse != [] and columns is None:
            columns = reverse

        if columns is None:
            columns = list(self.columns)

        if isinstance(columns, str):
            columns = [columns]

        if isinstance(reverse, str):
            reverse = [reverse]

        columns = list(columns)

        if reverse and not (set(columns) & set(reverse)):
            for c in reverse:
                if c in columns:
                    continue

                columns.append(c)

        dtypes = [(c, self.columns[c].dtype) for c in columns]
        data = numpy.array(self.columns[columns], dtype="O").T
        for c in reverse:
            index = columns.index(c)
            func = _reverse_num if array_is_num_type(self.columns[c]) else _reverse_str
            func = numpy.vectorize(func)
            data[:, index] = func(data[:, index])

        data = numpy.rec.fromarrays(data.copy().T, dtype=dtypes)
        indices = data.argsort()

        attr = self._get_persistent_attrs()
        attr |= kwargs
        result = Table(**attr)
        for c in self.columns:
            result.columns[c] = self.columns[c][indices]

        return result

    def _formatted_by_col(self, missing_data="", pad=True):
        """returns self as formatted strings

        Parameters
        ----------
        missing_data : str
            default str value for missing
        pad : bool
            if True, adds padding

        Returns
        -------
        dict of formatted columns, list of [(name, length),..]
        """
        missing_data = missing_data or self._missing_data
        formatted = {}
        col_widths = []
        for c in self.columns.order:
            data = self.columns[c]
            if len(data) == 0:
                continue

            format_spec = self._column_templates.get(c, None)
            frmt, c, width = table_format.formatted_array(
                data,
                c,
                format_spec=format_spec,
                missing_data=missing_data,
                precision=self._digits,
                pad=pad,
            )
            col_widths.append((c, width))
            formatted[c] = frmt

        return formatted, col_widths

    def _formatted(self, missing_data="", stripped=False):
        """returns self as formatted strings

        Parameters
        ----------
        missing_data : str
            default str value for missing
        stripped : bool
            if True, removes padding

        """
        formatted_cols, _ = self._formatted_by_col(
            missing_data=missing_data, pad=not stripped
        )
        ordered = [(self.columns.order.index(c.strip()), c) for c in formatted_cols]
        ordered.sort()
        formatted = [[c] + formatted_cols[c] for _, c in ordered]
        formatted = [list(e) for e in zip(*formatted)]
        if not formatted and self.header:
            formatted = [self.header]
        return formatted

    def to_csv(self, with_title=False, with_legend=False):
        """return table formatted as comma separated values

        Parameters
        ----------
        with_title : bool
            include the table title
        with_legend : bool
            include table legend

        Returns
        -------
        str
        """
        formatted_table = self._formatted(stripped=True)
        header = formatted_table.pop(0)
        title = self.title if with_title else None
        legend = self.legend if with_legend else None
        return table_format.separator_format(
            header, formatted_table, title=title, legend=legend, sep=","
        )

    def to_latex(
        self, concat_title_legend=True, justify=None, label=None, position=None
    ):
        """Returns the text a LaTeX table.

        Parameters
        ----------
        concat_title_legend : bool
            the table caption is formed by concatenating the table title and legend
        justify
            column justification, default is right aligned.
        label
            for cross referencing
        position
            table page position, default is here, top separate page

        Notes
        -----
        The \\caption*{} command is provided with the caption package. See
        https://ctan.org/pkg/caption for more details.
        """
        formatted_table = self._formatted()
        header = formatted_table.pop(0)
        caption = self.title or None
        legend = self.legend or None
        if concat_title_legend and (caption or legend):
            caption = " ".join([caption or "", legend or ""])
            caption = caption.strip()
            legend = None
        return table_format.latex(
            formatted_table,
            header,
            caption=caption,
            legend=legend,
            justify=justify,
            label=label,
            position=position,
        )

    def to_markdown(self, space=1, justify=None):
        """
        returns markdown formatted table

        Parameters
        ----------
        space
            number of spaces surrounding the cell contents, must be >= 1
        justify
            characters indicating alignment of columns

        Returns
        -------
        str
        """
        formatted_table = self._formatted()
        header = formatted_table.pop(0)
        return table_format.markdown(
            header, formatted_table, space=space, justify=justify
        )

    def to_rst(self, csv_table=False):
        """returns rst formatted table

        Parameters
        ----------
        csv_table : bool
            use csv-directive, grid table otherwise

        Returns
        -------
        str
        """
        stripped = csv_table
        formatted_table = self._formatted(stripped=stripped)
        header = formatted_table.pop(0)
        if csv_table:
            result = table_format.rst_csv_table(
                header, formatted_table, title=self.title, legend=self.legend
            )
        else:
            result = table_format.grid_table_format(
                header, formatted_table, title=self.title, legend=self.legend
            )
        return result

    def to_string(
        self,
        format="",
        borders=True,
        sep=None,
        center=False,
        concat_title_legend=True,
        **kwargs,
    ):
        """Return the table as a formatted string.

        Parameters
        ----------
        format
            possible formats are 'rest'/'rst', 'markdown'/'md',
            'latex', 'html', 'phylip', 'bedgraph', 'csv', 'tsv', or 'simple'
            (default).
        sep
            A string separator for delineating columns, e.g. ',' or
            '\t'. Overrides format.
        center : bool
            content is centered in the column, default is right
            justified
        concat_title_legend : bool
            Concat the title and legend.

        Notes
        -----
        If format is bedgraph, assumes that column headers are chrom, start,
        end, value. In that order!
        """
        if format == "bedgraph":
            # todo remove requirement for column order
            assert self.shape[1] == 4, "bedgraph format is for 4 column tables"
            # assuming that header order is chrom, start, end, val
            formatted_table = bedgraph.bedgraph(self.sorted().array.tolist(), **kwargs)
            return formatted_table

        if format.lower() in ("tsv", "csv"):
            sep = sep or {"tsv": "\t", "csv": ","}[format.lower()]
            format = ""

        if sep != "\t":
            sep = sep.strip() if sep else None

        if sep == ",":
            return self.to_csv(**kwargs)

        if sep == "\t":
            return self.to_tsv(**kwargs)

        if format in ("rest", "rst"):
            return self.to_rst(**kwargs)

        if format in ("markdown", "md"):
            return self.to_markdown(**kwargs)

        if format.endswith("tex"):
            return self.to_latex(concat_title_legend=concat_title_legend, **kwargs)

        if format == "html":
            return self.to_html(**kwargs)

        if format == "phylip":
            # need to eliminate row identifiers
            columns = [c for c in self.columns if c != self.index_name]
            table = self[:, columns]
            formatted_table = table._formatted(missing_data="0.0000")
            header = formatted_table.pop(0)
            return table_format.phylip_matrix(formatted_table, header)

        # convert self to a 2D list after caching current column templates
        col_formats = {}
        for c in self.columns:
            if c in self._column_templates:
                col_formats[c] = self._column_templates[c]
                continue

            col_formats[c] = ">" if array_is_num_type(self.columns[c]) else "<"

        orig_formats = self._column_templates
        self._column_templates = col_formats

        formatted_table = self._formatted(stripped=sep is not None)
        self._column_templates = orig_formats

        header = formatted_table.pop(0)
        args = (header, formatted_table, self.title, self.legend)

        if sep:
            return table_format.separator_format(*args, sep=sep)

        return table_format.simple_format(
            *args + (self._max_width, self.index_name, borders, self.space)
        )

    def to_tsv(self, with_title=False, with_legend=False):
        """return table formatted as tab separated values

        Parameters
        ----------
        with_title : bool
            include the table title
        with_legend : bool
            include table legend

        Returns
        -------
        str
        """
        formatted_table = self._formatted(stripped=True)
        header = formatted_table.pop(0)
        title = self.title if with_title else None
        legend = self.legend if with_legend else None
        return table_format.separator_format(
            header, formatted_table, title=title, legend=legend, sep="\t"
        )

    def to_html(self, column_alignment=None):
        """construct html table

        Parameters
        ----------
        column_alignment : dict
            {col_name: alignment character, ...} where alignment character
            can be one of 'l', 'c', 'r'. Defaults to 'r' for numeric columns,
            'l' for text columns.

        Returns
        -------
        string

        Notes
        -----
        Placed within c3table div element, embeds CSS style.
        """
        column_alignment = column_alignment or {}
        for c, v in column_alignment.items():
            v = v.lower()
            if v not in "clr":
                raise ValueError(f"invalid column alignment value {v} for '{c}'")

            column_alignment[c] = {
                "c": "c3col_center",
                "l": "c3col_left",
                "r": "c3col_right",
            }[v]

        base_colour = "rgba(161, 195, 209, {alpha})"
        # alpha applied to the index_name column background
        alpha = 0.0 if self.index_name is None else 0.25

        style = table_format.css_c3table_template % dict(
            colour=base_colour.format(alpha=alpha),
            head_colour=base_colour.format(alpha=0.75),
        )

        HtmlElement = table_format.HtmlElement
        tables = [str(HtmlElement(style, "style", newline=True))]

        cols, widths = self._formatted_by_col(pad=False)
        headers = table_format.get_continuation_tables_headers(
            widths,
            index_name=self.index_name,
            max_width=self._max_width,
        )

        for c in self.columns:
            c = c.strip()
            css_class = (
                "c3col_right" if array_is_num_type(self.columns[c]) else "c3col_left"
            )

            css_class = [column_alignment.get(c, css_class)]
            cols[c] = [
                HtmlElement(
                    HtmlElement(v, "span", css_classes=css_class),
                    "td",
                    css_classes=["index"] if c == self.index_name else None,
                )
                for v in cols[c]
            ]

        title = self.title or ""
        if title and not table_format.is_html_markup(title):
            title = escape(title)

        legend = self.legend or ""
        if legend and not table_format.is_html_markup(legend):
            legend = escape(legend)

        for i, header in enumerate(headers):
            if title and i == 0:
                st = HtmlElement(title, "span", css_classes=["cell_title"])
            elif title:
                st = HtmlElement("continuation", "span", css_classes=["cell_title"])
            else:
                st = None

            if legend and i == 0:
                legend = HtmlElement(legend, "span", css_classes=["cell_legend"])
                legend = HtmlElement(legend, "br") if title else legend
                st = f"{st} {legend}" if st else f"{legend}"

            caption = str(HtmlElement(st, "caption", newline=True)) if st else ""
            rows = []
            for i, row in enumerate(zip(*[cols[c] for c in header])):
                txt = HtmlElement("".join(str(e) for e in row), "tr")
                rows.append(str(txt))

            rows = str(HtmlElement("\n".join(rows), "tbody", newline=True))

            if column_alignment:
                header = [
                    HtmlElement(c, "span", css_classes=column_alignment.get(c, None))
                    for c in header
                ]

            header = "".join(str(HtmlElement(c, "th")) for c in header)
            header = str(
                HtmlElement(header, "thead", css_classes=["head_cell"], newline=True)
            )

            subtable = HtmlElement(
                "".join([caption, header, rows]), "table", newline=True
            )
            tables.append(str(subtable))

        return str(
            HtmlElement("\n".join(tables), "div", css_classes=["c3table"], newline=True)
        )

    def tolist(self, columns=None):
        """Returns raw data as a list

        Parameters
        ----------
        columns
            if None, all data are returned

        Notes
        -----
        If one column, a 1D list is returned.
        """
        if columns is None:
            columns = self.columns.order

        columns = [columns] if isinstance(columns, str) else columns
        if len(columns) == 1:
            result = self.columns[columns[0]].tolist()
            return result

        subtable = self.get_columns(columns)
        result = subtable.columns.array.tolist()

        return result

    @extend_docstring_from(DictArray.to_dict)
    def to_dict(self, flatten=False):
        index = self.columns[self.index_name] if self.index_name else self.shape[0]
        template = DictArrayTemplate(index, self.columns.order)
        darr = template.wrap(self.array)
        return darr.to_dict(flatten=flatten)

    def to_rich_dict(self):
        data = self.__getstate__()
        data["type"] = get_object_provenance(self)
        data["version"] = None  # todo
        return data

    def to_json(self):
        data = self.to_rich_dict()
        return json.dumps(data)

    def to_dataframe(self, categories=None):
        """returns pandas DataFrame instance

        Parameters
        ----------
        categories
            converts these columns to category dtype in the data
            frame. Note, categories are not ordered.
        """
        try:
            from pandas import DataFrame
        except ImportError:
            raise ImportError("pandas not installed")

        index = None if not self.index_name else self.columns[self.index_name]
        data = {c: self.columns[c] for c in self.columns if c != self.index_name}
        df = DataFrame(data=data, index=index)
        if categories is not None:
            categories = [categories] if type(categories) == str else categories
            df = df.astype({n: "category" for n in categories})

        return df

    def to_plotly(self, width=500, font_size=12, layout=None, **kwargs):
        """returns a Plotly Table"""
        from cogent3.draw.drawable import Drawable

        columns, _ = self._formatted_by_col(pad=False)
        header = [f"<b>{c}</b>" for c in self.header]
        if self.index_name:
            body_colour = ["white"] * self.shape[0]
            index_colour = ["rgba(161, 195, 209, 0.5)"] * self.shape[0]
            colours = [index_colour] + [body_colour[:] for i in range(self.shape[1])]
            columns[self.index_name] = [f"<b>{e}</b>" for e in columns[self.index_name]]
        else:
            colours = "white"

        tab = UnionDict(
            type="table",
            header=dict(
                values=header,
                fill=dict(color="rgba(161, 195, 209, 1)"),
                font=dict(size=font_size),
                align="center",
            ),
            cells=dict(
                values=[columns[c] for c in self.header], fill=dict(color=colours)
            ),
        )

        draw = Drawable()
        aspect_ratio = self.shape[0] / self.shape[1]
        layout = layout or {}
        default_layout = dict(
            width=width,
            height=aspect_ratio * width,
            autosize=False,
            title=self.title,
            margin=dict(l=10, r=10, t=30, b=10, pad=10),
        )
        default_layout.update(layout)
        draw.traces.append(tab)
        draw.layout |= default_layout
        return draw

    def to_categorical(self, columns=None, index_name=None):
        """construct object that can be used for statistical tests

        Parameters
        ----------
        columns
            columns to include. These correspond to contingency column
            labels. The row labels come from values under the index_name
            column. Defaults to all columns.

        Returns
        -------
        CategoryCounts, an object for performing statistical tests on
        contingency tables.

        Notes
        -----
        Only applies to cases where an index_name is defined. The selected columns
        must be int types and represent the counts of corresponding categories.
        """
        from cogent3.maths.stats.contingency import CategoryCounts
        from cogent3.util.dict_array import DictArrayTemplate

        self.index_name = index_name if index_name is not None else self.index_name
        if self.index_name is None:
            raise ValueError("requires index_name be set")

        columns = list(self.header) if columns is None else columns

        columns = [columns] if isinstance(columns, str) else columns
        if not set(columns) <= set(self.header):
            raise ValueError(f"unknown columns {columns}")

        if self.index_name in columns:
            columns.remove(self.index_name)
        row_cats = self.columns[self.index_name]
        # must be convertible to int
        for col in columns:
            if "int" not in self.columns[col].dtype.name:
                raise TypeError(f"{col} is not of int type")

        matrix = self.get_columns(columns, with_index=False).array.astype(int)

        data = DictArrayTemplate(row_cats, columns).wrap(matrix)
        return CategoryCounts(data)

    def transposed(self, new_column_name, select_as_header=None, **kwargs):
        """returns the transposed table.

        Parameters
        ----------
        new_column_name
            the existing header will become a column with
            this name
        select_as_header
            current column name containing data to be used
            as the header. Defaults to the first column.
        """
        select_as_header = select_as_header or self.columns.order[0]
        assert (
            select_as_header in self.columns
        ), f'"{select_as_header}" not in table header'

        if len(self.distinct_values(select_as_header)) != len(self):
            raise ValueError(f"not all '{select_as_header}' values unique")

        attr = self._get_persistent_attrs()
        # on transpose, a row index_name becomes a column, so pop
        del attr["index_name"]

        attr |= kwargs
        result = self.__class__(**attr)

        # place new column header first
        columns = [select_as_header] + [
            c for c in self.columns if c != select_as_header
        ]
        data = self[:, columns].array
        result.columns[new_column_name] = columns[1:]
        for row in data.tolist():
            c = str(row.pop(0))
            result.columns[c] = row
        return result

    def write(
        self,
        filename,
        mode=None,
        writer=None,
        format=None,
        sep=None,
        compress=None,
        **kwargs,
    ):
        """Write table to filename in the specified format.

        Parameters
        ----------
        mode
            file opening mode
        format
            Valid formats are those of the to_string method plus pickle. Will
            try and guess from filename if not specified.
        writer
            a function for formatting the data for output.
        sep
            a character delimiter for fields.
        compress
            if True, gzips the file and appends .gz to the filename (if not
            already added).

        Notes
        -----
        If a format is not specified, it attempts to use a filename suffix.
        Unformatted numerical values are written to file in order to preserve
        numerical accuracy.
        """
        filename = pathlib.Path(filename)
        file_suffix, compress_suffix = get_format_suffixes(filename)
        format = format or file_suffix
        compress = compress or compress_suffix is not None

        mode = mode or {"pickle": "wb"}.get(format, "w")

        if format == "json":
            with atomic_write(filename, mode="wt") as f:
                f.write(self.to_json())
            return

        if compress:
            if ".gz" not in filename.suffixes:
                filename = pathlib.Path(f"{filename}.gz")
            mode = "wt"

        outfile = atomic_write(filename, mode=mode)

        format = format if format else file_suffix

        if format == "csv":
            sep = sep or ","
        elif format == "tsv":
            sep = sep or "\t"

        if writer:
            rows = self.tolist()
            rows.insert(0, self.header[:])
            rows = writer(rows, has_header=True)
            outfile.write("\n".join(rows))
        elif format == "pickle":
            data = self.__getstate__()
            pickle.dump(data, outfile, protocol=1)
        elif sep is not None and format != "bedgraph":
            writer = csv.writer(outfile, delimiter=sep, lineterminator="\n")
            if self.title:
                writer.writerow([self.title])
            writer.writerow(self.header)
            writer.writerows(self.array)
            if self.legend:
                writer.writerow([self.legend])
        else:
            table = self.to_string(format=format, sep=sep, **kwargs)
            outfile.write(table + "\n")

        outfile.close()
