#!/usr/bin/env python
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
import pickle
import re
import warnings

from collections import defaultdict
from collections.abc import Callable, MutableMapping
from itertools import product
from xml.sax.saxutils import escape

import numpy

from cogent3.format import bedgraph
from cogent3.format import table as table_format
from cogent3.util.dict_array import DictArray, DictArrayTemplate
from cogent3.util.misc import (
    extend_docstring_from,
    get_format_suffixes,
    get_object_provenance,
    open_,
)
from cogent3.util.union_dict import UnionDict
from cogent3.util.warning import deprecated


try:
    from IPython.display import display
except ImportError:
    display = lambda x: print(repr(x))

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = ["Gavin Huttley", "Felix Schill"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
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


def formatted_array(
    series, title="", precision=4, format_spec=None, missing_data="", center=False,
):
    """converts elements in a numpy array series to an equal length string.

    Parameters
    ----------
    series
        the series of table rows
    title
        title of series
    precision
        number of decimal places. Can be overridden by following.
    format_spec
        format specification as per the python Format Specification, Mini-Language
        or a callable function.
    missing_data
        default missing data value.

    Returns
    -------
    list of formatted series, formatted title
    """
    if callable(format_spec):
        formatter = format_spec
        format_spec = base_format = ""
    else:
        formatter = None

    if isinstance(format_spec, str):
        format_spec = format_spec.replace("%", "")

    if format_spec:
        match = re.search("[<>^]", format_spec[:2])
        final_align = ">" if match is None else match.group()
        align = ""
    else:
        final_align = align = ">"

    base_format = format_spec if format_spec else ""
    assert isinstance(series, numpy.ndarray), "must be numpy array"
    if format_spec is None:
        type_name = series.dtype.name
        align = "^" if center else ">"
        if "int" in type_name:
            base_format = "d"
        elif "float" in type_name:
            base_format = f".{precision}f"
        elif "bool" == type_name:
            base_format = ""
        else:
            # handle mixed types with a custom formatter
            formatter = _MixedFormatter(
                align, len(title), precision, missing_data=missing_data
            )
            format_spec = base_format = ""

        format_spec = base_format

    formatted = []
    max_length = len(title)
    for i, v in enumerate(series):
        if formatter:
            v = formatter(v)
        else:
            try:
                v = format(v, format_spec)
            except (TypeError, ValueError):
                # could be a python object
                v = str(v)
        l = len(v)
        if l > max_length:
            max_length = l
            format_spec = f"{align}{max_length}{base_format}"
        formatted.append(v)

    # title is always right aligned, for now
    title = format(title, f">{max_length}")
    # now adjust to max_len
    format_spec = f"{final_align}{max_length}s"
    for i in range(len(series)):
        if len(formatted[i]) < max_length:
            formatted[i] = format(formatted[i].strip(), format_spec)
    return formatted, title


def cast_str_to_numeric(values):
    """converts a series of strings to numeric values"""
    if not isinstance(values, numpy.ndarray):
        values = numpy.array(values, dtype="U")

    for typ in (int, float, complex):
        try:
            values = values.astype(typ)
            break
        except Exception:
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

    result = {c: [data[c][r] for r in row_order] for c in data}
    return result


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


class _MixedFormatter:
    """handles formatting of mixed data types"""

    def __init__(
        self, alignment, length, precision=4, float_type="f", missing_data=None
    ):
        self.missing_data = missing_data
        self.length = length
        self.alignment = alignment
        self.precision = precision
        self.float_type = float_type

    def __call__(self, val):
        prefix = f"{self.alignment}{self.length}"
        float_spec = f"{prefix}.{self.precision}{self.float_type}"
        int_spec = f"{prefix}d"
        result = str(val)
        if self.missing_data is not None and not result:
            return self.missing_data

        for fspec in (int_spec, float_spec, prefix):
            try:
                result = format(val, fspec)
                break
            except (TypeError, ValueError):
                pass

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
            key = [self._get_key_(k) for k in key]
        elif isinstance(key, numpy.ndarray):
            # we try slicing by array
            cols = numpy.array(self.order, dtype="U")
            try:
                key = cols[key]
            except Exception:
                msg = f"{key} could not be used to slice columns"
                raise KeyError(msg)
        else:
            raise KeyError(f"{key}")

        return key

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
        key = str(key)
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
        txt = f"{self.__class__.__name__}({', '.join(v)})"
        return txt

    def __str__(self):
        return repr(self)

    def iter_rows(self):
        columns = [self[c] for c in self]
        for row in zip(*columns):
            yield self._template.wrap(row, dtype=object)

    @property
    def index_name(self):
        # check
        return self._index_name

    @index_name.setter
    def index_name(self, name):
        if name is None:
            return

        if name not in self:
            raise ValueError(f"'{name}' unknown, index must be an existing column")

        # make sure index has unique values
        unique = set(self[name])
        if len(unique) != self._num_rows:
            raise ValueError(f"cannot use '{name}' as index, not all values unique")

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
        arr = [self[k] for k in self.order]
        return numpy.array(arr, dtype="O").T

    @property
    def order(self):
        """column order"""
        # if index_name not first, we re-order
        if self._index_name and self._order[0] != self._index_name:
            order = [self._index_name] + [
                c for c in self._order if c != self._index_name
            ]
            self._order = tuple(order)
        return self._order

    def to_dict(self):
        """returns column based dict"""
        result = {c: self[c].tolist() for c in self}
        return result

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
        row_ids=None,
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
        attrs = {
            k: v
            for k, v in locals().items()
            if k not in ("self", "__class__", "data", "header", "kwargs")
        }
        rows = kwargs.pop("rows", None)

        assert not (rows and data), "rows is deprecated, use data"
        if rows:
            # todo deprecation warning
            data = rows

        attrs.update(kwargs)

        self._persistent_attrs = attrs

        self.columns = Columns()
        self._template = None

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

        if header is None and isinstance(data, dict):
            header = list(data)

        if row_ids:
            if row_ids == True:
                row_ids = header[0]
                deprecated("argument", "row_ids: bool", "row_ids: string", "2020.6")
                warnings.warn(
                    "support for row_ids as bool discontinued in "
                    "version 2020.6, use a column name instead",
                    DeprecationWarning,
                )

        if data:
            row_order = kwargs.get("row_order", None)
            data = cast_to_1d_dict(data, row_order=row_order)
            if row_ids:
                try:
                    self.columns[row_ids] = data[row_ids]
                except KeyError:
                    raise ValueError(f"'{row_ids}' not in data")

            for c in header:
                if c == row_ids:
                    continue
                self.columns[c] = data[c]

        elif header:
            # empty table
            for c in header:
                self.columns[c] = []

        if row_ids:
            self._index_name = row_ids
        else:
            self._index_name = None

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

        self.format = format

        # define the repr() display policy
        random = 0
        self._repr_policy = dict(head=None, tail=None, random=random)
        self.format = format
        self._missing_data = missing_data

    def __iter__(self):
        return iter(self.columns.iter_rows())

    def __len__(self):
        return self.columns._num_rows

    def __getitem__(self, names):
        # this is funky, but a side-effect of construction allowing setting
        # prior to having assigned the index column
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

        # if a row index has been specified, via row_ids, we need to interpret
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
        row_ids = attr.pop("row_ids")
        result = self.__class__(**attr)
        for c in columns:
            result.columns[c] = self.columns[c][rows]

        if row_ids in result.columns:
            result.index_name = row_ids

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
        row_ids = kwargs.pop("row_ids")
        table = self.__class__(**kwargs)
        table.columns.__setstate__(data["data"])
        table.index_name = row_ids
        self.__dict__.update(table.__dict__)

    def __repr__(self):
        table, shape_info = self._get_repr_()
        result = "\n".join([str(table), shape_info])
        return result

    def __str__(self):
        return self.to_string(self.format)

    def _get_repr_(self):
        """returns a table for __repr__"""
        rn = self._repr_policy["random"]
        head = self._repr_policy["head"]
        tail = self._repr_policy["tail"]
        if head is None and tail is None:
            if self.shape[0] < 50:
                head = self.shape[0]
                tail = None
            else:
                head, tail = 5, 5
            self._repr_policy["head"] = head
            self._repr_policy["tail"] = tail

        shape_info = ""
        ellipsis = [["..."] * len(self.header)]
        if rn:
            indices = numpy.random.choice(self.shape[0], size=rn, replace=False)
            indices = list(sorted(indices))
            rows = self.array.take(indices, axis=0).tolist()
            shape_info = f"Random selection of {rn} rows"
        elif all([head, tail]):
            top = self[:head].tolist()
            bottom = self[-tail:].tolist()
            rows = top + ellipsis + bottom
        elif head:
            rows = self[:head].tolist()
        elif tail:
            rows = self[-tail:].tolist()
        else:
            rows = self.tolist()

        shape_info += f"\n{self.shape[0]:,} rows x {self.shape[1]:,} columns"
        kwargs = self._get_persistent_attrs()
        table = self.__class__(header=self.header, data=rows, **kwargs)
        table._column_templates.update(self._column_templates)
        return table, shape_info

    def _repr_html_(self, include_shape=True):
        """returns html, used by Jupyter"""
        base_colour = "rgba(161, 195, 209, {alpha})"

        def row_cell_func(val, row, col):
            colour = base_colour.format(alpha=0.25)
            try:
                float(val)
            except ValueError:
                is_numeric = False
            else:
                is_numeric = True

            if self.index_name and col == 0:
                style = f' style="background: {colour}; font-weight: 600;"'
            elif is_numeric:
                style = f' style="font-family: monospace !important;"'
            else:
                style = ""
            val = f"<td{style}>{val}</td>"
            return val

        table, shape_info = self._get_repr_()
        shape_info = f"<p>{shape_info}</p>"
        if not include_shape:
            shape_info = ""

        title, legend = table.title, table.legend
        # current rich_html does not provide a good mechanism for custom
        # formatting of titles, legends
        table.title, table.legend = None, None
        head_colour = base_colour.format(alpha=0.75)
        element_format = dict(
            thead=f'<thead style="background: {head_colour}; '
            'font-weight: bold; text-align: center;">'
        )
        html = table.to_rich_html(
            row_cell_func=row_cell_func, element_formatters=element_format
        )
        if title or legend:
            title = title or ""
            legend = legend or ""
            caption = (
                '<caption style="color: rgb(250, 250, 250); '
                'background: rgba(30, 140, 200, 1); align=top;">'
                f'<span style="font-weight: bold;">{title}</span>'
                f"<span>{legend}</span></caption>"
            )
            html = html.splitlines()
            html.insert(1, caption)
            html = "\n".join(html)
        html = html.splitlines()
        html.insert(
            1,
            "\n".join(
                [
                    "<style>",
                    "tr:last-child {border-bottom: 1px solid #000;} "
                    "tr > th {text-align: center !important;} tr > td {text-align: left !important;}",
                    "</style>",
                ]
            ),
        )
        html = "\n".join(["\n".join(html), shape_info])
        return html

    def _get_persistent_attrs(self):
        attrs = UnionDict(self._persistent_attrs.copy())
        return attrs

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

    def set_repr_policy(self, head=None, tail=None, random=0):
        """specify policy for repr(self)

        Parameters
        ----------

        - head: number of top rows to included in represented display
        - tail: number of bottom rows to included in represented display
        - random: number of rows to sample randomly (supercedes head/tail)
        """
        if not any([head, tail, random]):
            return
        if random:
            assert (
                type(random) == int and random > 0
            ), "random must be a positive integer"
            head = tail = None
        self._repr_policy = dict(head=head, tail=tail, random=random)

    @property
    def format(self):
        """the str display format"""
        return self._format

    @format.setter
    def format(self, new):
        """the str display format"""
        # setting the default format for str(self)
        if new.lower() in table_format.known_formats:
            new = new.lower()
        else:
            new = "simple"
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
        self._column_templates[column_head] = format_template

    def head(self, nrows=5):
        """displays top nrows"""
        repr_policy = self._repr_policy
        self._repr_policy = dict(head=nrows, tail=None, random=None)
        display(self)
        self._repr_policy = repr_policy

    def tail(self, nrows=5):
        """displays bottom nrows"""
        repr_policy = self._repr_policy
        self._repr_policy = dict(head=None, tail=nrows, random=None)
        display(self)
        self._repr_policy = repr_policy

    @property
    def index_name(self):
        if self._index_name and not self._template:
            self.columns.index_name = self._index_name
            self.index_name = self._index_name

        return self._index_name

    @index_name.setter
    def index_name(self, value):
        if value is None:
            return

        self.columns.index_name = value
        self._index_name = value
        self._template = DictArrayTemplate(self.columns[value])

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
        self, other, columns_self=None, columns_other=None, use_index=True, **kwargs,
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
            index, this will be used.

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

        # key is a tuple made from specified columns; data is the row index
        other_row_index = defaultdict(list)
        # subtable = other.columns.take_columns(columns_other)
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
        self, other, columns_self=None, columns_other=None, inner_join=True, **kwargs,
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
        if not isinstance(callback, Callable):
            data = subset
        else:
            data = subset.array

        num_columns = len(columns)
        match = not negate
        indices = numpy.array(
            [
                True
                if _callback(callback, row=row, num_columns=num_columns) == match
                else False
                for row in data
            ]
        )
        return indices

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
            data = list(tuple(e) for e in data)

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
        All tables must have the same columns.
        """
        if new_column:
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
        for table in table_series:
            assert set(table.columns.order) == columns, "columns don't match"
            if new_column:
                new_col.extend([table.title] * table.shape[0])
            data = table.columns.to_dict()
            for c, v in data.items():
                raw_data[c].extend(v)

        dtypes = {c: self.columns[c].dtype for c in self.columns}
        if new_column:
            columns = (new_column,) + self.columns.order
            raw_data[new_column] = new_col
            dtypes[new_column] = "<U15"
        else:
            columns = self.columns.order
        for c in columns:
            result.columns[c] = numpy.array(raw_data[c], dtype=dtypes[c])
        return result

    def get_columns(self, columns):
        """Return a Table with just columns"""
        if self.index_name:
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
        row_ids = attr.pop("row_ids")
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
        if not isinstance(callback, Callable):
            data = subset
        else:
            data = subset.array

        num_columns = len(columns)
        values = numpy.array(
            [_callback(callback, row=row, num_columns=num_columns) for row in data]
        )

        if dtype:
            values = numpy.array(values, dtype=dtype)

        result.columns[new_column] = values

        if row_ids in result.columns:
            result.index_name = row_ids

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
        if indices is None:
            data = self.array
        else:
            data = self[indices, :].array

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

    def sorted(self, columns=None, reverse=False, **kwargs):
        """Returns a new table sorted according to columns order.

        Parameters
        ----------
        columns
            column headings, their order determines the sort order.
        reverse
            column headings, these columns will be reverse sorted.

            Either can be provided as just a single string, or a series of
            strings.

        Notes
        -----
        If only reverse is provided, that order is used.
        """
        reverse = reverse if reverse else []
        if reverse and columns is None:
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
            dtype = self.columns[c].dtype.name
            if "int" in dtype or "float" in dtype or "complex" in dtype:
                func = _reverse_num
            else:
                func = _reverse_str
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

    def _formatted(self, missing_data=""):
        missing_data = missing_data or self._missing_data
        formatted = []
        for c in self.columns.order:
            data = self.columns[c]
            format_spec = self._column_templates.get(c, None)
            frmt, c = formatted_array(
                data,
                c,
                format_spec=format_spec,
                missing_data=missing_data,
                precision=self._digits,
            )
            formatted.append([c] + frmt)

        formatted = list([list(e) for e in zip(*formatted)])
        return formatted

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

        if format.lower() == "phylip":
            missing_data = "0.0000"
        else:
            missing_data = self._missing_data

        if format.lower() in ("tsv", "csv"):
            sep = sep or {"tsv": "\t", "csv": ","}[format.lower()]
            format = ""

        # convert self to a 2D list
        if format != "phylip":
            formatted_table = self._formatted()
        else:
            columns = [c for c in self.columns if c != self.index_name]
            table = self[:, columns]
            formatted_table = table._formatted(missing_data=missing_data)

        header = formatted_table.pop(0)
        args = (header, formatted_table, self.title, self.legend)
        if sep and format != "bedgraph":
            args = (header, formatted_table, None, None)
            return table_format.separator_format(*args + (sep,))
        elif format in ("rest", "rst"):
            return table_format.grid_table_format(*args)
        elif format in ("markdown", "md"):
            return table_format.markdown(header, formatted_table, **kwargs)
        elif format.endswith("tex"):
            caption = self.title or None
            legend = self.legend or None
            if concat_title_legend and (caption or legend):
                caption = " ".join([caption or "", legend or ""])
                caption = caption.strip()
                legend = None
            return table_format.latex(
                formatted_table, header, caption=caption, legend=legend, **kwargs
            )
        elif format == "html":
            return self.to_rich_html(**kwargs)
        elif format == "phylip":
            # need to eliminate row identifiers
            return table_format.phylip_matrix(formatted_table, header)
        else:
            return table_format.simple_format(
                *args + (self._max_width, self.index_name, borders, self.space)
            )

    def to_rich_html(
        self,
        row_cell_func=None,
        header_cell_func=None,
        element_formatters=None,
        merge_identical=False,
        compact=False,
    ):
        """returns just the table as html.

        Parameters
        ----------
        row_cell_func
            callback function that formats the row values. Must
            take the row value and coordinates (row index, column index).
        header_cell_func
            callback function that formats the column headings
            must take the header label value and coordinate
        element_formatters
            a dictionary of specific callback funcs for
            formatting individual html table elements.
            e.g. {'table': lambda x: '<table border="1" class="docutils">'}
        merge_identical
            cells within a row are merged to one span.

        """
        element_formatters = element_formatters or {}
        formatted_table = self.array.tolist()
        header, formatted_table = table_format.formatted_cells(
            formatted_table,
            self.header,
            digits=self._digits,
            column_templates=self._column_templates,
            missing_data=self._missing_data,
        )
        subtables = table_format.get_continuation_tables(
            header,
            formatted_table,
            identifiers=self.index_name,
            max_width=self._max_width,
        )
        tables = []
        title = self.title if self.title else ""
        if title:
            title = escape(title)
        legend = self.legend if self.legend else ""
        if legend:
            legend = escape(legend)
        for i, (h, t) in enumerate(subtables):
            # but we strip the cell spacing
            sh = [v.strip() for v in h]
            t = [[c.strip() for c in r] for r in t]

            if title and i == 0:
                st = element_formatters.get(
                    "caption", f'<span style="font-weight:bold">{title}</span>'
                )
            elif title:
                st = element_formatters.get(
                    "caption", f'<span style="font-weight:bold">continuation</span>'
                )
            else:
                st = None

            if legend and i == 0:
                title = f"{st} {legend}" if st else legend

            caption = st if st else None
            subtable = table_format.rich_html(
                t,
                row_cell_func=row_cell_func,
                header=sh,
                header_cell_func=header_cell_func,
                element_formatters=element_formatters,
                merge_identical=merge_identical,
                compact=compact,
                caption=caption,
            )
            tables.append(subtable)
        return "\n".join(tables)

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
        if self.index_name:
            row_ids = self.columns[self.index_name]
        else:
            row_ids = self.shape[0]
        template = DictArrayTemplate(row_ids, self.columns.order)
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

        rows = self.array.tolist()
        header, rows = table_format.formatted_cells(
            rows,
            self.header,
            digits=self._digits,
            column_templates=self._column_templates,
            missing_data=self._missing_data,
            center=False,
        )
        # we strip white space padding from header and cells
        header = [e.strip() for e in header]
        rows = [[e.strip() for e in row] for row in rows]
        rows = list(zip(*rows))
        if self.index_name:
            body_colour = ["white"] * self.shape[0]
            index_colour = ["rgba(161, 195, 209, 0.5)"] * self.shape[0]
            colours = [index_colour] + [body_colour[:] for i in range(self.shape[1])]
            rows[0] = [f"<b>{e}</b>" for e in rows[0]]
        else:
            colours = "white"

        tab = UnionDict(
            type="table",
            header=dict(
                values=[f"<b>{c}</b>" for c in header],
                fill=dict(color="rgba(161, 195, 209, 1)"),
                font=dict(size=font_size),
                align="center",
            ),
            cells=dict(values=rows, fill=dict(color=colours)),
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
        assert select_as_header in self.columns, (
            '"%s" not in table header' % select_as_header
        )

        if len(self.distinct_values(select_as_header)) != len(self):
            raise ValueError(f"not all '{select_as_header}' values unique")

        attr = self._get_persistent_attrs()
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
        file_suffix, compress_suffix = get_format_suffixes(filename)
        format = format or file_suffix
        compress = compress or compress_suffix is not None

        mode = mode or {"pickle": "wb"}.get(format, "w")

        if compress:
            if not filename.endswith(".gz"):
                filename = "%s.gz" % filename
            mode = "wt"

        outfile = open_(filename, mode)

        if format is None:
            # try guessing from filename suffix
            if compress:
                index = -2
            else:
                index = -1
            suffix = filename.split(".")
            if len(suffix) > 1:
                format = suffix[index]

        if format == "csv":
            sep = sep or ","
        elif format == "tsv":
            sep = sep or "\t"

        if writer:
            rows = self.tolist()
            rows.insert(0, self.header[:])
            rows = writer(rows, has_header=True)
            outfile.writelines("\n".join(rows))
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
            outfile.writelines(table + "\n")
        outfile.close()
