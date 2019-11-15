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
import warnings

from collections.abc import Callable
from xml.sax.saxutils import escape

import numpy

from cogent3.format import bedgraph
from cogent3.format import table as table_format
from cogent3.util.dict_array import DictArray
from cogent3.util.misc import get_format_suffixes, get_object_provenance, open_
from cogent3.util.union_dict import UnionDict


try:
    from pandas import DataFrame

    _pandas_available = True
except ImportError:
    _pandas_available = False

try:
    from IPython.display import display
except ImportError:
    display = lambda x: print(repr(x))


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Felix Schill"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
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


def convert2DDict(twoDdict, header=None, row_order=None):
    """Returns a 2 dimensional list.

    Parameters
    ----------
    twoDdict
        a 2 dimensional dict with top level keys corresponding to
        column headings, lower level keys correspond to row headings but are
        not preserved.
    header
        series with column headings. If not provided, the sorted top
        level dict keys are used.
    row_order
        a specified order to generate the rows.

    """
    if not header:
        header = list(twoDdict.keys())
        header.sort()

    if not row_order:  # we assume rows consistent across dict
        row_order = list(twoDdict[header[0]].keys())
        row_order.sort()

    # make twoD list
    table = []
    for row in row_order:
        string_row = []
        for column in header:
            string_row.append(twoDdict[column][row])
        table.append(string_row)
    return table


class _Header(list):
    """a convenience class for storing the header"""

    def __new__(cls, arg):
        n = list.__new__(cls, list(arg))
        return n

    def __setslice__(self, *args):
        """disallowed"""
        raise RuntimeError("Table header is immutable, use with_new_header")

    def __setitem__(self, *args):
        """disallowed"""
        raise RuntimeError("Table header is immutable, use with_new_header")


class Table(DictArray):
    def __init__(
        self,
        header=None,
        rows=None,
        row_order=None,
        digits=4,
        space=4,
        title="",
        missing_data="",
        max_width=1e100,
        row_ids=None,
        legend="",
        column_templates=None,
        dtype=None,
        data_frame=None,
        format="simple",
    ):
        """

        Parameters
        ----------
        header
            column headings
        rows
            a 2D dict, list or tuple. If a dict, it must have column
            headings as top level keys, and common row labels as keys in each
            column.
        row_order
            the order in which rows will be pulled from the twoDdict
        digits
            floating point resolution
        space
            number of spaces between columns or a string
        title
            as implied
        missing_data
            character assigned if a row has no entry for a column
        max_width
            maximum column width for printing
        row_ids
            if True, the 0'th column is used as row identifiers and keys
            for slicing.
        legend
            table legend
        column_templates
            dict of column headings
            or a function that will handle the formatting.
        dtype
            optional numpy array typecode.
        data_frame
            pandas DataFrame, Table will be created from this
        format
            output format when using str(Table)

        """
        if data_frame is not None and not _pandas_available:
            raise ValueError("data_frame provided when pandas not installed")
        elif data_frame is not None:
            if rows or header:
                warnings.warn(
                    "provided rows/header will be over ridden by " "DataFrame"
                )

            rows = data_frame.to_records(index=False).tolist()
            header = data_frame.columns.tolist()

        if type(header) == numpy.ndarray:
            header = header.tolist()

        if not header:
            raise ValueError("header must be provided to Table")
        elif rows is None:
            raise ValueError("rows cannot be None")

        if len(rows) == 0:
            rows = numpy.empty((0, len(header)))

        try:
            num_cols = len(header)
            assert num_cols > 0
            if type(rows) == numpy.ndarray:
                assert num_cols == rows.shape[1]
            elif type(rows) == dict:
                assert num_cols == len(rows)
            else:
                assert num_cols == len(rows[0])
        except (IndexError, TypeError, AssertionError):
            raise RuntimeError("header and rows must be provided to Table")

        header = [str(head) for head in header]
        if isinstance(rows, dict):
            rows = convert2DDict(rows, header=header, row_order=row_order)

        # if row_ids, we select that column as the row identifiers
        if row_ids is not None:
            identifiers = [row[0] for row in rows]
        else:
            identifiers = len(rows)

        if not dtype:
            dtype = "O"
        DictArray.__init__(self, rows, identifiers, header, dtype=dtype)

        # forcing all column headings to be strings
        self._header = _Header([str(head) for head in header])
        self._missing_data = missing_data

        # default title / legend to be empty strings
        self.title = str(title) if title else ""
        self.legend = str(legend) if legend else ""
        try:
            self.space = " " * space
        except TypeError:
            self.space = space
        self._digits = digits
        self._row_ids = row_ids
        self._max_width = max_width

        # some attributes are not preserved in any file format, so always based
        # on args
        self._column_templates = column_templates or {}

        self.format = format

        # define the repr() display policy
        random = 0
        if self.shape[0] < 50:
            head = self.shape[0]
            tail = None
        else:
            head, tail = 5, 5

        self._repr_policy = dict(head=tail, tail=tail, random=random)

    def __repr__(self):
        table, shape_info = self._get_repr_()
        result = "\n".join([str(table), shape_info])
        return result

    def _get_repr_(self):
        """returns a table for __repr__"""
        rn = self._repr_policy["random"]
        head = self._repr_policy["head"]
        tail = self._repr_policy["tail"]
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
        table = self.__class__(header=self.header, rows=rows, **kwargs)
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

            if self._row_ids and col == 0:
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

    def __str__(self):
        return self.to_string(self.format)

    def __getitem__(self, names):
        (index, remaining) = self.template.interpret_index(names)
        # if we have two integers, return a single value
        ints = [isinstance(idx, int) for idx in index]
        if len(ints) == 2 and min(ints):
            return self.array[index]
        new_index = list(index)
        for i, idx in enumerate(new_index):
            if isinstance(idx, int):
                new_index[i] = slice(idx, idx + 1, None)

        index = tuple(new_index)
        rows = self.array[index]
        result = None
        if len(index) > 1:
            header = numpy.asarray(self.header, dtype="O")[index[1:]]
        else:
            header = self.header
        if remaining is not None:
            kwargs = self._get_persistent_attrs()
            result = self.__class__(header, rows, **kwargs)
        return result

    def __getstate__(self):
        data = self._get_persistent_attrs()
        del data["column_templates"]
        data.update(dict(header=self.header, rows=self.tolist()))
        return data

    def __setstate__(self, data):
        limit_ids = data.pop("limit_ids", None)
        if limit_ids is not None:
            data["row_ids"] = limit_ids or False
        new = Table(**data)
        self.__dict__.update(new.__dict__)
        return self

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
    def header(self):
        """returns header value"""
        return self._header

    @header.setter
    def header(self, data):
        """disallowed"""
        raise RuntimeError("not allowed to set the header")

    @property
    def format(self):
        """the display format"""
        return self._format

    @format.setter
    def format(self, new):
        """the display format"""
        # setting the default format for str(self)
        if new.lower() in table_format.known_formats:
            new = new.lower()
        else:
            new = "simple"
        self._format = new

    def with_new_header(self, old, new, **kwargs):
        """returns a new Table with old header labels replaced by new

        Parameters
        ----------
        old
            the old column header(s). Can be a string or series of
            them.
        new
            the new column header(s). Can be a string or series of
            them.

        """
        if type(old) == str:
            old = [old]
            new = [new]

        assert len(old) == len(new), "Mismatched number of old/new labels"
        indices = list(map(self.header.index, old))
        new_header = list(self.header)
        for i in range(len(old)):
            new_header[indices[i]] = new[i]

        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header=new_header, rows=self.tolist(), **kw)

    def _get_persistent_attrs(self):
        kws = dict(
            row_ids=self._row_ids,
            title=self.title,
            legend=self.legend,
            digits=self._digits,
            space=self.space,
            max_width=self._max_width,
            missing_data=self._missing_data,
            column_templates=self._column_templates or None,
            format=self._format,
        )
        return kws

    def format_column(self, column_head, format_template):
        """Provide a formatting template for a named column.

        Parameters
        ----------
        column_head
            the column label.
        format_template
            string formatting template or a function that
            will handle the formatting.

        """
        assert column_head in self.header, "Unknown column heading %s" % column_head

        self._column_templates[column_head] = format_template

    def to_string(self, format="", borders=True, sep=None, center=False, **kwargs):
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
        center
            content is centered in the column, default is right
            justified

        NOTE: If format is bedgraph, assumes that column headers are chrom,
        start, end, value. In that order!
        """
        if format.lower() == "phylip":
            missing_data = "%.4f" % 0.0
        else:
            missing_data = self._missing_data

        if format.lower() in ("tsv", "csv"):
            sep = sep or {"tsv": "\t", "csv": ","}[format.lower()]
            format = ""

        # convert self to a 2D list
        formatted_table = self.array.tolist()
        if format != "bedgraph":
            header, formatted_table = table_format.formatted_cells(
                formatted_table,
                self.header,
                digits=self._digits,
                column_templates=self._column_templates,
                missing_data=missing_data,
                center=center,
            )
            args = (header, formatted_table, self.title, self.legend)
        if sep and format != "bedgraph":
            return table_format.separator_format(*args + (sep,))
        elif format in ("rest", "rst"):
            return table_format.grid_table_format(*args)
        elif format in ("markdown", "md"):
            return table_format.markdown(header, formatted_table, **kwargs)
        elif format.endswith("tex"):
            caption = None
            if self.title or self.legend:
                caption = " ".join([self.title or "", self.legend or ""])
            return table_format.latex(
                formatted_table, header, caption=caption, **kwargs
            )
        elif format == "html":
            return self.to_rich_html(**kwargs)
        elif format == "phylip":
            # need to eliminate row identifiers
            formatted_table = [row[self._row_ids :] for row in formatted_table]
            header = header[self._row_ids :]
            return table_format.phylip_matrix(formatted_table, header)
        elif format == "bedgraph":
            assert self.shape[1] == 4, "bedgraph format is for 4 column tables"
            # assuming that header order is chrom, start, end, val
            formatted_table = bedgraph.bedgraph(self.sorted().array.tolist(), **kwargs)
            return formatted_table
        else:
            return table_format.simple_format(
                *args + (self._max_width, self._row_ids, borders, self.space)
            )

    def to_rich_html(
        self,
        row_cell_func=None,
        header_cell_func=None,
        element_formatters=None,
        merge_identical=False,
        compact=False,
    ):
        """returns just the table html code.

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
            identifiers=self._row_ids,
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

            legend = self.legend if self.legend else ""
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
        """Write table to filename in the specified format. If a format is not
        specified, it attempts to use a filename suffix. Note if a sep argument
        is provided, unformatted values are written to file in order to
        preserve numerical accuracy.

        Parameters
        ----------
        mode
            file opening mode
        format
            Valid formats are those of the to_string method plus
            pickle. Will try and guess from filename if not specified.
        writer
            a function for formatting the data for output.
        sep
            a character delimiter for fields.
        compress
            if True, gzips the file and appends .gz to the
            filename (if not already added).

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

    def appended(self, new_column, *tables, **kwargs):
        """Append an arbitrary number of tables to the end of this one.
        Returns a new table object. Optional keyword arguments to the new
        tables constructor may be passed.

        Parameters
        ----------
        new_column
            provide a heading for the new column, each tables
            title will be placed in it. If value is false, the result is no
            additional column.

        """
        # default title is no title
        kwargs["title"] = kwargs.get("title", "")
        # convert series of tables
        if isinstance(tables[0], tuple) or isinstance(tables[0], list):
            tables = tuple(tables[0])
        # for each table, determine it's number of rows and create an
        # equivalent length vector of its title
        if new_column:
            header = [new_column] + self.header
        else:
            header = self.header

        big_twoD = ()
        table_series = (self,) + tables
        for table in table_series:
            # check compatible tables
            assert (
                self.header == table.header
            ), "Inconsistent tables -- column headings are not the same."
            new_twoD = []
            for row in table:
                if new_column:
                    new_twoD.append([table.title] + row.to_array().tolist())
                else:
                    new_twoD.append(row.to_array().tolist())
            new_twoD = tuple(new_twoD)
            big_twoD += new_twoD
        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header, big_twoD, **kw)

    def tolist(self, columns=None):
        """Returns raw data as a 1D or 2D list of rows from columns. If one
        column, its a 1D list.

        Parameters
        ----------
        columns
            if None, all data are returned

        """

        if columns is None:
            return self.array.tolist()

        if isinstance(columns, str):
            # assumes all column headings are strings.
            columns = (columns,)

        column_indices = list(map(self.header.index, columns))
        result = self.array.take(column_indices, axis=1)

        if len(columns) == 1:
            result = result.flatten()

        return result.tolist()

    def to_dataframe(self, categories=None):
        """returns pandas DataFrame instance

        Parameters
        ----------
        categories
            converts these columns to category dtype in the data
            frame. Note, categories are not ordered.

        """
        if not _pandas_available:
            raise ImportError("pandas not installed")

        index = None if self._row_ids is None else self.template.names[0]
        data = dict(zip(self.header, self.to_array().T.tolist()))
        df = DataFrame(data=data, index=index)
        if categories is not None:
            categories = [categories] if type(categories) == str else categories
            df = df.astype({n: "category" for n in categories})

        return df

    def _callback(self, callback, row, columns=None, num_columns=None):
        if isinstance(callback, Callable):
            row_segment = row.take(columns)
            if num_columns == 1:
                row_segment = row_segment[0]
            return callback(row_segment)
        else:
            return eval(callback, {}, row)

    def filtered(self, callback, columns=None, **kwargs):
        """Returns a sub-table of rows for which the provided callback
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

        if columns:
            num_columns = len(columns)
        else:
            num_columns = None

        row_indexes = []
        if not isinstance(callback, Callable):
            data = self
            cols = columns
        else:
            data = self.array
            cols = list(map(self.header.index, columns))

        for rdex, row in enumerate(data):
            if self._callback(callback, row, cols, num_columns):
                row_indexes.append(rdex)

        sub_set = numpy.take(self, row_indexes, 0)

        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header=self.header, rows=sub_set, **kw)

    def filtered_by_column(self, callback, **kwargs):
        """Returns a table with columns identified by callback

        Parameters
        ----------
        callback
            A function which takes the columns delimited
            by columns and returns True/False, or a string representing valid
            python code to be evaluated.

        """
        data = self.array.transpose()
        column_indices = []
        append = column_indices.append
        for index, row in enumerate(data):
            if callback(row):
                append(index)
        columns = numpy.take(self.header, column_indices)
        return self.get_columns(columns, **kwargs)

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

        if columns:
            num_columns = len(columns)
        else:
            num_columns = None

        count = 0
        if not isinstance(callback, Callable):
            data = self
            cols = columns
        else:
            data = self.array
            cols = list(map(self.header.index, columns))

        for row in data:
            if self._callback(callback, row, cols, num_columns):
                count += 1
        return count

    def sorted(self, columns=None, reverse=None, **kwargs):
        """Returns a new table sorted according to columns order.

        If only reverse is provided, that order is used.

        Parameters
        ----------
        columns
            column headings, their order determines the sort order.
        reverse
            column headings, these columns will be reverse sorted.

            Either can be provided as just a single string, or a series of
            strings.
        """

        if reverse and columns is None:
            columns = reverse

        if columns is None:
            columns = self.header
        elif isinstance(columns, str):
            columns = [columns]

        indices = [self.header.index(col) for col in columns]

        if not reverse:
            reverse_indices = []
        else:
            if type(reverse) == str:
                reverse = [reverse]
            reverse_indices = []
            for index, header_index in enumerate(indices):
                col = self.header[header_index]
                if col in reverse:
                    reverse_indices += [index]

        reverse_indices = numpy.array(reverse_indices)

        dtypes = [(self.header[i], self.array.dtype) for i in indices]

        # applying the decorate-sort-undecorate approach
        aux_list = self.array.take(indices, axis=1)

        # we figure out the casting funcs for any reversed elements
        for index in reverse_indices:
            val = aux_list[0, index]
            try:
                val = val.translate(_reversed_chrs)
                func = _reverse_str
            except AttributeError:
                func = _reverse_num
            func = numpy.vectorize(func)
            aux_list[:, index] = func(aux_list[:, index])

        aux_list = numpy.rec.fromarrays(aux_list.copy().T, dtype=dtypes)
        indices = aux_list.argsort()
        new_twoD = self.array.take(indices, axis=0)

        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header=self.header, rows=new_twoD, **kw)

    def get_columns(self, columns, **kwargs):
        """Return a Table with just columns"""
        # check whether we have integer columns

        if isinstance(columns, str):
            columns = [columns]

        is_int = min([isinstance(val, int) for val in columns])
        indexes = []
        if is_int:
            indexes = columns
        else:
            indexes = [self.header.index(head) for head in columns]

        if self._row_ids:
            # we disallow reordering of identifiers, and ensure they are only
            # presented once
            for val in range(self._row_ids):
                try:
                    indexes.remove(val)
                except ValueError:
                    pass
            indexes = list(range(self._row_ids)) + indexes

        columns = numpy.take(numpy.asarray(self.header, dtype="O"), indexes)
        new = numpy.take(self.array, indexes, axis=1)

        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header=columns, rows=new, **kw)

    def with_new_column(self, new_column, callback, columns=None, **kwargs):
        """Returns a new table with an additional column, computed using
        callback.

        Parameters
        ----------
        new_column
            new column heading
        columns
            the columns whose values determine whether a row is to
            be included.
        callback
            Can be a function, which takes the sub
            by columns and returns True/False, or a string representing valid
            python code to be evaluated.

        """
        if new_column in self.header:
            raise AssertionError("column '%s' already exists" % new_column)

        if isinstance(columns, str):
            columns = (columns,)

        if columns is not None:
            num_columns = len(columns)
        else:
            num_columns = None

        if not isinstance(callback, Callable):
            data = self
            cols = columns
        else:
            data = self.array
            cols = list(map(self.header.index, columns))

        twoD = [
            list(row) + [self._callback(callback, row, cols, num_columns)]
            for row in data
        ]

        kw = self._get_persistent_attrs()
        kw.update(kwargs)
        return Table(header=self.header + [new_column], rows=twoD, **kw)

    def distinct_values(self, column):
        """returns the set of distinct values for the named column(s)"""
        columns = [column, [column]][type(column) == str]
        data = self.tolist(column)

        if len(columns) > 1:
            data = [tuple(row) for row in data]

        return set(data)

    def joined(
        self,
        other_table,
        columns_self=None,
        columns_other=None,
        inner_join=True,
        **kwargs,
    ):
        """returns a new table containing the join of this table and
        other_table. Default behaviour is the natural inner join. Checks for
        equality in the specified columns (if provided) or all columns; a
        combined row is included in the output if all indices match exactly. A
        combined row contains first the row of this table, and then columns
        from the other_table that are not key columns (i.e. not specified in
        columns_other). The order (of self, then other)
        is preserved. The column headers of the output are made unique by
        replacing the headers of other_table with
        <other_table.Title>_<other_table.header>.

        Parameters
        ----------
        other_table
            A table object which will be joined with this
            table. other_table must have a title.
        columns_self, columns_other
            indices of key columns that will
            be compared in the join operation. Can be either column index,
            or a string matching the column header. The order matters, and
            the dimensions of columns_self and columns_other have to match.
            A row will be included in the output iff
            self[row][columns_self[i]]==other_table[row][columns_other[i]]
            for all i
        inner_join
            if False, the outer join of the two tables is
            returned.

        """
        other_title = other_table.title if other_table.title else "right"
        if self.title == other_title:
            raise RuntimeError("Cannot join if a table.Title's are equal")

        columns_self = [columns_self, [columns_self]][type(columns_self) == str]
        columns_other = [columns_other, [columns_other]][type(columns_other) == str]
        if not inner_join:
            assert columns_self is None and columns_other is None, (
                "Cannot " "specify column indices for an outer join"
            )
            columns_self = []
            columns_other = []

        if columns_self is None and columns_other is None:
            # we do the natural inner join
            columns_self = []
            columns_other = []
            for col_head in self.header:
                if col_head in other_table.header:
                    columns_self.append(self.header.index(col_head))
                    columns_other.append(other_table.header.index(col_head))
        elif columns_self is None or columns_other is None:
            # the same column labels will be used for both tables
            columns_self = columns_self or columns_other
            columns_other = columns_self or columns_other
        elif len(columns_self) != len(columns_other):
            raise RuntimeError(
                "Error during table join: key columns have " "different dimensions!"
            )

        # create new 2d list for the output
        joined_table = []

        # resolve column indices from header, if necessary
        columns_self_indices = []
        columns_other_indices = []
        for col in columns_self:
            if type(col) == int:
                columns_self_indices.append(col)
            else:
                columns_self_indices.append(self.header.index(col))

        for col in columns_other:
            if type(col) == int:
                columns_other_indices.append(col)
            else:
                columns_other_indices.append(other_table.header.index(col))
        # create a mask of which columns of the other_table will end up in the
        # output
        output_mask_other = []
        for col in range(0, len(other_table.header)):
            if not (col in columns_other_indices):
                output_mask_other.append(col)
        # use a dictionary for the key lookup
        # key dictionary for other_table.
        # key is a tuple made from specified columns; data is the row index
        # for lookup...
        key_lookup = {}
        row_index = 0
        for row in other_table:
            # insert new entry for each row
            key = tuple([row[col] for col in columns_other_indices])
            if key in key_lookup:
                key_lookup[key].append(row_index)
            else:
                key_lookup[key] = [row_index]
            row_index = row_index + 1

        for this_row in self:
            # assemble key for query of other_table
            key = tuple([this_row[col] for col in columns_self_indices])
            if key in key_lookup:
                for output_row_index in key_lookup[key]:
                    other_row = [
                        other_table[output_row_index, c] for c in output_mask_other
                    ]
                    joined_table.append(list(this_row) + other_row)

        new_header = self.header + [
            other_title + "_" + other_table.header[c] for c in output_mask_other
        ]
        return Table(header=new_header, rows=joined_table, **kwargs)

    def summed(self, indices=None, col_sum=True, strict=True, **kwargs):
        """returns the sum of numerical values for column(s)/row(s)

        Parameters
        ----------
        indices
            column name(s) or indices or row indices
        col_sum
            sums values in the indicated column, the default. If
            False, returns the row sum.
        strict
            if False, ignores cells with non
            column/row.

        """

        all = indices is None

        if type(indices) == str:
            assert col_sum, "Must use row integer indices"
            indices = self.header.index(indices)
        elif type(indices) == int:  # a single index
            indices = [indices]
        elif not all:
            raise RuntimeError("unknown indices type: %s" % indices)

        if not all:
            vals = self.array.take([indices], axis=[0, 1][col_sum]).flatten()
            if strict:
                result = vals.sum()
            else:
                result = sum(v for v in vals if type(v) != str)
        else:
            # a multi-rowed result
            if col_sum:
                vals = self.array
            else:
                vals = self.array.transpose()

            if strict:
                result = vals.sum(axis=0).tolist()
            else:
                result = []
                append = result.append
                # we need to iterate over the elements to be summed, so we
                # have to transpose
                for row in vals.transpose():
                    try:
                        num = row.sum()
                    except TypeError:
                        num = sum(r for r in row if type(r) != str)
                    append(num)

        return result

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

        if denominator_func:
            data = self.array
            if not by_row:
                data = data.transpose()
            denominators = [denominator_func(row) for row in data]
        else:
            denominators = self.summed(col_sum=not by_row)

        if by_row:
            values = self.array
        else:
            values = self.array.transpose()

        rows = [values[i] / denom for i, denom in enumerate(denominators)]
        rows = numpy.array(rows)

        if not by_row:
            rows = rows.transpose()

        return Table(header=self.header, rows=rows, **kwargs)

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
        select_as_header = select_as_header or self.header[0]
        assert select_as_header in self.header, (
            '"%s" not in table header' % select_as_header
        )

        raw_data = self.tolist()
        raw_data.insert(0, self.header)
        transposed = numpy.array(raw_data, dtype="O")
        transposed = transposed.transpose()

        # indices for the header and non header rows
        header_index = self.header.index(select_as_header)

        data_indices = list(range(0, header_index)) + list(
            range(header_index + 1, len(transposed))
        )

        header = list(numpy.take(transposed, [header_index], axis=0)[0])
        # [1:] slice excludes old name
        header = [new_column_name] + header[1:]
        rows = numpy.take(transposed, data_indices, axis=0)
        return Table(header=header, rows=rows, **kwargs)

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
        if self._row_ids:
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

    def to_rich_dict(self):
        data = self.__getstate__()
        data["type"] = get_object_provenance(self)
        data["version"] = __version__
        return data

    def to_json(self):
        return json.dumps(self.to_rich_dict())
