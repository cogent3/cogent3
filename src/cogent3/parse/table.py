#!/usr/bin/env python

import csv

from collections.abc import Callable

from cogent3.util.misc import open_

from .record_finder import is_empty


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class ConvertFields(object):
    """converter for input data to Table"""

    def __init__(self, conversion, by_column=True):
        """handles conversions of columns or lines

        Parameters
        ----------
        by_column
            conversion will by done for each column, otherwise
            done by entire line

        """
        super(ConvertFields, self).__init__()
        self.conversion = conversion
        self.by_column = by_column

        self._func = self.convert_by_columns

        if not self.by_column:
            assert isinstance(
                conversion, Callable
            ), "conversion must be callable to convert by line"
            self._func = self.convert_by_line

    def convert_by_columns(self, line):
        """converts each column in a line"""
        for index, cast in self.conversion:
            line[index] = cast(line[index])
        return line

    def convert_by_line(self, line):
        """converts each column in a line"""
        return self.conversion(line)

    def __call__(self, *args, **kwargs):
        return self._func(*args, **kwargs)


def SeparatorFormatParser(
    with_header=True,
    converter=None,
    ignore=None,
    sep=",",
    strip_wspace=True,
    limit=None,
    **kw,
):
    """Returns a parser for a delimited tabular file.

    Parameters
    ----------
    with_header
        when True, first line is taken to be the header. Not
        passed to converter.
    converter
        a callable that returns a correctly formatted line.
    ignore
        lines for which ignore returns True are ignored. White
        lines are always skipped.
    sep
        the delimiter deparating fields.
    strip_wspace
        removes redundant white
    limit
        exits after this many lines

    """
    if ignore is None:  # keep all lines
        ignore = lambda x: False

    by_column = getattr(converter, "by_column", True)

    def callable(lines):
        num_lines = 0
        header = None
        for line in lines:
            if is_empty(line):
                continue

            line = line.strip("\n").split(sep)
            if strip_wspace and by_column:
                line = [field.strip() for field in line]

            if with_header and not header:
                header = True
                yield line
                continue

            if converter:
                line = converter(line)

            if ignore(line):
                continue

            yield line

            num_lines += 1
            if limit is not None and num_lines >= limit:
                break

    return callable


def load_delimited(
    filename,
    header=True,
    delimiter=",",
    with_title=False,
    with_legend=False,
    limit=None,
):
    if limit is not None:
        limit += 1  # don't count header line

    f = open_(filename)

    reader = csv.reader(f, dialect="excel", delimiter=delimiter)
    if with_title:
        title = "".join(next(reader))
    else:
        title = ""

    rows = []
    num_lines = 0
    for row in reader:
        rows.append(row)
        num_lines += 1
        if limit is not None and num_lines >= limit:
            break
    f.close()
    if header:
        header = rows.pop(0)
    else:
        header = None
    if with_legend:
        legend = "".join(rows.pop(-1))
    else:
        legend = ""
    # now do type casting in the order int, float, default is string
    return header, rows, title, legend
