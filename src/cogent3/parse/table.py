import csv
import pathlib

from cogent3.util.io import open_

from .record_finder import is_empty


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class FilteringParser:
    """A parser for a delimited tabular file that returns records matching a condition."""

    def __init__(
        self,
        row_condition=None,
        negate=False,
        with_header=True,
        columns=None,
        sep=",",
        limit=None,
    ):
        """
        Parameters
        ----------
        row_condition : callable
            callback that takes an entire line (except header) and returns True/False.
            A line is kept if condition(line) is True.
        negate : bool
            A line is kept if condition(line) is False.
        with_header : bool
            when True, first line is taken to be the header. Not
            passed to converter.
        columns
            series of indices, or column names, to return. If column names provided, with_header must be true.
        sep : str
            the delimiter separating fields.
        strip_wspace : bool
            removes redundant white
        limit : int
            exits after this many lines

        Notes
        -----
        The line elements are strings.
        """
        self.with_header = with_header
        columns = [columns] if isinstance(columns, (int, str)) else columns
        if columns is not None and isinstance(columns[0], str) and not with_header:
            raise ValueError("with_header must be True for columns with str values")

        self.columns = columns

        self.condition = row_condition
        self.negate = negate
        self.sep = sep
        self.limit = limit

    def _column_names_to_indices(self, header):
        if isinstance(self.columns[0], int):
            return

        absent = set(self.columns) - set(header)
        if absent:
            raise ValueError(f"columns not present in header {absent}")

        indices = [header.index(c) for c in self.columns]
        self.columns = indices

    def __call__(self, lines):
        """a generator that yields individual lines processed according to the
        provided conditions

        Parameters
        ----------
        lines: path or iterable
            If file path, handles file open and close. Will expand user
            component (i.e. '~/') of path.

        Notes
        -----
        Elements within a row are strings
        """
        input_from_path = False
        if isinstance(lines, (str, pathlib.Path)):
            path = pathlib.Path(lines).expanduser()
            input_from_path = path.exists()

            if input_from_path:
                lines = open_(path)

        num_lines = 0
        header = None
        match = not self.negate
        for line in lines:
            if is_empty(line):
                continue

            line = line.split(self.sep)
            line = [e.strip() for e in line]
            if header is None and self.with_header:
                header = True
                if self.columns:
                    self._column_names_to_indices(line)
                    line = [line[i] for i in self.columns]
                yield line
                continue

            if self.columns:
                line = [line[i] for i in self.columns]

            if self.condition and self.condition(line) != match:
                continue

            yield line

            num_lines += 1
            if self.limit is not None and num_lines >= self.limit:
                break

        if input_from_path:
            lines.close()


def load_delimited(
    filename,
    header=True,
    sep=",",
    with_title=False,
    with_legend=False,
    limit=None,
    **kwargs,
):
    """
    basic processing of tabular data

    Parameters
    ----------
    filename: Path
        path to delimited file (can begin with ~)
    header: bool
        whether the first line of the file (after the title, if present) is a header
    sep: str
        the character separating columns
    with_title: bool
        whether the first line of the file is a title
    with_legend: bool
        whether the last line of the file is a legend
    limit: int
        maximum number of lines to read from the file

    Returns
    -------
    header, rows, title, legend

    Notes
    -----
    All row values remain as strings.
    """
    if limit is not None and header:
        limit += 1  # don't count header line

    with open_(filename) as f:
        reader = csv.reader(f, dialect="excel", delimiter=sep)
        title = "".join(next(reader)) if with_title else ""
        rows = []
        num_lines = 0
        for row in reader:
            rows.append(row)
            num_lines += 1
            if limit is not None and num_lines >= limit:
                break

    header = rows.pop(0) if header else None
    legend = "".join(rows.pop(-1)) if with_legend else ""
    return header, rows, title, legend
