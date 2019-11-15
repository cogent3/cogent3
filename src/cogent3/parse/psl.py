#!/usr/bin/env python
"""Parser for PSL format (default output by blat).
   Compatible with blat v.34
"""

from cogent3.parse.table import ConvertFields
from cogent3.util.table import Table


__author__ = "Gavin Huttley, Anuj Pahwa"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Gavin Huttley", "Anuj Pahwa"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Development"


def make_header(lines):
    """returns one header line from multiple header lines"""
    lengths = list(map(len, lines))
    max_length = max(lengths)
    for index, line in enumerate(lines):
        if lengths[index] != max_length:
            for i in range(lengths[index], max_length):
                line.append("")

    header = []
    for t, b in zip(*lines):
        if t.strip().endswith("-"):
            c = t.strip() + b
        else:
            c = " ".join([t.strip(), b.strip()])
        header += [c.strip()]
    return header


int_series = lambda x: list(map(int, x.replace(",", " ").split()))

row_converter = ConvertFields(
    [(i, int) for i in range(8)]
    + [(i, int) for i in range(10, 13)]
    + [(i, int) for i in range(14, 18)]
    + [(i, int_series) for i in range(18, 21)]
)


def MinimalPslParser(data, row_converter=row_converter):
    """returns version, header and rows from data"""
    if type(data) == str:
        data = open(data)

    psl_version = None
    header = None
    rows = []

    for record in data:
        if psl_version is None:
            assert "psLayout version" in record
            psl_version = record.strip()
            yield psl_version
            continue

        if not record.strip():
            continue

        if header is None and record[0] == "-":
            header = make_header(rows)
            yield header
            rows = []
            continue

        rows += [record.rstrip().split("\t")]
        if header is not None:
            yield row_converter(rows[0])
            rows = []

    try:
        data.close()
    except AttributeError:
        pass


def PslToTable(data):
    """converts psl format to a table"""
    parser = MinimalPslParser(data)
    version = next(parser)
    header = next(parser)
    rows = [row for row in parser]
    table = Table(header=header, rows=rows, title=version)
    return table
