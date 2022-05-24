#!/usr/bin/env python
"""Parser for PSL format (default output by blat).
   Compatible with blat v.34
"""

from cogent3.util.table import Table


__author__ = "Gavin Huttley, Anuj Pahwa"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Gavin Huttley", "Anuj Pahwa"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
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


def MinimalPslParser(data):
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
            yield rows[0]
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
    return Table(header=header, data=rows, title=version)
