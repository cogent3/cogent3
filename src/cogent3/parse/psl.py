#!/usr/bin/env python
"""Parser for PSL format (default output by blat).
Compatible with blat v.34
"""

import contextlib

from cogent3.util.table import Table


def make_header(lines):
    """returns one header line from multiple header lines"""
    lengths = list(map(len, lines))
    max_length = max(lengths)
    for index, line in enumerate(lines):
        if lengths[index] != max_length:
            for _i in range(lengths[index], max_length):
                line.append("")

    header = []
    for t, b in zip(*lines, strict=False):
        c = t.strip() + b if t.strip().endswith("-") else f"{t.strip()} {b.strip()}"
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

    with contextlib.suppress(AttributeError):
        data.close()


def PslToTable(data):
    """converts psl format to a table"""
    parser = MinimalPslParser(data)
    version = next(parser)
    header = next(parser)
    rows = list(parser)
    return Table(header=header, data=rows, title=version)
