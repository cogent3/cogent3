#!/usr/bin/env python

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Matthew Wakefield",
    "Gavin Huttley",
    "Christopher Bradley",
]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

from io import StringIO
from pathlib import Path

from cogent3.util.misc import open_


def gff_parser(f):
    """delegates to the correct gff_parser based on the version"""
    f = f if not isinstance(f, Path) else str(f)
    if isinstance(f, str):
        with open_(f) as infile:
            yield from gff2_parser(infile)
    elif isinstance(f, StringIO):
        yield from gff2_parser(f)
    else:
        raise TypeError


def gff2_parser(f):
    assert not isinstance(f, str)
    for line in f:
        # comments and blank lines
        if "#" in line:
            (line, comments) = line.split("#", 1)
        else:
            comments = None
        line = line.strip()
        if not line:
            continue

        # parse columns
        cols = line.split("\t")
        if len(cols) == 8:
            cols.append("")
        assert len(cols) == 9, line
        (seqname, source, feature, start, end, score, strand, frame, attributes) = cols

        # adjust for python 0-based indexing etc.
        (start, end) = (int(start) - 1, int(end))
        # start is always meant to be less than end in GFF
        # and in v 2.0, features that extend beyond sequence have negative
        # indices
        if start < 0 or end < 0:
            start, end = abs(start), abs(end)
            if start > end:
                start, end = end, start

        # but we use reversal of indices when the feature is on the opposite
        # strand
        if strand == "-":
            (start, end) = (end, start)

        # should parse attributes too
        yield (
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            frame,
            attributes,
            comments,
        )


def parse_attributes(attribute_string):
    """Returns region of attribute string between first pair of double quotes"""
    attribute_string = attribute_string[attribute_string.find('"') + 1 :]
    if '"' in attribute_string:
        attribute_string = attribute_string[: attribute_string.find('"')]
    return attribute_string
