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
__version__ = "2019.10.24a"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

from pathlib import Path

from cogent3.util.misc import open_


def gff_parser(f):
    """parses a gff file
    Parameters
    -----------
    f
        accepts string path or pathlib.Path or file-like object (e.g. StringIO)
    """

    f = f if not isinstance(f, Path) else str(f)
    if isinstance(f, str):
        with open_(f) as infile:
            yield from _parse_gff(infile)
    else:
        yield from _parse_gff(f)


def _parse_gff(f):
    """parses a gff file"""

    gff3_header = "gff-version 3"
    if isinstance(f, list):
        gff3 = f and gff3_header in f[0]
    else:
        gff3 = gff3_header in f.readline()
        f.seek(0)

    for line in f:
        # comments and blank lines
        if "#" in line:
            (line, comments) = line.split("#", 1)
        else:
            comments = None
        line = line.strip()
        if not line:
            continue

        cols = line.split("\t")
        # the final column (attributes) may be empty
        if len(cols) == 8:
            cols.append("")
        assert len(cols) == 9, len(line)
        (seqid, source, type, start, end, score, strand, phase, attributes) = cols

        # adjust for 0-based indexing
        (start, end) = (int(start) - 1, int(end))
        # start is always meant to be less than end in GFF
        # features that extend beyond sequence have negative indices
        if start < 0 or end < 0:
            start, end = abs(start), abs(end)
            if start > end:
                start, end = end, start
        # reverse indices when the feature is on the opposite strand
        if strand == "-":
            (start, end) = (end, start)

        if gff3:
            attribute_parser = parse_attributes_gff3
        else:
            attribute_parser = parse_attributes_gff2
        attributes = attribute_parser(attributes)

        yield (
            seqid,
            source,
            type,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
            comments,
        )


def parse_attributes_gff2(attributes):
    """Returns a dict with name and info keys"""
    name = attributes[attributes.find('"') + 1 :]
    if '"' in name:
        name = name[: name.find('"')]
    attr_dict = {"Name": name, "Info": attributes}
    return attr_dict


def parse_attributes_gff3(attributes):
    """Returns a dictionary containing all the attributes"""
    attributes = attributes.strip(";")
    attributes = attributes.split(";")
    attributes = dict(t.split("=") for t in attributes)
    return attributes


def gff_label(attributes):
    """Returns an identifier from the attributes"""
    if "ID" in attributes.keys():
        return attributes["ID"]
    elif "Name" in attributes.keys():
        return attributes["Name"]
    else:
        return str(attributes)
