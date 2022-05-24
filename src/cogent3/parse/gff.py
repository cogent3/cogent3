#!/usr/bin/env python
from pathlib import Path

from cogent3.util.io import open_


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Matthew Wakefield",
    "Gavin Huttley",
    "Christopher Bradley",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def gff_parser(f):
    """parses a gff file
    Parameters
    -----------
    f
        accepts string path or pathlib.Path or file-like object (e.g. StringIO)

    Returns
    -------
    dict
        contains each of the 9 parameters specified by gff3, and comments.
    """

    # calling a separate function to ensure file closes correctly
    f = f if not isinstance(f, Path) else str(f)
    if isinstance(f, str):
        with open_(f) as infile:
            yield from _gff_parser(infile)
    else:
        yield from _gff_parser(f)


def _gff_parser(f):
    """parses a gff file"""

    gff3_header = "gff-version 3"
    if isinstance(f, list):
        gff3 = f and gff3_header in f[0]
    else:
        gff3 = gff3_header in f.readline()
        f.seek(0)

    if gff3:
        attribute_parser = parse_attributes_gff3
    else:
        attribute_parser = parse_attributes_gff2

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
        seqid, source, type_, start, end, score, strand, phase, attributes = cols

        # adjust for 0-based indexing
        start, end = int(start) - 1, int(end)
        # start is always meant to be less than end in GFF
        # features that extend beyond sequence have negative indices
        if start < 0 or end < 0:
            start, end = abs(start), abs(end)
            if start > end:
                start, end = end, start
        # reverse indices when the feature is on the opposite strand
        if strand == "-":
            (start, end) = (end, start)

        # all attributes have an "ID" but this may not be unique
        attributes = attribute_parser(attributes, (start, end))

        rtn = {
            "SeqID": seqid,
            "Source": source,
            "Type": type_,
            "Start": start,
            "End": end,
            "Score": score,
            "Strand": strand,
            "Phase": phase,
            "Attributes": attributes,
            "Comments": comments,
        }
        yield rtn


def parse_attributes_gff2(attributes, span):
    """Returns a dict with name and info keys"""
    name = attributes[attributes.find('"') + 1 :]
    if '"' in name:
        name = name[: name.find('"')]
    return {"ID": name, "Info": attributes}


def parse_attributes_gff3(attributes, span):
    """Returns a dictionary containing all the attributes"""
    attributes = attributes.strip(";")
    attributes = attributes.split(";")
    attributes = dict(t.split("=") for t in attributes) if attributes[0] else {}
    if "Parent" in attributes:
        # There may be multiple parents
        if "," in attributes["Parent"]:
            attributes["Parent"] = attributes["Parent"].split(",")
        else:
            attributes["Parent"] = [attributes["Parent"]]
    attributes["ID"] = attributes.get("ID", "")
    return attributes
