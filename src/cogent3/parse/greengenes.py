#!/usr/bin/env python

"""Parse the Greengenes formatted sequence data records

The script is intended to be used with the following input:
http://greengenes.lbl.gov/Download/Sequence_Data/Greengenes_format/greengenes16SrRNAgenes.txt.gz
"""

from cogent3.parse.record import DelimitedSplitter, GenericRecord
from cogent3.parse.record_finder import DelimitedRecordFinder


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"


def make_ignore_f(start_line):
    """Make an ignore function that ignores bad gg lines"""

    def ignore(line):
        """Return false if line is bad"""
        return not line or ["", ""] == line or [start_line, ""] == line

    return ignore


def DefaultDelimitedSplitter(delimiter):
    """Wraps delimited splitter to handle empty records"""
    parser = DelimitedSplitter(delimiter=delimiter)

    def f(line):
        parsed = parser(line)
        if len(parsed) == 1:
            parsed.append("")
        return parsed

    return f


def MinimalGreengenesParser(lines, LineDelim="=", RecStart="BEGIN", RecEnd="END"):
    """Parses raw Greengeens 16S rRNA Gene records

    lines  :  open records file
    LineDelim  :  individual line delimiter, eg foo=bar
    RecStart  :  start identifier for a record
    RecEnd  :  end identifier for a record
    """
    line_parser = DefaultDelimitedSplitter(delimiter=LineDelim)

    # parse what the ending record looks like so it can match after being split
    RecordDelim = line_parser(RecEnd)

    # make sure to ignore the starting record
    ignore = make_ignore_f(RecStart)

    parser = DelimitedRecordFinder(
        RecordDelim, constructor=line_parser, keep_delimiter=False, ignore=ignore
    )

    for record in parser(lines):
        yield GenericRecord(record)


all_ids = lambda x, y: True
specific_ids = lambda x, y: x in y


def SpecificGreengenesParser(lines, fields, ids=None, **kwargs):
    """Yield specific fields from successive Greengenes records

    If ids are specified, only the records for the set of ids passed in will
    be returned. Parser will silently ignore ids that are not present in the
    set of ids as well as silently ignore ids in the set that are not present
    in the records file.

    ids : must either test True or be an iterable with prokMSA_ids

    Returns tuples in 'fields' order
    """
    parser = MinimalGreengenesParser(lines, **kwargs)

    if ids:
        ids = set(ids)
        id_lookup = specific_ids
    else:
        id_lookup = all_ids

    for record in parser:
        if id_lookup(record["prokMSA_id"], ids):
            yield tuple([record[field] for field in fields])
