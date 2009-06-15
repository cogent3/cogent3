#!/usr/bin/env python

import cPickle, csv
from record_finder import is_empty

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def ConvertFields(conversions):
    """Factory function for converting indexed fields. Useful for the
    SeparatorFormatParser.
    
    Arguments:
        - conversions: a series consisting of index,converter callable pairs,
          eg [(0, int), (4, float)]"""
    def callable(line):
        for index, cast in conversions:
            line[index] = cast(line[index])
        return line
    
    return callable

def SeparatorFormatParser(with_header=True, converter = None, ignore = is_empty,
                sep=",", strip_wspace=True, **kw):
    """Returns a parser for a delimited tabular file.
    
    Arguments:
        - with_header: when True, first line is taken to be the header. Not
          passed to converter.
        - converter: a callable that returns a correctly formatted line.
        - ignore: lines for which ignore returns True are ignored
        - sep: the delimiter deparating fields.
        - strip_wspace: removes redundant white-space from strings."""
    sep = kw.get("delim", sep)
    def callable(lines):
        header = None
        for line in lines:
            if ignore(line):
                continue
            line = line.strip('\n').split(sep)
            if strip_wspace:
                line = [field.strip() for field in line]
            if with_header and not header:
                header = True
            elif converter:
                line = converter(line)
            yield line
    
    return callable

def load_delimited(filename, header = True, delimiter = ',',
        with_title = False, with_legend = False):
    f = file(filename, "U")
    reader = csv.reader(f, dialect = 'excel', delimiter = delimiter)
    rows = [row for row in reader]
    f.close()
    if with_title:
        title = ''.join(rows.pop(0))
    else:
        title = ''
    if header:
        header = rows.pop(0)
    else:
        header = None
    if with_legend:
        legend = ''.join(rows.pop(-1))
    else:
        legend = ''
    # now do type casting in the order int, float, default is string
    for row in rows:
        for cdex, cell in enumerate(row):
            try:
                cell = int(cell)
                row[cdex] = cell
            except ValueError:
                try:
                    cell = float(cell)
                    row[cdex] = cell
                except ValueError:
                    pass
                pass
    return header, rows, title, legend

