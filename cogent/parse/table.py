#!/usr/bin/env python

import cPickle, csv
from record_finder import is_empty
from gzip import GzipFile

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class ConvertFields(object):
    """converter for input data to Table"""
    
    def __init__(self, conversion, by_column=True):
        """handles conversions of columns or lines
        
        Arguments:
            - by_column: conversion will by done for each column, otherwise
              done by entire line
            - """
        super(ConvertFields, self).__init__()
        self.conversion = conversion
        self.by_column = by_column
        
        self._func = self.convertByColumns
        
        if not self.by_column:
            assert callable(conversion), \
                "conversion must be callable to convert by line"
            self._func = self.convertByLine
    
    def convertByColumns(self, line):
        """converts each column in a line"""
        for index, cast in self.conversion:
            line[index] = cast(line[index])
        return line
    
    def convertByLine(self, line):
        """converts each column in a line"""
        return self.conversion(line)
    
    def _call(self, *args, **kwargs):
        return self._func(*args, **kwargs)
    
    __call__ = _call
    

def SeparatorFormatParser(with_header=True, converter = None, ignore = None,
                sep=",", strip_wspace=True, limit=None, **kw):
    """Returns a parser for a delimited tabular file.
    
    Arguments:
        - with_header: when True, first line is taken to be the header. Not
          passed to converter.
        - converter: a callable that returns a correctly formatted line.
        - ignore: lines for which ignore returns True are ignored. White-space
          lines are always skipped.
        - sep: the delimiter deparating fields.
        - strip_wspace: removes redundant white-space from strings.
        - limit: exits after this many lines"""
    sep = kw.get("delim", sep)
    if ignore is None: # keep all lines
        ignore = lambda x: False
    
    by_column = getattr(converter, 'by_column', True)
    
    def callable(lines):
        num_lines = 0
        header = None
        for line in lines:
            if is_empty(line):
                continue
            
            line = line.strip('\n').split(sep)
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

def autogen_reader(infile, sep, with_title, limit=None):
    """returns a SeparatorFormatParser with field convertor for numeric column
    types."""
    seen_title_line = False
    for first_data_row in infile:
        if seen_title_line:
            break
        if sep in first_data_row and not seen_title_line:
            seen_title_line = True
    
    infile.seek(0) # reset to start of file
    
    numeric_fields = []
    for index, value in enumerate(first_data_row.strip().split(sep)):
        try:
            v = float(value)
        except ValueError:
            try:
                v = long(value)
            except ValueError:
                continue
        
        numeric_fields += [(index, eval(value).__class__)]
    
    return SeparatorFormatParser(converter=ConvertFields(numeric_fields),
                                 sep=sep, limit=limit)

def load_delimited(filename, header = True, delimiter = ',',
        with_title = False, with_legend = False, limit=None):
    if limit is not None:
        limit += 1 # don't count header line
    
    if filename.endswith('gz'):
        f = GzipFile(filename, 'rb')
    else:
        f = file(filename, "U")
    
    reader = csv.reader(f, dialect = 'excel', delimiter = delimiter)
    rows = []
    num_lines = 0
    for row in reader:
        rows.append(row)
        num_lines += 1
        if limit is not None and num_lines >= limit:
            break
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

