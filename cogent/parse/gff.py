#!/usr/bin/env python

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Matthew Wakefield", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def GffParser(f):
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
        cols = line.split('\t')
        if len(cols) == 8:
            cols.append('')
        assert len(cols) == 9, line
        (seqname, source, feature, start, end, score,
                strand, frame, attributes) = cols
        
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
        if strand == '-':
            (start, end) = (end, start)
        
        # should parse attributes too
        yield (seqname, source, feature, start, end, score,
                strand, frame, attributes, comments)

def parse_attributes(attribute_string):
    """Returns region of attribute string between first pair of double quotes"""
    attribute_string = attribute_string[attribute_string.find('"')+1:]
    if '"' in attribute_string:
        attribute_string = attribute_string[:attribute_string.find('"')]
    return attribute_string

