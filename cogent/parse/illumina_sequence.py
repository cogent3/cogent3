#!/usr/bin/env python
"""Tests of Illumina sequence file parser.
"""

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Production"

def MinimalIlluminaSequenceParser(data):
    """ Yields (header, sequence, quality strings) given an Illumina sequence file
    """
    if hasattr(data,'lower'):
        # passed a string so open filepath
        data = open(data,'U')
    
    for line in data:
        fields = line.strip().split(':')
        yield fields[:-2], fields[-2], fields[-1]
    
    try:
        # if data is a file, close it
        data.close()
    except AttributeError:
        # otherwise do nothing
        pass
    
    return
