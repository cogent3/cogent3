#!/usr/bin/env python
"""Retrieve information from web databases.
"""
from urllib import urlopen, urlretrieve, quote_plus

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class UrlGetter(object):
    Defaults = {}       #override in derived classes -- default values
    PrintedFields = {}  #override in derived classes -- fields to print
    BaseUrl = ''        #override in derived classes
    KeyValDelimiter = '='
    FieldDelimiter = '&'
    
    def __init__(self, **kwargs):
        """Returns new instance with arbitrary kwargs."""
        self.__dict__.update(self.Defaults)
        self.__dict__.update(kwargs)
        self._temp_args = {}

    def __str__(self):
        to_get = self.__dict__.copy()
        to_get.update(self._temp_args)
        return self.BaseUrl + self.FieldDelimiter.join(\
            [quote_plus(k)+self.KeyValDelimiter+quote_plus(str(v)) for k, v in to_get.items()\
            if k in self.PrintedFields])

    def open(self, **kwargs):
        """Returns a stream handle to URL result, temporarily overriding kwargs.."""
        self._temp_args = kwargs
        result = urlopen(str(self))
        self._temp_args = {}
        return result

    def read(self, **kwargs):
        """Gets URL and reads into memory, temporarily overriding kwargs."""
        result = self.open(**kwargs)
        data = result.read()
        result.close()
        return data

    def retrieve(self, fname, **kwargs):
        """Gets URL and writes to file fname, temporarily overriding kwargs. 
        
        Note: produces no return value."""
        self._temp_args = kwargs
        urlretrieve(str(self), fname)
        self._temp_args = None
        
def expand_slice(s):
    """Takes a start and end accession, and gets the whole range.

    WARNING: Unlike standard slices, includes the last item in the range.
    In other words, obj[AF1001:AF1010] will include AF1010.
    
    Both accessions must have the same non-numeric prefix.
    """
    start, step, end = s.start, s.step, s.stop
    #find where the number is
    start_index = last_nondigit_index(start)
    end_index = last_nondigit_index(end)
    prefix = start[:start_index]
    if prefix != end[:end_index]:
        raise TypeError, "Range start and end don't have same prefix"
    if not step:
        step = 1
    range_start = long(start[start_index:])
    range_end = long(end[end_index:])
    field_width = str(len(start) - start_index)
    format_string = '%'+field_width+'.'+field_width+'d'
    return [prefix + format_string % i \
        for i in range(range_start, range_end+1, step)]
    
def last_nondigit_index(s):
    """Returns the index of s such that s[i:] is numeric, or None."""
    for i in range(len(s)):
        if s[i:].isdigit():
            return i
    #if we get here, there weren't any trailing digits
    return None
