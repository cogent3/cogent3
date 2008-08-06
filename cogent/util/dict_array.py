#!/usr/bin/env python
"""Wrapper for numpy arrays so that they can be indexed by name
    
    >>> a = numpy.identity(3, int)
    >>> b = DictArrayTemplate('abc', 'ABC').wrap(a)
    >>> b[0]
    ===========
    A    B    C
    -----------
    1    0    0
    -----------
    >>> b['a']
    ===========
    A    B    C
    -----------
    1    0    0
    -----------
    >>> b.keys()
    ['a', 'b', 'c']
    >>> b['a'].keys()
    ['A', 'B', 'C']
"""

import numpy
from cogent.format import table

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class DictArrayTemplate(object):
    def __init__(self, *dimensions):
        self.names = []
        self.ordinals = []
        for names in dimensions:
            if isinstance(names, int):
                names = range(names)
            else:
                names = list(names)[:]
            self.names.append(names)
            self.ordinals.append(dict((c,i) for (i,c) in enumerate(names)))
    
    def wrap(self, array, dtype = None):
        # dtype is numpy
        array = numpy.asarray(array, dtype=dtype)
        for (dim, categories) in enumerate(self.names):
            assert len(categories) == numpy.shape(array)[dim], "cats=%s; dim=%s" % (categories, dim)
        return DictArray(array, self)
    
    def interpretIndex(self, names):
        if not isinstance(names, tuple):
            names = (names,)
        index = []
        remaining = []
        for (ordinals, allnames, name) in zip(self.ordinals, self.names, names):
            if not isinstance(name, (int, slice)):
                name = ordinals[name]
            elif isinstance(name, slice):
                start = name.start
                stop = name.stop
                if isinstance(name.start, str):
                    start = allnames.index(name.start)
                if isinstance(name.stop, str):
                    stop = allnames.index(name.stop)
                name = slice(start, stop, name.step)
                remaining.append(allnames.__getitem__(name))
            index.append(name)
        remaining.extend(self.names[len(index):])
        if remaining:
            klass = type(self)(*remaining)
        else:
            klass = None
        return (tuple(index), klass)
    
    def array_repr(self, a):
        if len(a.shape) == 1:
            heading = [str(n) for n in self.names[0]]
            a = a[numpy.newaxis, :]
        elif len(a.shape) == 2:
            heading = [''] + [str(n) for n in self.names[1]]
            a = [[str(name)] + list(row) for (name, row) in zip(self.names[0], a)]
        else:
            return '%s dimensional %s' % (
                len(self.names), type(self).__name__)
        
        formatted = table.formattedCells(rows=a, header=heading)
        return str(table.simpleFormat(formatted[0], formatted[1], space=4))
    

class DictArray(object):
    """Wraps a numpy array so that it can be indexed with strings like nested
    dictionaries (only ordered), for things like substitution matricies and
    bin probabilities."""
    
    def __init__(self, *args, **kwargs):
        """allow alternate ways of creating for time being"""
        if len(args) <= 2:
            self.array = args[0]
            self.template = args[1]
        else:
            if 'dtype' in kwargs or 'typecode' in kwargs:
                dtype = kwargs['dtype']
                kwargs.pop('dtype', None)
                kwargs.pop('typecode', None)
            else:
                dtype = None
            create_new = DictArrayTemplate(*args[1:]).wrap(args[0], dtype=dtype)
            self.__dict__ = create_new.__dict__
        self.Shape = self.array.shape
    
    def asarray(self):
        return self.array
    
    def __array__(self):
        return self.array
    
    def asdict(self):
        return dict(self.items())
    
    def __getitem__(self, names):
        (index, remaining) = self.template.interpretIndex(names)
        result = self.array[index]
        if remaining is not None:
            result = self.__class__(result, remaining)
        return result
    
    def __iter__(self):
        (index, remaining) = self.template.interpretIndex(0)
        for elt in self.array:
            if remaining is None:
                yield elt
            else:
                yield remaining.wrap(elt)
    
    def __len__(self):
        return len(self.template.names[0])
    
    def keys(self):
        return self.template.names[0][:]
    
    def items(self):
        return [(n,self[n]) for n in self.keys()]
    
    def __repr__(self):
        return self.template.array_repr(self.array)
    
    def __eq__(self, other):
        if not isinstance(other, dict):
            other = other.asdict()
        return self.asdict() == other
    
