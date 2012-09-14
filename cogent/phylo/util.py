#!/usr/bin/env python

import numpy
Float = numpy.core.numerictypes.sctype2char(float)
# Distance matricies are presently represented as simple dictionaries, which
# need to be converted into numpy arrays before being fed into phylogenetic
# reconstruction algorithms.

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "pm67nz@gmail.com"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def namesFromDistanceDict(dists):
    """Unique names from within the tuples which make up the keys of 'dists'"""
    names = []
    for key in dists:
        for name in key:
            if name not in names:
                names.append(name)
    return names

def lookupSymmetricDict(dists, a, b):
    """dists[a,b] or dists[b,a], whichever is present, so long as they
    don't contradict each other"""
    v1 = dists.get((a, b), None)
    v2 = dists.get((b, a), None)
    if v1 is None and v2 is None:
        raise KeyError((a,b))
    elif v1 is None or v2 is None or v1 == v2:
        return v1 or v2
    else:
        raise ValueError("d[%s,%s] != d[%s,%s]" % (a,b,b,a))

def distanceDictTo2D(dists):
    """(names, dists).  Distances converted into a straightforward distance
    matrix"""
    names = namesFromDistanceDict(dists)
    L = len(names)
    d = numpy.zeros([L, L], Float)
    for (i, a) in enumerate(names):
        for (j, b) in enumerate(names):
            if i != j:
                d[i, j] = lookupSymmetricDict(dists, a, b)
    return (names, d)

def triangularOrder(keys):
    """Indices for extracting a 1D representation of a triangular matrix
    where j > i and i is the inner dimension:
    Yields (0,1), (0,2), (1, 2), (0,3), (1,3), (2,3), (0,4)..."""
    N = len(keys)
    for j in range(1, N):
        for i in range(0, j):
            yield (keys[i], keys[j])
            
def distanceDictAndNamesTo1D(dists, names):
    """Distances converted into a triangular matrix implemented as a 1D array
    where j > i and i is the inner dimension:
    d[0,1], d[0, 2], d[1, 2], d[0, 3]..."""
    d = []
    for (name_i, name_j) in triangularOrder(names):
        d.append(lookupSymmetricDict(dists, name_i, name_j))
    return numpy.array(d)

def distanceDictTo1D(dists):
    """(names, dists).  Distances converted into a triangular matrix
    implemented as a 1D array where j > i and i is the inner dimension:
    d[0,1], d[0, 2], d[1, 2], d[0, 3]..."""
    names = namesFromDistanceDict(dists)
    d = distanceDictAndNamesTo1D(dists, names)
    return (names, d)
