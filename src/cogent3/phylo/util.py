#!/usr/bin/env python

import numpy


Float = numpy.core.numerictypes.sctype2char(float)
# Distance matricies are presently represented as simple dictionaries, which
# need to be converted into numpy arrays before being fed into phylogenetic
# reconstruction algorithms.

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "pm67nz@gmail.com"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def names_from_distance_dict(dists):
    """Unique names from within the tuples which make up the keys of 'dists'"""
    names = []
    for key in dists:
        for name in key:
            if name not in names:
                names.append(name)
    return names


def lookup_symmetric_dict(dists, a, b):
    """dists[a,b] or dists[b,a], whichever is present, so long as they
    don't contradict each other"""
    v1 = dists.get((a, b), None)
    v2 = dists.get((b, a), None)
    try:
        v1 = None if numpy.isnan(v1) else v1
    except TypeError:
        pass
    try:
        v2 = None if numpy.isnan(v2) else v2
    except TypeError:
        pass
    if v1 is None and v2 is None:
        raise KeyError((a, b))
    elif v1 is None or v2 is None or v1 == v2:
        return v1 or v2
    else:
        raise ValueError(f"d[{a},{b}] != d[{b},{a}]")


def distance_dict_to_2D(dists):
    """(names, dists).  Distances converted into a straightforward distance
    matrix"""
    names = names_from_distance_dict(dists)
    L = len(names)
    d = numpy.zeros([L, L], Float)
    for (i, a) in enumerate(names):
        for (j, b) in enumerate(names):
            if i != j:
                d[i, j] = lookup_symmetric_dict(dists, a, b)
    return (names, d)


def triangular_order(keys):
    """Indices for extracting a 1D representation of a triangular matrix
    where j > i and i is the inner dimension:
    Yields (0,1), (0,2), (1, 2), (0,3), (1,3), (2,3), (0,4)..."""
    N = len(keys)
    for j in range(1, N):
        for i in range(0, j):
            yield (keys[i], keys[j])


def distance_dict_and_names_to_1D(dists, names):
    """Distances converted into a triangular matrix implemented as a 1D array
    where j > i and i is the inner dimension:
    d[0,1], d[0, 2], d[1, 2], d[0, 3]..."""
    d = []
    for (name_i, name_j) in triangular_order(names):
        d.append(lookup_symmetric_dict(dists, name_i, name_j))
    return numpy.array(d)


def distance_dict_to_1D(dists):
    """(names, dists).  Distances converted into a triangular matrix
    implemented as a 1D array where j > i and i is the inner dimension:
    d[0,1], d[0, 2], d[1, 2], d[0, 3]..."""
    names = names_from_distance_dict(dists)
    d = distance_dict_and_names_to_1D(dists, names)
    return (names, d)
