#!/usr/bin/env python
"""Provides small utility functions for numpy arrays.
"""
from operator import __getitem__ as getitem
from operator import mul

import numpy

from numpy import (
    arange,
    argmin,
    argsort,
    array,
    clip,
    compress,
    concatenate,
    cumsum,
    identity,
    less,
    log,
    logical_not,
    maximum,
    min,
    newaxis,
    nonzero,
    pi,
    product,
    put,
    ravel,
    repeat,
    reshape,
    searchsorted,
    sort,
    sqrt,
    sum,
    take,
    trace,
    where,
    zeros,
)
from numpy.random import normal, randint


numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)
Int = numerictypes(int)
err = numpy.seterr(divide="raise")

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2020.12.14a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"


def safe_p_log_p(data):
    """Returns -(p*log2(p)) for every non-negative, nonzero p in a.

    WARNING: log2 is only defined on positive numbers, so make sure
    there are no negative numbers in the array.

    Always returns an array with floats in there to avoid unexpected
    results when applying it to an array with just integers.
    """

    non_zero = data != 0
    result = numpy.zeros(data.shape, dtype=float)
    with numpy.errstate(invalid="raise"):
        result[non_zero] = -data[non_zero] * numpy.log2(data[non_zero])
    return result


def safe_log(data):
    """Returns the log (base 2) of each nonzero item in a.

    WARNING: log2 is only defined on positive numbers, so make sure
    there are no negative numbers in the array. Will either return an
    array containing floating point exceptions or will raise an
    exception, depending on platform.

    Always returns an array with floats in there to avoid unexpected
    results when applying it to an array with just integers.
    """

    non_zero = data != 0
    result = numpy.zeros(data.shape, dtype=float)
    with numpy.errstate(invalid="raise"):
        result[non_zero] = numpy.log2(data[non_zero])
    return result


def row_uncertainty(a):
    """Returns uncertainty (Shannon's entropy) for each row in a IN BITS

    a: numpy array (has to be 2-dimensional!)

    The uncertainty is calculated in BITS not NATS!!!

    Will return 0 for every empty row, but an empty array for every empty column,
    thanks to this sum behavior:
    >>> sum(array([[]]),1)
    array([0])
    >>> sum(array([[]]))
    zeros((0,), 'l')
    """
    try:
        return sum(safe_p_log_p(a), 1)
    except ValueError:
        raise ValueError("Array has to be two-dimensional")


def column_uncertainty(a):
    """Returns uncertainty (Shannon's entropy) for each column in a in BITS

    a: numpy array (has to be 2-dimensional)

    The uncertainty is calculated in BITS not NATS!!!

    Will return 0 for every empty row, but an empty array for every empty column,
    thanks to this sum behavior:
    >>> sum(array([[]]),1)
    array([0])
    >>> sum(array([[]]))
    zeros((0,), 'l')

    """
    if len(a.shape) < 2:
        raise ValueError("Array has to be two-dimensional")
    return sum(safe_p_log_p(a), axis=0)


def row_degeneracy(a, cutoff=0.5):
    """Returns the number of characters that's needed to cover >= cutoff

    a: numpy array
    cutoff: number that should be covered in the array

    Example:
    [   [.1 .3  .4  .2],
        [.5 .3  0   .2],
        [.8 0   .1  .1]]
    if cutoff = .75: row_degeneracy -> [3,2,1]
    if cutoff = .95: row_degeneracy -> [4,3,3]

    WARNING: watch out with floating point numbers.
    if the cutoff= 0.9 and in the array is also 0.9, it might not be found
    >>> searchsorted(cumsum(array([.6,.3,.1])),.9)
    2
    >>> searchsorted(cumsum(array([.5,.4,.1])),.9)
    1

    If the cutoff value is not found, the result is clipped to the
    number of columns in the array.
    """
    if not a.any():
        return []
    try:
        b = cumsum(sort(a)[:, ::-1], 1)
    except IndexError:
        raise ValueError("Array has to be two dimensional")
    degen = [searchsorted(aln_pos, cutoff) for aln_pos in b]
    # degen contains now the indices at which the cutoff was hit
    # to change to the number of characters, add 1
    return clip(array(degen) + 1, 0, a.shape[1])


def column_degeneracy(a, cutoff=0.5):
    """Returns the number of characters that's needed to cover >= cutoff

    a: numpy array
    cutoff: number that should be covered in the array

    Example:
    [   [.1 .8  .3],
        [.3 .2  .3],
        [.6 0   .4]]
    if cutoff = .75: column_degeneracy -> [2,1,3]
    if cutoff = .45: column_degeneracy -> [1,1,2]

    WARNING: watch out with floating point numbers.
    if the cutoff= 0.9 and in the array is also 0.9, it might not be found
    >>> searchsorted(cumsum(array([.6,.3,.1])),.9)
    2
    >>> searchsorted(cumsum(array([.5,.4,.1])),.9)
    1

    If the cutoff value is not found, the result is clipped to the
    number of rows in the array.
    """
    if not a.any():
        return []
    b = cumsum(sort(a, 0)[::-1], axis=0)
    try:
        degen = [searchsorted(b[:, idx], cutoff) for idx in range(len(b[0]))]
    except TypeError:
        raise ValueError("Array has to be two dimensional")
    # degen contains now the indices at which the cutoff was hit
    # to change to the number of characters, add 1
    return clip(array(degen) + 1, 0, a.shape[0])


def validate_freqs_array(data, axis=None):
    """input data is a valid frequency array
    Parameters
    ----------
    data : ndarray
        numpy array
    axis : int or None
        indicates along which axis the array should sum to 1

    Notes
    -----
    Raises ValueError if any element < 0 or series do not sum to 1
    """
    if (data < 0).any():
        raise ValueError("negative frequency not allowed")

    # we explicitly ignore nan
    result = data.sum(axis=axis)
    if not numpy.allclose(result[numpy.isnan(result) == False], 1):
        raise ValueError("invalid frequencies, sum(axis=1) is not equal to 1")
