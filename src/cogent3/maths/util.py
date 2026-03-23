"""Provides small utility functions for numpy arrays."""

from collections.abc import Sequence as PySeq

import numpy
import numpy.typing as npt
from numpy import sum

NumpyArrayType = npt.NDArray[numpy.number]
NumpyIntArrayType = npt.NDArray[numpy.integer]
NumpyFloatArrayType = npt.NDArray[numpy.floating]


err = numpy.seterr(divide="raise")


def safe_p_log_p(data: NumpyArrayType) -> NumpyFloatArrayType:
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


def safe_log(data: NumpyArrayType) -> NumpyFloatArrayType:
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


def validate_freqs_array(data: NumpyFloatArrayType, axis: int | None = None) -> None:
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
        msg = "negative frequency not allowed"
        raise ValueError(msg)

    # we explicitly ignore nan
    result = data.sum(axis=axis)
    if not numpy.allclose(result[numpy.isnan(result) == False], 1):  # noqa
        msg = "invalid frequencies, sum(axis=1) is not equal to 1"
        raise ValueError(msg)


def ratios_to_proportions(total: float, params: PySeq[float]) -> list[float]:
    """Produces a list of N proportions from N-1 ratios and a total

    A recursive function that is the inverse of  proportions_to_ratios.

    Paramters
    ---------
    total: int
         The sum of `values` the array put into  proportions_to_ratios

    params: sequence
          The sequence output by  proportions_to_ratios, with int values
          between 0 and infinity

    Returns
    -------
    sequence: The `values` , the array that was put into
     proportions_to_ratios to get `params`

    Examples
    --------
    >>> ratios_to_proportions(1.0, [3, 1, 1])
    [0.125, 0.125, 0.375, 0.375]
    """

    if len(params) == 0:
        return [total]

    assert params[0] > 0, f"Ratios must be positive: {params[0]}"
    half = (len(params) + 1) // 2
    part = 1.0 / (params[0] + 1.0)  # ratio -> proportion
    return ratios_to_proportions(total * part, params[1:half]) + ratios_to_proportions(
        total * (1.0 - part),
        params[half:],
    )


def proportions_to_ratios(values: PySeq[float]) -> list[float]:
    """Produces a list of N-1 ratios from N proportions

    An invertible map that takes `values` an array of N numbers > 0
    whose sum is total and converts to an array of N-1 numbers between 0 and
    infinity.

    Parameters
    ----------
    values: sequence
        A sequence of N ints, where the ints
        have values between 0 and 1 (non inclusive).

    Raises
    ------
    AssertionError Exception
         Raises if there will be a negative or 0 value in the output array.

    Returns
    -------
    list: returns a list of size N-1 ints, where the unts take values
    between 0 and infinity.

    Notes:
    ------
    The function recursively halves the list into left side and right
    side (in that order) and for each halving divides the sum of the right half by the
    sum of the left half and adds the resulting value to the returned list.

    Examples
    --------

    >>>  proportions_to_ratios([0.125, 0.125, 0.375, 0.375])
    [3, 1, 1]

    >>>  proportions_to_ratios([0.1, 0.2, 0.9, 0])
    AssertionError
    """
    if len(values) == 1:
        return []
    half = len(values) // 2
    (num, denom) = (sum(values[half:]), sum(values[:half]))
    assert (num > 0) and (denom > 0)  # noqa
    ratio = num / denom
    return [
        ratio,
        *proportions_to_ratios(values[:half]),
        *proportions_to_ratios(values[half:]),
    ]
