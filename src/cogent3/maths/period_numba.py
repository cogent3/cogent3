import numba
import numpy
from numba import njit

# turn off code coverage as njit-ted code not accessible to coverage


@njit(numba.float64(numba.float64[:], numba.int64, numba.int64), cache=True)
def goertzel_inner(x, N, period):  # pragma: no cover
    """returns the power from series x for period"""
    coeff = 2.0 * numpy.cos(2 * numpy.pi / period)
    s_prev = 0.0
    s_prev2 = 0.0
    for n in range(N):
        s = x[n] + coeff * s_prev - s_prev2
        s_prev2 = s_prev
        s_prev = s

    return numpy.sqrt(s_prev2**2 + s_prev**2 - coeff * s_prev2 * s_prev)


@njit(
    numba.complex128[:](
        numba.float64[:],
        numba.complex128[:],
        numba.complex128[:],
        numba.int64,
        numba.int64,
    ),
    cache=True,
)
def ipdft_inner(x, X, W, ulim, N):  # pragma: no cover
    for p in range(ulim):
        w = 1.0
        for n in range(N):
            if n != 0:
                w = w * W[p]
            X[p] = X[p] + x[n] * w
    return X


@njit(numba.void(numba.float64[:], numba.float64[:], numba.int64), cache=True)
def autocorr_inner(x, xc, N) -> None:  # pragma: no cover
    for m in range(-N + 1, N):
        for n in range(N):
            if 0 <= n - m < N:
                xc[m + N - 1] += x[n] * x[n - m]
