#!/usr/bin/env python

"""computes probabilities for the kolmogorov distribution.

Translated from R 2.4 by Gavin Huttley
"""

from numpy import arange, array, asarray
from numpy import dot as matrixmultiply
from numpy import (
    exp,
    fabs,
    floor,
    log,
    ones,
    pi,
    ravel,
    reshape,
    sqrt,
    sum,
    zeros,
)

from cogent3.maths.stats.special import combinations


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

PIO4 = pi / 4
PIO2 = pi / 2
INVSQRT2PI = 1 / sqrt(2 * pi)


def mpower(A, exponent):
    """matrix power"""
    new = A
    for i in range(1, exponent):
        new = matrixmultiply(new, A)
    return new


def pkolmogorov1x(statistic, n):
    """Probability function for the one-sided one sample Kolmogorov statistics.

    Translated from R 2.4."""
    statistic = asarray(statistic)
    if statistic <= 0:
        return 0.0
    if statistic >= 1:
        return 1.0
    to = floor(n * (1 - statistic)) + 1
    j = arange(0, to)
    coeffs = asarray([log(combinations(n, i)) for i in j])
    p = sum(
        exp(
            coeffs
            + (n - j) * log(1 - statistic - j / n)
            + (j - 1) * (log(statistic + j / n))
        )
    )
    return 1 - statistic * p


def pkolmogorov2x(statistic, n):
    """Probability function for Kolmogorovs distribution."""
    k = int(n * statistic) + 1
    m = 2 * k - 1
    h = k - n * statistic
    H = ones(m ** 2, "d")
    for i in range(m):
        for j in range(m):
            if i - j + 1 < 0:
                H[i * m + j] = 0

    for i in range(m):
        H[i * m] -= h ** (i + 1)
        H[(m - 1) * m + i] -= h ** (m - i)
    H[(m - 1) * m] += [0, (2 * h - 1) ** m][2 * h - 1 > 0]
    for i in range(m):
        for j in range(m):
            if i - j + 1 > 0:
                for g in range(1, i - j + 2):
                    H[i * m + j] /= g
    Q = ravel(mpower(reshape(H, (m, m)), n))
    s = Q[(k - 1) * m + k - 1]
    for i in range(1, n + 1):
        s *= i / n
    return s


def pkstwo(x_vector, tolerance=1e-6):
    """Probability from the Kolmogorov asymptotic distribution."""
    # if isinstance(x_vector, float):
    #    x_vector = asarray(x_vector)
    x_vector = array(x_vector, ndmin=1)
    size = len(x_vector)
    k_max = int(sqrt(2 - log(tolerance)))
    for i in range(size):
        if x_vector[i] < 1:
            z = -(PIO2 * PIO4) / x_vector[i] ** 2
            w = log(x_vector[i])
            s = 0
            for k in range(1, k_max, 2):
                s += exp(k ** 2 * z - w)
            x_vector[i] = s / INVSQRT2PI
        else:
            z = -2 * x_vector[i] ** 2
            s = -1
            k = 1
            old = 0
            new = 1
            while fabs(old - new) > tolerance:
                old = new
                new += 2 * s * exp(z * k ** 2)
                s *= -1
                k += 1
            x_vector[i] = new

    return x_vector


def psmirnov2x(statistic, least, most):
    if least > most:
        least, most = most, least
    q = floor(statistic * most * least - 1e-7) / (least * most)
    u_vector = zeros(most + 1, "d")
    for j in range(most + 1):
        # SUPPORT2425
        u_vector[j] = [1, 0][int(j / most > q)]
    for i in range(1, least + 1):
        w = i / (i + most)
        if i / least > q:
            u_vector[0] = 0
        else:
            u_vector[0] = w * u_vector[0]
        for j in range(1, most + 1):
            if fabs(i / least - j / most) > q:
                u_vector[j] = 0
            else:
                u_vector[j] = w * u_vector[j] + u_vector[j - 1]

    return u_vector[most]
