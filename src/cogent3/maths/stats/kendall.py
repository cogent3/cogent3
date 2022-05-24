"""
Computes kendalls tau statistic and associated probabilities. A wrapper
function is provided in cogent3.maths.stats.kendall_correlation

Translated from R 2.5 by Gavin Huttley
"""


from numpy import array, floor, sqrt

from cogent3.maths.stats.distribution import zprob
from cogent3.maths.stats.number import CategoryCounter


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Daniel McDonald"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def as_paired_ranks(x, y):
    """return as matrix of paired ranks"""
    n = len(x)
    paired = list(zip(x, y))
    x = list(x)
    y = list(y)
    x.sort()
    y.sort()
    rank_val_map_x = dict(list(zip(x, list(range(n)))))
    rank_val_map_y = dict(list(zip(y, list(range(n)))))
    ranked = []
    for i in range(n):
        ranked += [[rank_val_map_x[paired[i][0]], rank_val_map_y[paired[i][1]]]]
    return ranked


def ckendall(k, n, w):
    # translated from R 2.5
    combin = n * (n - 1) / 2
    if k < 0 or k > combin:
        return 0
    if w[n][k] < 0:
        if n == 1:
            w[n][k] = k == 0
        else:
            s = 0
            for i in range(n):
                result = ckendall(k - i, n - 1, w)
                s += result
            w[n][k] = s
    return w[n][k]


def pkendall(x, n, divisor, working):
    # translated from R 2.5
    q = floor(x + 1e-7)
    if q < 0:
        x = 0
    elif q > n * (n - 1) / 2:
        x = 1
    else:
        p = 0
        for k in range(int(q) + 1):
            result = ckendall(k, n, working)
            p += result
        x = p / divisor
    return x


def kendalls_tau(x, y, return_p=True):
    """returns kendall's tau

    Parameters
    ----------
    return_p
        returns the probability from the normal approximation when
        True, otherwise just returns tau

    """
    ranked = as_paired_ranks(x, y)
    n = len(ranked)
    con = 0
    discor = 0
    x_tied = 0
    y_tied = 0
    for i in range(n - 1):
        x_1 = ranked[i][0]
        y_1 = ranked[i][1]
        for j in range(i + 1, n):
            x_2 = ranked[j][0]
            y_2 = ranked[j][1]
            x_diff = x_1 - x_2
            y_diff = y_1 - y_2
            if x_diff * y_diff > 0:
                con += 1
            elif x_diff and y_diff:
                discor += 1
            else:
                if x_diff:
                    y_tied += 1
                if y_diff:
                    x_tied += 1

    diff = con - discor
    total = con + discor
    denom = ((total + y_tied) * (total + x_tied)) ** 0.5
    variance = (4 * n + 10) / (9 * n * (n - 1))
    tau = diff / denom
    stat = tau

    if x_tied or y_tied:
        x_tied = array([v for v in CategoryCounter(x).values() if v > 1])
        y_tied = array([v for v in CategoryCounter(y).values() if v > 1])
        t0 = n * (n - 1) / 2
        t1 = sum(x_tied * (x_tied - 1)) / 2
        t2 = sum(y_tied * (y_tied - 1)) / 2
        stat = tau * sqrt((t0 - t1) * (t0 - t2))
        v0 = n * (n - 1) * (2 * n + 5)
        vt = sum(x_tied * (x_tied - 1) * (2 * x_tied + 5))
        vu = sum(y_tied * (y_tied - 1) * (2 * y_tied + 5))
        v1 = sum(x_tied * (x_tied - 1)) * sum(y_tied * (y_tied - 1))
        v2 = sum(x_tied * (x_tied - 1) * (x_tied - 2)) * sum(
            y_tied * (y_tied - 1) * (y_tied - 2)
        )
        variance = (
            (v0 - vt - vu) / 18
            + v1 / (2 * n * (n - 1))
            + v2 / (9 * n * (n - 1) * (n - 2))
        )
    if return_p:
        return tau, zprob(stat / variance ** 0.5)
    else:
        return tau
