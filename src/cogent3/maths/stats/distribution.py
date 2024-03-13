#!/usr/bin/env python
"""Translations of functions from Release 2.3 of the Cephes Math Library,
which is (c) Stephen L. Moshier 1984, 1995.
"""

from numpy import arctan as atan
from numpy import array, exp, sqrt
from scipy.stats import f, norm, t
from scipy.stats.distributions import chi2

from cogent3.maths.stats.special import (
    MACHEP,
    MAXNUM,
    PI,
    betai,
    expm1,
    fix_rounding_error,
    igam,
    igamc,
    igami,
    incbi,
    ln_binomial,
    log1p,
    ndtri,
)

from cogent3.util import warning as c3warn

# ndtri import b/c it should be available via this module


incbet = betai  # shouldn't have renamed it...

# Probability integrals: low gives left-hand tail, high gives right-hand tail.


def zprob(x):
    """Returns both tails of z distribution (-inf to -x, inf to x)."""
    return 2 * norm.sf(abs(x))


def tprob(x, df):
    """Returns both tails of t distribution (-infinity to -x, infinity to x)"""
    return 2 * t.sf(abs(x), df)


def poisson_high(successes, mean):
    """Returns right tail of Poission distribution, Pr(X > x).

    successes ranges from 0 to infinity. mean must be positive.
    """
    return pdtrc(successes, mean)


def poisson_low(successes, mean):
    """Returns left tail of Poisson distribution, Pr(X <= x).

    successes ranges from 0 to infinity. mean must be positive.
    """
    return pdtr(successes, mean)


@c3warn.deprecated_callable(version="2024.9", reason="use scipy.stats.poisson.pmf()instead", is_discontinued=True)
def poisson_exact(successes, mean): # pragma: no cover
    """being removed""" 
    if successes == 0:
        return pdtr(0, mean)
    elif successes < mean:  # use left tail
        return pdtr(successes, mean) - pdtr(successes - 1, mean)
    else:  # successes > mean: use right tail
        return pdtrc(successes - 1, mean) - pdtrc(successes, mean)


def binomial_exact(successes, trials, prob):
    """Returns binomial probability of exactly X successes.

    Works for integer and floating point values.

    Note: this function is only a probability mass function for integer
    values of 'trials' and 'successes', i.e. if you sum up non-integer
    values you probably won't get a sum of 1.
    """
    if (prob < 0) or (prob > 1):
        raise ValueError("Binomial prob must be between 0 and 1.")
    if (successes < 0) or (trials < successes):
        raise ValueError("Binomial successes must be between 0 and trials.")
    return exp(ln_binomial(successes, trials, prob))


def fprob(dfn, dfd, F, side="right"):
    """Returns both tails of F distribution (-inf to F and F to inf)

    Use in case of two-tailed test. Usually this method is called by
    f_two_sample, so you don't have to worry about choosing the right side.

    side: right means return twice the right-hand tail of the F-distribution.
        Use in case var(a) > var (b)
          left means return twice the left-hand tail of the F-distribution.
        Use in case var(a) < var(b)
    """
    if F < 0:
        raise ValueError(f"fprob: F must be >= 0 (got {F}).")
    if side == "right":
        return 2 * f.sf(F, dfn, dfd)
    elif side == "left":
        return 2 * f.cdf(F, dfn, dfd)
    else:
        raise ValueError(f"Not a valid value for side {side}")


def stdtr(k, t):
    """Student's t distribution, -infinity to t.

    See Cephes docs for details.
    """
    if k <= 0:
        raise ValueError("stdtr: df must be > 0.")
    if t == 0:
        return 0.5
    if t < -2:
        rk = k
        z = rk / (rk + t * t)
        return 0.5 * betai(0.5 * rk, 0.5, z)
    # compute integral from -t to + t
    if t < 0:
        x = -t
    else:
        x = t

    rk = k  # degrees of freedom
    z = 1 + (x * x) / rk
    # test if k is odd or even
    if (k & 1) != 0:
        # odd k
        xsqk = x / sqrt(rk)
        p = atan(xsqk)
        if k > 1:
            f = 1
            tz = 1
            j = 3
            while (j <= (k - 2)) and ((tz / f) > MACHEP):
                tz *= (j - 1) / (z * j)
                f += tz
                j += 2
            p += f * xsqk / z
        p *= 2 / PI
    else:
        # even k
        f = 1
        tz = 1
        j = 2
        while (j <= (k - 2)) and ((tz / f) > MACHEP):
            tz *= (j - 1) / (z * j)
            f += tz
            j += 2
        p = f * x / sqrt(z * rk)
    # common exit
    if t < 0:
        p = -p  # note destruction of relative accuracy
    p = 0.5 + 0.5 * p
    return p


def bdtr(k, n, p):
    """Binomial distribution, 0 through k.

    Uses formula bdtr(k, n, p) = betai(n-k, k+1, 1-p)

    See Cephes docs for details.
    """
    p = fix_rounding_error(p)
    if (p < 0) or (p > 1):
        raise ValueError("Binomial p must be between 0 and 1.")
    if (k < 0) or (n < k):
        raise ValueError("Binomial k must be between 0 and n.")
    if k == n:
        return 1
    dn = n - k
    if k == 0:
        return pow(1 - p, dn)
    else:
        return betai(dn, k + 1, 1 - p)


def bdtrc(k, n, p):
    """Complement of binomial distribution, k+1 through n.

    Uses formula bdtrc(k, n, p) = betai(k+1, n-k, p)

    See Cephes docs for details.
    """
    p = fix_rounding_error(p)
    if (p < 0) or (p > 1):
        raise ValueError("Binomial p must be between 0 and 1.")
    if (k < 0) or (n < k):
        raise ValueError("Binomial k must be between 0 and n.")
    if k == n:
        return 0
    dn = n - k
    if k == 0:
        if p < 0.01:
            dk = -expm1(dn * log1p(-p))
        else:
            dk = 1 - pow(1.0 - p, dn)
    else:
        dk = k + 1
        dk = betai(dk, dn, p)
    return dk


def pdtr(k, m):
    """Returns sum of left tail of Poisson distribution, 0 through k.

    See Cephes docs for details.
    """
    if k < 0:
        raise ValueError("Poisson k must be >= 0.")
    if m < 0:
        raise ValueError("Poisson m must be >= 0.")
    return igamc(k + 1, m)


def pdtrc(k, m):
    """Returns sum of right tail of Poisson distribution, k+1 through infinity.

    See Cephes docs for details.
    """
    if k < 0:
        raise ValueError("Poisson k must be >= 0.")
    if m < 0:
        raise ValueError("Poisson m must be >= 0.")
    return igam(k + 1, m)


def gdtr(a, b, x):
    """Returns integral from 0 to x of Gamma distribution with params a and b."""
    if x < 0.0:
        raise ZeroDivisionError("x must be at least 0.")
    return igam(b, a * x)


def gdtrc(a, b, x):
    """Returns integral from x to inf of Gamma distribution with params a and b."""
    if x < 0.0:
        raise ZeroDivisionError("x must be at least 0.")
    return igamc(b, a * x)


# note: ndtri for the normal distribution is already imported


def stdtri(k, p):
    """Returns inverse of Student's t distribution. k = df."""
    p = fix_rounding_error(p)
    # handle easy cases
    if k <= 0 or p < 0.0 or p > 1.0:
        raise ZeroDivisionError("k must be >= 1, p between 1 and 0.")
    rk = k
    # handle intermediate values
    if p > 0.25 and p < 0.75:
        if p == 0.5:
            return 0.0
        z = 1.0 - 2.0 * p
        z = incbi(0.5, 0.5 * rk, abs(z))
        t = sqrt(rk * z / (1.0 - z))
        if p < 0.5:
            t = -t
        return t
    # handle extreme values
    rflg = -1
    if p >= 0.5:
        p = 1.0 - p
        rflg = 1
    z = incbi(0.5 * rk, 0.5, 2.0 * p)

    if MAXNUM * z < rk:
        return rflg * MAXNUM
    t = sqrt(rk / z - rk)
    return rflg * t


def pdtri(k, p):
    """Inverse of Poisson distribution.

    Finds Poission mean such that integral from 0 to k is p.
    """
    p = fix_rounding_error(p)
    if k < 0 or p < 0.0 or p >= 1.0:
        raise ZeroDivisionError("k must be >=0, p between 1 and 0.")
    v = k + 1
    return igami(v, p)


def bdtri(k, n, y):
    """Inverse of binomial distribution.

    Finds binomial p such that sum of terms 0-k reaches cum probability y.
    """
    y = fix_rounding_error(y)
    if y < 0.0 or y > 1.0:
        raise ZeroDivisionError("y must be between 1 and 0.")
    if k < 0 or n <= k:
        raise ZeroDivisionError("k must be between 0 and n")
    dn = n - k
    if k == 0:
        if y > 0.8:
            p = -expm1(log1p(y - 1.0) / dn)
        else:
            p = 1.0 - y ** (1.0 / dn)
    else:
        dk = k + 1
        p = incbet(dn, dk, 0.5)
        if p > 0.5:
            p = incbi(dk, dn, 1.0 - y)
        else:
            p = 1.0 - incbi(dn, dk, y)
    return p


def gdtri(a, b, y):
    """Returns Gamma such that y is the probability in the integral.

    WARNING: if 1-y == 1, gives incorrect result. The scipy implementation
    gets around this by using cdflib, which is in Fortran. Until someone
    gets around to translating that, only use this function for values of
    p greater than 1e-15 or so!
    """
    y = fix_rounding_error(y)
    if y < 0.0 or y > 1.0 or a <= 0.0 or b < 0.0:
        raise ZeroDivisionError("a and b must be non-negative, y from 0 to 1.")
    return igami(b, 1.0 - y) / a


def fdtri(a, b, y):
    """Returns inverse of F distribution."""
    y = fix_rounding_error(y)
    if a < 1.0 or b < 1.0 or y <= 0.0 or y > 1.0:
        raise ZeroDivisionError("y must be between 0 and 1; a and b >= 1")
    y = 1.0 - y
    # Compute probability for x = 0.5
    w = incbet(0.5 * b, 0.5 * a, 0.5)
    # If that is greater than y, then the solution w < .5.
    # Otherwise, solve at 1-y to remove cancellation in (b - b*w).
    if w > y or y < 0.001:
        w = incbi(0.5 * b, 0.5 * a, y)
        x = (b - b * w) / (a * w)
    else:
        w = incbi(0.5 * a, 0.5 * b, 1.0 - y)
        x = b * w / (a * (1.0 - w))
    return x


def probability_points(n):
    """return series of n probabilities

    Returns
    -------
    Numpy array of probabilities

    Notes
    -----
    Useful for plotting probability distributions
    """
    assert n > 0, f"{n} must be > 0"
    adj = 0.5 if n > 10 else 3 / 8
    denom = n if n > 10 else n + 1 - 2 * adj
    return array([(i - adj) / denom for i in range(1, n + 1)])


def theoretical_quantiles(n, dist, *args):
    """returns theoretical quantiles from dist

    Parameters
    ----------
    n : int
        number of elements
    dist : str
        one of 'normal', 'chisq', 't', 'uniform'

    Returns
    -------
    Numpy array of quantiles
    """
    dist = dist.lower()
    funcs = dict(
        normal=ndtri,
        chisq=chi2.isf,
        t=stdtri,
    )

    if dist != "uniform" and dist not in funcs:
        raise ValueError(f"'{dist} not in {list(funcs)}")

    probs = probability_points(n)
    if dist == "uniform":
        return probs

    func = funcs[dist]

    if not args:
        return array([func(p) for p in probs])

    if (
        dist == "chisq"
    ):  # to use scipy.stats.distributions.chi2.isf we need to reverse the order of the arguments
        return array([func(*((p,) + args)) for p in probs])

    return array([func(*(args + (p,))) for p in probs])
