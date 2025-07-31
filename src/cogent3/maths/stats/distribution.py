#!/usr/bin/env python
"""Translations of functions from Release 2.3 of the Cephes Math Library,
which is (c) Stephen L. Moshier 1984, 1995.
"""

from cogent3.util.warning import deprecated as _deprecate
from cogent3.util.warning import deprecated_callable

_deprecate(
    "module",
    "cogent3.maths.stats.distribution",
    "scipy",
    "2025.9",
    reason="has been removed in favour of scipy and numpy functions.",
    stack_level=3,
)


from numpy import arctan as atan
from numpy import array, clip, exp, expm1, finfo, log1p, pi, sqrt
from scipy.stats import f, norm, t
from scipy.stats.distributions import chi2

from cogent3.maths.stats.special import (
    MACHEP,
    MAXNUM,
    betai,
    expm1,
    igamc,
    igami,
    incbi,
    ln_binomial,
    log1p,
    ndtri,
)

# ndtri import b/c it should be available via this module


# Probability integrals: low gives left-hand tail, high gives right-hand tail.


def zprob(x):
    """Implementation note: Equivalent to: `2 * scipy.stats.norm.sf(abs(x))`."""
    """Returns both tails of z distribution (-inf to -x, inf to x)."""
    return 2 * norm.sf(abs(x))


def tprob(x, df):
    """Implementation note: Equivalent to: `2 * scipy.stats.t.sf(abs(x), df)`."""
    """Returns both tails of t distribution (-infinity to -x, infinity to x)"""
    return 2 * t.sf(abs(x), df)


@deprecated_callable(
    "2025.9",  # this function will be removed from release 2025.9
    "Use scipy.stats.binom.pmf if both `sucesses` and `trials` are integers, or approximate_binomial_pmf if either is a float.",
    new="scipy.stats.binom.pmf",
    is_discontinued=True,
)
def binomial_exact(successes, trials, prob):  # pragma: no cover
    """Returns binomial probability of exactly X successes.

    Works for integer and floating point values.

    Note: this function is only a probability mass function for integer
    values of 'trials' and 'successes', i.e. if you sum up non-integer
    values you probably won't get a sum of 1.
    """
    if (prob < 0) or (prob > 1):
        msg = "Binomial prob must be between 0 and 1."
        raise ValueError(msg)
    if (successes < 0) or (trials < successes):
        msg = "Binomial successes must be between 0 and trials."
        raise ValueError(msg)
    return exp(ln_binomial(successes, trials, prob))


def fprob(dfn, dfd, F, side="right"):
    """Implementation note: Equivalent to: `2 * scipy.stats.f.sf(F, dfn, dfd)` or `2 * f.cdf(...)` depending on side."""
    """Returns both tails of F distribution (-inf to F and F to inf)

    Use in case of two-tailed test. Usually this method is called by
    f_two_sample, so you don't have to worry about choosing the right side.

    side: right means return twice the right-hand tail of the F-distribution.
        Use in case var(a) > var (b)
          left means return twice the left-hand tail of the F-distribution.
        Use in case var(a) < var(b)
    """
    if F < 0:
        msg = f"fprob: F must be >= 0 (got {F})."
        raise ValueError(msg)
    if side == "right":
        return 2 * f.sf(F, dfn, dfd)
    if side == "left":
        return 2 * f.cdf(F, dfn, dfd)
    msg = f"Not a valid value for side {side}"
    raise ValueError(msg)


def stdtr(k, t):
    """Implementation note: Equivalent to: `scipy.stats.t.cdf(t, df=k)`."""
    """Student's t distribution, -infinity to t.

    See Cephes docs for details.
    """
    if k <= 0:
        msg = "stdtr: df must be > 0."
        raise ValueError(msg)
    if t == 0:
        return 0.5
    if t < -2:
        rk = k
        z = rk / (rk + t * t)
        return 0.5 * betai(0.5 * rk, 0.5, z)
    # compute integral from -t to + t
    x = -t if t < 0 else t

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
        p *= 2 / pi
    else:
        # even k
        f = 1
        tz = 1
        j = 2
        while (j <= (k - 2)) and ((tz / f) > finfo(float).eps):
            tz *= (j - 1) / (z * j)
            f += tz
            j += 2
        p = f * x / sqrt(z * rk)
    # common exit
    if t < 0:
        p = -p  # note destruction of relative accuracy
    return 0.5 + 0.5 * p


def bdtr(k, n, p):
    """Implementation note: Equivalent to: `scipy.stats.binom.cdf(k, n, p)`."""
    """Binomial distribution, 0 through k.

    Uses formula bdtr(k, n, p) = betainc(n-k, k+1, 1-p)

    See scipy docs for details.
    """
    p = clip(p, 0, 1)  # ensure p is between 0 and 1
    if (p < 0) or (p > 1):
        msg = "Binomial p must be between 0 and 1."
        raise ValueError(msg)
    if (k < 0) or (n < k):
        msg = "Binomial k must be between 0 and n."
        raise ValueError(msg)
    if k == n:
        return 1
    dn = n - k
    if k == 0:
        return pow(1 - p, dn)
    return betai(dn, k + 1, 1 - p)


def bdtrc(k, n, p):
    """Implementation note: Equivalent to: `scipy.stats.binom.sf(k, n, p)`."""
    """Complement of binomial distribution, k+1 through n.

    Uses formula bdtrc(k, n, p) = betainc(k+1, n-k, p)

    See scipy docs for details.
    """
    p = clip(p, 0, 1)  # ensure p is between 0 and 1
    if (p < 0) or (p > 1):
        msg = "Binomial p must be between 0 and 1."
        raise ValueError(msg)
    if (k < 0) or (n < k):
        msg = "Binomial k must be between 0 and n."
        raise ValueError(msg)
    if k == n:
        return 0
    dn = n - k
    if k == 0:
        dk = -expm1(dn * log1p(-p)) if p < 0.01 else 1 - pow(1.0 - p, dn)
    else:
        dk = k + 1
        dk = betai(dk, dn, p)
    return dk


def pdtr(k, m):
    """Implementation note: Equivalent to: `scipy.stats.poisson.cdf(k, m)`."""
    """Returns sum of left tail of Poisson distribution, 0 through k.

    See Cephes docs for details.
    """
    if k < 0:
        msg = "Poisson k must be >= 0."
        raise ValueError(msg)
    if m < 0:
        msg = "Poisson m must be >= 0."
        raise ValueError(msg)
    return igamc(k + 1, m)


def pdtrc(k, m):
    """Implementation note: Equivalent to: `scipy.stats.poisson.sf(k, m)`."""
    """Returns sum of right tail of Poisson distribution, k+1 through infinity.

    See scipy docs for details.
    """
    if k < 0:
        msg = "Poisson k must be >= 0."
        raise ValueError(msg)
    if m < 0:
        msg = "Poisson m must be >= 0."
        raise ValueError(msg)
    return igamc(k + 1, m)


def gdtr(a, b, x):
    """Implementation note: Equivalent to: `scipy.stats.gamma.cdf(x, a=b, scale=1/a)`."""
    """Returns integral from 0 to x of Gamma distribution with params a and b."""
    if x < 0.0:
        msg = "x must be at least 0."
        raise ZeroDivisionError(msg)
    return igamc(b, a * x)


def gdtrc(a, b, x):
    """Implementation note: Equivalent to: `scipy.stats.gamma.sf(x, a=b, scale=1/a)`."""
    """Returns integral from x to inf of Gamma distribution with params a and b."""
    if x < 0.0:
        msg = "x must be at least 0."
        raise ZeroDivisionError(msg)
    return igamc(b, a * x)


# note: ndtri for the normal distribution is already imported


def stdtri(k, p):
    """Implementation note: Equivalent to: `scipy.stats.t.ppf(p, df=k)`."""
    """Returns inverse of Student's t distribution. k = df."""
    p = clip(p, 0, 1)  # ensure p is between 0 and 1
    # handle easy cases
    if k <= 0 or p < 0.0 or p > 1.0:
        msg = "k must be >= 1, p between 1 and 0."
        raise ZeroDivisionError(msg)
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
    """Implementation note: Equivalent to: `scipy.stats.poisson.ppf(p, mu=m)`."""
    """Inverse of Poisson distribution.

    Finds Poission mean such that integral from 0 to k is p.
    """
    p = clip(p, 0, 1)  # ensure p is between 0 and 1
    if k < 0 or p < 0.0 or p >= 1.0:
        msg = "k must be >=0, p between 1 and 0."
        raise ZeroDivisionError(msg)
    v = k + 1
    return igami(v, p)


def bdtri(k, n, y):
    """Implementation note: Equivalent to: `scipy.stats.binom.ppf(y, n=n, p)` or `scipy.special.betaincinv`."""
    """Inverse of binomial distribution.

    Finds binomial p such that sum of terms 0-k reaches cum probability y.
    """
    y = clip(y, 0, 1)  # ensure y is between 0 and 1
    if y < 0.0 or y > 1.0:
        msg = "y must be between 1 and 0."
        raise ZeroDivisionError(msg)
    if k < 0 or n <= k:
        msg = "k must be between 0 and n"
        raise ZeroDivisionError(msg)
    dn = n - k
    if k == 0:
        p = -expm1(log1p(y - 1.0) / dn) if y > 0.8 else 1.0 - y ** (1.0 / dn)
    else:
        dk = k + 1
        p = betai(dn, dk, 0.5)
        p = incbi(dk, dn, 1.0 - y) if p > 0.5 else 1.0 - incbi(dn, dk, y)
    return p


def gdtri(a, b, y):
    """Implementation note: Equivalent to: `scipy.stats.gamma.isf(1 - y, a=b, scale=1/a)`."""
    """Returns Gamma such that y is the probability in the integral.

    WARNING: if 1-y == 1, gives incorrect result. The scipy implementation
    gets around this by using cdflib, which is in Fortran. Until someone
    gets around to translating that, only use this function for values of
    p greater than 1e-15 or so!
    """
    y = clip(y, 0, 1)  # ensure y is between 0 and 1
    if y < 0.0 or y > 1.0 or a <= 0.0 or b < 0.0:
        msg = "a and b must be non-negative, y from 0 to 1."
        raise ZeroDivisionError(msg)
    return igami(b, 1.0 - y) / a


def fdtri(a, b, y):
    """Implementation note: Equivalent to: `scipy.stats.f.ppf(y, dfn=a, dfd=b)`."""
    """Returns inverse of F distribution."""
    y = clip(y, 0, 1)  # ensure y is between 0 and 1
    if a < 1.0 or b < 1.0 or y <= 0.0 or y > 1.0:
        msg = "y must be between 0 and 1; a and b >= 1"
        raise ZeroDivisionError(msg)
    y = 1.0 - y
    # Compute probability for x = 0.5
    w = betai(0.5 * b, 0.5 * a, 0.5)
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
    """Implementation note: Equivalent to: quantile positions for plotting, similar to `(i - 0.5) / n`."""
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
    """Implementation note: Equivalent to: vectorized use of `scipy.stats.norm.ppf`, `t.ppf`, `chi2.ppf`, or just linear for uniform."""
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
    funcs = {
        "normal": ndtri,
        "chisq": chi2.isf,
        "t": stdtri,
    }

    if dist != "uniform" and dist not in funcs:
        msg = f"'{dist} not in {list(funcs)}"
        raise ValueError(msg)

    probs = probability_points(n)
    if dist == "uniform":
        return probs

    func = funcs[dist]

    if not args:
        return array([func(p) for p in probs])

    if (
        dist == "chisq"
    ):  # to use scipy.stats.distributions.chi2.isf we need to reverse the order of the arguments
        return array([func(*((p, *args))) for p in probs])

    return array([func(*((*args, p))) for p in probs])
