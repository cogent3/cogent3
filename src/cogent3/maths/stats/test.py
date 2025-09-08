"""Provides standard statistical tests. Tests produce statistic and P-value."""

import numpy.typing as npt
from numpy import (
    array,
    finfo,
    float64,
    log,
    nonzero,
    ravel,
    take,
)
from numpy import sum as npsum
from numpy.random import randint
from scipy.special import ndtri
from scipy.stats import ks_2samp, mannwhitneyu, t
from scipy.stats.distributions import chi2

# defining globals for the alternate hypotheses
ALT_TWO_SIDED = "2"
ALT_LOW = "low"
ALT_HIGH = "high"


class IndexOrValueError(IndexError, ValueError):
    pass


class ZeroExpectedError(ValueError):
    """Class for handling tests where an expected value was zero."""


def G_2_by_2(a, b, c, d, williams=1, directional=1):
    """G test for independence in a 2 x 2 table.

    Usage: G, prob = G_2_by_2(a, b, c, d, willliams, directional)

    Cells are in the order:

        a b
        c d

    a, b, c, and d can be int, float, or long.
    williams is a boolean stating whether to do the Williams correction.
    directional is a boolean stating whether the test is 1-tailed.

    Briefly, computes sum(f ln f) for cells - sum(f ln f) for
    rows and columns + f ln f for the table.

    Always has 1 degree of freedom

    To generalize the test to r x c, use the same protocol:
    2*(cells - rows/cols + table), then with (r-1)(c-1) df.

    Note that G is always positive: to get a directional test,
    the appropriate ratio (e.g. a/b > c/d) must be tested
    as a separate procedure. Find the probability for the
    observed G, and then either halve or halve and subtract from
    one depending on whether the directional prediction was
    upheld.

    The default test is now one-tailed (Rob Knight 4/21/03).

    See Sokal & Rohlf (1995), ch. 17. Specifically, see box 17.6 (p731).
    """
    cells = [a, b, c, d]
    n = npsum(cells)
    # return 0 if table was empty
    if not n:
        return (0, 1)
    # raise error if any counts were negative
    if min(cells) < 0:
        msg = "G_2_by_2 got negative cell counts(s): must all be >= 0."
        raise ValueError(msg)

    G = 0
    # Add x ln x for items, adding zero for items whose counts are zero
    for i in [_f for _f in cells if _f]:
        G += i * log(i)
    # Find totals for rows and cols
    ab = a + b
    cd = c + d
    ac = a + c
    bd = b + d
    rows_cols = [ab, cd, ac, bd]
    # exit if we are missing a row or column entirely: result counts as
    # never significant
    if min(rows_cols) == 0:
        return (0, 1)
    # Subtract x ln x for rows and cols
    for i in [_f for _f in rows_cols if _f]:
        G -= i * log(i)
    # Add x ln x for table
    G += n * log(n)
    # Result needs to be multiplied by 2
    G *= 2

    # apply Williams correction
    if williams:
        q = 1 + ((((n / ab) + (n / cd)) - 1) * (((n / ac) + (n / bd)) - 1)) / (6 * n)
        G /= q

    p = chi2.sf(max(G, 0), 1)

    # find which tail we were in if the test was directional
    if directional:
        is_high = (b == 0) or (d != 0 and (a / b > c / d))
        p = p / 2 if is_high else 1 - p / 2
        if not is_high:
            G = -G
    return G, p


def safe_sum_p_log_p(a, base=None):
    """Calculates p * log(p) safely for an array that may contain zeros."""
    # TODO rewrite using numpy funcs to be more succinct
    flat = ravel(a)
    nz = take(flat, nonzero(flat)[0])
    logs = log(nz)
    if base:
        logs /= log(base)
    return npsum(nz * logs, 0)


def G_ind(m, williams=False):
    """Returns G test for independence in an r x c table.

    Requires input data as a numpy array. From Sokal and Rohlf p 738.
    """
    # TODO rewrite using numpy funcs to be more succinct
    f_ln_f_elements = safe_sum_p_log_p(m)
    f_ln_f_rows = safe_sum_p_log_p(npsum(m, 0))
    f_ln_f_cols = safe_sum_p_log_p(npsum(m, 1))
    tot = npsum(ravel(m))
    f_ln_f_table = tot * log(tot)

    df = (len(m) - 1) * (len(m[0]) - 1)
    G = 2 * (f_ln_f_elements - f_ln_f_rows - f_ln_f_cols + f_ln_f_table)
    if williams:
        q = 1 + (
            (tot * npsum(1.0 / npsum(m, 1)) - 1)
            * (tot * npsum(1.0 / npsum(m, 0)) - 1)
            / (6 * tot * df)
        )
        G = G / q
    return G, chi2.sf(max(G, 0), df)


def calc_contingency_expected(matrix):
    """Calculates expected frequencies from a table of observed frequencies

    The input matrix is a dict2D object and represents a frequency table
    with different variables in the rows and columns. (observed
    frequencies as values)

    The expected value is calculated with the following equation:
        Expected = row_total x column_total / overall_total

    The returned matrix (dict2D) has lists of the observed and the
    expected frequency as values
    """
    # transpose matrix for calculating column totals
    t_matrix = matrix.copy()
    t_matrix.transpose()

    overall_total = npsum(list(matrix.Items))
    # make new matrix for storing results
    result = matrix.copy()

    # populate result with expected values
    for row in matrix:
        row_sum = npsum(list(matrix[row].values()))
        for item in matrix[row]:
            column_sum = npsum(list(t_matrix[item].values()))
            # calculate expected frequency
            Expected = (row_sum * column_sum) / overall_total
            result[row][item] = [result[row][item]]
            result[row][item].append(Expected)
    return result


def G_fit(obs, exp, williams=1):
    """G test for fit between two lists of counts.

    Usage: test, prob = G_fit(obs, exp, williams)

    obs and exp are two lists of numbers.
    williams is a boolean stating whether to do the Williams correction.

    SUM(2 f(obs)ln (f(obs)/f(exp)))

    See Sokal and Rohlf chapter 17.
    """
    obs = array(obs)
    exp = array(exp)
    if obs.shape != exp.shape:
        msg = "requires data with equal dimensions."
        raise ValueError(msg)
    if (obs < 0).any():
        msg = "requires all observed values to be positive."
        raise ValueError(msg)
    if (exp == 0).any() or (exp < 0).any():
        msg = "requires all expected values to be positive"
        raise ZeroExpectedError(msg)
    non_zero = obs != 0
    G = 2 * (obs[non_zero] * (log(obs[non_zero]) - log(exp[non_zero]))).sum()
    k = len(obs)
    if williams:
        q = 1 + (k + 1) / (6 * obs.sum())
        G /= q

    return G, chi2.sf(G, k - 1)


def _get_bootstrap_sample(x, y, num_reps):
    combined = array(list(x) + list(y))
    total_obs = len(combined)
    num_x = len(x)
    for _ in range(num_reps):
        indices = randint(0, total_obs, total_obs)
        sampled = combined.take(indices)
        yield sampled[:num_x], sampled[num_x:]


def ks_boot(x, y, alt="two-sided", num_reps=1000):
    """Monte Carlo bootstrap variant of KS test (SciPy 1.16 compatible)."""
    if alt != "two-sided":
        raise NotImplementedError("Only 'two-sided' tests are supported in SciPy 1.16")

    tol = finfo(float64).eps * 100
    result = ks_2samp(x, y)
    observed_stat = result.statistic  # type: ignore[attr-defined]
    num_greater = 0
    for sampled_x, sampled_y in _get_bootstrap_sample(x, y, num_reps):
        result = ks_2samp(sampled_x, sampled_y)
        sample_stat = result.statistic  # type: ignore[attr-defined]
        if sample_stat >= (observed_stat - tol):
            num_greater += 1
    return observed_stat, num_greater / num_reps


def mw_boot(x, y, num_reps=1000):
    """Monte Carlo (bootstrap) variant of the Mann-Whitney test.

    Parameters
    ----------
    x, y
        vectors of numbers
    num_reps
        number of replicates for the  bootstrap

    Notes
    -----
    Uses the same Monte-Carlo resampling code as kw_boot
    """
    tol = finfo(float64).eps * 100

    # Get observed U statistic (taking the max to match previous behavior)
    result = mannwhitneyu(x, y, alternative="two-sided")
    observed_stat = result.statistic

    num_greater = 0
    for sampled_x, sampled_y in _get_bootstrap_sample(x, y, num_reps):
        sample_result = mannwhitneyu(sampled_x, sampled_y, alternative="two-sided")
        sample_stat = sample_result.statistic
        if sample_stat >= (observed_stat - tol):
            num_greater += 1

    return observed_stat, num_greater / num_reps


def probability_points(n) -> npt.NDArray:
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


def theoretical_quantiles(n, dist, **kwargs) -> npt.NDArray:
    """returns theoretical quantiles from dist

    Parameters
    ----------
    n
        number of elements
    dist
        one of 'normal', 'chisq', 't', 'uniform'
    kwargs
        additional keyword arguments (eg, df=2) to pass to the scipy distribution function

    Notes
    -----
    For details on kwargs see the documentation for the scipy functions
    `scipy.stats.norm.ppf`, `scipy.stats.t.ppf`, and `scipy.stats.chi2.ppf`.

    Returns
    -------
    Numpy array of quantiles
    """

    dist_name = dist.lower()
    funcs = {
        "normal": ndtri,
        "chisq": chi2.isf,
        "t": t.ppf,
    }

    if dist_name != "uniform" and dist_name not in funcs:
        msg = f"{dist_name!r} not in {list(funcs)}"
        raise ValueError(msg)

    probs = probability_points(n)
    if dist_name == "uniform":
        return probs

    func = funcs[dist_name]

    return array([func(p, **kwargs) for p in probs])
