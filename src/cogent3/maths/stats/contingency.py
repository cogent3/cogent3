from numpy import log, outer, sqrt, zeros
from numpy.random import shuffle
from numpy.testing import assert_allclose

from cogent3.maths.stats import chisqprob
from cogent3.maths.stats.test import G_fit
from cogent3.util.dict_array import DictArray


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


# todo this should probably go into different module
def shuffled_matrix(matrix):
    """returns a randomly sampled matrix with same marginals"""
    # SLOW algorithm
    expanded_row = []
    for i, row_total in enumerate(matrix.sum(axis=1)):
        expanded_row.extend([i] * row_total)
    expanded_col = []
    for i, col_total in enumerate(matrix.sum(axis=0)):
        expanded_col.extend([i] * col_total)

    shuffled = zeros(matrix.shape, dtype=int)
    shuffle(expanded_col)
    for i, j in zip(expanded_row, expanded_col):
        shuffled[i, j] += 1
    return shuffled


# todo following functions should be moved into stats.test and replace
# or merge with the older implementations
def calc_expected(observed):
    """returns the expected array from product of marginal frequencies"""
    if observed.ndim == 1 or (observed.ndim == 2 and 1 in observed.shape):
        expecteds = zeros(observed.shape, dtype=float)
        expecteds.fill(observed.mean())
    elif observed.ndim == 2:
        rsum = observed.sum(axis=1)
        rfreq = rsum / rsum.sum()
        csum = observed.sum(axis=0)
        cfreq = csum / csum.sum()
        expecteds = outer(rfreq, cfreq) * rsum.sum()
    else:
        raise NotImplementedError("too many dimensions")
    return expecteds


def calc_chisq(observed, expected):
    """returns the chisq statistic for the two numpy arrays"""
    stat = (observed - expected) ** 2
    stat = (stat / expected).sum()
    return stat


def calc_G(observed, expected, williams=True):
    """returns the G statistic for the two numpy arrays

    Parameters
    ----------
    observed : numpy.ndarray
        Observed counts
    expected : numpy.ndarray
        Expected values
    williams : bool
        Applies Williams correction for small sample size

    Returns
    -------
    G statistic
    """
    num_dim = len(observed.shape)
    df = observed.shape[0] - 1
    if num_dim == 2:
        df *= observed.shape[1] - 1

    non_zero = observed != 0
    if not non_zero.all():
        G = (
            2
            * (
                observed[non_zero] * (log(observed[non_zero]) - log(expected[non_zero]))
            ).sum()
        )
    else:
        G = 2 * (observed * (log(observed) - log(expected))).sum()
    if williams and num_dim > 1:
        total = observed.sum()
        denom = 6 * total * df
        q = (
            1
            + ((total / observed.sum(axis=0)).sum() - 1)
            * ((total / observed.sum(axis=1)).sum() - 1)
            / denom
        )
        G /= q
    return G


def estimate_pval(observed, stat_func, num_reps=1000):
    """returns p-value from resampling of observed

    valid for 2D categorical data only"""
    expected = calc_expected(observed)
    obs_stat = stat_func(observed, expected)
    num_gt = 0
    for i in range(num_reps):
        resamp_obs = shuffled_matrix(observed)
        resamp_exp = calc_expected(resamp_obs)
        resamp_stat = stat_func(resamp_obs, resamp_exp)
        if resamp_stat >= obs_stat:
            num_gt += 1
    return num_gt / num_reps


def _astype(data, dtype):
    """returns numpy array of correct type, raises TypeError if fails"""
    converted = data.astype(dtype)
    try:
        assert_allclose(
            converted.tolist(),
            data.tolist(),
        )
    except AssertionError:
        msg = f"could not reliably be converted to {dtype} from dtype={data.dtype}"
        raise TypeError(msg)

    return converted


class CategoryCounts:
    """CategoryCounts for performing contingency tests

    Has attributes for observed, expected and residuals.
    The latter is calculated using the G-test, for goodness-of-fit if expecteds
    are provided, G-test of independence if not provided.
    """

    def __init__(self, observed, expected=None):
        """Parameters
        -------------
        observed
            a DictArray instance, or something that can be converted to one.
            Values must be integers.
        expected
            provide in the case where you know the prior proportions, otherwise
            calculated from marginal frequencies
        """
        if not isinstance(observed, DictArray):
            observed = DictArray(observed)

        # make sure values are int
        observed.array = _astype(observed.array, int)

        if observed.array.sum() == 0:
            raise ValueError("at least one value must be > 0")

        if observed.array.min() < 0:
            raise ValueError("negative values encountered")

        if observed.array.ndim > 2:
            raise NotImplementedError("not designed for >2D")

        self._observed = observed
        self.expected = expected
        self._residuals = None
        self._df = None
        self.shape = observed.shape

    def _get_repr_(self, html=False):
        obs = self.observed.to_table()
        obs.title = "Observed"
        exp = self.expected.to_table()
        exp.title = "Expected"
        exp.digits = 2
        res = self.residuals.to_table()
        res.title = "Residuals"
        res.digits = 2

        ndim = self.observed.array.ndim
        if ndim == 1:
            result = obs.appended("", exp, res, title=None, digits=2)
            if html:
                result.set_repr_policy(show_shape=False)
                result = result._repr_html_()
            else:
                result, _, _ = result._get_repr_()
                result = str(result)
            return result

        result = []
        for t in (obs, exp, res):
            t.set_repr_policy(show_shape=False)
            if html:
                t = t._repr_html_()
            else:
                t, _, _ = t._get_repr_()
                t = str(t)

            result.append(t)

        joiner = "<br>" if html else "\n"
        return joiner.join(result)

    def _repr_html_(self):
        return self._get_repr_(html=True)

    def __repr__(self):
        return self._get_repr_(html=False)

    def __str__(self):
        return self._get_repr_(html=False)

    @property
    def observed(self):
        return self._observed

    @property
    def expected(self):
        if not self._expected:
            expecteds = calc_expected(self.observed.array)
            expecteds = self.observed.template.wrap(expecteds)
            self._expected = expecteds

        return self._expected

    @expected.setter
    def expected(self, expected):
        if expected is None:
            self._expected = None
            return

        expected = self.observed.template.wrap(expected)
        expected.array = _astype(expected.array, float)

        if expected.array.min() < 0:
            raise ValueError("negative values encountered")

        assert_allclose(
            self.observed.array.sum(), expected.array.sum()
        ), "unequal totals of observed and expected"

        self._expected = expected

    @property
    def residuals(self):
        if not self._residuals:
            r = self.observed.array - self.expected.array
            r /= sqrt(self.expected.array)
            self._residuals = self.observed.template.wrap(r)
        return self._residuals

    @property
    def df(self):
        if not self._df:
            self._df = self.shape[0] - 1
            if len(self.shape) == 2:
                self._df *= self.shape[1] - 1
        return self._df

    def chisq_test(self, shuffled=0):
        """performs the chisq test
        Parameters
        ----------
        shuffled
            pvalue is estimated via resampling from the observed data,
            preserving the marginals
        """
        stat = calc_chisq(self.observed.array, self.expected.array)
        if not shuffled:
            pval = chisqprob(stat, self.df)
        else:
            pval = estimate_pval(self.observed.array, calc_chisq, num_reps=shuffled)
        title = "Chisq-test for independence"
        return TestResult(
            self.observed,
            self.expected,
            self.residuals,
            "chisq",
            stat,
            self.df,
            pval,
            test_name=title,
        )

    def G_independence(self, pseudo_count=0, williams=True, shuffled=0):
        """performs the independence G test
        Parameters
        ----------
        pseudo_count : int
            added to observed to avoid zero division
        williams : bool
            Applies Williams correction for small sample size
        shuffled : int
            pvalue is estimated via resampling shuffled times from the observed
            data, preserving the marginals
        """
        assert type(pseudo_count) == int, f"{pseudo_count} not an integer"
        obs = self.observed
        exp = self.expected
        if pseudo_count and (obs.array == 0).any():
            obs = obs.template.wrap(obs.array + pseudo_count)
            exp = calc_expected(obs.array)
            exp = obs.template.wrap(exp)

        assert type(shuffled) == int, f"{shuffled} not an integer"
        G = calc_G(
            obs.array,
            exp.array,
            williams=williams,
        )
        if not shuffled:
            pval = chisqprob(G, self.df)
        else:
            pval = estimate_pval(obs.array, calc_G, num_reps=shuffled)

        title = "G-test for independence"
        amendments = ""
        if pseudo_count:
            amendments = f"pseudo_count={pseudo_count}, "

        if williams:
            amendments = f"{amendments}Williams correction"

        if amendments:
            title = f"{title} (with {amendments})"

        return TestResult(
            obs,
            exp,
            self.residuals,
            "G",
            G,
            self.df,
            pval,
            test_name=title,
        )

    def G_fit(self, williams=True):
        """performs the goodness-of-fit G test

        Parameters
        ----------
        williams : bool
            Applies Williams correction for small sample size
        """
        obs = self.observed.array
        G, pval = G_fit(obs.flatten(), self.expected.array.flatten(), williams=williams)
        title = "G-test goodness-of-fit"
        if williams:
            title = f"{title} (with Williams correction)"

        return TestResult(
            self.observed,
            self.expected,
            self.residuals,
            "G",
            G,
            self.df,
            pval,
            test_name=title,
        )

    def to_dict(self):
        return dict(
            observed=self.observed.to_dict(),
            expected=self.expected.to_dict(),
            residuals=self.residuals.to_dict(),
        )


class TestResult:
    """result of a contingency test"""

    def __init__(
        self, observed, expected, residuals, stat_name, stat, df, pvalue, test_name=""
    ):
        """
        Parameters
        ----------
        observed
            observed counts as a DictArray
        expected
            expected counts as a DictArray
        residuals
            Pearson residuals between observed and expected as a DictArray
        stat_name
            Name of the statistic, e.g. G, chisq.
        stat : float
            value of the statistic
        df : int
            degrees of freedom for the hypothesis test
        pvalue
            pvalue from the hypothesis test
        test_name
            name of the statistical test
        """
        self.pvalue = pvalue
        self.df = df
        self.stat = stat
        self.test_name = test_name
        self.stat_name = stat_name
        self.residuals = residuals
        self.expected = expected
        self.observed = observed
        setattr(self, stat_name, stat)

    def _get_repr_(self):
        from cogent3.util.table import Table

        header = [str(self.stat_name), "df", "pvalue"]
        col_templates = {
            str(self.stat_name): "%.3f",
            "df": "%s",
            "pvalue": "%.4f" if self.pvalue > 1e-3 else "%.2e",
        }
        table = Table(
            header,
            [[self.stat, self.df, self.pvalue]],
            title=self.test_name,
            column_templates=col_templates,
        )
        table.set_repr_policy(show_shape=False)
        return table

    def __repr__(self):
        result = str(self._get_repr_())
        components = CategoryCounts(
            self.observed.to_dict(), expected=self.expected.to_dict()
        )
        result = [result, str(components)]
        return "\n".join(result)

    def __str__(self):
        return repr(self)

    def _repr_html_(self):
        table = self._get_repr_()
        table.set_repr_policy(show_shape=False)
        components = CategoryCounts(
            self.observed.to_dict(), expected=self.expected.to_dict()
        )
        html = [table._repr_html_()]
        html.append(components._repr_html_())
        return "\n".join(html)

    @property
    def statistics(self):
        """returns Table of stat, df and p-value"""
        return self._get_repr_()
