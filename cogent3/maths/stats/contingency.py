from numpy import zeros, sqrt, outer, log
from numpy.random import shuffle

from cogent3.format.table import rich_html, formatted_cells, simple_format
from cogent3.maths.stats import chisqprob
from cogent3.maths.stats.test import G_ind, G_fit
from cogent3.util.dict_array import DictArray, DictArrayTemplate

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def _get_bin(bins, value):
    """returns bin index corresponding to value"""
    pass


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
def calc_expected(observed, pseudo_count=0):
    """returns the expected array from product of marginal frequencies"""
    if pseudo_count and (observed == 0).any():
        observed = observed.copy()
        observed += pseudo_count

    num_dim = len(observed.shape)
    if num_dim == 2:
        rsum = observed.sum(axis=1)
        rfreq = rsum / rsum.sum()
        csum = observed.sum(axis=0)
        cfreq = csum / csum.sum()
        expecteds = outer(rfreq, cfreq) * rsum.sum()
    elif num_dim == 1:
        expecteds = [observed.mean()] * observed.shape[0]
    else:
        raise NotImplementedError('too many dimensions')
    return expecteds


def calc_chisq(observed, expected):
    """returns the chisq statistic for the two numpy arrays"""
    stat = (observed - expected) ** 2
    stat = (stat / expected).sum()
    return stat


def calc_G(observed, expected, pseudo_count=0, williams=True):
    """returns the G statistic for the two numpy arrays"""
    num_dim = len(observed.shape)
    df = observed.shape[0] - 1
    if num_dim == 2:
        df *= (observed.shape[1] - 1)

    G = 2 * (observed * (log(observed) - log(expected))).sum()
    if williams and num_dim > 1:
        total = observed.sum()
        denom = 6 * total * df
        q = 1 + ((total / observed.sum(axis=0)).sum() - 1) * \
            ((total / observed.sum(axis=1)).sum() - 1) / denom
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
            a DictArray instance, or something that can be converted to one
        expected
            provide in the case where you know the prior proportions, otherwise
            calculated from marginal frequencies
        """
        if not isinstance(observed, DictArray):
            observed = DictArray(observed)

        if expected:
            expected = observed.template.wrap(expected)

        if observed.array.min() < 0 or expected and expected.array.min() < 0:
            raise ValueError('negative values encountered')

        self._observed = observed
        self._expected = expected
        self._residuals = None
        self._df = None
        self.shape = observed.shape

    def _get_repr_(self, html=False):
        def row_cell_func(val, row, col):
            if val in row_labels:
                result = f'<td><b>{val}<b></td>'
            else:
                result = f'<td style="text-align:right">{val}</td>'
            return result

        obs = self.observed.array.tolist()
        exp = self.expected.array.tolist()
        res = self.residuals.array.tolist()

        ndim = len(self.observed.shape)
        if ndim == 1:
            row_labels = 'Observed', 'Expected', 'Residuals'
            col_labels = [str(c) for c in self.observed.template.names[0]]
            rows = []
            # format floats for expecteds and resid
            for row_label, row in zip(row_labels, [obs, exp, res]):
                if row_label == 'Observed':
                    row = [row_label] + [f'{v:,}' for v in row]
                else:
                    row = [row_label] + [f'{v:,.2f}' for v in row]
                rows.append(row)

            if html:
                rows = rich_html(rows, header=[''] + col_labels,
                                 row_cell_func=row_cell_func,
                                 merge_identical=False)
            else:
                header, rows = formatted_cells(
                    rows, header=[''] + col_labels)
                rows = simple_format(header, rows)

        else:
            row_labels = self.observed.template.names[0]
            col_labels = self.observed.template.names[1]
            result = []
            for caption, table in zip(('Observed', 'Expected', 'Residuals'),
                                      (obs, exp, res)):
                rows = []
                for i, r in enumerate(table):
                    if caption == 'Observed':
                        r = [f'{v:,}' for v in r]
                    else:
                        r = [f'{v:,.2f}' for v in r]
                    rows.append([row_labels[i]] + r)
                if html:
                    result.append(rich_html(rows, header=[''] + col_labels,
                                            caption=f'<b>{caption}</b>',
                                            row_cell_func=row_cell_func,
                                            merge_identical=False))
                else:
                    header, rows = formatted_cells(rows,
                                                   header=[''] + col_labels)
                    result.append(simple_format(header, rows, title=caption))
            joiner = '<br>' if html else '\n'
            rows = joiner.join(result)
        return rows

    def _repr_html_(self):
        result = self._get_repr_(html=True)
        return result

    def __repr__(self):
        result = self._get_repr_(html=False)
        return result

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

    @property
    def residuals(self):
        if not self._residuals:
            r = (self.observed.array - self.expected.array)
            r /= sqrt(self.expected.array)
            self._residuals = self.observed.template.wrap(r)
        return self._residuals

    @property
    def df(self):
        if not self._df:
            self._df = self.shape[0] - 1
            if len(self.shape) == 2:
                self._df *= (self.shape[1] - 1)
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
            pval = estimate_pval(self.observed.array, calc_chisq,
                                 num_reps=shuffled)
        return stat, self.df, pval

    def G_independence(self, pseudo_count=0, williams=False, shuffled=0):
        """performs the independence G test
        Parameters
        ----------
        pseudo_count : int
            added to observed to avoid zero division
        shuffled : int
            pvalue is estimated via resampling shuffled times from the observed
            data, preserving the marginals
        """
        assert type(pseudo_count) == int, f'{pseudo_count} not an integer'
        assert type(shuffled) == int, f'{shuffled} not an integer'
        G = calc_G(self.observed.array, self.expected.array,
                   pseudo_count=pseudo_count, williams=williams)
        if not shuffled:
            pval = chisqprob(G, self.df)
        else:
            pval = estimate_pval(self.observed.array, calc_G,
                                 num_reps=shuffled)
        return G, self.df, pval

    def G_fit(self, pseudo_count=0, williams=True):
        """performs the goodness-of-fit G test"""
        assert type(pseudo_count) == int, f'{pseudo_count} not an integer'
        obs = self.observed.array
        if pseudo_count:
            obs += pseudo_count

        G, pval = G_fit(obs, self.expected.array, williams=williams)
        return G, self.df, pval

    def todict(self):
        result = dict(observed=self.observed.todict(),
                      expected=self.expected.todict(),
                      residuals=self.residuals.todict())
        return result

