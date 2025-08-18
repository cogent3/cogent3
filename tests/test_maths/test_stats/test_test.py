"""Unit tests for statistical tests and utility functions."""

from unittest import TestCase

import pytest
from numpy import (
    arange,
    array,
    asarray,
    concatenate,
    fill_diagonal,
    isfinite,
    isnan,
    logical_and,
    ones,
    reshape,
    testing,
    tril,
)
from numpy.testing import assert_allclose, assert_almost_equal, assert_equal

from cogent3.maths.stats.number import NumberCounter
from cogent3.maths.stats.test import (
    ALT_HIGH,
    ALT_LOW,
    ALT_TWO_SIDED,
    G_2_by_2,
    G_fit,
    G_ind,
    ZeroExpectedError,
    _flatten_lower_triangle,
    _get_alternate,
    _get_rank,
    _permute_observations,
    is_symmetric_and_hollow,
    multiple_n,
    mw_boot,
    probability_points,
    safe_sum_p_log_p,
    theoretical_quantiles,
)


def is_prob(value):
    """helper function to establish a 0 <= value <= 1"""
    value = asarray(value)
    return logical_and(value >= 0, value <= 1.0).all()


def similar_means(observed, expected, pvalue=0.01):
    """False if observed p is lower than pvalue"""

    observed, expected = asarray(observed), asarray(expected)

    t, p = t_two_sample(observed, expected)

    # handle case where all elements were the same
    if p is None or not isfinite(p):
        if not observed.shape:
            observed = observed.reshape((1,))

        if not expected.shape:
            expected = expected.reshape((1,))

        if observed[0] == expected[0]:
            return True

    return p > pvalue


class TestsHelper(TestCase):
    """Class with utility methods useful for other tests."""

    def setUp(self):
        """Sets up variables used in the tests."""
        # How many times a p-value should be tested to fall in a given range
        # before failing the test.
        self.p_val_tests = 10

    def assertCorrectPValue(
        self,
        exp_min,
        exp_max,
        fn,
        args=None,
        kwargs=None,
        p_val_idx=0,
    ):
        """Tests that the stochastic p-value falls in the specified range.

        Performs the test self.p_val_tests times and fails if the observed
        p-value does not fall into the specified range at least once. Each
        p-value is also tested that it falls in the range 0.0 to 1.0.

        This method assumes that fn is callable, and will unpack and pass args
        and kwargs to fn if they are provided. It also assumes that fn returns
        a single value (the p-value to be tested) or a tuple of results (any
        length greater than or equal to 1), with the p-value at position
        p_val_idx.

        This is primarily used for testing the Mantel and correlation_test
        functions.
        """
        found_match = False
        for _i in range(self.p_val_tests):
            if args is not None and kwargs is not None:
                obs = fn(*args, **kwargs)
            elif args is not None:
                obs = fn(*args)
            elif kwargs is not None:
                obs = fn(**kwargs)
            else:
                obs = fn()

            try:
                p_val = float(obs)
            except TypeError:
                p_val = obs[p_val_idx]

            assert is_prob(p_val)
            if p_val >= exp_min and p_val <= exp_max:
                found_match = True
                break
        assert found_match


class TestsTests(TestCase):
    """Tests miscellaneous functions."""

    def test_multiple_n(self):
        """multiple_n should swap parameters in multiple_comparisons"""
        assert_allclose(multiple_n(1e-7, 1 - 0.9990005), 10000, rtol=1e-6, atol=1e-6)
        assert_allclose(multiple_n(0.05, 0.4012631), 10, rtol=1e-6, atol=1e-6)
        assert_allclose(multiple_n(1e-20, 1e-20), 1)
        assert_allclose(multiple_n(1e-300, 1e-300), 1)
        assert_allclose(multiple_n(0.95, 0.99987499999999996), 3)
        assert_allclose(multiple_n(0.5, 0.96875), 5)
        assert_allclose(multiple_n(1e-20, 1e-19), 10)

    def test_get_alternate(self):
        """correctly identifies the specified alternate hypothesis"""
        alt = _get_alternate("lo")
        assert alt == ALT_LOW
        alt = _get_alternate("hi")
        assert alt == ALT_HIGH
        alt = _get_alternate("2")
        assert alt == ALT_TWO_SIDED
        with pytest.raises(ValueError):
            _get_alternate("22")


class GTests(TestCase):
    """Tests implementation of the G tests for fit and independence."""

    def test_G_2_by_2_2tailed_equal(self):
        """G_2_by_2 should return 0 if all cell counts are equal"""
        assert_allclose(0, G_2_by_2(1, 1, 1, 1, False, False)[0], atol=1e-12)
        assert_allclose(0, G_2_by_2(100, 100, 100, 100, False, False)[0], atol=1e-12)
        assert_allclose(0, G_2_by_2(100, 100, 100, 100, True, False)[0], atol=1e-12)

    def test_G_2_by_2_bad_data(self):
        """G_2_by_2 should raise ValueError if any counts are negative"""
        self.assertRaises(ValueError, G_2_by_2, 1, -1, 1, 1)

    def test_G_2_by_2_2tailed_examples(self):
        """G_2_by_2 values should match examples in Sokal & Rohlf"""
        # example from p 731, Sokal and Rohlf (1995)
        # without correction
        assert_allclose(G_2_by_2(12, 22, 16, 50, False, False)[0], 1.33249, 0.0001)
        assert_allclose(G_2_by_2(12, 22, 16, 50, False, False)[1], 0.24836, 0.0001)
        # with correction
        assert_allclose(G_2_by_2(12, 22, 16, 50, True, False)[0], 1.30277, 0.0001)
        assert_allclose(G_2_by_2(12, 22, 16, 50, True, False)[1], 0.25371, 0.0001)

    def test_G_2_by_2_1tailed_examples(self):
        """G_2_by_2 values should match values from codon_binding program"""
        # first up...the famous arginine case
        assert_allclose(
            G_2_by_2(36, 16, 38, 106),
            (29.111609, 0),
            rtol=0.00001,
            atol=1e-6,
        )
        # then some other miscellaneous positive and negative values
        assert_allclose(
            G_2_by_2(0, 52, 12, 132),
            (-7.259930, 0.996474),
            rtol=0.00001,
            atol=1e-6,
        )
        assert_allclose(
            G_2_by_2(5, 47, 14, 130),
            (-0.000481, 0.508751),
            rtol=0.00001,
            atol=1e-6,
        )
        assert_allclose(
            G_2_by_2(5, 47, 36, 108),
            (-6.065167, 0.993106),
            rtol=0.00001,
            atol=1e-6,
        )

    def test_Gfit_unequal_lists(self):
        """Gfit should raise errors if lists unequal"""
        # lists must be equal
        self.assertRaises(ValueError, G_fit, [1, 2, 3], [1, 2])

    def test_Gfit_negative_observeds(self):
        """Gfit should raise ValueError if any observeds are negative."""
        self.assertRaises(ValueError, G_fit, [-1, 2, 3], [1, 2, 3])

    def test_Gfit_nonpositive_expecteds(self):
        """Gfit should raise ZeroExpectedError if expecteds are zero/negative"""
        self.assertRaises(ZeroExpectedError, G_fit, [1, 2, 3], [0, 1, 2])
        self.assertRaises(ZeroExpectedError, G_fit, [1, 2, 3], [-1, 1, 2])

    def test_Gfit_good_data(self):
        """Gfit tests for fit should match examples in Sokal and Rohlf"""
        # example from p. 699, Sokal and Rohlf (1995)
        obs = [63, 31, 28, 12, 39, 16, 40, 12]
        exp = [
            67.78125,
            22.59375,
            22.59375,
            7.53125,
            45.18750,
            15.06250,
            45.18750,
            15.06250,
        ]
        # without correction
        assert_allclose(G_fit(obs, exp, False)[0], 8.82397, 0.00002)
        assert_allclose(G_fit(obs, exp, False)[1], 0.26554, 0.00002)
        # with correction
        assert_allclose(G_fit(obs, exp)[0], 8.76938, 0.00002)
        assert_allclose(G_fit(obs, exp)[1], 0.26964, 0.00002)

        # example from p. 700, Sokal and Rohlf (1995)
        obs = [130, 46]
        exp = [132, 44]
        # without correction
        assert_allclose(G_fit(obs, exp, False)[0], 0.12002, 0.00002)
        assert_allclose(G_fit(obs, exp, False)[1], 0.72901, 0.00002)
        # with correction
        assert_allclose(G_fit(obs, exp)[0], 0.11968, 0.00002)
        assert_allclose(G_fit(obs, exp)[1], 0.72938, 0.00002)

    def test_safe_sum_p_log_p(self):
        """safe_sum_p_log_p should ignore zero elements, not raise error"""
        m = array([2, 4, 0, 8])
        assert safe_sum_p_log_p(m, 2) == 2 * 1 + 4 * 2 + 8 * 3

    def test_G_ind(self):
        """G test for independence should match Sokal and Rohlf p 738 values"""
        a = array([[29, 11], [273, 191], [8, 31], [64, 64]])
        assert_allclose(G_ind(a)[0], 28.59642)
        assert_allclose(G_ind(a, True)[0], 28.31244, rtol=1e-5)


class BayesUpdateTests(TestCase):
    """Tests implementation of Bayes calculations"""

    def setUp(self):
        first = [0.25, 0.25, 0.25]
        second = [0.1, 0.75, 0.3]
        third = [0.95, 1e-10, 0.2]
        fourth = [0.01, 0.9, 0.1]
        bad = [1, 2, 1, 1, 1]
        self.bad = [first, bad, second, third]
        self.test = [first, second, third, fourth]
        self.permuted = [fourth, first, third, second]
        self.deleted = [second, fourth, third]
        self.extra = [first, second, first, third, first, fourth, first]

        # BEWARE: low precision in second item, so need to adjust threshold
        # for assertFloatEqual accordingly (and use assertFloatEqualAbs).
        self.result = [0.136690646154, 0.000000009712, 0.863309344133]


class StatTests(TestsHelper):
    """Tests that the t and z tests are implemented correctly"""

    def setUp(self):
        super().setUp()

        self.x = [
            7.33,
            7.49,
            7.27,
            7.93,
            7.56,
            7.81,
            7.46,
            6.94,
            7.49,
            7.44,
            7.95,
            7.47,
            7.04,
            7.10,
            7.64,
        ]

        self.y = [
            7.53,
            7.70,
            7.46,
            8.21,
            7.81,
            8.01,
            7.72,
            7.13,
            7.68,
            7.66,
            8.11,
            7.66,
            7.20,
            7.25,
            7.79,
        ]

    def test_permute_observations(self):
        """Test works correctly on small input dataset."""
        I = [10, 20.0, 1]
        II = [2, 4, 5, 7]
        obs = _permute_observations(I, II, 1)
        assert len(obs[0]) == 1
        assert len(obs[1]) == 1
        assert len(obs[0][0]) == len(I)
        assert len(obs[1][0]) == len(II)
        assert_allclose(sorted(concatenate((obs[0][0], obs[1][0]))), sorted(I + II))

class CorrelationTests(TestsHelper):
    """Tests of correlation coefficients and Mantel test."""

    def setUp(self):
        """Sets up variables used in the tests."""
        super().setUp()

        # For testing spearman and correlation_test using method='spearman'.
        # Taken from the Spearman wikipedia article. Also used for testing
        # Pearson (verified with R).
        self.data1 = [106, 86, 100, 101, 99, 103, 97, 113, 112, 110]
        self.data2 = [7, 0, 27, 50, 28, 29, 20, 12, 6, 17]

        # For testing spearman.
        self.a = [1, 2, 4, 3, 1, 6, 7, 8, 10, 4]
        self.b = [2, 10, 20, 1, 3, 7, 5, 11, 6, 13]
        self.c = [7, 1, 20, 13, 3, 57, 5, 121, 2, 9]
        self.r = (1.7, 10, 20, 1.7, 3, 7, 5, 11, 6.5, 13)
        self.x = (1, 2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)

        # Ranked copies for testing spearman.
        self.b_ranked = [2, 7, 10, 1, 3, 6, 4, 8, 5, 9]
        self.c_ranked = [5, 1, 8, 7, 3, 9, 4, 10, 2, 6]

    def test_is_symmetric_and_hollow(self):
        """Should correctly test for symmetry and hollowness of dist mats."""
        assert is_symmetric_and_hollow(array([[0, 1], [1, 0]]))
        assert is_symmetric_and_hollow(array([[0, 1], [1, 0]]))
        assert is_symmetric_and_hollow(array([[0.0, 0], [0.0, 0]]))
        assert not is_symmetric_and_hollow(array([[0.001, 1], [1, 0]]))
        assert not is_symmetric_and_hollow(array([[0, 1.1], [1, 0]]))
        assert not is_symmetric_and_hollow(array([[0.5, 1.1], [1, 0]]))

    def test_flatten_lower_triangle(self):
        """Test flattening various dms' lower triangulars."""
        assert _flatten_lower_triangle(array([[8]])) == []
        assert _flatten_lower_triangle(array([[1, 2], [3, 4]])) == [3]
        assert _flatten_lower_triangle(array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])) == [
            4,
            7,
            8,
        ]

    def test_get_rank(self):
        """Test the _get_rank function with valid input."""
        exp = (
            [1.5, 3.5, 7.5, 5.5, 1.5, 9.0, 10.0, 11.0, 12.0, 7.5, 14.0, 3.5, 5.5, 13.0],
            4,
        )
        obs = _get_rank(self.x)
        assert_allclose(obs[0], exp[0])
        assert obs[1] == exp[1]

        exp = ([1.5, 3.0, 5.5, 4.0, 1.5, 7.0, 8.0, 9.0, 10.0, 5.5], 2)
        obs = _get_rank(self.a)
        assert_allclose(obs[0], exp[0])
        assert obs[1] == exp[1]

        exp = ([2, 7, 10, 1, 3, 6, 4, 8, 5, 9], 0)
        obs = _get_rank(self.b)
        assert_allclose(obs[0], exp[0])
        assert obs[1] == exp[1]

        exp = ([1.5, 7.0, 10.0, 1.5, 3.0, 6.0, 4.0, 8.0, 5.0, 9.0], 1)
        obs = _get_rank(self.r)
        assert_allclose(obs[0], exp[0])
        assert obs[1] == exp[1]

        exp = ([], 0)
        obs = _get_rank([])
        assert_allclose(obs[0], exp[0])
        assert obs[1] == exp[1]

    def test_get_rank_invalid_input(self):
        """Test the _get_rank function with invalid input."""
        vec = [1, "a", 3, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, {1: 2}, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, [23, 1], 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, (1,), 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

class TestDistMatrixPermutationTest(TestCase):
    """Tests of distance_matrix_permutation_test"""

    def setUp(self):
        """sets up variables for testing"""
        self.matrix = array(
            [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]],
        )
        self.cells = [(0, 1), (1, 3)]
        self.cells2 = [(0, 2), (2, 3)]


    def test_probability_points(self):
        """generates evenly spaced probabilities"""
        expect = (
            0.1190476190476190,
            0.3095238095238095,
            0.5000000000000000,
            0.6904761904761905,
            0.8809523809523809,
        )
        got = probability_points(5)
        assert_almost_equal(got, expect)
        expect = (
            0.04545454545454546,
            0.13636363636363635,
            0.22727272727272727,
            0.31818181818181818,
            0.40909090909090912,
            0.50000000000000000,
            0.59090909090909094,
            0.68181818181818177,
            0.77272727272727271,
            0.86363636363636365,
            0.95454545454545459,
        )
        got = probability_points(11)
        assert_almost_equal(got, expect)

    def test_theoretical_quantiles(self):
        """correctly produce theoretical quantiles"""
        expect = probability_points(4)
        got = theoretical_quantiles(4, dist="uniform")
        assert_almost_equal(got, expect)
        expect = (
            -1.049131397963971,
            -0.299306910465667,
            0.299306910465667,
            1.049131397963971,
        )
        probability_points(4)
        got = theoretical_quantiles(len(expect), dist="normal")
        assert_almost_equal(got, expect)

        # for gamma with shape 2, scale 1/3
        expect = [
            3.833845224364122,
            1.922822334309249,
            0.9636761737854768,
            0.3181293892593747,
        ]
        got = theoretical_quantiles(4, "chisq", df=2)
        assert_almost_equal(got, expect)

        expect = (
            -1.2064470985524887,
            -0.3203979544794824,
            0.3203979544794824,
            1.2064470985524887,
        )
        got = theoretical_quantiles(4, "t", df=4)
        assert_almost_equal(got, expect)


# execute tests if called from command line
