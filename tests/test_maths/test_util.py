"""Provides tests for array.py"""

# SUPPORT2425
# from __future__ import with_statement

from unittest import TestCase
from warnings import filterwarnings

import pytest
from numpy import array, transpose
from numpy.testing import assert_allclose, assert_equal

from cogent3.maths.util import (
    safe_log,
    safe_p_log_p,
)

filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)


class ArrayMathTests(TestCase):
    def test_safe_p_log_p(self):
        """safe_p_log_p: should handle pos/neg/zero/empty arrays"""
        # normal valid array
        a = array([[4, 0, 8], [2, 16, 4]])
        assert_equal(safe_p_log_p(a), array([[-8, 0, -24], [-2, -64, -8]]))
        # just zeros
        a = array([[0, 0], [0, 0]])
        assert_equal(safe_p_log_p(a), array([[0, 0], [0, 0]]))
        # negative number -- throw error
        with pytest.raises(FloatingPointError):
            safe_p_log_p(array([-4]))
        # integer input, float output
        assert_allclose(safe_p_log_p(array([3])), array([-4.75488750]))
        # empty array
        assert_equal(safe_p_log_p(array([])), array([]))

    def test_safe_log(self):
        """safe_log: should handle pos/neg/zero/empty arrays"""
        # normal valid array
        a = array([[4, 0, 8], [2, 16, 4]])
        assert_equal(safe_log(a), array([[2, 0, 3], [1, 4, 2]]))
        # input integers, output floats
        assert_allclose(safe_log(array([1, 2, 3])), array([0, 1, 1.5849625]))
        # just zeros
        a = array([[0, 0], [0, 0]])
        assert_equal(safe_log(a), array([[0, 0], [0, 0]]))
        # negative number

        with pytest.raises(FloatingPointError):
            safe_log(array([0, 3, -4]))

        # empty array
        assert_equal(safe_log(array([])), array([]))
        # double empty array
        assert_equal(safe_log(array([[]])), array([[]]))

class TestUtils(TestCase):
    def test_proportions_and_ratios(self):
        """interconverts proportions and ratios"""
        from cogent3.maths.util import (
            proportions_to_ratios,
            ratios_to_proportions,
        )

        probs = array([0.3, 0.1, 0.1, 0.5])
        ratios = proportions_to_ratios(probs)
        assert_allclose(ratios, [0.6 / 0.4, 0.1 / 0.3, 0.5 / 0.1])

        probs = array([0.3, 0.1, 0.6])
        ratios = proportions_to_ratios(probs)
        assert_allclose(ratios, [0.7 / 0.3, 0.6 / 0.1])

        got = ratios_to_proportions(1, ratios)
        assert_allclose(got, probs)

        probs = array([0.3, 0.1, -0.1, 0.5])
        with pytest.raises(AssertionError):
            proportions_to_ratios(probs)

        probs = array([0.3, 0.1, 0.0, 0.5])
        with pytest.raises(AssertionError):
            proportions_to_ratios(probs)

        with pytest.raises(AssertionError):
            ratios_to_proportions(1.0, [2.3, 1.1, -0.3])
