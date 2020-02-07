#!/usr/bin/env python
"""Provides tests for array.py
"""
# SUPPORT2425
# from __future__ import with_statement

from warnings import filterwarnings

import numpy

from numpy import (
    arange,
    array,
    ones,
    ravel,
    reshape,
    sqrt,
    trace,
    transpose,
    zeros,
)

from cogent3.maths.util import (
    column_degeneracy,
    column_uncertainty,
    row_degeneracy,
    row_uncertainty,
    safe_log,
    safe_p_log_p,
)
from cogent3.util.unit_test import TestCase, main  # , numpy_err


filterwarnings("ignore", "invalid value encountered in", category=RuntimeWarning)


Float = numpy.core.numerictypes.sctype2char(float)

__author__ = "Rob Knight and Jeremy Widmann"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Sandra Smit"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"


class ArrayMathTests(TestCase):
    def test_safe_p_log_p(self):
        """safe_p_log_p: should handle pos/neg/zero/empty arrays
        """
        # normal valid array
        a = array([[4, 0, 8], [2, 16, 4]])
        self.assertEqual(safe_p_log_p(a), array([[-8, 0, -24], [-2, -64, -8]]))
        # just zeros
        a = array([[0, 0], [0, 0]])
        self.assertEqual(safe_p_log_p(a), array([[0, 0], [0, 0]]))
        # negative number -- throw error
        with self.assertRaises(FloatingPointError):
            safe_p_log_p(array([-4]))
        # integer input, float output
        self.assertFloatEqual(safe_p_log_p(array([3])), array([-4.75488750]))
        # empty array
        self.assertEqual(safe_p_log_p(array([])), array([]))

    def test_safe_log(self):
        """safe_log: should handle pos/neg/zero/empty arrays
        """
        # normal valid array
        a = array([[4, 0, 8], [2, 16, 4]])
        self.assertEqual(safe_log(a), array([[2, 0, 3], [1, 4, 2]]))
        # input integers, output floats
        self.assertFloatEqual(safe_log(array([1, 2, 3])), array([0, 1, 1.5849625]))
        # just zeros
        a = array([[0, 0], [0, 0]])
        self.assertEqual(safe_log(a), array([[0, 0], [0, 0]]))
        # negative number

        with self.assertRaises(FloatingPointError):
            safe_log(array([0, 3, -4]))

        # empty array
        self.assertEqual(safe_log(array([])), array([]))
        # double empty array
        self.assertEqual(safe_log(array([[]])), array([[]]))

    def test_row_uncertainty(self):
        """row_uncertainty: should handle pos/neg/zero/empty arrays
        """
        # normal valid array
        b = transpose(
            array(
                [
                    [0.25, 0.2, 0.45, 0.25, 1],
                    [0.25, 0.2, 0.45, 0, 0],
                    [0.25, 0.3, 0.05, 0.75, 0],
                    [0.25, 0.3, 0.05, 0, 0],
                ]
            )
        )
        self.assertFloatEqual(row_uncertainty(b), [2, 1.97, 1.47, 0.81, 0], 1e-3)
        # one-dimensional array
        self.assertRaises(ValueError, row_uncertainty, array([0.25, 0.25, 0.25, 0.25]))
        # zeros
        self.assertEqual(row_uncertainty(array([[0, 0]])), array([0]))
        # empty 2D array
        self.assertEqual(row_uncertainty(array([[]])), array([0]))
        self.assertEqual(row_uncertainty(array([[], []])), array([0, 0]))
        # negative number -- throw error
        with self.assertRaises(FloatingPointError):
            row_uncertainty(array([[-2]]))

    def test_col_uncertainty(self):
        """column_uncertainty: should handle pos/neg/zero/empty arrays
        """
        b = array(
            [
                [0.25, 0.2, 0.45, 0.25, 1],
                [0.25, 0.2, 0.45, 0, 0],
                [0.25, 0.3, 0.05, 0.75, 0],
                [0.25, 0.3, 0.05, 0, 0],
            ]
        )
        self.assertFloatEqual(column_uncertainty(b), [2, 1.97, 1.47, 0.81, 0], 1e-3)
        # one-dimensional array
        self.assertRaises(
            ValueError, column_uncertainty, array([0.25, 0.25, 0.25, 0.25])
        )
        # zeros
        self.assertEqual(column_uncertainty(array([[0, 0]])), array([0, 0]))
        # empty 2D array
        self.assertEqual(column_uncertainty(array([[]])), array([]))
        self.assertEqual(column_uncertainty(array([[], []])), array([]))
        # negative number -- throw error
        with self.assertRaises(FloatingPointError):
            column_uncertainty(array([[-2]]))

    def test_row_degeneracy(self):
        """row_degeneracy: should work with different cutoff values and arrays
        """
        a = array([[0.1, 0.3, 0.4, 0.2], [0.5, 0.3, 0, 0.2], [0.8, 0, 0.1, 0.1]])
        self.assertEqual(row_degeneracy(a, cutoff=0.75), [3, 2, 1])
        self.assertEqual(row_degeneracy(a, cutoff=0.95), [4, 3, 3])
        # one-dimensional array
        self.assertRaises(ValueError, row_degeneracy, array([0.25, 0.25, 0.25, 0.25]))
        # if cutoff value is not found, results are clipped to the
        # number of columns in the array
        self.assertEqual(row_degeneracy(a, cutoff=2), [4, 4, 4])
        # same behavior on empty array
        self.assertEqual(row_degeneracy(array([[]])), [])

    def test_column_degeneracy(self):
        """column_degeneracy: should work with different cutoff values
        """
        a = array([[0.1, 0.8, 0.3], [0.3, 0.2, 0.3], [0.6, 0, 0.4]])
        self.assertEqual(column_degeneracy(a, cutoff=0.75), [2, 1, 3])
        self.assertEqual(column_degeneracy(a, cutoff=0.45), [1, 1, 2])
        # one-dimensional array
        self.assertRaises(
            ValueError, column_degeneracy, array([0.25, 0.25, 0.25, 0.25])
        )
        # if cutoff value is not found, results are clipped to the
        # number of rows in the array
        self.assertEqual(column_degeneracy(a, cutoff=2), [3, 3, 3])
        # same behavior on empty array
        self.assertEqual(column_degeneracy(array([[]])), [])


if __name__ == "__main__":
    main()
