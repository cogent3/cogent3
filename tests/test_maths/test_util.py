#!/usr/bin/env python
"""Provides tests for array.py
"""
# SUPPORT2425
#from __future__ import with_statement

from warnings import filterwarnings
filterwarnings("ignore", "invalid value encountered in",
               category=RuntimeWarning)

from cogent3.util.unit_test import main, TestCase  # , numpy_err
from cogent3.maths.util import safe_p_log_p, \
    safe_log, row_uncertainty, column_uncertainty,\
    row_degeneracy, column_degeneracy

import numpy
Float = numpy.core.numerictypes.sctype2char(float)
from numpy import array, zeros, transpose, sqrt, reshape, arange, \
    ravel, trace, ones

__author__ = "Rob Knight and Jeremy Widmann"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Sandra Smit"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"


class ArrayMathTests(TestCase):

    def test_safe_p_log_p(self):
        """safe_p_log_p: should handle pos/neg/zero/empty arrays as expected
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
        """safe_log: should handle pos/neg/zero/empty arrays as expected
        """
        # normal valid array
        a = array([[4, 0, 8], [2, 16, 4]])
        self.assertEqual(safe_log(a), array([[2, 0, 3], [1, 4, 2]]))
        # input integers, output floats
        self.assertFloatEqual(
            safe_log(array([1, 2, 3])), array([0, 1, 1.5849625]))
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
        """row_uncertainty: should handle pos/neg/zero/empty arrays as expected
        """
        # normal valid array
        b = transpose(array([[.25, .2, .45, .25, 1], [.25, .2, .45, 0, 0],
                             [.25, .3, .05, .75, 0], [.25, .3, .05, 0, 0]]))
        self.assertFloatEqual(row_uncertainty(
            b), [2, 1.97, 1.47, 0.81, 0], 1e-3)
        # one-dimensional array
        self.assertRaises(ValueError, row_uncertainty,
                          array([.25, .25, .25, .25]))
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
        b = array([[.25, .2, .45, .25, 1], [.25, .2, .45, 0, 0], [.25, .3, .05, .75, 0],
                   [.25, .3, .05, 0, 0]])
        self.assertFloatEqual(column_uncertainty(
            b), [2, 1.97, 1.47, 0.81, 0], 1e-3)
        # one-dimensional array
        self.assertRaises(ValueError, column_uncertainty,
                          array([.25, .25, .25, .25]))
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
        a = array([[.1, .3, .4, .2], [.5, .3, 0, .2], [.8, 0, .1, .1]])
        self.assertEqual(row_degeneracy(a, cutoff=.75), [3, 2, 1])
        self.assertEqual(row_degeneracy(a, cutoff=.95), [4, 3, 3])
        # one-dimensional array
        self.assertRaises(ValueError, row_degeneracy,
                          array([.25, .25, .25, .25]))
        # if cutoff value is not found, results are clipped to the
        # number of columns in the array
        self.assertEqual(row_degeneracy(a, cutoff=2), [4, 4, 4])
        # same behavior on empty array
        self.assertEqual(row_degeneracy(array([[]])), [])

    def test_column_degeneracy(self):
        """column_degeneracy: should work with different cutoff values
        """
        a = array([[.1, .8, .3], [.3, .2, .3], [.6, 0, .4]])
        self.assertEqual(column_degeneracy(a, cutoff=.75), [2, 1, 3])
        self.assertEqual(column_degeneracy(a, cutoff=.45), [1, 1, 2])
        # one-dimensional array
        self.assertRaises(ValueError, column_degeneracy,
                          array([.25, .25, .25, .25]))
        # if cutoff value is not found, results are clipped to the
        # number of rows in the array
        self.assertEqual(column_degeneracy(a, cutoff=2), [3, 3, 3])
        # same behavior on empty array
        self.assertEqual(column_degeneracy(array([[]])), [])


if __name__ == '__main__':
    main()
