#!/usr/bin/env python
"""Unit tests for matrix logarithm."""
from unittest import TestCase, main

from numpy import array

from cogent3.maths.matrix_logarithm import (
    is_generator_unique,
    logm,
    logm_taylor,
)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley", "Ben Kaehler"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose


class logarithm_tests(TestCase):
    """Tests of top-level matrix logarithm functions."""

    def test_logm(self):
        """logm results should match scipy's"""
        p = array(
            [
                [0.86758487, 0.05575623, 0.0196798, 0.0569791],
                [0.01827347, 0.93312148, 0.02109664, 0.02750842],
                [0.04782582, 0.1375742, 0.80046869, 0.01413129],
                [0.23022035, 0.22306947, 0.06995306, 0.47675713],
            ]
        )

        q = logm(p)
        assert_allclose(
            q,
            array(
                [
                    [-0.15572053, 0.04947485, 0.01918653, 0.08705915],
                    [0.01405019, -0.07652296, 0.02252941, 0.03994336],
                    [0.05365208, 0.15569116, -0.22588966, 0.01654642],
                    [0.35144866, 0.31279003, 0.10478999, -0.76902868],
                ]
            ),
            rtol=1e-6,
        )

    def test_logm_taylor(self):
        """logm_taylor should return same result as logm"""
        q_eig = logm(
            [
                [0.86758487, 0.05575623, 0.0196798, 0.0569791],
                [0.01827347, 0.93312148, 0.02109664, 0.02750842],
                [0.04782582, 0.1375742, 0.80046869, 0.01413129],
                [0.23022035, 0.22306947, 0.06995306, 0.47675713],
            ]
        )
        q_taylor = logm_taylor(
            [
                [0.86758487, 0.05575623, 0.0196798, 0.0569791],
                [0.01827347, 0.93312148, 0.02109664, 0.02750842],
                [0.04782582, 0.1375742, 0.80046869, 0.01413129],
                [0.23022035, 0.22306947, 0.06995306, 0.47675713],
            ]
        )
        assert_allclose(q_taylor, q_eig)

    def test_is_generator_unique(self):
        """is_generator_unique should identify non-unique primary roots or
        raise a NotImplementedError for non-primary roots"""
        q_fail = array(
            [
                [-9.08955989, 4.50419008, 2.93729567, 1.64807414],
                [3.0820101213, -11.867582855, 3.0196380713, 5.7659346624],
                [3.1293061336, 0.6470353007, -4.340134955, 0.5637935206],
                [9.95494662, 0.63789574, 1.39069539, -11.98353775],
            ]
        )
        self.assertFalse(is_generator_unique(q_fail))

        q_pass = array(
            [
                [-3.764760594, 1.1273556812, 1.3310122018, 1.3063927109],
                [0.920950736, -2.6797373188, 0.4269374722, 1.3318491106],
                [1.1327752022, 1.1606494551, -2.7789984239, 0.4855737666],
                [1.2387180594, 0.1873997167, 0.9202488686, -2.3463666447],
            ]
        )
        self.assertTrue(is_generator_unique(q_pass))

        q_raise = array(
            [
                [-2.77453845e-03, 2.77453817e-03, 1.00809110e-10, 1.74200370e-10],
                [6.81595242e-11, -3.43166912e-10, 1.00810050e-10, 1.74197338e-10],
                [6.81605094e-11, 2.77453817e-03, -2.77453841e-03, 1.74198726e-10],
                [6.81594144e-11, 1.38727059e-10, 1.00808021e-04, -1.00808228e-04],
            ]
        )
        self.assertRaises(NotImplementedError, is_generator_unique, q_raise)

        # currently only support 3x3 or 4x4 matrices
        q_raise = array(
            [
                [1, 2, 3, 4, 5],
                [1, 2, 3, 4, 5],
                [1, 2, 3, 4, 5],
                [1, 2, 3, 4, 5],
                [1, 2, 3, 4, 5],
            ]
        )
        self.assertRaises(NotImplementedError, is_generator_unique, q_raise)


if __name__ == "__main__":
    main()
