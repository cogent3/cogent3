from unittest import TestCase, main

import numpy as np

from numpy import array, diag, dot, exp

import cogent3.maths.matrix_exponentiation as cmme

from cogent3.maths import matrix_exponential_integration as expm


__author__ = "Ben Kaehler"
__copyright__ = "Copyright 2007-2014, The Cogent Project"
__credits__ = ["Ben Kaehler", "Ananias Iliadis", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Ben Kaehler"
__email__ = "benjamin.kaehler@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose


class TestIntegratingExponentiator(TestCase):
    def setUp(self) -> None:
        self.result = 0.7295333
        q = array([[0.5, 0.2, 0.1, 0.2]] * 4)
        for i in range(4):
            q[i, i] = 0.0
            q[i, i] = -sum(
                q[
                    i,
                ]
            )
        self.q = q
        self.p0 = array([0.2, 0.3, 0.3, 0.2])

    def test_van_loan_integrating_exponentiator(self):
        """VanLoanIntegratingExponentiator should reproduce Felsenstein
        analytic result, should throw if we pass it a defected matrix and ask
        it to use CheckedExponentiator, will work with a defective matrix (that
        we can integrate by hand) if we use the default RobustExponentiator,
        and should work for different choices of R and exponentiatior."""
        # Result from Von Bing's R code
        I = expm.VanLoanIntegratingExponentiator(self.q, -diag(self.q))(1.0)
        assert_allclose(dot(self.p0, I), self.result)

        self.assertRaises(
            ArithmeticError,
            expm.VanLoanIntegratingExponentiator,
            self.q,
            -diag(self.q),
            cmme.CheckedExponentiator,
        )

        Q = array([[1.0, 1.0], [0.0, 1.0]])

        def integral(t):
            return array(
                [[exp(t) - 1.0, exp(t) * (t - 1.0) + 1.0], [0.0, exp(t) - 1.0]]
            )

        assert_allclose(expm.VanLoanIntegratingExponentiator(Q)(1.0), integral(1.0))
        assert_allclose(expm.VanLoanIntegratingExponentiator(Q)(2.0), integral(2.0))

        R = array([[1.0], [1.0]])
        assert_allclose(
            expm.VanLoanIntegratingExponentiator(Q, R, cmme.TaylorExponentiator)(1.0),
            dot(integral(1.0), R),
        )

    def test_von_bing_integrating_exponentiator(self):
        """VonBingIntegratingExponentiator should reproduce Felsenstein
        analytic result, should throw if we pass it a defective matrix, and
        should match results obtained from VanLoanIntegratingExponentiator for
        a diagonisable matrix."""
        # Result from Von Bing's R code.
        I = expm.VonBingIntegratingExponentiator(self.q)(1.0)
        assert_allclose(dot(dot(self.p0, I), -diag(self.q)), self.result)

        self.assertRaises(
            ArithmeticError,
            expm.VonBingIntegratingExponentiator,
            array([[1.0, 1.0], [0.0, 1.0]]),
        )

        p = array(
            [
                [0.86758487, 0.05575623, 0.0196798, 0.0569791],
                [0.01827347, 0.93312148, 0.02109664, 0.02750842],
                [0.04782582, 0.1375742, 0.80046869, 0.01413129],
                [0.23022035, 0.22306947, 0.06995306, 0.47675713],
            ]
        )

        assert_allclose(
            expm.VonBingIntegratingExponentiator(p)(1.0),
            expm.VanLoanIntegratingExponentiator(
                p, exponentiator=cmme.FastExponentiator
            )(1.0),
        )
        assert_allclose(
            expm.VonBingIntegratingExponentiator(p)(2.0),
            expm.VanLoanIntegratingExponentiator(
                p, exponentiator=cmme.FastExponentiator
            )(2.0),
        )

    def test_repr(self):
        """repr() works for the integrating exponentiators"""
        for klass in (
            expm.VanLoanIntegratingExponentiator,
            expm.VonBingIntegratingExponentiator,
        ):
            i = klass(self.q)
            g = repr(i)
            self.assertIsInstance(g, str)

    def test_calc_number_subs(self):
        """correctly compute ENS"""
        mprobs = diag([0.1, 0.2, 0.3, 0.4])
        moprobs = array([0.1, 0.2, 0.3, 0.4])

        def get_calibrated_Q(R):
            Q = dot(R, mprobs)
            diag_add = diag(np.sum(Q, axis=1))
            to_divide = np.dot(moprobs, np.sum(Q, axis=1))
            Q -= diag_add
            Q /= to_divide
            return Q

        R = array([[0, 2, 1, 1], [2, 0, 1, 1], [1, 1, 0, 2], [1, 1, 2, 0]], dtype=float)

        Q = get_calibrated_Q(R)
        length = 0.1
        got = expm.expected_number_subs(moprobs, Q, length)
        assert_allclose(got, length)
        # case 2, length != ENS

        A = array(
            [[0, 1, 1, 1], [2, 0, 1, 1], [1, 1, 0, 40], [1, 1, 1, 0]], dtype=float
        )
        Q = get_calibrated_Q(A)
        length = 0.2
        got = expm.expected_number_subs(moprobs, Q, length)
        self.assertNotAlmostEqual(got, length)


if __name__ == "__main__":
    main()
