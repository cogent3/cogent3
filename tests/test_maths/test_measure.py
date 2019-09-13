from unittest import TestCase, main

from numpy import diag_indices, dot
from numpy.random import random
from numpy.testing import assert_allclose

from cogent3.maths.matrix_exponentiation import PadeExponentiator
from cogent3.maths.matrix_logarithm import logm
from cogent3.maths.measure import (
    jsd,
    jsm,
    paralinear_continuous_time,
    paralinear_discrete_time,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Stephen Ka-Wah Ma"]
__license__ = "BSD-3"
__version__ = "2019.9.13a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


def gen_q_p():
    q1 = random((4, 4))
    indices = diag_indices(4)
    q1[indices] = 0
    q1[indices] = -q1.sum(axis=1)
    p1 = PadeExponentiator(q1)()
    return q1, p1


def gen_qs_ps():
    q1, p1 = gen_q_p()
    q2, p2 = gen_q_p()
    p3 = dot(p1, p2)
    q3 = logm(p3)
    return (q1, p1), (q2, p2), (q3, p3)


def next_pi(pi, p):
    return dot(pi, p)


class ParalinearTest(TestCase):
    def test_paralinear_discrete_time(self):
        """tests paralinear_discrete_time to compare it with the output of paralinear_continuous_time"""
        qp1, qp2, qp3 = gen_qs_ps()
        pi1 = random(4)
        pi1 /= pi1.sum()
        pi2 = next_pi(pi1, qp1[1])
        pi3 = next_pi(pi2, qp2[1])

        con_time_pl1 = paralinear_continuous_time(qp1[1], pi1, qp1[0])
        dis_time_pl1 = paralinear_discrete_time(qp1[1], pi1)
        assert_allclose(con_time_pl1, dis_time_pl1)

        con_time_pl2 = paralinear_continuous_time(qp2[1], pi2, qp2[0])
        dis_time_pl2 = paralinear_discrete_time(qp2[1], pi2)
        assert_allclose(con_time_pl2, dis_time_pl2)

        con_time_pl3 = paralinear_continuous_time(qp3[1], pi3, qp3[0])
        dis_time_pl3 = paralinear_discrete_time(qp3[1], pi3)
        assert_allclose(con_time_pl3, dis_time_pl3)

    def test_paralinear_continuous_time(self):
        """paralinear_continuous_time is additive from random matrices"""
        qp1, qp2, qp3 = gen_qs_ps()
        pi1 = random(4)
        pi1 /= pi1.sum()
        pi2 = next_pi(pi1, qp1[1])

        pl1 = paralinear_continuous_time(qp1[1], pi1, qp1[0])
        pl2 = paralinear_continuous_time(qp2[1], pi2, qp2[0])
        pl3 = paralinear_continuous_time(qp3[1], pi1, qp3[0])

        assert_allclose(pl1 + pl2, pl3)

    def test_paralinear_continuous_time_validate(self):
        """paralinear_continuous_time validate check consistency"""
        qp1, qp2, qp3 = gen_qs_ps()
        pi1 = random(4)

        with self.assertRaises(AssertionError):
            paralinear_continuous_time(
                qp1[1], qp1[0], qp1[0], validate=True
            )  # pi invalid shape

        with self.assertRaises(AssertionError):
            paralinear_continuous_time(
                qp1[1], pi1, qp1[0], validate=True
            )  # pi invalid values

        pi1 /= pi1.sum()
        with self.assertRaises(AssertionError):
            paralinear_continuous_time(qp1[1], pi1, qp1[1], validate=True)  # invalid Q

        with self.assertRaises(AssertionError):
            paralinear_continuous_time(qp1[0], pi1, qp1[0], validate=True)  # invalid P

        qp2[0][0, 0] = 9
        with self.assertRaises(AssertionError):
            paralinear_continuous_time(qp1[1], pi1, qp2[0], validate=True)  # invalid Q

        qp2[1][0, 3] = 9
        with self.assertRaises(AssertionError):
            paralinear_continuous_time(qp2[1], pi1, qp1[0], validate=True)  # invalid P


class TestJensenShannon(TestCase):
    def test_jsd_validation(self):
        """jsd fails with malformed data"""
        freqs1 = random(5)
        normalised_freqs1 = freqs1 / freqs1.sum()
        two_dimensional_freqs1 = [freqs1, freqs1]
        shorter_freqs1 = freqs1[:4]

        freqs2 = random(5)
        normalised_freqs2 = freqs2 / freqs2.sum()
        two_dimensional_freqs2 = [freqs2, freqs2]
        shorter_freqs2 = freqs2[:4]

        with self.assertRaises(AssertionError):
            jsd(
                freqs1, two_dimensional_freqs2, validate=True
            )  # freqs1/freqs2 mismatched shape

        with self.assertRaises(AssertionError):
            jsd(
                two_dimensional_freqs1, freqs2, validate=True
            )  # freqs1/freqs2 mismatched shape

        with self.assertRaises(AssertionError):
            jsd(freqs1, shorter_freqs2, validate=True)  # freqs1/freqs2 mismatched shape

        with self.assertRaises(AssertionError):
            jsd(shorter_freqs1, freqs2, validate=True)  # freqs1/freqs2 mismatched shape

        with self.assertRaises(AssertionError):
            jsd(
                two_dimensional_freqs1, freqs2, validate=True
            )  # freqs1 has incorrect dimension

        with self.assertRaises(AssertionError):
            jsd(
                two_dimensional_freqs1, two_dimensional_freqs2, validate=True
            )  # freqs1 has incorrect dimension

        with self.assertRaises(AssertionError):
            jsd(
                freqs1, two_dimensional_freqs2, validate=True
            )  # freqs2 has incorrect dimension

        with self.assertRaises(AssertionError):
            jsd(freqs1, freqs2, validate=True)  # invalid freqs1

        with self.assertRaises(AssertionError):
            jsd(freqs1, normalised_freqs2, validate=True)  # invalid freqs1

        with self.assertRaises(AssertionError):
            jsd(normalised_freqs1, freqs2, validate=True)  # invalid freqs2

    def test_jsd(self):
        """case1 is testing if the jsd between two identical distributions is 0.0"""
        case1 = [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        for pointer in range(10):
            case1[0][pointer] = 1.0
            case1[1][pointer] = 1.0
            assert_allclose(
                jsd(case1[0], case1[1], validate=True),
                0.0,
                err_msg="Testing case1 for jsd failed",
            )
            case1[0][pointer] = 0.0
            case1[1][pointer] = 0.0
        """case2 is testing the numerical output of jsd between two random distributions"""
        case2 = [[1.0 / 10, 9.0 / 10, 0], [0, 1.0 / 10, 9.0 / 10]]
        assert_allclose(
            jsd(case2[0], case2[1], validate=True),
            0.7655022032053593,
            err_msg="Testing case2 for jsd failed",
        )
        """case3 is testing the numerical output of jsd between two random distributions"""
        case3 = [[1.0, 0.0], [0.5, 0.5]]
        assert_allclose(
            jsd(case3[0], case3[1], validate=True),
            0.3112781244591328,
            err_msg="Testing case3 for jsd failed",
        )
        """case4 is testing if the jsd between two identical uniform distributions is 0.0"""
        case4 = [
            [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        ]
        assert_allclose(
            jsd(case4[0], case4[1], validate=True),
            0.0,
            err_msg="Testing case4 for jsd failed",
        )

    def test_jsm(self):
        """case1 is testing if the jsm between two identical distributions is 0.0"""
        case1 = [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
        for pointer in range(10):
            case1[0][pointer] = 1.0
            case1[1][pointer] = 1.0
            assert_allclose(
                jsm(case1[0], case1[1], validate=True),
                0.0,
                err_msg="Testing case1 for jsm failed",
            )
            case1[0][pointer] = 0.0
            case1[1][pointer] = 0.0
        """case2 is testing the numerical output of jsm between two random distributions"""
        case2 = [[1.0 / 10, 9.0 / 10, 0], [0, 1.0 / 10, 9.0 / 10]]
        assert_allclose(
            jsm(case2[0], case2[1], validate=True),
            0.8749298275892526,
            err_msg="Testing case2 for jsm failed",
        )
        """case3 is testing the numerical output of jsm between two random distributions"""
        case3 = [[1.0, 0.0], [0.5, 0.5]]
        assert_allclose(
            jsm(case3[0], case3[1], validate=True),
            0.5579230452841438,
            err_msg="Testing case3 for jsm failed",
        )
        """case4 is testing if the jsm between two identical uniform distributions is 0.0"""
        case4 = [
            [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        ]
        assert_allclose(
            jsm(case4[0], case4[1], validate=True),
            0.0,
            err_msg="Testing case4 for jsm failed",
        )


if __name__ == "__main__":
    main()
