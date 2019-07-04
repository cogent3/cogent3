from unittest import TestCase, main

from numpy import diag_indices, dot
from numpy.random import random
from numpy.testing import assert_allclose

from cogent3.maths.matrix_exponentiation import PadeExponentiator
from cogent3.maths.matrix_logarithm import logm
from cogent3.maths.measure import paralinear


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
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
    def test_paralinear(self):
        """paralinear is additive from random matrices"""
        qp1, qp2, qp3 = gen_qs_ps()
        pi1 = random(4)
        pi1 /= pi1.sum()
        pi2 = next_pi(pi1, qp1[1])

        pl1 = paralinear(qp1[0], qp1[1], pi1)
        pl2 = paralinear(qp2[0], qp2[1], pi2)
        pl3 = paralinear(qp3[0], qp3[1], pi1)

        assert_allclose(pl1 + pl2, pl3)

    def test_paralinear_validate(self):
        """paralinear validate check consistency"""
        qp1, qp2, qp3 = gen_qs_ps()
        pi1 = random(4)
        with self.assertRaises(AssertionError):
            paralinear(qp1[0], qp1[1], qp1[0], validate=True)  # pi invalid shape

        with self.assertRaises(AssertionError):
            paralinear(qp1[0], qp1[1], pi1, validate=True)  # pi invalid values

        pi1 /= pi1.sum()
        with self.assertRaises(AssertionError):
            paralinear(qp1[1], qp1[1], pi1, validate=True)  # invalid Q

        with self.assertRaises(AssertionError):
            paralinear(qp1[0], qp1[0], pi1, validate=True)  # invalid P

        qp2[0][0, 0] = 9
        with self.assertRaises(AssertionError):
            paralinear(qp2[0], qp1[1], pi1, validate=True)  # invalid Q

        qp2[1][0, 3] = 9
        with self.assertRaises(AssertionError):
            paralinear(qp1[0], qp2[1], pi1, validate=True)  # invalid P


if __name__ == "__main__":
    main()
