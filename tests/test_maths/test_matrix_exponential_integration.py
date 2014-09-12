from cogent.util.unit_test import TestCase, main

from numpy import array, dot, diag, exp

import cogent.maths.matrix_exponentiation as cmme

from cogent.maths import matrix_exponential_integration as expm

__author__ = 'Ben Kaehler'
__copyright__ = "Copyright 2007-2014, The Cogent Project"
__credits__ = ['Ben Kaehler']
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = "Production"

class TestIntegratingExponentiator(TestCase):
    def test_van_loan_integrating_exponentiator(self):
        """VanLoanIntegratingExponentiator should reproduce Felsenstein
        analytic result, should throw if we pass it a defected matrix and ask
        it to use CheckedExponentiator, will work with a defective matrix (that
        we can integrate by hand) if we use the default RobustExponentiator,
        and should work for different choices of R and exponentiatior."""
        # Result from Von Bing's R code
        result = 0.7295333
        q = array([[0.5, 0.2, 0.1, 0.2]]*4)
        for i in range(4):
            q[i, i] = 0.
            q[i, i] = -sum(q[i,])
        p0 = array([0.2, 0.3, 0.3, 0.2])

        I = expm.VanLoanIntegratingExponentiator(q, -diag(q))(1.0)
        self.assertFloatEqual(dot(p0, I), result)

        self.assertRaises(ArithmeticError,
                expm.VanLoanIntegratingExponentiator, 
                q, -diag(q), cmme.CheckedExponentiator)

        Q = array([[1., 1.], [0., 1.]])
        def integral(t):
            return array([[exp(t)-1., exp(t)*(t-1.)+1.], [0., exp(t)-1.]])
        
        self.assertFloatEqual(expm.VanLoanIntegratingExponentiator(Q)(1.),
            integral(1.))
        self.assertFloatEqual(expm.VanLoanIntegratingExponentiator(Q)(2.),
            integral(2.))

        R = array([[1.],[1.]])
        self.assertFloatEqual(expm.VanLoanIntegratingExponentiator(Q, R,
            cmme.TaylorExponentiator)(1.), dot(integral(1.), R))

    def test_von_bing_integrating_exponentiator(self):
        """VonBingIntegratingExponentiator should reproduce Felsenstein 
        analytic result, should throw if we pass it a defective matrix, and
        should match results obtained from VanLoanIntegratingExponentiator for
        a diagonisable matrix."""
        # Result from Von Bing's R code.
        result = 0.7295333
        q = array([[0.5, 0.2, 0.1, 0.2]]*4)
        for i in range(4):
            q[i, i] = 0.
            q[i, i] = -sum(q[i,])
        p0 = array([0.2, 0.3, 0.3, 0.2])

        I = expm.VonBingIntegratingExponentiator(q)(1.0)
        self.assertFloatEqual(dot(dot(p0, I), -diag(q)), result)

        self.assertRaises(ArithmeticError,
                expm.VonBingIntegratingExponentiator,
                array([[1., 1.], [0., 1.]]))

        p = array([[ 0.86758487,  0.05575623,  0.0196798 ,  0.0569791 ],
        [ 0.01827347,  0.93312148,  0.02109664,  0.02750842],
        [ 0.04782582,  0.1375742 ,  0.80046869,  0.01413129],
        [ 0.23022035,  0.22306947,  0.06995306,  0.47675713]])

        self.assertFloatEqual(expm.VonBingIntegratingExponentiator(p)(1.),
                expm.VanLoanIntegratingExponentiator(p,
                    exponentiator=cmme.FastExponentiator)(1.))
        self.assertFloatEqual(expm.VonBingIntegratingExponentiator(p)(2.),
                expm.VanLoanIntegratingExponentiator(p, 
                    exponentiator=cmme.FastExponentiator)(2.))



if __name__ == '__main__':
    main()
