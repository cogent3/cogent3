from unittest import TestCase, main

from numpy import arange, array, convolve, exp, float64, pi, random, sin, zeros

from cogent3.maths.period import _autocorr_inner2 as py_autocorr_inner
from cogent3.maths.period import _goertzel_inner as py_goertzel_inner
from cogent3.maths.period import _ipdft_inner2 as py_ipdft_inner
from cogent3.maths.period import auto_corr, dft, goertzel, hybrid, ipdft
from cogent3.maths.period_numba import autocorr_inner as numba_autocorr_inner
from cogent3.maths.period_numba import goertzel_inner as numba_goertzel_inner
from cogent3.maths.period_numba import ipdft_inner as numba_ipdft_inner


__author__ = "Hua Ying, Julien Epps and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Julien Epps", "Hua Ying", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose, assert_almost_equal, assert_equal


class TestPeriod(TestCase):
    def setUp(self):
        t = arange(0, 10, 0.1)
        n = random.randn(len(t))
        nse = convolve(n, exp(-t / 0.05)) * 0.1
        nse = nse[: len(t)]
        self.sig = sin(2 * pi * t) + nse
        self.p = 10

    def test_inner_funcs(self):
        """python and pyrexed implementation should be the same"""
        x = array(
            [
                0.04874203,
                0.56831373,
                0.94267804,
                0.95664485,
                0.60719478,
                -0.09037356,
                -0.69897319,
                -1.11239811,
                -0.84127485,
                -0.56281126,
                0.02301213,
                0.56250284,
                1.0258557,
                1.03906527,
                0.69885916,
                0.10103556,
                -0.43248024,
                -1.03160503,
                -0.84901545,
                -0.84934356,
                0.00323728,
                0.44344594,
                0.97736748,
                1.01635433,
                0.38538423,
                0.09869918,
                -0.60441861,
                -0.90175391,
                -1.00166887,
                -0.66303249,
                -0.02070569,
                0.76520328,
                0.93462426,
                0.97011673,
                0.63199999,
                0.0764678,
                -0.55680168,
                -0.92028808,
                -0.98481451,
                -0.57600588,
                0.0482667,
                0.57572519,
                1.02077883,
                0.93271663,
                0.41581696,
                -0.07639671,
                -0.71426286,
                -0.97730119,
                -1.0370596,
                -0.67919572,
                0.03779302,
                0.60408759,
                0.87826068,
                0.79126442,
                0.69769622,
                0.01419442,
                -0.42917556,
                -1.00100485,
                -0.83945546,
                -0.55746313,
                0.12730859,
                0.60057659,
                0.98059721,
                0.83275501,
                0.69031804,
                0.02277554,
                -0.63982729,
                -1.23680355,
                -0.79477887,
                -0.67773375,
                -0.05204714,
                0.51765381,
                0.77691955,
                0.8996709,
                0.5153137,
                0.01840839,
                -0.65124866,
                -1.13269058,
                -0.92342177,
                -0.45673709,
                0.11212881,
                0.50153941,
                1.09329507,
                0.96457193,
                0.80271578,
                -0.0041043,
                -0.81750772,
                -0.99259986,
                -0.92343788,
                -0.57694955,
                0.13982059,
                0.56653375,
                0.82217563,
                0.85162513,
                0.3984116,
                -0.18937514,
                -0.65304629,
                -1.0067146,
                -1.0037422,
                -0.68011283,
            ]
        )
        N = 100
        period = 10
        assert_allclose(
            py_goertzel_inner(x, N, period), numba_goertzel_inner(x, N, period)
        )

        ulim = 8
        N = 8
        x = array([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0])
        X = zeros(8, dtype="complex128")
        W = array(
            [
                1.00000000e00 + 2.44929360e-16j,
                -1.00000000e00 - 1.22464680e-16j,
                -5.00000000e-01 - 8.66025404e-01j,
                6.12323400e-17 - 1.00000000e00j,
                3.09016994e-01 - 9.51056516e-01j,
                5.00000000e-01 - 8.66025404e-01j,
                6.23489802e-01 - 7.81831482e-01j,
                7.07106781e-01 - 7.07106781e-01j,
            ]
        )
        py_result = py_ipdft_inner(x, X, W, ulim, N)
        numba_result = numba_ipdft_inner(x, X, W, ulim, N)
        for i, j in zip(py_result, numba_result):
            assert_allclose(abs(i), abs(j), rtol=1e-6)

        x = array(
            [
                -0.07827614,
                0.56637551,
                1.01320526,
                1.01536245,
                0.63548361,
                0.08560101,
                -0.46094955,
                -0.78065656,
                -0.8893556,
                -0.56514145,
                0.02325272,
                0.63660719,
                0.86291302,
                0.82953598,
                0.5706848,
                0.11655242,
                -0.6472655,
                -0.86178218,
                -0.96495057,
                -0.76098445,
                -0.18911517,
                0.59280646,
                1.00248693,
                0.89241423,
                0.52475111,
                -0.01620599,
                -0.60199278,
                -0.98279829,
                -1.12469771,
                -0.61355799,
                0.04321191,
                0.52784788,
                0.68508784,
                0.86015123,
                0.66825756,
                -0.0802846,
                -0.63626753,
                -0.93023345,
                -0.99129547,
                -0.46891033,
                0.04145813,
                0.71226518,
                1.01499246,
                0.94726778,
                0.63598143,
                -0.21920589,
                -0.48071702,
                -0.86041579,
                -0.9046141,
                -0.55714746,
                -0.10052384,
                0.69708969,
                1.02575789,
                1.16524031,
                0.49895282,
                -0.13068573,
                -0.45770419,
                -0.86155787,
                -0.9230734,
                -0.6590525,
                -0.05072955,
                0.52380317,
                1.02674335,
                0.87778499,
                0.4303284,
                -0.01855665,
                -0.62858193,
                -0.93954774,
                -0.94257301,
                -0.49692951,
                0.00699347,
                0.69049074,
                0.93906549,
                1.06339809,
                0.69337543,
                0.00252569,
                -0.57825881,
                -0.88460603,
                -0.99259672,
                -0.73535697,
                0.12064751,
                0.91159174,
                0.88966993,
                1.02159917,
                0.43479926,
                -0.06159005,
                -0.61782651,
                -0.95284676,
                -0.8218889,
                -0.52166419,
                0.021961,
                0.52268762,
                0.79428288,
                1.01642697,
                0.49060377,
                -0.02183994,
                -0.52743836,
                -0.99363909,
                -1.02963821,
                -0.64249996,
            ]
        )
        py_xc = zeros(2 * len(x) - 1, dtype=float64)
        numba_xc = py_xc.copy()
        N = 100
        py_autocorr_inner(x, py_xc, N)
        numba_autocorr_inner(x, numba_xc, N)
        for i, j in zip(py_xc, numba_xc):
            assert_allclose(i, j)

    def test_autocorr(self):
        """correctly compute autocorrelation"""
        s = [1, 1, 1, 1]
        X, periods = auto_corr(s, llim=-3, ulim=None)
        exp_X = array([1, 2, 3, 4, 3, 2, 1], dtype=float)
        assert_equal(X, exp_X)

        auto_x, auto_periods = auto_corr(self.sig, llim=2, ulim=50)
        max_idx = list(auto_x).index(max(auto_x))
        auto_p = auto_periods[max_idx]
        self.assertEqual(auto_p, self.p)

    def test_dft(self):
        """correctly compute discrete fourier transform"""
        dft_x, dft_periods = dft(self.sig)
        dft_x = abs(dft_x)
        max_idx = list(dft_x).index(max(dft_x))
        dft_p = dft_periods[max_idx]
        self.assertEqual(int(dft_p), self.p)

    def test_ipdft(self):
        """correctly compute integer discrete fourier transform"""
        s = [0, 1, 0, -1, 0, 1, 0, -1]
        X, periods = ipdft(s, llim=1, ulim=len(s))
        exp_X = abs(
            array(
                [
                    0,
                    0,
                    -1.5 + 0.866j,
                    -4j,
                    2.927 - 0.951j,
                    1.5 + 0.866j,
                    0.302 + 0.627j,
                    0,
                ]
            )
        )
        X = abs(X)
        assert_almost_equal(X, exp_X, decimal=4)

        ipdft_x, ipdft_periods = ipdft(self.sig, llim=2, ulim=50)
        ipdft_x = abs(ipdft_x)
        max_idx = list(ipdft_x).index(max(ipdft_x))
        ipdft_p = ipdft_periods[max_idx]
        self.assertEqual(ipdft_p, self.p)

    def test_goertzel(self):
        """goertzel and ipdft should be the same"""
        ipdft_pwr, ipdft_prd = ipdft(self.sig, llim=10, ulim=10)
        assert_allclose(goertzel(self.sig, 10), ipdft_pwr)

    def test_hybrid(self):
        """correctly compute hybrid statistic"""
        hybrid_x, hybrid_periods = hybrid(self.sig, llim=None, ulim=50)
        hybrid_x = abs(hybrid_x)
        max_idx = list(hybrid_x).index(max(hybrid_x))
        hybrid_p = hybrid_periods[max_idx]
        self.assertEqual(hybrid_p, self.p)

    def test_hybrid_returns_all(self):
        """correctly returns hybrid, ipdft and autocorr statistics"""
        ipdft_pwr, ipdft_prd = ipdft(self.sig, llim=2, ulim=50)
        auto_x, auto_periods = auto_corr(self.sig, llim=2, ulim=50)
        hybrid_x, hybrid_periods = hybrid(self.sig, llim=None, ulim=50)
        hybrid_ipdft_autocorr_stats, hybrid_periods = hybrid(
            self.sig, llim=None, ulim=50, return_all=True
        )
        assert_equal(hybrid_ipdft_autocorr_stats[0], hybrid_x)
        assert_equal(hybrid_ipdft_autocorr_stats[1], ipdft_pwr)
        assert_equal(hybrid_ipdft_autocorr_stats[2], auto_x)

        ipdft_pwr, ipdft_prd = ipdft(self.sig, llim=10, ulim=10)
        auto_x, auto_periods = auto_corr(self.sig, llim=10, ulim=10)
        hybrid_x, hybrid_periods = hybrid(self.sig, llim=10, ulim=10)
        hybrid_ipdft_autocorr_stats, hybrid_periods = hybrid(
            self.sig, llim=10, ulim=10, return_all=True
        )
        self.assertEqual(hybrid_ipdft_autocorr_stats[0], hybrid_x)
        self.assertEqual(hybrid_ipdft_autocorr_stats[1], ipdft_pwr)
        self.assertEqual(hybrid_ipdft_autocorr_stats[2], auto_x)


if __name__ == "__main__":
    main()
