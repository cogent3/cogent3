#!/usr/bin/env python
"""Unit tests for special functions used in statistics.
"""
import math

from unittest import TestCase, main

from cogent3.maths.stats.special import (
    combinations,
    combinations_exact,
    igami,
    incbi,
    ln_binomial,
    ln_combinations,
    ln_permutations,
    log1p,
    log_one_minus,
    ndtri,
    one_minus_exp,
    permutations,
    permutations_exact,
)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Sandra Smit"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose


class SpecialTests(TestCase):
    """Tests miscellaneous functions."""

    def test_permutations(self):
        """permutations should return expected results"""
        self.assertEqual(permutations(1, 1), 1)
        self.assertEqual(permutations(2, 1), 2)
        self.assertEqual(permutations(3, 1), 3)
        self.assertEqual(permutations(4, 1), 4)
        self.assertEqual(permutations(4, 2), 12)
        self.assertEqual(permutations(4, 3), 24)
        self.assertEqual(permutations(4, 4), 24)
        assert_allclose(permutations(300, 100), 3.8807387193009318e239)

    def test_permutations_errors(self):
        """permutations should raise errors on invalid input"""
        self.assertRaises(IndexError, permutations, 10, 50)
        self.assertRaises(IndexError, permutations, -1, 50)
        self.assertRaises(IndexError, permutations, 10, -5)

    def test_permutations_float(self):
        """permutations should use gamma function when floats as input"""
        assert_allclose(permutations(1.0, 1), 1)
        assert_allclose(permutations(2, 1.0), 2)
        assert_allclose(permutations(3.0, 1.0), 3)
        assert_allclose(permutations(4.0, 1), 4)
        assert_allclose(permutations(4.0, 2.0), 12)
        assert_allclose(permutations(4.0, 3.0), 24)
        assert_allclose(permutations(4, 4.0), 24)
        assert_allclose(permutations(300, 100), 3.8807387193009318e239)

    def test_permutations_range(self):
        """permutations should increase gradually with increasing k"""
        start = 5  # permuations(10,5) = 30240
        end = 6  # permutations(10,6) = 151200
        step = 0.1
        lower_lim = 30240
        upper_lim = 151200
        previous_value = 30239.9999
        while start <= end:
            obs = permutations(10, start)
            assert lower_lim <= obs <= upper_lim
            assert obs > previous_value
            previous_value = obs
            start += step

    def test_permutations_exact(self):
        """permutations_exact should return expected results"""
        assert_allclose(permutations_exact(1, 1), 1)
        assert_allclose(permutations_exact(2, 1), 2)
        assert_allclose(permutations_exact(3, 1), 3)
        assert_allclose(permutations_exact(4, 1), 4)
        assert_allclose(permutations_exact(4, 2), 12)
        assert_allclose(permutations_exact(4, 3), 24)
        assert_allclose(permutations_exact(4, 4), 24)
        assert_allclose(permutations_exact(300, 100) / 3.8807387193009318e239, 1.0)

    def test_ln_permutations(self):
        """ln_permutations should return expected results"""
        assert_allclose(ln_permutations(1, 1), math.log(1))
        assert_allclose(ln_permutations(2, 1), math.log(2))
        assert_allclose(ln_permutations(3, 1.0), math.log(3))
        assert_allclose(ln_permutations(4, 1), math.log(4))
        assert_allclose(ln_permutations(4.0, 2), math.log(12))
        assert_allclose(ln_permutations(4, 3.0), math.log(24))
        assert_allclose(ln_permutations(4, 4), math.log(24))
        assert_allclose(ln_permutations(300.0, 100), math.log(3.8807387193009318e239))

    def test_combinations(self):
        """combinations should return expected results when int as input"""
        self.assertEqual(combinations(1, 1), 1)
        self.assertEqual(combinations(2, 1), 2)
        self.assertEqual(combinations(3, 1), 3)
        self.assertEqual(combinations(4, 1), 4)
        self.assertEqual(combinations(4, 2), 6)
        self.assertEqual(combinations(4, 3), 4)
        self.assertEqual(combinations(4, 4), 1)
        self.assertEqual(combinations(20, 4), 19 * 17 * 15)
        assert_allclose(combinations(300, 100), 4.1582514632578812e81)

    def test_combinations_errors(self):
        """combinations should raise errors on invalid input"""
        self.assertRaises(IndexError, combinations, 10, 50)
        self.assertRaises(IndexError, combinations, -1, 50)
        self.assertRaises(IndexError, combinations, 10, -5)

    def test_combinations_float(self):
        """combinations should use gamma function when floats as input"""
        assert_allclose(combinations(1.0, 1.0), 1)
        assert_allclose(combinations(2.0, 1.0), 2)
        assert_allclose(combinations(3.0, 1.0), 3)
        assert_allclose(combinations(4.0, 1.0), 4)
        assert_allclose(combinations(4.0, 2), 6)
        assert_allclose(combinations(4, 3.0), 4)
        assert_allclose(combinations(4.0, 4.0), 1)
        assert_allclose(combinations(20.0, 4.0), 19 * 17 * 15)
        assert_allclose(combinations(300, 100.0), 4.1582514632578812e81)

    def test_combinations_range(self):
        """combinations should decrease gradually with increasing k"""
        start = 5  # combinations(10,5) = 252
        end = 6  # combinations(10,6) = 210
        step = 0.1
        lower_lim = 210
        upper_lim = 252
        previous_value = 252.00001
        while start <= end:
            obs = combinations(10, start)
            assert lower_lim <= obs <= upper_lim
            assert obs < previous_value
            previous_value = obs
            start += step

    def test_combinations_exact(self):
        """combinations_exact should return expected results"""
        self.assertEqual(combinations_exact(1, 1), 1)
        self.assertEqual(combinations_exact(2, 1), 2)
        self.assertEqual(combinations_exact(3, 1), 3)
        self.assertEqual(combinations_exact(4, 1), 4)
        self.assertEqual(combinations_exact(4, 2), 6)
        self.assertEqual(combinations_exact(4, 3), 4)
        self.assertEqual(combinations_exact(4, 4), 1)
        self.assertEqual(combinations_exact(20, 4), 19 * 17 * 15)
        assert_allclose(combinations_exact(300, 100), 4.1582514632578812e81)

    def test_ln_combinations(self):
        """ln_combinations should return expected results"""
        assert_allclose(ln_combinations(1, 1), math.log(1))
        assert_allclose(ln_combinations(2, 1), math.log(2))
        assert_allclose(ln_combinations(3, 1), math.log(3))
        assert_allclose(ln_combinations(4.0, 1), math.log(4))
        assert_allclose(ln_combinations(4, 2.0), math.log(6))
        assert_allclose(ln_combinations(4, 3), math.log(4))
        assert_allclose(ln_combinations(4, 4.0), math.log(1))
        assert_allclose(ln_combinations(20, 4), math.log(19 * 17 * 15))
        assert_allclose(ln_combinations(300, 100), math.log(4.1582514632578812e81))

    def test_ln_binomial_integer(self):
        """ln_binomial should match R results for integer values"""
        expected = {
            (10, 60, 0.1): -3.247883,
            (1, 1, 0.5): math.log(0.5),
            (1, 1, 0.0000001): math.log(1e-07),
            (1, 1, 0.9999999): math.log(0.9999999),
            (3, 5, 0.75): math.log(0.2636719),
            (0, 60, 0.5): math.log(8.673617e-19),
            (129, 130, 0.5): math.log(9.550892e-38),
            (299, 300, 0.099): math.log(1.338965e-298),
            (9, 27, 0.0003): math.log(9.175389e-26),
            (1032, 2050, 0.5): math.log(0.01679804),
        }
        for (key, value) in list(expected.items()):
            assert_allclose(ln_binomial(*key), value, 1e-4)

    def test_ln_binomial_floats(self):
        """Binomial exact should match values from R for integer successes"""
        expected = {
            (18.3, 100, 0.2): (math.log(0.09089812), math.log(0.09807429)),
            (2.7, 1050, 0.006): (math.log(0.03615498), math.log(0.07623827)),
            (2.7, 1050, 0.06): (math.log(1.365299e-25), math.log(3.044327e-24)),
            (2, 100.5, 0.6): (math.log(7.303533e-37), math.log(1.789727e-36)),
            (0.2, 60, 0.5): (math.log(8.673617e-19), math.log(5.20417e-17)),
            (0.5, 5, 0.3): (math.log(0.16807), math.log(0.36015)),
            (10, 100.5, 0.5): (math.log(7.578011e-18), math.log(1.365543e-17)),
        }

        for (key, value) in list(expected.items()):
            min_val, max_val = value
            assert min_val < ln_binomial(*key) < max_val
            # self.assert_allclose(binomial_exact(*key), value, 1e-4)

    def test_ln_binomial_range(self):
        """ln_binomial should increase in a monotonically increasing region."""
        start = 0
        end = 1
        step = 0.1
        lower_lim = -1.783375 - 1e-4
        upper_lim = -1.021235 + 1e-4
        previous_value = -1.784
        while start <= end:
            obs = ln_binomial(start, 5, 0.3)
            assert lower_lim <= obs <= upper_lim
            assert obs > previous_value
            previous_value = obs
            start += step

    def test_log_one_minus_large(self):
        """log_one_minus_x should return math.log(1-x) if x is large"""
        assert_allclose(log_one_minus(0.2), math.log(1 - 0.2))

    def test_log_one_minus_small(self):
        """log_one_minus_x should return -x if x is small"""
        assert_allclose(log_one_minus(1e-30), -1e-30, rtol=1e-30, atol=1e-30)

    def test_one_minus_exp_large(self):
        """one_minus_exp_x should return 1 - math.exp(x) if x is large"""
        assert_allclose(one_minus_exp(0.2), 1 - (math.exp(0.2)))

    def test_one_minus_exp_small(self):
        """one_minus_exp_x should return -x if x is small"""
        assert_allclose(one_minus_exp(1e-30), -1e-30)

    def test_log1p(self):
        """log1p should give same results as cephes"""
        p_s = [
            1e-10,
            1e-5,
            0.1,
            0.8,
            0.9,
            0.95,
            0.999,
            0.9999999,
            1,
            1.000000001,
            1.01,
            2,
        ]
        exp = [
            9.9999999995e-11,
            9.99995000033e-06,
            0.0953101798043,
            0.587786664902,
            0.641853886172,
            0.667829372576,
            0.692647055518,
            0.69314713056,
            0.69314718056,
            0.69314718106,
            0.698134722071,
            1.09861228867,
        ]
        for p, e in zip(p_s, exp):
            assert_allclose(log1p(p), e)

    def test_igami(self):
        """igami should give same result as cephes implementation"""
        a_vals = [1e-10, 1e-5, 0.5, 1, 10, 200]
        y_vals = list(range(0, 10, 2))
        obs = [igami(a, y / 10.0) for a in a_vals for y in y_vals]
        exp = [
            1.79769313486e308,
            0.0,
            0.0,
            0.0,
            0.0,
            1.79769313486e308,
            0.0,
            0.0,
            0.0,
            0.0,
            1.79769313486e308,
            0.821187207575,
            0.3541631504,
            0.137497948864,
            0.0320923773337,
            1.79769313486e308,
            1.60943791243,
            0.916290731874,
            0.510825623766,
            0.223143551314,
            1.79769313486e308,
            12.5187528198,
            10.4756841889,
            8.9044147366,
            7.28921960854,
            1.79769313486e308,
            211.794753362,
            203.267574402,
            196.108740945,
            188.010915412,
        ]
        for o, e in zip(obs, exp):
            assert_allclose(o, e)

    def test_ndtri(self):
        """ndtri should give same result as implementation in cephes"""
        exp = [
            -1.79769313486e308,
            -2.32634787404,
            -2.05374891063,
            -1.88079360815,
            -1.75068607125,
            -1.64485362695,
            -1.5547735946,
            -1.47579102818,
            -1.40507156031,
            -1.34075503369,
            -1.28155156554,
            -1.22652812004,
            -1.17498679207,
            -1.12639112904,
            -1.08031934081,
            -1.03643338949,
            -0.99445788321,
            -0.954165253146,
            -0.915365087843,
            -0.877896295051,
            -0.841621233573,
            -0.806421247018,
            -0.772193214189,
            -0.738846849185,
            -0.70630256284,
            -0.674489750196,
            -0.643345405393,
            -0.612812991017,
            -0.582841507271,
            -0.553384719556,
            -0.524400512708,
            -0.495850347347,
            -0.467698799115,
            -0.439913165673,
            -0.412463129441,
            -0.385320466408,
            -0.358458793251,
            -0.331853346437,
            -0.305480788099,
            -0.279319034447,
            -0.253347103136,
            -0.227544976641,
            -0.201893479142,
            -0.176374164781,
            -0.150969215497,
            -0.125661346855,
            -0.100433720511,
            -0.0752698620998,
            -0.0501535834647,
            -0.0250689082587,
            0.0,
            0.0250689082587,
            0.0501535834647,
            0.0752698620998,
            0.100433720511,
            0.125661346855,
            0.150969215497,
            0.176374164781,
            0.201893479142,
            0.227544976641,
            0.253347103136,
            0.279319034447,
            0.305480788099,
            0.331853346437,
            0.358458793251,
            0.385320466408,
            0.412463129441,
            0.439913165673,
            0.467698799115,
            0.495850347347,
            0.524400512708,
            0.553384719556,
            0.582841507271,
            0.612812991017,
            0.643345405393,
            0.674489750196,
            0.70630256284,
            0.738846849185,
            0.772193214189,
            0.806421247018,
            0.841621233573,
            0.877896295051,
            0.915365087843,
            0.954165253146,
            0.99445788321,
            1.03643338949,
            1.08031934081,
            1.12639112904,
            1.17498679207,
            1.22652812004,
            1.28155156554,
            1.34075503369,
            1.40507156031,
            1.47579102818,
            1.5547735946,
            1.64485362695,
            1.75068607125,
            1.88079360815,
            2.05374891063,
            2.32634787404,
        ]
        obs = [ndtri(i / 100.0) for i in range(100)]
        assert_allclose(obs, exp)

    def test_incbi(self):
        """incbi results should match cephes libraries"""
        aa_range = [0.1, 0.2, 0.5, 1, 2, 5]
        bb_range = aa_range
        yy_range = [0.1, 0.2, 0.5, 0.9]
        exp = [
            8.86928001193e-08,
            9.08146855855e-05,
            0.5,
            0.999999911307,
            4.39887474012e-09,
            4.50443299194e-06,
            0.0416524955556,
            0.997881005025,
            3.46456275553e-10,
            3.54771169012e-07,
            0.00337816430373,
            0.732777808689,
            1e-10,
            1.024e-07,
            0.0009765625,
            0.3486784401,
            3.85543289443e-11,
            3.94796342545e-08,
            0.000376636057552,
            0.154915841005,
            1.33210087225e-11,
            1.36407136078e-08,
            0.000130149552409,
            0.056682323296,
            0.00211899497509,
            0.0646097657259,
            0.958347504444,
            0.999999995601,
            0.000247764691908,
            0.00788804962659,
            0.5,
            0.999752235308,
            3.09753032747e-05,
            0.000990813218262,
            0.092990311753,
            0.906714634947,
            1e-05,
            0.00032,
            0.03125,
            0.59049,
            4.01878917904e-06,
            0.000128614607219,
            0.0126923538971,
            0.309157452156,
            1.41593162013e-06,
            4.5316442592e-05,
            0.00449136140034,
            0.122896698096,
            0.267222191311,
            0.684264602461,
            0.996621835696,
            0.999999999654,
            0.0932853650529,
            0.321847764104,
            0.907009688247,
            0.999969024697,
            0.0244717418524,
            0.0954915028125,
            0.5,
            0.975528258148,
            0.01,
            0.04,
            0.25,
            0.81,
            0.00445768188762,
            0.0179929616503,
            0.120614758428,
            0.531877433474,
            0.00165851285512,
            0.00672409501831,
            0.046687245337,
            0.247272226803,
            0.6513215599,
            0.8926258176,
            0.9990234375,
            0.9999999999,
            0.40951,
            0.67232,
            0.96875,
            0.99999,
            0.19,
            0.36,
            0.75,
            0.99,
            0.1,
            0.2,
            0.5,
            0.9,
            0.0513167019495,
            0.105572809,
            0.292893218813,
            0.683772233983,
            0.020851637639,
            0.04364750021,
            0.129449436704,
            0.36904265552,
            0.845084158995,
            0.956946913164,
            0.999623363942,
            0.999999999961,
            0.690842547844,
            0.850620771098,
            0.987307646103,
            0.999995981211,
            0.468122566526,
            0.629849697132,
            0.879385241572,
            0.995542318112,
            0.316227766017,
            0.4472135955,
            0.707106781187,
            0.948683298051,
            0.195800105659,
            0.287140725417,
            0.5,
            0.804199894341,
            0.0925952589131,
            0.13988068827,
            0.264449983296,
            0.510316306551,
            0.943317676704,
            0.984896695084,
            0.999869850448,
            0.999999999987,
            0.877103301904,
            0.944441767096,
            0.9955086386,
            0.999998584068,
            0.752727773197,
            0.841546267738,
            0.953312754663,
            0.998341487145,
            0.63095734448,
            0.724779663678,
            0.870550563296,
            0.979148362361,
            0.489683693449,
            0.577552475154,
            0.735550016704,
            0.907404741087,
            0.300968763593,
            0.366086516536,
            0.5,
            0.699031236407,
        ]
        i = 0
        for a in aa_range:
            for b in bb_range:
                for y in yy_range:
                    result = incbi(a, b, y)
                    e = exp[i]
                    assert_allclose(e, result)
                    i += 1
        # specific cases that failed elsewhere
        assert_allclose(incbi(999, 2, 1e-10), 0.97399698104554944)


# execute tests if called from command line
if __name__ == "__main__":
    main()
