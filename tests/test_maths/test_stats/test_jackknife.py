from unittest import TestCase, main

import numpy as np

from cogent3.maths.stats.jackknife import JackknifeStats


__author__ = "Anuj Pahwa, Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Anuj Pahwa", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose


def pmcc(data, axis=1):
    """Compute the Product-moment correlation coefficient.
    Expression 15.3 from Biometry by Sokal/Rohlf
    This code implementation is on the proviso that the data that is provided
    is two dimensional: [[Y1], [Y2]] (trying to determine the correlation
    coefficient between data sets Y1 and Y2"""

    if axis == 0:
        data = data.transpose()
        axis = 1

    other_axis = 0
    mean = data.mean(axis=axis)
    data_less_mean = np.array([data[0] - mean[0], data[1] - mean[1]])
    sum_squares = np.sum(np.square(data_less_mean), axis=axis)
    sum_products = np.sum(np.prod(data_less_mean, axis=other_axis))
    pmcc = np.divide(sum_products, np.sqrt(np.prod(sum_squares)))
    z_trans = np.arctanh(pmcc)
    return z_trans


# test data from Box 15.2; Biometry by Sokal/Rohlf
data = np.array(
    [
        [159, 179, 100, 45, 384, 230, 100, 320, 80, 220, 320, 210],
        [
            14.40,
            15.20,
            11.30,
            2.50,
            22.70,
            14.90,
            1.41,
            15.81,
            4.19,
            15.39,
            17.25,
            9.52,
        ],
    ]
)

# factory function generator for the statistical function of interest


def stat_maker(func, data, axis):
    def calc_stat(coords):
        subset_data = data.take(coords, axis)
        return func(subset_data, axis)

    return calc_stat


# function to compute mean of a np array


def mean(data, axis):
    return data.mean(axis=axis)


class JackknifeTests(TestCase):
    def test_proper_initialise(self):
        """jackknife should initialise correctly"""
        # Scalar
        pmcc_stat = stat_maker(pmcc, data, 1)
        test_knife = JackknifeStats(data.shape[1], pmcc_stat)
        self.assertEqual(test_knife.n, data.shape[1])
        self.assertEqual(test_knife._jackknifed_stat, None)

        # Vector
        mean_stat = stat_maker(mean, data, 1)
        test_knife = JackknifeStats(data.shape[1], mean_stat)
        self.assertEqual(test_knife.n, data.shape[1])
        self.assertEqual(test_knife._jackknifed_stat, None)

    def test_jackknife_stats(self):
        """jackknife results should match Sokal & Rolf example"""
        # Scalar
        pmcc_stat = stat_maker(pmcc, data, 1)
        test_knife = JackknifeStats(data.shape[1], pmcc_stat)
        assert_allclose(test_knife.jackknifed_stat, 1.2905845)
        assert_allclose(test_knife.standard_error, 0.2884490)
        self.assertTrue(test_knife._jackknifed_stat != None)

        # Vector
        mean_stat = stat_maker(mean, data, 1)
        test_knife = JackknifeStats(data.shape[1], mean_stat)
        expected_jk_stat = data.mean(axis=1)
        got_jk_stat = test_knife.jackknifed_stat
        expected_standard_err = [30.69509346, 1.87179671]
        got_standard_err = test_knife.standard_error

        for index in [0, 1]:
            assert_allclose(got_jk_stat[index], expected_jk_stat[index])
            assert_allclose(got_standard_err[index], expected_standard_err[index])

    def test_tables(self):
        """jackknife should work for calculators return scalars or vectors"""
        # Scalar
        pmcc_stat = stat_maker(pmcc, data, 1)
        test_knife = JackknifeStats(data.shape[1], pmcc_stat)
        expected_subsample_stats = [
            1.4151,
            1.3946,
            1.4314,
            1.1889,
            1.1323,
            1.3083,
            1.3561,
            1.3453,
            1.2412,
            1.3216,
            1.2871,
            1.3664,
        ]
        expected_pseudovalues = [
            0.1968,
            0.4224,
            0.0176,
            2.6852,
            3.3084,
            1.3718,
            0.8461,
            0.9650,
            2.1103,
            1.2253,
            1.6049,
            0.7333,
        ]
        test_knife.jackknife()
        got_subsample_stats = test_knife._subset_statistics
        got_pseudovalues = test_knife._pseudovalues
        for index in range(data.shape[1]):
            np.testing.assert_almost_equal(
                got_subsample_stats[index], expected_subsample_stats[index], 4
            )
            np.testing.assert_approx_equal(
                got_pseudovalues[index], expected_pseudovalues[index], 4
            )

        # Vector
        mean_stat = stat_maker(mean, data, 1)
        test_knife = JackknifeStats(data.shape[1], mean_stat)

        test_knife.jackknife()
        expected_pseudovalues = data.transpose()
        expected_subsample_stats = [
            [198.9091, 11.8336],
            [197.0909, 11.7609],
            [204.2727, 12.1155],
            [209.2727, 12.9155],
            [178.4545, 11.0791],
            [192.4545, 11.7882],
            [204.2727, 13.0145],
            [184.2727, 11.7055],
            [206.0909, 12.7618],
            [193.3636, 11.7436],
            [184.2727, 11.5745],
            [194.2727, 12.2773],
        ]
        got_subsample_stats = test_knife._subset_statistics
        got_pseudovalues = test_knife._pseudovalues

        for index1 in range(data.shape[1]):
            for index2 in range(data.shape[0]):
                np.testing.assert_almost_equal(
                    got_subsample_stats[index1][index2],
                    expected_subsample_stats[index1][index2],
                    4,
                )
                np.testing.assert_almost_equal(
                    got_pseudovalues[index1][index2],
                    expected_pseudovalues[index1][index2],
                    4,
                )

    def test_tabular_properties(self):
        """constructs tabular properties"""
        pmcc_stat = stat_maker(pmcc, data, 1)
        test_knife = JackknifeStats(data.shape[1], pmcc_stat)
        ss = test_knife.sub_sample_stats
        self.assertEqual(ss.shape, (12, 2))
        ss = test_knife.sample_stat
        pvs = test_knife.pseudovalues
        self.assertEqual(pvs.shape, (12, 2))
        ss = test_knife.summary_stats
        self.assertEqual(ss.shape, (1, 3))


if __name__ == "__main__":
    main()
