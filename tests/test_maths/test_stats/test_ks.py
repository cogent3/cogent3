#!/usr/bin/env python
from unittest import TestCase, main

from cogent3.maths.stats.ks import (
    pkolmogorov1x,
    pkolmogorov2x,
    pkstwo,
    psmirnov2x,
)
from cogent3.maths.stats.test import ks_boot, ks_test


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose


class KSTests(TestCase):
    """Tests Kolmogorov-Smirnov."""

    def setUp(self):
        self.x1 = [
            0.09916191,
            0.29732882,
            0.41475044,
            0.68816838,
            0.20841367,
            0.46129887,
            0.22074544,
            0.06889561,
            0.88264852,
            0.87726406,
            0.76905072,
            0.86178033,
            0.42596777,
            0.59443782,
            0.68852176,
            0.66032130,
            0.72683791,
            0.02363118,
            0.82384762,
            0.32759965,
            0.69231127,
            0.50848596,
            0.67500888,
            0.84919139,
            0.70774136,
            0.97847465,
            0.59784714,
            0.82033663,
            0.45640039,
            0.13054766,
            0.01227875,
            0.21229238,
            0.37054602,
            0.80905622,
            0.26056527,
            0.01662457,
            0.76277188,
            0.76892495,
            0.39186350,
            0.61468789,
            0.83247770,
            0.69946238,
            0.80550609,
            0.22336814,
            0.62491296,
            0.03413056,
            0.74500251,
            0.36008309,
            0.19443889,
            0.06808133,
        ]
        self.x2 = [
            1.1177760,
            0.9984325,
            0.8113576,
            0.7247507,
            0.9473543,
            1.1192222,
            1.2577115,
            0.6168244,
            0.9616475,
            1.0677138,
            0.5106196,
            1.2334833,
            0.3750225,
            0.9788191,
            1.1366872,
            0.8212352,
            0.7665240,
            0.4409294,
            0.4447418,
            1.1381901,
            0.7299300,
            1.1307991,
            0.5356031,
            0.3193794,
            1.2476867,
            0.7909454,
            0.7781800,
            0.8438637,
            1.1814135,
            1.0117055,
            0.7433708,
            0.7917239,
            0.5080752,
            0.9014003,
            0.5960710,
            0.9646521,
            0.9263595,
            0.7969784,
            1.2847108,
            0.6393015,
            0.6828791,
            1.0817340,
            0.6586887,
            0.7314203,
            0.3998812,
            0.9988478,
            1.0225579,
            1.2721428,
            0.6465969,
            0.9133413,
        ]

    def test_pk1x(self):
        """1 sample 1-sided should match answers from R"""
        assert_allclose(pkolmogorov1x(0.06, 30), 0.2248113)

    def test_pk2x(self):
        """1 sample 2-sided should match answers from R"""
        assert_allclose(pkolmogorov2x(0.7199, 50), (1 - 6.661e-16), rtol=1e-5)
        assert_allclose(pkolmogorov2x(0.08, 30), 0.01754027, rtol=1e-5)
        assert_allclose(pkolmogorov2x(0.03, 300), 0.05753413, rtol=1e-5)

    def test_ps2x(self):
        """2 sample 2-sided smirnov should match answers from R"""
        assert_allclose(psmirnov2x(0.48, 20, 50), 0.9982277)
        assert_allclose(psmirnov2x(0.28, 20, 50), 0.8161612)
        assert_allclose(psmirnov2x(0.28, 50, 20), 0.8161612)

    def tes_pk2x(self):
        """2 sample  2-sided kolmogorov should match answers from R"""
        assert_allclose(pkolmogorov1x(0.058, 50), 0.007530237)
        assert_allclose(pkolmogorov1x(0.018, 50), 4.887356e-26)
        assert_allclose(pkolmogorov1x(0.018, 5000), 0.922618)

    def test_pkstwo(self):
        """kolmogorov asymptotic should match answers from R"""
        assert_allclose(pkstwo(2.3), [1 - 5.084e-05], rtol=1e-5)

    def test_ks2x(self):
        """KS two-sample, 2-sided should match answers from R"""
        D, Pval = ks_test(self.x1, self.x2)
        assert_allclose((D, Pval), (0.46, 3.801e-05), rtol=1e-4)
        D, Pval = ks_test(self.x1, self.x2, exact=False)
        assert_allclose((D, Pval), (0.46, 5.084e-05), rtol=1e-4)
        D, Pval = ks_test(self.x1, self.x2[:20])
        assert_allclose((D, Pval), (0.53, 0.0003576), rtol=1e-4)
        D, Pval = ks_test(self.x2[:20], self.x1)
        assert_allclose((D, Pval), (0.53, 0.0003576), rtol=1e-4)
        D, Pval = ks_test(self.x1[:20], self.x2)
        assert_allclose((D, Pval), (0.48, 0.001772), rtol=1e-3)
        D, Pval = ks_test(self.x1, self.x2, alt="greater")
        assert_allclose((D, Pval), (0.46, 2.542e-05), rtol=1e-4)
        D, Pval = ks_test(self.x1, self.x2, alt="g")
        assert_allclose((D, Pval), (0.46, 2.542e-05), rtol=1e-4)
        D, Pval = ks_test(self.x1, self.x2, alt="less")
        assert_allclose((D, Pval), (6.9388939039072284e-18, 1.0), rtol=1e-4)
        D, Pval = ks_test(self.x2, self.x1, alt="l")
        assert_allclose((D, Pval), (0.46, 2.542e-05), rtol=1e-4)

    def test_ks_boot(self):
        """excercising the bootstrapped version of KS"""
        D, Pval = ks_boot(self.x1[:10], self.x2[:10], num_reps=10)


if __name__ == "__main__":
    main()
