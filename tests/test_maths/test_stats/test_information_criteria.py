from unittest import TestCase

from cogent3.maths.stats.information_criteria import aic, bic
from numpy.testing import assert_allclose


class InformationCriteria(TestCase):
    """Tests calculation of AIC and BIC measures."""

    def test_aic(self):
        """correctly compute AIC from Burnham & Anderson 2002, p102"""
        assert_allclose(aic(-9.7039, 4), 27.4078)

    def test_aic_corrected(self):
        """correctly compute AIC corrected for small sample size"""
        # from Burnham & Anderson 2002, p102
        assert_allclose(aic(-9.7039, 4, sample_size=13), 32.4078)

    def test_bic(self):
        """correctly compute BIC"""
        # against hand calculated
        assert_allclose(bic(-9.7039, 4, 13), 29.6675974298)
