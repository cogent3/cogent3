#!/usr/bin/env python
#file test_alpha_diversity.py
from numpy import array, log, sqrt, exp
from math import e
from cogent.util.unit_test import TestCase, main
from cogent.maths.stats.alpha_diversity import expand_counts, counts, observed_species, singles, \
    doubles, osd, margalef, menhinick, dominance, simpson, reciprocal_simpson,\
    shannon, equitability, berger_parker_d, mcintosh_d, brillouin_d, \
    strong, kempton_taylor_q, fisher_alpha, \
    mcintosh_e, heip_e, simpson_e, robbins, robbins_confidence, \
    chao1_uncorrected, chao1_bias_corrected, chao1, chao1_var, \
    chao1_confidence, ACE

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class diversity_tests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Set up shared variables"""
        self.TestData = array([0,1,1,4,2,5,2,4,1,2])
        self.NoSingles = array([0,2,2,4,5,0,0,0,0,0])
        self.NoDoubles = array([0,1,1,4,5,0,0,0,0,0])

    def test_expand_counts(self):
        """expand_counts should return correct expanded array"""
        c = array([2,0,1,2])
        self.assertEqual(expand_counts(c), array([0,0,2,3,3]))

    def test_counts(self):
        """counts should return correct array"""
        c = array([5,0,1,1,5,5])
        obs = counts(c)
        exp = array([1,2,0,0,0,3])
        self.assertEqual(obs, exp)
        d = array([2,2,1,0])
        obs = counts(d, obs)
        exp = array([2,3,2,0,0,3])
        self.assertEqual(obs, exp)

    def test_singles(self):
        """singles should return correct # of singles"""
        self.assertEqual(singles(self.TestData), 3)
        self.assertEqual(singles(array([0,3,4])), 0)
        self.assertEqual(singles(array([1])), 1)

    def test_doubles(self):
        """doubles should return correct # of doubles"""
        self.assertEqual(doubles(self.TestData), 3)
        self.assertEqual(doubles(array([0,3,4])), 0)
        self.assertEqual(doubles(array([2])), 1)

    def test_osd(self):
        """osd should return correct # of observeds, singles, doubles"""
        self.assertEqual(osd(self.TestData), (9,3,3))

    def test_margalef(self):
        """margalef should match hand-calculated values"""
        self.assertEqual(margalef(self.TestData), 8/log(22))

    def test_menhinick(self):
        """menhinick should match hand-calculated values"""
        self.assertEqual(menhinick(self.TestData), 9/sqrt(22))

    def test_dominance(self):
        """dominance should match hand-calculated values"""
        c = array([1,0,2,5,2])
        self.assertFloatEqual(dominance(c), .34)
        d = array([5])
        self.assertEqual(dominance(d), 1)

    def test_simpson(self):
        """simpson should match hand-calculated values"""
        c = array([1,0,2,5,2])
        self.assertFloatEqual(simpson(c), .66)
        d = array([5])
        self.assertFloatEqual(simpson(d), 0)

    def test_reciprocal_simpson(self):
        """reciprocal_simpson should match hand-calculated results"""
        c = array([1,0,2,5,2])
        self.assertFloatEqual(reciprocal_simpson(c), 1/.66)

    def test_shannon(self):
        """shannon should match hand-calculated values"""
        c = array([5])
        self.assertFloatEqual(shannon(c), 0)
        c = array([5,5])
        self.assertFloatEqual(shannon(c), 1)
        c = array([1,1,1,1,0])
        self.assertEqual(shannon(c), 2)

    def test_equitability(self):
        """equitability should match hand-calculated values"""
        c = array([5])
        self.assertFloatEqual(equitability(c), 0)
        c = array([5,5])
        self.assertFloatEqual(equitability(c), 1)
        c = array([1,1,1,1,0])
        self.assertEqual(equitability(c), 1)

    def test_berger_parker_d(self):
        """berger-parker_d should match hand-calculated values"""
        c = array([5])
        self.assertFloatEqual(berger_parker_d(c), 1)
        c = array([5,5])
        self.assertFloatEqual(berger_parker_d(c), 0.5)
        c = array([1,1,1,1,0])
        self.assertEqual(berger_parker_d(c), 0.25)

    def test_mcintosh_d(self):
        """mcintosh_d should match hand-calculated values"""
        c = array([1,2,3])
        self.assertFloatEqual(mcintosh_d(c), 0.636061424871458)

    def test_brillouin_d(self):
        """brillouin_d should match hand-calculated values"""
        c = array([1,2,3,1])
        self.assertFloatEqual(brillouin_d(c), 0.86289353018248782)

    def test_strong(self):
        """strong's dominance index should match hand-calculated values"""
        c = array([1,2,3,1])
        self.assertFloatEqual(strong(c), 0.214285714)

    def test_kempton_taylor_q(self):
        """kempton_taylor_q should approximate Magurran 1998 calculation p143"""
        c = array([2,3,3,3,3,3,4,4,4,6,6,7,7,9,9,11,14,15,15,20,29,33,34,
            36,37,53,57,138,146,170])
        self.assertFloatEqual(kempton_taylor_q(c), 14/log(34/4))

    def test_fisher_alpha(self):
        """fisher alpha should match hand-calculated value."""
        c = array([4,3,4,0,1,0,2])
        obs = fisher_alpha(c)
        self.assertFloatEqual(obs, 2.7823795367398798)

    def test_mcintosh_e(self):
        """mcintosh e should match hand-calculated value."""
        c = array([1,2,3,1])
        num = sqrt(15)
        den = sqrt(19)
        exp = num/den
        self.assertEqual(mcintosh_e(c), exp)

    def test_heip_e(self):
        """heip e should match hand-calculated value"""
        c = array([1,2,3,1])
        h = shannon(c, base=e)
        expected = exp(h-1)/3
        self.assertEqual(heip_e(c), expected)

    def test_simpson_e(self):
        """simpson e should match hand-calculated value"""
        c = array([1,2,3,1])
        s = simpson(c)
        self.assertEqual((1/s)/4, simpson_e(c))

    def test_robbins(self):
        """robbins metric should match hand-calculated value"""
        c = array([1,2,3,0,1])
        r = robbins(c)
        self.assertEqual(r,2./7) 

    def test_robbins_confidence(self):
        """robbins CI should match hand-calculated value"""
        c = array([1,2,3,0,1])
        r = robbins_confidence(c, 0.05)
        n = 7
        s = 2
        k = sqrt(8/0.05)
        self.assertEqual(r, ((s-k)/(n+1), (s+k)/(n+1))) 


    def test_observed_species(self):
        """observed_species should return # observed species"""
        c = array([4,3,4,0,1,0,2])
        obs = observed_species(c)
        exp = 5
        self.assertEqual(obs, exp)
        c = array([0,0,0])
        obs = observed_species(c)
        exp = 0
        self.assertEqual(obs, exp)
        self.assertEqual(observed_species(self.TestData), 9)

    def test_chao1_bias_corrected(self):
        """chao1_bias_corrected should return same result as EstimateS"""
        obs = chao1_bias_corrected(*osd(self.TestData))
        self.assertEqual(obs, 9.75)

    def test_chao1_uncorrected(self):
        """chao1_uncorrected should return same result as EstimateS"""
        obs = chao1_uncorrected(*osd(self.TestData))
        self.assertEqual(obs, 10.5)

    def test_chao1(self):
        """chao1 should use right decision rules"""
        self.assertEqual(chao1(self.TestData), 9.75)
        self.assertEqual(chao1(self.TestData,bias_corrected=False),10.5)
        self.assertEqual(chao1(self.NoSingles), 4)
        self.assertEqual(chao1(self.NoSingles,bias_corrected=False),4)
        self.assertEqual(chao1(self.NoDoubles), 5)
        self.assertEqual(chao1(self.NoDoubles,bias_corrected=False),5)

    def test_chao1_var(self):
        """chao1_var should match observed results from EstimateS"""
        #NOTE: EstimateS reports sd, not var, and rounds to 2 dp
        self.assertFloatEqual(chao1_var(self.TestData), 1.42**2, eps=0.01)
        self.assertFloatEqual(chao1_var(self.TestData,bias_corrected=False),\
            2.29**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoSingles), 0.39**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoSingles, \
            bias_corrected=False), 0.39**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoDoubles), 2.17**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoDoubles, \
            bias_corrected=False), 2.17**2, eps=0.01)

    def test_chao1_confidence(self):
        """chao1_confidence should match observed results from EstimateS""" 
        #NOTE: EstimateS rounds to 2 dp
        self.assertFloatEqual(chao1_confidence(self.TestData), (9.07,17.45), \
            eps=0.01)
        self.assertFloatEqual(chao1_confidence(self.TestData, \
            bias_corrected=False), (9.17,21.89), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoSingles),\
            (4, 4.95), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoSingles, \
            bias_corrected=False), (4,4.95), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoDoubles), \
            (4.08,17.27), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoDoubles, \
            bias_corrected=False), (4.08,17.27), eps=0.01)
    
    def test_ACE(self):
        """ACE should match values calculated by hand""" 
        self.assertFloatEqual(ACE(array([2,0])), 1.0, eps=0.001)
        # next: just returns the number of species when all are abundant
        self.assertFloatEqual(ACE(array([12,0,9])), 2.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,2,8])), 3.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,2,1])), 4.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,1,2,1])), 7.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,3,2,1])), 4.6, eps=0.001)
        self.assertFloatEqual(ACE(array([12,3,6,1,10])), 5.62749672, eps=0.001)



if __name__ == '__main__':
    main()
