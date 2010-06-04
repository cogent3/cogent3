#!/usr/bin/env python
#file test_fisher_exact.py

__author__ = "AdeeB NooR"
__copyright__ = "Copyright 2010, The GenomeDB project"
__credits__ = ["Adeeb Noor", "Rob Knight" , "Jesse Zaneveld"]
__license__ = "GPL"
_version__ = "1.0-dev"
__maintainer__ = "Adeeb Noor"
__email__ = "adeeb.noor@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from fisher_exact_test import fact,n_choose_k,hypergeometric_pmf,hypergeometric_cdf,hypergeometric_sf,fisher_exact_test

class FisherExactTests(TestCase):
    """ Test code for fisher_exact_test.py """

    def test_fact(self): 
        """ fact should calculate factorial """
        x=10
        obs = fact(x)
        exp = 3628800 
    	self.assertEqual(obs,exp)

    def test_n_choose_k(self): 
        """ n_choose_k should calculate combinations """
        n= 10
        k= 10
        obs = n_choose_k(n,k)
        exp = 1
    	self.assertEqual(obs,exp)

    def test_hypergeometric_pmf(self):
        """ hypergeometric_pmf should calculate probability mass """
        k = 4	
        N = 50
        m = 5	
        n = 10
        obs = hypergeometric_pmf(k,N,m,n)
        exp = 0.0039645830580150
        self.assertFloatEqual(obs,exp)


    def test_hypergeometric_cdf(self):
        """ hypergeometric_cdf should calculate cummulative distribution """
        M = 100	
        x = 2
        K = 20	
        N = 10
        obs = hypergeometric_cdf(x,M,K,N)
        exp = 0.68122006381768885
        self.assertFloatEqual(obs,exp)

    def test_hypergeometric_sf(self):
        """ hypergeometric_pmf should calculate probability mass """
        M = 2	
        x = 100
        K = 20	
        N = 10
        obs = hypergeometric_sf(M,x,K,N)
        exp = 0.31877993618231115
        self.assertFloatEqual(obs,exp)

    def test_fisher_exact_test(self):
        """ fisher_exact_test should calculate Fisher exact test on a 2x2 contingency table  """
        input_array = ([[2, 5], [10, 4]])
        obs = fisher_exact_test(input_array)
        exp = 0.16, 0.15882352941176558
        self.assertFloatEqual(obs,exp)


if __name__ == '__main__':
    main()
