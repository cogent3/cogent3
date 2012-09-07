import numpy
from cogent.util.unit_test import TestCase, main

from cogent.maths.stats.period import chi_square, factorial, g_statistic, \
        circular_indices, _seq_to_symbols, seq_to_symbols, blockwise_bootstrap, \
        SeqToSymbols
from cogent.maths.period import ipdft, hybrid, auto_corr, Hybrid, Ipdft, \
        AutoCorrelation

__author__ = "Hua Ying, Julien Epps and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Julien Epps", "Hua Ying", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.2-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class TestPeriodStat(TestCase):
    def setUp(self):
        x = [1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 
            0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 
            1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 
            1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 
            1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 
            0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 
            0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 
            0, 0, 0, 0, 1, 1, 1, 0, 0, 0]
        self.x = numpy.array(x)
        self.sig = numpy.array(self.x, numpy.float64)
        self.motifs = ['AA', 'TT', 'TA']
    
    def test_chi_square(self):
        D, cs_p_val = chi_square(self.x, 10)
        self.assertEqual('%.4f'%D, '0.4786')
        self.assertEqual('%.4f'%cs_p_val, '0.4891')
    
    def test_factorial(self):
        self.assertEqual(factorial(1), 1)
        self.assertEqual(factorial(4), 24)
        self.assertEqual(factorial(0), 1)
    
    def test_g_statitic(self):
        """calc g-stat correctly"""
        X, periods = ipdft(self.sig, llim=2, ulim=39)
        g_obs, p_val = g_statistic(X)
        self.assertFloatEqual(p_val, 0.9997, eps=1e-3)
        self.assertFloatEqual(g_obs, 0.0577, eps=1e-3)
    
    def test_circular_indices(self):
        v = range(10)
        self.assertEqual(circular_indices(v, 8, 10, 4), [8,9,0,1])
        self.assertEqual(circular_indices(v, 9, 10, 4), [9,0,1,2])
        self.assertEqual(circular_indices(v, 4, 10, 4), [4,5,6,7])
    
    def test_seq_to_symbol(self):
        """both py and pyx seq_to_symbol versions correctly convert a sequence"""
        motifs = ['AA', 'AT', 'TT']
        symbols = _seq_to_symbols('AATGGTTA', motifs, 2)
        self.assertEqual(symbols, numpy.array([1,1,0,0,0,1,0,0]))
        symbols = seq_to_symbols('AAGATT', motifs, 2, numpy.zeros(6, numpy.uint8))
        self.assertEqual(symbols, numpy.array([1,0,0,1,1,0]))
    
    def test_seq_to_symbol_factory(self):
        """checks factory function for conversion works"""
        motifs = ['AA', 'AT', 'TT']
        seq_to_symbols = SeqToSymbols(motifs)
        self.assertEqual(seq_to_symbols('AATGGTTA'),
                        numpy.array([1,1,0,0,0,1,0,0]))
        self.assertEqual(seq_to_symbols('AAGATT'),
                        numpy.array([1,0, 0, 1, 1, 0], numpy.uint8))
    
    def test_permutation(self):
        s = 'ATCGTTGGGACCGGTTCAAGTTTTGGAACTCGCAAGGGGTGAATGGTCTTCGTCTAACGCTGG'\
            'GGAACCCTGAATCGTTGTAACGCTGGGGTCTTTAACCGTTCTAATTTAACGCTGGGGGGTTCT'\
            'AATTTTTAACCGCGGAATTGCGTC'
        seq_to_symbol = SeqToSymbols(self.motifs, length=len(s))
        hybrid_calc = Hybrid(len(s), llim=2, period = 4)
        ipdft_calc = Ipdft(len(s), llim=2, period = 4)
        stat, p = blockwise_bootstrap(s, hybrid_calc, block_size=10,
                    num_reps=1000, seq_to_symbols=seq_to_symbol)
        # print 's=%.4f; p=%.3f' % (stat, p)
        stat, p = blockwise_bootstrap(s, ipdft_calc, block_size=10,
                    num_reps=1000, seq_to_symbols=seq_to_symbol)
        # print 's=%.4f; p=%.3f' % (stat, p)
    
    def test_permutation_all(self):
        """performs permutation test of Hybrid, but considers all stats"""
        s = 'ATCGTTGGGACCGGTTCAAGTTTTGGAACTCGCAAGGGGTGAATGGTCTTCGTCTAACGCTGG'\
            'GGAACCCTGAATCGTTGTAACGCTGGGGTCTTTAACCGTTCTAATTTAACGCTGGGGGGTTCT'\
            'AATTTTTAACCGCGGAATTGCGTC'
        seq_to_symbol = SeqToSymbols(self.motifs, length=len(s))
        hybrid_calc = Hybrid(len(s), period = 4, return_all=True)
        stat, p = blockwise_bootstrap(s, hybrid_calc, block_size=10,
                    num_reps=1000, seq_to_symbols=seq_to_symbol)
        # print 's=%s; p=%s' % (stat, p)
    
    def test_get_num_stats(self):
        """calculators should return correct num stats"""
        hybrid_calc = Hybrid(150, llim=2, period = 4)
        ipdft_calc = Ipdft(150, llim=2, period = 4)
        autocorr_calc = AutoCorrelation(150, llim=2, period = 4)
        self.assertEqual(hybrid_calc.getNumStats(), 1)
        self.assertEqual(ipdft_calc.getNumStats(), 1)
        self.assertEqual(autocorr_calc.getNumStats(), 1)
        hybrid_calc = Hybrid(150, llim=2, period = 4, return_all=True)
        self.assertEqual(hybrid_calc.getNumStats(), 3)
    
    def test_permutation_skips(self):
        """permutation test correctly handles data without symbols"""
        s = 'N' * 150
        seq_to_symbol = SeqToSymbols(self.motifs, length=len(s))
        ipdft_calc = Ipdft(len(s), llim=2, period = 4)
        stat, p = blockwise_bootstrap(s, ipdft_calc, block_size=10,
                    num_reps=1000, seq_to_symbols=seq_to_symbol, num_stats=1)
        self.assertEqual(stat, 0.0)
        self.assertEqual(p, 1.0)
    

if __name__ == '__main__':
    main()
