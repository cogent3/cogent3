#!/usr/bin/env python
"""Unit tests of the basic CAI calculations."""
from cogent.util.unit_test import TestCase, main
from math import log, exp
from operator import mul
from cogent.maths.stats.cai.util import cu, as_rna, synonyms_to_rna, \
    get_synonyms, sum_codon_freqs, norm_to_max, arithmetic_mean, \
    geometric_mean, codon_adaptiveness_all, codon_adaptiveness_blocks, \
    valid_codons, set_min, cai_1, cai_2, cai_3

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"


def product(x): return reduce(mul, x)
def amean(x): return sum(x)/float(len(x))
def gmean(x): return (product(x))**(1./len(x))

class cai_tests(TestCase):
    """Tests of top-level functionality."""

    def test_as_rna(self):
        """as_rna should do correct conversion to RNA"""
        self.assertEqual(as_rna('TCGT'), 'UCGU')

    def test_synonyms_to_rna(self):
        """synonyms_to_rna should convert as expected"""
        s = {'*':['TAA','TAG'], 'F':['TTT','TTC']}
        self.assertEqual(synonyms_to_rna(s), {'*':['UAA','UAG'],'F':['UUU','UUC']})

    def test_synonyms(self):
        """synonyms should produce correct results for standard genetic code.
        
        NOTE: for standard genetic code, expect the following:
        - Stop codons are UGA, UAA, UAG
        - Single-codon blocs are M = 'AUG', W = 'UGG'
        """
        result = get_synonyms()
        self.assertEqual(len(result), 18)
        self.assertEqual(''.join(sorted(result)), 'ACDEFGHIKLNPQRSTVY')
        self.assertEqual(sorted(result['I']), ['AUA','AUC','AUU'])

        #check that we can do it without eliminating single-codon blocks
        result = get_synonyms(singles_removed=False)
        self.assertEqual(len(result), 20)
        self.assertEqual(''.join(sorted(result)), 'ACDEFGHIKLMNPQRSTVWY')
        self.assertEqual(sorted(result['I']), ['AUA','AUC','AUU'])
        self.assertEqual(result['W'], ['UGG'])

    def test_sum_codon_freqs(self):
        """sum_codon_freqs should add list of codon freqs together, incl, missing keys"""
        d = {'x':3, 'UUU':5, 'UAC':3}
        d2 = {'y':5, 'UUU':1, 'AGG':2}
        result = sum_codon_freqs([d,d2])
        self.assertEqual(len(result), 64)
        assert 'x' not in result    #should exclude bad keys
        self.assertEqual(result['UUU'], 6.0)
        self.assertEqual(result['AGG'], 2.0)
        self.assertEqual(result['UAC'], 3.0)
        self.assertEqual(sum(result.values()), 11.0)

    def test_norm_to_max(self):
        """norm_to_max should normalize vals in list to best val"""
        a = [1,2,3,4]
        self.assertEqual(norm_to_max(a), [.25,.5,.75,1])

    def test_arithmetic_mean(self):
        """arithmetic_mean should average a list of means with freqs"""
        obs = arithmetic_mean([1,2,3],[2,1,3])
        exp = sum([1,1,2,3,3,3])/6.0
        self.assertEqual(obs, exp)
        #should also work without freqs
        self.assertFloatEqual(arithmetic_mean([1,3,7]), 11/3.)

    def test_geometric_mean(self):
        """geometric_mean should average a list of means with freqs"""
        obs = geometric_mean([1,2,3],[2,1,3])
        exp = (1*1*2*3*3*3)**(1/6.)
        self.assertEqual(obs, exp)
        obs = geometric_mean([0.01, 0.2, 0.5], [5, 2, 3])
        exp = (.01*.01*.01*.01*.01*.2*.2*.5*.5*.5)**(0.1)
        self.assertFloatEqual(obs,exp)
        #should also work without freqs
        self.assertFloatEqual(geometric_mean([0.01,0.2,0.5]), (.01*.2*.5)**(1/3.))

    def test_codon_adaptiveness_all(self):
        """codon_adaptiveness_all should normalize all codons relative to the best one."""
        codons = {'x':4, 'y':3, 'z':2, 'zz':0, 'zzz':4}
        result = codon_adaptiveness_all(codons)
        self.assertEqual(result, {'x':1., 'y':.75, 'z':.5, 'zz':0, 'zzz':1.})

    def test_codon_adaptiveness_blocks(self):
        """codon_adaptiveness_blocks should normalize codons by the best in each block"""
        codons = {'x':4, 'y':1, 'z':2, 'zz':0, 'zzz':2, 'zzzz':1}
        blocks = {'A': ['x','y','z'], 'B':['zz','zzz','zzzz']}
        result = codon_adaptiveness_blocks(codons, blocks)
        self.assertEqual(result, {'x':1., 'y':.25, 'z':.5, 'zz':0, 'zzz':1., 'zzzz':.5})

    def test_set_min(self):
        """set_min should set minimum value to specified threshold."""
        codons = {'x':4, 'y':1e-5, 'z':0}
        set_min(codons, 1)
        self.assertEqual(codons, {'x':4, 'y':1, 'z':1})

    def test_valid_codons(self):
        """valid_codons should extract all valid codons from blocks"""
        blocks = {'A':['GCA','GCG'], 'C':['UGU','UGC']}
        self.assertEqual(list(sorted(valid_codons(blocks))), ['GCA','GCG','UGC','UGU'])

    def test_cai_1(self):
        """cai_1 should produce expected results"""
        ref_freqs = cu.copy()
        ref_freqs.update({'AGA':4, 'AGG':2, 'CCC':4, 'CCA':1, 'UGG':1})
        #tests with arithmetic mean
        gene_freqs = {'AGA':1}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs = {'AGA':5}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs = {'AGG':5}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=arithmetic_mean), 0.5)
        gene_freqs = {'AGG':5,'AGA':5}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=arithmetic_mean), 0.75)
        gene_freqs={'AGA':5,'CCC':1}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs={'AGA':5,'CCA':5}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=arithmetic_mean), 0.625)
        ref_freqs_2 = cu.copy()
        ref_freqs_2.update({'AGA':4, 'AGG':2, 'CCC':5, 'CCA':1, 'UGG':1})
        ref_freqs_2.update({'UUU':2,'UUC':1})
        gene_freqs = {'AGA':3,'AGG':1,'CCC':2,'CCA':1,'UUU':1, 'UUC':2}
        obs = cai_1(ref_freqs_2, gene_freqs, average=arithmetic_mean)
        vals = [.8,.8,.8,.4,1,1,.2,.4,.2,.2]
        expect = sum(vals)/len(vals)
        self.assertFloatEqual(obs, expect)
        #tests with geometric mean
        gene_freqs = {'AGA':1}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs = {'AGA':5}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs = {'AGG':5}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=geometric_mean), 0.5)
        gene_freqs = {'AGG':5,'AGA':5}
        self.assertFloatEqual(cai_1(ref_freqs, gene_freqs, average=geometric_mean), \
            (1**5 * 0.5**5)**(0.1))
        gene_freqs={'AGA':5,'CCC':1}
        self.assertEqual(cai_1(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs={'AGA':5,'CCA':5}
        self.assertFloatEqual(cai_1(ref_freqs, gene_freqs, average=geometric_mean), \
            (1**5 * 0.25**5)**0.1)
        ref_freqs_2 = cu.copy()
        ref_freqs_2.update({'AGA':4, 'AGG':2, 'CCC':5, 'CCA':1, 'UGG':1})
        ref_freqs_2.update({'UUU':2,'UUC':1})
        gene_freqs = {'AGA':3,'AGG':1,'CCC':2,'CCA':1,'UUU':1, 'UUC':2}
        obs = cai_1(ref_freqs_2, gene_freqs, average=geometric_mean)
        vals = [.8,.8,.8,.4,1,1,.2,.4,.2,.2]
        expect = (product(vals))**(1./len(vals))
        self.assertFloatEqual(obs, expect)

    def test_cai_2(self):
        """cai_2 should produce expected results"""
        ref_freqs = cu.copy()
        ref_freqs.update({'AGA':4, 'AGG':2, 'CCC':5, 'CCA':1, 'UGG':1})
        #tests with arithmetic mean
        gene_freqs = {'AGA':1}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs = {'AGA':5}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs = {'AGG':5}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=arithmetic_mean), 0.5)
        gene_freqs = {'AGG':5,'AGA':5}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=arithmetic_mean), 0.75)
        gene_freqs={'AGA':5,'CCC':1}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs={'AGA':5,'CCA':5}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=arithmetic_mean), 0.6)
        ref_freqs_2 = ref_freqs.copy()
        ref_freqs_2.update({'UUU':2,'UUC':1})
        gene_freqs = {'AGA':3,'AGG':1,'CCC':2,'CCA':1,'UUU':1, 'UUC':2}
        obs = cai_2(ref_freqs_2, gene_freqs, average=arithmetic_mean)
        vals = [1,1,1,.5,1,1,.2,1,.5,.5]
        expect = sum(vals)/len(vals)
        self.assertEqual(obs, expect)
        #tests with geometric mean
        gene_freqs = {'AGA':1}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs = {'AGA':5}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs = {'AGG':5}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=geometric_mean), 0.5)
        gene_freqs = {'AGG':5,'AGA':5}
        self.assertFloatEqual(cai_2(ref_freqs, gene_freqs, average=geometric_mean), \
            (1**5 * 0.5**5)**(0.1))
        gene_freqs={'AGA':5,'CCC':1}
        self.assertEqual(cai_2(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs={'AGA':5,'CCA':5}
        self.assertFloatEqual(cai_2(ref_freqs, gene_freqs, average=geometric_mean), \
            (1**5 * 0.2**5)**0.1)
        ref_freqs_2 = ref_freqs.copy()
        ref_freqs_2.update({'UUU':2,'UUC':1})
        gene_freqs = {'AGA':3,'AGG':1,'CCC':2,'CCA':1,'UUU':1, 'UUC':2}
        obs = cai_2(ref_freqs_2, gene_freqs, average=geometric_mean)
        vals = [1,1,1,.5,1,1,.2,1,.5,.5]
        expect = (product(vals))**(1./len(vals))
        self.assertEqual(obs, expect)
        #test that results match example on Gang Wu's CAI calculator page
        ref_freqs = cu.copy()
        ref_freqs.update({'UUU':78743, 'UUC':56591, 'UUA':51320, 'UUG':45581, \
            'CUU':42704, 'CUC':35873, 'CUA':15275, 'CUG':168885})
        gene_freqs={'UUU':6, 'UUC':3, 'CUU':3, 'CUC':2, 'CUG':8}
        self.assertFloatEqual(cai_2(ref_freqs, gene_freqs, average=geometric_mean), \
            exp((6*log(1) + 3*log(56591./78743) + 3*log(42704./168885) + \
            2*log(35873./168885)+8*log(1))/22.))

    def test_cai_3(self):
        """cai_3 should produce expected results"""
        ref_freqs = cu.copy()
        ref_freqs.update({'AGA':4, 'AGG':2, 'CCC':5, 'CCA':1, 'UGG':1})
        #tests with arithmetic mean
        gene_freqs = {'AGA':1}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs = {'AGA':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs = {'AGG':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=arithmetic_mean), 0.5)
        gene_freqs = {'AGG':5,'AGA':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=arithmetic_mean), 0.75)
        gene_freqs={'AGA':5,'CCC':1}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=arithmetic_mean), 1)
        gene_freqs={'AGA':5,'CCA':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=arithmetic_mean), 0.6)
        ref_freqs_2 = ref_freqs.copy()
        ref_freqs_2.update({'UUU':2,'UUC':1})
        gene_freqs = {'AGA':3,'AGG':1,'CCC':2,'CCA':1,'UUU':1, 'UUC':2}
        obs = cai_3(ref_freqs_2, gene_freqs, average=arithmetic_mean)
        family_vals = [[1,1,1,.5],[1,1,.2],[1,.5,.5]]
        family_averages = map(amean, family_vals)
        expect = amean(family_averages)
        self.assertEqual(obs, expect)
        #tests with geometric mean
        gene_freqs = {'AGA':1}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs = {'AGA':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs = {'AGG':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=geometric_mean), 0.5)
        gene_freqs = {'AGG':5,'AGA':5}
        self.assertFloatEqual(cai_3(ref_freqs, gene_freqs, average=geometric_mean), \
            (1**5 * 0.5**5)**(0.1))
        gene_freqs={'AGA':5,'CCC':1}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average=geometric_mean), 1)
        gene_freqs={'AGA':5,'CCA':5}
        self.assertFloatEqual(cai_3(ref_freqs, gene_freqs, average=geometric_mean), \
            (1**5 * 0.2**5)**0.1)
        ref_freqs_2 = ref_freqs.copy()
        ref_freqs_2.update({'UUU':2,'UUC':1})
        gene_freqs = {'AGA':3,'AGG':1,'CCC':2,'CCA':1,'UUU':1, 'UUC':2}
        obs = cai_3(ref_freqs_2, gene_freqs, average=geometric_mean)
        family_vals = [[1,1,1,.5],[1,1,.2],[1,.5,.5]]
        family_averages = map(gmean, family_vals)
        expect = gmean(family_averages)
        self.assertEqual(obs, expect)
        #tests with Eyre-Walker's variant -- should be same as geometric mean
        gene_freqs = {'AGA':1}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average='eyre_walker'), 1)
        gene_freqs = {'AGA':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average='eyre_walker'), 1)
        gene_freqs = {'AGG':5}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average='eyre_walker'), 0.5)
        gene_freqs = {'AGG':5,'AGA':5}
        self.assertFloatEqual(cai_3(ref_freqs, gene_freqs, average='eyre_walker'), \
            (1**5 * 0.5**5)**(0.1))
        gene_freqs={'AGA':5,'CCC':1}
        self.assertEqual(cai_3(ref_freqs, gene_freqs, average='eyre_walker'), 1)
        gene_freqs={'AGA':5,'CCA':5}
        self.assertFloatEqual(cai_3(ref_freqs, gene_freqs, average='eyre_walker'), \
            (1**5 * 0.2**5)**0.1)
        ref_freqs_2 = ref_freqs.copy()
        ref_freqs_2.update({'UUU':2,'UUC':1})
        gene_freqs = {'AGA':3,'AGG':1,'CCC':2,'CCA':1,'UUU':1, 'UUC':2}
        obs = cai_3(ref_freqs_2, gene_freqs, average='eyre_walker')
        family_vals = [[1,1,1,.5],[1,1,.2],[1,.5,.5]]
        family_averages = map(gmean, family_vals)
        expect = gmean(family_averages)
        self.assertEqual(obs, expect)
        #test results for Gang Wu's example (unfortunately, no worked example for
        #this model)
        ref_freqs = cu.copy()
        ref_freqs.update({'UUU':78743, 'UUC':56591, 'UUA':51320, 'UUG':45581, \
            'CUU':42704, 'CUC':35873, 'CUA':15275, 'CUG':168885})
        gene_freqs={'UUU':6, 'UUC':3, 'CUU':3, 'CUC':2, 'CUG':8}
        obs = cai_3(ref_freqs, gene_freqs, average=geometric_mean)
        family_vals =  [6*[1]+3*[56591./78743],\
            3*[42704./168885] + 2*[35873./168885]+8*[1]]
        family_averages = map(gmean, family_vals)
        expect = gmean(family_averages)
        self.assertFloatEqual(obs, expect)

if __name__ == '__main__':
    main()
