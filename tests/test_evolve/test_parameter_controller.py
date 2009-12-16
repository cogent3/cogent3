#! /usr/bin/env python
# Matthew Wakefield Feb 2004

import unittest
import os

from cogent import LoadSeqs, LoadTree
import cogent.evolve.parameter_controller, cogent.evolve.substitution_model
from cogent.maths import optimisers

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, 'data')

good_rule_sets = [
    [
    {'par_name' : 'length','is_independent':True},
    ],
    [
    {'par_name' : 'length','is_independent':True},
    ],
    [
    {'par_name' : 'length','is_clade' :True, 'is_independent':True, 'edges' : ['a','b']},
    ],
    [
    {'par_name' : 'length','is_independent':True, 'edges' : ['a','c','e']},
    ],
    [
    {'par_name' : 'length','is_independent':True, 'edge' : 'a'},
    ],
]
bad_rule_sets = [
    [
    {'par_name' : 'length','is_clade' :True, 'edges' : ['b','f'],},
    ],
]
class test_parameter_controller(unittest.TestCase):
    """Tesing Parameter Controller"""
    def setUp(self):
        #length all edges 1 except c=2.  b&d transitions all other transverions
        self.al = LoadSeqs(
            data={'a':'tata', 'b':'tgtc', 'c':'gcga', 'd':'gaac', 'e':'gagc',})
        self.tree = LoadTree(treestring='((a,b),(c,d),e);')
        self.model = cogent.evolve.substitution_model.Nucleotide(
            do_scaling=True, equal_motif_probs=True, model_gaps=True)
        
    def test_scoped_local(self):
        model = cogent.evolve.substitution_model.Nucleotide(
                do_scaling=True, equal_motif_probs=True, model_gaps=True,
                predicates = {'kappa':'transition'})
        lf = model.makeLikelihoodFunction(self.tree)
        lf.setConstantLengths()
        lf.setAlignment(self.al)
        null = lf.getNumFreeParams()
        lf.setParamRule(par_name='kappa',
                            is_independent=True,
                            edges=['b','d'])
        self.assertEqual(null+2, lf.getNumFreeParams())
    
    def test_setMotifProbs(self):
        """Mprobs supplied to the parameter controller"""
        model = cogent.evolve.substitution_model.Nucleotide(
            model_gaps=True, motif_probs=None)
        lf = model.makeLikelihoodFunction(self.tree, 
                motif_probs_from_align=False)
                
        mprobs = {'A':0.1,'C':0.2,'G':0.2,'T':0.5,'-':0.0}
        lf.setMotifProbs(mprobs)
        self.assertEqual(lf.getMotifProbs(), mprobs)
        
        lf.setMotifProbsFromData(self.al[:1], is_const=True)
        self.assertEqual(lf.getMotifProbs()['G'], 0.6)
        
        lf.setMotifProbsFromData(self.al[:1], pseudocount=1)
        self.assertNotEqual(lf.getMotifProbs()['G'], 0.6)
        
        # test with consideration of ambiguous states
        al = LoadSeqs(data = {'seq1': 'ACGTAAGNA', 'seq2': 'ACGTANGTC',
                               'seq3': 'ACGTACGTG'})
        lf.setMotifProbsFromData(al, include_ambiguity=True, is_const=True)
        motif_probs = dict(lf.getMotifProbs())
        correct_probs = {'A': 8.5/27, 'C': 5.5/27, '-': 0.0, 'T': 5.5/27,
                         'G': 7.5/27}
        self.assertEqual(motif_probs, correct_probs)
        self.assertEqual(sum(motif_probs.values()), 1.0)

    def test_setMultiLocus(self):
        """2 loci each with own mprobs"""
        model = cogent.evolve.substitution_model.Nucleotide(motif_probs=None)
        lf = model.makeLikelihoodFunction(self.tree, 
                motif_probs_from_align=False, loci=["a", "b"])
                
        mprobs_a = dict(A=.2, T=.2, C=.3, G=.3)
        mprobs_b = dict(A=.1, T=.2, C=.3, G=.4)
        
        for is_const in [False, True]:
            lf.setMotifProbs(mprobs_a, is_const=is_const)
            s = str(lf)
            lf.setMotifProbs(mprobs_b, locus="b")
            self.assertEqual(lf.getMotifProbs(locus="a"), mprobs_a)
            self.assertEqual(lf.getMotifProbs(locus="b"), mprobs_b)
            s = str(lf)
            #lf.setParamRule('mprobs', is_independent=False)

    def test_setParamRules(self):
        lf = self.model.makeLikelihoodFunction(self.tree)
        def do_rules(rule_set):
            for rule in rule_set:
                lf.setParamRule(**rule)
        for rule_set in good_rule_sets:
            lf.setDefaultParamRules()
            do_rules(rule_set)
        for rule_set in bad_rule_sets:
            lf.setDefaultParamRules()
            self.assertRaises((KeyError, TypeError,
                    AssertionError, ValueError), do_rules, rule_set)
    
    def test_setLocalClock(self):
        pass

    def test_setConstantLengths(self):
        t = LoadTree(treestring='((a:1,b:2):3,(c:4,d:5):6,e:7);')
        lf = self.model.makeLikelihoodFunction(t)#self.tree)
        lf.setParamRule('length', is_const=True)
        # lf.setConstantLengths(t)
        lf.setAlignment(self.al)
        self.assertEqual(lf.getParamValue('length', 'b'), 2)
        self.assertEqual(lf.getParamValue('length', 'd'), 5)

    def test_pairwise_clock(self):
        al = LoadSeqs(data={'a':'agct','b':'ggct'})
        tree = LoadTree(treestring='(a,b);')
        model = cogent.evolve.substitution_model.Dinucleotide(
                do_scaling=True, equal_motif_probs=True, model_gaps=True,
                mprob_model='tuple')
        lf = model.makeLikelihoodFunction(tree)
        lf.setLocalClock('a','b')
        lf.setAlignment(al)
        lf.optimise(local=True, show_progress=False)
        rd = lf.getStatisticsAsDict()
        self.assertAlmostEquals(lf.getLogLikelihood(),-10.1774488956)
        self.assertEqual(rd['length']['a'],rd['length']['b'])

    def test_local_clock(self):
        lf = self.model.makeLikelihoodFunction(self.tree)
        lf.setLocalClock('c','d')
        lf.setAlignment(self.al)
        lf.optimise(local=True, show_progress=False, 
                tolerance=1e-8, max_restarts=2)
        rd = lf.getStatisticsAsDict()
        self.assertAlmostEquals(lf.getLogLikelihood(),-27.84254174)
        self.assertEqual(rd['length']['c'],rd['length']['d'])
        self.assertNotEqual(rd['length']['a'],rd['length']['e'])
        
    def test_complex_parameter_rules(self):
            # This test has many local minima and so does not cope
            # with changes to optimiser details.
        model = cogent.evolve.substitution_model.Nucleotide(
                do_scaling=True, equal_motif_probs=True, model_gaps=True,
                predicates = {'kappa':'transition'})
        lf = model.makeLikelihoodFunction(self.tree)
        lf.setParamRule(par_name='kappa',
                                is_independent=True)
        lf.setParamRule(par_name='kappa',
                                is_independent=False,
                                edges=['b','d'])
        lf.setConstantLengths(LoadTree(
                            treestring='((a:1,b:1):1,(c:2,d:1):1,e:1);'))
        #print self.pc
        lf.setAlignment(self.al)
        lf.optimise(show_progress=False, local=True)
        rd = lf.getStatisticsAsDict()
        self.assertAlmostEquals(lf.getLogLikelihood(),-27.3252, 3)
        self.assertEqual(rd['kappa']['b'],rd['kappa']['d'])
        self.assertNotEqual(rd['kappa']['a'],rd['kappa']['b'])
        
        
if __name__ == '__main__':
    unittest.main()

    
    
        
