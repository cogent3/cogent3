#!/usr/bin/env python

import warnings

warnings.filterwarnings("ignore", "Motif probs overspecified")
warnings.filterwarnings("ignore", "Model not reversible")

from numpy import ones, dot

from cogent import LoadSeqs, DNA, LoadTree, LoadTable
from cogent.evolve.substitution_model import Nucleotide
from cogent.evolve.predicate import MotifChange
from cogent.util.unit_test import TestCase, main
from cogent.maths.matrix_exponentiation import PadeExponentiator as expm

__author__ = "Peter Maxwell and  Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def _dinuc_root_probs(x,y=None):
    if y is None:
        y = x
    return dict([(n1+n2, p1*p2) 
            for n1,p1 in x.items() for n2,p2 in y.items()])

def _trinuc_root_probs(x,y,z):
    return dict([(n1+n2+n3, p1*p2*p3) 
            for n1,p1 in x.items() for n2,p2 in y.items()
            for n3,p3 in z.items()])

def make_p(length, coord, val):
    """returns a probability matrix with value set at coordinate in
    instantaneous rate matrix"""
    Q = ones((4,4), float)*0.25 # assumes equi-frequent mprobs at root
    for i in range(4):
        Q[i,i] = 0.0
    Q[coord] *= val
    row_sum = Q.sum(axis=1)
    scale = 1/(.25*row_sum).sum()
    for i in range(4):
        Q[i,i] -= row_sum[i]
    Q *= scale
    return expm(Q)(length)

class NewQ(TestCase):
    aln = LoadSeqs(data={
    'seq1': 'TGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT',
    'seq2': 'TGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT'},
    moltype=DNA)
    tree = LoadTree(tip_names=['seq1', 'seq2'])
    
    symm_nuc_probs = dict(A=0.25,T=0.25,C=0.25,G=0.25)
    symm_root_probs = _dinuc_root_probs(symm_nuc_probs)
    asymm_nuc_probs = dict(A=0.1,T=0.1,C=0.4,G=0.4)
    asymm_root_probs = _dinuc_root_probs(asymm_nuc_probs)
    posn_root_probs = _dinuc_root_probs(symm_nuc_probs, asymm_nuc_probs)
    cond_root_probs = dict([(n1+n2, p1*[.1, .7][n1==n2]) 
            for n1,p1 in asymm_nuc_probs.items() for n2 in 'ATCG'])
    
    # Each of these (data, model) pairs should give a result different 
    # from any of the simpler models applied to the same data.  
    ordered_by_complexity = [
            # P(AA) == P(GG) == P(AG)
            [symm_root_probs, 'tuple'],
            
            # P(GA) == P(AG) but P(AA) != P(GG)
            [asymm_root_probs, 'monomer'],
            
            # P(AG) == P(A?)*P(?G) but P(A?) != P(?A)
            [posn_root_probs, 'monomers'],
            
            # P(AG) != P(A?)*P(?G)
            [cond_root_probs, 'conditional'],
            ]
    
    def test_newQ_is_nuc_process(self):
        """newQ is an extension of an independent nucleotide process"""
        nuc = Nucleotide(motif_probs = self.asymm_nuc_probs)
        new_di = Nucleotide(motif_length=2, mprob_model='monomer',
            motif_probs = self.asymm_root_probs)
        
        nuc_lf = nuc.makeLikelihoodFunction(self.tree)
        new_di_lf = new_di.makeLikelihoodFunction(self.tree)
        # newQ branch length is exactly motif_length*nuc branch length
        nuc_lf.setParamRule('length', is_independent=False, init=0.2)
        new_di_lf.setParamRule('length', is_independent=False, init=0.4)
        
        nuc_lf.setAlignment(self.aln)
        new_di_lf.setAlignment(self.aln)
        self.assertFloatEqual(nuc_lf.getLogLikelihood(),
                                new_di_lf.getLogLikelihood())
    
    def test_lf_display(self):
        """str of likelihood functions should not fail"""
        for (dummy, model) in self.ordered_by_complexity:
            di = Nucleotide(motif_length=2, mprob_model=model)
            di.adaptMotifProbs(self.cond_root_probs, auto=True)
            lf = di.makeLikelihoodFunction(self.tree)
            s = str(lf)
    
    def test_get_statistics(self):
        """get statistics should correctly apply arguments"""
        for (mprobs, model) in self.ordered_by_complexity:
            di = Nucleotide(motif_length=2, motif_probs=mprobs, 
                    mprob_model=model)
            lf = di.makeLikelihoodFunction(self.tree)
            for wm, wt in [(True, True), (True, False), (False, True),
                           (False, False)]:
                stats = lf.getStatistics(with_motif_probs=wm, with_titles=wt)
    
    def test_sim_alignment(self):
        """should be able to simulate an alignment under all models"""
        for (mprobs, model) in self.ordered_by_complexity:
            di = Nucleotide(motif_length=2, motif_probs=mprobs, 
                    mprob_model=model)
            lf = di.makeLikelihoodFunction(self.tree)
            lf.setParamRule('length', is_independent=False, init=0.4)
            lf.setAlignment(self.aln)
            sim = lf.simulateAlignment()
    
    def test_reconstruct_ancestor(self):
        """should be able to reconstruct ancestral sequences under all
        models"""
        for (mprobs, model) in self.ordered_by_complexity:
            di = Nucleotide(motif_length=2, mprob_model=model)
            di.adaptMotifProbs(mprobs, auto=True)
            lf = di.makeLikelihoodFunction(self.tree)
            lf.setParamRule('length', is_independent=False, init=0.4)
            lf.setAlignment(self.aln)
            ancestor = lf.reconstructAncestralSeqs()
    
    def test_results_different(self):
        for (i, (mprobs, dummy)) in enumerate(self.ordered_by_complexity):
            results = []
            for (dummy, model) in self.ordered_by_complexity:
                di = Nucleotide(motif_length=2, motif_probs=mprobs, 
                        mprob_model=model)
                lf = di.makeLikelihoodFunction(self.tree)
                lf.setParamRule('length', is_independent=False, init=0.4)
                lf.setAlignment(self.aln)
                lh = lf.getLogLikelihood()
                for other in results[:i]:
                    self.failIfAlmostEqual(other, lh, places=2)
                for other in results[i:]:
                    self.assertFloatEqual(other, lh)
                results.append(lh)
    
    def test_position_specific_mprobs(self):
        """correctly compute likelihood when positions have distinct
        probabilities"""
        aln_len = len(self.aln)
        posn1 = []
        posn2 = []
        for name, seq in self.aln.todict().items():
            p1 = [seq[i] for i in range(0,aln_len,2)]
            p2 = [seq[i] for i in range(1,aln_len,2)]
            posn1.append([name, ''.join(p1)])
            posn2.append([name, ''.join(p2)])
        
        # the position specific alignments
        posn1 = LoadSeqs(data=posn1)
        posn2 = LoadSeqs(data=posn2)
        
        # a newQ dinucleotide model
        sm = Nucleotide(motif_length=2, mprob_model='monomer', do_scaling=False)
        lf = sm.makeLikelihoodFunction(self.tree)
        lf.setAlignment(posn1)
        posn1_lnL = lf.getLogLikelihood()
        lf.setAlignment(posn2)
        posn2_lnL = lf.getLogLikelihood()
        expect_lnL = posn1_lnL+posn2_lnL
        
        # the joint model
        lf.setAlignment(self.aln)
        aln_lnL = lf.getLogLikelihood()
        
        # setting the full alignment, which has different motif probs, should
        # produce a different lnL
        self.failIfAlmostEqual(expect_lnL, aln_lnL)
        
        # set the arguments for taking position specific mprobs
        sm = Nucleotide(motif_length=2, mprob_model='monomers',
                        do_scaling=False)
        lf = sm.makeLikelihoodFunction(self.tree)
        lf.setAlignment(self.aln)
        posn12_lnL = lf.getLogLikelihood()
        self.assertFloatEqual(expect_lnL, posn12_lnL)
    
    def test_compute_conditional_mprobs(self):
        """equal likelihood from position specific and conditional mprobs"""
        def compare_models(motif_probs, motif_length):
            # if the 1st and 2nd position motifs are independent of each other
            # then conditional is the same as positional
            ps = Nucleotide(motif_length=motif_length, motif_probs=motif_probs,
                mprob_model='monomers')
            cd = Nucleotide(motif_length=motif_length,motif_probs=motif_probs,
                            mprob_model='conditional')
            
            ps_lf = ps.makeLikelihoodFunction(self.tree)
            ps_lf.setParamRule('length', is_independent=False, init=0.4)
            ps_lf.setAlignment(self.aln)
            
            cd_lf = cd.makeLikelihoodFunction(self.tree)
            cd_lf.setParamRule('length', is_independent=False, init=0.4)
            cd_lf.setAlignment(self.aln)
            self.assertFloatEqual(cd_lf.getLogLikelihood(),
                    ps_lf.getLogLikelihood())
        
        compare_models(self.posn_root_probs, 2)
        # trinucleotide
        trinuc_mprobs = _trinuc_root_probs(self.asymm_nuc_probs,
                            self.asymm_nuc_probs, self.asymm_nuc_probs)
        compare_models(trinuc_mprobs, 3)
    
    def test_cond_pos_differ(self):
        """lnL should differ when motif probs are not multiplicative"""
        dinuc_probs = {'AA': 0.088506666666666664, 'AC': 0.044746666666666664,
            'GT': 0.056693333333333332, 'AG': 0.070199999999999999,
            'CC': 0.048653333333333333, 'TT': 0.10678666666666667,
            'CG': 0.0093600000000000003, 'GG': 0.049853333333333333,
            'GC': 0.040253333333333335, 'AT': 0.078880000000000006,
            'GA': 0.058639999999999998, 'TG': 0.081626666666666667,
            'TA': 0.068573333333333333, 'CA': 0.06661333333333333,
            'TC': 0.060866666666666666, 'CT': 0.069746666666666665}
        
        mg = Nucleotide(motif_length=2, motif_probs=dinuc_probs,
                        mprob_model='monomer')
        mg_lf = mg.makeLikelihoodFunction(self.tree)
        mg_lf.setParamRule('length', is_independent=False, init=0.4)
        mg_lf.setAlignment(self.aln)
        
        cd = Nucleotide(motif_length=2, motif_probs=dinuc_probs,
                        mprob_model='conditional')
        
        cd_lf = cd.makeLikelihoodFunction(self.tree)
        cd_lf.setParamRule('length', is_independent=False, init=0.4)
        cd_lf.setAlignment(self.aln)
        self.assertNotAlmostEqual(mg_lf.getLogLikelihood(),
                                    cd_lf.getLogLikelihood())
    
    def test_getting_node_mprobs(self):
        """return correct motif probability vector for tree nodes"""
        tree = LoadTree(treestring='(a:.2,b:.2,(c:.1,d:.1):.1)')
        aln = LoadSeqs(data={
        'a': 'TGTG',
        'b': 'TGTG',
        'c': 'TGTG',
        'd': 'TGTG',
        })
        
        motifs = ['T', 'C', 'A', 'G']
        aX = MotifChange(motifs[0], motifs[3], forward_only=True).aliased('aX')
        bX = MotifChange(motifs[3], motifs[0], forward_only=True).aliased('bX')
        edX = MotifChange(motifs[1], motifs[2], forward_only=True).aliased('edX')
        cX = MotifChange(motifs[2], motifs[1], forward_only=True).aliased('cX')
        sm = Nucleotide(predicates=[aX, bX, edX, cX], equal_motif_probs=True)
        
        lf = sm.makeLikelihoodFunction(tree)
        lf.setParamRule('aX', edge='a', value=8.0)
        lf.setParamRule('bX', edge='b', value=8.0)
        lf.setParamRule('edX', edge='edge.0', value=2.0)
        lf.setParamRule('cX', edge='c', value=0.5)
        lf.setParamRule('edX', edge='d', value=4.0)
        lf.setAlignment(aln)
        
        # we construct the hand calc variants
        mprobs = ones(4, float) * .25
        a = make_p(.2, (0,3), 8)
        a = dot(mprobs, a)
        
        b = make_p(.2, (3, 0), 8)
        b = dot(mprobs, b)
        
        e = make_p(.1, (1, 2), 2)
        e = dot(mprobs, e)
        
        c = make_p(.1, (2, 1), 0.5)
        c = dot(e, c)
        
        d = make_p(.1, (1, 2), 4)
        d = dot(e, d)
        
        prob_vectors = lf.getMotifProbsByNode()
        self.assertFloatEqual(prob_vectors['a'].array, a)
        self.assertFloatEqual(prob_vectors['b'].array, b)
        self.assertFloatEqual(prob_vectors['c'].array, c)
        self.assertFloatEqual(prob_vectors['d'].array, d)
        self.assertFloatEqual(prob_vectors['edge.0'].array, e)
    

if __name__ == "__main__":
    main()

