#!/usr/bin/env python

from cogent import DNA
from cogent.align.align import classic_align_pairwise, make_dna_scoring_dict
from cogent.evolve.models import HKY85
import cogent.evolve.substitution_model
dna_model = cogent.evolve.substitution_model.Nucleotide(
        model_gaps=False, equal_motif_probs=True)

import cogent.align.progressive

import unittest

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def matchedColumns(align):
    """Count the matched columns in an alignment"""
    def all_same(column):
        consensus = None
        for motif in column:
            if consensus is None:
                consensus = motif
            elif motif != consensus:
                return False
        return True
    
    return len(align.filtered(all_same))

seq1 = DNA.makeSequence('aaaccggacattacgtgcgta', Name='FAKE01')
seq2 = DNA.makeSequence( 'ccggtcaggttacgtacgtt', Name= 'FAKE02')

class AlignmentTestCase(unittest.TestCase):
    def _aligned_both_ways(self, seq1, seq2, **kw):
        S = make_dna_scoring_dict(10, -1, -8)
        a1 = classic_align_pairwise(seq1, seq2, S, 10, 2, **kw)
        a2 = classic_align_pairwise(seq2, seq1, S, 10, 2, **kw)
        return [a1, a2]
    
    def test_local(self):
        for a in self._aligned_both_ways(seq1, seq2, local=True):
            self.assertEqual(matchedColumns(a), 15)
            self.assertEqual(len(a), 19)
    
    def test_gap_at_one_end(self):
        for a in self._aligned_both_ways(seq1, seq2, local=False):
            self.assertEqual(matchedColumns(a), 15)
            self.assertEqual(len(a), 23)
    
    def test_gaps_at_both_ends(self):
        s = 'aaaccggttt'
        s1 = DNA.makeSequence(s[:-2], Name="A")
        s2 = DNA.makeSequence(s[2:], Name="B")
        for a in self._aligned_both_ways(s1, s2, local=False):
            self.assertEqual(matchedColumns(a), 6)
            self.assertEqual(len(a), 10)
    
    def test_short(self):
        s1 = DNA.makeSequence('tacagta', Name="A")
        s2 = DNA.makeSequence('tacgtc', Name="B")
        for a in self._aligned_both_ways(s1, s2, local=False):
            self.assertEqual(matchedColumns(a), 5)
            self.assertEqual(len(a), 7)
    
    def test_codon(self):
        s1 = DNA.makeSequence('tacgccgta', Name="A")
        s2 = DNA.makeSequence('tacgta', Name="B")
        codon_model = cogent.evolve.substitution_model.Codon(
                                 model_gaps=False, equal_motif_probs=True,
                                 mprob_model='conditional')
        tree = cogent.LoadTree(tip_names=['A', 'B'])
        lf = codon_model.makeLikelihoodFunction(tree, aligned=False)
        lf.setSequences(dict(A=s1, B=s2))
        (score, a) = lf.getLogLikelihood().edge.getViterbiScoreAndAlignment()
        self.assertEqual(matchedColumns(a), 6)
        self.assertEqual(len(a), 9)
    

class UnalignedPairTestCase(unittest.TestCase):
    def test_forward(self):
        tree = cogent.LoadTree(tip_names='AB')
        pc = dna_model.makeLikelihoodFunction(tree, aligned=False)  
        pc.setSequences({'A':seq1, 'B':seq2})
        LnL = pc.getLogLikelihood()
        assert isinstance(LnL, float)
    

class MultipleAlignmentTestCase(unittest.TestCase):
    def _make_aln(self, orig, model=dna_model, param_vals=None, 
            indel_rate=0.1, indel_length=0.5, **kw):
        kw['indel_rate'] = indel_rate
        kw['indel_length'] = indel_length
        seqs = dict((key, DNA.makeSequence(value)) 
                for (key, value) in orig.items())
        if len(seqs) == 2:
            tree = cogent.LoadTree(tip_names=seqs.keys())
            tree = cogent.LoadTree(treestring="(A:.1,B:.1)")
        else:
            tree = cogent.LoadTree(treestring="(((A:.1,B:.1):.1,C:.1):.1,D:.1)")
        aln, tree = cogent.align.progressive.TreeAlign(model, seqs,
                tree=tree, show_progress=False, param_vals=param_vals, **kw)
        return aln
    
    def _test_aln(self, seqs, model=dna_model, param_vals=None, **kw):
        orig = dict((n,s.replace('-', '')) for (n,s) in seqs.items())
        aln = self._make_aln(orig, model=model, param_vals=param_vals, **kw)
        result = dict((n,s.lower()) for (n,s) in aln.todict().items())
        # assert the alignment result is correct
        self.assertEqual(seqs, result)
        # assert the returned alignment has the correct parameter values in the
        # align.Info object.
        if param_vals:
            for param, val in param_vals:
                self.assertEqual(aln.Info.AlignParams[param], val)
    
    def test_progressive1(self):
        """test progressive alignment, gaps in middle"""
        self._test_aln({
                'A': 'tacagta', 
                'B': 'tac-gtc',
                'C': 'ta---ta', 
                'D': 'tac-gtc',
                })
         
    def test_progressive2(self):
        """test progressive alignment, gaps in middle"""
        self._test_aln({
                'A': 'ac-ttgt', 
                'B': 'ac---gt',
                'C': 'aca--gt', 
                'D': 'ac---gt',
                })
    
    def test_progressive_params(self):
        """excercise progressive alignment providing model params"""
        self._test_aln({
                'A': 'tacagta', 
                'B': 'tac-gtc',
                'C': 'ta---ta', 
                'D': 'cac-cta',
                }, model=HKY85(), param_vals=[('kappa',2.0)])
    
    def test_TreeAlign_does_pairs(self):
        """test TreeAlign handles pairs of sequences"""
        self._test_aln({
                'A': 'acttgtac', 
                'B': 'ac--gtac',
                })
    
    def test_gap_at_start(self):
        """test progressive alignment, gaps at start"""
        self._test_aln({
                'A': '-ac', 
                'B': '-ac',
                'C': '-ac', 
                'D': 'gac',
                })
    
    def test_gap_at_end(self):
        """test progressive alignment, gaps at end"""
        self._test_aln({
                'A': 'gt-', 
                'B': 'gt-',
                'C': 'gt-', 
                'D': 'gta',
                })
    
    def test_gaps2(self):
        """Gaps have real costs, even end gaps"""
        self._test_aln({
                'A': 'g-', 
                'B': 'g-',
                'C': 'ga', 
                'D': 'a-',
                })
    
        self._test_aln({
                'A': '-g', 
                'B': '-g',
                'C': 'ag', 
                'D': '-a',
                })

    def test_difficult_end_gaps(self):
        self._test_aln({
                'A': '--cctc', 
                'B': '--cctc',
                'C': 'gacctc', 
                'D': 'ga----',
                })  
        return  

        self._test_aln({
                'A': 'gcctcgg------', 
                'B': 'gcctcgg------',
                'C': 'gcctcggaaacgt', 
                'D': '-------aaacgt',
                })    



class HirschbergTestCase(MultipleAlignmentTestCase):
    # Force use of linear space algorithm
    
    def _test_aln(self, seqs, **kw):
        tmp = cogent.align.pairwise.HIRSCHBERG_LIMIT
        try:
            cogent.align.pairwise.HIRSCHBERG_LIMIT = 100
            result = MultipleAlignmentTestCase._test_aln(self, seqs, **kw)
        finally:
            cogent.align.pairwise.HIRSCHBERG_LIMIT = tmp
        return result


    
if __name__ == '__main__':
    unittest.main()

