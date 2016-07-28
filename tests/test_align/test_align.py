#!/usr/bin/env python

from cogent3 import DNA, LoadSeqs
from cogent3.align.align import classic_align_pairwise, make_dna_scoring_dict,\
    local_pairwise, global_pairwise
from cogent3.evolve.models import HKY85
import cogent3.evolve.substitution_model
dna_model = cogent3.evolve.substitution_model.Nucleotide(
    model_gaps=False, equal_motif_probs=True)

import cogent3.align.progressive

import unittest

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "3.0.alpha"
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

seq1 = DNA.make_sequence('aaaccggacattacgtgcgta', name='FAKE01')
seq2 = DNA.make_sequence('ccggtcaggttacgtacgtt', name='FAKE02')


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
        s1 = DNA.make_sequence(s[:-2], name="A")
        s2 = DNA.make_sequence(s[2:], name="B")
        for a in self._aligned_both_ways(s1, s2, local=False):
            self.assertEqual(matchedColumns(a), 6)
            self.assertEqual(len(a), 10)

    def test_short(self):
        s1 = DNA.make_sequence('tacagta', name="A")
        s2 = DNA.make_sequence('tacgtc', name="B")
        for a in self._aligned_both_ways(s1, s2, local=False):
            self.assertEqual(matchedColumns(a), 5)
            self.assertEqual(len(a), 7)

    def test_pairwise_returns_score(self):
        """exercise pairwise local/global returns alignment score"""
        S = make_dna_scoring_dict(10, -1, -8)
        aln, score = local_pairwise(seq1, seq2, S, 10, 2, return_score=True)
        self.assertTrue(score > 100)
        aln, score = global_pairwise(seq1, seq2, S, 10, 2, return_score=True)
        self.assertTrue(score > 100)

    def test_codon(self):
        s1 = DNA.make_sequence('tacgccgta', name="A")
        s2 = DNA.make_sequence('tacgta', name="B")
        codon_model = cogent3.evolve.substitution_model.Codon(
            model_gaps=False, equal_motif_probs=True,
            mprob_model='conditional')
        tree = cogent3.LoadTree(tip_names=['A', 'B'])
        lf = codon_model.make_likelihood_function(tree, aligned=False)
        lf.set_sequences(dict(A=s1, B=s2))
        a = lf.get_log_likelihood().edge.get_viterbi_path().get_alignment()
        self.assertEqual(matchedColumns(a), 6)
        self.assertEqual(len(a), 9)

    def test_local_tiebreak(self):
        """Should pick the first best-equal hit rather than the last one"""
        # so that the Pyrex and Python versions give the same result.
        score_matrix = make_dna_scoring_dict(match=1, transition=-1,
                                             transversion=-1)
        pattern = DNA.make_sequence('cwc', name='pattern')
        two_hit = DNA.make_sequence('cactc', name='target')
        aln = local_pairwise(pattern, two_hit, score_matrix, 5, 2)
        hit = aln.named_seqs['target']
        self.assertEqual(str(hit).lower(), 'cac')


class UnalignedPairTestCase(unittest.TestCase):

    def test_forward(self):
        tree = cogent3.LoadTree(tip_names='AB')
        pc = dna_model.make_likelihood_function(tree, aligned=False)
        pc.set_sequences({'A': seq1, 'B': seq2})
        LnL = pc.get_log_likelihood()
        assert isinstance(LnL, float)


class MultipleAlignmentTestCase(unittest.TestCase):

    def _make_aln(self, orig, model=dna_model, param_vals=None,
                  indel_rate=0.1, indel_length=0.5, **kw):
        kw['indel_rate'] = indel_rate
        kw['indel_length'] = indel_length
        seqs = dict((key, DNA.make_sequence(value))
                    for (key, value) in list(orig.items()))
        if len(seqs) == 2:
            tree = cogent3.LoadTree(tip_names=list(seqs.keys()))
            tree = cogent3.LoadTree(treestring="(A:.1,B:.1)")
        else:
            tree = cogent3.LoadTree(
                treestring="(((A:.1,B:.1):.1,C:.1):.1,D:.1)")
        aln, tree = cogent3.align.progressive.TreeAlign(model, seqs,
                                                        tree=tree, param_vals=param_vals, show_progress=False, **kw)
        return aln

    def _test_aln(self, seqs, model=dna_model, param_vals=None, **kw):
        orig = dict((n, s.replace('-', '')) for (n, s) in list(seqs.items()))
        aln = self._make_aln(orig, model=model, param_vals=param_vals, **kw)
        result = dict((n, s.lower()) for (n, s) in list(aln.todict().items()))
        # assert the alignment result is correct
        self.assertEqual(seqs, result)
        # assert the returned alignment has the correct parameter values in the
        # align.info object.
        if param_vals:
            for param, val in param_vals:
                self.assertEqual(aln.info.AlignParams[param], val)

    def test_progressive1(self):
        """test progressive alignment, gaps in middle"""
        self._test_aln({
            'A': 'tacagta',
            'B': 'tac-gtc',
            'C': 'ta---ta',
            'D': 'tac-gtc',
            })

    def test_progressive_est_tree(self):
        """excercise progressive alignment without a guide tree"""
        seqs = LoadSeqs(data={'A': "TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA",
                              'B': "TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT",
                              'C': "TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT"}, aligned=False)
        aln, tree = cogent3.align.progressive.TreeAlign(HKY85(), seqs, show_progress=False,
                                                        param_vals={'kappa': 4.0})

        expect = {'A': 'TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA-------',
                  'C': 'TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT',
                  'B': 'TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT'}
        self.assertEqual(aln.todict(), expect)

    def test_progressive_params(self):
        """excercise progressive alignment providing model params"""
        self._test_aln({
            'A': 'tacagta',
            'B': 'tac-gtc',
            'C': 'ta---ta',
            'D': 'cac-cta',
            }, model=HKY85(), param_vals=[('kappa', 2.0)])

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
        """gaps have real costs, even end gaps"""
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
        tmp = cogent3.align.pairwise.HIRSCHBERG_LIMIT
        try:
            cogent3.align.pairwise.HIRSCHBERG_LIMIT = 100
            result = MultipleAlignmentTestCase._test_aln(self, seqs, **kw)
        finally:
            cogent3.align.pairwise.HIRSCHBERG_LIMIT = tmp
        return result


if __name__ == '__main__':
    unittest.main()
