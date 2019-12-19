from unittest import TestCase

from cogent3 import load_aligned_seqs, load_tree
from cogent3.app.evo import model


class LikelihoodTree2Tests(TestCase):
    def test_likelihood_calc_consistency(self):
        aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
        aln = aln.take_seqs(aln.names[:5])
        tree = load_tree("data/murphy.tree", "murphy.tree")
        tree = tree.get_sub_tree(aln.names)

        codon = model(
            "GNC",
            tree=tree,
            show_progress=True,
            opt_args=dict(max_evaluations=1000, limit_action="ignore"),
        )
        result = codon(aln)
        self.assertEqual("%.4f" % result.lnL, "-8022.5806")
        self.assertEqual(result.nfp, 79)
        self.assertTrue(result.DLC)
