import itertools
import pathlib
import unittest

import numpy
import pytest

import cogent3.align.progressive
import cogent3.evolve.substitution_model

from cogent3 import DNA, load_aligned_seqs, load_tree, make_unaligned_seqs
from cogent3.align.align import (
    classic_align_pairwise,
    global_pairwise,
    local_pairwise,
    make_dna_scoring_dict,
    make_generic_scoring_dict,
)
from cogent3.evolve.models import HKY85, get_model


DATADIR = pathlib.Path(__file__).parent.parent / "data"

dna_model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
    model_gaps=False, equal_motif_probs=True
)


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


seq1 = DNA.make_seq("aaaccggacattacgtgcgta", name="FAKE01")
seq2 = DNA.make_seq("ccggtcaggttacgtacgtt", name="FAKE02")


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

    def test_pwise_protein(self):
        """works for pairwise protein alignment"""
        from cogent3 import PROTEIN

        S = make_generic_scoring_dict(1, PROTEIN)
        seq1 = PROTEIN.make_seq("MAYPFQLGLQD", "seq1")
        seq2 = PROTEIN.make_seq("MAYPFGLQD", "seq2")
        a1 = classic_align_pairwise(seq1, seq2, S, 10, 2, local=False)
        self.assertEqual(a1.to_dict(), dict(seq1="MAYPFQLGLQD", seq2="MAYPF--GLQD"))

    def test_gap_at_one_end(self):
        for a in self._aligned_both_ways(seq1, seq2, local=False):
            self.assertEqual(matchedColumns(a), 15)
            self.assertEqual(len(a), 23)

    def test_gaps_at_both_ends(self):
        s = "aaaccggttt"
        s1 = DNA.make_seq(s[:-2], name="A")
        s2 = DNA.make_seq(s[2:], name="B")
        for a in self._aligned_both_ways(s1, s2, local=False):
            self.assertEqual(matchedColumns(a), 6)
            self.assertEqual(len(a), 10)

    def test_short(self):
        s1 = DNA.make_seq("tacagta", name="A")
        s2 = DNA.make_seq("tacgtc", name="B")
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
        s1 = DNA.make_seq("tacgccgta", name="A")
        s2 = DNA.make_seq("tacgta", name="B")
        codon_model = cogent3.evolve.substitution_model.TimeReversibleCodon(
            model_gaps=False, equal_motif_probs=True, mprob_model="conditional"
        )
        tree = cogent3.make_tree(tip_names=["A", "B"])
        lf = codon_model.make_likelihood_function(tree, aligned=False)
        lf.set_sequences(dict(A=s1, B=s2))
        a = lf.get_log_likelihood().edge.get_viterbi_path().get_alignment()
        self.assertEqual(matchedColumns(a), 6)
        self.assertEqual(len(a), 9)

    def test_local_tiebreak(self):
        """Should pick the first best-equal hit rather than the last one"""
        # so that the Pyrex and Python versions give the same result.
        score_matrix = make_dna_scoring_dict(match=1, transition=-1, transversion=-1)
        pattern = DNA.make_seq("cwc", name="pattern")
        two_hit = DNA.make_seq("cactc", name="target")
        aln = local_pairwise(pattern, two_hit, score_matrix, 5, 2)
        hit = aln.named_seqs["target"]
        self.assertEqual(str(hit).lower(), "cac")


class UnalignedPairTestCase(unittest.TestCase):
    def test_forward(self):
        tree = cogent3.make_tree(tip_names="AB")
        pc = dna_model.make_likelihood_function(tree, aligned=False)
        pc.set_sequences({"A": seq1, "B": seq2})
        LnL = pc.get_log_likelihood()
        assert isinstance(LnL, float)


class MultipleAlignmentTestCase(unittest.TestCase):
    def _make_aln(
        self,
        orig,
        model=dna_model,
        param_vals=None,
        indel_rate=0.1,
        indel_length=0.5,
        **kw,
    ):
        kw["indel_rate"] = indel_rate
        kw["indel_length"] = indel_length
        seqs = {key: DNA.make_seq(value) for (key, value) in list(orig.items())}
        if len(seqs) == 2:
            tree = cogent3.make_tree(treestring="(A:.1,B:.1)")
        else:
            tree = cogent3.make_tree(treestring="(((A:.1,B:.1):.1,C:.1):.1,D:.1)")
        aln, tree = cogent3.align.progressive.tree_align(
            model, seqs, tree=tree, param_vals=param_vals, show_progress=False, **kw
        )
        return aln

    def _test_aln(self, seqs, model=dna_model, param_vals=None, **kw):
        orig = {n: s.replace("-", "") for (n, s) in list(seqs.items())}
        aln = self._make_aln(orig, model=model, param_vals=param_vals, **kw)
        result = {n: s.lower() for (n, s) in list(aln.to_dict().items())}
        # assert the alignment result is correct
        self.assertEqual(seqs, result)
        # and the moltype matches the model
        model = get_model(model)
        self.assertIs(aln.moltype, model.moltype)

        # assert the returned alignment has the correct parameter values in the
        # align.info object.
        if param_vals:
            for param, val in param_vals:
                self.assertEqual(aln.info.align_params[param], val)

    def test_progressive1(self):
        """test progressive alignment, gaps in middle"""
        self._test_aln(
            {"A": "tacagta", "B": "tac-gtc", "C": "ta---ta", "D": "tac-gtc"},
            model="F81",
        )

    def test_progessive_model_name(self):
        """tree_align handles models specified by name"""
        self._test_aln({"A": "tacagta", "B": "tac-gtc", "C": "ta---ta", "D": "tac-gtc"})

    def test_progressive_est_tree(self):
        """exercise progressive alignment without a guide tree"""
        seqs = make_unaligned_seqs(
            data={
                "A": "TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA",
                "B": "TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT",
                "C": "TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT",
            }
        )
        aln, tree = cogent3.align.progressive.tree_align(
            HKY85(), seqs, show_progress=False, param_vals={"kappa": 4.0}
        )

        expect = {
            "A": "TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA-------",
            "C": "TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT",
            "B": "TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT",
        }
        self.assertEqual(aln.to_dict(), expect)

        aln, tree = cogent3.align.progressive.tree_align(
            HKY85(),
            seqs,
            show_progress=False,
            params_from_pairwise=True,
        )

        expect = {
            "A": "TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA-------",
            "C": "TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT",
            "B": "TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT",
        }
        self.assertEqual(aln.to_dict(), expect)

    def test_align_info(self):
        """alignment info object has parameter values"""
        aln = self._make_aln(
            {"A": "gcctcgg", "B": "gcctcgg", "C": "gcctcggaaacgt", "D": "aaacgt"}
        )
        self.assertTrue(aln.info["align_params"]["lnL"] < 0)

    def test_progressive_params(self):
        """excercise progressive alignment providing model params"""
        self._test_aln(
            {"A": "tacagta", "B": "tac-gtc", "C": "ta---ta", "D": "cac-cta"},
            model=HKY85(),
            param_vals=[("kappa", 2.0)],
        )

    def test_tree_align_does_pairs(self):
        """test tree_align handles pairs of sequences"""
        self._test_aln({"A": "acttgtac", "B": "ac--gtac"})

    def test_gap_at_start(self):
        """test progressive alignment, gaps at start"""
        self._test_aln({"A": "-ac", "B": "-ac", "C": "-ac", "D": "gac"})

    def test_gap_at_end(self):
        """test progressive alignment, gaps at end"""
        self._test_aln({"A": "gt-", "B": "gt-", "C": "gt-", "D": "gta"})

    def test_gaps2(self):
        """gaps have real costs, even end gaps"""
        self._test_aln({"A": "g-", "B": "g-", "C": "ga", "D": "a-"})

        self._test_aln({"A": "-g", "B": "-g", "C": "ag", "D": "-a"})

    def test_difficult_end_gaps(self):
        self._test_aln({"A": "--cctc", "B": "--cctc", "C": "gacctc", "D": "ga----"})
        self._test_aln(
            {
                "A": "gcctcgg------",
                "B": "gcctcgg------",
                "C": "gcctcggaaacgt",
                "D": "-------aaacgt",
            }
        )

    def test_tree_align_handles_zero_lengths(self):
        seqs = make_unaligned_seqs(
            data={
                "A": "TTAATTTTAGTAGTGCTATCCCC",
                "B": "TTAATTTTAGTAGTGCTATCCCA",
                "C": "TTAATTTTAGTAGTGCTATCC",
            },
            moltype="dna",
        )

        expected = {
            "A": "TTAATTTTAGTAGTGCTATCCCC",
            "B": "TTAATTTTAGTAGTGCTATCCCA",
            "C": "TTAATTTTAGTAGTGCTATCC--",
        }

        tree_mapping = [
            "A: 0.0225400070648391",
            "B: 0.0225400070648391",
            "C: 0.0",
        ]
        tree_variants = itertools.permutations(tree_mapping, r=3)

        for tree_encoding in tree_variants:
            aln, _ = cogent3.align.progressive.tree_align(
                model="F81",
                seqs=seqs,
                tree=cogent3.make_tree("({},{},{})".format(*tree_encoding)),
                show_progress=False,
            )
            self.assertEqual(aln.to_dict(), expected)


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


@pytest.fixture(scope="session")
def seqs():
    tree = load_tree(DATADIR / "brca1_5.tree")
    aln = load_aligned_seqs(DATADIR / "brca1.fasta", moltype="dna")
    seqs = aln[200:1200].take_seqs(tree.get_tip_names()).degap()
    return seqs


def test_tree_align_pwise_iter(seqs):
    kwargs = dict(
        model="F81", seqs=seqs, show_progress=False, indel_rate=1e-3, indel_length=1e-1
    )
    aln, tree = cogent3.align.progressive.tree_align(iters=None, **kwargs)
    one = aln.alignment_quality(app_name="sp_score", calc="pdist")
    for _ in range(10):
        aln, tree = cogent3.align.progressive.tree_align(
            iters=1, approx_dists=True, **kwargs
        )
        two = aln.alignment_quality(app_name="sp_score", calc="pdist")
        # the quality scores will differ, but they're not deterministic
        # because the alignments are not deterministic
        if not numpy.allclose(two, one):
            break
    else:
        raise AssertionError("all attempts produced alignments with identical quality")


def test_tree_align_dists_from_pairwise_align(seqs):
    # difficult to test reliably so only exercising use of option
    aln, tree = cogent3.align.progressive.tree_align(
        model="F81", seqs=seqs, show_progress=False, approx_dists=False
    )
    assert aln


def test_tree_align_two(seqs):
    seqs = seqs.take_seqs(["Human", "NineBande"])
    aln, tree = cogent3.align.progressive.tree_align(
        model="F81", seqs=seqs, show_progress=False, iters=1, approx_dists=True
    )
    # the tree should have equal branch lengths
    dist = set(tree.get_distances().values())
    assert len(dist) == 1
    assert isinstance(list(dist)[0], float)
    assert len(aln) >= seqs.get_lengths().array.min()
