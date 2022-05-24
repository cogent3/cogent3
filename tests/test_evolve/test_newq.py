#!/usr/bin/env python

import warnings

from unittest import TestCase, main

from numpy import dot, ones
from numpy.testing import assert_allclose

from cogent3 import (
    DNA,
    load_aligned_seqs,
    load_tree,
    make_aligned_seqs,
    make_tree,
)
from cogent3.evolve.ns_substitution_model import (
    NonReversibleCodon,
    NonReversibleNucleotide,
)
from cogent3.evolve.predicate import MotifChange
from cogent3.evolve.substitution_model import (
    TimeReversibleCodon,
    TimeReversibleNucleotide,
)
from cogent3.maths.matrix_exponentiation import PadeExponentiator as expm


warnings.filterwarnings("ignore", "Motif probs overspecified")
warnings.filterwarnings("ignore", "Model not reversible")


__author__ = "Peter Maxwell and  Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def _dinuc_root_probs(x, y=None):
    if y is None:
        y = x
    return dict(
        [(n1 + n2, p1 * p2) for n1, p1 in list(x.items()) for n2, p2 in list(y.items())]
    )


def _trinuc_root_probs(x, y, z):
    return dict(
        [
            (n1 + n2 + n3, p1 * p2 * p3)
            for n1, p1 in list(x.items())
            for n2, p2 in list(y.items())
            for n3, p3 in list(z.items())
        ]
    )


def make_p(length, coord, val):
    """returns a probability matrix with value set at coordinate in
    instantaneous rate matrix"""
    Q = ones((4, 4), float) * 0.25  # assumes equi-frequent mprobs at root
    for i in range(4):
        Q[i, i] = 0.0
    Q[coord] *= val
    row_sum = Q.sum(axis=1)
    scale = 1 / (0.25 * row_sum).sum()
    for i in range(4):
        Q[i, i] -= row_sum[i]
    Q *= scale
    return expm(Q)(length)


class NewQ(TestCase):
    aln = make_aligned_seqs(
        data={
            "seq1": "TGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
            "seq2": "TGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
        },
        moltype=DNA,
    )
    tree = make_tree(tip_names=["seq1", "seq2"])

    symm_nuc_probs = dict(A=0.25, T=0.25, C=0.25, G=0.25)
    symm_root_probs = _dinuc_root_probs(symm_nuc_probs)
    asymm_nuc_probs = dict(A=0.1, T=0.1, C=0.4, G=0.4)
    asymm_root_probs = _dinuc_root_probs(asymm_nuc_probs)
    posn_root_probs = _dinuc_root_probs(symm_nuc_probs, asymm_nuc_probs)
    cond_root_probs = dict(
        [
            (n1 + n2, p1 * [0.1, 0.7][n1 == n2])
            for n1, p1 in list(asymm_nuc_probs.items())
            for n2 in "ATCG"
        ]
    )

    # Each of these (data, model) pairs should give a result different
    # from any of the simpler models applied to the same data.
    ordered_by_complexity = [
        # P(AA) == P(GG) == P(AG)
        [symm_root_probs, "tuple"],
        # P(GA) == P(AG) but P(AA) != P(GG)
        [asymm_root_probs, "monomer"],
        # P(AG) == P(A?)*P(?G) but P(A?) != P(?A)
        [posn_root_probs, "monomers"],
        # P(AG) != P(A?)*P(?G)
        [cond_root_probs, "conditional"],
    ]

    def test_newQ_is_nuc_process(self):
        """newQ is an extension of an independent nucleotide process"""
        nuc = TimeReversibleNucleotide(motif_probs=self.asymm_nuc_probs)
        new_di = TimeReversibleNucleotide(
            motif_length=2, mprob_model="monomer", motif_probs=self.asymm_root_probs
        )

        nuc_lf = nuc.make_likelihood_function(self.tree)
        new_di_lf = new_di.make_likelihood_function(self.tree)
        # newQ branch length is exactly motif_length*nuc branch length
        nuc_lf.set_param_rule("length", is_independent=False, init=0.2)
        new_di_lf.set_param_rule("length", is_independent=False, init=0.4)

        nuc_lf.set_alignment(self.aln)
        new_di_lf.set_alignment(self.aln)
        assert_allclose(nuc_lf.get_log_likelihood(), new_di_lf.get_log_likelihood())

    def test_lf_display(self):
        """str of likelihood functions should not fail"""
        for (dummy, model) in self.ordered_by_complexity:
            di = TimeReversibleNucleotide(motif_length=2, mprob_model=model)
            di.adapt_motif_probs(self.cond_root_probs, auto=True)
            lf = di.make_likelihood_function(self.tree)
            str(lf)

    def test_get_statistics(self):
        """get statistics should correctly apply arguments"""
        for (mprobs, model) in self.ordered_by_complexity:
            di = TimeReversibleNucleotide(
                motif_length=2, motif_probs=mprobs, mprob_model=model
            )
            lf = di.make_likelihood_function(self.tree)
            for wm, wt in [(True, True), (True, False), (False, True), (False, False)]:
                stats = lf.get_statistics(with_motif_probs=wm, with_titles=wt)

    def test_get_statistics_mprobs(self):
        """get_statistics motif probs table has motifs as title"""
        sm = NonReversibleCodon()
        lf = sm.make_likelihood_function(self.tree)
        stats = lf.get_statistics(with_motif_probs=True, with_titles=True)
        mprobs = stats[-1]
        self.assertEqual(set(mprobs.header), set(sm.get_motifs()))

    def test_get_motif_probs(self):
        """exercise getting motif probs under all models"""
        for (mprobs, model) in self.ordered_by_complexity:
            di = TimeReversibleNucleotide(
                motif_length=2, motif_probs=mprobs, mprob_model=model
            )
            lf = di.make_likelihood_function(self.tree)
            lf.set_alignment(self.aln)
            if model == "monomers":
                _ = lf.get_motif_probs(position=0)

    def test_sim_alignment(self):
        """should be able to simulate an alignment under all models"""
        for (mprobs, model) in self.ordered_by_complexity:
            di = TimeReversibleNucleotide(
                motif_length=2, motif_probs=mprobs, mprob_model=model
            )
            lf = di.make_likelihood_function(self.tree)
            lf.set_param_rule("length", is_independent=False, init=0.4)
            lf.set_alignment(self.aln)
            lf.simulate_alignment()

    def test_reconstruct_ancestor(self):
        """should be able to reconstruct ancestral sequences under all
        models"""
        for (mprobs, model) in self.ordered_by_complexity:
            di = TimeReversibleNucleotide(motif_length=2, mprob_model=model)
            di.adapt_motif_probs(mprobs, auto=True)
            lf = di.make_likelihood_function(self.tree)
            lf.set_param_rule("length", is_independent=False, init=0.4)
            lf.set_alignment(self.aln)
            lf.reconstruct_ancestral_seqs()

    def test_results_different(self):
        for (i, (mprobs, dummy)) in enumerate(self.ordered_by_complexity):
            results = []
            for (dummy, model) in self.ordered_by_complexity:
                di = TimeReversibleNucleotide(
                    motif_length=2, motif_probs=mprobs, mprob_model=model
                )
                lf = di.make_likelihood_function(self.tree)
                lf.set_param_rule("length", is_independent=False, init=0.4)
                lf.set_alignment(self.aln)
                lh = lf.get_log_likelihood()
                for other in results[:i]:
                    self.assertNotAlmostEqual(other, lh, places=2)
                for other in results[i:]:
                    assert_allclose(other, lh)
                results.append(lh)

    def test_position_specific_mprobs(self):
        """correctly compute likelihood when positions have distinct
        probabilities"""
        aln_len = len(self.aln)
        posn1 = []
        posn2 = []
        for name, seq in list(self.aln.to_dict().items()):
            p1 = [seq[i] for i in range(0, aln_len, 2)]
            p2 = [seq[i] for i in range(1, aln_len, 2)]
            posn1.append([name, "".join(p1)])
            posn2.append([name, "".join(p2)])

        # the position specific alignments
        posn1 = make_aligned_seqs(data=posn1)
        posn2 = make_aligned_seqs(data=posn2)

        # a newQ dinucleotide model
        sm = TimeReversibleNucleotide(motif_length=2, mprob_model="monomer")
        lf = sm.make_likelihood_function(self.tree)
        lf.set_alignment(posn1)
        posn1_lnL = lf.get_log_likelihood()
        lf.set_alignment(posn2)
        posn2_lnL = lf.get_log_likelihood()
        expect_lnL = posn1_lnL + posn2_lnL

        # the joint model
        lf.set_alignment(self.aln)
        aln_lnL = lf.get_log_likelihood()

        # setting the full alignment, which has different motif probs, should
        # produce a different lnL
        self.assertNotAlmostEqual(aln_lnL, expect_lnL)

        # set the arguments for taking position specific mprobs
        sm = TimeReversibleNucleotide(motif_length=2, mprob_model="monomers")
        lf = sm.make_likelihood_function(self.tree)
        lf.set_alignment(self.aln)
        lf.get_motif_probs()
        posn12_lnL = lf.get_log_likelihood()
        assert_allclose(posn12_lnL, expect_lnL, rtol=1e-4)

    def test_compute_conditional_mprobs(self):
        """equal likelihood from position specific and conditional mprobs"""

        def compare_models(motif_probs, motif_length):
            # if the 1st and 2nd position motifs are independent of each other
            # then conditional is the same as positional
            ps = TimeReversibleNucleotide(
                motif_length=motif_length,
                motif_probs=motif_probs,
                mprob_model="monomers",
            )
            cd = TimeReversibleNucleotide(
                motif_length=motif_length,
                motif_probs=motif_probs,
                mprob_model="conditional",
            )

            ps_lf = ps.make_likelihood_function(self.tree)
            ps_lf.set_param_rule("length", is_independent=False, init=0.4)
            ps_lf.set_alignment(self.aln)

            cd_lf = cd.make_likelihood_function(self.tree)
            cd_lf.set_param_rule("length", is_independent=False, init=0.4)
            cd_lf.set_alignment(self.aln)
            assert_allclose(cd_lf.get_log_likelihood(), ps_lf.get_log_likelihood())

        compare_models(self.posn_root_probs, 2)
        # trinucleotide
        trinuc_mprobs = _trinuc_root_probs(
            self.asymm_nuc_probs, self.asymm_nuc_probs, self.asymm_nuc_probs
        )
        compare_models(trinuc_mprobs, 3)

    def test_cond_pos_differ(self):
        """lnL should differ when motif probs are not multiplicative"""
        dinuc_probs = {
            "AA": 0.088506666666666664,
            "AC": 0.044746666666666664,
            "GT": 0.056693333333333332,
            "AG": 0.070199999999999999,
            "CC": 0.048653333333333333,
            "TT": 0.10678666666666667,
            "CG": 0.0093600000000000003,
            "GG": 0.049853333333333333,
            "GC": 0.040253333333333335,
            "AT": 0.078880000000000006,
            "GA": 0.058639999999999998,
            "TG": 0.081626666666666667,
            "TA": 0.068573333333333333,
            "CA": 0.06661333333333333,
            "TC": 0.060866666666666666,
            "CT": 0.069746666666666665,
        }

        mg = TimeReversibleNucleotide(
            motif_length=2, motif_probs=dinuc_probs, mprob_model="monomer"
        )
        mg_lf = mg.make_likelihood_function(self.tree)
        mg_lf.set_param_rule("length", is_independent=False, init=0.4)
        mg_lf.set_alignment(self.aln)

        cd = TimeReversibleNucleotide(
            motif_length=2, motif_probs=dinuc_probs, mprob_model="conditional"
        )

        cd_lf = cd.make_likelihood_function(self.tree)
        cd_lf.set_param_rule("length", is_independent=False, init=0.4)
        cd_lf.set_alignment(self.aln)
        self.assertNotAlmostEqual(
            mg_lf.get_log_likelihood(), cd_lf.get_log_likelihood()
        )

    def test_getting_node_mprobs(self):
        """return correct motif probability vector for tree nodes"""
        tree = make_tree(treestring="(a:.2,b:.2,(c:.1,d:.1):.1)")
        aln = make_aligned_seqs(
            data={"a": "TGTG", "b": "TGTG", "c": "TGTG", "d": "TGTG"}
        )

        motifs = ["T", "C", "A", "G"]
        aX = MotifChange(motifs[0], motifs[3], forward_only=True).aliased("aX")
        bX = MotifChange(motifs[3], motifs[0], forward_only=True).aliased("bX")
        edX = MotifChange(motifs[1], motifs[2], forward_only=True).aliased("edX")
        cX = MotifChange(motifs[2], motifs[1], forward_only=True).aliased("cX")
        sm = NonReversibleNucleotide(
            predicates=[aX, bX, edX, cX], equal_motif_probs=True
        )

        lf = sm.make_likelihood_function(tree)
        lf.set_param_rule("aX", edge="a", value=8.0)
        lf.set_param_rule("bX", edge="b", value=8.0)
        lf.set_param_rule("edX", edge="edge.0", value=2.0)
        lf.set_param_rule("cX", edge="c", value=0.5)
        lf.set_param_rule("edX", edge="d", value=4.0)
        lf.set_alignment(aln)

        # we construct the hand calc variants
        mprobs = ones(4, float) * 0.25
        a = make_p(0.2, (0, 3), 8)
        a = dot(mprobs, a)

        b = make_p(0.2, (3, 0), 8)
        b = dot(mprobs, b)

        e = make_p(0.1, (1, 2), 2)
        e = dot(mprobs, e)

        c = make_p(0.1, (2, 1), 0.5)
        c = dot(e, c)

        d = make_p(0.1, (1, 2), 4)
        d = dot(e, d)

        prob_vectors = lf.get_motif_probs_by_node()
        assert_allclose(prob_vectors["a"].array, a)
        assert_allclose(prob_vectors["b"].array, b)
        assert_allclose(prob_vectors["c"].array, c)
        assert_allclose(prob_vectors["d"].array, d)
        assert_allclose(prob_vectors["edge.0"].array, e)

    def test_get_motif_probs_by_node_mg94(self):
        """handles different statespace dimensions from process and stationary distribution"""
        from cogent3.evolve.models import get_model

        aln = load_aligned_seqs("data/primates_brca1.fasta", moltype="dna")
        aln = aln.no_degenerates(motif_length=3)

        tree = load_tree("data/primates_brca1.tree")

        # root mprobs are constant
        sm = get_model("MG94HKY")
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        mprobs = lf.get_motif_probs()

        mprobs = lf.get_motif_probs_by_node()
        self.assertEqual(mprobs.shape, (len(tree.get_edge_vector()), 61))

        # root mprobs are variable
        sm = get_model("MG94HKY", optimise_motif_probs=True)
        sm = get_model("MG94HKY")
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        mprobs = lf.get_motif_probs_by_node()
        self.assertEqual(mprobs.shape, (len(tree.get_edge_vector()), 61))

        # not imlemented for monomers variant
        sm = TimeReversibleCodon(
            mprob_model="monomers", model_gaps=False, recode_gaps=True
        )
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        with self.assertRaises(NotImplementedError):
            _ = lf.get_motif_probs_by_node()


if __name__ == "__main__":
    main()
