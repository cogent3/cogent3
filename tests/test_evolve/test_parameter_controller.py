import os
import warnings

from unittest import TestCase, main

import cogent3.evolve.parameter_controller
import cogent3.evolve.substitution_model

from cogent3 import make_aligned_seqs, make_tree


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose, assert_almost_equal


base_path = os.getcwd()
data_path = os.path.join(base_path, "data")

good_rule_sets = [
    [{"par_name": "length", "is_independent": True}],
    [{"par_name": "length", "is_independent": True}],
    [
        {
            "par_name": "length",
            "clade": True,
            "is_independent": True,
            "edges": ["a", "b"],
        }
    ],
    [{"par_name": "length", "is_independent": True, "edges": ["a", "c", "e"]}],
    [{"par_name": "length", "is_independent": True, "edge": "a"}],
]
bad_rule_sets = [[{"par_name": "length", "clade": True, "edges": ["b", "f"]}]]


class test_parameter_controller(TestCase):
    """Tesing Parameter Controller"""

    def setUp(self):
        # length all edges 1 except c=2.  b&d transitions all other
        # transverions
        self.al = make_aligned_seqs(
            data={"a": "tata", "b": "tgtc", "c": "gcga", "d": "gaac", "e": "gagc"}
        )
        self.tree = make_tree(treestring="((a,b),(c,d),e);")
        self.model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
            equal_motif_probs=True, model_gaps=True
        )

    def test_scoped_local(self):
        model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
            equal_motif_probs=True, model_gaps=True, predicates={"kappa": "transition"}
        )
        lf = model.make_likelihood_function(self.tree)
        lf.set_constant_lengths()
        lf.set_alignment(self.al)
        null = lf.get_num_free_params()
        lf.set_param_rule(par_name="kappa", is_independent=True, edges=["b", "d"])
        self.assertEqual(null + 2, lf.get_num_free_params())

    def test_set_get_motif_probs_nstat(self):
        from cogent3 import get_model

        aln = make_aligned_seqs(
            data=dict(
                a="AACGAAGCAGAGTCACGGCA",
                b="ACGGAAGTTGAGTCACCCCA",
                c="TGCATCGAAAAGTCACGCTG",
            ),
            moltype="dna",
        )
        bases = "ACGT"
        expect = aln.get_motif_probs()
        expect = [expect[b] for b in bases]
        tree = make_tree("(a,b,c)")
        gn = get_model("GN")
        lf = gn.make_likelihood_function(tree)
        lf.set_alignment(aln)
        got = lf.get_motif_probs().to_dict()
        got = [got[b] for b in bases]
        assert_allclose(got, expect)

    def test_set_motif_probs(self):
        """Mprobs supplied to the parameter controller"""

        def compare_mprobs(got, exp):
            # handle min val
            motifs = list(got)
            assert_almost_equal(
                [got[m] for m in motifs], [exp[m] for m in motifs], decimal=5
            )

        model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
            model_gaps=True, motif_probs=None
        )
        lf = model.make_likelihood_function(self.tree, motif_probs_from_align=False)

        mprobs = {"A": 0.1, "C": 0.2, "G": 0.2, "T": 0.5, "-": 0.0}
        lf.set_motif_probs(mprobs)
        # node the LF adjust motif probs so they are all >= 1e-6
        got = lf.get_motif_probs().to_dict()
        compare_mprobs(got, mprobs)

        lf.set_motif_probs_from_data(self.al[:1], is_constant=True)
        assert_almost_equal(lf.get_motif_probs()["G"], 0.6, decimal=4)

        lf.set_motif_probs_from_data(self.al[:1], pseudocount=1)
        self.assertNotEqual(lf.get_motif_probs()["G"], 0.6)

        # test with consideration of ambiguous states
        al = make_aligned_seqs(
            data={"seq1": "ACGTAAGNA", "seq2": "ACGTANGTC", "seq3": "ACGTACGTG"}
        )
        lf.set_motif_probs_from_data(al, include_ambiguity=True, is_constant=True)
        motif_probs = dict(lf.get_motif_probs())
        correct_probs = {
            "A": 8.5 / 27,
            "C": 5.5 / 27,
            "-": 0.0,
            "T": 5.5 / 27,
            "G": 7.5 / 27,
        }
        compare_mprobs(motif_probs, correct_probs)
        assert_allclose(sum(motif_probs.values()), 1.0)

    def test_set_multilocus(self):
        """2 loci each with own mprobs"""
        model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
            motif_probs=None
        )
        lf = model.make_likelihood_function(
            self.tree, motif_probs_from_align=False, loci=["a", "b"]
        )

        mprobs_a = dict(A=0.2, T=0.2, C=0.3, G=0.3)
        mprobs_b = dict(A=0.1, T=0.2, C=0.3, G=0.4)

        for is_constant in [False, True]:
            lf.set_motif_probs(mprobs_a, is_constant=is_constant)
            lf.set_motif_probs(mprobs_b, locus="b")
            self.assertEqual(lf.get_motif_probs(locus="a"), mprobs_a)
            self.assertEqual(lf.get_motif_probs(locus="b"), mprobs_b)

    def test_set_param_rules(self):
        lf = self.model.make_likelihood_function(self.tree)

        def do_rules(rule_set):
            for rule in rule_set:
                lf.set_param_rule(**rule)

        for rule_set in good_rule_sets:
            lf.set_default_param_rules()
            do_rules(rule_set)
        for rule_set in bad_rule_sets:
            lf.set_default_param_rules()
            self.assertRaises(
                (KeyError, TypeError, AssertionError, ValueError), do_rules, rule_set
            )

    def test_set_constant_lengths(self):
        t = make_tree(treestring="((a:1,b:2):3,(c:4,d:5):6,e:7);")
        lf = self.model.make_likelihood_function(t)  # self.tree)
        lf.set_param_rule("length", is_constant=True)
        # lf.set_constant_lengths(t)
        lf.set_alignment(self.al)
        self.assertEqual(lf.get_param_value("length", "b"), 2)
        self.assertEqual(lf.get_param_value("length", "d"), 5)

    def test_pairwise_clock(self):
        al = make_aligned_seqs(data={"a": "agct", "b": "ggct"})
        tree = make_tree(treestring="(a,b);")
        model = cogent3.evolve.substitution_model.TimeReversibleDinucleotide(
            equal_motif_probs=True, model_gaps=True, mprob_model="tuple"
        )
        lf = model.make_likelihood_function(tree)
        lf.set_local_clock("a", "b")
        lf.set_alignment(al)
        lf.optimise(local=True)
        rd = lf.get_param_value_dict(["edge"], params=["length"])
        self.assertAlmostEqual(lf.get_log_likelihood(), -10.1774488956)
        self.assertEqual(rd["length"]["a"], rd["length"]["b"])

    def test_local_clock(self):
        lf = self.model.make_likelihood_function(self.tree)
        lf.set_local_clock("c", "d")
        lf.set_alignment(self.al)
        lf.optimise(local=True, tolerance=1e-8, max_restarts=2)
        rd = lf.get_param_value_dict(["edge"], params=["length"])
        self.assertAlmostEqual(lf.get_log_likelihood(), -27.84254174)
        self.assertEqual(rd["length"]["c"], rd["length"]["d"])
        self.assertNotEqual(rd["length"]["a"], rd["length"]["e"])

    def test_complex_parameter_rules(self):
        # This test has many local minima and so does not cope
        # with changes to optimiser details.
        model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
            equal_motif_probs=True, model_gaps=True, predicates={"kappa": "transition"}
        )
        lf = model.make_likelihood_function(self.tree)
        lf.set_param_rule(par_name="kappa", is_independent=True)
        lf.set_param_rule(par_name="kappa", is_independent=False, edges=["b", "d"])
        lf.set_constant_lengths(make_tree(treestring="((a:1,b:1):1,(c:2,d:1):1,e:1);"))
        # print self.pc
        lf.set_alignment(self.al)
        lf.optimise(local=True)
        rd = lf.get_param_value_dict(["edge"], params=["kappa"])
        self.assertAlmostEqual(lf.get_log_likelihood(), -27.3252, 3)
        self.assertEqual(rd["kappa"]["b"], rd["kappa"]["d"])
        self.assertNotEqual(rd["kappa"]["a"], rd["kappa"]["b"])

    def test_bounds(self):
        """Test setting upper and lower bounds for parameters"""
        lf = self.model.make_likelihood_function(self.tree)
        lf.set_param_rule("length", value=3, lower=0, upper=5)

        # Out of bounds value should warn and keep bounded
        with warnings.catch_warnings(record=True) as w:
            lf.set_param_rule("length", lower=0, upper=2)
            self.assertTrue(len(w), "No warning issued")
        self.assertEqual(lf.get_param_value("length", edge="a"), 2)

        # upper < lower bounds should fail
        self.assertRaises(ValueError, lf.set_param_rule, "length", lower=2, upper=0)


if __name__ == "__main__":
    main()
