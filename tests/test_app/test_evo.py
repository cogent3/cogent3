from unittest import TestCase, main
from unittest.mock import MagicMock

from numpy.testing import assert_allclose, assert_raises

from cogent3 import load_aligned_seqs, make_aligned_seqs, make_tree
from cogent3.app import evo as evo_app
from cogent3.app.result import hypothesis_result
from cogent3.evolve.models import get_model


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.9.13a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestModel(TestCase):
    basedir = "data"

    def test_model_str(self):
        """correct str representation"""
        model = evo_app.model("HKY85", time_het="max")
        got = str(model)
        self.assertEqual(
            got,
            (
                "model(type='model', sm='HKY85', tree=None, "
                "name=None, sm_args=None, lf_args=None, "
                "time_het='max', param_rules=None, "
                "opt_args=None, split_codons=False, "
                "show_progress=False, verbose=False)"
            ),
        )

    def test_model_tree(self):
        """allows tree to be string, None or tree"""
        treestring = "(a,b,c)"
        for tree in (treestring, make_tree(treestring=treestring), None):
            mod = evo_app.model("HKY85", tree=tree)
            expect = None if tree is None else make_tree(treestring=treestring)
            self.assertIsInstance(mod._tree, expect.__class__)

    def test_unique_models(self):
        """hypothesis raises ValueError if models not unique"""
        model1 = evo_app.model("HKY85")
        model2 = evo_app.model("HKY85", time_het="max")
        with self.assertRaises(ValueError):
            hyp = evo_app.hypothesis(model1, model2)

    def test_model_time_het(self):
        """support lf time-het argument edge_sets"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "GN",
            time_het=[dict(edges=["Mouse", "Human"], is_independent=False)],
            opt_args=dict(max_evaluations=25, limit_action="ignore"),
        )
        result = mod(aln)
        # 11 free params per calibrated GN matrix, there are 2
        # 3 params for root motif probs, 3 branch lengths
        expect_nfp = 11 * 2 + 3 + 3
        self.assertEqual(result.lf.nfp, expect_nfp)

    def test_model_param_rules(self):
        """applies upper bound if sensible"""
        mod = evo_app.model(
            "GN",
            param_rules=[dict(par_name="length", edge="Mouse", is_independent=False)],
        )
        self.assertEqual(mod._param_rules[0].get("upper"), 50)
        mod = evo_app.model(
            "GN", param_rules=[dict(par_name="length", edge="Mouse", is_constant=True)]
        )
        self.assertEqual(mod._param_rules[0].get("upper", None), None)

    def test_discrete_time_model(self):
        """works with discrete-time submodel"""
        from cogent3.app.composable import NotCompleted

        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "BH", opt_args=dict(max_evaluations=100, limit_action="ignore")
        )
        r = mod(aln)
        self.assertNotIsInstance(r, NotCompleted)

    def test_model_hypothesis_result_repr(self):
        """result objects __repr__ and _repr_html_ methods work correctly"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        model1 = evo_app.model(
            "F81", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        model2 = evo_app.model(
            "HKY85", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        hyp = evo_app.hypothesis(model1, model2)
        result = hyp(aln)
        self.assertIsInstance(result.__repr__(), str)
        self.assertIsInstance(result._repr_html_(), str)
        self.assertIsInstance(result.null.__repr__(), str)
        self.assertIsInstance(result.null._repr_html_(), str)

    def test_hypothesis_str(self):
        """correct str representation"""
        model1 = evo_app.model("HKY85")
        model2 = evo_app.model("HKY85", name="hky85-max-het", time_het="max")
        hyp = evo_app.hypothesis(model1, model2)
        got = str(hyp)
        expect = (
            "hypothesis(type='hypothesis', null='HKY85', "
            "alternates=(model(type='model', sm='HKY85', tree=None, "
            "name='hky85-max-het', sm_args=None, lf_args=None, "
            "time_het='max', param_rules=None, opt_args=None,"
            " split_codons=False, show_progress=False, verbose=False),),"
            " init_alt=None)"
        )
        self.assertEqual(got, expect)

    def test_split_pos_model(self):
        """model with split codons, access .lf using codon position int"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        tree = make_tree(tip_names=aln.names)
        mod = evo_app.model(
            "F81",
            tree=tree,
            split_codons=True,
            opt_args=dict(max_evaluations=5, limit_action="ignore"),
        )
        result = mod(aln)
        aln1 = result.lf[1].get_param_value("alignment").to_dict()
        self.assertEqual(aln1, aln[::3].to_dict())

    def test_model_summed_branch_lengths(self):
        """returns summed branch lengths"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        model1 = evo_app.model(
            "F81", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        result = model1(aln)
        tree = result.lf.get_annotated_tree()
        assert_allclose(result.total_length(), tree.total_length())
        tree = result.lf.get_annotated_tree(length_as="paralinear")
        assert_allclose(
            result.total_length(length_as="paralinear"), tree.total_length()
        )

    def test_model_result_total_length(self):
        """returns summed branch lengths"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        model1 = evo_app.model(
            "GN", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        result = model1(aln)
        expect_tree = result.lf.get_annotated_tree(length_as="ENS")
        assert_allclose(result.tree.total_length(), expect_tree.total_length())
        # it will be different to the standard length values
        expect_tree = result.lf.get_annotated_tree()
        assert_raises(
            AssertionError,
            assert_allclose,
            result.tree.total_length(),
            expect_tree.total_length(),
        )

    def test_model_result_total_length_split_codon(self):
        """returns summed branch lengths across positions when split_codons True"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        model1 = evo_app.model(
            "GN",
            split_codons=True,
            opt_args=dict(max_evaluations=25, limit_action="ignore"),
        )
        result = model1(aln)
        expect = 0.0
        for lf in result.lf.values():
            tree = lf.get_annotated_tree(length_as="ENS")
            expect += tree.total_length()

        got = result.total_length(length_as="ENS")
        assert_allclose(got, expect)


def _make_getter(val):
    def call(**kwargs):
        return val

    return call


def _make_hyp(aic1, aic2, aic3, nfp1, nfp2, nfp3):
    null = MagicMock()
    null.name = "unrooted"
    null.lf.get_aic = _make_getter(aic1)
    null.nfp = nfp1
    alt1 = MagicMock()
    alt1.name = "alt1"
    alt1.lf.get_aic = _make_getter(aic2)
    alt1.nfp = nfp2
    alt2 = MagicMock()
    alt2.name = "alt2"
    alt2.lf.get_aic = _make_getter(aic3)
    alt2.nfp = nfp3
    hyp = hypothesis_result("unrooted", source="something")
    for m in (null, alt1, alt2):
        hyp[m.name] = m
    return hyp


class TestHypothesisResult(TestCase):
    def test_get_best_model(self):
        """should correctly identify the best model"""
        # aic order is null, alt1, alt2
        # in this case, no substantial diff, so should return smaller nfp, ie null
        hyp = _make_hyp(112, 110, 111, 10, 11, 12)
        got = hyp.get_best_model(threshold=0.05)
        self.assertIs(got, hyp.null)
        # here alt2 is winner
        hyp = _make_hyp(110, 111, 104, 10, 11, 12)
        got = hyp.get_best_model(threshold=0.05)
        self.assertIs(got, hyp["alt2"])
        # but if we set threshold more permissive, it will return null
        got = hyp.get_best_model(threshold=0.03)
        self.assertIs(got, hyp.null)

    def test_select_model(self):
        """correctly identify models"""
        hyp = _make_hyp(112, 110, 111, 10, 11, 12)
        got = set(hyp.select_models(threshold=0.05))
        expect = set(hyp.values())
        self.assertEqual(got, expect)
        # single model
        hyp = _make_hyp(110, 111, 104, 10, 11, 12)
        got = hyp.select_models(threshold=0.05)
        self.assertEqual(len(got), 1)
        self.assertIs(got[0], hyp["alt2"])
        # but if we set threshold more permissive, it will return all
        got = hyp.select_models(threshold=0.03)
        self.assertEqual(len(got), 3)
        expect = set(hyp.values())
        self.assertEqual(set(got), expect)


class TestAncestralStates(TestCase):
    def test_ancestral(self):
        """recon ancestral states works"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "GN", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        anc = evo_app.ancestral_states()
        result = anc(mod(aln))
        self.assertEqual(result["root"].shape, (len(aln), 4))
        assert_allclose(result["root"].row_sum(), 1)


class TestNatSel(TestCase):
    # needs to work for single edge, just two edges, and all combos of clae,
    # stem etc..
    def test_zhang(self):
        """natsel_zhang correctly configured and should not fail"""
        opt = dict(max_evaluations=20, limit_action="ignore")
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
        natsel = evo_app.natsel_zhang(
            "CNFGTR",
            tree="data/primate_brca1.tree",
            tip1="Human",
            tip2="Chimpanzee",
            opt_args=opt,
        )
        result = natsel(aln)
        self.assertEqual(result.df, 3)
        self.assertEqual(result.alt.nfp, 21)
        # the naming scheme is model name followed by null/alt
        self.assertTrue("CNFGTR-null" in result)
        self.assertTrue("CNFGTR-alt" in result)

        # result keys correct when given a model
        Y98 = get_model("Y98")
        natsel = evo_app.natsel_zhang(
            Y98,
            tree="data/primate_brca1.tree",
            tip1="Human",
            tip2="Chimpanzee",
            opt_args=opt,
        )
        result = natsel(aln)
        self.assertEqual(result.df, 3)
        self.assertTrue("Y98-null" in result)
        self.assertTrue("Y98-alt" in result)

        # fails if not a codon model
        with self.assertRaises(ValueError):
            _ = evo_app.natsel_zhang(
                "F81",
                tree="data/primate_brca1.tree",
                tip1="Human",
                tip2="Chimpanzee",
                opt_args=opt,
            )

        # fails if no tip names provided
        with self.assertRaises(ValueError):
            _ = evo_app.natsel_zhang(
                "Y98", tree="data/primate_brca1.tree", opt_args=opt
            )

    def test_zhang_mtseq(self):
        """genetic code setting should work"""
        from cogent3.app.composable import NotCompleted

        opt = dict(max_evaluations=20, limit_action="ignore")
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        natsel = evo_app.natsel_zhang("CNFGTR", tip1="Human", opt_args=opt, gc=2)
        result = natsel(aln)
        self.assertEqual(result.df, 3)
        # but if provide wrong gc, get NotCompleted
        natsel = evo_app.natsel_zhang("CNFGTR", tip1="Human", opt_args=opt, gc=1)
        result = natsel(aln)
        self.assertIsInstance(result, NotCompleted)

    def test_zhang_mprobs(self):
        """natsel_zhang optimise_motif_probs setting should work"""
        opt = dict(max_evaluations=2, limit_action="ignore")
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        # default, not optimising root probs
        natsel = evo_app.natsel_zhang("MG94HKY", tip1="Human", opt_args=opt, gc=2)
        result = natsel(aln)
        self.assertEqual(result.null.lf.nfp, 6)

        # optimising root probs
        natsel = evo_app.natsel_zhang(
            "MG94HKY", tip1="Human", opt_args=opt, gc=2, optimise_motif_probs=True
        )
        result = natsel(aln)
        self.assertEqual(result.null.lf.nfp, 9)

    def test_neutral(self):
        """test of neutrality, one omega != 1"""
        opt = dict(max_evaluations=20, limit_action="ignore")
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
        neutral = evo_app.natsel_neutral(
            "MG94HKY", tree="data/primate_brca1.tree", opt_args=opt
        )
        result = neutral(aln)
        self.assertEqual(result.df, 1)
        self.assertTrue("MG94HKY-null" in result)
        self.assertTrue("MG94HKY-alt" in result)
        # fails if not a codon model
        with self.assertRaises(ValueError):
            _ = evo_app.natsel_neutral("F81", tree="data/primate_brca1.tree")

    def test_neutral_mtdna(self):
        """test of neutrality, different genetic code"""
        from cogent3.app.composable import NotCompleted

        opt = dict(max_evaluations=2, limit_action="ignore")
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        neutral = evo_app.natsel_neutral("MG94HKY", opt_args=opt, gc=2)
        result = neutral(aln)
        self.assertEqual(result.df, 1)
        # not completed if wrong gc
        neutral = evo_app.natsel_neutral("MG94HKY", opt_args=opt, gc=1)
        result = neutral(aln)
        self.assertIsInstance(result, NotCompleted)

    def test_neutral_mprobs(self):
        """test of neutrality, optimise_motif_probs setting should work"""
        opt = dict(max_evaluations=2, limit_action="ignore")
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        # default, not optimising root probs
        natsel = evo_app.natsel_neutral("MG94HKY", opt_args=opt, gc=2)
        result = natsel(aln)
        self.assertEqual(result.null.lf.nfp, 4)

        # optimising root probs
        natsel = evo_app.natsel_neutral(
            "MG94HKY", opt_args=opt, gc=2, optimise_motif_probs=True
        )
        result = natsel(aln)
        self.assertEqual(result.null.lf.nfp, 7)


class TestTabulateStats(TestCase):
    def test_tabulate(self):
        """call returns tabular_result with Tables"""
        from cogent3.util.table import Table

        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "GN", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        result = mod(aln)
        tabulator = evo_app.tabulate_stats()
        tabulated = tabulator(result)
        self.assertEqual(len(tabulated), 3)
        for title in ("motif params", "global params", "edge params"):
            self.assertTrue(title in tabulated)
            self.assertIsInstance(tabulated[title], Table)


if __name__ == "__main__":
    main()
