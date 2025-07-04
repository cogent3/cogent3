import pathlib
from os.path import dirname, join
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import MagicMock

import pytest
from numpy.testing import assert_allclose, assert_raises

from cogent3 import (
    get_app,
    get_dataset,
    load_aligned_seqs,
    make_aligned_seqs,
    make_tree,
    open_data_store,
)
from cogent3.app import evo as evo_app
from cogent3.app import io
from cogent3.app.composable import NotCompleted
from cogent3.app.result import (
    hypothesis_result,
    model_collection_result,
    model_result,
)
from cogent3.evolve.models import get_model
from cogent3.util.deserialise import deserialise_object

data_dir = join(dirname(dirname(__file__)), "data")


class TestModel(TestCase):
    basedir = "data"

    def test_model_opt_mprob_arg(self):
        """argument controls optimisability of motif prob settings"""
        for mn in ("HKY85", "GN", "CNFGTR"):
            for value in (True, False):
                # check setting via sm_args is overridden
                with pytest.raises(ValueError):
                    _ = evo_app.model(
                        mn,
                        optimise_motif_probs=value,
                        sm_args={"optimise_motif_probs": not value},
                    )
                model = evo_app.model(
                    mn,
                    optimise_motif_probs=value,
                )
                assert model._sm._optimise_motif_probs == value
                # check picking a different value for constructor get's overriden
                model = evo_app.model(
                    get_model(mn, optimise_motif_probs=not value),
                    optimise_motif_probs=value,
                )
                assert model._sm._optimise_motif_probs == value

    def test_model_tree(self):
        """allows tree to be string, None or tree"""
        treestring = "(a,b,c)"
        for tree in (treestring, make_tree(treestring=treestring), None):
            mod = evo_app.model("HKY85", tree=tree)
            expect = None if tree is None else make_tree(treestring=treestring)
            assert isinstance(mod._tree, expect.__class__)

    def test_unique_models(self):
        """hypothesis raises ValueError if models not unique"""
        model1 = evo_app.model("HKY85")
        model2 = evo_app.model("HKY85", time_het="max")
        with pytest.raises(ValueError):
            evo_app.hypothesis(model1, model2)

    def test_hyp_init(self):
        """uses user specified init_alt function, or not"""
        opt_args = {"max_evaluations": 25, "limit_action": "ignore"}
        model1 = evo_app.model("F81", opt_args=opt_args)
        model2 = evo_app.model("HKY85", opt_args=opt_args)
        # defaults to using null for init
        hyp = evo_app.hypothesis(model1, model2)
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        result = hyp(aln)
        assert result.df == 1

        # user specified function
        hyp = evo_app.hypothesis(model1, model2, init_alt=lambda x, y: x)
        result = hyp(aln)
        assert result.df == 1

    def test_hyp_init_sequential(self):
        """uses preceding model to initialise function"""
        opt_args = {"max_evaluations": 15, "limit_action": "ignore"}
        model1 = evo_app.model("F81", opt_args=opt_args)
        model2 = evo_app.model("HKY85", opt_args=opt_args)
        model3 = evo_app.model("GTR", opt_args=opt_args)
        # defaults to initialise model3 from model 2 from model1
        hyp = evo_app.hypothesis(model1, model2, model3, sequential=True)
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        result = hyp(aln)
        assert result["F81"].lf.lnL < result["HKY85"].lf.lnL < result["GTR"].lf.lnL

        # can be set to False, in which case all models start at defaults
        hyp = evo_app.hypothesis(model1, model2, model3, sequential=False)
        result = hyp(aln)
        assert not result["F81"].lf.lnL < result["HKY85"].lf.lnL < result["GTR"].lf.lnL

    def test_model_collection_init_sequential(self):
        """model collection uses preceding model to initialise function"""
        opt_args = {"max_evaluations": 15, "limit_action": "ignore"}
        model1 = evo_app.model("F81", opt_args=opt_args)
        model2 = evo_app.model("HKY85", opt_args=opt_args)
        model3 = evo_app.model("GTR", opt_args=opt_args)
        # defaults to initialise model3 from model 2 from model1
        mod_coll = evo_app.model_collection(model1, model2, model3, sequential=True)
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna", info={"source": "blah"})
        result = mod_coll(aln)
        assert result["F81"].lf.lnL < result["HKY85"].lf.lnL < result["GTR"].lf.lnL

        # can be set to False, in which case all models start at defaults
        mod_coll = evo_app.hypothesis(model1, model2, model3, sequential=False)
        result = mod_coll(aln)
        assert not result["F81"].lf.lnL < result["HKY85"].lf.lnL < result["GTR"].lf.lnL

        assert isinstance(result, model_collection_result)

        # now with a single discrete edge
        lf_args = {"discrete_edges": ["Opossum"]}
        model2 = evo_app.model("HKY85", opt_args=opt_args, lf_args=lf_args)
        model3 = evo_app.model("GTR", opt_args=opt_args, lf_args=lf_args)
        # defaults to initialise model3 from model 2 from model1
        mod_coll = evo_app.model_collection(model2, model3, sequential=True)
        result = mod_coll(aln)
        assert isinstance(result, model_collection_result)

    def test_model_time_het(self):
        """support lf time-het argument edge_sets"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        mod = evo_app.model(
            "GN",
            optimise_motif_probs=True,
            time_het=[{"edges": ["Mouse", "Human"], "is_independent": False}],
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        result = mod(aln)
        # 11 free params per calibrated GN matrix, there are 2
        # 3 params for root motif probs, 3 branch lengths
        expect_nfp = 11 * 2 + 3 + 3
        assert result.lf.nfp == expect_nfp

    def test_setting_model_bounds(self):
        upper = 10.0
        lower = 0.5
        app = evo_app.model(
            "HKY85",
            optimise_motif_probs=True,
            show_progress=False,
            unique_trees=True,
            time_het="max",
            lower=lower,
            upper=upper,
        )

        aln = make_aligned_seqs(
            {"s1": "ACGT", "s2": "ACGC", "s3": "AAGT"},
            moltype="dna",
        )
        result = app(aln)
        rules = result.lf.get_param_rules()
        kappa_bounds = {
            (r["lower"], r["upper"]) for r in rules if r["par_name"] == "kappa"
        }
        assert kappa_bounds == {(lower, upper)}

    def test_model_param_rules(self):
        """applies upper bound if sensible"""
        mod = evo_app.model(
            "GN",
            param_rules=[
                {"par_name": "length", "edge": "Mouse", "is_independent": False},
            ],
        )
        assert mod._param_rules[0].get("upper") == 50
        mod = evo_app.model(
            "GN",
            param_rules=[{"par_name": "length", "edge": "Mouse", "is_constant": True}],
        )
        assert mod._param_rules[0].get("upper", None) is None

    def test_discrete_time_model(self):
        """works with discrete-time submodel"""
        from cogent3.app.composable import NotCompleted

        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        mod = evo_app.model(
            "BH",
            opt_args={"max_evaluations": 100, "limit_action": "ignore"},
        )
        r = mod(aln)
        assert not isinstance(r, NotCompleted)

    def test_model_hypothesis_result_repr(self):
        """result objects __repr__ and _repr_html_ methods work correctly"""
        import re

        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        model1 = evo_app.model(
            "F81",
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        model2 = evo_app.model(
            "HKY85",
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        hyp = evo_app.hypothesis(model1, model2)
        result = hyp(aln)
        # check the p-val formatted as %.4f
        pval = str(result).splitlines()[4].split()[-1]
        assert re.search("\\d\\.\\d+", pval) is not None
        assert isinstance(result.__repr__(), str)
        assert isinstance(result._repr_html_(), str)
        assert isinstance(result.null.__repr__(), str)
        assert isinstance(result.null._repr_html_(), str)
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
        aln = aln.take_seqs(["Human", "Rhesus", "Galago"])[2::3].omit_gap_pos()
        model1 = evo_app.model(
            "F81",
            optimise_motif_probs=False,
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        model2 = evo_app.model(
            "HKY85",
            optimise_motif_probs=False,
            opt_args={"max_evaluations": 100, "limit_action": "ignore"},
        )
        hyp = evo_app.hypothesis(model1, model2)
        result = hyp(aln)
        pval = str(result).splitlines()[4].split()[-1]
        assert re.search("[0-9\\.]+e-\\d+", pval) is not None

    def test_split_pos_model(self):
        """model with split codons, access .lf using codon position int"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        tree = make_tree(tip_names=aln.names)
        mod = evo_app.model(
            "F81",
            tree=tree,
            split_codons=True,
            opt_args={"max_evaluations": 5, "limit_action": "ignore"},
        )
        result = mod(aln)
        aln1 = result.lf[1].get_param_value("alignment").to_dict()
        assert aln1 == aln[::3].to_dict()

    def test_split_codon_model_result_json(self):
        """round trip split_codon result"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        tree = make_tree(tip_names=aln.names)
        mod = evo_app.model(
            "F81",
            tree=tree,
            split_codons=True,
            opt_args={"max_evaluations": 5, "limit_action": "ignore"},
        )
        result = mod(aln)
        lf1 = result.lf[1]
        json = result.to_json()
        deser = deserialise_object(json)
        assert_allclose(deser.lf[1].lnL, lf1.lnL)

    def test_model_summed_branch_lengths(self):
        """returns summed branch lengths"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        model1 = evo_app.model(
            "F81",
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        result = model1(aln)
        tree = result.lf.get_annotated_tree()
        assert_allclose(result.total_length(), tree.total_length())
        tree = result.lf.get_annotated_tree(length_as="paralinear")
        assert_allclose(
            result.total_length(length_as="paralinear"),
            tree.total_length(),
        )

    def test_model_result_total_length(self):
        """returns summed branch lengths"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        model1 = evo_app.model(
            "GN",
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
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
        aln = make_aligned_seqs(_data, moltype="dna")
        model1 = evo_app.model(
            "GN",
            split_codons=True,
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        result = model1(aln)
        expect = 0.0
        for lf in result.lf.values():
            tree = lf.get_annotated_tree(length_as="ENS")
            expect += tree.total_length()

        got = result.total_length(length_as="ENS")
        assert_allclose(got, expect)

    def test_model_tree_unique_trees(self):
        """handles case of using unique trees for each alignment"""
        with pytest.raises(AssertionError):
            _ = evo_app.model("GN", tree="(a,b,c)", unique_trees=True)
        _data1 = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        _data2 = {
            "Dog": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }

        aln1 = make_aligned_seqs(_data1, moltype="dna")
        aln2 = make_aligned_seqs(_data2, moltype="dna")
        model = evo_app.model(
            "GN",
            unique_trees=True,
            opt_args={"max_evaluations": 2, "limit_action": "ignore"},
        )
        for aln in (aln1, aln2):
            result = model(aln)
            assert isinstance(result, model_result)

        # but the second one fails if unique_trees=False
        model = evo_app.model(
            "GN",
            unique_trees=False,
            opt_args={"max_evaluations": 2, "limit_action": "ignore"},
        )
        for aln, expect_type in ((aln1, model_result), (aln2, NotCompleted)):
            result = model(aln)
            assert isinstance(result, expect_type)


def _make_getter(val):
    def call(**kwargs):
        return val

    return call


def _make_hyp(aic1, aic2, aic3, nfp1, nfp2, nfp3):
    null = _make_mock_result("unrooted", aic1, nfp1)
    alt1 = _make_mock_result("alt1", aic2, nfp2)
    alt2 = _make_mock_result("alt2", aic3, nfp3)
    # this is a really ugly hack to address type validation on result setitem!
    hypothesis_result._item_types = ("model_result", "MagicMock")
    hyp = hypothesis_result("unrooted", source="something")
    for m in (null, alt1, alt2):
        hyp[m.name] = m
    return hyp


def _make_mock_result(arg0, arg1, arg2):
    result = MagicMock()
    result.name = arg0
    result.lf.get_aic = _make_getter(arg1)
    result.nfp = arg2
    return result


class TestHypothesisResult(TestCase):
    def test_get_best_model(self):
        """should correctly identify the best model"""
        # aic order is null, alt1, alt2
        # in this case, no substantial diff, so should return smaller nfp, ie null
        hyp = _make_hyp(112, 110, 111, 10, 11, 12)
        got = hyp.get_best_model(threshold=0.05)
        assert got is hyp.null
        # here alt2 is winner
        hyp = _make_hyp(110, 111, 104, 10, 11, 12)
        got = hyp.get_best_model(threshold=0.05)
        assert got is hyp["alt2"]
        # but if we set threshold more permissive, it will return null
        got = hyp.get_best_model(threshold=0.03)
        assert got is hyp.null

    def test_select_model(self):
        """correctly identify models"""
        hyp = _make_hyp(112, 110, 111, 10, 11, 12)
        got = set(hyp.select_models(threshold=0.05))
        expect = set(hyp.values())
        assert got == expect
        # single model
        hyp = _make_hyp(110, 111, 104, 10, 11, 12)
        got = hyp.select_models(threshold=0.05)
        assert len(got) == 1
        assert got[0] is hyp["alt2"]
        # but if we set threshold more permissive, it will return all
        got = hyp.select_models(threshold=0.03)
        assert len(got) == 3
        expect = set(hyp.values())
        assert set(got) == expect

    def test_null_hyp_fail_error(self):
        """if null fails NotCompleted.origin should be model"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        tree = "((Mouse,Rat),Human,Opossum)"
        m1 = evo_app.model("F81", tree=tree)
        m2 = evo_app.model("GTR", tree=tree)
        hyp = evo_app.hypothesis(m1, m2)
        r = hyp(aln)
        assert r.origin == "model"

    def test_hyp_split_codon_select_models(self):
        """hypothesis_result identifies selects best model when split_codon"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        opt_args = {"max_evaluations": 10, "limit_action": "ignore"}
        m1 = evo_app.model(
            "F81",
            optimise_motif_probs=False,
            split_codons=True,
            opt_args=opt_args,
            param_rules=[
                {"par_name": "length", "upper": 10.0, "lower": 1e-09},
            ],
            lower=1e-06,
            upper=1000000.0,
        )
        m2 = evo_app.model(
            "GTR",
            optimise_motif_probs=False,
            split_codons=True,
            opt_args=opt_args,
            param_rules=[
                {"par_name": "length", "upper": 10.0, "lower": 1e-09},
            ],
            lower=1e-06,
            upper=1000000.0,
        )
        hyp = evo_app.hypothesis(m1, m2)
        r = hyp(aln)
        bm = r.select_models()
        assert_allclose(bm[0].lnL, -85.00043312185628)

    def test_alt_hyp_fail_error(self):
        """if alt fails NotCompleted.origin should be model"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGA",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGA",
            "Opossum": "TGACCAGTGAAAGTGGCGGCGGTGGCTGA",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        tree = "(Mouse,Human,Opossum)"
        m1 = evo_app.model("F81", tree=tree)
        m2 = evo_app.model("MG94HKY", tree=tree)
        hyp = evo_app.hypothesis(m1, m2)
        r = hyp(aln)
        assert r.origin == "model"

    def test_model_moltype_mismatch(self):
        """if model and alignment moltypes incompatible"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGA",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGA",
            "Opossum": "TGACCAGTGAAAGTGGCGGCGGTGGCTGA",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        tree = "(Mouse,Human,Opossum)"
        m1 = evo_app.model("JTT92", tree=tree)
        r = m1(aln)
        assert r.origin == "model"


class TestAncestralStates(TestCase):
    def test_ancestral(self):
        """recon ancestral states works"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        mod = evo_app.model(
            "GN",
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        anc = evo_app.ancestral_states()
        result = anc(mod(aln))
        assert result["root"].shape == (len(aln), 4)
        assert_allclose(result["root"].row_sum(), 1)


class TestNatSel(TestCase):
    # needs to work for single edge, just two edges, and all combos of clae,
    # stem etc..
    def test_zhang(self):
        """natsel_zhang correctly configured and should not fail"""
        opt = {"max_evaluations": 20, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
        natsel = evo_app.natsel_zhang(
            "CNFGTR",
            tree="data/primate_brca1.tree",
            tip1="Human",
            tip2="Chimpanzee",
            opt_args=opt,
        )
        result = natsel(aln)
        assert result.df == 3
        assert result.alt.nfp == 21
        # the naming scheme is model name followed by null/alt
        assert "CNFGTR-null" in result
        assert "CNFGTR-alt" in result

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
        assert result.df == 3
        assert "Y98-null" in result
        assert "Y98-alt" in result

        # fails if not a codon model
        with pytest.raises(ValueError):
            _ = evo_app.natsel_zhang(
                "F81",
                tree="data/primate_brca1.tree",
                tip1="Human",
                tip2="Chimpanzee",
                opt_args=opt,
            )

        # fails if no tip names provided
        with pytest.raises(ValueError):
            _ = evo_app.natsel_zhang(
                "Y98",
                tree="data/primate_brca1.tree",
                opt_args=opt,
            )

    def test_zhang_mtseq(self):
        """genetic code setting should work"""
        from cogent3.app.composable import NotCompleted

        opt = {"max_evaluations": 20, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        natsel = evo_app.natsel_zhang("CNFGTR", tip1="Human", opt_args=opt, gc=2)
        result = natsel(aln)
        assert result.df == 3
        # but if provide wrong gc, get NotCompleted
        natsel = evo_app.natsel_zhang("CNFGTR", tip1="Human", opt_args=opt, gc=1)
        result = natsel(aln)
        assert isinstance(result, NotCompleted)

    def test_zhang_mprobs(self):
        """natsel_zhang optimise_motif_probs setting should work"""
        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        # default, not optimising root probs
        natsel = evo_app.natsel_zhang("MG94HKY", tip1="Human", opt_args=opt, gc=2)
        result = natsel(aln)
        assert result.null.lf.nfp == 6

        # optimising root probs
        natsel = evo_app.natsel_zhang(
            "MG94HKY",
            tip1="Human",
            opt_args=opt,
            gc=2,
            optimise_motif_probs=True,
        )
        result = natsel(aln)
        assert result.null.lf.nfp == 9

    def test_neutral(self):
        """test of neutrality, one omega != 1"""
        opt = {"max_evaluations": 20, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
        neutral = evo_app.natsel_neutral(
            "MG94HKY",
            tree="data/primate_brca1.tree",
            opt_args=opt,
        )
        result = neutral(aln)
        assert result.df == 1
        assert "MG94HKY-null" in result
        assert "MG94HKY-alt" in result
        # fails if not a codon model
        with pytest.raises(ValueError):
            _ = evo_app.natsel_neutral("F81", tree="data/primate_brca1.tree")

    def test_neutral_mtdna(self):
        """test of neutrality, different genetic code"""
        from cogent3.app.composable import NotCompleted

        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        neutral = evo_app.natsel_neutral("MG94HKY", opt_args=opt, gc=2)
        result = neutral(aln)
        assert result.df == 1
        # not completed if wrong gc
        neutral = evo_app.natsel_neutral("MG94HKY", opt_args=opt, gc=1)
        result = neutral(aln)
        assert isinstance(result, NotCompleted)

    def test_neutral_mprobs(self):
        """test of neutrality, optimise_motif_probs setting should work"""
        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        # default, not optimising root probs
        natsel = evo_app.natsel_neutral("MG94HKY", opt_args=opt, gc=2)
        result = natsel(aln)
        assert result.null.lf.nfp == 4

        # optimising root probs
        natsel = evo_app.natsel_neutral(
            "MG94HKY",
            opt_args=opt,
            gc=2,
            optimise_motif_probs=True,
        )
        result = natsel(aln)
        assert result.null.lf.nfp == 7

    def test_neutral_nstat_model(self):
        """test of neutrality, non-stationary codon model"""
        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        neutral = evo_app.natsel_neutral("GNC", opt_args=opt, gc=2)
        result = neutral(aln)
        # 11 rate matrix params for GNC (omega omitted in null), 3 edges
        assert result.null.lf.nfp == 3 + 11

    def test_natsel_sitehet(self):
        """site-het natsel hypothesis test"""
        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
        # default, not optimising root probs
        natsel = evo_app.natsel_sitehet(
            "MG94HKY",
            tree="data/primate_brca1.tree",
            opt_args=opt,
        )
        result = natsel(aln)
        # one free param for each edge, 1 for kappa, 1 for omega, 1 for bprobs
        assert result.null.lf.nfp == 14
        # plus one extra bprob and one extra omega
        assert result.alt.lf.nfp == 16
        # fails if not a codon model
        with pytest.raises(ValueError):
            _ = evo_app.natsel_sitehet("F81", tree="data/primate_brca1.tree")

    def test_natsel_sitehet_mprob(self):
        """natsel_sitehet correctly applies genetic code and optimise_motif_probs args"""
        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        # optimising root probs
        natsel = evo_app.natsel_sitehet(
            "MG94HKY",
            opt_args=opt,
            gc=2,
            optimise_motif_probs=True,
        )
        # test of genetic code is implicit, if not correct, the following
        # call would return NotCompleted (for this mtDNA gene), which does not
        # have a .null attribute
        result = natsel(aln)
        # 3 edges, 1 kappa, 1 omega, 1 bprob, 3 mprob
        assert result.null.lf.nfp == 9

    def test_natsel_timehet(self):
        """natsel_timehet works"""
        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
        natsel = evo_app.natsel_timehet(
            "MG94HKY",
            tree="data/primate_brca1.tree",
            tip1="Human",
            tip2="Chimpanzee",
            opt_args=opt,
        )
        result = natsel(aln)
        assert result.df == 1
        # the naming scheme is model name followed by null/alt
        assert "MG94HKY-null" in result
        assert "MG94HKY-alt" in result
        # that is_independent works
        natsel = evo_app.natsel_timehet(
            "MG94HKY",
            tree="data/primate_brca1.tree",
            tip1="Human",
            tip2="Chimpanzee",
            is_independent=True,
            opt_args=opt,
        )
        result = natsel(aln)
        assert result.df == 2

        # handle specifying just single edge
        natsel = evo_app.natsel_timehet(
            "MG94HKY",
            tree="data/primate_brca1.tree",
            tip1="Human",
            opt_args=opt,
        )
        result = natsel(aln)
        assert result.df == 1

        # fails if not a codon model
        with pytest.raises(ValueError):
            _ = evo_app.natsel_timehet("F81", tip1="Human")

    def test_natsel_timehet_mprobs(self):
        """natsel_timehet works with gc and mprobs settings"""
        opt = {"max_evaluations": 2, "limit_action": "ignore"}
        aln = load_aligned_seqs("data/ENSG00000198712.fa", moltype="dna")
        natsel = evo_app.natsel_timehet(
            "MG94HKY",
            tip1="Human",
            tip2="Chimp",
            opt_args=opt,
            gc=2,
            optimise_motif_probs=True,
        )
        result = natsel(aln)
        assert result.df == 1
        assert result.null.lf.nfp == 3 + 3 + 1 + 1
        assert result.alt.lf.nfp == 3 + 3 + 1 + 2
        # the naming scheme is model name followed by null/alt
        assert "MG94HKY-null" in result
        assert "MG94HKY-alt" in result


class TestTabulateStats(TestCase):
    def test_tabulate(self):
        """call returns tabular_result with Tables"""
        from cogent3.core.table import Table

        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(_data, moltype="dna")
        mod = evo_app.model(
            "GN",
            opt_args={"max_evaluations": 25, "limit_action": "ignore"},
        )
        result = mod(aln)
        tabulator = evo_app.tabulate_stats()
        tabulated = tabulator(result)
        assert len(tabulated) == 3
        for title in ("motif params", "global params", "edge params"):
            assert title in tabulated
            assert isinstance(tabulated[title], Table)


class TestBootstrap(TestCase):
    """testing the bootstrap app"""

    def test_bstrap(self):
        """exercising bootstrap with simple hypothesis"""
        aln = load_aligned_seqs(join(data_dir, "brca1.fasta"), moltype="dna")
        aln = aln.take_seqs(aln.names[:3])
        aln = aln.omit_gap_pos(allowed_gap_frac=0)
        opt_args = {"max_evaluations": 20, "limit_action": "ignore"}
        m1 = evo_app.model("F81", opt_args=opt_args)
        m2 = evo_app.model("HKY85", opt_args=opt_args)
        hyp = evo_app.hypothesis(m1, m2)
        strapper = evo_app.bootstrap(hyp, num_reps=2, parallel=False)
        result = strapper(aln)
        nd = result.null_dist
        assert {type(v) for v in nd}, {float}
        json = result.to_json()
        got = deserialise_object(json)
        assert isinstance(got, evo_app.bootstrap_result)

    def test_bstrap_fail(self):
        """invalid data returns meaningful error"""
        aln = load_aligned_seqs(join(data_dir, "brca1.fasta"), moltype="dna")
        aln = aln.take_seqs(aln.names[:3])
        opt_args = {"max_evaluations": 20, "limit_action": "ignore"}
        m1 = evo_app.model("F81", opt_args=opt_args)
        # we've retained gaps, so this should fail at first call as incompatible with model
        m2 = evo_app.model("GTR", opt_args=opt_args, sm_args={"recode_gaps": False})
        hyp = evo_app.hypothesis(m1, m2)
        strapper = evo_app.bootstrap(hyp, num_reps=2, parallel=False)
        result = strapper(aln)
        # correct message being relayed
        assert "ValueError: '-' at" in result.message

    def test_bootstrap_composability(self):
        """can be composed with load_db and write_db"""
        m1 = evo_app.model("F81")
        m2 = evo_app.model("HKY85")
        hyp = evo_app.hypothesis(m1, m2)
        with TemporaryDirectory(dir=".") as dirname:
            dirname = pathlib.Path(dirname)
            out_dstore = open_data_store(dirname / "delme.sqlitedb", mode="w")
            writer = io.write_db(out_dstore)
            _ = io.load_db() + evo_app.bootstrap(hyp, num_reps=2) + writer


def test_bstrap_parallel():
    """exercising bootstrap with parallel"""
    aln = load_aligned_seqs(join(data_dir, "brca1.fasta"), moltype="dna")
    aln = aln.take_seqs(aln.names[:3])
    aln = aln.omit_gap_pos(allowed_gap_frac=0)
    opt_args = {"max_evaluations": 20, "limit_action": "ignore"}
    m1 = evo_app.model("F81", opt_args=opt_args)
    m2 = evo_app.model("HKY85", opt_args=opt_args)
    hyp = evo_app.hypothesis(m1, m2)
    strapper = evo_app.bootstrap(hyp, num_reps=2, parallel=True)
    result = strapper(aln)
    assert isinstance(result, evo_app.bootstrap_result)


def test_model_opt_args():
    opt_args = {"max_restarts": 10, "tolerance": 1e-8}

    model = evo_app.model(
        "GN",
        opt_args=opt_args,
    )
    assert model._opt_args == {**opt_args, "show_progress": False}


@pytest.mark.internet
def test_get_app_tree_is_url():
    """A check that the model app can use a url for the tree"""
    tree_url = "https://github.com/user-attachments/files/20329020/brca1_5.tree.gz"
    mod = get_app("model", "F81", tree=tree_url)
    assert isinstance(mod, evo_app.model)


def test_model_str():
    """correct str representation"""
    model = evo_app.model("HKY85", time_het="max")
    got = " ".join(str(model).splitlines())
    expect_start = "model(sm='HKY85', tree=None, unique_trees=False, "
    assert got.startswith(expect_start)


def test_hypothesis_str():
    """correct str representation"""
    model1 = evo_app.model("HKY85")
    model2 = evo_app.model("HKY85", name="hky85-max-het", time_het="max")
    hyp = evo_app.hypothesis(model1, model2)
    got = " ".join(str(hyp).splitlines())
    expect_start = "hypothesis(null=model(sm='HKY85', tree=None, unique_trees=False, "
    assert got.startswith(expect_start)


def test_model_tree_func(DATA_DIR):
    dist_cal = get_app("fast_slow_dist", fast_calc="paralinear", moltype="dna")
    est_tree = get_app("quick_tree", drop_invalid=True)
    app = dist_cal + est_tree
    model = get_app(
        "model",
        "HKY85",
        tree_func=app,
        opt_args={"max_evaluations": 25, "limit_action": "ignore"},
    )

    aln = load_aligned_seqs(DATA_DIR / "brca1_5.paml", moltype="dna")
    result = model(aln)
    assert set(result.tree.get_tip_names()) == set(aln.names)


def test_model_invalid_tree_func():
    with pytest.raises(AssertionError):
        get_app(
            "model",
            "HKY85",
            tree_func="123",
        )


def test_model_bounds_allpar():
    upper = 10.0
    lower = 0.5
    app = get_app(
        "model",
        "HKY85",
        optimise_motif_probs=True,
        show_progress=False,
        unique_trees=True,
        lower=lower,
        upper=upper,
    )

    aln = make_aligned_seqs(
        {"s1": "ACGT", "s2": "ACGC", "s3": "AAGT"},
        moltype="dna",
    )
    result = app(aln)
    rules = result.lf.get_param_rules()
    par_bounds = {(r["lower"], r["upper"]) for r in rules if r["par_name"] == ("kappa")}
    assert par_bounds == {(lower, upper)}


def test_model_bounds_kappa():
    upper_kappa = 99
    lower_kappa = 9
    app = get_app(
        "model",
        "HKY85",
        optimise_motif_probs=True,
        show_progress=False,
        unique_trees=True,
        param_rules=[{"par_name": "kappa", "upper": upper_kappa, "lower": lower_kappa}],
    )

    aln = make_aligned_seqs(
        {"s1": "ACGT", "s2": "ACGC", "s3": "AAGT"},
        moltype="dna",
    )
    result = app(aln)
    rules = result.lf.get_param_rules()
    kappa_bounds = {(r["lower"], r["upper"]) for r in rules if r["par_name"] == "kappa"}
    assert kappa_bounds == {(lower_kappa, upper_kappa)}


def test_source_propagated():
    aln = get_dataset("brca1")
    three = get_app("take_n_seqs", 4)
    dcalc = get_app("fast_slow_dist", fast_calc="pdist", moltype="dna")
    qtree = get_app("quick_tree")
    app = three + dcalc + qtree
    result = app(aln)
    assert result.source == pathlib.Path(aln.source).name
