from unittest import TestCase, main

from cogent3 import make_aligned_seqs
from cogent3.app import evo as evo_app
from cogent3.app.result import (
    generic_result,
    model_collection_result,
    model_result,
)
from cogent3.util.deserialise import deserialise_object


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2021, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2021.04.20a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestGenericResult(TestCase):
    def test_deserialised_values(self):
        """correctly deserialises values"""
        from cogent3 import DNA

        data = {"type": "cogent3.core.moltype.MolType", "moltype": "dna"}
        result = generic_result(source="blah.json")
        result["key"] = data
        result.deserialised_values()
        got = result["key"]
        self.assertEqual(got, DNA)
        # if we have a type value without "cogent3", leaves as is
        data = {"type": "core.moltype.MolType", "moltype": "dna"}
        result = generic_result(source="blah.json")
        result["key"] = data
        result.deserialised_values()
        got = result["key"]
        self.assertEqual(got, data)
        # or if no "type" entry, leaves as is
        data = {"moltype": "dna"}
        result = generic_result(source="blah.json")
        result["key"] = data
        result.deserialised_values()
        got = result["key"]
        self.assertEqual(got, data)

    def test_repr_str(self):
        """it works"""
        data = {"type": "cogent3.core.moltype.MolType", "moltype": "dna"}
        result = generic_result(source="blah.json")
        result["key"] = data
        r = repr(result)
        s = str(result)

    def test_keys(self):
        """it works"""
        data = {"type": "cogent3.core.moltype.MolType", "moltype": "dna"}
        result = generic_result(source="blah.json")
        result["key"] = data
        keys = result.keys()
        self.assertEqual(keys, ["key"])


class TestModelResult(TestCase):
    def test_model_result_alignment(self):
        """returns alignment from lf"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "F81",
            show_progress=False,
            opt_args=dict(max_evaluations=5, limit_action="ignore"),
        )
        result = mod(aln)
        got = result.alignment
        self.assertEqual(got.to_dict(), _data)

    def test_model_name_lf_name(self):
        """model_result.name is set as lf.name"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "F81",
            name="blah",
            show_progress=False,
            opt_args=dict(max_evaluations=5, limit_action="ignore"),
        )
        result = mod(aln)
        self.assertEqual(result.name, result.lf.name)
        print(result)

    def test_model_result_alignment_split_pos_model(self):
        """returns alignment from lf with split codon positions"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "F81",
            split_codons=True,
            show_progress=False,
            opt_args=dict(max_evaluations=5, limit_action="ignore"),
        )
        result = mod(aln)
        for i in range(1, 4):
            got = result.alignment[i]
            expect = aln[i - 1 :: 3]
            self.assertEqual(got.to_dict(), expect.to_dict())

    def test_model_result_repr_split_pos_model(self):
        """repr works for model_result of split codon positions"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "F81",
            split_codons=True,
            show_progress=False,
            opt_args=dict(max_evaluations=55, limit_action="ignore"),
        )
        result = mod(aln)
        s = repr(result)

    def test_model_result_tree_split_pos_model(self):
        """returns tree from lf with split codon positions"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "F81",
            split_codons=True,
            show_progress=False,
            opt_args=dict(max_evaluations=55, limit_action="ignore"),
        )
        result = mod(aln)
        self.assertTrue(len(result.tree), 3)
        # check the trees are different by summing lengths
        lengths = set()
        for i, t in result.tree.items():
            lengths.add(t.total_length())
        self.assertTrue(len(lengths) > 1)

    def test_model_result_simulate_alignment(self):
        """returns tree from lf with split codon positions"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "F81",
            split_codons=True,
            show_progress=False,
            opt_args=dict(max_evaluations=55, limit_action="ignore"),
        )
        result = mod(aln)
        got = result.simulate_alignment()
        self.assertEqual(len(aln), len(got))
        self.assertNotEqual(aln.to_dict(), got.to_dict())

    def test_model_result_tree_discrete_time(self):
        """returns paralinear lengths"""

        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        model1 = evo_app.model(
            "BH", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        result = model1(aln)
        got = result.tree
        self.assertEqual(
            got.children[0].params["length"], got.children[0].params["paralinear"]
        )

    def test_model_result_setitem(self):
        """TypeError if value a likelihood function, or a dict with correct type"""
        v = dict(type="arbitrary")
        r = model_result(name="one", source="two")
        with self.assertRaises(TypeError):
            r["name"] = v

        with self.assertRaises(TypeError):
            r["name"] = 4

        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        with self.assertRaises(TypeError):
            r["name"] = aln


class TestModelCollectionResult(TestCase):
    _model_results = {}

    def setUp(self):
        """constructs _model_results if they don't already exist"""
        if self._model_results:
            return

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
        mr1 = model1(aln)
        mr2 = model2(aln)
        self._model_results[mr1.name] = mr1
        self._model_results[mr2.name] = mr2

    def test_get_best_model(self):
        """should correctly identify the best model"""
        coll = model_collection_result(None)
        coll.update(self._model_results)
        got = coll.get_best_model()
        # we ensure a model_result instance is returned from the possible set
        self.assertIn(got, self._model_results.values())

    def test_select_model(self):
        """correctly select models"""
        # we ensure a series of model_result instances is returned
        coll = model_collection_result(None)
        coll.update(self._model_results)
        got = coll.select_models()
        self.assertTrue(len(got) > 0)
        possible = list(self._model_results.values())
        for m in got:
            self.assertIn(m, possible)

    def test_model_collection_result_repr(self):
        """constructed result can do the different repr"""
        result = model_collection_result(None)
        coll = model_collection_result(None)
        coll.update(self._model_results)
        got = result.__repr__()
        self.assertIsInstance(got, str)
        got = result._repr_html_()
        self.assertIsInstance(got, str)

    def test_json_roundtrip(self):
        """roundtrip from json correct"""
        coll = model_collection_result(name="blah", source="blah2")
        coll.update(self._model_results)
        self.assertEqual(coll.name, "blah")
        self.assertEqual(coll.source, "blah2")
        orig = coll.__repr__()
        got = deserialise_object(coll.to_json())
        self.assertEqual(got.__repr__(), orig)
        self.assertIsInstance(got, model_collection_result)
        self.assertEqual(got.name, coll.name)
        self.assertEqual(got.source, coll.source)
        # select_models() should not fail
        got = deserialise_object(coll.to_json())
        m = got.select_models()
        self.assertIsInstance(m[0], model_result)


class TestHypothesisResult(TestCase):
    def test_pvalue(self):
        """hypothesis test p-value property"""
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
        self.assertTrue(0 <= result.pvalue <= 1)


if __name__ == "__main__":
    main()
