import pathlib

from unittest import TestCase, main

from cogent3 import make_aligned_seqs, make_table
from cogent3.app import evo as evo_app
from cogent3.app.data_store import DataStoreMember
from cogent3.app.result import (
    generic_result,
    hypothesis_result,
    model_collection_result,
    model_result,
    tabular_result,
)
from cogent3.util.deserialise import deserialise_object
from cogent3.util.dict_array import DictArray


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
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
        repr(result)
        str(result)

    def test_keys(self):
        """it works"""
        data = {"type": "cogent3.core.moltype.MolType", "moltype": "dna"}
        result = generic_result(source="blah.json")
        result["key"] = data
        keys = result.keys()
        self.assertEqual(keys, ["key"])

    def test_invalid_setitem(self):
        """generic_result raise TypeError if trying to set invalid item type for json"""
        gr = generic_result("null")
        with self.assertRaises(TypeError):
            gr["null"] = {0, 23}

    def test_infers_source(self):
        """flexible handling of data source"""
        # works for string
        source = "path/blah.fasta"
        aln = make_aligned_seqs(
            {"A": "ACGT"}, info=dict(source=source, random_key=1234)
        )
        gr = generic_result(aln)
        self.assertEqual(gr.source, "path/blah.fasta")

        # or Path
        aln.info.source = pathlib.Path(source)
        gr = generic_result(aln)
        self.assertEqual(str(gr.source), str(pathlib.Path("path/blah.fasta")))

        # or DataStoreMember
        aln.info.source = DataStoreMember(source)
        gr = generic_result(aln)
        self.assertEqual(str(gr.source), "path/blah.fasta")

        aln.info = {}
        with self.assertRaises(ValueError):
            generic_result(aln)


class TestModelResult(TestCase):
    def test_repr(self):
        """does not fail"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = make_aligned_seqs(data=_data, moltype="dna")
        mod = evo_app.model(
            "F81",
            show_progress=False,
            opt_args=dict(max_evaluations=1, limit_action="ignore"),
        )
        result = mod(aln)
        self.assertIsInstance(repr(result), str)
        # no values set
        self.assertIsInstance(repr(model_result(source="blah")), str)

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
        repr(result)

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
        lengths = {t.total_length() for _, t in result.tree.items()}
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

    def test_repr_str(self):
        """it works even when no values"""
        mr = model_result(source="blah")
        self.assertIsInstance(repr(mr), str)

    def test_model_result_invalid_setitem(self):
        """model_result raise TypeError if trying to set incorrect item type"""
        mr = model_result(source="blah")
        with self.assertRaises(TypeError):
            mr["null"] = 23


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
        model3 = evo_app.model(
            "GTR", opt_args=dict(max_evaluations=25, limit_action="ignore")
        )
        mr1 = model1(aln)
        mr2 = model2(aln)
        mr3 = model3(aln)
        self._model_results[mr1.name] = mr1
        self._model_results[mr2.name] = mr2
        self._model_results[mr3.name] = mr3

    def test_get_best_model(self):
        """should correctly identify the best model"""
        coll = model_collection_result(source="blah")
        coll.update(self._model_results)
        got = coll.get_best_model()
        # we ensure a model_result instance is returned from the possible set
        self.assertIn(got, self._model_results.values())

    def test_select_model(self):
        """correctly select models"""
        # we ensure a series of model_result instances is returned
        coll = model_collection_result(source="blah")
        coll.update(self._model_results)
        got = coll.select_models()
        self.assertTrue(len(got) > 0)
        possible = list(self._model_results.values())
        for m in got:
            self.assertIn(m, possible)

    def test_model_collection_result_repr(self):
        """constructed result can do the different repr"""
        result = model_collection_result(source="blah")
        coll = model_collection_result(source="blah")
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

    def test_to_hypothesis(self):
        """creates a hypothesis_result from two model results"""
        mr = model_collection_result(source="blah")
        mr.update(self._model_results)
        hyp = mr.get_hypothesis_result("F81", "HKY85")
        self.assertIsInstance(hyp, hypothesis_result)
        self.assertEqual(hyp.null.name, "F81")

    def test_repr_str(self):
        """it works even when no values"""
        mr = model_collection_result(source="blah")
        self.assertIsInstance(repr(mr), str)

    def test_model_collection_result_invalid_setitem(self):
        """model_collection_result raise TypeError if trying to set incorrect item type"""
        mcr = model_collection_result(source="blah")
        with self.assertRaises(TypeError):
            mcr["null"] = 23


class TestHypothesisResult(TestCase):
    def test_repr_str(self):
        """it works even when no values"""
        hr = hypothesis_result(name_of_null="null", source="blah")
        self.assertIsInstance(repr(hr), str)

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

    def test_invalid_setitem(self):
        """hypothesis_result raise TypeError if trying to set incorrect item type"""
        hr = hypothesis_result(name_of_null="null", source="blah")
        with self.assertRaises(TypeError):
            hr["null"] = {0, 23}


class TestTabularResult(TestCase):
    def test_valid_setitem(self):
        """tabular_result works when set correct item type"""
        tr = tabular_result("null")
        tr["result"] = make_table(data={"A": [0, 1]})
        darr = DictArray({"A": [0, 1]})
        tr["result2"] = darr
        js = tr.to_json()
        self.assertIsInstance(js, str)

    def test_invalid_setitem(self):
        """tabular_result raise TypeError if trying to set incorrect item type"""
        tr = tabular_result("null")
        with self.assertRaises(TypeError):
            tr["null"] = {0, 23}


if __name__ == "__main__":
    main()
