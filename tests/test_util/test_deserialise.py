import json

from tempfile import TemporaryDirectory

from numpy.testing import assert_allclose

from cogent3 import LoadSeqs, LoadTree
from cogent3.app.result import model_result
from cogent3.core import alignment, moltype
from cogent3.evolve.models import get_model
from cogent3.util.deserialise import deserialise_object
from cogent3.util.unit_test import TestCase, main


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.20a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestDeserialising(TestCase):
    def test_roundtrip_codon_alphabet(self):
        """codon alphabet to_json enables roundtrip"""
        data = moltype.STANDARD_CODON.to_json()
        got = deserialise_object(data)
        self.assertEqual(type(got), type(moltype.STANDARD_CODON))
        self.assertEqual(list(got), list(moltype.STANDARD_CODON))

    def test_roundtrip_alphabet(self):
        """alphabet to_json enables roundtrip"""
        dna = moltype.get_moltype("dna")
        data = dna.alphabet.to_json()
        got = deserialise_object(data)
        self.assertEqual(type(got), type(dna.alphabet))
        self.assertEqual(list(got), list(dna.alphabet))

    def test_roundtrip_moltype(self):
        """moltype to_json enables roundtrip"""
        dna = moltype.get_moltype("dna")
        data = dna.to_json()
        got = deserialise_object(data)
        self.assertEqual(type(got), type(dna))
        self.assertEqual(list(got), list(dna))
        self.assertEqual(dna, got)

    def test_roundtrip_seq(self):
        """seq to_json enables roundtrip"""
        for mtype in ("dna", "protein"):
            mtype = moltype.get_moltype(mtype)
            seq = mtype.make_seq("ACGGTCGG", "label", info={"something": 3})
            got = deserialise_object(seq.to_json())
            self.assertEqual(got.info.something, 3)
            self.assertEqual(got.name, "label")
            self.assertEqual(got.moltype, seq.moltype)
            self.assertEqual(str(got), str(seq))

    def test_roundtrip_seqcoll(self):
        """SequenceCollection to_json enables roundtrip"""
        data = dict(A="TTGT", B="GGCT")
        seqcoll = LoadSeqs(data=data, moltype="dna", aligned=False)
        got = deserialise_object(seqcoll.to_json())
        self.assertEqual(got.rc().to_dict(), seqcoll.rc().to_dict())
        self.assertIsInstance(got, alignment.SequenceCollection)

    def test_roundtrip_arrayalign(self):
        """ArrayAlignment to_json enables roundtrip"""
        data = dict(A="TTGTA", B="GGCT-")
        arrayalign = LoadSeqs(data=data, moltype="dna")
        got = deserialise_object(arrayalign.to_json())
        self.assertEqual(got.rc().to_dict(), arrayalign.rc().to_dict())
        self.assertIsInstance(got, alignment.ArrayAlignment)

    def test_roundtrip_align(self):
        """Alignment to_json enables roundtrip"""
        data = dict(A="TTGTA", B="GGCT-")
        align = LoadSeqs(data=data, moltype="dna", array_align=False)
        got = deserialise_object(align.to_json())
        self.assertEqual(got.rc().to_dict(), align.rc().to_dict())
        self.assertIsInstance(got, alignment.Alignment)

    def test_roundtrip_tree(self):
        """Tree to_json enables roundtrip"""
        tree = LoadTree(treestring="(c:01,d:0.3,(a:0.05,b:0.08)xx:0.2)")
        got = deserialise_object(tree.to_json())
        self.assertFloatEqual(got.get_node_matching_name("a").length, 0.05)
        self.assertFloatEqual(got.get_node_matching_name("xx").length, 0.2)

    def test_roundtrip_submod(self):
        """substitution model to_json enables roundtrip"""
        sm = get_model("HKY85")
        data = sm.to_json()
        got = deserialise_object(data)
        self.assertEqual(got.to_rich_dict(), sm.to_rich_dict())
        sm = get_model("GN")
        data = sm.to_json()
        got = deserialise_object(data)
        self.assertEqual(got.to_rich_dict(), sm.to_rich_dict())
        sm = get_model("CNFGTR")
        data = sm.to_json()
        got = deserialise_object(data)
        self.assertEqual(got.to_rich_dict(), sm.to_rich_dict())

    def test_roundtrip_likelihood_function(self):
        """likelihood function.to_json enables roundtrip"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = LoadSeqs(data=_data, moltype="dna")
        tree = LoadTree(tip_names=aln.names)
        sm = get_model("HKY85")
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        edge_vals = zip(aln.names, (2, 3, 4))
        for edge, val in edge_vals:
            lf.set_param_rule("kappa", edge=edge, init=val)
        lnL = lf.get_log_likelihood()
        data = lf.to_json()
        got_obj = deserialise_object(data)
        self.assertFloatEqual(got_obj.get_log_likelihood(), lnL)

    def test_roundtrip_het_lf(self):
        """correctly round trips a site-het model"""
        with open("data/site-het-param-rules.json") as infile:
            rules = json.load(infile)

        aln = LoadSeqs("data/primates_brca1.fasta", moltype="dna")
        tree = LoadTree("data/primates_brca1.tree")
        rule_lnL = rules.pop("phylohmm-gamma-kappa")
        sm = get_model("HKY85", ordered_param="rate", distribution="gamma")
        lf1 = sm.make_likelihood_function(tree, bins=4, sites_independent=False)
        lf1.set_alignment(aln)
        lf1.apply_param_rules(rule_lnL["rules"])
        data = lf1.to_json()
        got_lf = deserialise_object(data)
        assert_allclose(lf1.lnL, got_lf.lnL)

    def test_roundtrip_from_file(self):
        """correctly roundtrips a likelihood function fro json file"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = LoadSeqs(data=_data, moltype="dna")
        tree = LoadTree(tip_names=aln.names)
        sm = get_model("HKY85")
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        edge_vals = zip(aln.names, (2, 3, 4))
        for edge, val in edge_vals:
            lf.set_param_rule("kappa", edge=edge, init=val)
        lnL = lf.get_log_likelihood()
        data = lf.to_json()
        with TemporaryDirectory(dir=".") as dirname:
            outpath = dirname + "/delme.json"
            with open(outpath, "w") as outfile:
                outfile.write(data)

            got = deserialise_object(outpath)
            self.assertFloatEqual(got.get_log_likelihood(), lnL)

    def test_roundtrip_model_result(self):
        """mode_result.to_json enables roundtrip and lazy evaluation"""
        _data = {
            "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
            "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
            "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
        }
        aln = LoadSeqs(data=_data, moltype="dna")
        tree = LoadTree(tip_names=aln.names)
        sm = get_model("HKY85")
        lf = sm.make_likelihood_function(tree)
        lf.set_alignment(aln)
        edge_vals = zip(aln.names, (2, 3, 4))
        for edge, val in edge_vals:
            lf.set_param_rule("kappa", edge=edge, init=val)
        result = model_result(name="test")
        result[1] = lf
        self.assertIs(result[1], lf)
        self.assertEqual(result.nfp, lf.nfp)
        self.assertEqual(result.lnL, lf.lnL)

        data = result.to_json()
        got_obj = deserialise_object(data)
        # lazy evaluation means initially, the value is a dict
        self.assertIsInstance(got_obj[1], dict)
        # and properties match original
        self.assertEqual(got_obj.lnL, result.lnL)
        self.assertEqual(got_obj.nfp, result.nfp)
        self.assertEqual(got_obj.DLC, result.DLC)
        # when we ask for the lf attribute, it's no longer a dict
        self.assertNotIsInstance(got_obj.lf, dict)
        self.assertEqual(got_obj.lf.nfp, got_obj.nfp)

    def test_not_completed_result(self):
        """correctly reconstructs a NotCompletedResult object"""
        from cogent3.app.composable import NotCompleted

        val = NotCompleted("ERROR", "nothing", "some error", source="here")
        expect = val.to_rich_dict()
        json = val.to_json()
        got = deserialise_object(json)
        self.assertEqual(got.to_rich_dict(), expect)

    def test_deserialise_tabular_table(self):
        """correctly deserialises Table"""
        from cogent3 import LoadTable

        table = LoadTable(
            header=["id", "foo", "bar"],
            rows=[
                [1, "abc", 11],
                [2, "bca", 22],
                [3, "cab", 33],
                [4, "abc", 44],
                [5, "bca", 55],
            ],
        )
        json = table.to_json()
        got = deserialise_object(json)
        self.assertEqual(got.todict(), table.todict())

    def test_deserialise_tabular_dictarray(self):
        """correctly deserialises DictArray"""
        from cogent3.util.dict_array import DictArrayTemplate

        template = DictArrayTemplate(5, ["id", "foo", "bar"])
        data = [
            [1, "abc", 11],
            [2, "bca", 22],
            [3, "cab", 33],
            [4, "abc", 44],
            [5, "bca", 55],
        ]
        darr = template.wrap(data)
        json = darr.to_json()
        got = deserialise_object(json)
        self.assertEqual(got.todict(), darr.todict())

    def test_deserialise_tabular_distancematrix(self):
        """correctly deserialises DistanceMatrix"""
        from cogent3.evolve.fast_distance import DistanceMatrix

        data = {
            ("ABAYE2984", "Atu3667"): None,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): None,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): None,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): None,
            ("BAA10469", "Atu3667"): None,
            ("BAA10469", "Avin_42730"): 1.85,
        }

        dm = DistanceMatrix(data)
        json = dm.to_json()
        got = deserialise_object(json)
        self.assertEqual(dm.todict(), got.todict())


if __name__ == "__main__":
    main()
