#!/usr/bin/env python

import os

from unittest import TestCase, main

from cogent3.evolve import predicate, substitution_model
from cogent3.evolve.models import F81, GN, HKY85
from cogent3.evolve.predicate import MotifChange


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, "data")


class NucleotideModelTestMethods(TestCase):
    def setUp(self):
        self.submodel = substitution_model.TimeReversibleNucleotide(model_gaps=False)

    def test_isTransition(self):
        """testing isTransition"""
        isTransition = self.submodel.get_predefined_predicate("transition")
        assert isTransition("A", "G")
        assert isTransition("C", "T")
        assert not isTransition("A", "T")
        assert not isTransition("C", "G")

        isTransition = self.submodel.get_predefined_predicate("kappa")
        assert isTransition("A", "G")
        assert isTransition("C", "T")
        assert not isTransition("A", "T")
        assert not isTransition("C", "G")

    def test_isTransversion(self):
        """testing isTransversion"""
        isTransversion = self.submodel.get_predefined_predicate("transversion")
        assert not isTransversion("A", "G")
        assert not isTransversion("C", "T")
        assert isTransversion("A", "T")
        assert isTransversion("C", "G")

    def test_isIndel(self):
        """testing indel comparison nucleotide model"""
        model = substitution_model.TimeReversibleNucleotide(model_gaps=True)
        isIndel = model.get_predefined_predicate("indel")
        assert isIndel("A", "-")
        assert isIndel("-", "G")
        # assert not self.submodel.isIndel('-', '-')
        assert not isIndel("a", "t")

    def test_PredicateChecks(self):
        # overparameterisation
        self.assertRaises(
            ValueError,
            substitution_model.TimeReversibleNucleotide,
            model_gaps=False,
            predicates=["transition", "transversion"],
        )

    def test_nonrev_exception(self):
        """constructing a Nucleotide model with non-reversible preds raises exception"""
        preds = predicate.MotifChange("A", "G", forward_only=True)
        with self.assertRaises(ValueError):
            sm = substitution_model.TimeReversibleNucleotide(predicates=[preds])

    def test_get_param_matrix_coords(self):
        """return correct coords for continuous-time models"""
        f81 = F81()
        self.assertEqual(f81.get_param_matrix_coords(), {})
        self.assertTrue(
            len(f81.get_param_matrix_coords(include_ref_cell=True)["ref_cell"]) == 12
        )
        hky85 = HKY85()
        coords = hky85.get_param_matrix_coords()
        self.assertEqual(set(coords), set(["kappa"]))
        coords = hky85.get_param_matrix_coords(include_ref_cell=True)
        self.assertEqual(set(coords), set(["kappa", "ref_cell"]))
        gn = GN()
        coords = gn.get_param_matrix_coords(include_ref_cell=True)
        self.assertTrue(len(coords) == 12)
        self.assertTrue(len(coords["ref_cell"]) == 1)

    def test_to_rich_dict(self):
        """returns complete dict of attributes"""
        F81().to_rich_dict()
        HKY85().to_rich_dict()
        GN().to_rich_dict()
        # TODO need to assess ability to reconstruct from this


class MultiLetterMotifSubstModelTests(TestCase):
    def setUp(self):
        self.submodel = substitution_model.TimeReversibleDinucleotide(
            model_gaps=True, mprob_model="tuple"
        )

    def test_ascii_art(self):
        model = substitution_model.TimeReversibleDinucleotide(
            mprob_model="tuple", predicates=["k:transition"]
        )
        model.ascii_art()
        model = substitution_model.TimeReversibleDinucleotide(mprob_model="tuple")
        model.ascii_art()

    def test_isIndel(self):
        """testing indel comparison for dinucleotide model"""
        # these are non-instantaneous
        isIndel = self.submodel.get_predefined_predicate("indel")
        assert not isIndel("AA", "--")
        assert not isIndel("--", "CT")

        # assert not self.submodel.isIndel('--', '--')
        assert not isIndel("AT", "AA")

        assert isIndel("AA", "A-")
        assert isIndel("-A", "CA")
        assert isIndel("TA", "-A")

        # isIndel can now assume it won't get any non-instantaneous pairs
        # assert self.submodel.isIndel('-a', 'a-') == 0

    def test_to_rich_dict(self):
        """returns complete dict of attributes"""
        got = self.submodel.to_rich_dict()
        self.assertEqual(got["motif_length"], 2)


nuc_probs = [("T", 0.1), ("C", 0.2), ("A", 0.3), ("G", 0.4)]


class TupleModelMotifProbFuncs(TestCase):
    dinucs = (
        "TT",
        "CT",
        "AT",
        "GT",
        "TC",
        "CC",
        "AC",
        "GC",
        "TA",
        "CA",
        "AA",
        "GA",
        "TG",
        "CG",
        "AG",
        "GG",
    )
    nuc_probs = nuc_probs
    dinuc_probs = [(m2 + m1, p1 * p2) for m1, p1 in nuc_probs for m2, p2 in nuc_probs]
    mat_indices = dict(
        C=set(
            [
                (0, 1),
                (0, 4),
                (1, 5),
                (2, 1),
                (2, 6),
                (3, 1),
                (3, 7),
                (4, 5),
                (6, 5),
                (7, 5),
                (8, 4),
                (8, 9),
                (9, 5),
                (10, 6),
                (10, 9),
                (11, 7),
                (11, 9),
                (12, 4),
                (12, 13),
                (13, 5),
                (14, 6),
                (14, 13),
                (15, 7),
                (15, 13),
            ]
        ),
        A=set(
            [
                (0, 2),
                (0, 8),
                (1, 2),
                (1, 9),
                (2, 10),
                (3, 2),
                (3, 11),
                (4, 6),
                (4, 8),
                (5, 6),
                (5, 9),
                (6, 10),
                (7, 6),
                (7, 11),
                (8, 10),
                (9, 10),
                (11, 10),
                (12, 8),
                (12, 14),
                (13, 9),
                (13, 14),
                (14, 10),
                (15, 11),
                (15, 14),
            ]
        ),
        G=set(
            [
                (0, 3),
                (0, 12),
                (1, 3),
                (1, 13),
                (2, 3),
                (2, 14),
                (3, 15),
                (4, 7),
                (4, 12),
                (5, 7),
                (5, 13),
                (6, 7),
                (6, 14),
                (7, 15),
                (8, 11),
                (8, 12),
                (9, 11),
                (9, 13),
                (10, 11),
                (10, 14),
                (11, 15),
                (12, 15),
                (13, 15),
                (14, 15),
            ]
        ),
        T=set(
            [
                (1, 0),
                (2, 0),
                (3, 0),
                (4, 0),
                (5, 1),
                (5, 4),
                (6, 2),
                (6, 4),
                (7, 3),
                (7, 4),
                (8, 0),
                (9, 1),
                (9, 8),
                (10, 2),
                (10, 8),
                (11, 3),
                (11, 8),
                (12, 0),
                (13, 1),
                (13, 12),
                (14, 2),
                (14, 12),
                (15, 3),
                (15, 12),
            ]
        ),
    )


class ThreeLetterMotifSubstModelTests(TestCase):
    def setUp(self):
        self.submodel = substitution_model.TimeReversibleNucleotide(
            motif_length=3, mprob_model="tuple"
        )

    def test_isIndel(self):
        """testing indel comparison for trinucleotide model"""
        isIndel = self.submodel.get_predefined_predicate("indel")
        assert isIndel("AAA", "AA-")
        assert isIndel("-CA", "CCA")
        assert isIndel("TAC", "T-C")

        # isIndel can now assume it won't get any non-instantaneous pairs
        assert not isIndel("AAA", "---")
        assert not isIndel("---", "CTT")
        assert not isIndel("AAA", "--A")
        assert not isIndel("C--", "CTT")

    def test_to_rich_dict(self):
        """returns complete dict of attributes"""
        got = self.submodel.to_rich_dict()
        self.assertEqual(got["motif_length"], 3)

    def test_nr_trinuc(self):
        """This is exercising a TimeReversibleTriucleotide"""
        preds = [
            MotifChange("A", "C"),
            MotifChange("G", "A"),
            MotifChange("CGA", "TGA"),
        ]
        sm = substitution_model.TimeReversibleTrinucleotide(predicates=preds)
        got = sm.get_param_list()
        self.assertEqual(got, ["A/C", "G/A", "CGA/TGA"])
        self.assertEqual(len(sm.get_motifs()), 64)


class CodonSubstModelTests(TestCase):
    def setUp(self):
        self.standardcode = substitution_model.TimeReversibleCodon(
            model_gaps=True, gc=1, mprob_model="tuple"
        )
        self.mitocode = substitution_model.TimeReversibleCodon(
            model_gaps=False, gc=2, mprob_model="tuple"
        )

    def test_isTransition(self):
        """testing codon isTransition"""
        isTransition = self.standardcode.get_predefined_predicate("transition")
        # first position
        assert isTransition("TGC", "CGC")
        assert isTransition("GGC", "AGC")
        # second position
        assert isTransition("CTT", "CCT")
        assert isTransition("CAT", "CGT")
        # thirs position
        assert isTransition("CTT", "CTC")
        assert isTransition("CTA", "CTG")
        # mito code
        assert isTransition("CTT", "CTC")
        assert isTransition("CTA", "CTG")

        assert not isTransition("GAG", "GTG")
        assert not isTransition("CCC", "CGC")

    def test_isReplacement(self):
        """test isReplacement for the two major genetic codes"""
        isReplacement = self.standardcode.get_predefined_predicate("replacement")
        # for the standard code, a replacement
        assert isReplacement("CTG", "ATG")
        assert not isReplacement("AGT", "TCC")
        assert not isReplacement("CTG", "---")
        assert not isReplacement("---", "CTA")
        # for the vert mitocho code, instantaneous replacement
        isReplacement = self.mitocode.get_predefined_predicate("replacement")
        assert isReplacement("AAA", "AAC")

        # check using 'omega'
        isReplacement = self.standardcode.get_predefined_predicate("omega")
        assert isReplacement("CTG", "ATG")
        assert not isReplacement("AGT", "TCC")
        assert not isReplacement("CTG", "---")
        assert not isReplacement("---", "CTA")
        isReplacement = self.mitocode.get_predefined_predicate("omega")
        assert isReplacement("AAA", "AAC")

    def test_isSilent(self):
        """testing isSilent for the two major genetic codes"""
        isSilent = self.standardcode.get_predefined_predicate("silent")
        assert isSilent("CTA", "CTG")
        assert not isSilent("AGT", "AAG")
        assert not isSilent("CTG", "---")
        assert not isSilent("---", "CTG")
        # for vert mito code
        isSilent = self.mitocode.get_predefined_predicate("silent")
        assert isSilent("TCC", "TCA")

    def test_isIndel(self):
        """test isIndel for codon model"""
        isIndel = self.standardcode.get_predefined_predicate("indel")
        assert isIndel("CTA", "---")
        assert not isIndel("---", "---")
        assert isIndel("---", "TTC")

    def test_str_(self):
        """str() and repr() of a substitution model"""
        str(self.standardcode)
        repr(self.standardcode)


class ModelDataInteractionTestMethods(TestCase):
    def test_excludeinggaps(self):
        """testing excluding gaps from model"""
        model = substitution_model.TimeReversibleNucleotide(model_gaps=False)
        assert len(model.get_alphabet()) == 4

    def test_includinggaps(self):
        """testing excluding gaps from model"""
        model = substitution_model.TimeReversibleNucleotide(model_gaps=True)
        assert len(model.get_alphabet()) == 5

    def test_getMotifs(self):
        """testing return of motifs"""
        substitution_model.TimeReversibleNucleotide().get_motifs()

    def test_get_param_list(self):
        """testing getting the parameter list"""
        model = substitution_model.TimeReversibleNucleotide()
        self.assertEqual(model.get_param_list(), [])

        model = substitution_model.TimeReversibleNucleotide(
            predicates=["beta:transition"]
        )
        self.assertEqual(model.get_param_list(), ["beta"])

    # need to ensure entering motif probs that sum to 1, that motif sets are
    # the same


if __name__ == "__main__":
    main()
