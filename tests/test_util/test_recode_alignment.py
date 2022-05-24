from unittest import TestCase, main

from numpy.testing import assert_allclose

from cogent3 import make_aligned_seqs
from cogent3.core.alignment import ArrayAlignment
from cogent3.evolve.models import DSO78_freqs, DSO78_matrix
from cogent3.util.recode_alignment import (
    alphabets,
    build_alphabet_map,
    recode_count_matrix,
    recode_counts_and_freqs,
    recode_dense_alignment,
    recode_freq_vector,
)


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Greg Caporaso"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Beta"


class RecodeAlignmentTests(TestCase):
    """Tests of functions in recode_alphabet.py

    These functions will probably move at some point, and the unit tests
        will move with them.
    """

    def setUp(self):
        """Initialize some variables for the tests"""
        self.canonical_abbrevs = "ACDEFGHIKLMNPQRSTVWY"
        self.ambiguous_abbrevs = "BXZ"

        self.all_to_a = [("A", self.canonical_abbrevs + self.ambiguous_abbrevs)]
        self.charge_2 = alphabets["charge_2"]
        self.hydropathy_3 = alphabets["hydropathy_3"]
        self.orig = alphabets["orig"]
        self.aln = ArrayAlignment(data={"1": "CDDFBXZ", "2": "CDD-BXZ", "3": "AAAASS-"})
        self.aln2 = make_aligned_seqs(
            data={"1": "CDDFBXZ", "2": "CDD-BXZ", "3": "AAAASS-"}
        )

    def test_build_alphabet_map_handles_bad_data(self):
        """build_alphabet_map:  bad data raises error"""
        self.assertRaises(ValueError, build_alphabet_map)
        self.assertRaises(ValueError, build_alphabet_map, "not_a_valid_id")
        self.assertRaises(
            ValueError, build_alphabet_map, alphabet_def=["A", "BCD", "B", "EFG"]
        )

    def test_build_alphabet_map_w_alphabet_id(self):
        """build_alphabet_map: returns correct dict when given alphabet_id"""
        expected = dict(
            [
                ("G", "G"),
                ("A", "G"),
                ("V", "G"),
                ("L", "G"),
                ("I", "G"),
                ("S", "G"),
                ("P", "G"),
                ("T", "G"),
                ("C", "G"),
                ("N", "G"),
                ("D", "G"),
                ("X", "G"),
                ("B", "G"),
                ("M", "M"),
                ("F", "M"),
                ("Y", "M"),
                ("W", "M"),
                ("Q", "M"),
                ("K", "M"),
                ("H", "M"),
                ("R", "M"),
                ("E", "M"),
                ("Z", "M"),
            ]
        )
        self.assertEqual(build_alphabet_map("size_2"), expected)
        self.assertEqual(build_alphabet_map("charge_3")["E"], "D")
        self.assertEqual(build_alphabet_map("charge_3")["B"], "A")
        self.assertEqual(build_alphabet_map("charge_3")["K"], "K")

    def test_build_alphabet_map_w_alphabet_def(self):
        """build_alphabet_map: returns correct dict when given alphabet_def"""
        expected = dict(
            [
                ("G", "S"),
                ("A", "S"),
                ("V", "S"),
                ("L", "S"),
                ("I", "S"),
                ("S", "S"),
                ("P", "S"),
                ("T", "S"),
                ("C", "S"),
                ("N", "S"),
                ("D", "S"),
                ("X", "S"),
                ("B", "S"),
                ("M", "L"),
                ("F", "L"),
                ("Y", "L"),
                ("W", "L"),
                ("Q", "L"),
                ("K", "L"),
                ("H", "L"),
                ("R", "L"),
                ("E", "L"),
                ("Z", "L"),
            ]
        )
        self.assertEqual(
            build_alphabet_map(
                alphabet_def=[("S", "GAVLISPTCNDXB"), ("L", "MFYWQKHREZ")]
            ),
            expected,
        )

    def test_build_alphabet_map_handles_all_ids_and_defs_wo_error(self):
        """build_alphabet_map: handles all pre-defined alphabets w/o error"""
        for alphabet_id, alphabet_def in list(alphabets.items()):
            try:
                build_alphabet_map(alphabet_id=alphabet_id)
            except ValueError:
                raise AssertionError(f"Failed on id: {alphabet_id}")
            try:
                build_alphabet_map(alphabet_def=alphabet_def)
            except ValueError:
                raise AssertionError(f"Failed on def: {str(alphabet_def)}")

    def test_recode_dense_alignment_handles_all_ids_and_defs_wo_error(self):
        """recode_dense_alignment: handles pre-defined alphabets w/o error"""
        for alphabet_id, alphabet_def in list(alphabets.items()):
            try:
                recode_dense_alignment(self.aln, alphabet_id=alphabet_id)
            except ValueError:
                raise AssertionError(f"Failed on id: {alphabet_id}")
            try:
                recode_dense_alignment(self.aln, alphabet_def=alphabet_def)
            except ValueError:
                raise AssertionError(f"Failed on def: {str(alphabet_def)}")

    def test_recode_dense_alignment_leaves_original_alignment_intact(self):
        """recode_dense_alignment: leaves input alignment intact"""
        # provided with alphabet_id
        actual = recode_dense_alignment(self.aln, alphabet_id="charge_2")
        self.assertNotEqual(actual, self.aln)
        # provided with alphabet_def
        actual = recode_dense_alignment(self.aln, alphabet_def=self.charge_2)
        self.assertNotEqual(actual, self.aln)

    def test_recode_dense_alignment(self):
        """recode_dense_alignment: recode alignment to charge_2 alpha works"""
        expected_c2 = ArrayAlignment(
            data={"1": "AKKAKAK", "2": "AKK-KAK", "3": "AAAAAA-"}
        )
        expected_h3 = ArrayAlignment(
            data={"1": "PRRPRPR", "2": "PRR-RPR", "3": "PPPPYY-"}
        )
        expected_aa = ArrayAlignment(
            data={"1": "AAAAAAA", "2": "AAA-AAA", "3": "AAAAAA-"}
        )

        # provided with alphabet_id
        actual = recode_dense_alignment(self.aln, alphabet_id="charge_2")
        self.assertEqual(actual, expected_c2)
        # provided with alphabet_def
        actual = recode_dense_alignment(self.aln, alphabet_def=self.charge_2)
        self.assertEqual(actual, expected_c2)

        # different alphabet
        actual = recode_dense_alignment(self.aln, alphabet_id="hydropathy_3")
        self.assertEqual(actual, expected_h3)
        actual = recode_dense_alignment(self.aln, alphabet_def=self.hydropathy_3)
        self.assertEqual(actual, expected_h3)

        # different alphabet
        actual = recode_dense_alignment(self.aln, alphabet_def=self.all_to_a)
        self.assertEqual(actual, expected_aa)

        # original charactars which aren't remapped are let in original state
        actual = recode_dense_alignment(self.aln, alphabet_def=[("a", "b")])
        self.assertEqual(actual, self.aln)

        # non-alphabetic character mapped same as alphabetic characters
        actual = recode_dense_alignment(self.aln, alphabet_def=[(".", "-")])
        expected = ArrayAlignment(data={"1": "CDDFBXZ", "2": "CDD.BXZ", "3": "AAAASS."})
        self.assertEqual(actual, expected)

    def test_recode_dense_alignment_to_orig(self):
        """recode_dense_alignment: recode aln to orig returns original aln"""
        # provided with alphabet_id
        self.assertEqual(recode_dense_alignment(self.aln, alphabet_id="orig"), self.aln)
        # provided with alphabet_def
        self.assertEqual(
            recode_dense_alignment(self.aln, alphabet_def=self.orig), self.aln
        )

    def test_recode_freq_vector(self):
        """recode_freq_vector: bg freqs updated to reflect recoded alphabet"""

        freqs = {"A": 0.21, "E": 0.29, "C": 0.05, "D": 0.45}
        a_def = [("A", "AEC"), ("E", "D")]
        expected = {"A": 0.55, "E": 0.45}
        self.assertEqual(recode_freq_vector(a_def, freqs), expected)
        # reversal of alphabet
        freqs = {"A": 0.21, "E": 0.29, "C": 0.05, "D": 0.45}
        a_def = [("A", "D"), ("E", "C"), ("C", "E"), ("D", "A")]
        expected = {"A": 0.45, "E": 0.05, "C": 0.29, "D": 0.21}
        self.assertEqual(recode_freq_vector(a_def, freqs), expected)

        # no change in freqs (old alphabet = new alphabet)
        freqs = {"A": 0.21, "E": 0.29, "C": 0.05, "D": 0.45}
        a_def = [("A", "A"), ("E", "E"), ("C", "C"), ("D", "D")]
        self.assertEqual(recode_freq_vector(a_def, freqs), freqs)

        freqs = {"A": 0.21, "E": 0.29, "C": 0.05, "D": 0.45}
        a_def = [("X", "AEC"), ("Y", "D")]
        expected = {"X": 0.55, "Y": 0.45}
        self.assertEqual(recode_freq_vector(a_def, freqs), expected)

    def test_recode_freq_vector_ignores(self):
        """recode_freq_vector: ignored chars are ignored"""
        freqs = {"A": 0.21, "B": 0.29, "C": 0.05, "D": 0.45, "X": 0.22, "Z": 0.5}
        a_def = [("A", "A"), ("B", "B"), ("C", "C"), ("D", "D"), ("X", "X"), ("Z", "Z")]
        expected = {"A": 0.21, "C": 0.05, "D": 0.45}
        self.assertEqual(recode_freq_vector(a_def, freqs), expected)

        freqs = {
            "D": 0.21,
            "E": 0.29,
            "N": 0.05,
            "Q": 0.45,
            "B": 0.26,
            "Z": 0.74,
            "X": 1.0,
        }
        a_def = [("D", "DEN"), ("Q", "Q")]
        expected = {"D": 0.55, "Q": 0.45}
        self.assertEqual(recode_freq_vector(a_def, freqs), expected)


class RecodeMatrixTests(TestCase):
    """Tests of substitution matrix recoding."""

    def setUp(self):
        """Create variables for use in the tests"""
        self.m1 = [
            [0, 4, 1, 3, 5],
            [4, 0, 2, 4, 6],
            [1, 2, 0, 7, 8],
            [3, 4, 7, 0, 9],
            [5, 6, 8, 9, 0],
        ]
        self.recoded_m1 = [
            [0, 0, 21, 0, 0],
            [0, 0, 0, 0, 0],
            [21, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ]
        self.aa_order1 = "DELIV"
        self.input_freqs1 = dict(list(zip(self.aa_order1, [0.2] * 5)))
        self.alphabet1 = [("D", "DE"), ("L", "LIV")]

        # create_recoded_rate_matrix(alphabets['a1_4'])
        self.m2 = [
            [0, 8, 6, 5, 1],
            [8, 0, 7, 3, 0],
            [6, 7, 0, 4, 2],
            [5, 3, 4, 0, 0],
            [1, 0, 2, 0, 0],
        ]
        self.recoded_m2 = [
            [0, 0, 21, 0, 1],
            [0, 0, 0, 0, 0],
            [21, 0, 0, 0, 2],
            [0, 0, 0, 0, 0],
            [1, 0, 2, 0, 0],
        ]
        self.aa_order2 = "DELIC"
        self.input_freqs2 = dict(list(zip(self.aa_order2, [0.2] * 5)))
        self.alphabet2 = [("D", "DE"), ("L", "LI"), ("C", "C")]
        self.alphabet2_w_ambig = [("D", "DEX"), ("L", "LIB"), ("C", "CZ")]

    def test_recode_counts_and_freqs(self):
        """recode_counts_and_freqs: functions as expected"""
        alphabet = alphabets["charge_his_3"]
        aa_order = "ACDEFGHIKLMNPQRSTVWY"
        actual = recode_counts_and_freqs(alphabet)
        expected_matrix = recode_count_matrix(
            alphabet, count_matrix=DSO78_matrix, aa_order=aa_order
        )
        expected_freqs = {}.fromkeys(aa_order, 0.0)
        expected_freqs.update(recode_freq_vector(alphabet, DSO78_freqs))
        expected = (expected_matrix, expected_freqs)
        assert_allclose(actual[0], expected[0])
        self.assertEqual(actual[1], expected[1])

    def test_recode_count_matrix_2_states(self):
        """recode_count_matrix: returns correct result with 2-state alphabet"""
        actual = recode_count_matrix(self.alphabet1, self.m1, self.aa_order1)
        expected = self.recoded_m1
        assert_allclose(actual, expected)

    def test_recode_count_matrix_3_states(self):
        """recode_count_matrix: returns correct result with 3-state alphabet"""
        actual = recode_count_matrix(self.alphabet2, self.m2, self.aa_order2)
        expected = self.recoded_m2
        assert_allclose(actual, expected)

    def test_recode_count_matrix_3_states_ambig_ignored(self):
        """recode_count_matrix: correct result w 3-state alphabet w ambig chars"""
        actual = recode_count_matrix(self.alphabet2_w_ambig, self.m2, self.aa_order2)
        expected = self.recoded_m2
        assert_allclose(actual, expected)

    def test_recode_count_matrix_no_change(self):
        """recode_count_matrix: no changes applied when they shouldn't be"""
        # recoding recoded matrices
        actual = recode_count_matrix(self.alphabet1, self.recoded_m1, self.aa_order1)
        expected = self.recoded_m1
        assert_allclose(actual, expected)

        actual = recode_count_matrix(self.alphabet2, self.recoded_m2, self.aa_order2)
        expected = self.recoded_m2
        assert_allclose(actual, expected)


if __name__ == "__main__":
    main()
