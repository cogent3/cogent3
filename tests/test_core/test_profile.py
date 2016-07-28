#!/usr/bin/env python
"""Provides tests for classes and functions in profile.py
"""
from numpy import array, sum, sqrt, transpose, add, subtract, multiply,\
    divide, zeros
from numpy.random import random

from cogent3.util.unit_test import TestCase, main  # , numpy_err
from cogent3.core.moltype import DNA
from cogent3.core.sequence import ModelSequence
from cogent3.core.profile import Profile, ProfileError, CharMeaningProfile
from cogent3.core.alignment import DenseAlignment as Alignment

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley", "Rob Knight",
               "Peter Maxwell"]
__license__ = "GPL"
__version__ = "3.0.prealpha"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

translate = str.translate


class ProfileTests(TestCase):
    """Tests for Profile object"""

    def setUp(self):
        """setUp method for all Profile tests"""
        self.full = Profile(array([[2, 4], [3, 5], [4, 8]]), "AB")
        self.empty = Profile(array([[]]), "AB")
        self.empty_row = Profile(array([[1, 1], [0, 0]]), "AB")
        self.empty_col = Profile(array([[0, 1], [0, 1]]), "AB")
        self.consensus = Profile(array([[.2, 0, .8, 0], [0, .1, .2, .7], [0, 0, 0, 1],
                                        [.2, .3, .4, .1], [.5, .5, 0, 0]]),
                                 alphabet=DNA, CharOrder="TCAG")
        self.not_same_value = Profile(array([[.3, .5, .1, .1], [.4, .6, 0, .7],
                                             [.3, .2, 0, 0], [0, 0, 4, 0]]), alphabet=DNA, CharOrder="TCAG")
        self.zero_entry = Profile(array([[.3, .2, 0, .5], [0, 0, .8, .2]]),
                                  alphabet="UCAG")
        self.score1 = Profile(Data=array([[-1, 0, 1, 2], [-2, 2, 0, 0], [-3, 5, 1, 0]]),
                              alphabet=DNA, CharOrder="ATGC")
        self.score2 = Profile(array([[.2, .4, .4, 0], [.1, 0, .9, 0], [.1, .2, .3, .4]]),
                              alphabet="TCAG")
        self.oned = Profile(array([.25, .25, .25, .25]), "ABCD")
        self.pp = Profile(
            array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]), "ABCD")

    def test_init(self):
        """__init__: should set all attributed correctly"""
        self.assertRaises(TypeError, Profile)
        self.assertRaises(TypeError, Profile, array([[2, 3]]))
        # only alphabet
        p = Profile(array([[.2, .8], [.7, .3]]), "AB")
        self.assertEqual(p.Data, [[.2, .8], [.7, .3]])
        self.assertEqual(p.alphabet, "AB")
        self.assertEqual(p.CharOrder, list("AB"))
        self.assertEqual(translate("ABBA", p._translation_table),
                         "\x00\x01\x01\x00")
        # alphabet and char order
        p = Profile(array([[.1, .2], [.4, .3]]), alphabet=DNA,
                    CharOrder="AG")
        self.assertEqual(p.CharOrder, "AG")
        assert p.alphabet is DNA
        # non-character alphabet
        p = Profile(array([[.1, .2], [.4, .3]]), alphabet=[7, 3],
                    CharOrder=[3, 7])
        self.assertEqual(p.CharOrder, [3, 7])
        self.assertEqual(p.alphabet, [7, 3])
        self.assertEqual(p.Data, [[.1, .2], [.4, .3]])

    def test_str(self):
        """__str__: should return string representation of data in profile
        """
        self.assertEqual(str(self.empty_row), str(array([[1, 1], [0, 0]])))

    def test_make_translation_table(self):
        """_make_translation_table: should return correct table from char order
        """
        p = Profile(array([[.2, .8], [.7, .3]]), "ABCDE", "AB")
        self.assertEqual(translate("ABBA", p._translation_table),
                         "\x00\x01\x01\x00")

    def test_has_valid_data(self):
        """has_valid_data: should work on full and empty profiles"""
        full = self.full.copy()
        full.normalize_positions()
        self.assertEqual(full.has_valid_data(), True)
        self.assertEqual(self.empty_row.has_valid_data(), False)
        self.assertEqual(self.empty.has_valid_data(), False)

    def test_has_valid_attributes(self):
        """has_valid_attributes: should work for different alphabets/char orders
        """
        p = Profile(array([[1, 2], [3, 4]]), alphabet="ABCD", CharOrder="BAC")
        # self.Data doesn't match len(CharOrder)
        self.assertEqual(p.has_valid_attributes(), False)
        p = Profile(array([[1, 2], [3, 4]]), alphabet="ABCD", CharOrder="AX")
        # not all chars in CharOrder in Alphabet
        self.assertEqual(p.has_valid_attributes(), False)
        p = Profile(array([[1, 2], [3, 4]]), alphabet="ABCD", CharOrder="CB")
        # should be fine
        self.assertEqual(p.has_valid_attributes(), True)

    def test_is_valid(self):
        """is_valid: should work as expected"""
        # everything valid
        p1 = Profile(array([[.3, .7], [.8, .2]]),
                     alphabet="AB", CharOrder="AB")
        # invalid data, valid attributes
        p2 = Profile(array([[1, 2], [3, 4]]), alphabet="ABCD", CharOrder="BA")
        # invalid attributes, valid data
        p3 = Profile(array([[.3, .7], [.8, .2]]),
                     alphabet="ABCD", CharOrder="AF")

        self.assertEqual(p1.is_valid(), True)
        self.assertEqual(p2.is_valid(), False)
        self.assertEqual(p3.is_valid(), False)

    def test_data_at(self):
        """data_at: should work on valid position and character"""
        p = Profile(array([[.2, .4, .4, 0], [.1, 0, .9, 0], [.1, .2, .3, .4]]),
                    alphabet="TCAG")
        self.assertEqual(p.data_at(0, 'C'), .4)
        self.assertEqual(p.data_at(1, 'T'), .1)
        self.assertRaises(ProfileError, p.data_at, 1, 'U')
        self.assertRaises(ProfileError, p.data_at, -2, 'T')
        self.assertRaises(ProfileError, p.data_at, 5, 'T')

    def test_copy(self):
        """copy: should act as expected while rebinding/modifying attributes
        """
        p = Profile(array([[1, 1], [.7, .3]]), {
                    'A': 'A', 'G': 'G', 'R': 'AG'}, "AG")
        p_copy = p.copy()
        assert p.Data is p_copy.Data
        assert p.alphabet is p_copy.alphabet
        assert p.CharOrder is p_copy.CharOrder

        # modifying p.Data modifies p_copy.Data
        p.Data[1, 1] = 100
        assert p.alphabet is p_copy.alphabet

        # normalizing p.Data rebinds it, so p_copy.Data is unchanged
        p.normalize_positions()
        assert not p.Data is p_copy.Data

        # Adding something to the alphabet changes both p and p_copy
        p.alphabet['Y'] = 'TC'
        assert p.alphabet is p_copy.alphabet

        # Rebinding the CharOrder does only change the original
        p.CharOrder = 'XX'
        assert not p.CharOrder is p_copy.CharOrder

    def test_normalizePositions(self):
        """normalizePositions: should normalize or raise appropriate error
        """
        p = self.full.copy()
        p.normalize_positions()
        self.assertEqual(p.Data, array(
            [[2 / 6, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]))
        self.assertEqual(sum(p.Data, 1), [1, 1, 1])
        p = self.empty_col.copy()
        p.normalize_positions()
        self.assertEqual(p.Data, array([[0, 1], [0, 1]]))
        p = self.empty_row.copy()
        self.assertRaises(ProfileError, p.normalize_positions)
        p = Profile(array([[0.0, 0.0]]), "AB")
        self.assertRaises(ProfileError, p.normalize_positions)

        # negative numbers!!!!!!
        p1 = Profile(array([[3, -2], [4, -3]]), "AB")
        p1.normalize_positions()
        self.assertEqual(p1.Data, array([[3, -2], [4, -3]]))
        p2 = Profile(array([[3, -3], [4, -3]]), "AB")
        self.assertRaises(ProfileError, p2.normalize_positions)

    def test_normalize_seqs(self):
        """normalize_seqs: should normalize or raise appropriate error
        """
        p = self.full.copy()
        p.normalize_seqs()
        self.assertEqual(p.Data, array(
            [[2 / 9, 4 / 17], [3 / 9, 5 / 17], [4 / 9, 8 / 17]]))
        self.assertEqual(sum(p.Data, axis=0), [1, 1])
        p = self.empty_row.copy()
        p.normalize_seqs()
        self.assertEqual(p.Data, array([[1, 1], [0, 0]]))
        p = self.empty_col.copy()
        self.assertRaises(ProfileError, p.normalize_seqs)
        p = Profile(array([[0.0], [0.0]]), "AB")
        self.assertRaises(ProfileError, p.normalize_seqs)

        # negative numbers!!!!!!
        p1 = Profile(array([[3, 4], [-2, -3]]), "AB")
        p1.normalize_seqs()
        self.assertEqual(p1.Data, array([[3, 4], [-2, -3]]))
        p2 = Profile(array([[3, 4], [-3, -3]]), "AB")
        self.assertRaises(ProfileError, p2.normalize_seqs)

    def test_pretty_print_without_parameters(self):
        """pretty_print: should work without parameters passed in"""
        p = self.full
        self.assertEqual(p.pretty_print(), "2\t4\n3\t5\n4\t8")
        self.assertEqual(p.pretty_print(include_header=True),
                         "A\tB\n2\t4\n3\t5\n4\t8")
        self.assertEqual(p.pretty_print(transpose_data=True),
                         "2\t3\t4\n4\t5\t8")
        self.assertEqual(p.pretty_print(include_header=True,
                                       transpose_data=True), "A\t2\t3\t4\nB\t4\t5\t8")
        # empty
        self.assertEqual(self.empty.pretty_print(), "")
        self.assertEqual(self.empty.pretty_print(transpose_data=True), "")

        # it will still print with invalid data (e.g if len(CharOrder)
        # doesn't match the data
        p = self.full.copy()
        p.CharOrder = "ABC"

        self.assertEqual(p.pretty_print(include_header=True),
                         "A\tB\tC\n2\t4\t \n3\t5\t \n4\t8\t ")
        # it will truncate the CharOrder if data is transposed
        # and CharOrder is longer then the number of rows in the
        # transposed data
        self.assertEqual(p.pretty_print(include_header=True,
                                       transpose_data=True), "A\t2\t3\t4\nB\t4\t5\t8")

    def test_pretty_print_four_cases(self):
        """pretty_print: with/without header/transpose/limit"""
        p = self.full
        p = self.pp
        self.assertEqual(p.pretty_print(),
                         "1\t 2\t 3\t 4\n5\t 6\t 7\t 8\n9\t10\t11\t12")
        self.assertEqual(p.pretty_print(column_limit=3),
                         "1\t 2\t 3\n5\t 6\t 7\n9\t10\t11")
        self.assertEqual(p.pretty_print(column_limit=3, include_header=True),
                         "A\t B\t C\n1\t 2\t 3\n5\t 6\t 7\n9\t10\t11")
        self.assertEqual(p.pretty_print(column_limit=3, include_header=False,
                                       transpose_data=True),
                         "1\t5\t 9\n2\t6\t10\n3\t7\t11\n4\t8\t12")
        self.assertEqual(p.pretty_print(column_limit=2, include_header=False,
                                       transpose_data=True),
                         "1\t5\n2\t6\n3\t7\n4\t8")
        self.assertEqual(p.pretty_print(column_limit=3, include_header=True,
                                       transpose_data=True),
                         "A\t1\t5\nB\t2\t6\nC\t3\t7\nD\t4\t8")

    def test_reduce_wrong_size(self):
        """reduce: should fail when profiles have different sizes"""
        p1 = Profile(array([[1, 0], [0, 1]]), alphabet="AB")
        p2 = Profile(array([[1, 0, 0], [1, 0, 0]]), alphabet="ABC")
        self.assertRaises(ProfileError, p1.reduce, p2)

    def test_reduce_normalization_error(self):
        """reduce: fails when input or output can't be normalized"""
        # Will raise errors when input data can't be normalized
        self.assertRaises(ProfileError, self.empty.reduce, self.empty, add)
        self.assertRaises(ProfileError, self.full.reduce, self.empty_row, add)

        # don't normalize input, but do normalize output
        # fails when one row adds up to zero
        p1 = Profile(array([[3, 3], [4, 4]]), "AB")
        p2 = Profile(array([[3, 3], [-4, -4]]), "AB")
        self.assertRaises(ProfileError, p1.reduce, p2, add, False, True)

    def test_reduce_operators(self):
        """reduce: should work fine with different operators
        """
        # different operators, normalize input, don't normalize output
        p1 = Profile(array([[1, 0, 0], [0, 1, 0]]), alphabet="ABC")
        p2 = Profile(array([[1, 0, 0], [0, 0, 1]]), alphabet="ABC")

        self.assertEqual(p1.reduce(p2).Data, array([[1, 0, 0], [0, .5, .5]]))
        self.assertEqual(p1.reduce(p2, add, normalize_input=True,
                                   normalize_output=False).Data, array([[2, 0, 0], [0, 1, 1]]))
        self.assertEqual(p1.reduce(p2, subtract, normalize_input=True,
                                   normalize_output=False).Data, array([[0, 0, 0], [0, 1, -1]]))
        self.assertEqual(p1.reduce(p2, multiply, normalize_input=True,
                                   normalize_output=False).Data, array([[1, 0, 0], [0, 0, 0]]))

        self.assertRaises(ProfileError, p1.reduce, p2, divide,
                          normalize_input=True, normalize_output=False)

        # don't normalize and normalize only input
        p3 = Profile(array([[1, 2], [3, 4]]), alphabet="AB")
        p4 = Profile(array([[4, 3], [2, 1]]), alphabet="AB")

        self.assertEqual(p3.reduce(p4, add, normalize_input=False,
                                   normalize_output=False).Data, array([[5, 5], [5, 5]]))
        self.assertFloatEqual(p3.reduce(p4, add, normalize_input=True,
                                        normalize_output=False).Data, array([[19 / 21, 23 / 21], [23 / 21, 19 / 21]]))

        # normalize input and output
        p5 = Profile(array([[1, 1, 0, 0], [1, 1, 1, 1]]), alphabet="ABCD")
        p6 = Profile(array([[1, 0, 0, 0], [1, 0, 0, 1]]), alphabet="ABCD")

        self.assertEqual(p5.reduce(p6, add, normalize_input=True,
                                   normalize_output=True).Data, array([[.75, .25, 0, 0],
                                                                      [.375, .125, .125, .375]]))

        # it can collapse empty profiles when normalizing is turned off
        self.assertEqual(self.empty.reduce(self.empty,
                                           normalize_input=False, normalize_output=False).Data.tolist(), [[]])

        # more specific tests of the operators will be in the
        # separate functions

    def test__add_(self):
        """__add__: should not normalize input or output, just add"""
        p1 = Profile(
            array([[.3, .4, .1, 0], [.1, .1, .1, .7]]), alphabet="ABCD")
        p2 = Profile(array([[1, 0, 0, 0], [1, 0, 0, 1]]), alphabet="ABCD")
        self.assertEqual(
            (p1 + p2).Data, array([[1.3, .4, .1, 0], [1.1, .1, .1, 1.7]]))
        self.assertRaises(ProfileError, self.empty.__add__, p1)
        self.assertEqual((self.empty + self.empty).Data.tolist(), [[]])

    def test__sub_(self):
        """__sub__: should subtract two profiles, no normalization"""
        p1 = Profile(
            array([[.3, .4, .1, 0], [.1, .1, .1, .7]]), alphabet="ABCD")
        p2 = Profile(array([[1, 0, 0, 0], [1, 0, 0, 1]]), alphabet="ABCD")
        self.assertFloatEqual((p1 - p2).Data, array([[-.7, .4, .1, 0],
                                                   [-.9, .1, .1, -.3]]))

    def test__mul_(self):
        """__mul__: should multiply two profiles, no normalization"""
        p1 = Profile(array([[1, -2, 3, 0], [1, 1, 1, .5]]), alphabet="ABCD")
        p2 = Profile(array([[1, 0, 0, 0], [1, 0, 3, 2]]), alphabet="ABCD")
        self.assertEqual((p1 * p2).Data, array([[1, 0, 0, 0],
                                              [1, 0, 3, 1]]))

    def test__div_(self):
        """__div__ and __truediv__: always true division b/c __future__.division
        """
        p1 = Profile(array([[2, 3], [4, 5]]), "AB")
        p2 = Profile(array([[1, 0], [4, 5]]), "AB")  # Int 0
        p3 = Profile(array([[1, 0.0], [4, 5]]), "AB")  # Float 0.0
        p4 = Profile(array([[1, 2], [8.0, 5]]), "AB")  # Float 0.0

        self.assertRaises(ProfileError, p1.__truediv__, p2)
        # infinity in result data
        self.assertRaises(ProfileError, p1.__div__, p3)
        self.assertFloatEqual((p1.__div__(p4)).Data,
                              array([[2, 1.5], [0.5, 1]]))

    def test_distance(self):
        """distance: should return correct distance between the profiles
        """
        p1 = Profile(array([[2, 4], [3, 1]]), "AB")
        p2 = Profile(array([[4, 6], [5, 3]]), "AB")
        p3 = Profile(array([[4, 6], [5, 3], [1, 1]]), "AB")
        p4 = Profile(array([2, 2]), "AB")
        p5 = Profile(array([2, 2, 2]), "AB")
        p6 = Profile(array([[]]), "AB")

        self.assertEqual(p1.distance(p2), 4)
        self.assertEqual(p2.distance(p1), 4)
        self.assertEqual(p1.distance(p4), sqrt(6))
        self.assertEqual(p6.distance(p6), 0)

        # Raises error when frames are not aligned
        self.assertRaises(ProfileError, p1.distance, p3)
        self.assertRaises(ProfileError, p1.distance, p5)

    def test_to_odds_matrix(self):
        """to_odds_matrix: should work on valid data or raise an error
        """
        p = Profile(array([[.1, .3, .5, .1], [.25, .25, .25, .25],
                           [.05, .8, .05, .1], [.7, .1, .1, .1], [.6, .15, .05, .2]]),
                    alphabet="ACTG")
        p_exp = Profile(array([[.4, 1.2, 2, .4], [1, 1, 1, 1], [.2, 3.2, .2, .4],
                               [2.8, .4, .4, .4], [2.4, .6, .2, .8]]), alphabet="ACTG")
        self.assertEqual(p.to_odds_matrix().Data, p_exp.Data)
        assert p.alphabet is p.to_odds_matrix().alphabet
        self.assertEqual(p.to_odds_matrix([.25, .25, .25, .25]).Data, p_exp.Data)

        # fails if symbol_freqs has wrong size
        self.assertRaises(ProfileError, p.to_odds_matrix,
                          [.25, .25, .25, .25, .25, .25])
        self.assertRaises(ProfileError, self.zero_entry.to_odds_matrix,
                          [.1, .2, .3])
        # works on empty profile
        self.assertEqual(self.empty.to_odds_matrix().Data.tolist(), [[]])
        # works with different input
        self.assertEqual(self.zero_entry.to_odds_matrix().Data,
                         array([[1.2, .8, 0, 2], [0, 0, 3.2, .8]]))
        self.assertFloatEqual(self.zero_entry.to_odds_matrix([.1, .2, .3, .4]).Data,
                              array([[3, 1, 0, 1.25], [0, 0, 2.667, .5]]), 1e-3)
        # fails when one of the background frequencies is 0
        self.assertRaises(ProfileError, self.zero_entry.to_odds_matrix,
                          [.1, .2, .3, 0])

    def test_to_log_odds_matrix(self):
        """to_log_odds_matrix: should work as expected"""
        # This test can be short, because it mainly depends on to_odds_matrix
        # for which everything has been tested
        p = Profile(array([[.1, .3, .5, .1], [.25, .25, .25, .25],
                           [.05, .8, .05, .1], [.7, .1, .1, .1], [.6, .15, .05, .2]]),
                    alphabet="ACTG")
        p_exp = Profile(array(
            [[-1.322, 0.263, 1., -1.322],
             [ 0., 0., 0., 0.],
             [-2.322,  1.678, -2.322, -1.322],
             [ 1.485, -1.322, -1.322, -1.322],
             [ 1.263, -0.737, -2.322, -0.322]]),
            alphabet="ACTG")
        self.assertFloatEqual(p.to_log_odds_matrix().Data, p_exp.Data, eps=1e-3)
        # works on empty matrix
        self.assertEqual(self.empty.to_log_odds_matrix().Data.tolist(), [[]])

    def test__score_indices(self):
        """_score_indices: should work on valid input"""
        self.assertEqual(self.score1._score_indices(array([0, 1, 1, 3, 0, 3]),
                                                    offset=0), [6, 2, -3, 0])
        self.assertFloatEqual(self.score2._score_indices(
            array([3, 1, 2, 0, 2, 2, 3]), offset=0), [.3, 1.4, .8, 1.4, 1.7])
        self.assertFloatEqual(self.score2._score_indices(
            array([3, 1, 2, 0, 2, 2, 3]), offset=3), [1.4, 1.7])
        # Errors will be raised on invalid input. Errors are not handled
        # in this method. Validation of the input is done elsewhere
        self.assertRaises(IndexError, self.score2._score_indices,
                          array([3, 1, 63, 0, 4, 2, 3]), offset=3)

    def test__score_profile(self):
        """_score_profile: should work on valid input"""
        p1 = Profile(array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, .5, .5], [0, 0, 0, 1],
                            [.25, .25, .25, .25]]), "TCAG")
        p2 = Profile(array([[0, 1, 0, 0], [.2, 0, .8, 0], [0, 0, .5, .5], [1 / 3, 1 / 3, 0, 1 / 3],
                            [.25, .25, .25, .25]]), "TCAG")

        self.assertFloatEqual(self.score2._score_profile(p1, offset=0),
                              [.55, 1.25, .45])
        self.assertFloatEqual(self.score2._score_profile(p1, offset=2),
                              [.45])
        self.assertFloatEqual(self.score2._score_profile(p2, offset=0),
                              [1.49, 1.043, .483], 1e-3)

        # Errors will be raised on invalid input. Errors are not handled
        # in this method. Validation of the input is done elsewhere
        # In this case you don't get an error, but for sure an unexpected
        # result
        self.assertFloatEqual(self.score2._score_profile(p1, offset=3).tolist(),
                              [])

    def test_score_sequence(self):
        """score: should work correctly for Sequence as input
        """
       # works on normal valid data
        s1 = self.score1.score("ATTCAC", offset=0)
        self.assertEqual(s1,
                         [6, 2, -3, 0])
        self.assertFloatEqual(self.score2.score("TCAAGT", offset=0),
                              [.5, 1.6, 1.7, 0.5])
        # works with different offset
        self.assertFloatEqual(self.score2.score("TCAAGT", offset=2),
                              [1.7, 0.5])
        self.assertFloatEqual(self.score2.score("TCAAGT", offset=3),
                              [0.5])
        # raises error on invalid offset
        self.assertRaises(ProfileError, self.score2.score,
                          "TCAAGT", offset=4)
        # works on seq of minimal length
        self.assertFloatEqual(self.score2.score("AGT", offset=0),
                              [0.5])
        # raises error when sequence is too short
        self.assertRaises(ProfileError, self.score2.score, "", offset=0)
        # raises error on empty profile
        self.assertRaises(ProfileError, self.empty.score, "ACGT")
        # raises error when sequence contains characters that
        # are not in the characterorder
        self.assertRaises(ProfileError, self.score2.score, "ACBRT")

    def test_score_sequence_object(self):
        """score: should work correctly on Sequence object as input
        """
        # DnaSequence object
        ds = self.score1.score(DNA.make_sequence("ATTCAC"), offset=0)
        self.assertEqual(ds, [6, 2, -3, 0])
        # ModelSequence object
        ms = self.score1.score(ModelSequence("ATTCAC", alphabet=DNA.alphabet),
                               offset=0)
        self.assertEqual(ms, [6, 2, -3, 0])

    def test_score_no_trans_table(self):
        """score: should work when no translation table is present
        """
        p = Profile(Data=array([[-1, 0, 1, 2], [-2, 2, 0, 0], [-3, 5, 1, 0]]),
                    alphabet=DNA, CharOrder="ATGC")
        # remove translation table
        del p.__dict__['_translation_table']
        # then score the profile
        s1 = p.score(DNA.make_sequence("ATTCAC"), offset=0)
        self.assertEqual(s1, [6, 2, -3, 0])

    def test_score_profile(self):
        """score: should work correctly for Profile as input
        """
        p1 = Profile(array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, .5, .5], [0, 0, 0, 1],
                            [.25, .25, .25, .25]]), "TCAG")
        p2 = Profile(array([[0, 1, 0, 0], [.2, 0, .8, 0], [0, 0, .5, .5], [1 / 3, 1 / 3, 0, 1 / 3],
                            [.25, .25, .25, .25]]), "TCAG")
        p3 = Profile(array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]), "TCAG")
        p4 = Profile(array([[1, 0, 0, 0], [0, 1, 0, 0]]), "TCAG")
        p5 = Profile(array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]), "AGTC")

        # works on normal valid data
        self.assertFloatEqual(self.score2.score(p1, offset=0),
                              [.55, 1.25, .45])
        self.assertFloatEqual(self.score2.score(p2, offset=0),
                              [1.49, 1.043, .483], 1e-3)
        # works with different offset
        self.assertFloatEqual(self.score2.score(p1, offset=1),
                              [1.25, 0.45])
        self.assertFloatEqual(self.score2.score(p1, offset=2),
                              [0.45])
        # raises error on invalid offset
        self.assertRaises(ProfileError, self.score2.score,
                          p1, offset=3)
        # works on profile of minimal length
        self.assertFloatEqual(self.score2.score(p3, offset=0),
                              [0.6])
        # raises error when profile is too short
        self.assertRaises(ProfileError, self.score2.score, p4, offset=0)
        # raises error on empty profile
        self.assertRaises(ProfileError, self.empty.score, p1)
        # raises error when character order doesn't match
        self.assertRaises(ProfileError, self.score2.score, p5)

    def test_row_uncertainty(self):
        """row_uncertainty: should handle full and empty profiles
        """
        p = Profile(array([[.25, .25, .25, .25], [.5, .5, 0, 0]]), "ABCD")
        self.assertEqual(p.row_uncertainty(), [2, 1])

        # for empty rows 0 is returned as the uncertainty
        self.assertEqual(self.empty.row_uncertainty().tolist(), [])
        p = Profile(array([[], [], []]), "")
        self.assertEqual(p.row_uncertainty().tolist(), [])
        # doesn't work on 1D array
        self.assertRaises(ProfileError, self.oned.row_uncertainty)

    def test_columnUncertainty(self):
        """columnUncertainty: should handle full and empty profiles
        """
        p = Profile(array([[.25, .5], [.25, .5], [.25, 0], [.25, 0]]), "AB")
        self.assertEqual(p.column_uncertainty(), [2, 1])
        # for empty cols nothing is returned as the uncertainty
        self.assertEqual(self.empty.column_uncertainty().tolist(), [])
        p = Profile(array([[], [], []]), "")
        self.assertEqual(p.column_uncertainty().tolist(), [])
        # doesn't work on 1D array
        self.assertRaises(ProfileError, self.oned.column_uncertainty)

    def test_row_degeneracy(self):
        """rowDegneracy: should work as expected"""
        p1 = self.consensus
        p2 = self.not_same_value

        self.assertEqual(p1.row_degeneracy(), [1, 1, 1, 2, 1])
        self.assertEqual(p1.row_degeneracy(cutoff=.5), [1, 1, 1, 2, 1])
        self.assertEqual(p1.row_degeneracy(cutoff=.75), [1, 2, 1, 3, 2])
        # when a row seems to add up to the cutoff value, it's not
        # always found because of floating point error. E.g. second row
        # in this example
        self.assertEqual(p1.row_degeneracy(cutoff=1), [2, 4, 1, 4, 2])
        # when the cutoff can't be found, the number of columns in the
        # profile is returned (for each row)
        self.assertEqual(p1.row_degeneracy(cutoff=1.5), [4, 4, 4, 4, 4])

        self.assertEqual(p2.row_degeneracy(cutoff=.95), [4, 2, 4, 1])
        self.assertEqual(p2.row_degeneracy(cutoff=1.4), [4, 3, 4, 1])

        self.assertEqual(self.empty.row_degeneracy(), [])

    def test_column_degeneracy(self):
        """column_degeneracy: shoudl work as expected"""
        p1 = self.consensus
        p1.Data = transpose(p1.Data)
        p2 = self.not_same_value
        p2.Data = transpose(p2.Data)
        p1d = p1.column_degeneracy()
        self.assertEqual(p1d, [1, 1, 1, 2, 1])
        self.assertEqual(p1.column_degeneracy(cutoff=.5), [1, 1, 1, 2, 1])
        self.assertEqual(p1.column_degeneracy(cutoff=.75), [1, 2, 1, 3, 2])
        # when a row seems to add up to the cutoff value, it's not
        # always found because of floating point error. E.g. second row
        # in this example
        self.assertEqual(p1.column_degeneracy(cutoff=1), [2, 4, 1, 4, 2])
        # when the cutoff can't be found, the number of rows in the
        # profile is returned (for each column)
        self.assertEqual(p1.column_degeneracy(cutoff=1.5), [4, 4, 4, 4, 4])

        self.assertEqual(p2.column_degeneracy(cutoff=.95), [4, 2, 4, 1])
        self.assertEqual(p2.column_degeneracy(cutoff=1.4), [4, 3, 4, 1])

        self.assertEqual(self.empty.column_degeneracy(), [])

    def test_row_max(self):
        """row_max should return max value in each row"""
        p1 = self.consensus
        obs = p1.row_max()
        self.assertEqual(obs, array([.8, .7, 1, .4, .5]))

    def test_to_consensus(self):
        """to_consensus: should work with all the different options
        """
        p = self.consensus
        self.assertEqual(p.to_consensus(fully_degenerate=False), "AGGAT")
        self.assertEqual(p.to_consensus(fully_degenerate=True), "WVGNY")
        self.assertEqual(p.to_consensus(cutoff=0.75), "ARGHY")
        self.assertEqual(p.to_consensus(cutoff=0.95), "WVGNY")
        self.assertEqual(p.to_consensus(cutoff=2), "WVGNY")

        p = self.not_same_value
        self.assertEqual(p.to_consensus(fully_degenerate=False), "CGTA")
        self.assertEqual(p.to_consensus(fully_degenerate=True), "NBYA")
        self.assertEqual(p.to_consensus(cutoff=0.75), "YSYA")
        self.assertEqual(p.to_consensus(cutoff=2), "NBYA")
        self.assertEqual(p.to_consensus(cutoff=5), "NBYA")

        # when you specify both fully_generate and a cutoff value
        # the cutoff takes priority and is used in the calculation
        self.assertEqual(p.to_consensus(cutoff=0.75, fully_degenerate=True),
                         "YSYA")

        # raises AttributeError when Alphabet doens't have degenerates
        p = Profile(array([[.2, .8], [.7, .3]]), "AB")
        self.assertRaises(AttributeError, p.to_consensus, cutoff=.5)

    def test_to_consensus_include_all(self):
        """to_consensus: Should include all possibilities when include_all=True
        """
        p1 = Profile(array([[.2, 0, .8, 0], [0, .1, .2, .7], [0, 0, 0, 1],
                            [.2, .3, .4, .1], [.5, .5, 0, 0]]),
                     alphabet=DNA, CharOrder="TCAG")
        self.assertEqual(p1.to_consensus(cutoff=0.4, include_all=True),
                         "AGGAY")
        p2 = Profile(array([[.25, 0.25, .25, 0.25], [0.1, .1, .1, 0],
                            [.4, 0, .4, 0], [0, .2, 0.2, 0.3]]),
                     alphabet=DNA, CharOrder="TCAG")
        self.assertEqual(p2.to_consensus(cutoff=0.4,
                                        include_all=True), "NHWV")

    def test_random_indices(self):
        """random_indices: 99% of new frequencies should be within 3*SD
        """
        r_num, c_num = 100, 20
        num_elements = r_num * c_num
        r = random([r_num, c_num])
        p = Profile(r, "A" * c_num)
        p.normalize_positions()
        d = p.Data
        n = 1000

        # Test only works on normalized profile, b/c of 1-d below
        means = n * d
        three_stds = sqrt(d * (1 - d) * n) * 3
        result = [p.random_indices() for x in range(n)]
        a = Alignment(transpose(result))

        def absoluteProfile(alignment, char_order):
            f = a.column_freqs()
            res = zeros([len(f), len(char_order)])
            for row, freq in enumerate(f):
                for i in freq:
                    res[row, ord(i)] = freq[i]
            return res

        ap = absoluteProfile(a, p.CharOrder)
        failure = abs(ap - means) > three_stds
        assert sum(sum(failure)) / num_elements <= 0.01

    def test_random_sequence(self):
        """random_sequence: 99% of new frequencies should be within 3*SD"""
        r_num, c_num = 100, 20
        num_elements = r_num * c_num
        alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        r = random([r_num, c_num])
        p = Profile(r, alpha[:c_num])
        p.normalize_positions()
        d = p.Data
        n = 1000

        # Test only works on normalized profile, b/c of 1-d below
        means = n * d
        three_stds = sqrt(d * (1 - d) * n) * 3

        a = Alignment([p.random_sequence() for x in range(n)])

        def absoluteProfile(alignment, char_order):
            f = a.column_freqs()
            res = zeros([len(f), len(char_order)])
            for row, freq in enumerate(f):
                for i in freq:
                    col = char_order.index(i)
                    res[row, col] = freq[i]
            return res

        ap = absoluteProfile(a, p.CharOrder)
        failure = abs(ap - means) > three_stds
        assert sum(sum(failure)) / num_elements <= 0.01


class ModuleLevelFunctionsTest(TestCase):
    """Contains tests for the module level functions in profile.py"""

    def setUp(self):
        """setUp to change the alphabet for testing general CharMeaningProfile
        """
        self.alt_dna = DNA
        DnaDegenerateSymbols = {'R': 'AG',
            'N': 'TCAG', 'Y': 'TC', '?': 'TCAG-'}
        self.alt_dna.degenerates = DnaDegenerateSymbols

    def test_CharMeaningProfile(self):
        """CharMeaningProfile: should work as expected
        """
        p1 = CharMeaningProfile(self.alt_dna, "AGCT")
        p1_exp = [('A', [1, 0, 0, 0]), ('G', [0, 1, 0, 0]), ('C', [0, 0, 1, 0]),
                  ('T', [0, 0, 0, 1])]
        p2 = CharMeaningProfile(self.alt_dna, "TCAG")
        p2_exp = [('A', [0, 0, 1, 0]), ('G', [0, 0, 0, 1]), ('C', [0, 1, 0, 0]),
                  ('T', [1, 0, 0, 0])]
        # split_degen, but only whose chars are all in char order
        # so ? is ignored right now
        p3 = CharMeaningProfile(self.alt_dna, "TCAG", split_degenerates=True)
        p3_exp = [('A', [0, 0, 1, 0]), ('G', [0, 0, 0, 1]), ('C', [0, 1, 0, 0]),
                  ('T', [1, 0, 0, 0]), ('R', [
                   0, 0, .5, .5]), ('Y', [.5, .5, 0, 0]),
                  ('N', [.25, .25, .25, .25])]
        # if we add '-' to the character order, ? is split up as well
        p4 = CharMeaningProfile(self.alt_dna, "TCAG-", split_degenerates=True)
        p4_exp = [('A', [0, 0, 1, 0, 0]), ('G', [0, 0, 0, 1, 0]), ('C', [0, 1, 0, 0, 0]),
                  ('T', [1, 0, 0, 0, 0]), ('R', [
                   0, 0, .5, .5, 0]), ('Y', [.5, .5, 0, 0, 0]),
                  ('N', [.25, .25, .25, .25, 0]), ('-', [0, 0, 0, 0, 1]), ('?', [.2, .2, .2, .2, .2])]
        # Degenerate characters in the character order, when split_degenerates
        # is True, won't be split up, they'll get a 1 in their own column.
        p5 = CharMeaningProfile(self.alt_dna, "AGN", split_degenerates=True)
        p5_exp = [('A', [1, 0, 0]), ('G', [0, 1, 0]), ('N', [0, 0, 1]),
                  ('R', [.5, .5, 0])]
        # defaults char_order to list(alphabet)
        p6 = CharMeaningProfile(self.alt_dna)
        p6_exp = [('A', [0, 0, 1, 0]), ('G', [0, 0, 0, 1]), ('C', [0, 1, 0, 0]),
                  ('T', [1, 0, 0, 0])]
        # also accepts empty char_order -> set to list(alphabet)
        p7 = CharMeaningProfile(self.alt_dna, "")
        p7_exp = [('A', [0, 0, 1, 0]), ('G', [0, 0, 0, 1]), ('C', [0, 1, 0, 0]),
                  ('T', [1, 0, 0, 0])]

        for obs, exp in [(p1, p1_exp), (p2, p2_exp), (p3, p3_exp), (p4, p4_exp),
                        (p5, p5_exp), (p6, p6_exp), (p7, p7_exp)]:
            nz = [(chr(i), r.tolist())
                   for i, r in enumerate(obs.Data) if r.any()]
            self.assertEqualItems(nz, exp)

        self.assertRaises(ValueError, CharMeaningProfile, self.alt_dna,
                          "AGNX", split_degenerates=True)

if __name__ == "__main__":
    main()
