#!/usr/bin/env python
"""Provides tests for classes and functions in profile.py
"""
from collections import Counter
from numpy import array, sum, sqrt, transpose, add, subtract, multiply,\
    divide, zeros, vstack
from numpy.testing import assert_allclose
from numpy.random import random

from cogent3.util.unit_test import TestCase, main  # , numpy_err
from cogent3.core.moltype import DNA
from cogent3.core.sequence import ArraySequence
from cogent3.core.profile import PSSM, MotifCountsArray, MotifFreqsArray
from cogent3.core.alignment import ArrayAlignment as Alignment

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class MotifCountsArrayTests(TestCase):
    def test_construct_succeeds(self):
        """construct from int array or list"""
        from cogent3.maths.stats.number import CategoryCounter
        states = 'ACGT'
        rows = [CategoryCounter([b] * 20) for b in 'ACGT']
        rows = [r.tolist(states) for r in rows]
        pwm = MotifCountsArray(rows, states)

        data = [[2, 4], [3, 5], [4, 8]]
        got = MotifCountsArray(array(data), "AB")
        self.assertEqual(got.array.tolist(), data)

        got = MotifCountsArray(data, "AB")
        self.assertEqual(got.array.tolist(), data)

    def test_construct_fails(self):
        """fails if given wrong data type or no data"""
        # can't use a string
        data = [['A', 'A'], ['A', 'A'], ['A', 'A']]
        with self.assertRaises(ValueError):
            got = MotifCountsArray(data, "AB")

        # or a float
        data = [[1.1, 2.1], [0.0, 2.1], [3., 4.5]]
        with self.assertRaises(ValueError):
            got = MotifCountsArray(data, "AB")
        # or be empty
        with self.assertRaises(ValueError):
            got = MotifCountsArray([], "AB")

        with self.assertRaises(ValueError):
            got = MotifCountsArray([[], []], "AB")

        data = [[2, 4], [3, 5], [4, 8]]
        with self.assertRaises(ValueError):
            pssm = PSSM(data, "ACGT")

    def test_str_repr(self):
        """exercise str and repr"""
        data = array([[2, 4],
                      [3, 5],
                      [4, 8]])
        marr = MotifCountsArray(array(data), "AB")
        str(marr)
        repr(marr)

    def test_getitem(self):
        """slicing should return correct class"""
        data = array([[2, 4],
                      [3, 5],
                      [4, 8]])
        marr = MotifCountsArray(array(data), "AB")
        # print(marr[[1, 2], :])
        self.assertEqual(marr[0].array.tolist(), [2, 4])
        self.assertEqual(marr[0, 'B'], 4)
        self.assertEqual(marr[0, :].array.tolist(), [2, 4])
        self.assertEqual(marr[:, 'A'].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr[:, 'A':'B'].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr[1, 'A'], 3)
        marr = MotifCountsArray(array(data), "AB", row_indices=['a', 'b', 'c'])
        self.assertEqual(marr['a'].array.tolist(), [2, 4])
        self.assertEqual(marr['a', 'B'], 4)
        self.assertEqual(marr['a', :].array.tolist(), [2, 4])
        self.assertEqual(marr[:, 'A'].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr[:, 'A':'B'].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr['b', 'A'], 3)

    def test_todict(self):
        """correctly converts to a dict"""
        motifs = ['A', 'C', 'D']
        counts =  [[4, 0, 0]]
        marr = MotifCountsArray(counts, motifs)
        self.assertEqual(marr.todict(), {0: {'A': 4, 'C': 0, 'D': 0}})

    def test_to_freqs(self):
        """produces a freqs array"""
        data = array([[2, 4],
                      [3, 5],
                      [4, 8]])
        marr = MotifCountsArray(array(data), "AB")
        expect = data / vstack(data.sum(axis=1))
        got = marr.to_freq_array()
        assert_allclose(got.array, expect)

    def test_to_pssm(self):
        """produces a PSSM array"""
        data = array([[10, 30, 50, 10],
                      [25, 25, 25, 25],
                      [5, 80, 5, 10],
                      [70, 10, 10, 10],
                      [60, 15, 5, 20]])
        marr = MotifCountsArray(array(data), "ACGT")
        got = marr.to_pssm()
        expect = array(
            [[-1.322, 0.263, 1., -1.322],
             [0., 0., 0., 0.],
             [-2.322, 1.678, -2.322, -1.322],
             [1.485, -1.322, -1.322, -1.322],
             [1.263, -0.737, -2.322, -0.322]])
        assert_allclose(got.array, expect, atol=1e-3)

    def test_iter(self):
        """iter count array traverses positions"""
        data = [[2, 4], [3, 5], [4, 8]]
        got = MotifCountsArray(array(data), "AB")
        for row in got:
            self.assertEqual(row.shape, (2,))

    def test_take(self):
        """take works like numpy take, supporting negation"""
        data = array([[2, 4, 9, 2],
                      [3, 5, 8, 0],
                      [4, 8, 25, 13]])
        marr = MotifCountsArray(data, ['A', 'B', 'C', 'D'])
        # fails if don't provide an indexable indices
        with self.assertRaises(ValueError):
            marr.take(1, axis=1)

        # indexing columns using keys
        cols = marr.take(['A', 'D'], axis=1)
        assert_allclose(cols.array, data.take([0, 3], axis=1))
        cols = marr.take(['A', 'D'], negate=True, axis=1)
        assert_allclose(cols.array, data.take([1, 2], axis=1))
        # indexing columns using indexs
        cols = marr.take([0, 3], axis=1)
        assert_allclose(cols.array, data.take([0, 3], axis=1))
        cols = marr.take([0, 3], negate=True, axis=1)
        assert_allclose(cols.array, data.take([1, 2], axis=1))

        marr = MotifCountsArray(data, ['A', 'B', 'C', 'D'],
                                row_indices=['a', 'b', 'c'])
        # rows using keys
        rows = marr.take(['a', 'c'], axis=0)
        assert_allclose(rows.array, data.take([0, 2], axis=0))
        rows = marr.take(['a'], negate=True, axis=0)
        assert_allclose(rows.array, data.take([1, 2], axis=0))
        # rows using indexes
        rows = marr.take([0, 2], axis=0)
        assert_allclose(rows.array, data.take([0, 2], axis=0))
        rows = marr.take([0], negate=True, axis=0)
        assert_allclose(rows.array, data.take([1, 2], axis=0))

        # 1D profile
        marr = MotifCountsArray(data[0], ['A', 'B', 'C', 'D'])
        cols = marr.take([0], negate=True, axis=1)
        assert_allclose(cols.array, data[0].take([1, 2, 3]))


class MotifFreqsArrayTests(TestCase):
    def test_construct_succeeds(self):
        """construct from float array or list"""
        data = [[2 / 6, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        got = MotifFreqsArray(array(data), "AB")
        data = [[2 / 6, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        got = MotifFreqsArray(data, "AB")

    def test_construct_fails(self):
        """valid freqs only"""
        # no negatives
        data = [[-2 / 6, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        with self.assertRaises(ValueError):
            got = MotifFreqsArray(data, "AB")

        # must sum to 1 on axis=1
        data = [[2 / 5, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        with self.assertRaises(ValueError):
            got = MotifFreqsArray(data, "AB")

        data = [['A', 'A'], ['A', 'A'], ['A', 'A']]
        with self.assertRaises(ValueError):
            got = MotifFreqsArray(data, "AB")

        # int's not allowed
        data = [[2, 4], [3, 5], [4, 8]]
        with self.assertRaises(ValueError):
            got = MotifFreqsArray(data, "AB")

    def test_entropy(self):
        """calculates entripies correctly"""
        data = [[.25, .25, .25, .25],
                [.5, .5, 0, 0]]
        got = MotifFreqsArray(array(data), "ABCD")
        entropy = got.entropy()
        assert_allclose(entropy, [2, 1])

    def test_information(self):
        """calculates entr0pies correctly"""
        data = [[.25, .25, .25, .25],
                [.5, .5, 0, 0]]
        got = MotifFreqsArray(array(data), "ABCD")
        entropy = got.information()
        assert_allclose(entropy, [0, 1])

    def test_sim_seq(self):
        """exercising simulating a sequence"""
        data = [[.25, .25, .25, .25],
                [.5, .5, 0, 0]]
        farr = MotifFreqsArray(array(data), "ACGT")
        pos1 = Counter()
        pos2 = Counter()
        num = 1000
        for i in range(num):
            seq = farr.simulate_seq()
            self.assertEqual(len(seq), 2)
            pos1[seq[0]] += 1
            pos2[seq[1]] += 1
        self.assertEqual(len(pos1), 4)
        self.assertEqual(len(pos2), 2)
        self.assertTrue(min(pos1.values()) > 0)
        assert_allclose(pos2['C'] + pos2['A'], num)
        self.assertTrue(0 < pos2['C'] / num < 1)


class PSSMTests(TestCase):
    def test_construct_succeeds(self):
        """correctly construction from freqs"""
        data = [[.1, .3, .5, .1],
                [.25, .25, .25, .25],
                [.05, .8, .05, .1],
                [.7, .1, .1, .1],
                [.6, .15, .05, .2]]
        pssm = PSSM(data, "ACTG")
        expect = array(
            [[-1.322, 0.263, 1., -1.322],
             [0., 0., 0., 0.],
             [-2.322, 1.678, -2.322, -1.322],
             [1.485, -1.322, -1.322, -1.322],
             [1.263, -0.737, -2.322, -0.322]])
        assert_allclose(pssm.array, expect, atol=1e-3)

    def test_construct_fails(self):
        """construction fails for non-freq based input"""
        data = [[2, 4], [3, 5], [4, 8]]
        with self.assertRaises(ValueError):
            pssm = PSSM(data, "AB")

        # fails for negative numbers
        data = [[-1.322, 0.263, 1., -1.322],
                [0., 0., 0., 0.],
                [-2.322, 1.678, -2.322, -1.322],
                [1.485, -1.322, -1.322, -1.322],
                [1.263, -0.737, -2.322, -0.322]]
        with self.assertRaises(ValueError):
            pssm = PSSM(data, "ACTG")

    def test_score_indices(self):
        """produce correct score from indexed seq"""
        data = [[.1, .3, .5, .1],
                [.25, .25, .25, .25],
                [.05, .8, .05, .1],
                [.7, .1, .1, .1],
                [.6, .15, .05, .2]]
        pssm = PSSM(data, "ACTG")
        indices = [3, 1, 2, 0, 2, 2, 3]
        scores = pssm.score_indexed_seq(indices)
        assert_allclose(scores, [-4.481, -5.703, -2.966], atol=1e-3)

    def test_score_str(self):
        """produce correct score from seq"""
        data = [[.1, .3, .5, .1],
                [.25, .25, .25, .25],
                [.05, .8, .05, .1],
                [.7, .1, .1, .1],
                [.6, .15, .05, .2]]
        pssm = PSSM(data, "ACTG")
        seq = ''.join("ACTG"[i] for i in [3, 1, 2, 0, 2, 2, 3])
        scores = pssm.score_seq(seq)
        assert_allclose(scores, [-4.481, -5.703, -2.966], atol=1e-3)

    def test_score_seq_obj(self):
        """produce correct score from seq"""
        from cogent3 import DNA
        data = [[.1, .3, .5, .1],
                [.25, .25, .25, .25],
                [.05, .8, .05, .1],
                [.7, .1, .1, .1],
                [.6, .15, .05, .2]]
        pssm = PSSM(data, "ACTG")
        seq = DNA.make_seq(''.join("ACTG"[i] for i in [3, 1, 2, 0, 2, 2, 3]))
        scores = pssm.score_seq(seq)
        assert_allclose(scores, [-4.481, -5.703, -2.966], atol=1e-3)


if __name__ == "__main__":
    main()
