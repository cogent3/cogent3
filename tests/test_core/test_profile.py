#!/usr/bin/env python
"""Provides tests for classes and functions in profile.py
"""
from collections import Counter
from unittest import TestCase, main

from numpy import array, log2, nan, vstack
from numpy.testing import assert_allclose

from cogent3.core.profile import PSSM, MotifCountsArray, MotifFreqsArray


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Sandra Smit", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class MotifCountsArrayTests(TestCase):
    def test_construct_succeeds(self):
        """construct from int array or list"""
        from cogent3.maths.stats.number import CategoryCounter

        states = "ACGT"
        rows = [CategoryCounter([b] * 20) for b in "ACGT"]
        rows = [r.tolist(states) for r in rows]
        MotifCountsArray(rows, states)

        data = [[2, 4], [3, 5], [4, 8]]
        got = MotifCountsArray(array(data), "AB")
        self.assertEqual(got.array.tolist(), data)

        got = MotifCountsArray(data, "AB")
        self.assertEqual(got.array.tolist(), data)

    def test_construct_fails(self):
        """fails if given wrong data type or no data"""
        # can't use a string
        data = [["A", "A"], ["A", "A"], ["A", "A"]]
        with self.assertRaises(ValueError):
            MotifCountsArray(data, "AB")

        # or a float
        data = [[1.1, 2.1], [0.0, 2.1], [3.0, 4.5]]
        with self.assertRaises(ValueError):
            MotifCountsArray(data, "AB")
        # or be empty
        with self.assertRaises(ValueError):
            MotifCountsArray([], "AB")

        with self.assertRaises(ValueError):
            MotifCountsArray([[], []], "AB")

        data = [[2, 4], [3, 5], [4, 8]]
        with self.assertRaises(ValueError):
            PSSM(data, "ACGT")

    def test_str_repr(self):
        """exercise str and repr"""
        data = array([[2, 4], [3, 5], [4, 8]])
        marr = MotifCountsArray(array(data), "AB")
        str(marr)
        repr(marr)

    def test_getitem(self):
        """slicing should return correct class"""
        data = array([[2, 4], [3, 5], [4, 8]])
        marr = MotifCountsArray(array(data), "AB")
        # print(marr[[1, 2], :])
        self.assertEqual(marr[0].array.tolist(), [2, 4])
        self.assertEqual(marr[0, "B"], 4)
        self.assertEqual(marr[0, :].array.tolist(), [2, 4])
        self.assertEqual(marr[:, "A"].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr[:, "A":"B"].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr[1, "A"], 3)
        marr = MotifCountsArray(array(data), "AB", row_indices=["a", "b", "c"])
        self.assertEqual(marr["a"].array.tolist(), [2, 4])
        self.assertEqual(marr["a", "B"], 4)
        self.assertEqual(marr["a", :].array.tolist(), [2, 4])
        self.assertEqual(marr[:, "A"].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr[:, "A":"B"].array.tolist(), [[2], [3], [4]])
        self.assertEqual(marr["b", "A"], 3)

    def test_sliced_range(self):
        """a sliced range should preserve row indices"""
        motifs = ("A", "C", "G", "T")
        names = ["FlyingFox", "DogFaced", "FreeTaile"]
        data = [[316, 134, 133, 317], [321, 136, 123, 314], [331, 143, 127, 315]]
        counts = MotifCountsArray(data, motifs, row_indices=names)
        self.assertEqual(counts.keys(), names)
        subset = counts[:2]
        self.assertEqual(subset.keys(), names[:2])

    def test_to_dict(self):
        """correctly converts to a dict"""
        motifs = ["A", "C", "D"]
        counts = [[4, 0, 0]]
        marr = MotifCountsArray(counts, motifs)
        self.assertEqual(marr.to_dict(), {0: {"A": 4, "C": 0, "D": 0}})

    def test_to_freqs(self):
        """produces a freqs array"""
        data = array([[2, 4], [3, 5], [4, 8]])
        marr = MotifCountsArray(array(data), "AB")
        expect = data / vstack(data.sum(axis=1))
        got = marr.to_freq_array()
        assert_allclose(got.array, expect)

    def test_to_freqs_pseudocount(self):
        """produces a freqs array with pseudocount"""
        data = array([[2, 4], [3, 5], [0, 8]])
        marr = MotifCountsArray(array(data), "AB")
        got = marr.to_freq_array(pseudocount=1)
        adj = data + 1
        expect = adj / vstack(adj.sum(axis=1))
        assert_allclose(got.array, expect)

        got = marr.to_freq_array(pseudocount=0.5)
        adj = data + 0.5
        expect = adj / vstack(adj.sum(axis=1))
        assert_allclose(got.array, expect)

    def test_to_freqs_1d(self):
        """produce a freqs array from 1D counts"""
        data = [43, 48, 114, 95]
        total = sum(data)
        a = MotifCountsArray([43, 48, 114, 95], motifs=("T", "C", "A", "G"))
        f = a.to_freq_array()
        assert_allclose(f.array, array([v / total for v in data], dtype=float))

    def test_to_pssm(self):
        """produces a PSSM array"""
        data = array(
            [
                [10, 30, 50, 10],
                [25, 25, 25, 25],
                [5, 80, 5, 10],
                [70, 10, 10, 10],
                [60, 15, 5, 20],
            ]
        )
        marr = MotifCountsArray(array(data), "ACGT")
        got = marr.to_pssm()
        expect = array(
            [
                [-1.322, 0.263, 1.0, -1.322],
                [0.0, 0.0, 0.0, 0.0],
                [-2.322, 1.678, -2.322, -1.322],
                [1.485, -1.322, -1.322, -1.322],
                [1.263, -0.737, -2.322, -0.322],
            ]
        )
        assert_allclose(got.array, expect, atol=1e-3)

    def test_to_pssm_pseudocount(self):
        """produces a PSSM array with pseudocount"""
        data = array(
            [
                [10, 30, 50, 10],
                [25, 25, 25, 25],
                [5, 80, 0, 10],
                [70, 10, 10, 10],
                [60, 15, 0, 20],
            ]
        )
        marr = MotifCountsArray(array(data), "ACGT")
        got = marr.to_pssm(pseudocount=1)
        freqs = marr._to_freqs(pseudocount=1)
        expect = log2(freqs / 0.25)
        assert_allclose(got.array, expect, atol=1e-3)

    def test_iter(self):
        """iter count array traverses positions"""
        data = [[2, 4], [3, 5], [4, 8]]
        got = MotifCountsArray(array(data), "AB")
        for row in got:
            self.assertEqual(row.shape, (2,))

    def test_take(self):
        """take works like numpy take, supporting negation"""
        data = array([[2, 4, 9, 2], [3, 5, 8, 0], [4, 8, 25, 13]])
        marr = MotifCountsArray(data, ["A", "B", "C", "D"])
        # fails if don't provide an indexable indices
        with self.assertRaises(ValueError):
            marr.take(1, axis=1)

        # indexing columns using keys
        cols = marr.take(["A", "D"], axis=1)
        assert_allclose(cols.array, data.take([0, 3], axis=1))
        cols = marr.take(["A", "D"], negate=True, axis=1)
        assert_allclose(cols.array, data.take([1, 2], axis=1))
        # indexing columns using indexs
        cols = marr.take([0, 3], axis=1)
        assert_allclose(cols.array, data.take([0, 3], axis=1))
        cols = marr.take([0, 3], negate=True, axis=1)
        assert_allclose(cols.array, data.take([1, 2], axis=1))

        marr = MotifCountsArray(data, ["A", "B", "C", "D"], row_indices=["a", "b", "c"])
        # rows using keys
        rows = marr.take(["a", "c"], axis=0)
        assert_allclose(rows.array, data.take([0, 2], axis=0))
        rows = marr.take(["a"], negate=True, axis=0)
        assert_allclose(rows.array, data.take([1, 2], axis=0))
        # rows using indexes
        rows = marr.take([0, 2], axis=0)
        assert_allclose(rows.array, data.take([0, 2], axis=0))
        rows = marr.take([0], negate=True, axis=0)
        assert_allclose(rows.array, data.take([1, 2], axis=0))

        # 1D profile
        marr = MotifCountsArray(data[0], ["A", "B", "C", "D"])
        cols = marr.take([0], negate=True, axis=1)
        assert_allclose(cols.array, data[0].take([1, 2, 3]))


class MotifFreqsArrayTests(TestCase):
    def test_construct_succeeds(self):
        """construct from float array or list"""
        data = [[2 / 6, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        MotifFreqsArray(array(data), "AB")
        data = [[2 / 6, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        MotifFreqsArray(data, "AB")

    def test_construct_fails(self):
        """valid freqs only"""
        # no negatives
        data = [[-2 / 6, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        with self.assertRaises(ValueError):
            MotifFreqsArray(data, "AB")

        # must sum to 1 on axis=1
        data = [[2 / 5, 4 / 6], [3 / 8, 5 / 8], [4 / 12, 8 / 12]]
        with self.assertRaises(ValueError):
            MotifFreqsArray(data, "AB")

        data = [["A", "A"], ["A", "A"], ["A", "A"]]
        with self.assertRaises(ValueError):
            MotifFreqsArray(data, "AB")

        # int's not allowed
        data = [[2, 4], [3, 5], [4, 8]]
        with self.assertRaises(ValueError):
            MotifFreqsArray(data, "AB")

    def test_entropy_terms(self):
        """Checks entropy_terms works correctly"""
        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        got = MotifFreqsArray(array(data), "ABCD")
        entropy_terms = got.entropy_terms()
        expect = [[0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0, 0]]
        assert_allclose(entropy_terms.array, expect)

    def test_entropy(self):
        """calculates entripies correctly"""
        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        got = MotifFreqsArray(array(data), "ABCD")
        entropy = got.entropy()
        assert_allclose(entropy, [2, 1])

    def test_relative_entropy_terms(self):
        """Check that relative_entropy_terms works for different background distributions"""
        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        got = MotifFreqsArray(array(data), "ABCD")
        rel_entropy = got.relative_entropy_terms(background=None)
        expected = [[0, 0, 0, 0], [-0.25, -0.25, -0.5, -0.5]]
        assert_allclose(rel_entropy, expected)

        background = {"A": 0.5, "B": 0.25, "C": 0.125, "D": 0.125}
        rel_entropy = got.relative_entropy_terms(background=background)
        expected = [[0.5, 0, -0.125, -0.125], [0, -0.25, -0.375, -0.375]]
        assert_allclose(rel_entropy, expected)

        with self.assertRaises(ValueError):
            got.relative_entropy_terms(background=dict(A=-0.5, B=1.5))

        with self.assertRaises(ValueError):
            got.relative_entropy_terms(background={"A": 0.5, "B": 0.25, "C": 0.125})

    def test_relative_entropy(self):
        """calculates relative entropy correctly"""
        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        got = MotifFreqsArray(array(data), "ABCD")
        rel_entropy = got.relative_entropy(background=None)
        assert_allclose(rel_entropy, [0, -1.5])

        background = {"A": 0.5, "B": 0.25, "C": 0.125, "D": 0.125}
        rel_entropy = got.relative_entropy(background=background)
        expected = [0.25, -1]
        assert_allclose(rel_entropy, expected)

    def test_pairwise_jsd(self):
        """correctly constructs pairwise JSD dict"""
        from numpy.random import random

        from cogent3.maths.measure import jsd

        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        expect = jsd(data[0], data[1])
        freqs = MotifFreqsArray(array(data), "ACGT")
        got = freqs.pairwise_jsd()
        assert_allclose(list(got.values())[0], expect)

        data = []
        for _ in range(6):
            freqs = random(4)
            freqs = freqs / freqs.sum()
            data.append(freqs)

        freqs = MotifFreqsArray(array(data), "ACGT")
        pwise = freqs.pairwise_jsd()
        self.assertEqual(len(pwise), 6 * 6 - 6)

    def test_pairwise_jsm(self):
        """correctly constructs pairwise JS metric dict"""
        from numpy.random import random

        from cogent3.maths.measure import jsm

        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        expect = jsm(data[0], data[1])
        freqs = MotifFreqsArray(array(data), "ACGT")
        got = freqs.pairwise_jsm()
        assert_allclose(list(got.values())[0], expect)

        data = []
        for _ in range(6):
            freqs = random(4)
            freqs = freqs / freqs.sum()
            data.append(freqs)

        freqs = MotifFreqsArray(array(data), "ACGT")
        pwise = freqs.pairwise_jsm()
        self.assertEqual(len(pwise), 6 * 6 - 6)

    def test_pairwise_(self):
        """returns None when single row"""
        # ndim=1
        data = [0.25, 0.25, 0.25, 0.25]
        freqs = MotifFreqsArray(array(data), "ACGT")
        got = freqs.pairwise_jsm()
        self.assertEqual(got, None)

        # ndim=2
        data = array([[0.25, 0.25, 0.25, 0.25]])
        freqs = MotifFreqsArray(data, "ACGT")
        got = freqs.pairwise_jsm()
        self.assertEqual(got, None)

    def test_information(self):
        """calculates entropies correctly"""
        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        got = MotifFreqsArray(array(data), "ABCD")
        entropy = got.information()
        assert_allclose(entropy, [0, 1])

    def test_sim_seq(self):
        """exercising simulating a sequence"""
        data = [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0, 0]]
        farr = MotifFreqsArray(array(data), "ACGT")
        pos1 = Counter()
        pos2 = Counter()
        num = 1000
        for _ in range(num):
            seq = farr.simulate_seq()
            self.assertEqual(len(seq), 2)
            pos1[seq[0]] += 1
            pos2[seq[1]] += 1
        self.assertEqual(len(pos1), 4)
        self.assertEqual(len(pos2), 2)
        self.assertTrue(min(pos1.values()) > 0)
        assert_allclose(pos2["C"] + pos2["A"], num)
        self.assertTrue(0 < pos2["C"] / num < 1)

    def test_slicing(self):
        """slice by keys should work"""
        counts = MotifCountsArray(
            [[3, 2, 3, 2], [3, 2, 3, 2]],
            ["A", "C", "G", "T"],
            row_indices=["DogFaced", "FlyingFox"],
        )
        freqs = counts.to_freq_array()
        got = freqs["FlyingFox"].to_array()
        assert_allclose(got, [0.3, 0.2, 0.3, 0.2])

    def test_to_pssm(self):
        """construct PSSM from freqs array"""
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        farr = MotifFreqsArray(data, "ACTG")
        pssm = farr.to_pssm()
        expect = array(
            [
                [-1.322, 0.263, 1.0, -1.322],
                [0.0, 0.0, 0.0, 0.0],
                [-2.322, 1.678, -2.322, -1.322],
                [1.485, -1.322, -1.322, -1.322],
                [1.263, -0.737, -2.322, -0.322],
            ]
        )
        assert_allclose(pssm.array, expect, atol=1e-3)

    def test_logo(self):
        """produces a Drawable with correct layout elements"""
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        farr = MotifFreqsArray(data, "ACTG")
        # with defaults, has a single x/y axes and number of shapes
        logo = farr.logo(ylim=0.5)
        fig = logo.figure
        self.assertEqual(fig.data, [{}])
        self.assertTrue(len(fig.layout.xaxis) > 10)
        self.assertTrue(len(fig.layout.yaxis) > 10)
        self.assertEqual(fig.layout.yaxis.range, [0, 0.5])
        # since the second row are equi-frequent, their information is 0 so
        # we substract 4 shapes from that column
        self.assertEqual(len(fig.layout.shapes), farr.shape[0] * farr.shape[1] - 4)
        # wrapping across multiple rows should produce multiple axes
        logo = farr.logo(ylim=0.5, wrap=3)
        fig = logo.figure
        for axis in ("axis", "axis2"):
            self.assertIn(f"x{axis}", fig.layout)
            self.assertIn(f"y{axis}", fig.layout)

        # fails if vspace not in range 0-1
        with self.assertRaises(AssertionError):
            farr.logo(vspace=20)


class PSSMTests(TestCase):
    def test_construct_succeeds(self):
        """correctly construction from freqs"""
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        pssm = PSSM(data, "ACTG")
        expect = array(
            [
                [-1.322, 0.263, 1.0, -1.322],
                [0.0, 0.0, 0.0, 0.0],
                [-2.322, 1.678, -2.322, -1.322],
                [1.485, -1.322, -1.322, -1.322],
                [1.263, -0.737, -2.322, -0.322],
            ]
        )
        assert_allclose(pssm.array, expect, atol=1e-3)

    def test_construct_fails(self):
        """construction fails for invalid input"""

        # fails for entries all zero
        data_all_zero = [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
        with self.assertRaises(ValueError):
            PSSM(data_all_zero, "ACTG")

        # fails for numpy.nan
        data_nan = [
            [nan, -0.263, -1.0, -1.322],
            [-2.322, -1.678, -2.322, -1.322],
            [-1.485, -1.322, -1.322, -1.322],
            [-1.263, -0.737, -2.322, -0.322],
        ]
        with self.assertRaises(ValueError):
            PSSM(data_nan, "ACTG")

        # fails for entries all negative numbers
        data = [
            [-1.322, -0.263, -1.0, -1.322],
            [-2.322, -1.678, -2.322, -1.322],
            [-1.485, -1.322, -1.322, -1.322],
            [-1.263, -0.737, -2.322, -0.322],
        ]
        with self.assertRaises(ValueError):
            PSSM(data, "ACTG")

    def test_score_indices(self):
        """produce correct score from indexed seq"""
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        pssm = PSSM(data, "ACTG")
        indices = [3, 1, 2, 0, 2, 2, 3]
        scores = pssm.score_indexed_seq(indices)
        assert_allclose(scores, [-4.481, -5.703, -2.966], atol=1e-3)

        indices = [4, 1, 2, 0, 2, 2, 3]
        scores = pssm.score_indexed_seq(indices)
        # log2 of (0.25 * 0.05 * 0.7 * 0.05) / .25**4 = -3.158...
        assert_allclose(scores, [-3.158, -5.703, -2.966], atol=1e-3)

        # fails if sequence too short
        with self.assertRaises(ValueError):
            pssm.score_indexed_seq(indices[:3])

    def test_score_str(self):
        """produce correct score from seq"""
        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        pssm = PSSM(data, "ACTG")
        seq = "".join("ACTG"[i] for i in [3, 1, 2, 0, 2, 2, 3])
        scores = pssm.score_seq(seq)
        assert_allclose(scores, [-4.481, -5.703, -2.966], atol=1e-3)
        with self.assertRaises(ValueError):
            pssm.score_seq(seq[:3])

    def test_score_seq_obj(self):
        """produce correct score from seq"""
        from cogent3 import DNA

        data = [
            [0.1, 0.3, 0.5, 0.1],
            [0.25, 0.25, 0.25, 0.25],
            [0.05, 0.8, 0.05, 0.1],
            [0.7, 0.1, 0.1, 0.1],
            [0.6, 0.15, 0.05, 0.2],
        ]
        pssm = PSSM(data, "ACTG")
        seq = DNA.make_seq("".join("ACTG"[i] for i in [3, 1, 2, 0, 2, 2, 3]))
        scores = pssm.score_seq(seq)
        assert_allclose(scores, [-4.481, -5.703, -2.966], atol=1e-3)


if __name__ == "__main__":
    main()
