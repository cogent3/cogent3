#!/usr/bin/env python
"""Unit tests for distance_transform.py functions.
"""

from unittest import TestCase, main

from numpy import array, ones, shape, sqrt

from cogent3.maths.distance_transform import (
    binary_dist_chisq,
    binary_dist_chord,
    binary_dist_euclidean,
    binary_dist_hamming,
    binary_dist_jaccard,
    binary_dist_lennon,
    binary_dist_ochiai,
    binary_dist_otu_gain,
    binary_dist_pearson,
    binary_dist_sorensen_dice,
    dist_abund_jaccard,
    dist_bray_curtis,
    dist_bray_curtis_faith,
    dist_bray_curtis_magurran,
    dist_canberra,
    dist_chisq,
    dist_chord,
    dist_euclidean,
    dist_gower,
    dist_hellinger,
    dist_kulczynski,
    dist_manhattan,
    dist_morisita_horn,
    dist_pearson,
    dist_soergel,
    dist_spearman_approx,
    dist_specprof,
    numpy,
    trans_chisq,
    trans_chord,
    trans_hellinger,
    trans_specprof,
    zeros,
)


__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__contributors__ = ["Justin Kuczynski", "Zongzhi Liu", "Greg Caporaso"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

from numpy.testing import assert_allclose, assert_equal


class functionTests(TestCase):
    """Tests of top-level functions."""

    def setUp(self):
        self.mat_test = array([[10, 10, 20], [10, 15, 10], [15, 5, 5]], "float")

        self.emptyarray = array([], "d")
        self.mtx1 = array([[1, 3], [0.0, 23.1]], "d")
        self.dense1 = array([[1, 3], [5, 2], [0.1, 22]], "d")

        self.zeromtx = array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], "d"
        )
        self.sparse1 = array(
            [[0.0, 0.0, 5.33], [0.0, 0.0, 0.4], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], "d"
        )
        self.input_binary_dist_otu_gain1 = array(
            [[2, 1, 0, 0], [1, 0, 0, 1], [0, 0, 3, 0], [0, 0, 0, 1]]
        )

    def get_sym_mtx_from_uptri(self, mtx):
        """helper fn, only for square matrices"""
        numrows, numcols = shape(mtx)
        for i in range(numrows):
            for j in range(i):
                if i == j:
                    break
                mtx[i, j] = mtx[j, i]  # j < i, so row<col => upper triangle
        return mtx

    def test_dist_canberra(self):
        """tests dist_canberra

        tests inputs of empty mtx, zeros, and results compared with calcs done
        by hand"""

        assert_allclose(dist_canberra(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array([[0.0, 46.2 / 52.2], [46.2 / 52.2, 0.0]], "d")
        assert_allclose(dist_canberra(self.mtx1), mtx1expected)

        sparse1exp = ones((self.sparse1.shape[0], self.sparse1.shape[0]))
        # remove diagonal
        sparse1exp[0, 0] = sparse1exp[1, 1] = sparse1exp[2, 2] = sparse1exp[3, 3] = 0.0

        sparse1exp[0, 1] = sparse1exp[1, 0] = (5.33 - 0.4) / (5.33 + 0.4)
        assert_allclose(dist_canberra(self.sparse1), sparse1exp)

    def test_dist_canberra_bug(self):
        i = array([[0, 0, 1], [0, 1, 1]])
        d = (1.0 / 2.0) * sum(
            [abs(0.0 - 1.0) / (0.0 + 1.0), abs(1.0 - 1.0) / (1.0 + 1.0)]
        )
        expected = array([[0.0, d], [d, 0.0]])
        actual = dist_canberra(i)
        assert_allclose(expected, actual)

    def test_dist_euclidean(self):
        """tests dist_euclidean

        tests inputs of empty mtx, zeros, and dense1 compared with calcs done
        by hand"""

        assert_allclose(dist_euclidean(self.zeromtx), zeros((4, 4), "d"))

        dense1expected = array(
            [
                [0.0, sqrt(17.0), sqrt(0.9 ** 2 + 19 ** 2)],
                [sqrt(17.0), 0.0, sqrt(4.9 ** 2 + 20 ** 2)],
                [sqrt(0.9 ** 2 + 19 ** 2), sqrt(4.9 ** 2 + 20 ** 2), 0.0],
            ],
            "d",
        )
        assert_allclose(dist_euclidean(self.dense1), dense1expected)

    def test_dist_gower(self):
        """tests dist_gower

        tests inputs of empty mtx, zeros, and results compared with calcs done
        by hand"""

        assert_allclose(dist_gower(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array([[0.0, 2.0], [2.0, 0.0]], "d")
        assert_allclose(dist_gower(self.mtx1), mtx1expected)

        sparse1expected = array(
            [
                [0.0, 4.93 / 5.33, 2, 1],
                [4.93 / 5.33, 0.0, 1 + 0.4 / 5.33, 0.4 / 5.33],
                [2, 1 + 0.4 / 5.33, 0, 1],
                [1, 0.4 / 5.33, 1, 0.0],
            ],
            "d",
        )
        assert_allclose(dist_gower(self.sparse1), sparse1expected)

    def test_dist_manhattan(self):
        """tests dist_manhattan

        tests inputs of empty mtx, zeros, and dense1 compared with calcs done
        by hand"""

        assert_allclose(dist_manhattan(self.zeromtx), zeros((4, 4), "d"))

        dense1expected = array(
            [[0.0, 5.0, 019.9], [5.0, 0.0, 24.9], [19.9, 24.90, 0.0]], "d"
        )
        assert_allclose(dist_manhattan(self.dense1), dense1expected)

    def test_dist_abund_jaccard(self):
        """dist_abund_jaccard should compute distances for dense1 and mtx1"""
        mtx1_expected = array([[0, 0.25], [0.25, 0]], "d")
        assert_equal(dist_abund_jaccard(self.mtx1), mtx1_expected)

        dense1_expected = zeros((3, 3), "d")
        assert_equal(dist_abund_jaccard(self.dense1), dense1_expected)

        sparse1_expected = array(
            [
                [0.0, 0.0, 1.0, 1.0],
                [0.0, 0.0, 1.0, 1.0],
                [1.0, 1.0, 0.0, 1.0],
                [1.0, 1.0, 1.0, 0.0],
            ],
            "d",
        )
        assert_equal(dist_abund_jaccard(self.sparse1), sparse1_expected)

    def test_dist_morisita_horn(self):
        """tests dist_morisita_horn

        tests inputs of empty mtx, zeros, and dense1 compared with calcs done
        by hand"""

        assert_allclose(dist_morisita_horn(self.zeromtx), zeros((4, 4), "d"))

        a = 1 - 2 * 69.3 / (26 / 16.0 * 23.1 * 4)
        mtx1expected = array([[0, a], [a, 0]], "d")
        assert_allclose(dist_morisita_horn(self.mtx1), mtx1expected)

    def test_dist_bray_curtis(self):
        """tests dist_bray_curtis

        tests inputs of empty mtx, zeros, and mtx1 compared with calcs done
        by hand"""

        assert_allclose(dist_manhattan(self.zeromtx), zeros((4, 4) * 1, "d"))

        mtx1expected = array([[0, 21.1 / 27.1], [21.1 / 27.1, 0]], "d")
        assert_allclose(dist_bray_curtis(self.mtx1), mtx1expected)

    def test_dist_bray_curtis_faith(self):
        """tests dist_bray_curtis_faith

        tests inputs of empty mtx, zeros, and mtx1 compared with calcs done
        by hand"""

        assert_allclose(dist_manhattan(self.zeromtx), zeros((4, 4) * 1, "d"))

        mtx1expected = array([[0, 21.1 / 27.1], [21.1 / 27.1, 0]], "d")
        assert_allclose(dist_bray_curtis_faith(self.mtx1), mtx1expected)

    def test_dist_soergel(self):
        """tests dist_soergel

        tests inputs of empty mtx, zeros, and dense1 compared with calcs done
        by hand/manhattan dist"""

        assert_allclose(dist_soergel(self.zeromtx), zeros((4, 4) * 1, "d"))

        dense1expected = dist_manhattan(self.dense1)
        dense1norm = array([[1, 8, 23], [8, 1, 27], [23, 27, 1]], "d")
        dense1expected /= dense1norm

        assert_allclose(dist_soergel(self.dense1), dense1expected)

    def test_dist_kulczynski(self):
        """tests dist_kulczynski

        tests inputs of empty mtx, zeros, and mtx1 compared with calcs done
        by hand"""

        assert_allclose(dist_kulczynski(self.zeromtx), zeros((4, 4) * 1, "d"))

        mtx1expected = array(
            [
                [0, 1.0 - 1.0 / 2.0 * (3.0 / 4.0 + 3.0 / 23.1)],
                [1.0 - 1.0 / 2.0 * (3.0 / 4.0 + 3.0 / 23.1), 0],
            ],
            "d",
        )

        assert_allclose(dist_kulczynski(self.mtx1), mtx1expected)

    def test_dist_pearson(self):
        """tests dist_pearson

        tests inputs of empty mtx, zeros, mtx compared with calcs done
        by hand, and an example from
        http://davidmlane.com/hyperstat/A56626.html
        """

        assert_allclose(dist_pearson(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array([[0, 0], [0, 0]], "d")
        assert_allclose(dist_pearson(self.mtx1), mtx1expected)

        # example 1 from http://davidmlane.com/hyperstat/A56626.html
        ex1 = array([[1, 2, 3], [2, 5, 6]], "d")
        ex1res = 1 - 4.0 / sqrt(2.0 * (8 + 2.0 / 3.0))
        ex1expected = array([[0, ex1res], [ex1res, 0]], "d")

        assert_allclose(dist_pearson(ex1), ex1expected)

    def test_dist_spearman_approx(self):
        """tests dist_spearman_approx

        tests inputs of empty mtx, zeros, and an example from wikipedia
        """

        assert_allclose(dist_spearman_approx(self.zeromtx), zeros((4, 4) * 1, "d"))

        # ex1 from wikipedia Spearman's_rank_correlation_coefficient 20jan2009
        ex1 = array(
            [
                [106, 86, 100, 101, 99, 103, 97, 113, 112, 110],
                [7, 0, 27, 50, 28, 29, 20, 12, 6, 17],
            ],
            "d",
        )
        ex1res = 6.0 * 194.0 / (10.0 * 99.0)
        ex1expected = array([[0, ex1res], [ex1res, 0]], "d")
        assert_allclose(dist_spearman_approx(ex1), ex1expected)

    # now binary fns
    def test_binary_dist_otu_gain(self):
        """binary OTU gain functions as expected"""
        actual = binary_dist_otu_gain(self.input_binary_dist_otu_gain1)
        expected = array([[0, 1, 2, 2], [1, 0, 2, 1], [1, 1, 0, 1], [1, 0, 1, 0]])
        assert_equal(actual, expected)

    def test_binary_dist_chisq(self):
        """tests binary_dist_chisq

        tests inputs of empty mtx, zeros, and mtx1 compared with calcs done
        by hand"""

        assert_allclose(binary_dist_chisq(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array([[0, sqrt(9 / 8.0)], [sqrt(9 / 8.0), 0]], "d")
        assert_allclose(binary_dist_chisq(self.mtx1), mtx1expected)

    def test_binary_dist_chord(self):
        """tests binary_dist_chord

        tests inputs of empty mtx, zeros, and results compared with calcs done
        by hand"""

        assert_allclose(binary_dist_chord(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array(
            [
                [0, sqrt(1 / 2.0 + (1.0 / sqrt(2.0) - 1.0) ** 2)],
                [sqrt(1 / 2.0 + (1.0 / sqrt(2.0) - 1.0) ** 2), 0],
            ],
            "d",
        )
        assert_allclose(binary_dist_chord(self.mtx1), mtx1expected)

    def test_binary_dist_lennon(self):
        """tests binary_dist_lennon

        tests inputs of empty mtx, zeros, and results compared with calcs done
        by hand"""

        assert_allclose(binary_dist_lennon(self.zeromtx), zeros((4, 4), "d"))

        mtxa = array([[5.2, 9, 0.2], [0, 99, 1], [0, 0.0, 8233.1]], "d")
        assert_allclose(binary_dist_lennon(mtxa), zeros((3, 3), "d"))

        mtxb = array([[5.2, 0, 0.2, 9.2], [0, 0, 0, 1], [0, 3.2, 0, 8233.1]], "d")
        mtxbexpected = array([[0, 0, 0.5], [0, 0, 0], [0.5, 0, 0]], "d")
        assert_allclose(binary_dist_lennon(mtxb), mtxbexpected)

    def test_binary_dist_pearson(self):
        """tests binary_dist_pearson

        tests inputs of empty mtx, zeros, and dense1 compared with calcs done
        by hand"""

        assert_allclose(binary_dist_pearson(self.zeromtx), zeros((4, 4), "d"))

        assert_allclose(binary_dist_pearson(self.dense1), zeros((3, 3)))

    def test_binary_dist_jaccard(self):
        """tests binary_dist_jaccard

        tests inputs of empty mtx, zeros, and sparse1 compared with calcs done
        by hand"""

        assert_allclose(binary_dist_jaccard(self.zeromtx), zeros((4, 4), "d"))

        sparse1expected = array(
            [[0, 0, 1.0, 1.0], [0, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]], "d"
        )
        assert_allclose(binary_dist_jaccard(self.sparse1), sparse1expected)

        sparse1expected = dist_manhattan(self.sparse1.astype(bool))
        sparse1norm = array(
            [[1, 1, 2, 1], [1, 1, 2, 1], [2, 2, 1, 1], [1, 1, 1, 100]], "d"
        )
        sparse1expected /= sparse1norm
        assert_allclose(binary_dist_jaccard(self.sparse1), sparse1expected)

    def test_binary_dist_ochiai(self):
        """tests binary_dist_ochiai

        tests inputs of empty mtx, zeros, and mtx1 compared with calcs done
        by hand"""

        assert_allclose(binary_dist_ochiai(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array([[0, 1 - 1 / sqrt(2.0)], [1 - 1 / sqrt(2.0), 0]], "d")
        assert_allclose(binary_dist_ochiai(self.mtx1), mtx1expected)

    def test_binary_dist_hamming(self):
        """tests binary_dist_hamming

        tests inputs of empty mtx, zeros, and mtx1 compared with calcs done
        by hand"""

        assert_allclose(binary_dist_hamming(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array([[0, 1], [1, 0]], "d")
        assert_allclose(binary_dist_hamming(self.mtx1), mtx1expected)

    def test_binary_dist_sorensen_dice(self):
        """tests binary_dist_sorensen_dice

        tests inputs of empty mtx, zeros, and mtx1 compared with calcs done
        by hand"""

        assert_allclose(binary_dist_sorensen_dice(self.zeromtx), zeros((4, 4), "d"))

        mtx1expected = array([[0, 1 / 3.0], [1 / 3.0, 0]], "d")
        assert_allclose(binary_dist_sorensen_dice(self.mtx1), mtx1expected)

        sparse1expected = array(
            [[0, 0, 1.0, 1.0], [0, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]], "d"
        )

        assert_allclose(binary_dist_sorensen_dice(self.sparse1), sparse1expected)

    def test_binary_dist_euclidean(self):
        """tests binary_dist_euclidean

        tests two inputs compared with calculations by hand, and runs zeros
        and an empty input"""
        dense1expected = array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], "d")
        sparse1expected = zeros((4, 4), "d")
        sparse1expected[0, 2] = sqrt(2)
        sparse1expected[0, 3] = 1.0
        sparse1expected[1, 2] = sqrt(2)
        sparse1expected[1, 3] = 1.0
        sparse1expected[2, 3] = 1.0
        sparse1expected = self.get_sym_mtx_from_uptri(sparse1expected)

        assert_allclose(binary_dist_euclidean(self.dense1), dense1expected)
        assert_allclose(binary_dist_euclidean(self.sparse1), sparse1expected)
        assert_allclose(binary_dist_euclidean(self.zeromtx), zeros((4, 4), "d"))

    # zj's stuff
    def test_chord_transform(self):
        """trans_chord should return the exp result in the ref paper."""

        exp = [
            [0.40824829, 0.40824829, 0.81649658],
            [0.48507125, 0.72760688, 0.48507125],
            [0.90453403, 0.30151134, 0.30151134],
        ]
        res = trans_chord(self.mat_test)
        assert_allclose(res, exp)

    def test_chord_dist(self):
        """dist_chord should return the exp result."""

        assert_allclose(dist_chord(self.zeromtx), zeros((4, 4), "d"))

        exp = [
            [0.0, 0.46662021, 0.72311971],
            [0.46662021, 0.0, 0.62546036],
            [0.72311971, 0.62546036, 0.0],
        ]
        dist = dist_chord(self.mat_test)
        assert_allclose(dist, exp)

    def test_chisq_transform(self):
        """trans_chisq should return the exp result in the ref paper."""
        exp_m = [
            [0.42257713, 0.45643546, 0.84515425],
            [0.48294529, 0.7824608, 0.48294529],
            [1.01418511, 0.36514837, 0.3380617],
        ]
        res_m = trans_chisq(self.mat_test)
        assert_allclose(res_m, exp_m)

    def test_chisq_distance(self):
        """dist_chisq should return the exp result."""

        assert_allclose(dist_chisq(self.zeromtx), zeros((4, 4), "d"))

        exp_d = [
            [0.0, 0.4910521, 0.78452291],
            [0.4910521, 0.0, 0.69091002],
            [0.78452291, 0.69091002, 0.0],
        ]
        res_d = dist_chisq(self.mat_test)
        assert_allclose(res_d, exp_d)

    def test_hellinger_transform(self):
        """dist_hellinger should return the exp result in the ref paper."""
        exp = [
            [0.5, 0.5, 0.70710678],
            [0.53452248, 0.65465367, 0.53452248],
            [0.77459667, 0.4472136, 0.4472136],
        ]
        res = trans_hellinger(self.mat_test)
        assert_allclose(res, exp)

    def test_hellinger_distance(self):
        """dist_hellinger should return the exp result."""

        assert_allclose(dist_hellinger(self.zeromtx), zeros((4, 4), "d"))

        exp = [
            [0.0, 0.23429661, 0.38175149],
            [0.23429661, 0.0, 0.32907422],
            [0.38175149, 0.32907422, 0.0],
        ]
        dist = dist_hellinger(self.mat_test)
        assert_allclose(dist, exp)

    def test_species_profile_transform(self):
        """trans_specprof should return the exp result."""
        exp = [[0.25, 0.25, 0.5], [0.28571429, 0.42857143, 0.28571429], [0.6, 0.2, 0.2]]
        res = trans_specprof(self.mat_test)
        assert_allclose(res, exp)

    def test_species_profile_distance(self):
        """dist_specprof should return the exp result."""

        assert_allclose(dist_specprof(self.zeromtx), zeros((4, 4), "d"))

        exp = [
            [0.0, 0.28121457, 0.46368092],
            [0.28121457, 0.0, 0.39795395],
            [0.46368092, 0.39795395, 0.0],
        ]
        dist = dist_specprof(self.mat_test)
        assert_allclose(dist, exp)

    def test_dist_bray_curtis_magurran1(self):
        """zero values should return zero dist, or 1 with nonzero samples"""
        res = dist_bray_curtis_magurran(numpy.array([[0, 0, 0], [0, 0, 0], [1, 1, 1]]))
        assert_allclose(res, numpy.array([[0, 0, 1], [0, 0, 1], [1, 1, 0]]))

    def test_dist_bray_curtis_magurran2(self):
        """should match hand-calculated values"""
        res = dist_bray_curtis_magurran(numpy.array([[1, 4, 3], [1, 3, 5], [0, 2, 0]]))
        assert_allclose(
            res,
            numpy.array(
                [
                    [0, 1 - 14 / 17, 1 - (0.4)],
                    [1 - 14 / 17, 0, 1 - 4 / 11],
                    [1 - 0.4, 1 - 4 / 11, 0],
                ]
            ),
        )

    # def test_no_dupes(self):
    # """ here we check all distance functions in distance_transform for
    # duplicate
    # results.  Uses an unsafe hack to get all distance functions,
    # thus disabled by default
    # The dataset is from Legendre 2001, Ecologically Meaningful...
    # also, doesn't actually raise an error on failing, just prints to
    # stdout
    # """
    # import distance_transform
    # L19 dataset
    # L19data = array(
    # [[7,1,0,0,0,0,0,0,0],
    # [4,2,0,0,0,1,0,0,0],
    # [2,4,0,0,0,1,0,0,0],
    # [1,7,0,0,0,0,0,0,0],
    # [0,8,0,0,0,0,0,0,0],
    # [0,7,1,0,0,0,0,0,0],
    # [0,4,2,0,0,0,2,0,0],
    # [0,2,4,0,0,0,1,0,0],
    # [0,1,7,0,0,0,0,0,0],
    # [0,0,8,0,0,0,0,0,0],
    # [0,0,7,1,0,0,0,0,0],
    # [0,0,4,2,0,0,0,3,0],
    # [0,0,2,4,0,0,0,1,0],
    # [0,0,1,7,0,0,0,0,0],
    # [0,0,0,8,0,0,0,0,0],
    # [0,0,0,7,1,0,0,0,0],
    # [0,0,0,4,2,0,0,0,4],
    # [0,0,0,2,4,0,0,0,1],
    # [0,0,0,1,7,0,0,0,0]], 'd')

    # distfns = []
    # distfn_strs = dir(distance_transform)
    # warning: dangerous eval, and might catch bad or not functions
    # for fnstr in distfn_strs:
    # if fnstr.find('dist') != -1:
    # distfns.append(eval('%s' % fnstr))

    # dist_results = []
    # for distfn in distfns:
    # dist_results.append(distfn(L19data))
    # for i in range(len(dist_results)):
    # for j in range(i):
    # try:
    # assert_allclose(dist_results[i], dist_results[j])
    # except:
    # pass # should not be equal, so catch error and proceed
    # else:
    # print "duplicates found: ", distfns[i], distfns[j]


if __name__ == "__main__":
    main()
