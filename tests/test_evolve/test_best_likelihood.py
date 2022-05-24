#!/usr/bin/env python

import math

from unittest import TestCase, main

from cogent3 import DNA, make_aligned_seqs
from cogent3.evolve.best_likelihood import (
    BestLogLikelihood,
    _take,
    _transpose,
    aligned_columns_to_rows,
    count_column_freqs,
    get_G93_lnL_from_array,
    get_ML_probs,
)


__author__ = "Helen Lindsay"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Helen Lindsay"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Helen Lindsay"
__email__ = "helen.lindsay@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose


IUPAC_DNA_ambiguities = "NRYWSKMBDHV"


def makeSampleAlignment(gaps=False, ambiguities=False):
    if gaps:
        seqs_list = ["AAA--CTTTGG-T", "CCCCC-TATG-GT", "-AACCCTTTGGGT"]
    elif ambiguities:
        seqs_list = ["AARNCCTTTGGC", "CCNYCCTTTGSG", "CAACCCTGWGGG"]
    else:
        seqs_list = ["AAACCCGGGTTTA", "CCCGGGTTTAAAC", "GGGTTTAAACCCG"]
    seqs = list(zip("abc", seqs_list))
    return make_aligned_seqs(data=seqs)


class TestGoldman93(TestCase):
    def setUp(self):
        self.aln = makeSampleAlignment()
        self.gapped_aln = makeSampleAlignment(gaps=True)
        self.ambig_aln = makeSampleAlignment(ambiguities=True)

    def test_aligned_columns_to_rows(self):
        obs = aligned_columns_to_rows(self.aln[:-1], 3)
        expect = [
            ["AAA", "CCC", "GGG"],
            ["CCC", "GGG", "TTT"],
            ["GGG", "TTT", "AAA"],
            ["TTT", "AAA", "CCC"],
        ]
        assert obs == expect, (obs, expect)

        obs = aligned_columns_to_rows(self.aln, 1)
        expect = [
            ["A", "C", "G"],
            ["A", "C", "G"],
            ["A", "C", "G"],
            ["C", "G", "T"],
            ["C", "G", "T"],
            ["C", "G", "T"],
            ["G", "T", "A"],
            ["G", "T", "A"],
            ["G", "T", "A"],
            ["T", "A", "C"],
            ["T", "A", "C"],
            ["T", "A", "C"],
            ["A", "C", "G"],
        ]
        self.assertEqual(obs, expect)
        obs = aligned_columns_to_rows(self.gapped_aln[:-1], 3, allowed_chars="ACGT")
        expect = [["TTT", "TAT", "TTT"]]
        self.assertEqual(obs, expect)

        obs = aligned_columns_to_rows(
            self.ambig_aln, 2, exclude_chars=IUPAC_DNA_ambiguities
        )
        expect = [["AA", "CC", "CA"], ["CC", "CC", "CC"], ["TT", "TT", "TG"]]
        self.assertEqual(obs, expect)

    def test_count_column_freqs(self):
        columns = aligned_columns_to_rows(self.aln, 1)
        obs = count_column_freqs(columns)
        expect = {"A C G": 4, "C G T": 3, "G T A": 3, "T A C": 3}
        self.assertEqual(obs, expect)

        columns = aligned_columns_to_rows(self.aln[:-1], 2)
        obs = count_column_freqs(columns)
        expect = {
            "AA CC GG": 1,
            "AC CG GT": 1,
            "CC GG TT": 1,
            "GG TT AA": 1,
            "GT TA AC": 1,
            "TT AA CC": 1,
        }
        self.assertEqual(obs, expect)

    def test__transpose(self):
        """test transposing an array"""
        a = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
        e = [[0, 3, 6, 9], [1, 4, 7, 10], [2, 5, 8, 11]]
        self.assertEqual(_transpose(a), e)

    def test__take(self):
        """test taking selected rows from an array"""
        e = [[0, 3, 6, 9], [1, 4, 7, 10], [2, 5, 8, 11]]
        self.assertEqual(_take(e, [0, 1]), [[0, 3, 6, 9], [1, 4, 7, 10]])
        self.assertEqual(_take(e, [1, 2]), [[1, 4, 7, 10], [2, 5, 8, 11]])
        self.assertEqual(_take(e, [0, 2]), [[0, 3, 6, 9], [2, 5, 8, 11]])

    def test_get_ML_probs(self):
        columns = aligned_columns_to_rows(self.aln, 1)
        obs = get_ML_probs(columns, with_patterns=True)
        expect = {
            "A C G": 4 / 13.0,
            "C G T": 3 / 13.0,
            "G T A": 3 / 13.0,
            "T A C": 3 / 13.0,
        }
        sum = 0
        for pattern, lnL, freq in obs:
            assert_allclose(lnL, expect[pattern])
            sum += lnL
            self.assertTrue(lnL >= 0)
        assert_allclose(sum, 1)

    def test_get_G93_lnL_from_array(self):
        columns = aligned_columns_to_rows(self.aln, 1)
        obs = get_G93_lnL_from_array(columns)
        expect = math.log(math.pow(4 / 13.0, 4)) + 3 * math.log(math.pow(3 / 13.0, 3))
        assert_allclose(obs, expect)

    def test_BestLogLikelihood(self):
        obs = BestLogLikelihood(self.aln, DNA.alphabet)
        expect = math.log(math.pow(4 / 13.0, 4)) + 3 * math.log(math.pow(3 / 13.0, 3))
        assert_allclose(obs, expect)
        lnL, l = BestLogLikelihood(self.aln, DNA.alphabet, return_length=True)
        self.assertEqual(l, len(self.aln))


if __name__ == "__main__":
    main()
