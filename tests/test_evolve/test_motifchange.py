#!/usr/bin/env python

import unittest

from cogent3.core.moltype import CodonAlphabet
from cogent3.evolve.predicate import MotifChange, parse


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Rob Knight",
    "Matthew Wakefield",
    "Brett Easton",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class FakeModel(object):
    def __init__(self, alphabet):
        self.alphabet = alphabet
        self.moltype = alphabet.moltype

    def get_alphabet(self):
        return self.alphabet


class TestPredicates(unittest.TestCase):
    def setUp(self):
        self.alphabet = CodonAlphabet()
        self.model = FakeModel(self.alphabet)

    def _makeMotifChange(self, *args, **kw):
        pred = MotifChange(*args, **kw)
        return pred.interpret(self.model)

    def test_parse(self):
        """correctly construction"""
        ag = MotifChange("A", "G")
        got = parse(str(ag))
        self.assertEqual(str(got), "A/G")
        ts = MotifChange("A", "G") | MotifChange("C", "T")
        got = parse(str(ts))
        self.assertEqual(str(got), "(A/G | C/T)")
        a_g = MotifChange("A", "G", forward_only=True)
        t_c = MotifChange("T", "C", forward_only=True)
        sym = a_g | t_c
        got = parse(str(sym))
        self.assertEqual(str(got), "(A>G | T>C)")

    def assertMatch(self, pred, seq1, seq2):
        assert pred(seq1, seq2), (pred.__doc__, (seq1, seq2))

    def assertNoMatch(self, pred, seq1, seq2):
        assert not pred(seq1, seq2), ("not " + pred.__doc__, (seq1, seq2))

    def test_indels(self):
        indel = self._makeMotifChange("---", "NNN")
        self.assertMatch(indel, "---", "AAA")

    def test_impossible_change(self):
        self.assertRaises(Exception, self._makeMotifChange, "----", "NNNN")

    def test_isfromcpg(self):
        isFromCpG = self._makeMotifChange("CG", forward_only=True)
        self.assertMatch(isFromCpG, "CG", "CA")
        self.assertMatch(isFromCpG, "CG", "TG")
        self.assertMatch(isFromCpG, "ACG", "ATG")
        self.assertMatch(isFromCpG, "CGT", "CTT")

        self.assertNoMatch(isFromCpG, "CTT", "CGT")
        self.assertNoMatch(isFromCpG, "C", "G")

    def test_isfromtocpg(self):
        isFromToCpG = self._makeMotifChange("CG")
        self.assertMatch(isFromToCpG, "CG", "CA")
        self.assertMatch(isFromToCpG, "CG", "TG")
        self.assertMatch(isFromToCpG, "ACG", "ATG")
        self.assertMatch(isFromToCpG, "CGT", "CTT")
        self.assertMatch(isFromToCpG, "CTT", "CGT")

    def test_isFromToCpA_C_only(self):
        isFromToCpA_C_only = self._makeMotifChange("CA", diff_at=0)
        self.assertMatch(isFromToCpA_C_only, "CA", "TA")
        self.assertMatch(isFromToCpA_C_only, "TCA", "TTA")
        self.assertMatch(isFromToCpA_C_only, "TAA", "CAA")
        self.assertNoMatch(isFromToCpA_C_only, "TCA", "TCT")

    def test_isFromCpA_C_only(self):
        isFromCpA_C_only = self._makeMotifChange("CA", forward_only=True, diff_at=0)
        self.assertMatch(isFromCpA_C_only, "CA", "TA")
        self.assertMatch(isFromCpA_C_only, "TCA", "TTA")
        self.assertNoMatch(isFromCpA_C_only, "TAA", "CAA")

    def test_isCpT_T_only(self):
        isCpT_T_only = self._makeMotifChange("CT", diff_at=1)
        self.assertMatch(isCpT_T_only, "CT", "CA")
        self.assertMatch(isCpT_T_only, "TCA", "TCT")
        self.assertNoMatch(isCpT_T_only, "TTA", "TCA")
        self.assertNoMatch(isCpT_T_only, "TA", "CT")

    def test_isCCC(self):
        isCCC = self._makeMotifChange("CCC")
        self.assertNoMatch(isCCC, "CC", "CT")

    def test_isC(self):
        isC = self._makeMotifChange("C")
        self.assertMatch(isC, "C", "T")
        self.assertNoMatch(isC, "CA", "CT")
        self.assertMatch(isC, "CA", "CC")
        self.assertMatch(isC, "CAT", "GAT")
        self.assertMatch(isC, "CAT", "CCT")
        self.assertMatch(isC, "CAT", "CAC")
        self.assertNoMatch(isC, "CAT", "CAA")
        self.assertNoMatch(isC, "C", "C")

    def test_isCtoT(self):
        isCtoT = self._makeMotifChange("C", "T")
        self.assertMatch(isCtoT, "C", "T")
        self.assertMatch(isCtoT, "T", "C")
        self.assertNoMatch(isCtoT, "T", "A")
        isCtoT = self._makeMotifChange("C", "T", forward_only=True)
        self.assertMatch(isCtoT, "C", "T")
        self.assertNoMatch(isCtoT, "T", "C")

    def test_isCGtoCA(self):
        isCG_CA = self._makeMotifChange("CG", "CA")
        self.assertMatch(isCG_CA, "CG", "CA")
        self.assertMatch(isCG_CA, "CA", "CG")
        self.assertMatch(isCG_CA, "CAT", "CGT")
        self.assertMatch(isCG_CA, "CGT", "CAT")
        self.assertMatch(isCG_CA, "TCA", "TCG")
        self.assertNoMatch(isCG_CA, "TCT", "TCG")
        self.assertMatch(isCG_CA, "CGTT", "CATT")
        self.assertMatch(isCG_CA, "TCGT", "TCAT")
        self.assertMatch(isCG_CA, "TTCG", "TTCA")
        self.assertMatch(isCG_CA, "CATT", "CGTT")
        self.assertMatch(isCG_CA, "TCAT", "TCGT")
        self.assertMatch(isCG_CA, "TTCA", "TTCG")
        isCG_CA = self._makeMotifChange("CG", "CA", forward_only=True)
        self.assertMatch(isCG_CA, "CGTT", "CATT")
        self.assertMatch(isCG_CA, "TCGT", "TCAT")
        self.assertMatch(isCG_CA, "TTCG", "TTCA")
        self.assertNoMatch(isCG_CA, "CATT", "CGTT")
        self.assertNoMatch(isCG_CA, "TCAT", "TCGT")
        self.assertNoMatch(isCG_CA, "TTCA", "TTCG")

        isCG = self._makeMotifChange("CG", diff_at=1)
        self.assertMatch(isCG, "CGTT", "CATT")
        self.assertMatch(isCG, "TCGT", "TCAT")
        self.assertMatch(isCG, "TTCG", "TTCA")
        self.assertNoMatch(isCG, "CGTT", "TGTT")
        self.assertNoMatch(isCG, "TCGT", "TAGT")
        self.assertNoMatch(isCG, "TTCG", "--GG")

    def test_wildcards(self):
        isCG_CN = self._makeMotifChange("CG", "CN")
        self.assertMatch(isCG_CN, "CG", "CA")
        self.assertNoMatch(isCG_CN, "CG", "CG")
        self.assertNoMatch(isCG_CN, "CG", "C-")


if __name__ == "__main__":
    unittest.main()
