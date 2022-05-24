#!/usr/bin/env python
"""Tests of the Enumeration and Alphabet objects.

Note: individual alphabets are typically in MolType and are tested there.
"""
import pickle

from unittest import TestCase, main

from numpy.testing import assert_equal

from cogent3.core.alphabet import (
    CharAlphabet,
    Enumeration,
    JointEnumeration,
    _make_complement_array,
    _make_translation_tables,
    array,
    get_array_type,
    uint8,
    uint16,
    uint32,
)
from cogent3.core.moltype import RNA, get_moltype


DnaBases = CharAlphabet("TCAG")
RnaBases = CharAlphabet("UCAG")
AminoAcids = CharAlphabet("ACDEFGHIKLMNPQRSTVWY")

__author__ = "Rob Knight, Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Rob Knight", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class translation_table_tests(TestCase):
    """Tests of top-level translation table functions"""

    def test_make_translation_tables(self):
        """_make_translation_tables should translate from chars to indices"""
        a = "ucag"
        itoa, atoi = _make_translation_tables(a)
        s = "ggacu"
        obs = s.translate(atoi)
        self.assertEqual(obs, "\x03\x03\x02\x01\x00")
        orig = obs.translate(itoa)
        self.assertEqual(orig, s)

    def test_make_complement_array(self):
        """_make_complement_array should identify complements correctly"""
        complement_array = _make_complement_array(RNA.alphabet, RNA.complements)
        test = "UCAG"
        test_array = [RNA.alphabet.index(i) for i in test]
        complements = complement_array.take(test_array)
        result = "".join([RNA.alphabet[i] for i in complements])
        self.assertEqual(result, "AGUC")


class get_array_type_tests(TestCase):
    """Tests of the get_array_type top-level function."""

    def test_get_array_type(self):
        """get_array_type should return unsigned type that fits elements."""
        self.assertEqual(get_array_type(0), uint8)
        self.assertEqual(get_array_type(100), uint8)
        self.assertEqual(get_array_type(256), uint8)  # boundary case
        self.assertEqual(get_array_type(257), uint16)  # boundary case
        self.assertEqual(get_array_type(10000), uint16)
        self.assertEqual(get_array_type(65536), uint16)
        self.assertEqual(get_array_type(65537), uint32)


class EnumerationTests(TestCase):
    """Tests of the Enumeration object."""

    def test_init(self):
        """Enumeration init should work from any sequence"""
        a = Enumeration("abc")
        self.assertEqual(a.index("a"), 0)
        self.assertEqual(a.index("b"), 1)
        self.assertEqual(a.index("c"), 2)
        self.assertEqual(a[0], "a")
        self.assertEqual(a[1], "b")
        self.assertEqual(a[2], "c")
        self.assertEqual(a.array_type, uint8)

        a = Enumeration("bca")
        self.assertEqual(a.index("b"), 0)
        self.assertEqual(a.index("c"), 1)
        self.assertEqual(a.index("a"), 2)
        self.assertEqual(a[0], "b")
        self.assertEqual(a[1], "c")
        self.assertEqual(a[2], "a")

        a = Enumeration([1, "2"])
        self.assertEqual(a.index(1), 0)
        self.assertEqual(a.index("2"), 1)
        self.assertRaises(KeyError, a.index, "1")

        # check that it works with gaps
        a = Enumeration("ab-", "-")
        self.assertEqual(a.gap, "-")
        self.assertEqual(a.gap_index, 2)

        a = Enumeration(list(range(257)))  # too big to fit in uint8
        self.assertEqual(a.array_type, uint16)

    def test_index(self):
        """Enumeration index should return first index of item"""
        a = Enumeration("bca")
        self.assertEqual(a.index("b"), 0)
        self.assertEqual(a.index("c"), 1)
        self.assertEqual(a.index("a"), 2)

    def test_getitem(self):
        """Enumeration[i] should return character at i"""
        a = Enumeration("bca")
        self.assertEqual(a[0], "b")
        self.assertEqual(a[1], "c")
        self.assertEqual(a[2], "a")

    def test_to_indices(self):
        """Enumeration to_indices should return indices from elements"""
        a = Enumeration("bca")
        self.assertEqual(a.to_indices(""), [])
        self.assertEqual(a.to_indices("ccabac"), [1, 1, 2, 0, 2, 1])

    def test_is_valid(self):
        """Enumeration is_valid should return True for valid sequence"""
        a = Enumeration("bca")
        self.assertEqual(a.is_valid(""), True)
        self.assertEqual(a.is_valid("bbb"), True)
        self.assertEqual(a.is_valid("bbbaac"), True)
        self.assertEqual(a.is_valid("bbd"), False)
        self.assertEqual(a.is_valid("d"), False)
        self.assertEqual(a.is_valid(["a", "b"]), True)
        self.assertEqual(a.is_valid(["a", None]), False)

    def test_from_indices(self):
        """Enumeration from_indices should return elements from indices"""
        a = Enumeration("bca")
        self.assertEqual(a.from_indices([]), [])
        self.assertEqual(a.from_indices([1, 1, 2, 0, 2, 1]), list("ccabac"))

    def test_pow(self):
        """Enumeration pow should produce JointEnumeration with n copies"""
        a = AminoAcids ** 3
        self.assertEqual(a[0], (AminoAcids[0],) * 3)
        self.assertEqual(a[-1], (AminoAcids[-1],) * 3)
        self.assertEqual(len(a), len(AminoAcids) ** 3)
        self.assertEqual(a.array_type, uint16)

        # check that it works with gaps
        a = Enumeration("a-b", "-")
        b = a ** 3
        self.assertEqual(len(b), 27)
        self.assertEqual(b.gap, ("-", "-", "-"))
        self.assertEqual(b.gap_index, 13)
        self.assertEqual(b.array_type, uint8)

        # check that array type is set correctly if needed
        b = a ** 6  # too big to fit in char
        self.assertEqual(b.array_type, uint16)

    def test_mul(self):
        """Enumeration mul should produce correct JointEnumeration"""
        a = DnaBases * RnaBases
        self.assertEqual(len(a), 16)
        self.assertEqual(a[0], ("T", "U"))
        self.assertEqual(a[-1], ("G", "G"))

        # check that it works with gaps
        a = Enumeration("ab-", "-")
        b = Enumeration("xz", "z")
        x = a * b
        self.assertEqual(x.gap, ("-", "z"))
        self.assertEqual(x.gap_index, 5)
        self.assertEqual(len(x), 6)
        self.assertEqual(
            x, (("a", "x"), ("a", "z"), ("b", "x"), ("b", "z"), ("-", "x"), ("-", "z"))
        )
        # check that it doesn't work when only one seq has gaps
        c = Enumeration("c")
        x = a * c
        self.assertEqual(x.gap, None)

    def test_counts(self):
        """Enumeration counts should count freqs in array"""
        a = DnaBases
        f = array([[0, 0, 1, 0, 0, 3]])
        assert_equal(a.counts(f), array([4, 1, 0, 1]))
        # check that it works with byte array
        f = array([[0, 0, 1, 0, 0, 3]], "B")
        assert_equal(a.counts(f), array([4, 1, 0, 1]))
        # should ignore out-of-bounds items
        g = [0, 4]
        assert_equal(a.counts(g), array([1, 0, 0, 0]))
        # make sure it works for long sequences, i.e. no wraparound at 255
        h = [0, 3] * 70000
        assert_equal(a.counts(h), array([70000, 0, 0, 70000]))
        h2 = array(h).astype("B")
        assert_equal(a.counts(h2), array([70000, 0, 0, 70000]))
        i = array([0, 3] * 75000)
        assert_equal(a.counts(i), array([75000, 0, 0, 75000]))
        # make sure it works for long _binary_ sequences, e.g. the results
        # of array comparisons.
        a = array([0, 1, 2, 3] * 10000)
        b = array([0, 0, 0, 0] * 10000)
        same = a == b


class CharAlphabetTests(TestCase):
    """Tests of CharAlphabets."""

    def test_init(self):
        """CharAlphabet init should make correct translation tables"""
        r = CharAlphabet("UCAG")
        i2c, c2i = r._indices_nums_to_chars, r._chars_to_indices
        s = array([0, 0, 1, 0, 3, 2], "b").tobytes()
        self.assertEqual(s.translate(i2c), b"UUCUGA")
        self.assertEqual("UUCUGA".translate(c2i), "\000\000\001\000\003\002")

    def test_pickling(self):
        r = CharAlphabet("UCAG")
        wa = r.get_word_alphabet(2)
        pkl = pickle.dumps(r)
        got = pickle.loads(pkl)
        self.assertIsInstance(got, type(r))
        self.assertEqual(got.get_word_alphabet(2), wa)

    def test_from_string(self):
        """CharAlphabet from_string should return correct array"""
        r = CharAlphabet("UCAG")
        assert_equal(r.from_string("UUCUGA"), array([0, 0, 1, 0, 3, 2], "B"))

    def test_is_valid(self):
        """CharAlphabet is_valid should return True for valid sequence"""
        a = CharAlphabet("bca")
        self.assertEqual(a.is_valid(""), True)
        self.assertEqual(a.is_valid("bbb"), True)
        self.assertEqual(a.is_valid("bbbaac"), True)
        self.assertEqual(a.is_valid("bbd"), False)
        self.assertEqual(a.is_valid("d"), False)
        self.assertEqual(a.is_valid(["a", "b"]), True)
        self.assertEqual(a.is_valid(["a", None]), False)

    def test_from_array(self):
        """CharAlphabet from_array should return correct array"""
        r = CharAlphabet("UCAG")
        got = r.from_array(array(["UUC", "UGA"], "c"))
        assert_equal(got, array([[0, 0, 1], [0, 3, 2]], "B"))

    def test_to_chars(self):
        """CharAlphabet to_chars should convert an input array to chars"""
        r = CharAlphabet("UCAG")
        c = r.to_chars(array([[0, 0, 1], [0, 3, 2]], "B"))
        assert_equal(c, array(["UUC", "UGA"], "c"))

    def test_to_string(self):
        """CharAlphabet to_string should convert an input array to string"""
        r = CharAlphabet("UCAG")
        self.assertEqual(r.to_string(array([[0, 0, 1], [0, 3, 2]], "B")), "UUC\nUGA")
        # should work with single seq
        self.assertEqual(r.to_string(array([[0, 0, 1, 0, 3, 2]], "B")), "UUCUGA")
        # should work with single seq
        self.assertEqual(r.to_string(array([0, 0, 1, 0, 3, 2], "B")), "UUCUGA")
        # should work with empty seq
        self.assertEqual(r.to_string(array([], "B")), "")

    def test_pairs(self):
        """pairs should cache the same object."""
        r = CharAlphabet("UCAG")
        rp = r.pairs
        self.assertEqual(len(rp), 16)
        rp2 = r.pairs
        self.assertIs(rp, rp2)

    def test_triples(self):
        """triples should cache the same object."""
        r = CharAlphabet("UCAG")
        rt = r.Triples
        self.assertEqual(len(rt), 64)
        rt2 = r.Triples
        self.assertIs(rt, rt2)

    def test_from_seq_to_array(self):
        """convert a sequence into indices"""
        dna = get_moltype("dna")
        seq = dna.make_seq("ACGG")
        got = dna.alphabet.from_seq_to_array(seq)
        assert_equal(got, array([dna.alphabet.index(b) for b in seq]))

    def test_from_ordinals_to_seq(self):
        """check indices convert to a sequence"""
        indices = [2, 1, 3, 3]
        dna = get_moltype("dna")
        got = dna.alphabet.from_ordinals_to_seq(indices)
        self.assertEqual(got, dna.make_seq("ACGG"))


class JointEnumerationTests(TestCase):
    """Tests of JointEnumerations."""

    def test_init(self):
        """JointEnumeration init should work as expected"""
        # should work for alphabet object
        a = JointEnumeration([DnaBases, RnaBases])
        self.assertEqual(len(a), 16)
        self.assertEqual(a.shape, (4, 4))
        self.assertEqual(a[0], ("T", "U"))
        self.assertEqual(a[-1], ("G", "G"))
        assert_equal(a._sub_enum_factors, array([[4], [1]]))

        # should work for arbitrary sequences
        a = JointEnumeration(["TCAG", "UCAG"])
        self.assertEqual(len(a), 16)
        self.assertEqual(a[0], ("T", "U"))
        self.assertEqual(a[-1], ("G", "G"))
        assert_equal(a._sub_enum_factors, array([[4], [1]]))

        # should work for different length sequences
        a = JointEnumeration(["TCA", "UCAG"])
        self.assertEqual(a.shape, (3, 4))
        self.assertEqual(len(a), 12)
        self.assertEqual(a[0], ("T", "U"))
        self.assertEqual(a[-1], ("A", "G"))
        assert_equal(a._sub_enum_factors, array([[4], [1]]))  # note: _not_ [3,1]

    def test_to_indices(self):
        """JointEnumeration to_indices should convert tuples correctly"""
        a = JointEnumeration(["TCAG", "UCAG"])
        i = a.to_indices([("T", "U"), ("G", "G"), ("G", "G")])
        self.assertEqual(i, [0, 15, 15])

    def test_from_indices(self):
        """JointEnumeration from_indices should return correct tuples"""
        a = JointEnumeration(["TCAG", "UCAG"])
        i = a.from_indices([0, 15, 15])
        self.assertEqual(i, [("T", "U"), ("G", "G"), ("G", "G")])

    def test_pack_arrays(self):
        """JointEnumeration pack_arrays should return correct array."""
        a = JointEnumeration(["xyz", "abcd", "ef"])
        v = [[0, 1, 2, 0], [3, 3, 1, 0], [1, 1, 0, 0]]
        result = a.pack_arrays(v)
        assert_equal(result, array([7, 15, 18, 0]))

    def test_unpack_arrays(self):
        """JointEnumeration unpack_arrays should return correct arrays."""
        a = JointEnumeration(["xyz", "abcd", "ef"])
        v = [7, 15, 18, 0]
        result = a.unpack_arrays(v)
        assert_equal(result, array([[0, 1, 2, 0], [3, 3, 1, 0], [1, 1, 0, 0]]))


if __name__ == "__main__":
    main()
