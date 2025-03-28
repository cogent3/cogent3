"""Tests of the Enumeration and Alphabet objects.

Note: individual alphabets are typically in MolType and are tested there.
"""

import os
import pickle
from unittest import TestCase

import numpy
import pytest
from numpy import unravel_index
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
    uint64,
)
from cogent3.core.moltype import RNA

if "COGENT3_NEW_TYPE" in os.environ:
    pytest.skip(
        "Tests skipped because COGENT3_TYPE environment variable is defined",
        allow_module_level=True,
    )

DnaBases = CharAlphabet("TCAG")
RnaBases = CharAlphabet("UCAG")
AminoAcids = CharAlphabet("ACDEFGHIKLMNPQRSTVWY")


class translation_table_tests(TestCase):
    """Tests of top-level translation table functions"""

    def test_make_translation_tables(self):
        """_make_translation_tables should translate from chars to indices"""
        a = "ucag"
        itoa, atoi = _make_translation_tables(a)
        s = "ggacu"
        obs = s.translate(atoi)
        assert obs == "\x03\x03\x02\x01\x00"
        orig = obs.translate(itoa)
        assert orig == s

    def test_make_complement_array(self):
        """_make_complement_array should identify complements correctly"""
        complement_array = _make_complement_array(RNA.alphabet, RNA.complements)
        test = "UCAG"
        test_array = [RNA.alphabet.index(i) for i in test]
        complements = complement_array.take(test_array)
        result = "".join([RNA.alphabet[i] for i in complements])
        assert result == "AGUC"


class get_array_type_tests(TestCase):
    """Tests of the get_array_type top-level function."""

    def test_get_array_type(self):
        """get_array_type should return unsigned type that fits elements."""
        assert get_array_type(0) == uint8
        assert get_array_type(100) == uint8
        assert get_array_type(255) == uint8  # boundary case
        assert get_array_type(257) == uint16  # boundary case
        assert get_array_type(10000) == uint16
        assert get_array_type(65535) == uint16
        assert get_array_type(65537) == uint32
        assert get_array_type(1 + 2**32) == uint64

    def test_get_array_type_fail(self):
        """get_array_type should return unsigned type that fits elements."""
        with pytest.raises(NotImplementedError):
            assert get_array_type(2**64) == uint64


class EnumerationTests(TestCase):
    """Tests of the Enumeration object."""

    def test_init(self):
        """Enumeration init should work from any sequence"""
        a = Enumeration("abc")
        assert a.index("a") == 0
        assert a.index("b") == 1
        assert a.index("c") == 2
        assert a[0] == "a"
        assert a[1] == "b"
        assert a[2] == "c"
        assert a.array_type == uint8

        a = Enumeration("bca")
        assert a.index("b") == 0
        assert a.index("c") == 1
        assert a.index("a") == 2
        assert a[0] == "b"
        assert a[1] == "c"
        assert a[2] == "a"

        a = Enumeration([1, "2"])
        assert a.index(1) == 0
        assert a.index("2") == 1
        self.assertRaises(KeyError, a.index, "1")

        # check that it works with gaps
        a = Enumeration("ab-", "-")
        assert a.gap == "-"
        assert a.gap_index == 2

        a = Enumeration(list(range(257)))  # too big to fit in uint8
        assert a.array_type == uint16

    def test_index(self):
        """Enumeration index should return first index of item"""
        a = Enumeration("bca")
        assert a.index("b") == 0
        assert a.index("c") == 1
        assert a.index("a") == 2

    def test_getitem(self):
        """Enumeration[i] should return character at i"""
        a = Enumeration("bca")
        assert a[0] == "b"
        assert a[1] == "c"
        assert a[2] == "a"

    def test_to_indices(self):
        """Enumeration to_indices should return indices from elements"""
        a = Enumeration("bca")
        assert a.to_indices("") == []
        assert a.to_indices("ccabac") == [1, 1, 2, 0, 2, 1]

    def test_is_valid(self):
        """Enumeration is_valid should return True for valid sequence"""
        a = Enumeration("bca")
        assert a.is_valid("") is True
        assert a.is_valid("bbb") is True
        assert a.is_valid("bbbaac") is True
        assert a.is_valid("bbd") is False
        assert a.is_valid("d") is False
        assert a.is_valid(["a", "b"]) is True
        assert a.is_valid(["a", None]) is False

    def test_from_indices(self):
        """Enumeration from_indices should return elements from indices"""
        a = Enumeration("bca")
        assert a.from_indices([]) == []
        assert a.from_indices([1, 1, 2, 0, 2, 1]) == list("ccabac")

    def test_pow(self):
        """Enumeration pow should produce JointEnumeration with n copies"""
        a = AminoAcids**3
        assert a[0] == (AminoAcids[0],) * 3
        assert a[-1] == (AminoAcids[-1],) * 3
        assert len(a) == len(AminoAcids) ** 3
        assert a.array_type == uint16

        # check that it works with gaps
        a = Enumeration("a-b", "-")
        b = a**3
        assert len(b) == 27
        assert b.gap == ("-", "-", "-")
        assert b.gap_index == 13
        assert b.array_type == uint8

        # check that array type is set correctly if needed
        b = a**6  # too big to fit in char
        assert b.array_type == uint16

    def test_mul(self):
        """Enumeration mul should produce correct JointEnumeration"""
        a = DnaBases * RnaBases
        assert len(a) == 16
        assert a[0] == ("T", "U")
        assert a[-1] == ("G", "G")

        # check that it works with gaps
        a = Enumeration("ab-", "-")
        b = Enumeration("xz", "z")
        x = a * b
        assert x.gap == ("-", "z")
        assert x.gap_index == 5
        assert len(x) == 6
        assert x == (
            ("a", "x"),
            ("a", "z"),
            ("b", "x"),
            ("b", "z"),
            ("-", "x"),
            ("-", "z"),
        )
        # check that it doesn't work when only one seq has gaps
        c = Enumeration("c")
        x = a * c
        assert x.gap is None


class CharAlphabetTests(TestCase):
    """Tests of CharAlphabets."""

    def test_init(self):
        """CharAlphabet init should make correct translation tables"""
        r = CharAlphabet("UCAG")
        i2c, c2i = r._indices_nums_to_chars, r._chars_to_indices
        s = array([0, 0, 1, 0, 3, 2], "b").tobytes()
        assert s.translate(i2c) == b"UUCUGA"
        assert "UUCUGA".translate(c2i) == "\x00\x00\x01\x00\x03\x02"

    def test_pickling(self):
        r = CharAlphabet("UCAG")
        wa = r.get_word_alphabet(2)
        pkl = pickle.dumps(r)
        got = pickle.loads(pkl)
        assert isinstance(got, type(r))
        assert got.get_word_alphabet(2) == wa

    def test_pickling_moltype(self):
        from cogent3.core.moltype import DNA

        a = DNA.alphabet
        got = pickle.loads(pickle.dumps(a))
        assert got.moltype is not None

    def test_word_alphabet_order(self):
        bases = "TCAG"
        r = CharAlphabet(bases)
        indices = [unravel_index(i, shape=(4, 4, 4)) for i in range(64)]
        expect = tuple("".join([bases[b] for b in coord]) for coord in indices)
        got = tuple(r.get_word_alphabet(3))
        assert got == expect

    def test_is_valid(self):
        """CharAlphabet is_valid should return True for valid sequence"""
        a = CharAlphabet("bca")
        assert a.is_valid("")
        assert a.is_valid("bbb")
        assert a.is_valid("bbbaac")
        assert not a.is_valid("bbd")
        assert not a.is_valid("d")
        assert a.is_valid(["a", "b"])
        assert not a.is_valid(["a", None])

    def test_to_chars(self):
        """CharAlphabet to_chars should convert an input array to chars"""
        r = CharAlphabet("UCAG")
        c = r.to_chars(array([[0, 0, 1], [0, 3, 2]], "B"))
        assert_equal(c, array(["UUC", "UGA"], "c"))

    def test_to_string(self):
        """CharAlphabet to_string should convert an input array to string"""
        r = CharAlphabet("UCAG")
        assert r.to_string(array([[0, 0, 1], [0, 3, 2]], "B")) == "UUC\nUGA"
        # should work with single seq
        assert r.to_string(array([[0, 0, 1, 0, 3, 2]], "B")) == "UUCUGA"
        # should work with single seq
        assert r.to_string(array([0, 0, 1, 0, 3, 2], "B")) == "UUCUGA"
        # should work with empty seq
        assert r.to_string(array([], "B")) == ""


class JointEnumerationTests(TestCase):
    """Tests of JointEnumerations."""

    def test_init(self):
        """JointEnumeration init should work as expected"""
        # should work for alphabet object
        a = JointEnumeration([DnaBases, RnaBases])
        assert len(a) == 16
        assert a.shape == (4, 4)
        assert a[0] == ("T", "U")
        assert a[-1] == ("G", "G")
        assert_equal(a._sub_enum_factors, array([[4], [1]]))

        # should work for arbitrary sequences
        a = JointEnumeration(["TCAG", "UCAG"])
        assert len(a) == 16
        assert a[0] == ("T", "U")
        assert a[-1] == ("G", "G")
        assert_equal(a._sub_enum_factors, array([[4], [1]]))

        # should work for different length sequences
        a = JointEnumeration(["TCA", "UCAG"])
        assert a.shape == (3, 4)
        assert len(a) == 12
        assert a[0] == ("T", "U")
        assert a[-1] == ("A", "G")
        assert_equal(a._sub_enum_factors, array([[4], [1]]))  # note: _not_ [3,1]

    def test_to_indices(self):
        """JointEnumeration to_indices should convert tuples correctly"""
        a = JointEnumeration(["TCAG", "UCAG"])
        i = a.to_indices([("T", "U"), ("G", "G"), ("G", "G")])
        assert i == [0, 15, 15]

    def test_from_indices(self):
        """JointEnumeration from_indices should return correct tuples"""
        a = JointEnumeration(["TCAG", "UCAG"])
        i = a.from_indices([0, 15, 15])
        assert i == [("T", "U"), ("G", "G"), ("G", "G")]

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


@pytest.mark.parametrize("data", ["AUGGA", b"AUGGA"])
def test_to_indices(data):
    alphabet = RNA.alphabets.degen_gapped
    got = alphabet.to_indices(data)
    assert got.tolist() == [2, 0, 3, 3, 2]


@pytest.mark.parametrize("cls", [list, numpy.array, tuple])
def test_from_indices(cls):
    data = cls([2, 0, 3, 3, 2])
    alphabet = RNA.alphabets.degen_gapped
    got = alphabet.from_indices(data)
    assert got == "AUGGA"
