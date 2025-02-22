#!/usr/bin/env python
"""Tests of transformation and composition functions ."""

from unittest import TestCase

from numpy.testing import assert_allclose

from cogent3.util.transform import (
    KeepChars,
    first_index_in_set,
    for_seq,
    per_longest,
    per_shortest,
)


class has_x:
    # convenience class for has_field and related functions

    def __init__(self, x):
        self.x = x

    def __hash__(self):
        return hash(self.x)

    def __str__(self):
        return str(self.x)


class has_y:
    # convenience class for has_field and related functions

    def __init__(self, y):
        self.y = y

    def __hash__(self):
        return hash(self.y)

    def __str__(self):
        return str(self.y)


class metafunctionsTests(TestCase):
    """Tests of standalone functions."""

    def setUp(self):
        """Define some standard functions and data."""
        self.Numbers = list(range(20))
        self.SmallNumbers = list(range(3))
        self.SmallNumbersRepeated = list(range(5)) * 4
        self.Letters = "abcde"
        self.Mixed = list(self.Letters) + list(range(5))
        self.firsts = "ab2"
        self.seconds = "0bc"

        self.is_char = lambda x: isinstance(x, str) and len(x) == 1
        self.is_vowel = lambda x: x in "aeiou"
        self.is_consonant = lambda x: x not in "aeiuo"
        self.is_number = lambda x: isinstance(x, int)
        self.is_odd_number = lambda x: x % 2
        self.is_odd_letter = lambda x: x in "acegikmoqs"
        self.is_zero = lambda x: x == 0
        self.is_small = lambda x: x < 3
        self.double = lambda x: x * 2
        self.minusone = lambda x: x - 1

        # function to test *args, **kwargs)
        self.is_alpha_digit = lambda first, second: first.isalpha() and second.isdigit()
        self.is_digit_alpha = lambda first, second: first.isdigit() and second.isalpha()


class SequenceFunctionsTests(TestCase):
    """Tests of standalone functions for dealing with sequences."""

    def test_per_shortest(self):
        """per_shortest should divide by min(len(x), len(y))"""
        assert per_shortest(20, "aaaaaa", "bbbb") == 5
        assert per_shortest(20, "aaaaaa", "b") == 20
        assert per_shortest(20, "a", "bbbbb") == 20
        assert per_shortest(20, "", "b") == 0
        assert per_shortest(20, "", "") == 0
        # check that it does it in floating-point
        assert per_shortest(1, "aaaaaa", "bbbb") == 0.25
        # check that it raises TypeError on non-seq
        self.assertRaises(TypeError, per_shortest, 1, 2, 3)

    def test_per_longest(self):
        """per_longest should divide by max(len(x), len(y))"""
        assert per_longest(20, "aaaaaa", "bbbb") == 20 / 6.0
        assert per_longest(20, "aaaaaa", "b") == 20 / 6.0
        assert per_longest(20, "a", "bbbbb") == 20 / 5.0
        assert per_longest(20, "", "b") == 20
        assert per_longest(20, "", "") == 0
        # check that it does it in floating-point
        assert per_longest(1, "aaaaaa", "bbbb") == 1 / 6.0
        # check that it raises TypeError on non-seq
        self.assertRaises(TypeError, per_longest, 1, 2, 3)

    def test_for_seq(self):
        """for_seq should return the correct function"""

        def is_eq(x, y):
            return x == y

        def is_ne(x, y):
            return x != y

        def lt_5(x, y):
            return x + y < 5

        def diff(x, y):
            return x - y

        def sumsq(x):
            return sum([i * i for i in x])

        def long_norm(s, x, y):
            return (s + 0.0) / max(len(x), len(y))

        def times_two(s, x, y):
            return 2 * s

        s1 = [1, 2, 3, 4, 5]
        s2 = [1, 3, 2, 4, 5]
        s3 = [1, 1, 1, 1, 1]
        s4 = [5, 5, 5, 5, 5]
        s5 = [3, 3, 3, 3, 3]
        short = [1]

        # test behavior with default aggregator and normalizer
        f = for_seq(is_eq)
        assert_allclose(f(s1, s1), 1.0)
        assert_allclose(f(s1, short), 1.0)
        assert_allclose(f(short, s1), 1.0)
        assert_allclose(f(short, s4), 0.0)
        assert_allclose(f(s4, short), 0.0)
        assert_allclose(f(s1, s2), 0.6)

        f = for_seq(is_ne)
        assert_allclose(f(s1, s1), 0.0)
        assert_allclose(f(s1, short), 0.0)
        assert_allclose(f(short, s1), 0.0)
        assert_allclose(f(short, s4), 1.0)
        assert_allclose(f(s4, short), 1.0)
        assert_allclose(f(s1, s2), 0.4)

        f = for_seq(lt_5)
        assert_allclose(f(s3, s3), 1.0)
        assert_allclose(f(s3, s4), 0.0)
        assert_allclose(f(s2, s3), 0.6)

        f = for_seq(diff)
        assert_allclose(f(s1, s1), 0.0)
        assert_allclose(f(s4, s1), 2.0)
        assert_allclose(f(s1, s4), -2.0)

        # test behavior with different aggregator
        f = for_seq(diff)
        assert_allclose(f(s1, s5), 0)
        f = for_seq(diff, aggregator=sum)
        assert_allclose(f(s1, s5), 0)
        f = for_seq(diff, aggregator=sumsq)
        assert_allclose(f(s1, s5), 2.0)

        # test behavior with different normalizer
        f = for_seq(diff, aggregator=sumsq, normalizer=None)
        assert_allclose(f(s1, s5), 10)
        f = for_seq(diff, aggregator=sumsq)
        assert_allclose(f(s1, s5), 2.0)
        f = for_seq(diff, aggregator=sumsq, normalizer=times_two)
        assert_allclose(f(s1, s5), 20)
        f = for_seq(diff, aggregator=sumsq)
        assert_allclose(f(s5, short), 4)
        f = for_seq(diff, aggregator=sumsq, normalizer=long_norm)
        assert_allclose(f(s5, short), 0.8)


class Filter_Criteria_Tests(TestCase):
    """Tests of standalone functions used as filter criteria"""

    def test_KeepChars(self):
        """KeepChars returns a string containing only chars in keep"""
        f = KeepChars("ab c3*[")
        assert f("") == ""  # empty
        self.assertRaises(TypeError, f, None)  # None

        # one character, case sensitive
        assert f("b") == "b"
        assert f("g") == ""
        assert f("xyz123") == "3"
        assert f("xyz  123") == "  3"

        # more characters, case sensitive
        assert f("kjbwherzcagebcujrkcs") == "bcabcc"
        assert f("f[ffff*ff*fff3fff") == "[**3"

        # case insensitive
        f = KeepChars("AbC", False)
        assert f("abcdef") == "abc"
        assert f("ABCDEF") == "ABC"
        assert f("aBcDeF") == "aBc"

    def test_first_index_in_set(self):
        """first_index_in_set should return index of first occurrence"""
        vowels = "aeiou"
        s1 = "ebcua"
        s2 = "bcbae"
        s3 = ""
        s4 = "cbd"
        assert first_index_in_set(s1, vowels) == 0
        assert first_index_in_set(s2, vowels) == 3
        assert first_index_in_set(s3, vowels) is None
        assert first_index_in_set(s4, vowels) is None
