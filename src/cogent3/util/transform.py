#!/usr/bin/env python
"""Provides transformations of functions and other objects.

Includes:

Standard combinatorial higher-order functions adapted from David Mertz (2003),
"Text Processing in Python", Chapter 1.

Functions for performing complex tests on strings, e.g. includes_any or
includes_all.

Functions for generating combinations, permutations, or cartesian products
of lists.
"""

maketrans = str.maketrans

# standard combinatorial HOF's from Mertz


def per_shortest(total, x, y):
    """Divides total by min(len(x), len(y)).

    Useful for normalizing per-item results from sequences that are zipped
    together. Always returns 0 if one of the sequences is empty (to
    avoid divide by zero error).
    """
    shortest = min(len(x), len(y))
    if not shortest:
        return 0
    return total / shortest


def per_longest(total, x, y):
    """Divides total by max(len(x), len(y)).

    Useful for normalizing per-item results from sequences that are zipped
    together. Always returns 0 if one of the sequences is empty (to
    avoid divide by zero error).
    """
    longest = max(len(x), len(y))
    if not longest:
        return 0
    return total / longest


class for_seq:
    """Returns function that applies f(i,j) to i,j in zip(first, second).

    f: f(i,j) applying to elements of the sequence.

    aggregator: method to reduce the list of results to a scalar. Default: sum.

    normalizer: f(total, i, j) that normalizes the total as a function of
    i and j. Default is length_normalizer (divides by the length of the shorter
    of i and j). If normalizer is None, no normalization is performed.

    Will always truncate to length of the shorter sequence (because of the use
    of zip).
    """

    def __init__(self, f, aggregator=sum, normalizer=per_shortest):
        self.f = f
        self.aggregator = aggregator
        self.normalizer = normalizer

    def __call__(self, first, second):
        f = self.f
        if self.normalizer is None:
            return self.aggregator(
                [f(i, j) for i, j in zip(first, second, strict=False)],
            )
        return self.normalizer(
            self.aggregator([f(i, j) for i, j in zip(first, second, strict=False)]),
            first,
            second,
        )


# convenience functions for modifying objects


class KeepChars:
    """Returns a filter object o(s): call to return a filtered string.

    Specifically, strips out everything in s that is not in keep.
    This filter is case sensitive by default.
    """

    allchars = bytes(range(256))

    def __init__(self, keep, case_sens=True):
        """Returns a new KeepChars object, based on string keep"""
        if not case_sens:
            low = keep.lower()
            up = keep.upper()
            keep = low + up

        keep = keep.encode("utf-8")
        self._strip_table = dict([(c, None) for c in self.allchars if c not in keep])

    def __call__(self, s):
        """f(s) -> s, translates using self.allchars and self.delchars"""
        if s is None:
            raise TypeError
        if isinstance(s, bytes):
            s = s.decode("utf8")
        s = str(s)
        return s.translate(self._strip_table)


def first_index_in_set(seq, items):
    """Returns index of first occurrence of any of items in seq, or None."""
    for i, s in enumerate(seq):
        if s in items:
            return i
