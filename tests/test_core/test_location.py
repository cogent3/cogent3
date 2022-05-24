#!/usr/bin/env python

"""Unit tests for Range, Span and Point classes.
"""
from unittest import TestCase, main

from cogent3 import DNA
from cogent3.core.location import (
    Map,
    Range,
    RangeFromString,
    Span,
    gap_coords_to_map,
)


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class SpanTests(TestCase):
    """Tests of the Span object."""

    def setUp(self):
        """Define some standard Spans"""
        self.empty = Span(0, 0)
        self.full = Span(35, 30)  # will convert to (30, 35) internally
        self.overlapping = Span(32, 36)
        self.inside = Span(31, 32)
        self.before = Span(25, 30)
        self.after = Span(35, 40)
        self.reverse = Span(30, 35, reverse=True)
        self.spans_zero = Span(-5, 5)

    def test_init(self):
        """Span object should init with start, end, and Length"""
        s = Span(0)
        self.assertEqual(s.start, 0)
        self.assertEqual(s.end, 1)
        self.assertEqual(s.reverse, False)
        # to get an empty interval, must specify start and end explicitly
        t = Span(0, 0)
        self.assertEqual(t.start, 0)
        self.assertEqual(t.end, 0)
        self.assertEqual(t.reverse, False)
        # should be able to specify direction also
        u = Span(5, 15, reverse=True)
        self.assertEqual(u.start, 5)
        self.assertEqual(u.end, 15)
        self.assertEqual(u.reverse, True)
        # should be able to init from another span
        v = Span(u)
        self.assertEqual(v.start, 5)
        self.assertEqual(v.end, 15)
        self.assertEqual(v.reverse, True)

    def test_contains(self):
        """Span object contains its start but not its end"""
        self.assertNotIn(0, self.empty)
        self.assertIn(30, self.full)
        self.assertIn(34, self.full)
        self.assertNotIn(35, self.full)
        self.assertIn(self.inside, self.full)
        self.assertNotIn(self.overlapping, self.full)
        self.assertIn(0, self.spans_zero)
        self.assertIn(-5, self.spans_zero)
        self.assertNotIn(5, self.spans_zero)

    def test_overlaps(self):
        """Span objects should be able to overlap points or spans"""
        self.assertTrue(self.full.overlaps(self.overlapping))
        self.assertFalse(self.full.overlaps(self.before))
        self.assertFalse(self.before.overlaps(self.overlapping))
        self.assertFalse(self.full.overlaps(self.after))
        self.assertFalse(self.after.overlaps(self.before))
        self.assertTrue(self.full.overlaps(self.inside))
        self.assertTrue(self.spans_zero.overlaps(self.empty))
        self.assertTrue(self.empty.overlaps(self.spans_zero))

    def test_reverses(self):
        """Span.reverses should change direction"""
        self.assertFalse(self.empty.reverse)
        self.empty.reverses()
        self.assertTrue(self.empty.reverse)
        self.empty.reverses()
        self.assertFalse(self.empty.reverse)
        self.assertTrue(self.reverse.reverse)
        self.reverse.reverses()
        self.assertFalse(self.reverse.reverse)

    def test_iter(self):
        """Span iter should loop through (integer) contents"""
        self.assertEqual(list(iter(self.empty)), [])
        self.assertEqual(list(iter(self.full)), [30, 31, 32, 33, 34])
        self.assertEqual(
            list(iter(self.spans_zero)), [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4]
        )
        self.assertEqual(list(iter(self.inside)), [31])
        self.assertEqual(list(self.reverse), [34, 33, 32, 31, 30])

    def test_str(self):
        """Span str should print start, stop, reverse"""
        self.assertEqual(str(self.empty), "(0,0,False)")
        self.assertEqual(str(self.full), "(30,35,False)")
        self.assertEqual(str(self.reverse), "(30,35,True)")

    def test_len(self):
        """Span len should return difference between start and end"""
        self.assertEqual(len(self.empty), 0)
        self.assertEqual(len(self.full), 5)
        self.assertEqual(len(self.inside), 1)
        self.assertEqual(len(self.spans_zero), 10)

    def test_cmp(self):
        """Span cmp should support sort by 1st/2nd index and direction"""
        s, e, f, r, i, o = (
            self.spans_zero,
            self.empty,
            self.full,
            self.reverse,
            self.inside,
            self.overlapping,
        )

        n = Span(30, 36)

        expected_order = [s, e, f, r, n, i, o]
        first = expected_order[:]
        first.sort()
        second = [r, o, f, s, e, i, n]
        second.sort()
        for i, j in zip(first, second):
            self.assertIs(i, j)

        for i, j in zip(first, expected_order):
            self.assertIs(i, j)

    def test_sort(self):
        """Span should support sort by 1st/2nd index and direction"""
        s, e, f, r, i, o = (
            self.spans_zero,
            self.empty,
            self.full,
            self.reverse,
            self.inside,
            self.overlapping,
        )

        expected_order = [s, e]
        first = expected_order[:]
        first.sort()

        for i, j in zip(first, expected_order):
            self.assertIs(i, j)

    def test_starts_before(self):
        """Span starts_before should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertTrue(e.starts_before(f))
        self.assertFalse(f.starts_before(e))
        self.assertTrue(e.starts_before(1))
        self.assertTrue(e.starts_before(1000))
        self.assertFalse(e.starts_before(0))
        self.assertFalse(e.starts_before(-1))
        self.assertFalse(f.starts_before(30))
        self.assertTrue(f.starts_before(31))
        self.assertTrue(f.starts_before(1000))

    def test_starts_after(self):
        """Span starts_after should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertFalse(e.starts_after(f))
        self.assertTrue(f.starts_after(e))
        self.assertFalse(e.starts_after(1))
        self.assertFalse(e.starts_after(1000))
        self.assertFalse(e.starts_after(0))
        self.assertTrue(e.starts_after(-1))
        self.assertTrue(f.starts_after(29))
        self.assertFalse(f.starts_after(30))
        self.assertFalse(f.starts_after(31))
        self.assertFalse(f.starts_after(1000))

    def test_startsAt(self):
        """Span startsAt should return True if input matches"""
        e, f = self.empty, self.full
        s = Span(30, 1000)
        self.assertTrue(e.starts_at(0))
        self.assertTrue(f.starts_at(30))
        self.assertTrue(s.starts_at(30))
        self.assertTrue(f.starts_at(s))
        self.assertTrue(s.starts_at(f))
        self.assertFalse(e.starts_at(f))
        self.assertFalse(e.starts_at(-1))
        self.assertFalse(e.starts_at(1))
        self.assertFalse(f.starts_at(29))

    def test_startsInside(self):
        """Span startsInside should return True if input starts inside span"""
        e, f, i, o = self.empty, self.full, self.inside, self.overlapping
        self.assertFalse(e.starts_inside(0))
        self.assertFalse(f.starts_inside(30))
        self.assertFalse(e.starts_inside(f))
        self.assertTrue(i.starts_inside(f))
        self.assertFalse(f.starts_inside(i))
        self.assertTrue(o.starts_inside(f))
        self.assertFalse(o.ends_inside(i))

    def test_endsBefore(self):
        """Span endsBefore should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertTrue(e.ends_before(f))
        self.assertFalse(f.ends_before(e))
        self.assertTrue(e.ends_before(1))
        self.assertTrue(e.ends_before(1000))
        self.assertFalse(e.ends_before(0))
        self.assertFalse(e.ends_before(-1))
        self.assertFalse(f.ends_before(30))
        self.assertFalse(f.ends_before(31))
        self.assertTrue(f.ends_before(1000))

    def test_endsAfter(self):
        """Span endsAfter should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertFalse(e.ends_after(f))
        self.assertTrue(f.ends_after(e))
        self.assertFalse(e.ends_after(1))
        self.assertFalse(e.ends_after(1000))
        self.assertFalse(e.ends_after(0))
        self.assertTrue(e.ends_after(-1))
        self.assertTrue(f.ends_after(29))
        self.assertTrue(f.ends_after(30))
        self.assertTrue(f.ends_after(34))
        self.assertFalse(f.ends_after(35))
        self.assertFalse(f.ends_after(1000))

    def test_endsAt(self):
        """Span endsAt should return True if input matches"""
        e, f = self.empty, self.full
        s = Span(30, 1000)
        t = Span(-100, 35)
        self.assertTrue(e.ends_at(0))
        self.assertTrue(f.ends_at(35))
        self.assertTrue(s.ends_at(1000))
        self.assertFalse(f.ends_at(s))
        self.assertFalse(s.ends_at(f))
        self.assertTrue(f.ends_at(t))
        self.assertTrue(t.ends_at(f))

    def test_ends_inside(self):
        """Span ends_inside should return True if input ends inside span"""
        e, f, i, o = self.empty, self.full, self.inside, self.overlapping
        self.assertFalse(e.ends_inside(0))
        self.assertFalse(f.ends_inside(30))
        self.assertFalse(f.ends_inside(34))
        self.assertFalse(f.ends_inside(35))
        self.assertFalse(e.ends_inside(f))
        self.assertTrue(i.ends_inside(f))
        self.assertFalse(f.ends_inside(i))
        self.assertFalse(o.ends_inside(f))
        self.assertFalse(o.ends_inside(i))
        self.assertTrue(e.ends_inside(Span(-1, 1)))
        self.assertTrue(e.ends_inside(Span(0, 1)))
        self.assertFalse(e.ends_inside(Span(-1, 0)))


class RangeInterfaceTests(object):  # SpanTests):
    """A single-element Range should behave like the corresponding Span."""

    def setUp(self):
        """Define some standard Spans"""
        self.empty = Range(Span(0, 0))
        self.full = Range(Span(30, 35))
        self.overlapping = Range(Span(32, 36))
        self.inside = Range(Span(31, 32))
        self.before = Range(Span(25, 30))
        self.after = Range(Span(35, 40))
        self.reverse = Range(Span(30, 35, reverse=True))
        self.spans_zero = Range(Span(-5, 5))

    def test_str(self):
        """Range str should print start, stop, reverse for each Span"""
        # note that the Range adds an extra level of parens, since it can
        # contain more than one Span.
        self.assertEqual(str(self.empty), "((0,0,False))")
        self.assertEqual(str(self.full), "((30,35,False))")
        self.assertEqual(str(self.reverse), "((30,35,True))")


class RangeTests(TestCase):
    """Tests of the Range object."""

    def setUp(self):
        """Set up a few standard ranges."""
        self.one = Range(Span(0, 100))
        self.two = Range([Span(3, 5), Span(8, 11)])
        self.three = Range([Span(6, 7), Span(15, 17), Span(30, 35)])
        self.overlapping = Range([Span(6, 10), Span(7, 3)])
        self.single = Range(0)
        self.singles = Range([3, 11])
        self.twocopy = Range(self.two)
        self.twothree = Range([self.two, self.three])
        self.empty = Range([Span(6, 6), Span(8, 8)])

    def test_init(self):
        """Range init from Spans, numbers, or Ranges should work OK."""
        # single span
        self.assertEqual(self.one, Span(0, 100))
        # list of spans
        self.assertEqual(self.two.spans, [Span(3, 5), Span(8, 11)])
        # another range
        self.assertEqual(self.two, self.twocopy)
        # list of ranges
        self.assertEqual(
            self.twothree.spans,
            [Span(3, 5), Span(8, 11), Span(6, 7), Span(15, 17), Span(30, 35)],
        )
        # list of numbers
        self.assertEqual(self.singles.spans, [Span(3, 4), Span(11, 12)])
        # single number
        self.assertEqual(self.single.spans, [Span(0, 1)])
        # nothing
        self.assertEqual(Range().spans, [])

    def test_str(self):
        """Range str should print nested with parens"""
        self.assertEqual(str(self.one), "((0,100,False))")
        self.assertEqual(
            str(self.twothree),
            "((3,5,False),(8,11,False),(6,7,False),(15,17,False),(30,35,False))",
        )
        self.assertEqual(str(self.single), "((0,1,False))")

    def test_len(self):
        """Range len should sum span lengths"""
        self.assertEqual(len(self.one), 100)
        self.assertEqual(len(self.single), 1)
        self.assertEqual(len(self.empty), 0)
        self.assertEqual(len(self.three), 8)

    def test_cmp(self):
        """Ranges should compare equal if they have the same spans"""
        self.assertEqual(
            self.twothree,
            Range([Span(3, 5), Span(8, 11), Span(6, 7), Span(15, 17), Span(30, 35)]),
        )
        self.assertEqual(Range(), Range())

    def test_start_end(self):
        """Range start and end should behave as expected"""
        self.assertEqual(self.one.start, 0)
        self.assertEqual(self.one.end, 100)
        self.assertEqual(self.overlapping.start, 3)
        self.assertEqual(self.overlapping.end, 10)
        self.assertEqual(self.three.start, 6)
        self.assertEqual(self.three.end, 35)

    def test_reverse(self):
        """Range reverse method should reverse each span"""
        for s in self.overlapping.spans:
            self.assertFalse(s.reverse)
        self.overlapping.reverses()
        for s in self.overlapping.spans:
            self.assertTrue(s.reverse)
        self.overlapping.spans.append(Span(0, 100))
        self.overlapping.reverses()
        for s in self.overlapping.spans[0:1]:
            self.assertFalse(s.reverse)
        self.assertTrue(self.overlapping.spans[-1].reverse)

    def test_Reverse(self):
        """Range reverse property should return True if any span reversed"""
        self.assertFalse(self.one.reverse)
        self.one.reverses()
        self.assertTrue(self.one.reverse)
        self.assertFalse(self.two.reverse)
        self.two.spans.append(Span(0, 100, reverse=True))
        self.assertTrue(self.two.reverse)
        self.two.reverses()
        self.assertTrue(self.two.reverse)

    def test_contains(self):
        """Range contains an item if any span contains it"""
        self.assertIn(50, self.one)
        self.assertIn(0, self.one)
        self.assertIn(99, self.one)
        self.assertNotIn(100, self.one)
        self.assertIn(6, self.three)
        self.assertNotIn(7, self.three)
        self.assertNotIn(8, self.three)
        self.assertNotIn(14, self.three)
        self.assertIn(15, self.three)
        self.assertNotIn(29, self.three)
        self.assertIn(30, self.three)
        self.assertIn(34, self.three)
        self.assertNotIn(35, self.three)
        self.assertNotIn(40, self.three)
        # should work if a span is added
        self.three.spans.append(40)
        self.assertIn(40, self.three)
        # should work for spans
        self.assertIn(Span(31, 33), self.three)
        self.assertNotIn(Span(31, 37), self.three)
        # span contains itself
        self.assertIn(self.twocopy, self.two)
        # should work for ranges
        self.assertIn(Range([6, Span(15, 16), Span(30, 33)]), self.three)
        # should work for copy, except when extra piece added
        threecopy = Range(self.three)
        self.assertIn(threecopy, self.three)
        threecopy.spans.append(1000)
        self.assertNotIn(threecopy, self.three)
        self.three.spans.append(Span(950, 1050))
        self.assertIn(threecopy, self.three)
        self.assertNotIn(self.three, threecopy)

    def test_overlaps(self):
        """Range overlaps should return true if any component overlapping"""
        self.assertTrue(self.two.overlaps(self.one))
        self.assertTrue(self.one.overlaps(self.two))
        self.assertTrue(self.three.overlaps(self.one))
        # two and three are interleaved but not overlapping
        self.assertFalse(self.two.overlaps(self.three))
        self.assertFalse(self.three.overlaps(self.two))
        self.assertTrue(self.one.overlaps(self.empty))
        self.assertTrue(self.empty.overlaps(self.one))
        self.assertTrue(self.singles.overlaps(self.two))

    def test_overlaps_extent(self):
        """Range overlaps_extent should return true for interleaved ranges"""
        self.assertTrue(self.two.overlaps_extent(self.three))
        self.assertTrue(self.three.overlaps_extent(self.two))
        self.assertFalse(self.single.overlaps_extent(self.two))
        self.assertFalse(self.single.overlaps_extent(self.three))
        self.assertTrue(self.one.overlaps_extent(self.three))

    def test_sort(self):
        """Range sort should sort component spans"""
        one = self.one
        one.sort()
        self.assertEqual(one.spans, [Span(100, 0)])
        one.spans.append(Span(-20, -10))
        self.assertEqual(one.spans, [Span(0, 100), Span(-20, -10)])
        one.sort()
        self.assertEqual(one.spans, [Span(-20, -10), Span(0, 100)])
        one.spans.append(Span(-20, -10, reverse=True))
        self.assertEqual(
            one.spans, [Span(-20, -10), Span(0, 100), Span(-20, -10, reverse=True)]
        )
        one.sort()
        self.assertEqual(
            one.spans, [Span(-20, -10), Span(-20, -10, reverse=True), Span(0, 100)]
        )

    def test_iter(self):
        """Range iter should iterate through each span in turn"""
        self.assertEqual(list(iter(self.two)), [3, 4, 8, 9, 10])
        self.two.spans.insert(1, Span(103, 101, reverse=True))
        self.assertEqual(list(iter(self.two)), [3, 4, 102, 101, 8, 9, 10])

    def test_extent(self):
        """Range extent should span limits of range"""
        self.assertEqual(self.one.extent, Span(0, 100))
        self.assertEqual(self.three.extent, Span(6, 35))
        self.assertEqual(self.singles.extent, Span(3, 12))
        self.assertEqual(self.single.extent, Span(0, 1))
        self.three.spans.append(Span(100, 105, reverse=True))
        self.assertEqual(self.three.extent, Span(6, 105))
        self.three.spans.append(Span(-100, -1000))
        self.assertEqual(self.three.extent, Span(-1000, 105))

    def test_simplify(self):
        """Range reduce should group overlapping ranges"""
        # consolidate should have no effect when no overlap
        r = self.two
        r.simplify()
        self.assertEqual(r.spans, [Span(3, 5), Span(8, 11)])
        # should consolidate an overlap of the same direction
        r.spans.append(Span(-1, 4))
        r.simplify()
        self.assertEqual(r.spans, [Span(-1, 5), Span(8, 11)])
        # should also consolidate _adjacent_ spans of the same direction
        r.spans.append(Span(11, 14))
        r.simplify()
        self.assertEqual(r.spans, [Span(-1, 5), Span(8, 14)])
        # bridge should cause consolidations
        s = Range(r)
        s.spans.append(Span(5, 8))
        s.simplify()
        self.assertEqual(s.spans, [Span(-1, 14)])
        # ditto for bridge that overlaps everything
        s = Range(r)
        s.spans.append(Span(-100, 100))
        s.simplify()
        self.assertEqual(s.spans, [Span(-100, 100)])
        # however, can't consolidate span in other orientation
        s = Range(r)
        s.spans.append(Span(-100, 100, reverse=True))
        self.assertEqual(
            s.spans, [Span(-1, 5), Span(8, 14), Span(-100, 100, reverse=True)]
        )


class RangeFromStringTests(TestCase):
    """Tests of the RangeFromString factory function."""

    def test_init(self):
        self.assertEqual(RangeFromString(""), Range())
        self.assertEqual(RangeFromString("  3  , 4\t, ,, 10  ,"), Range([3, 4, 10]))
        self.assertEqual(
            RangeFromString("3,4-10,1-5"), Range([Span(3), Span(4, 10), Span(1, 5)])
        )


class MapTests(TestCase):
    """tests of the Map class"""

    def test_get_coords(self):
        """get_coordinates should return raw coordinates matching input"""
        spans = [(0, 9), (20, 32)]
        map = Map(spans, parent_length=100)
        coords = map.get_coordinates()
        self.assertEqual(coords, spans)

        # should work for reversed Maps too
        spans = [(32, 20), (9, 0)]
        map = Map(spans, parent_length=100)
        coords = map.get_coordinates()
        self.assertEqual(coords, spans)

    def test_get_gap_coords(self):
        """returns gap start and lengths"""
        m, seq = DNA.make_seq("-AC--GT-TTA--").parse_out_gaps()
        got = m.get_gap_coordinates()
        self.assertEqual(dict(got), {0: 1, 2: 2, 4: 1, 7: 2})

    def test_gap_coords_to_map(self):
        """construct a Map from coordinates of gap alone"""
        m, seq = DNA.make_seq("-AC--GT-TTA--").parse_out_gaps()
        gap_coords = {0: 1, 2: 2, 4: 1, 7: 2}
        seqlen = 70
        got = gap_coords_to_map(gap_coords, seqlen)
        self.assertEqual(len(got), seqlen + sum(gap_coords.values()))

        gap_coords = {5: 2, 17: 3, 10: 2}
        seqlen = 20
        got = gap_coords_to_map(gap_coords, seqlen)
        self.assertEqual(len(got), sum(gap_coords.values()) + seqlen)

        # roundtrip from Map.get_gap_coordinates()
        self.assertEqual(dict(got.get_gap_coordinates()), gap_coords)

        # and no gaps
        m, seq = DNA.make_seq("ACGTTTA").parse_out_gaps()
        got = gap_coords_to_map({}, len(seq))
        self.assertEqual(len(got), len(m))
        self.assertEqual(got.get_coordinates(), m.get_coordinates())

        # and gaps outside sequence
        with self.assertRaises(ValueError):
            got = gap_coords_to_map({20: 1}, len(seq))


# run the following if invoked from command-line
if __name__ == "__main__":
    main()
