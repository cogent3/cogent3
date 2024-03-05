"""Unit tests for Span classes.
"""
from unittest import TestCase

import pytest

from cogent3 import DNA
from cogent3.core.location import (
    FeatureMap,
    IndelMap,
    LostSpan,
    Span,
    TerminalPadding,
    gap_coords_to_map,
)


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


class MapTests(TestCase):
    """tests of the Map class"""

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


@pytest.mark.parametrize("cls", (IndelMap, FeatureMap))
def test_map_plus_position(cls):
    # seq is 9 long
    # plus coords  012345678
    # +slice         **
    # plus seq     AAACCCTGG

    # orig = Aligned(*DNA.make_seq("AAACCCTGG", name="a").parse_out_gaps())
    orig = cls(locations=[(0, 9)], parent_length=9)
    assert orig.absolute_position(2) == 2
    assert orig.absolute_position(6) == 6

    assert orig.relative_position(2) == 2
    assert orig.relative_position(6) == 6

    subseq = orig[2:4]
    assert subseq.absolute_position(0) == 2
    assert subseq.absolute_position(4) == 6

    assert subseq.relative_position(2) == 0
    assert subseq.relative_position(6) == 4

    # minus coords 012345678
    # rel coord      01234
    # -slice         *****
    # minus seq    CCAGGGTTT
    # plus coords  876543210
    rc = orig.nucleic_reversed()
    rcsubseq = rc[2:7]
    abs_coord = rcsubseq.absolute_position(0)


def test_indel_map_useful_complete():
    im = IndelMap(spans=[LostSpan(3)], parent_length=0)
    assert not im.useful
    assert not im.complete
    assert len(im) == im.length == 3


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_map_nucleic_reversed(cls):
    expect = [(0, 9)]
    # seq is 9 long
    # plus coords  012345678
    # +slice         **
    # plus seq     AAACCCTGG

    # orig = Aligned(*DNA.make_seq("AAACCCTGG", name="a").parse_out_gaps())
    orig = cls(locations=[(0, 9)], parent_length=9)
    # minus coords 012345678
    # rel coord      01234
    # -slice         *****
    # minus seq    CCAGGGTTT
    # plus coords  876543210
    rc = orig.nucleic_reversed()
    coords = rc.get_coordinates()
    assert coords == expect


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_coordinate(cls):
    # coordinates are for ungapped segments in underlying sequence
    #                   01   2 345
    seq = DNA.make_seq("AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = cls(spans=tuple(m.spans), parent_length=m.parent_length)
    got = m.get_coordinates()
    assert got == [(0, 2), (2, 3), (3, 6)]


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_gap_coordinate(cls):
    seq = DNA.make_seq("AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = cls(spans=tuple(m.spans), parent_length=m.parent_length)
    got = m.get_gap_coordinates()
    assert got == [(2, 3), (3, 1), (6, 2)]


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_gaps(cls):
    # returns spans corresponding to position on "aligned" seq of gaps
    #                   000000000011
    #                   012345678901
    seq = DNA.make_seq("AC---G-TAA--")
    m, s = seq.parse_out_gaps()
    m = cls(spans=tuple(m.spans), parent_length=m.parent_length)
    got = [(g.start, g.end) for g in tuple(m.gaps().spans)]
    assert got == [(2, 5), (6, 7), (10, 12)]


@pytest.mark.parametrize("cls", (IndelMap, FeatureMap))
def test_nongap(cls):
    # returns spans corresponding to position on "aligned" seq of nongaps
    #                   000000000011
    #                   012345678901
    seq = DNA.make_seq("AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = cls(spans=tuple(m.spans), parent_length=m.parent_length)

    got = [(g.start, g.end) for g in m.nongap().spans]
    assert got == [(0, 2), (5, 6), (7, 10)]


def test_round_trip_rich_dict():
    seq = DNA.make_seq("AC---G-TAA--")
    m, s = seq.parse_out_gaps()
    # reversed() reverses the order of spans, but keeps their coordinates
    # differs from nucleic reversed, which computes a new relative position
    im = IndelMap(spans=m.spans, parent_length=m.parent_length)
    rd = im.to_rich_dict()
    got = IndelMap.from_rich_dict(rd)
    assert im is not got
    assert got.to_rich_dict() == im.to_rich_dict()


def test_serialisable_attr():
    im = IndelMap(locations=[(0, 2)], parent_length=20)
    set_vals = {"locations": [(0, 2)], "parent_length": 20}
    got = {k: im._serialisable[k] for k in set_vals}
    assert got == set_vals


def test_terminal_unknown():
    # span idx          01 2  345  6
    seq = DNA.make_seq("-AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    # not unknown, by default
    m_spans = tuple(m.spans)
    assert m_spans[0].lost and not isinstance(m_spans[0], TerminalPadding)
    # use the constructor arg
    m = IndelMap(
        spans=tuple(m.spans), parent_length=m.parent_length, termini_unknown=True
    )
    m_spans = tuple(m.spans)
    assert isinstance(m_spans[0], TerminalPadding)
    assert isinstance(m_spans[-1], TerminalPadding)
    assert m_spans[2].lost and not isinstance(m_spans[1], TerminalPadding)
    assert m_spans[4].lost and not isinstance(m_spans[2], TerminalPadding)

    # use the method
    m, _ = seq.parse_out_gaps()
    m = IndelMap(
        spans=tuple(m.spans), parent_length=m.parent_length
    ).with_termini_unknown()
    m_spans = tuple(m.spans)
    assert isinstance(m_spans[0], TerminalPadding)
    assert isinstance(m_spans[-1], TerminalPadding)
    # middle gap is not terminal, so...
    assert not isinstance(m_spans[2], TerminalPadding)

    # no gaps, no effect
    seq = DNA.make_seq("ACGTAA")
    m, _ = seq.parse_out_gaps()
    # use the constructor arg
    m = IndelMap(
        spans=tuple(m.spans), parent_length=m.parent_length, termini_unknown=True
    )
    m_spans = tuple(m.spans)
    assert not isinstance(m_spans[0], TerminalPadding)
    # use the method
    m = IndelMap(
        spans=tuple(m.spans), parent_length=m.parent_length
    ).with_termini_unknown()
    m_spans = tuple(m.spans)
    assert not isinstance(m_spans[0], TerminalPadding)


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_map_inverse(cls):
    m = cls(locations=[(0, 2), (4, 6)], parent_length=6)
    assert len(m) == 4
    mi = m.inverse()
    assert len(mi) == 6
    mi_spans = tuple(mi.spans)
    assert mi_spans[1].lost and len(mi_spans[1]) == 2
    # invert the inversion, should give us back the original
    re_inv = mi.inverse()
    expect = m.to_rich_dict()
    got = re_inv.to_rich_dict()
    assert got == expect


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_map_offsets(cls):
    # offsets are absolute starts of spans
    #                              1
    #                   01 3  678  1
    seq = DNA.make_seq("-AC---G-TAA--")
    m, s = seq.parse_out_gaps()
    got = m.offsets
    assert got == [0, 1, 3, 6, 7, 8, 11]


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_map_indexed(cls):
    m = cls(locations=[(0, 2), (4, 6)], parent_length=6).inverse()
    indexed = m[2]
    assert len(indexed) == 1


def test_compare_map_indexed():
    from cogent3.core.alignment import Aligned

    seq = DNA.make_seq("--AC-GTAA--".replace("-", ""))
    spans = [LostSpan(2), Span(0, 2), LostSpan(2), Span(2, 6), LostSpan(2)]
    kwargs = dict(spans=spans, parent_length=len(seq))
    mm = FeatureMap(**kwargs)
    im = IndelMap(**kwargs)
    ma = Aligned(mm, seq)
    ia = Aligned(im, seq)
    first = ia[0]
    assert first == "-"
    length = len(ma)
    got = [str(ia[i]) for i in range(length)]
    expect = [str(ma[i]) for i in range(length)]
    assert got == expect


@pytest.mark.parametrize("slice_it", (True, False))
def test_indel_map_zeroed(slice_it):
    spans = [LostSpan(2), Span(2, 4), LostSpan(2), Span(4, 8), LostSpan(2)]
    kwargs = dict(spans=spans, parent_length=6)

    mm = FeatureMap(**kwargs)
    if slice_it:
        mm = mm[:6]
    mm_zeroed = mm.zeroed()

    im = IndelMap(**kwargs)
    if slice_it:
        im = im[:6]

    im_zeroed = im.zeroed()
    assert im_zeroed.get_coordinates() == mm_zeroed.get_coordinates()
    assert im_zeroed.parent_length == mm_zeroed.parent_length


def test_indelmap_to_feature_map():
    spans = [LostSpan(2), Span(2, 4), LostSpan(2), Span(4, 8), LostSpan(2)]
    kwargs = dict(spans=spans, parent_length=6)
    im = IndelMap(**kwargs)
    mm = im.to_feature_map()
    assert mm.get_coordinates() == im.get_coordinates()


def test_indelmap_nucleic_reversed():
    spans = [LostSpan(2), Span(2, 4), LostSpan(2), Span(4, 8), LostSpan(2)]
    kwargs = dict(spans=spans, parent_length=12)
    orig = IndelMap(**kwargs)
    rev = orig.nucleic_reversed()
    rev_spans = tuple(rev.spans)
    assert rev_spans[1].reverse == rev_spans[3].reverse == False
    old = orig.nucleic_reversed()
    old_spans = tuple(old.spans)
    assert rev.get_coordinates() == old.get_coordinates()


def test_indelmap_with_reverse_span():
    spans = [
        LostSpan(2),
        Span(8, 4, reverse=True),
        LostSpan(2),
        Span(4, 2, reverse=True),
        LostSpan(2),
    ]
    imap = IndelMap(spans=spans, parent_length=12)
    assert (imap.start, imap.end) == (8, 2)


def test_indelmap_no_gaps():
    imap = IndelMap(locations=(), parent_length=6)
    gaps = imap.gaps()
    assert not gaps


def test_get_coords():
    """get_coordinates should return raw coordinates matching input"""
    spans = [(0, 9), (20, 32)]
    map = FeatureMap(locations=spans, parent_length=100)
    coords = map.get_coordinates()
    assert coords == spans


def test_get_coords_invalid_order():
    """get_coordinates should return raw coordinates matching input"""

    # should work for reversed Maps too
    spans = [(32, 20), (9, 0)]
    with pytest.raises(ValueError):
        FeatureMap(locations=spans, parent_length=100)
