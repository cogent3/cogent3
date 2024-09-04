"""Unit tests for Span classes."""

from itertools import combinations
from unittest import TestCase

import numpy
import pytest
from cogent3 import DNA, make_seq
from cogent3.core import new_moltype
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
        m, _ = DNA.make_seq(seq="-AC--GT-TTA--").parse_out_gaps()
        got = m.get_gap_coordinates()
        self.assertEqual(dict(got), {0: 1, 2: 2, 4: 1, 7: 2})


def test_map_plus_position():
    # seq is 9 long
    # plus coords  012345678
    # +slice         **
    # plus seq     AAACCCTGG

    orig = FeatureMap.from_locations(locations=[(0, 9)], parent_length=9)
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
    im = IndelMap.from_spans(spans=[LostSpan(3)], parent_length=0)
    assert not im.useful
    assert not im.complete
    assert len(im) == 3


def test_map_nucleic_reversed():
    expect = [(0, 9)]
    # seq is 9 long
    # plus coords  012345678
    # plus seq     AAACCCTGG

    orig = FeatureMap.from_locations(locations=[(0, 9)], parent_length=9)
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
    #                       01   2 345
    seq = DNA.make_seq(seq="AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = cls.from_spans(spans=tuple(m.spans), parent_length=m.parent_length)
    got = m.get_coordinates()
    assert got == [(0, 2), (2, 3), (3, 6)]


@pytest.mark.parametrize("cls", (FeatureMap, IndelMap))
def test_gap_coordinate(cls):
    seq = DNA.make_seq(seq="AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = cls.from_spans(spans=tuple(m.spans), parent_length=m.parent_length)
    got = m.get_gap_coordinates()
    assert [tuple(c) for c in got] == [(2, 3), (3, 1), (6, 2)]


def test_gaps():
    # returns spans corresponding to position on "aligned" seq of gaps
    #                       000000000011
    #                       012345678901
    seq = DNA.make_seq(seq="AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = FeatureMap.from_spans(spans=tuple(m.spans), parent_length=m.parent_length)
    got = [(g.start, g.end) for g in tuple(m.gaps().spans)]
    assert got == [(2, 5), (6, 7), (10, 12)]


@pytest.mark.parametrize("cls", (IndelMap, FeatureMap))
def test_nongap(cls):
    # returns spans corresponding to position on "aligned" seq of nongaps
    #                       000000000011
    #                       012345678901
    seq = DNA.make_seq(seq="AC---G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = cls.from_spans(spans=tuple(m.spans), parent_length=m.parent_length)

    got = [(s.start, s.end) for s in m.nongap()]
    assert got == [(0, 2), (5, 6), (7, 10)]


@pytest.mark.parametrize("cls", (IndelMap, FeatureMap))
def test_nongap_startswith(cls):
    # returns spans corresponding to position on "aligned" seq of nongaps
    #                   012345678
    seq = DNA.make_seq(seq="--G-TAA--")
    m, _ = seq.parse_out_gaps()
    m = cls.from_spans(spans=tuple(m.spans), parent_length=m.parent_length)

    got = [(s.start, s.end) for s in m.nongap()]
    assert got == [(2, 3), (4, 7)]


@pytest.mark.parametrize("cls", (IndelMap, FeatureMap))
def test_nongap_not_endswith(cls):
    # returns spans corresponding to position on "aligned" seq of nongaps
    #                   0123456
    seq = DNA.make_seq(seq="--G-TAA")
    m, _ = seq.parse_out_gaps()
    m = cls.from_spans(spans=tuple(m.spans), parent_length=m.parent_length)

    got = [(s.start, s.end) for s in m.nongap()]
    assert got == [(2, 3), (4, 7)]


def test_spans_gen():
    # returns spans corresponding to position on "aligned" seq of nongaps
    #                       000000000011
    #                       012345678901
    seq = DNA.make_seq(seq="AC---G-TAA--")
    expect = [Span(0, 2), LostSpan(3), Span(2, 3), LostSpan(1), Span(3, 6), LostSpan(2)]
    m, _ = seq.parse_out_gaps()
    gap_data = numpy.array([(2, 3), (3, 1), (6, 2)])
    pos, lengths = gap_data.T
    im = IndelMap(gap_pos=pos, gap_lengths=lengths, parent_length=m.parent_length)
    got = list(im.spans)
    assert got == expect


def test_spans_gap_start():
    # returns spans corresponding to position on "aligned" seq of nongaps
    #                       000000000011
    #                       012345678901
    seq = DNA.make_seq(seq="---TAA")
    expect = [LostSpan(3), Span(0, 3)]
    im, _ = seq.parse_out_gaps()
    got = list(im.spans)
    assert got == expect


def test_spans_gen_nogaps():
    # returns spans corresponding to position on "aligned" seq of nongaps
    #                       000000000011
    #                       012345678901
    seq = DNA.make_seq(seq="ACGTAA")
    m, _ = seq.parse_out_gaps()
    spans = list(m.spans)
    assert len(spans) == 1
    assert len(spans[0]) == len(seq)


def test_round_trip_rich_dict():
    seq = DNA.make_seq(seq="AC---G-TAA--")
    im, _ = seq.parse_out_gaps()
    rd = im.to_rich_dict()
    got = IndelMap.from_rich_dict(rd)
    assert im is not got
    assert got.to_rich_dict() == im.to_rich_dict()


def test_terminal_unknown():
    # seq  idx   01   2 345  6
    #           -AC---G-TAA--
    # aligned seq length is 13
    gap_data = numpy.array([[0, 1], [2, 3], [3, 1], [6, 2]])
    gap_pos, gap_lengths = gap_data.T

    m = IndelMap(
        gap_pos=gap_pos.copy(), gap_lengths=gap_lengths.copy(), parent_length=6
    )
    # not unknown, by default
    m_spans = tuple(m.spans)
    assert m_spans[0].lost and not isinstance(m_spans[0], TerminalPadding)
    # use the constructor arg
    m = IndelMap(
        gap_pos=gap_pos.copy(),
        gap_lengths=gap_lengths.copy(),
        parent_length=6,
        termini_unknown=True,
    )
    m_spans = tuple(m.spans)
    assert isinstance(m_spans[0], TerminalPadding)
    assert isinstance(m_spans[-1], TerminalPadding)
    assert m_spans[2].lost and not isinstance(m_spans[1], TerminalPadding)
    assert m_spans[4].lost and not isinstance(m_spans[2], TerminalPadding)

    # use the method
    m = IndelMap(
        gap_pos=gap_pos.copy(), gap_lengths=gap_lengths.copy(), parent_length=6
    ).with_termini_unknown()
    m_spans = tuple(m.spans)
    assert isinstance(m_spans[0], TerminalPadding)
    assert isinstance(m_spans[-1], TerminalPadding)
    # middle gap is not terminal, so...
    assert not isinstance(m_spans[2], TerminalPadding)

    # no gaps, no effect
    # use the constructor arg
    empty = numpy.array([], dtype=int)
    m = IndelMap(
        gap_pos=empty.copy(),
        gap_lengths=empty.copy(),
        parent_length=6,
        termini_unknown=True,
    )
    m_spans = tuple(m.spans)
    assert len(m_spans) == 1 and not m_spans[0].lost
    # use the method
    m = IndelMap(
        gap_pos=empty.copy(), gap_lengths=empty.copy(), parent_length=6
    ).with_termini_unknown()
    m_spans = tuple(m.spans)
    assert len(m_spans) == 1 and not m_spans[0].lost


def test_featuremap_inverse():
    m = FeatureMap.from_locations(locations=[(0, 2), (4, 6)], parent_length=6)
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


def test_indelmap_from_aligned_segments():
    locations = [(0, 2), (4, 6)]
    im = IndelMap.from_aligned_segments(locations=locations, aligned_length=6)
    assert len(im) == 6
    expected_length = sum(e - s for s, e in locations)
    assert im.parent_length == expected_length
    im_spans = tuple(im.spans)
    assert im_spans[1].lost and len(im_spans[1]) == 2


@pytest.mark.parametrize("swap", (False, True))
def test_indelmap_inconsistent_input(swap):
    a = numpy.array([0, 1], dtype=int)
    b = numpy.array([3], dtype=int)
    a, b = (b, a) if swap else (a, b)
    with pytest.raises(ValueError):
        IndelMap(gap_pos=a, gap_lengths=b, parent_length=10)

    with pytest.raises(ValueError):
        IndelMap(gap_pos=a, cum_gap_lengths=b, parent_length=10)


def test_indelmap_from_aligned_segments2():
    locations = [(0, 5), (7, 12), (14, 21), (24, 27)]
    im = IndelMap.from_aligned_segments(locations=locations, aligned_length=27)
    expected_length = sum(e - s for s, e in locations)
    assert im.parent_length == expected_length


def test_indelmap_merge():
    im1 = IndelMap(
        gap_pos=numpy.array([], dtype=int),
        gap_lengths=numpy.array([], dtype=int),
        parent_length=4,
    )
    im2 = IndelMap(
        gap_pos=numpy.array([0], dtype=int),
        gap_lengths=numpy.array([2], dtype=int),
        parent_length=4,
    )

    assert len(im1) == 4
    assert im1.parent_length == 4
    inv = im1.merge_maps(im2)
    assert inv.parent_length == 4
    assert len(inv) == 6

    fmap1 = im1.to_feature_map()
    fmap2 = im2.to_feature_map()
    got = fmap2.inverse()[fmap1.inverse()].inverse()
    assert inv.get_gap_coordinates() == [
        list(coords) for coords in got.get_gap_coordinates()
    ]


def test_indelmap_merge_parent_length():
    im1 = IndelMap(
        gap_pos=numpy.array([], dtype=int),
        gap_lengths=numpy.array([], dtype=int),
        parent_length=4,
    )
    im2 = IndelMap(
        gap_pos=numpy.array([0], dtype=int),
        gap_lengths=numpy.array([2], dtype=int),
        parent_length=4,
    )
    # providing a value for parent_length overrides standard
    ov = im1.merge_maps(im2, parent_length=20)
    assert ov.parent_length != im2.parent_length
    assert ov.parent_length == 20


def test_map_indexed():
    m = FeatureMap.from_locations(locations=[(0, 2), (4, 6)], parent_length=6).inverse()
    indexed = m[2]
    assert len(indexed) == 1


def test_compare_map_indexed():
    from cogent3.core.alignment import Aligned

    raw_seq = "--AC-GTAA--"
    im, seq = DNA.make_seq(seq=raw_seq).parse_out_gaps()
    ia = Aligned(im, seq)
    length = len(raw_seq)
    got = [str(ia[i]) for i in range(length)]
    expect = list("--AC-GTAA--")
    assert got == expect


@pytest.mark.parametrize("slice_it", (True, False))
def test_feature_map_zeroed(slice_it):
    spans = [LostSpan(2), Span(2, 4), LostSpan(2), Span(4, 8), LostSpan(2)]
    kwargs = dict(spans=spans, parent_length=6)

    fm = FeatureMap(**kwargs)
    if slice_it:
        fm = fm[3:6]

    fm_zeroed = fm.zeroed()
    assert fm.start > 0
    assert fm_zeroed.start == 0


def test_indelmap_to_feature_map():
    spans = [LostSpan(2), Span(2, 4), LostSpan(2), Span(4, 8), LostSpan(2)]
    kwargs = dict(spans=spans, parent_length=8)
    im = IndelMap.from_spans(**kwargs)
    mm = im.to_feature_map()
    assert mm.get_coordinates() == im.get_coordinates()


@pytest.mark.parametrize("raw", ("--AC--GGGG--", "A-A-A", "-A-AA----A"))
def test_indelmap_nucleic_reversed(raw):
    from cogent3.core.alignment import Aligned

    plus = DNA.make_seq(seq=raw)
    minus = plus.rc()
    plus_imap, _ = DNA.make_seq(seq=raw).parse_out_gaps()
    minus_imap, minus_seq = minus.parse_out_gaps()
    got = plus_imap.nucleic_reversed()
    assert got.get_coordinates() == minus_imap.get_coordinates()
    assert (got.gap_pos == minus_imap.gap_pos).all()
    assert (got.cum_gap_lengths == minus_imap.cum_gap_lengths).all()
    assert got.parent_length == minus_imap.parent_length
    assert str(Aligned(got, minus_seq)) == str(minus)


def test_get_coords():
    """get_coordinates should return raw coordinates matching input"""
    spans = [(0, 9), (20, 32)]
    fmap = FeatureMap.from_locations(locations=spans, parent_length=100)
    coords = fmap.get_coordinates()
    assert coords == spans


def test_indel_map_get_coords():
    """get_coordinates should return raw coordinates matching input"""
    imap = IndelMap(
        gap_pos=numpy.array([9]), gap_lengths=numpy.array([10]), parent_length=20
    )
    locations = [(0, 9), (9, 20)]
    coords = imap.get_coordinates()
    assert coords == locations


def test_indel_map_get_gap_coords():
    gap_data = numpy.array([(2, 3), (3, 1), (6, 2)])
    gap_pos, lengths = gap_data[:, 0], gap_data[:, 1]
    m = IndelMap(gap_pos=gap_pos, gap_lengths=lengths, parent_length=10)
    got = m.get_gap_coordinates()
    assert got == gap_data.tolist()


@pytest.mark.parametrize(
    "coords", ([(0, 3), (7, 11)], [(0, 3)], [(2, 4), (6, 10)], [(4, 6)])
)
def test_indelmap_joined_segments(coords):
    raw = "--AC--GGGG--"
    expect, _ = DNA.make_seq(seq="".join(raw[s:e] for s, e in coords)).parse_out_gaps()
    imap, _ = DNA.make_seq(seq=raw).parse_out_gaps()
    got = imap.joined_segments(coords)
    assert got.gap_pos.tolist() == expect.gap_pos.tolist()
    assert got.cum_gap_lengths.tolist() == expect.cum_gap_lengths.tolist()


def test_indelmap_slice_terminating():
    raw = "-CB-A--"
    start, end = 4, 6
    expect, _ = DNA.make_seq(seq=raw[start:end]).parse_out_gaps()
    imap, _ = DNA.make_seq(seq=raw).parse_out_gaps()
    got = imap[start:end]
    assert got.gap_pos.tolist() == expect.gap_pos.tolist()
    assert got.cum_gap_lengths.tolist() == expect.cum_gap_lengths.tolist()


def test_indelmap_slice_zero():
    raw = "-CB-A--"
    start, end = 0, 0
    expect, s = DNA.make_seq(seq=raw[start:end]).parse_out_gaps()
    imap, _ = DNA.make_seq(seq=raw).parse_out_gaps()
    got = imap[start:end]
    assert got.parent_length == len(s)
    assert got.gap_pos.tolist() == expect.gap_pos.tolist()
    assert got.cum_gap_lengths.tolist() == expect.cum_gap_lengths.tolist()


def test_indelmap_invalid_slice_range():
    # If we mirror python slicing, an invalid slice should return an empty map
    imap = IndelMap(
        gap_pos=numpy.array([10], dtype=int),
        gap_lengths=numpy.array([2], dtype=int),
        parent_length=10,
    )

    expect = numpy.array([], dtype=int)

    got = imap[-100]
    assert (got.gap_pos == expect).all()
    assert (got.cum_gap_lengths == expect).all()

    got = imap[-100:]
    assert (got.gap_pos == expect).all()
    assert (got.cum_gap_lengths == expect).all()

    got = imap[:-99]
    assert (got.gap_pos == expect).all()
    assert (got.cum_gap_lengths == expect).all()


def test_indelmap_get_indices_errors():
    imap = IndelMap(
        gap_pos=numpy.array([10], dtype=int),
        gap_lengths=numpy.array([2], dtype=int),
        parent_length=10,
    )
    with pytest.raises(IndexError):
        imap.get_align_index(-1000)


def test_indelmap_slice_innner_gap():
    start, end = 4, 6
    raw = "--AC--GGGG--"
    expect, _ = DNA.make_seq(seq=raw[start:end]).parse_out_gaps()
    imap, _ = DNA.make_seq(seq=raw).parse_out_gaps()
    imap = imap[start:end]
    assert imap.gap_pos.tolist() == expect.gap_pos.tolist()


def test_indelmap_slice_cum_length():
    start, end = 7, 11
    raw = "--AC--GGGG--"
    expect, _ = DNA.make_seq(seq=raw[start:end]).parse_out_gaps()
    imap, _ = DNA.make_seq(raw).parse_out_gaps()
    imap = imap[start:end]
    assert imap.gap_pos.tolist() == expect.gap_pos.tolist()
    assert imap.cum_gap_lengths.tolist() == expect.cum_gap_lengths.tolist()


def test_get_coords_invalid_order():
    """get_coordinates should return raw coordinates matching input"""

    # should work for reversed Maps too
    spans = [(32, 20), (9, 0)]
    with pytest.raises(ValueError):
        FeatureMap.from_locations(locations=spans, parent_length=100)


def test_indelmap_slice_gap_into_seq():
    pos, lengths = numpy.array([[3, 1], [7, 1]]).T
    gaps = IndelMap(gap_pos=pos, gap_lengths=lengths, parent_length=7)
    # example seq
    # 012345678
    # 012 3456
    # ACA-CGAC-
    # slicing from aligned coords 3-5, gives '-C'
    got = gaps[3:5]
    assert got.gap_pos.tolist() == [0]
    assert got.cum_gap_lengths.tolist() == [1]


def test_indelmap_slice_one_gap():
    pos, lengths = numpy.array([[3, 6]]).T
    gaps = IndelMap(gap_pos=pos, gap_lengths=lengths, parent_length=24)
    # example seq
    #           1
    # 01234567890
    # 012      3456789
    # GTG------GTAGAAGTTCCAAATAATGAA
    # slicing from aligned coords 3-5, gives 'TG------G'
    # in sliced seq, gap at 2
    got = gaps[1:10]
    assert got.gap_pos.tolist() == [2]
    assert got.cum_gap_lengths.tolist() == [6]


@pytest.mark.parametrize(
    "data",
    (
        "AB---CD--EF",
        "---ABCD--EF",
        "ABCD---EF--",
        "-----ABCDEF",
        "ABCDEF-----",
        "-ABCDEF----",
        "-A-B-C-D-EF",
        "A-B-C-D-EF-",
    ),
)
@pytest.mark.parametrize("index", range(4, 6))  # the ungapped sequence is 6 long
def test_gapped_convert_seq2aln(data, index):
    # converting a sequence index to alignment index
    ungapped = data.replace("-", "")
    seq = make_seq(data, moltype="text")
    gaps, _ = seq.parse_out_gaps()
    idx = gaps.get_align_index(index)
    assert data[idx] == ungapped[index]


@pytest.mark.parametrize(
    "data",
    (
        "AB---CD--EF",
        "---ABCD--EF",
        "ABCD---EF--",
        "-----ABCDEF",
        "ABCDEF-----",
        "-ABCDEF----",
        "-A-B-C-D-EF",
        "A-B-C-D-EF-",
    ),
)
@pytest.mark.parametrize(
    "start,end", list(combinations(range(6), 2))
)  # the ungapped sequence is 6 long
def test_indelmap_align_index_slice_stop(data, start, end):
    # converting a sequence index to alignment index
    ungapped = data.replace("-", "")
    seq = make_seq(data, moltype="text")
    gaps, _ = seq.parse_out_gaps()
    stop = gaps.get_align_index(end, slice_stop=True)
    assert data[stop - 1] == ungapped[end - 1]


@pytest.mark.parametrize(
    "data",
    (
        "AB---CD--EF",
        "---ABCD--EF",
        "ABCD---EF--",
        "-----ABCDEF",
        "ABCDEF-----",
        "-ABCDEF----",
        "-A-B-C-D-EF",
        "A-B-C-D-EF-",
    ),
)
@pytest.mark.parametrize("index", range(6))  # the ungapped sequence is 6 long
def test_gapped_convert_seq2aln2seq(data, index):
    # round tripping seq to alignment to seq indices
    gaps, _ = make_seq(data, moltype="text").parse_out_gaps()
    align_index = gaps.get_align_index(index)
    got = gaps.get_seq_index(align_index)
    assert got == index


@pytest.mark.parametrize("expect", range(10))
def test_indelmap_get_seq_index_negative(expect):
    parent_length = 10
    gap_pos = [0]
    gap_lengths = [3]
    imap = IndelMap(
        gap_pos=numpy.array(gap_pos),
        gap_lengths=numpy.array(gap_lengths),
        parent_length=parent_length,
    )
    neg_index = expect - parent_length
    got = imap.get_seq_index(neg_index)
    assert got == expect


@pytest.mark.parametrize("expect", range(10))
def test_indelmap_get_align_index_negative(expect):
    parent_length = 10
    gap_pos = [0]
    gap_lengths = [3]
    gap_length = gap_lengths[0]
    imap = IndelMap(
        gap_pos=numpy.array(gap_pos),
        gap_lengths=numpy.array(gap_lengths),
        parent_length=parent_length,
    )
    neg_index = expect + gap_length - len(imap)
    got = imap.get_seq_index(neg_index)
    assert got == expect


@pytest.mark.parametrize(
    "data",
    (
        "AB--CDE-FG",
        "--ABC-DEFG",
        "AB--CDE-FG--",
        "ABCDE--FG---",
        "-----ABCDEFG",
        "-A-B-C-D-E-F-G-",
    ),
)
@pytest.mark.parametrize("seq_index", range(7))
def test_gapped_convert_aln2seq_nongap_char(data, seq_index):
    # test alignment indexes when aligned position is NOT a gap
    ungapped = "ABCDEFG"
    align_index = data.find(ungapped[seq_index])
    gaps, _ = make_seq(data, moltype="text").parse_out_gaps()
    idx = gaps.get_seq_index(align_index)
    assert idx == seq_index


def _find_nth_gap_index(data: str, n: int) -> int:
    num = -1
    for i, c in enumerate(data):
        if c == "-":
            num += 1
        if num == n:
            return i
    raise ValueError(f"{data=}, {n=}")


def _get_expected_seqindex(data: str, align_index: int) -> int:
    # compute the expected seqindex
    refseq = data.replace("-", "")
    got = data[align_index:].lstrip("-")
    return refseq.find(got[0]) if got else len(refseq)


@pytest.mark.parametrize(
    "data",
    (
        "AB-----CDE-F--G",
        "----ABC-DEFG---",
        "AB--CDE-FG-----",
        "ABCDE--FG------",
        "--------ABCDEFG",
        "-A-B-C-D-E-F-G-",
    ),
)
@pytest.mark.parametrize("gap_number", range(8))
def test_gapped_convert_aln2seq_gapchar(data, gap_number):
    # test alignment indexes when aligned position IS a gap
    # in this case we expect the position of the next non-gap
    # to be the result
    # find the alignment index corresponding to the
    align_index = _find_nth_gap_index(data, gap_number)
    assert data[align_index] == "-", (data, gap_number)
    # find nearest non-gap
    seq_index = _get_expected_seqindex(data, align_index)
    gaps, _ = make_seq(data, moltype="text").parse_out_gaps()
    idx = gaps.get_seq_index(align_index)
    assert idx == seq_index


def test_gapped_convert_aln2seq_invalid():
    gaps, _ = make_seq(seq="AC--GTA-TG", moltype="dna").parse_out_gaps()
    with pytest.raises(IndexError):
        # absolute value of negative indices must be < seq length
        gaps.get_seq_index(-100)


@pytest.mark.parametrize(
    "aslice",
    (
        slice(3, 7),
        slice(20, None),
    ),
)
def test_no_gaps_in_slice(aslice):
    # aligned length is 25
    seq_length = 20
    gap_length = 5
    pos, lengths = numpy.array([[10, gap_length]], dtype=numpy.int32).T
    gp = IndelMap(gap_pos=pos, gap_lengths=lengths, parent_length=seq_length)
    got = gp[aslice]
    assert not got.num_gaps
    start = aslice.start or 0
    stop = aslice.stop or (seq_length + gap_length)
    assert len(got) == stop - start


def test_len_gapped():
    seq_length = 20
    gap_length = 5
    pos, lengths = numpy.array([[10, gap_length]], dtype=numpy.int32).T
    gp = IndelMap(gap_pos=pos, gap_lengths=lengths, parent_length=seq_length)
    assert len(gp) == (seq_length + gap_length)


def test_all_gaps_in_slice():
    # slicing GapPositions
    # sample seq 1
    data = "AC--GTA-TG"
    gp, _ = make_seq(data, moltype="dna").parse_out_gaps()
    sl = slice(1, 9)

    got = gp[sl]
    expect_gaps, _ = make_seq(data[sl], moltype="dna").parse_out_gaps()
    assert got.get_gap_coordinates() == expect_gaps.get_gap_coordinates()
    assert got.parent_length == 5


@pytest.mark.parametrize(
    "data",
    (
        "----GTA-TG",
        "AC--GTA---",
        "AC--GTA-TG",
        "A-C-G-T-A-",
        "-A-C-G-T-A",
        "ACGTAACGTA",
        "----------",
    ),
)
@pytest.mark.parametrize(
    "aslice",
    [slice(i, j) for i, j in combinations(range(10), 2)],
)
def test_variant_slices(data, aslice):
    gaps, _ = make_seq(data, moltype="dna").parse_out_gaps()
    got = gaps[aslice]

    expect_gaps, expect_seq = make_seq(data[aslice], moltype="dna").parse_out_gaps()
    assert got.parent_length == len(expect_seq)
    assert (got.gap_pos == expect_gaps.gap_pos).all()
    assert (got.cum_gap_lengths == expect_gaps.cum_gap_lengths).all()


def test_indelmap_joined():
    pos = numpy.array([0, 1])
    cum_len = numpy.array([1, 5])
    imap = IndelMap(gap_pos=pos, cum_gap_lengths=cum_len, parent_length=5)
    fmap = FeatureMap.from_locations(locations=[(0, 1), (2, 5)], parent_length=10)
    coords = fmap.get_coordinates()
    got = imap.joined_segments(coords)
    assert got.num_gaps == 1
    assert got.gap_pos[0] == 0
    assert got.cum_gap_lengths[0] == (1 + 5 - 2)


def test_indel_map_sliced_with_negative():
    imap = IndelMap(
        gap_pos=numpy.array([0]), cum_gap_lengths=numpy.array([1]), parent_length=14
    )
    got = imap[:-3]
    assert got.parent_length == 14 - 3


def test_indelmap_roundtrip_json():
    from cogent3.util.deserialise import deserialise_object

    imap = IndelMap(
        gap_pos=numpy.array([0, 9]),
        cum_gap_lengths=numpy.array([1, 3]),
        parent_length=14,
    )
    got = deserialise_object(imap.to_json())
    assert (got.gap_pos == imap.gap_pos).all()
    assert (got.cum_gap_lengths == imap.cum_gap_lengths).all()
    assert got.parent_length == imap.parent_length


def test_featuremap_roundtrip_json():
    from cogent3.util.deserialise import deserialise_object

    fmap = FeatureMap.from_locations(
        locations=[[0, 9], [20, 23]],
        parent_length=140,
    )
    got = deserialise_object(fmap.to_json())
    coords = fmap.get_coordinates()
    assert coords == [(0, 9), (20, 23)]
    assert got.parent_length == fmap.parent_length == 140


@pytest.mark.parametrize(
    "error_type,locations", ((ValueError, ((-2, 2),)), (RuntimeError, ((20, 25),)))
)
def test_invalid_locations(error_type, locations):
    with pytest.raises(error_type):
        FeatureMap.from_locations(locations=locations, parent_length=10)


def test_gap_coords_to_map():
    """construct a Map from coordinates of gap alone"""
    gap_coords = {0: 1, 2: 2, 4: 1, 7: 2}
    seqlen = 70
    got = gap_coords_to_map(gap_coords, seqlen)
    assert len(got) == seqlen + sum(gap_coords.values())

    gap_coords = {5: 2, 17: 3, 10: 2}
    seqlen = 20
    got = gap_coords_to_map(gap_coords, seqlen)
    assert len(got) == sum(gap_coords.values()) + seqlen

    # roundtrip from Map.get_gap_coordinates()
    assert dict(got.get_gap_coordinates()) == gap_coords

    # and no gaps
    m, seq = DNA.make_seq(seq="ACGTTTA").parse_out_gaps()
    got = gap_coords_to_map({}, len(seq))
    assert len(got) == len(m)
    assert got.get_coordinates() == m.get_coordinates()

    # and gaps outside sequence
    with pytest.raises(ValueError):
        gap_coords_to_map({20: 1}, len(seq))


@pytest.mark.parametrize(
    "seq, expected",
    [
        ("ACGTACGT", []),
        ("----ACGTACGT----", [(0, 4), (12, 16)]),
        ("ACGT----ACGT", [4, 8]),
        ("----", [[0, 4]]),
        ("----ACGT", [[0, 4]]),
        ("ACGTACGT----", [(8, 12)]),
    ],
)
def test_get_gap_align_coordinates1(seq, expected):
    expected = numpy.array(expected, dtype=int)
    if not len(expected):
        expected = expected.reshape((0, 2))
    im, _ = make_seq(seq).parse_out_gaps()
    result = im.get_gap_align_coordinates()
    assert (result == expected).all(), f"{expected=}, {result=}"


@pytest.mark.parametrize(
    "seq, expected",
    [
        ("", []),
        ("A", []),
        ("-", [[0, 1]]),
    ],
)
def test_get_gap_align_coordinates_edge_cases(seq, expected):
    expected = numpy.array(expected, dtype=int)
    if not len(expected):
        expected = expected.reshape((0, 2))

    im, _ = make_seq(seq).parse_out_gaps()
    result = im.get_gap_align_coordinates()
    assert (result == expected).all(), f"{expected=}, {result=}"


def test_featuremap_add():
    spans_a = [LostSpan(2), Span(2, 4)]
    spans_b = [LostSpan(2), Span(4, 8), LostSpan(2)]
    kwargs = dict(parent_length=6)

    fm_a = FeatureMap(spans=spans_a, **kwargs)
    fm_b = FeatureMap(spans=spans_b, **kwargs)
    fm_ab = fm_a + fm_b
    assert list(fm_ab.spans) == (spans_a + spans_b)


def test_featuremap_mul():
    spans = [LostSpan(2), Span(2, 4)]
    fm = FeatureMap(spans=spans, parent_length=6)
    fm_3 = fm * 3
    assert list(fm_3.spans) == [sp * 3 for sp in spans]
    assert fm_3.parent_length == 6 * 3


def test_featuremap_div():
    spans = [LostSpan(3), Span(3, 6)]
    fm_3 = FeatureMap(spans=spans, parent_length=6)
    fm_1 = fm_3 / 3
    assert list(fm_1.spans) == [sp / 3 for sp in spans]
    assert fm_1.parent_length == 6 / 3


def test_indelmap_make_seq_feature_map():
    #           1
    # 01234567890
    # AC--GTA-TAA
    # 01    234 567
    im = IndelMap(
        gap_pos=numpy.array([2, 5], dtype=int),
        gap_lengths=numpy.array([2, 1], dtype=int),
        parent_length=8,
    )
    orig_spans = [Span(1, 5)]
    align_map = FeatureMap(spans=orig_spans, parent_length=11)
    spans = [Span(1, 3)]
    expect = FeatureMap(spans=spans, parent_length=8)
    got = im.make_seq_feature_map(align_map)
    assert got.get_coordinates() == expect.get_coordinates()
    assert got.parent_length == expect.parent_length

    # ignoring lost spans
    align_map = FeatureMap(spans=orig_spans + [LostSpan(4)], parent_length=11)
    got = im.make_seq_feature_map(align_map)
    assert got.get_coordinates() == expect.get_coordinates()
    assert got.parent_length == expect.parent_length


@pytest.mark.parametrize("start", range(10))
@pytest.mark.parametrize("stop", range(11))
@pytest.mark.parametrize("step", range(1, 5))
@pytest.mark.parametrize(
    "data",
    [
        "TCAGTCAGTC",
        "AAAAA-----",
        "-----AAAAA",
        "--AA--AA--",
        "AA--AA--AA",
        "A-A-A-A-A-",
        "---TTTT---",
        "C--------C",
        "----------",
    ],
)
def test_indelmap_positive_step_variant_slices(start, stop, step, data):
    imap, _ = new_moltype.DNA.make_seq(seq=data).parse_out_gaps()
    got = imap[start:stop:step]
    expect, _ = new_moltype.DNA.make_seq(seq=data[start:stop:step]).parse_out_gaps()
    assert (got.gap_pos == expect.gap_pos).all()
    assert (got.cum_gap_lengths == expect.cum_gap_lengths).all()


@pytest.mark.parametrize("start", range(11))
@pytest.mark.parametrize("stop", range(10))
@pytest.mark.parametrize("step", [-1, -2, -3])
@pytest.mark.parametrize(
    "data",
    [
        "TCAGTCAGTC",
        "AAAAA-----",
        "AAAAAAA--A",
        "-----AAAAA",
        "--AA--AA--",
        "AA--AA--AA",
        "A-A-A-A-A-",
        "---TTTT---",
        "C--------C",
        "----------",
    ],
)
def test_indelmap_negative_step_variant_slices(start, stop, step, data):
    imap, _ = new_moltype.DNA.make_seq(seq=data).parse_out_gaps()
    got = imap[start:stop:step]
    expect, _ = new_moltype.DNA.make_seq(seq=data[start:stop:step]).parse_out_gaps()
    assert (got.gap_pos == expect.gap_pos).all()
    assert (got.cum_gap_lengths == expect.cum_gap_lengths).all()
