#!/usr/bin/env python

"""Unit tests for Range, Span and Point classes.

NOTE: Map not currently tested.
"""
from cogent.util.unit_test import TestCase, main
from cogent.core.location import Span, Range, Point, RangeFromString

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class SpanTests(TestCase):
    """Tests of the Span object."""
    def setUp(self):
        """Define some standard Spans"""
        self.empty = Span(0, 0)
        self.full = Span(35, 30)    #will convert to (30, 35) internally
        self.overlapping = Span(32, 36)
        self.inside = Span(31, 32)
        self.before = Span(25, 30)
        self.after = Span(35, 40)
        self.reverse = Span(30, 35, Reverse=True)
        self.spans_zero = Span(-5, 5)
        
    def test_init(self):
        """Span object should init with Start, End, and Length"""
        s = Span(0)
        self.assertEqual(s.Start, 0)
        self.assertEqual(s.End, 1)
        self.assertEqual(s.Reverse, False)
        #to get an empty interval, must specify start and end explicitly
        t = Span(0, 0)
        self.assertEqual(t.Start, 0)
        self.assertEqual(t.End, 0)
        self.assertEqual(t.Reverse, False)
        #should be able to specify direction also
        u = Span(5, 15, Reverse=True)
        self.assertEqual(u.Start, 5)
        self.assertEqual(u.End, 15)
        self.assertEqual(u.Reverse, True)
        #should be able to init from another span
        v = Span(u)
        self.assertEqual(v.Start, 5)
        self.assertEqual(v.End, 15)
        self.assertEqual(v.Reverse, True)

    def test_contains(self):
        """Span object contains its start but not its end"""
        self.assertNotContains(self.empty, 0)
        self.assertContains(self.full, 30)
        self.assertContains(self.full, 34)
        self.assertNotContains(self.full, 35)
        self.assertContains(self.full, self.inside)
        self.assertNotContains(self.full, self.overlapping)
        self.assertContains(self.spans_zero, 0)
        self.assertContains(self.spans_zero, -5)
        self.assertNotContains(self.spans_zero, 5)
        
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

    def test_reverse(self):
        """Span reverse should change direction"""
        self.assertFalse(self.empty.Reverse)
        self.empty.reverse()
        self.assertTrue(self.empty.Reverse)
        self.empty.reverse()
        self.assertFalse(self.empty.Reverse)
        self.assertTrue(self.reverse.Reverse)
        self.reverse.reverse()
        self.assertFalse(self.reverse.Reverse)

    def test_iter(self):
        """Span iter should loop through (integer) contents"""
        self.assertEqual(list(iter(self.empty)), [])
        self.assertEqual(list(iter(self.full)), [30,31,32,33,34])
        self.assertEqual(list(iter(self.spans_zero)),[-5,-4,-3,-2,-1,0,1,2,3,4])
        self.assertEqual(list(iter(self.inside)), [31])
        self.assertEqual(list(self.reverse), [34,33,32,31,30])

    def test_str(self):
        """Span str should print start, stop, reverse"""
        self.assertEqual(str(self.empty), '(0,0,False)')
        self.assertEqual(str(self.full), '(30,35,False)')
        self.assertEqual(str(self.reverse), '(30,35,True)')

    def test_len(self):
        """Span len should return difference between start and end"""
        self.assertEqual(len(self.empty), 0)
        self.assertEqual(len(self.full), 5)
        self.assertEqual(len(self.inside),1)
        self.assertEqual(len(self.spans_zero), 10)

    def test_cmp(self):
        """Span cmp should support sort by 1st/2nd index and direction"""
        s, e, f, r, i, o = self.spans_zero, self.empty, self.full, \
            self.reverse, self.inside, self.overlapping

        n = Span(30, 36)

        expected_order = [s, e, f, r, n, i, o]
        first = expected_order[:]
        first.sort()
        second = [r, o, f, s, e, i, n]
        second.sort()
        for i, j in zip(first, second):
            self.assertSameObj(i, j)

        for i, j in zip(first, expected_order):
            self.assertSameObj(i, j)

    def test_startsBefore(self):
        """Span startsBefore should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertTrue(e.startsBefore(f))
        self.assertFalse(f.startsBefore(e))
        self.assertTrue(e.startsBefore(1))
        self.assertTrue(e.startsBefore(1000))
        self.assertFalse(e.startsBefore(0))
        self.assertFalse(e.startsBefore(-1))
        self.assertFalse(f.startsBefore(30))
        self.assertTrue(f.startsBefore(31))
        self.assertTrue(f.startsBefore(1000))
        
    def test_startsAfter(self):
        """Span startsAfter should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertFalse(e.startsAfter(f))
        self.assertTrue(f.startsAfter(e))
        self.assertFalse(e.startsAfter(1))
        self.assertFalse(e.startsAfter(1000))
        self.assertFalse(e.startsAfter(0))
        self.assertTrue(e.startsAfter(-1))
        self.assertTrue(f.startsAfter(29))
        self.assertFalse(f.startsAfter(30))
        self.assertFalse(f.startsAfter(31))
        self.assertFalse(f.startsAfter(1000))

    def test_startsAt(self):
        """Span startsAt should return True if input matches"""
        e, f = self.empty, self.full
        s = Span(30, 1000)
        self.assertTrue(e.startsAt(0))
        self.assertTrue(f.startsAt(30))
        self.assertTrue(s.startsAt(30))
        self.assertTrue(f.startsAt(s))
        self.assertTrue(s.startsAt(f))
        self.assertFalse(e.startsAt(f))
        self.assertFalse(e.startsAt(-1))
        self.assertFalse(e.startsAt(1))
        self.assertFalse(f.startsAt(29))

    def test_startsInside(self):
        """Span startsInside should return True if input starts inside span"""
        e, f, i, o = self.empty, self.full, self.inside, self.overlapping
        self.assertFalse(e.startsInside(0))
        self.assertFalse(f.startsInside(30))
        self.assertFalse(e.startsInside(f))
        self.assertTrue(i.startsInside(f))
        self.assertFalse(f.startsInside(i))
        self.assertTrue(o.startsInside(f))
        self.assertFalse(o.endsInside(i))
        
    def test_endsBefore(self):
        """Span endsBefore should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertTrue(e.endsBefore(f))
        self.assertFalse(f.endsBefore(e))
        self.assertTrue(e.endsBefore(1))
        self.assertTrue(e.endsBefore(1000))
        self.assertFalse(e.endsBefore(0))
        self.assertFalse(e.endsBefore(-1))
        self.assertFalse(f.endsBefore(30))
        self.assertFalse(f.endsBefore(31))
        self.assertTrue(f.endsBefore(1000))
        
    def test_endsAfter(self):
        """Span endsAfter should match hand-calculated results"""
        e, f = self.empty, self.full
        self.assertFalse(e.endsAfter(f))
        self.assertTrue(f.endsAfter(e))
        self.assertFalse(e.endsAfter(1))
        self.assertFalse(e.endsAfter(1000))
        self.assertFalse(e.endsAfter(0))
        self.assertTrue(e.endsAfter(-1))
        self.assertTrue(f.endsAfter(29))
        self.assertTrue(f.endsAfter(30))
        self.assertTrue(f.endsAfter(34))
        self.assertFalse(f.endsAfter(35))
        self.assertFalse(f.endsAfter(1000))

    def test_endsAt(self):
        """Span endsAt should return True if input matches"""
        e, f = self.empty, self.full
        s = Span(30, 1000)
        t = Span(-100, 35)
        self.assertTrue(e.endsAt(0))
        self.assertTrue(f.endsAt(35))
        self.assertTrue(s.endsAt(1000))
        self.assertFalse(f.endsAt(s))
        self.assertFalse(s.endsAt(f))
        self.assertTrue(f.endsAt(t))
        self.assertTrue(t.endsAt(f))

    def test_endsInside(self):
        """Span endsInside should return True if input ends inside span"""
        e, f, i, o = self.empty, self.full, self.inside, self.overlapping
        self.assertFalse(e.endsInside(0))
        self.assertFalse(f.endsInside(30))
        self.assertFalse(f.endsInside(34))
        self.assertFalse(f.endsInside(35))
        self.assertFalse(e.endsInside(f))
        self.assertTrue(i.endsInside(f))
        self.assertFalse(f.endsInside(i))
        self.assertFalse(o.endsInside(f))
        self.assertFalse(o.endsInside(i))
        self.assertTrue(e.endsInside(Span(-1,1)))
        self.assertTrue(e.endsInside(Span(0,1)))
        self.assertFalse(e.endsInside(Span(-1,0)))

class RangeInterfaceTests(object): #SpanTests):
    """A single-element Range should behave like the corresponding Span."""
    def setUp(self):
        """Define some standard Spans"""
        self.empty = Range(Span(0, 0))
        self.full = Range(Span(30, 35))
        self.overlapping = Range(Span(32, 36))
        self.inside = Range(Span(31, 32))
        self.before = Range(Span(25, 30))
        self.after = Range(Span(35, 40))
        self.reverse = Range(Span(30, 35, Reverse=True))
        self.spans_zero = Range(Span(-5, 5))

    def test_str(self):
        """Range str should print start, stop, reverse for each Span"""
        #note that the Range adds an extra level of parens, since it can
        #contain more than one Span.
        self.assertEqual(str(self.empty), '((0,0,False))')
        self.assertEqual(str(self.full), '((30,35,False))')
        self.assertEqual(str(self.reverse), '((30,35,True))')

 
class RangeTests(TestCase):
    """Tests of the Range object."""
    def setUp(self):
        """Set up a few standard ranges."""
        self.one = Range(Span(0,100))
        self.two = Range([Span(3,5), Span(8, 11)])
        self.three = Range([Span(6,7), Span(15, 17), Span(30, 35)])
        self.overlapping = Range([Span(6, 10), Span(7,3)])
        self.single = Range(0)
        self.singles = Range([3, 11])
        self.twocopy = Range(self.two)
        self.twothree = Range([self.two, self.three])
        self.empty = Range([Span(6,6), Span(8,8)])

    def test_init(self):
        """Range init from Spans, numbers, or Ranges should work OK."""
        #single span
        self.assertEqual(self.one, Span(0,100))
        #list of spans
        self.assertEqual(self.two.Spans, [Span(3,5), Span(8,11)])
        #another range
        self.assertEqual(self.two, self.twocopy)
        #list of ranges
        self.assertEqual(self.twothree.Spans, [Span(3,5), Span(8,11),
            Span(6,7), Span(15,17), Span(30,35)])
        #list of numbers
        self.assertEqual(self.singles.Spans, [Span(3,4), Span(11,12)])
        #single number
        self.assertEqual(self.single.Spans, [Span(0,1)])
        #nothing
        self.assertEqual(Range().Spans, [])
        
    def test_str(self):
        """Range str should print nested with parens"""
        self.assertEqual(str(self.one), '((0,100,False))')
        self.assertEqual(str(self.twothree), 
'((3,5,False),(8,11,False),(6,7,False),(15,17,False),(30,35,False))')
        self.assertEqual(str(self.single), '((0,1,False))')

    def test_len(self):
        """Range len should sum span lengths"""
        self.assertEqual(len(self.one), 100)
        self.assertEqual(len(self.single), 1)
        self.assertEqual(len(self.empty), 0)
        self.assertEqual(len(self.three), 8)

    def test_cmp(self):
        """Ranges should compare equal if they have the same spans"""
        self.assertEqual(self.twothree, Range([Span(3,5), Span(8, 11),
            Span(6,7), Span(15, 17), Span(30, 35)]))
        self.assertEqual(Range(), Range())

    def test_start_end(self):
        """Range Start and End should behave as expected"""
        self.assertEqual(self.one.Start, 0)
        self.assertEqual(self.one.End, 100)
        self.assertEqual(self.overlapping.Start, 3)
        self.assertEqual(self.overlapping.End, 10)
        self.assertEqual(self.three.Start, 6)
        self.assertEqual(self.three.End, 35)

    def test_reverse(self):
        """Range reverse method should reverse each span"""
        for s in self.overlapping.Spans:
            self.assertFalse(s.Reverse)
        self.overlapping.reverse()
        for s in self.overlapping.Spans:
            self.assertTrue(s.Reverse)
        self.overlapping.Spans.append(Span(0, 100))
        self.overlapping.reverse()
        for s in self.overlapping.Spans[0:1]:
            self.assertFalse(s.Reverse)
        self.assertTrue(self.overlapping.Spans[-1].Reverse)

    def test_Reverse(self):
        """Range Reverse property should return True if any span reversed"""
        self.assertFalse(self.one.Reverse)
        self.one.reverse()
        self.assertTrue(self.one.Reverse)
        self.assertFalse(self.two.Reverse)
        self.two.Spans.append(Span(0,100,Reverse=True))
        self.assertTrue(self.two.Reverse)
        self.two.reverse()
        self.assertTrue(self.two.Reverse)

    def test_contains(self):
        """Range contains an item if any span contains it"""
        self.assertContains(self.one, 50)
        self.assertContains(self.one, 0)
        self.assertContains(self.one, 99)
        self.assertNotContains(self.one, 100)
        self.assertContains(self.three, 6)
        self.assertNotContains(self.three, 7)
        self.assertNotContains(self.three, 8)
        self.assertNotContains(self.three, 14)
        self.assertContains(self.three, 15)
        self.assertNotContains(self.three, 29)
        self.assertContains(self.three, 30)
        self.assertContains(self.three, 34)
        self.assertNotContains(self.three, 35)
        self.assertNotContains(self.three, 40)
        #should work if a span is added
        self.three.Spans.append(40)
        self.assertContains(self.three, 40)
        #should work for spans
        self.assertContains(self.three, Span(31, 33))
        self.assertNotContains(self.three, Span(31, 37))
        #span contains itself
        self.assertContains(self.two, self.twocopy)
        #should work for ranges
        self.assertContains(self.three, Range([6, Span(15,16), Span(30,33)]))
        #should work for copy, except when extra piece added
        threecopy = Range(self.three)
        self.assertContains(self.three, threecopy)
        threecopy.Spans.append(1000)
        self.assertNotContains(self.three, threecopy)
        self.three.Spans.append(Span(950, 1050))
        self.assertContains(self.three, threecopy)
        self.assertNotContains(threecopy, self.three)

    def test_overlaps(self):
        """Range overlaps should return true if any component overlapping"""
        self.assertTrue(self.two.overlaps(self.one))
        self.assertTrue(self.one.overlaps(self.two))
        self.assertTrue(self.three.overlaps(self.one))
        #two and three are interleaved but not overlapping
        self.assertFalse(self.two.overlaps(self.three))
        self.assertFalse(self.three.overlaps(self.two))
        self.assertTrue(self.one.overlaps(self.empty))
        self.assertTrue(self.empty.overlaps(self.one))
        self.assertTrue(self.singles.overlaps(self.two))

    def test_overlapsExtent(self):
        """Range overlapsExtent should return true for interleaved ranges"""
        self.assertTrue(self.two.overlapsExtent(self.three))
        self.assertTrue(self.three.overlapsExtent(self.two))
        self.assertFalse(self.single.overlapsExtent(self.two))
        self.assertFalse(self.single.overlapsExtent(self.three))
        self.assertTrue(self.one.overlapsExtent(self.three))

    def test_sort(self):
        """Range sort should sort component spans"""
        one = self.one
        one.sort()
        self.assertEqual(one.Spans, [Span(100,0)])
        one.Spans.append(Span(-20,-10))
        self.assertEqual(one.Spans, [Span(0,100),Span(-20,-10)])
        one.sort()
        self.assertEqual(one.Spans, [Span(-20,-10),Span(0,100)])
        one.Spans.append(Span(-20, -10, Reverse=True))
        self.assertEqual(one.Spans, [Span(-20,-10),Span(0,100),
            Span(-20,-10,Reverse=True)])
        one.sort()
        self.assertEqual(one.Spans, [Span(-20,-10),Span(-20,-10,Reverse=True),
            Span(0,100)])

    def test_iter(self):
        """Range iter should iterate through each span in turn"""
        self.assertEqual(list(iter(self.two)), [3,4,8,9,10])
        self.two.Spans.insert(1, Span(103, 101, Reverse=True))
        self.assertEqual(list(iter(self.two)), [3,4,102,101,8,9,10])

    def test_Extent(self):
        """Range extent should span limits of range"""
        self.assertEqual(self.one.Extent, Span(0,100))
        self.assertEqual(self.three.Extent, Span(6,35))
        self.assertEqual(self.singles.Extent, Span(3, 12))
        self.assertEqual(self.single.Extent, Span(0,1))
        self.three.Spans.append(Span(100, 105, Reverse=True))
        self.assertEqual(self.three.Extent, Span(6,105))
        self.three.Spans.append(Span(-100, -1000))
        self.assertEqual(self.three.Extent, Span(-1000,105))

    def test_simplify(self):
        """Range reduce should group overlapping ranges"""
        #consolidate should have no effect when no overlap
        r = self.two
        r.simplify()
        self.assertEqual(r.Spans, [Span(3,5), Span(8,11)])
        #should consolidate an overlap of the same direction
        r.Spans.append(Span(-1, 4))
        r.simplify()
        self.assertEqual(r.Spans, [Span(-1,5), Span(8,11)])
        #should also consolidate _adjacent_ spans of the same direction
        r.Spans.append(Span(11,14))
        r.simplify()
        self.assertEqual(r.Spans, [Span(-1,5), Span(8,14)])
        #bridge should cause consolidations
        s = Range(r)
        s.Spans.append(Span(5,8))
        s.simplify()
        self.assertEqual(s.Spans, [Span(-1,14)])
        #ditto for bridge that overlaps everything
        s = Range(r)
        s.Spans.append(Span(-100, 100))
        s.simplify()
        self.assertEqual(s.Spans, [Span(-100,100)])
        #however, can't consolidate span in other orientation
        s = Range(r)
        s.Spans.append(Span(-100, 100, Reverse=True))
        self.assertEqual(s.Spans, [Span(-1,5), Span(8,14), \
            Span(-100,100,Reverse=True)])

class RangeFromStringTests(TestCase):
    """Tests of the RangeFromString factory function."""
    def test_init(self):
        self.assertEqual(RangeFromString(''), Range())
        self.assertEqual(RangeFromString('  3  , 4\t, ,, 10  ,'),
            Range([3,4,10]))
        self.assertEqual(RangeFromString('3,4-10,1-5'),
            Range([Span(3), Span(4,10), Span(1,5)]))

#run the following if invoked from command-line
if __name__ == "__main__":
    main()


