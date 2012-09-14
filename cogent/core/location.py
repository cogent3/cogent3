#!/usr/bin/env python
"""Alignments and Sequences are _Annotatables
_Annotatables hold a list of Maps.
Maps can be Features, Variables or AlignedSequences.
Maps have a list of Spans.

Also provides Range and Point classes for dealing with parts of sequences.

Span is a region with a start, an end, and a direction. Range is an ordered
collection of Spans (note: Range does _not_ support the list interface, but
you can always access Range.Spans directly). Map is like a Range but is
immutable and is able to be nested, i.e. Maps can be defined relative to
other Maps.

Implementation Notes

Span and Range behave much like Python's slices: a Span contains the element
after its Start but does not contain the element after its End. It may help to
think of the Span indices occurring _between_ the list elements:
    
    a b c d e
   | | | | | |
   0 1 2 3 4 5

...so that a Span whose Start is its End contains no elements (e.g. 2:2), and
a Span whose End is 2 more than its start contains 2 elements (e.g. 2:4 has c
and d), etc. Similarly, Span(0,2) does _not_ overlap Span(2,3), since the
former contains a and b while the latter contains c.

A Point is a Span whose Start and End refer to the same object, i.e. the same
position in the sequence. A Point occurs between elements in the sequence,
and so does not contain any elements itself.

WARNING: this differs from the way e.g. NCBI handles sequence indices, where
the sequence is 1-based, a single index is treated as containing one element,
the point 3 contains exactly one element, 3, rather than no elements, and a
range from 2:4 contains 2, 3 and 4, _not_ just 2 and 3.

"""
from cogent.util.misc import FunctionWrapper, ClassChecker, ConstrainedList, \
    iterable
from itertools import chain
from string import strip

from bisect import bisect_right, bisect_left
import copy

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Matthew Wakefield",
                    "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

def _norm_index(i, length, default):
    """For converting s[:3] to s[0:3], s[-1] to s[len(s)-1] and s[0:lots] to s[0:len(s)]"""
    if i is None:
        i = default
    elif i < 0:
        i += length
    return min(max(i,0),length)

def _norm_slice(index, length):
    """_norm_slice(slice(1, -2, 3), 10) -> (1,8,3)"""
    if isinstance(index, slice):
        start = _norm_index(index.start, length, 0)
        end = _norm_index(index.stop, length, length)
        return (start, end, index.step)
    else:
        start = index
        if start < 0: start += length
        if start >= length: raise IndexError(index)
        return (start, start+1, 1)

def as_map(slice, length):
    """Take anything that might be used as a subscript: Integer, Slice,
    or Map, and return a Map."""
    
    if isinstance(slice, (list, tuple)):
        spans = []
        for i in slice:
            spans.extend(as_map(i, length).spans)
        map = Map(spans=spans, parent_length=length)
    elif isinstance(slice, Map):
        map = slice
        # reasons for failure when the following is not commented out
        # should be checked further
        #assert map.parent_length == length, (map, length)
    else:
        (lo, hi, step) = _norm_slice(slice, length)
        assert (step or 1) == 1
        map = Map([(lo, hi)], parent_length=length)
    return map

class SpanI(object):
    """Abstract interface for Span and Range objects.
    
    Required properties: Start, End (must both be numbers)
    """
    __slots__ = []      #override in subclass
    
    def __contains__(self, other):
        """Returns True if other entirely contained in self."""
        raise NotImplementedError
    
    def overlaps(self, other):
        """Returns True if any positions in self are also in other."""
        raise NotImplementedError
    
    def reverse(self):
        """Reverses self."""
        raise NotImplementedError
    
    def __iter__(self):
        """Iterates over indices contained in self."""
        raise NotImplementedError
    
    def __str__(self):
        """Returns string representation of self."""
        return '(%s,%s)' % (self.Start, self.End)
    
    def __len__(self):
        """Returns length of self."""
        raise NotImplementedError
    
    def __cmp__(self):
        """Compares indices of self with indices of other."""
        raise NotImplementedError
    
    def startsBefore(self, other):
        """Returns True if self starts before other or other.Start."""
        try:
            return self.Start < other.Start
        except AttributeError:
            return self.Start < other
    
    def startsAfter(self, other):
        """Returns True if self starts after other or after other.Start."""
        try:
            return self.Start > other.Start
        except AttributeError:
            return self.Start > other
    
    def startsAt(self, other):
        """Returns True if self starts at the same place as other."""
        try:
            return self.Start == other.Start
        except AttributeError:
            return self.Start == other
    
    def startsInside(self, other):
        """Returns True if self's start in other or equal to other."""
        try:
            return self.Start in other
        except (AttributeError, TypeError):  #count other as empty span
            return False
    
    def endsBefore(self, other):
        """Returns True if self ends before other or other.End."""
        try:
            return self.End < other.End
        except AttributeError:
            return self.End < other
    
    def endsAfter(self, other):
        """Returns True if self ends after other or after other.End."""
        try:
            return self.End > other.End
        except AttributeError:
            return self.End > other
    
    def endsAt(self, other):
        """Returns True if self ends at the same place as other."""
        try:
            return self.End == other.End
        except AttributeError:
            return self.End == other
    
    def endsInside(self, other):
        """Returns True if self's end in other or equal to other."""
        try:
            return self.End in other
        except (AttributeError, TypeError):  #count other as empty span
            return False
    

class Span(SpanI):
    """A contiguous location, not much more than (start, end)
    
    Spans don't even know what map they are on.  The only smarts the class
    has is the ability to slice correctly.  Spans do not expect to be
    reverse-sliced (sl[5,3]) and treat positions as relative to themselves,
    not an underlying sequence (eg sl[:n] == sl[0:n]), so this slicing is
    very different to feature slicing.
    
    Spans may optionaly have a value, which gets preserved when they are remapped etc."""
    
    lost = False
    
    __slots__ = ( 'tidy_start', 'tidy_end', 'length', 'value',
            'Start', 'End', 'Reverse')
    
    def __init__(self, Start, End=None, tidy_start=False, tidy_end=False,
            value=None, Reverse=False):
        self._new_init(Start, End, Reverse)
        self.tidy_start = tidy_start
        self.tidy_end = tidy_end
        self.value = value
        self.length = self.End - self.Start
        assert self.length >= 0
    
    def _new_init(self, Start, End=None, Reverse=False):
        """Returns a new Span object, with Start, End, and Reverse properties.
        
        If End is not supplied, it is set to Start + 1 (providing a 1-element
        range).
        Reverse defaults to False.
        
        This should replace the current __init__ method when deprecated vars
        are removed.
        """
        #special handling in case we were passed another Span
        if isinstance(Start, Span):
            assert End is None
            self.Start, self.End, self.Reverse = Start.Start, Start.End, \
                Start.Reverse
        else:
            #reverse start and end so that start is always first
            if End is None:
                End = Start + 1
            elif Start > End:
                Start, End = End, Start
            
            self.Start = Start
            self.End = End
            self.Reverse = Reverse
    
    def __setstate__(self, args):
         self.__init__(*args)
    
    def __getstate__(self):
        return (self.Start, self.End, self.tidy_start, self.tidy_end, \
        self.value, self.Reverse)
    
    def __repr__(self):
        (start, end) = (self.Start, self.End)
        if self.Reverse:
            (end, start) = (start, end)
        return '%s:%s' % (start, end)
    
    def reversed(self):
        return self.__class__(self.Start, self.End, self.tidy_end, self.tidy_start, self.value, Reverse=not self.Reverse)
    
    def __getitem__(self, slice):
        start,end,step = _norm_slice(slice, self.length)
        assert (step or 1) == 1, slice
        assert start <= end, slice
        tidy_start = self.tidy_start and start==0
        tidy_end = self.tidy_end and end == self.length
        if self.Reverse:
            (Start, End, Reverse) = (self.End-end, self.End-start, True)
        else:
            (Start, End, Reverse) = (self.Start+start, self.Start+end, False)
        return type(self)(Start, End, tidy_start, tidy_end, self.value, Reverse)
    
    def __mul__(self, scale):
        return Span(self.Start * scale, self.End * scale,
                self.tidy_start, self.tidy_end, self.value, self.Reverse)
    
    def __div__(self, scale):
        assert not self.Start % scale or self.End % scale
        return Span(self.Start // scale, self.End // scale,
                self.tidy_start, self.tidy_end, self.value, self.Reverse)
    
    def remapWith(self, map):
        """The list of spans corresponding to this span on its grandparent, ie:
        C is a span of a feature on B which itself is a feature on A, so to
        place C on A return that part of B (map) covered by C (self)"""
        
        (offsets, spans) = (map.offsets, map.spans)
        map_length = offsets[-1] + spans[-1].length
        
        # don't try to remap any non-corresponding end region(s)
        # this won't matter if all spans lie properly within their
        # parent maps, but that might not be true of Display slices.
        (zlo, zhi) = (max(0, self.Start), min(map_length, self.End))
        
        # Find the right span(s) of the map
        first = bisect_right(offsets, zlo) - 1
        last = bisect_left(offsets, zhi, first) -1
        result = spans[first:last+1]
        
        # Cut off something at either end to get
        # the same position and length as 'self'
        if result:
            end_trim = offsets[last] + spans[last].length - zhi
            start_trim = zlo - offsets[first]
            if end_trim > 0:
                result[-1] = result[-1][:result[-1].length-end_trim]
            if start_trim > 0:
                result[0] = result[0][start_trim:]
        
        # May need to add a bit at either end if the span didn't lie entirely
        # within its parent map (eg: Display slice, inverse of feature map).
        if self.Start < 0:
            result.insert(0, LostSpan(-self.Start))
        if self.End > map_length:
            result.append(LostSpan(self.End-map_length))
        
        # If the ends of self are meaningful then so are the new ends,
        # but not any new internal breaks.
        if result:
            if self.tidy_start:
                result[0].tidy_start = True
            if self.tidy_end:
                result[-1].tidy_end = True
        
        # Deal with case where self is a reverse slice.
        if self.Reverse:
            result = [part.reversed() for part in result]
            result.reverse()
        
        if self.value is not None:
            result = [copy.copy(s) for s in result]
            for s in result: s.value = self.value
        
        return result
    
    def __contains__(self, other):
        """Returns True if other completely contained in self.
        
        other must either be a number or have Start and End properties.
        """
        try:
            return other.Start >= self.Start and other.End <= self.End
        except AttributeError:
            #other is scalar: must be _less_ than self.End,
            #for the same reason that 3 is not in range(3).
            return other >= self.Start and other < self.End
    
    def overlaps(self, other):
        """Returns True if any positions in self are also in other."""
        #remember to subtract 1 from the Ends, since self.End isn't really
        #in self...
        try:
            return (self.Start in other) or (other.Start in self)
        except AttributeError:  #other was probably a number?
            return other in self
    
    def reverse(self):
        """Reverses self."""
        self.Reverse = not self.Reverse
    
    def reversedRelativeTo(self, length):
        """Returns a new span with positions adjusted relative to length. For
        use in reverse complementing of nucleic acids"""
        
        # if reverse complementing, the start becomes the length minus the end
        # position
        start = length - self.End
        assert start >= 0
        end = start + self.length
        return self.__class__(start, end, value = self.value,
                Reverse = not self.Reverse)
    
    def __iter__(self):
        """Iterates over indices contained in self.
        
        NOTE: to make sure that the same items are contained whether going
        through the range in forward or reverse, need to adjust the indices
        by 1 if going backwards.
        """
        if self.Reverse:
            return iter(xrange(self.End-1, self.Start-1, -1))
        else:
            return iter(xrange(self.Start, self.End, 1))
    
    def __str__(self):
        """Returns string representation of self."""
        return '(%s,%s,%s)' % (self.Start, self.End, bool(self.Reverse))
    
    def __len__(self):
        """Returns length of self."""
        return self.End - self.Start
    
    def __cmp__(self, other):
        """Compares indices of self with indices of other."""
        if hasattr(other, 'Start') and hasattr(other, 'End'):
            return cmp(self.Start, other.Start) or cmp(self.End, other.End) \
                or cmp(self.Reverse, other.Reverse)
        else:
            return cmp(type(self), type(other))
    

class _LostSpan(object):
    """A placeholder span which doesn't exist in the underlying sequence"""
    
    __slots__ = ['length', 'value']
    lost = True
    terminal = False
    
    def __init__(self, length, value=None):
        self.length = length
        self.value = value
    
    def __len__(self):
        return self.length
    
    def __setstate__(self, args):
         self.__init__(*args)
    
    def __getstate__(self):
        return (self.length, self.value)
    
    def __repr__(self):
        return '-%s-' % (self.length)
    
    def where(self, index):
        return None
    
    def reversed(self):
        return self
    
    def __getitem__(self, slice):
        (start,end,step) = _norm_slice(slice, self.length)
        assert (step or 1) == 1, slice
        return self.__class__(abs(end-start), self.value)
    
    def __mul__(self, scale):
        return LostSpan(self.length * scale, self.value)
    
    def __div__(self, scale):
        assert not self.length % 3
        return LostSpan(self.length // scale, self.value)
    
    def remapWith(self, map):
        return [self]
    
    def reversedRelativeTo(self, length):
        return self
    

# Save memory by only making one of each small gap
_lost_span_cache = {}
def LostSpan(length, value=None):
    global _lost_span_cache
    if value is None and length < 1000:
        if length not in _lost_span_cache:
            _lost_span_cache[length] = _LostSpan(length, value)
        return _lost_span_cache[length]
    else:
        return _LostSpan(length, value)

class TerminalPadding(_LostSpan):
    terminal = True
    def __repr__(self):
        return '?%s?' % (self.length)
    

class Map(object):
    """A map holds a list of spans.  """
    
    def __init__(self, locations=None, spans=None, tidy=False,
            parent_length=None, termini_unknown=False):
        assert parent_length is not None
        if spans is None:
            spans = []
            for (start, end) in locations:
                diff = 0
                reverse = start > end
                if max(start, end) < 0 or min(start, end) > parent_length:
                    raise RuntimeError("located outside sequence: %s" % \
                                    str((start, end, parent_length)))
                elif max(start, end) < 0:
                    diff = min(start, end)
                    start = [start, 0][start < 0]
                    end = [end, 0][end < 0]
                elif min(start, end) > parent_length:
                    diff = max(start, end) - parent_length
                    start = [start, parent_length][start > parent_length]
                    end = [end, parent_length][end > parent_length]
                
                span = Span(start, end, tidy, tidy, Reverse=reverse)
                if diff < 0:
                    spans += [LostSpan(-diff), span]
                elif diff > 0:
                    spans += [span, LostSpan(diff)]
                else:
                    spans += [span]
        
        self.offsets = []
        self.useful = False
        self.complete = True
        self.Reverse = None
        posn = 0
        for span in spans:
            self.offsets.append(posn)
            posn += span.length
            if span.lost:
                self.complete = False
            elif not self.useful:
                self.useful = True
                (self.Start, self.End) = (span.Start, span.End)
                self.Reverse = span.Reverse
            else:
                self.Start = min(self.Start, span.Start)
                self.End = max(self.End, span.End)
                if self.Reverse is not None and (span.Reverse != self.Reverse):
                    self.Reverse = None
        
        if termini_unknown:
            if spans[0].lost:
                spans[0] = TerminalPadding(spans[0].length)
            if spans[-1].lost:
                spans[-1] = TerminalPadding(spans[-1].length)
        
        self.spans = spans
        self.length = posn
        self.parent_length = parent_length
        self.__inverse = None
    
    def __len__(self):
        return self.length
    
    def __repr__(self):
        return repr(self.spans) + '/%s' % self.parent_length
    
    def __getitem__(self, slice):
        # A possible shorter map at the same level
        slice = as_map(slice, len(self))
        new_parts = []
        for span in slice.spans:
            new_parts.extend(span.remapWith(self))
        return Map(spans=new_parts, parent_length=self.parent_length)
    
    def __mul__(self, scale):
        # For Protein -> DNA
        new_parts = []
        for span in self.spans:
            new_parts.append(span * scale)
        return Map(spans=new_parts, parent_length=self.parent_length*scale)
    
    def __div__(self, scale):
        # For DNA -> Protein
        new_parts = []
        for span in self.spans:
            new_parts.append(span / scale)
        return Map(spans=new_parts, parent_length=self.parent_length // scale)
    
    def __add__(self, other):
        if other.parent_length != self.parent_length:
            raise ValueError("Those maps belong to different sequences")
        return Map(spans=self.spans + other.spans, parent_length=self.parent_length)
        
    def withTerminiUnknown(self):
        return Map(self, spans=self.spans[:],
                parent_length=self.parent_length,
                termini_unknown = True)
    
    def getCoveringSpan(self):
        if self.Reverse:
            span = (self.End, self.Start)
        else:
            span = (self.Start, self.End)
        return Map([span], parent_length=self.parent_length)
    
    def covered(self):
        """>>> Map([(10,20), (15, 25), (80, 90)]).covered().spans
        [Span(10,25), Span(80, 90)]"""
        
        delta = {}
        for span in self.spans:
            if span.lost:
                continue
            delta[span.Start] = delta.get(span.Start, 0) + 1
            delta[span.End] = delta.get(span.End, 0) - 1
        positions = delta.keys()
        positions.sort()
        last_y = y = 0
        last_x = start = None
        result = []
        for x in positions:
            y += delta[x]
            if x == last_x:
                continue
            if y and not last_y:
                assert start is None
                start = x
            elif last_y and not y:
                result.append((start, x))
                start = None
            last_x = x
            last_y = y
        assert y == 0
        return Map(result, parent_length=self.parent_length)
    
    def reversed(self):
        """Reversed location on same parent"""
        spans = [s.reversed() for s in self.spans]
        spans.reverse()
        return Map(spans=spans, parent_length=self.parent_length)
    
    def nucleicReversed(self):
        """Same location on reversed parent"""
        spans = [s.reversedRelativeTo(self.parent_length) for s in self.spans]
        return Map(spans=spans, parent_length=self.parent_length)
    
    def gaps(self):
        """The gaps (lost spans) in this map"""
        locations = []
        offset = 0
        for s in self.spans:
            if s.lost:
                locations.append((offset, offset+s.length))
            offset += s.length
        return Map(locations, parent_length=len(self))
    
    def shadow(self):
        """The 'negative' map of the spans not included in this map"""
        return self.inverse().gaps()
    
    def nongap(self):
        locations = []
        offset = 0
        for s in self.spans:
            if not s.lost:
                locations.append((offset, offset+s.length))
            offset += s.length
        return Map(locations, parent_length=len(self))
    
    def withoutGaps(self):
        return Map(
                spans = [s for s in self.spans if not s.lost],
                parent_length = self.parent_length)
    
    def inverse(self):
        if self.__inverse is None:
            self.__inverse = self._inverse()
        return self.__inverse
    
    def _inverse(self):
        # can't work if there are overlaps in the map
        # tidy ends don't survive inversion
        if self.parent_length is None:
            raise ValueError("Uninvertable. Parent length not known")
        posn = 0
        temp = []
        for span in self.spans:
            if not span.lost:
                if span.Reverse:
                    temp.append((span.Start, span.End, posn+span.length, posn))
                else:
                    temp.append((span.Start, span.End, posn, posn+span.length))
            posn += span.length
        
        temp.sort()
        new_spans = []
        last_hi = 0
        for (lo, hi, start, end) in temp:
            if lo > last_hi:
                new_spans.append(LostSpan(lo-last_hi))
            elif lo < last_hi:
                raise ValueError, "Uninvertable. Overlap: %s < %s" % (lo, last_hi)
            new_spans.append(Span(start, end, Reverse=start>end))
            last_hi = hi
        if self.parent_length > last_hi:
            new_spans.append(LostSpan(self.parent_length-last_hi))
        return Map(spans=new_spans, parent_length=len(self))
    

class SpansOnly(ConstrainedList):
    """List that converts elements to Spans on addition."""
    Mask = FunctionWrapper(Span)
    _constraint = ClassChecker(Span)

class Range(SpanI):
    """Complex object consisting of many spans."""
    
    def __init__(self, Spans=[]):
        """Returns a new Range object with data in Spans.
        """
        result = SpansOnly()
        #need to check if we got a single Span, since they define __iter__.
        if isinstance(Spans, Span):
            result.append(Spans)
        elif hasattr(Spans, 'Spans'):   #probably a single range object?
            result.extend(Spans.Spans)
        else:
            for s in iterable(Spans):
                if hasattr(s, 'Spans'):
                  result.extend(s.Spans)
                else:
                    result.append(s)
        self.Spans = result
    
    def __str__(self):
        """Returns string representation of self."""
        return '(%s)' % ','.join(map(str, self.Spans))
    
    def __len__(self):
        """Returns sum of span lengths.
        
        NOTE: if spans overlap, will count multiple times. Use reduce() to
        get rid of overlaps.
        """
        return sum(map(len, self.Spans))
    
    def __cmp__(self, other):
        """Compares spans of self with indices of other."""
        if hasattr(other, 'Spans'):
            return cmp(self.Spans, other.Spans)
        elif len(self.Spans) == 1 and hasattr(other, 'Start') and \
            hasattr(other, 'End'):
            return cmp(self.Spans[0].Start, other.Start) or \
                cmp(self.Spans[0].End, other.End)
        else:
            return object.__cmp__(self, other)
    
    def _get_start(self):
        """Finds earliest start of items in self.Spans."""
        return min([i.Start for i in self.Spans])
    Start = property(_get_start)
    
    def _get_end(self):
        """Finds latest end of items in self.Spans."""
        return max([i.End for i in self.Spans])
    End = property(_get_end)
    
    def _get_reverse(self):
        """Reverse is True if any piece is reversed."""
        for i in self.Spans:
            if i.Reverse:
                return True
        return False
    Reverse = property(_get_reverse)
    
    def reverse(self):
        """Reverses all spans in self."""
        for i in self.Spans:
            i.reverse()
    
    def __contains__(self, other):
        """Returns True if other completely contained in self.
        
        other must either be a number or have Start and End properties.
        """
        if hasattr(other, 'Spans'):
            for curr in other.Spans:
                found = False
                for i in self.Spans:
                    if curr in i:
                        found = True
                        break
                if not found:
                    return False
            return True
        else:
            for i in self.Spans:
                if other in i:
                    return True
            return False
    
    def overlaps(self, other):
        """Returns True if any positions in self are also in other."""
        if hasattr(other, 'Spans'):
            for i in self.Spans:
                for j in other.Spans:
                    if i.overlaps(j):
                        return True
        else:
            for i in self.Spans:
                if i.overlaps(other):
                    return True
        return False
    
    def overlapsExtent(self, other):
        """Returns True if any positions in self's extent also in other's."""
        if hasattr(other, 'Extent'):
            return self.Extent.overlaps(other.Extent)
        else:
            return self.Extent.overlaps(other)
    
    def sort(self):
        """Sorts the spans in self."""
        self.Spans.sort()
    
    def __iter__(self):
        """Iterates over indices contained in self."""
        return chain(*[iter(i) for i in self.Spans])
    
    def _get_extent(self):
        """Returns Span object representing the extent of self."""
        return Span(self.Start, self.End)
    Extent = property(_get_extent)
    
    def simplify(self):
        """Reduces the spans in self in-place to get fewest spans.
        
        Will not condense spans with opposite directions.
        
        Will condense adjacent but nonoverlapping spans (e.g. (1,3) and (4,5)).
        """
        forward = []
        reverse = []
        spans = self.Spans[:]
        spans.sort()
        for span in spans:
            if span.Reverse:
                direction = reverse
            else:
                direction = forward
            
            found_overlap = False
            for other in direction:
                if span.overlaps(other) or (span.Start == other.End) or \
                    (other.Start == span.End):  #handle adjacent spans also
                    other.Start = min(span.Start, other.Start)
                    other.End = max(span.End, other.End)
                    found_overlap = True
                    break
            if not found_overlap:
                direction.append(span)
        self.Spans[:] = forward + reverse
    

class Point(Span):
    """Point is a special case of Span, where Start always equals End.
    
    Note that, as per Python standard, a point is _between_ two elements
    in a sequence. In other words, a point does not contain any elements.
    If you want a single element, use a Span where End = Start + 1.
    
    A Point does have a direction (i.e. a Reverse property) to indicate
    where successive items would go if it were expanded.
    """
    
    def __init__(self, Start, Reverse=False):
        """Returns new Point object."""
        self.Reverse = Reverse
        self._start = Start
    
    def _get_start(self):
        """Returns self.Start."""
        return self._start
    
    def _set_start(self, Start):
        """Sets self.Start and self.End."""
        self._start = Start
    
    Start = property(_get_start, _set_start)
    End = Start     #start and end are synonyms for the same property

def RangeFromString(string, delimiter=','):
    """Returns Range object from string of the form 1-5,11,20,30-50.
    
    Ignores whitespace; expects values to be comma-delimited and positive.
    """
    result = Range()
    pairs = map(strip, string.split(delimiter))
    for p in pairs:
        if not p:   #adjacent delimiters?
            continue
        if '-' in p:    #treat as pair
            first, second = p.split('-')
            result.Spans.append(Span(int(first), int(second)))
        else:
            result.Spans.append(Span(int(p)))
    return result

