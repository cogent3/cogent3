"""
Maps have a list of Spans.

Span is a region with a start, an end, and a direction.

Notes
-----

Spans behave much like Python's slices: a Span contains the element
after its start but does not contain the element after its end. It may help to
think of the Span indices occurring _between_ the list elements:

    a b c d e
   | | | | | |
   0 1 2 3 4 5

...so that a Span whose start is its end contains no elements (e.g. 2:2), and
a Span whose end is 2 more than its start contains 2 elements (e.g. 2:4 has c
and d), etc. Similarly, Span(0,2) does _not_ overlap Span(2,3), since the
former contains a and b while the latter contains c.
"""
import copy
import dataclasses
import inspect

from abc import ABC, abstractmethod
from bisect import bisect_left, bisect_right
from functools import total_ordering
from itertools import chain
from typing import Iterator, List, Optional, Sequence, Tuple, Union

from numpy import array, ndarray

from cogent3._version import __version__
from cogent3.util import warning as c3warn
from cogent3.util.misc import (
    ClassChecker,
    ConstrainedList,
    FunctionWrapper,
    get_object_provenance,
    iterable,
)


strip = str.strip


def _norm_index(i, length, default):
    """For converting s[:3] to s[0:3], s[-1] to s[len(s)-1] and s[0:lots] to s[0:len(s)]"""
    if i is None:
        i = default
    elif i < 0:
        i += length
    return min(max(i, 0), length)


def _norm_slice(index, length):
    """_norm_slice(slice(1, -2, 3), 10) -> (1,8,3)"""
    if isinstance(index, slice):
        start = _norm_index(index.start, length, 0)
        end = _norm_index(index.stop, length, length)
        return (start, end, index.step)
    else:
        start = index
        if start < 0:
            start += length
        if start >= length:
            raise IndexError(index)
        return (start, start + 1, 1)


def as_map(slice, length, cls):
    """Take anything that might be used as a subscript: Integer, Slice,
    or MapABC, and return cls."""

    if isinstance(slice, (list, tuple)):
        spans = []
        for i in slice:
            spans.extend(as_map(i, length, cls).spans)
        return cls(spans=spans, parent_length=length)
    elif isinstance(slice, (FeatureMap, IndelMap)):
        return slice
    else:
        lo, hi, step = _norm_slice(slice, length)
        assert (step or 1) == 1
        return cls(locations=[(lo, hi)], parent_length=length)


class SpanI(object):
    """Abstract interface for Span and Range objects.

    Required properties: start, end (must both be numbers)
    """

    __slots__ = []  # override in subclass

    def __contains__(self, other):
        """Returns True if other entirely contained in self."""
        raise NotImplementedError

    def overlaps(self, other):
        """Returns True if any positions in self are also in other."""
        raise NotImplementedError

    def reverses(self):
        """Reverses self."""
        raise NotImplementedError

    def __iter__(self):
        """Iterates over indices contained in self."""
        raise NotImplementedError

    def __str__(self):
        """Returns string representation of self."""
        return f"({self.start},{self.end})"

    def __len__(self):
        """Returns length of self."""
        raise NotImplementedError

    def __lt__(self, other):
        """Compares indices of self with indices of other."""
        raise NotImplementedError

    def starts_before(self, other):
        """Returns True if self starts before other or other.start."""
        try:
            return self.start < other.start
        except AttributeError:
            return self.start < other

    def starts_after(self, other):
        """Returns True if self starts after other or after other.start."""
        try:
            return self.start > other.start
        except AttributeError:
            return self.start > other

    def starts_at(self, other):
        """Returns True if self starts at the same place as other."""
        try:
            return self.start == other.start
        except AttributeError:
            return self.start == other

    def starts_inside(self, other):
        """Returns True if self's start in other or equal to other."""
        try:
            return self.start in other
        except (AttributeError, TypeError):  # count other as empty span
            return False

    def ends_before(self, other):
        """Returns True if self ends before other or other.end."""
        try:
            return self.end < other.end
        except AttributeError:
            return self.end < other

    def ends_after(self, other):
        """Returns True if self ends after other or after other.end."""
        try:
            return self.end > other.end
        except AttributeError:
            return self.end > other

    def ends_at(self, other):
        """Returns True if self ends at the same place as other."""
        try:
            return self.end == other.end
        except AttributeError:
            return self.end == other

    def ends_inside(self, other):
        """Returns True if self's end in other or equal to other."""
        try:
            return self.end in other
        except (AttributeError, TypeError):  # count other as empty span
            return False


@total_ordering
class Span(SpanI):
    """A contiguous location, not much more than (start, end)

    Spans don't even know what map they are on.  The only smarts the class
    has is the ability to slice correctly.  Spans do not expect to be
    reverse-sliced (sl[5,3]) and treat positions as relative to themselves,
    not an underlying sequence (eg sl[:n] == sl[0:n]), so this slicing is
    very different to feature slicing.

    Spans may optionaly have a value, which gets preserved when they are remapped etc.
    """

    lost = False

    __slots__ = (
        "tidy_start",
        "tidy_end",
        "length",
        "value",
        "start",
        "end",
        "reverse",
        "_serialisable",
    )

    def __init__(
        self,
        start,
        end=None,
        tidy_start=False,
        tidy_end=False,
        value=None,
        reverse=False,
    ):
        d = locals()
        x = ("self", "__class__", "__slots__")
        self._serialisable = {k: v for k, v in d.items() if k not in x}

        self._new_init(start, end, reverse)
        self.tidy_start = tidy_start
        self.tidy_end = tidy_end
        self.value = value
        self.length = self.end - self.start
        assert self.length >= 0

    def _new_init(self, start, end=None, reverse=False):
        """Returns a new Span object, with start, end, and reverse properties.

        If end is not supplied, it is set to start + 1 (providing a 1-element
        range).
        reverse defaults to False.
        """
        # This should replace the current __init__ method when deprecated vars
        # are removed.
        # special handling in case we were passed another Span
        if isinstance(start, Span):
            assert end is None
            self.start, self.end, self.reverse = start.start, start.end, start.reverse
        else:
            # reverse start and end so that start is always first
            if end is None:
                end = start + 1
            elif start > end:
                start, end = end, start

            self.start = start
            self.end = end
            self.reverse = reverse

    def to_rich_dict(self):
        attribs = copy.deepcopy(self._serialisable)
        attribs["type"] = get_object_provenance(self)
        attribs["version"] = __version__
        return attribs

    def __setstate__(self, args):
        self.__init__(*args)

    def __getstate__(self):
        return (
            self.start,
            self.end,
            self.tidy_start,
            self.tidy_end,
            self.value,
            self.reverse,
        )

    def __repr__(self):
        (start, end) = (self.start, self.end)
        if self.reverse:
            (end, start) = (start, end)
        return f"{start}:{end}"

    def reversed(self):
        return self.__class__(
            self.start,
            self.end,
            self.tidy_end,
            self.tidy_start,
            self.value,
            reverse=not self.reverse,
        )

    def __getitem__(self, slice):
        start, end, step = _norm_slice(slice, self.length)
        assert (step or 1) == 1, slice
        assert start <= end, slice
        tidy_start = self.tidy_start and start == 0
        tidy_end = self.tidy_end and end == self.length
        if self.reverse:
            (start, end, reverse) = (self.end - end, self.end - start, True)
        else:
            (start, end, reverse) = (self.start + start, self.start + end, False)
        return type(self)(start, end, tidy_start, tidy_end, self.value, reverse)

    def __mul__(self, scale):
        return Span(
            self.start * scale,
            self.end * scale,
            self.tidy_start,
            self.tidy_end,
            self.value,
            self.reverse,
        )

    def __div__(self, scale):
        assert not self.start % scale or self.end % scale
        return Span(
            self.start // scale,
            self.end // scale,
            self.tidy_start,
            self.tidy_end,
            self.value,
            self.reverse,
        )

    def remap_with(self, map):
        """The list of spans corresponding to this span on its grandparent, ie:
        C is a span of a feature on B which itself is a feature on A, so to
        place C on A return that part of B (map) covered by C (self)"""

        offsets, spans = map.offsets, list(map.spans)
        map_length = offsets[-1] + spans[-1].length

        # don't try to remap any non-corresponding end region(s)
        # this won't matter if all spans lie properly within their
        # parent maps, but that might not be true of Display slices.
        zlo, zhi = max(0, self.start), min(map_length, self.end)

        # Find the right span(s) of the map
        first = bisect_right(offsets, zlo) - 1
        last = bisect_left(offsets, zhi, first) - 1
        result = spans[first : last + 1]

        # Cut off something at either end to get
        # the same position and length as 'self'
        if result:
            end_trim = offsets[last] + spans[last].length - zhi
            start_trim = zlo - offsets[first]
            if end_trim > 0:
                result[-1] = result[-1][: result[-1].length - end_trim]
            if start_trim > 0:
                result[0] = result[0][start_trim:]

        # May need to add a bit at either end if the span didn't lie entirely
        # within its parent map (eg: Display slice, inverse of feature map).
        if self.start < 0:
            result.insert(0, LostSpan(-self.start))
        if self.end > map_length:
            result.append(LostSpan(self.end - map_length))

        # If the ends of self are meaningful then so are the new ends,
        # but not any new internal breaks.
        if result:
            if self.tidy_start:
                result[0].tidy_start = True
            if self.tidy_end:
                result[-1].tidy_end = True

        # Deal with case where self is a reverse slice.
        if self.reverse:
            result = [part.reversed() for part in result]
            result.reverse()

        if self.value is not None:
            result = [copy.copy(s) for s in result]
            for s in result:
                s.value = self.value

        return result

    def __contains__(self, other):
        """Returns True if other completely contained in self.

        other must either be a number or have start and end properties.
        """
        try:
            return other.start >= self.start and other.end <= self.end
        except AttributeError:
            # other is scalar: must be _less_ than self.end,
            # for the same reason that 3 is not in range(3).
            return other >= self.start and other < self.end

    def overlaps(self, other):
        """Returns True if any positions in self are also in other."""
        # remember to subtract 1 from the Ends, since self.end isn't really
        # in self...
        try:
            return (self.start in other) or (other.start in self)
        except AttributeError:  # other was probably a number?
            return other in self

    def reverses(self):
        """Reverses self."""
        self.reverse = not self.reverse

    def reversed_relative_to(self, length):
        """Returns a new span with positions adjusted relative to length. For
        use in reverse complementing of nucleic acids"""

        # if reverse complementing, the start becomes the length minus the end
        # position
        start = length - self.end
        assert start >= 0
        end = start + self.length
        return self.__class__(start, end, value=self.value, reverse=not self.reverse)

    def __iter__(self):
        """Iterates over indices contained in self.

        NOTE: to make sure that the same items are contained whether going
        through the range in forward or reverse, need to adjust the indices
        by 1 if going backwards.
        """
        if self.reverse:
            return iter(range(self.end - 1, self.start - 1, -1))
        else:
            return iter(range(self.start, self.end, 1))

    def __str__(self):
        """Returns string representation of self."""
        return f"({self.start},{self.end},{bool(self.reverse)})"

    def __len__(self):
        """Returns length of self."""
        return self.end - self.start

    def __lt__(self, other):
        """Compares indices of self with indices of other."""
        if hasattr(other, "start") and hasattr(other, "end"):
            s = (self.start, self.end, self.reverse)
            o = (other.start, other.end, other.reverse)
            return s < o
        else:
            return type(self) < type(other)

    def __eq__(self, other):
        """Compares indices of self with indices of other."""
        if hasattr(other, "start") and hasattr(other, "end"):
            return (
                self.start == other.start
                and self.end == other.end
                and self.reverse == other.reverse
            )
        else:
            return type(self) == type(other)


class _LostSpan:
    """A placeholder span which doesn't exist in the underlying sequence"""

    __slots__ = ["length", "value", "_serialisable"]
    lost = True
    terminal = False

    def __init__(self, length, value=None):
        d = locals()
        exclude = ("self", "__class__", "__slots__")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}

        self.length = length
        self.value = value

    def to_rich_dict(self):
        attribs = copy.deepcopy(self._serialisable)
        attribs["type"] = get_object_provenance(self)
        attribs["version"] = __version__
        return attribs

    def __len__(self):
        return self.length

    def __setstate__(self, args):
        self.__init__(*args)

    def __getstate__(self):
        return (self.length, self.value)

    def __repr__(self):
        return f"-{self.length}-"

    def where(self, index):
        return None

    def reversed(self):
        return self

    def __getitem__(self, slice):
        (start, end, step) = _norm_slice(slice, self.length)
        assert (step or 1) == 1, slice
        return self.__class__(abs(end - start), self.value)

    def __mul__(self, scale):
        return LostSpan(self.length * scale, self.value)

    def __div__(self, scale):
        assert not self.length % 3
        return LostSpan(self.length // scale, self.value)

    def remap_with(self, map):
        return [self]

    def reversed_relative_to(self, length):
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
        return f"?{self.length}?"


class Map:  # pragma: no cover
    """A map holds a list of spans."""

    def __init__(
        self,
        locations=None,
        spans=None,
        tidy=False,
        parent_length=None,
        termini_unknown=False,
    ):
        assert parent_length is not None
        d = locals()
        exclude = ("self", "__class__", "__slots__")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}

        if spans is None:
            spans = []
            for start, end in locations:
                diff = 0
                reverse = start > end
                if max(start, end) < 0 or min(start, end) > parent_length:
                    raise RuntimeError(
                        f"located outside sequence: {(start, end, parent_length)}"
                    )
                if max(start, end) > parent_length and min(start, end) < 0:
                    l_diff = min(start, end)
                    r_diff = max(start, end) - parent_length
                    start, end = (
                        (0, parent_length) if start < end else (parent_length, 0)
                    )
                    spans += [
                        LostSpan(abs(l_diff)),
                        Span(start, end, tidy, tidy, reverse=reverse),
                        LostSpan(abs(r_diff)),
                    ]
                elif min(start, end) < 0:
                    diff = min(start, end)
                    start = 0 if start < 0 else start
                    end = 0 if end < 0 else end
                    spans += [
                        LostSpan(abs(diff)),
                        Span(start, end, tidy, tidy, reverse=reverse),
                    ]
                elif max(start, end) > parent_length:
                    diff = max(start, end) - parent_length
                    start = parent_length if start > parent_length else start
                    end = parent_length if end > parent_length else end
                    spans += [
                        Span(start, end, tidy, tidy, reverse=reverse),
                        LostSpan(abs(diff)),
                    ]
                else:
                    spans += [Span(start, end, tidy, tidy, reverse=reverse)]

        self.offsets = []
        self.useful = False
        self.complete = True
        self.reverse = None
        posn = 0
        for span in spans:
            self.offsets.append(posn)
            posn += span.length
            if span.lost:
                self.complete = False
            elif not self.useful:
                self.useful = True
                (self.start, self.end) = (span.start, span.end)
                self.reverse = span.reverse
            else:
                self.start = min(self.start, span.start)
                self.end = max(self.end, span.end)
                if self.reverse is not None and (span.reverse != self.reverse):
                    self.reverse = None

        if termini_unknown:
            spans = list(spans)
            if spans[0].lost:
                spans[0] = TerminalPadding(spans[0].length)
            if spans[-1].lost:
                spans[-1] = TerminalPadding(spans[-1].length)

        self.spans = tuple(spans)
        self.length = posn
        self.parent_length = parent_length
        self.__inverse = None

    def __len__(self):
        return self.length

    def __repr__(self):
        return repr(list(self.spans)) + f"/{self.parent_length}"

    def __getitem__(self, slice):
        # A possible shorter map at the same level
        slice = as_map(slice, len(self), self.__class__)
        new_parts = []
        for span in slice.spans:
            new_parts.extend(span.remap_with(self))
        return Map(spans=new_parts, parent_length=self.parent_length)

    def __mul__(self, scale):
        # For Protein -> DNA
        new_parts = []
        for span in self.spans:
            new_parts.append(span * scale)
        return Map(spans=new_parts, parent_length=self.parent_length * scale)

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

    def with_termini_unknown(self):
        return Map(
            self,
            spans=self.spans[:],
            parent_length=self.parent_length,
            termini_unknown=True,
        )

    def get_covering_span(self):
        if self.reverse:
            span = (self.end, self.start)
        else:
            span = (self.start, self.end)
        return Map([span], parent_length=self.parent_length)

    def covered(self):
        """>>> Map([(10,20), (15, 25), (80, 90)]).covered().spans
        [Span(10,25), Span(80, 90)]"""

        delta = {}
        for span in self.spans:
            if span.lost:
                continue
            delta[span.start] = delta.get(span.start, 0) + 1
            delta[span.end] = delta.get(span.end, 0) - 1
        positions = list(delta.keys())
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
        return Map(locations=result, parent_length=self.parent_length)

    def reversed(self):
        """Reversed location on same parent"""
        spans = [s.reversed() for s in self.spans]
        spans.reverse()
        return Map(spans=spans, parent_length=self.parent_length)

    def nucleic_reversed(self):
        """Same location on reversed parent"""
        spans = [s.reversed_relative_to(self.parent_length) for s in self.spans]
        return Map(spans=spans, parent_length=self.parent_length)

    def get_gap_coordinates(self):
        """returns [(gap pos, gap length), ...]"""
        gap_pos = []
        for i, span in enumerate(self.spans):
            if not span.lost:
                continue

            pos = self.spans[i - 1].end if i else 0
            gap_pos.append((pos, len(span)))

        return gap_pos

    def gaps(self):
        """The gaps (lost spans) in this map"""
        locations = []
        offset = 0
        for s in self.spans:
            if s.lost:
                locations.append((offset, offset + s.length))
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
                locations.append((offset, offset + s.length))
            offset += s.length
        return Map(locations, parent_length=len(self))

    def without_gaps(self):
        return Map(
            spans=[s for s in self.spans if not s.lost],
            parent_length=self.parent_length,
        )

    def inverse(self):
        if self.__inverse is None:
            self.__inverse = self._inverse()
        return self.__inverse

    def _inverse(self):
        # can't work if there are overlaps in the map
        # tidy ends don't survive inversion
        if self.parent_length is None:
            raise ValueError("Uninvertable. parent length not known")
        posn = 0
        temp = []
        for span in self.spans:
            if not span.lost:
                if span.reverse:
                    temp.append((span.start, span.end, posn + span.length, posn))
                else:
                    temp.append((span.start, span.end, posn, posn + span.length))
            posn += span.length

        temp.sort()
        new_spans = []
        last_hi = 0
        for lo, hi, start, end in temp:
            if lo > last_hi:
                new_spans.append(LostSpan(lo - last_hi))
            elif lo < last_hi:
                raise ValueError(f"Uninvertable. Overlap: {lo} < {last_hi}")
            new_spans.append(Span(start, end, reverse=start > end))
            last_hi = hi
        if self.parent_length > last_hi:
            new_spans.append(LostSpan(self.parent_length - last_hi))
        return Map(spans=new_spans, parent_length=len(self))

    def get_coordinates(self):
        """returns span coordinates as [(v1, v2), ...]

        v1/v2 are (start, end) unless the map is reversed, in which case it will
        be (end, start)"""

        if self.reverse:
            order_func = lambda x: (max(x), min(x))
        else:
            order_func = lambda x: x

        coords = list(
            map(order_func, [(s.start, s.end) for s in self.spans if not s.lost])
        )

        return coords

    def to_rich_dict(self):
        """returns dicts for contained spans [dict(), ..]"""
        spans = [s.to_rich_dict() for s in self.spans]
        data = copy.deepcopy(self._serialisable)
        data.pop("locations")
        data["spans"] = spans
        data["type"] = get_object_provenance(self)
        data["version"] = __version__
        return data

    def zeroed(self):
        """returns a new instance with the first span starting at 0

        Note
        ----

        Useful when an Annotatable object is sliced, but the connection to
        the original parent is being deliberately broken as in the
        Sequence.deepcopy(sliced=True) case.
        """
        # todo there's probably a more efficient way to do this
        # create the new instance
        from cogent3.util.deserialise import deserialise_map_spans

        data = self.to_rich_dict()
        zeroed = deserialise_map_spans(data)
        zeroed.parent_length = len(zeroed.get_covering_span())
        shift = min(zeroed.start, zeroed.end)
        new_end = 0
        for span in zeroed.spans:
            if span.lost:
                continue
            span.start -= shift
            span.end -= shift
            new_end = max(new_end, span.end)

        zeroed.start = 0
        zeroed.end = new_end

        return zeroed

    T = Union[ndarray, int]

    def absolute_position(self, rel_pos: T) -> T:
        """converts rel_pos into an absolute position

        Raises
        ------
        raises ValueError if rel_pos < 0
        """
        check = array([rel_pos], dtype=int) if isinstance(rel_pos, int) else rel_pos
        if check.min() < 0:
            raise ValueError(f"must positive, not {rel_pos=}")

        if len(self) == self.parent_length:
            # handle case of reversed here?
            return rel_pos

        return self.start + rel_pos

    def relative_position(self, abs_pos: T) -> T:
        """converts abs_pos into an relative position

        Raises
        ------
        raises ValueError if abs_pos < 0
        """
        check = array([abs_pos], dtype=int) if isinstance(abs_pos, int) else abs_pos
        if check.min() < 0:
            raise ValueError(f"must positive, not {abs_pos=}")
        return abs_pos - self.start


@c3warn.deprecated_callable(
    version="2024.3", reason="No being used.", is_discontinued=True
)
class SpansOnly(ConstrainedList):  # pragma: no cover
    """List that converts elements to Spans on addition."""

    mask = FunctionWrapper(Span)
    _constraint = ClassChecker(Span)


@c3warn.deprecated_callable(
    version="2024.3", reason="No being used.", is_discontinued=True
)
@total_ordering
class Range(SpanI):  # pragma: no cover
    """Complex object consisting of many spans."""

    def __init__(self, spans=None):
        """Returns a new Range object with data in spans."""
        spans = [] if spans is None else spans
        result = SpansOnly()
        # need to check if we got a single Span, since they define __iter__.
        if isinstance(spans, Span):
            result.append(spans)
        elif hasattr(spans, "spans"):  # probably a single range object?
            result.extend(spans.spans)
        else:
            for s in iterable(spans):
                if hasattr(s, "spans"):
                    result.extend(s.spans)
                else:
                    result.append(s)
        self.spans = result

    def __str__(self):
        """Returns string representation of self."""
        return f"({','.join(map(str, self.spans))})"

    def __len__(self):
        """Returns sum of span lengths.

        NOTE: if spans overlap, will count multiple times. Use reduce() to
        get rid of overlaps.
        """
        return sum(map(len, self.spans))

    def __lt__(self, other):
        """Compares spans of self with indices of other."""
        if hasattr(other, "spans"):
            return self.spans < other.spans
        elif len(self.spans) == 1 and hasattr(other, "start") and hasattr(other, "end"):
            return self.spans[0].start < other.start or self.spans[0].end < other.end
        else:
            return object < other

    def __eq__(self, other):
        """Compares spans of self with indices of other."""
        if hasattr(other, "spans"):
            return self.spans == other.spans
        elif len(self.spans) == 1 and hasattr(other, "start") and hasattr(other, "end"):
            return self.spans[0].start == other.start and self.spans[0].end == other.end
        else:
            return object == other

    def _get_start(self):
        """Finds earliest start of items in self.spans."""
        return min([i.start for i in self.spans])

    start = property(_get_start)

    def _get_end(self):
        """Finds latest end of items in self.spans."""
        return max([i.end for i in self.spans])

    end = property(_get_end)

    def _get_reverse(self):
        """reverse is True if any piece is reversed."""
        for i in self.spans:
            if i.reverse:
                return True
        return False

    reverse = property(_get_reverse)

    def reverses(self):
        """Reverses all spans in self."""
        for i in self.spans:
            i.reverses()

    def __contains__(self, other):
        """Returns True if other completely contained in self.

        other must either be a number or have start and end properties.
        """
        if hasattr(other, "spans"):
            for curr in other.spans:
                found = False
                for i in self.spans:
                    if curr in i:
                        found = True
                        break
                if not found:
                    return False
            return True
        else:
            for i in self.spans:
                if other in i:
                    return True
            return False

    def overlaps(self, other):
        """Returns True if any positions in self are also in other."""
        if hasattr(other, "spans"):
            for i in self.spans:
                for j in other.spans:
                    if i.overlaps(j):
                        return True
        else:
            for i in self.spans:
                if i.overlaps(other):
                    return True
        return False

    def overlaps_extent(self, other):
        """Returns True if any positions in self's extent also in other's."""
        if hasattr(other, "extent"):
            return self.extent.overlaps(other.extent)
        else:
            return self.extent.overlaps(other)

    def sort(self):
        """Sorts the spans in self."""
        self.spans.sort()

    def __iter__(self):
        """Iterates over indices contained in self."""
        return chain(*[iter(i) for i in self.spans])

    def _get_extent(self):
        """Returns Span object representing the extent of self."""
        return Span(self.start, self.end)

    extent = property(_get_extent)

    def simplify(self):
        """Reduces the spans in self in-place to get fewest spans.

        Will not condense spans with opposite directions.

        Will condense adjacent but nonoverlapping spans (e.g. (1,3) and (4,5)).
        """
        forward = []
        reverse = []
        spans = self.spans[:]
        spans.sort()
        for span in spans:
            if span.reverse:
                direction = reverse
            else:
                direction = forward

            found_overlap = False
            for other in direction:
                if (
                    span.overlaps(other)
                    or (span.start == other.end)
                    or (other.start == span.end)
                ):  # handle adjacent spans also
                    other.start = min(span.start, other.start)
                    other.end = max(span.end, other.end)
                    found_overlap = True
                    break
            if not found_overlap:
                direction.append(span)
        self.spans[:] = forward + reverse


@c3warn.deprecated_callable(
    version="2024.3", reason="No being used.", is_discontinued=True
)
class Point(Span):  # pragma: no cover
    """Point is a special case of Span, where start always equals end.

    Note that, as per Python standard, a point is _between_ two elements
    in a sequence. In other words, a point does not contain any elements.
    If you want a single element, use a Span where end = start + 1.

    A Point does have a direction (i.e. a reverse property) to indicate
    where successive items would go if it were expanded.
    """

    def __init__(self, start, reverse=False):
        """Returns new Point object."""
        self.reverse = reverse
        self._start = start

    def _get_start(self):
        """Returns self.start."""
        return self._start

    def _set_start(self, start):
        """Sets self.start and self.end."""
        self._start = start

    start = property(_get_start, _set_start)
    end = start  # start and end are synonyms for the same property


@c3warn.deprecated_callable(
    version="2024.3", reason="No being used.", is_discontinued=True
)
def RangeFromString(string, delimiter=","):  # pragma: no cover
    """Returns Range object from string of the form 1-5,11,20,30-50.

    Ignores whitespace; expects values to be comma-delimited and positive.
    """
    result = Range()
    pairs = list(map(strip, string.split(delimiter)))
    for p in pairs:
        if not p:  # adjacent delimiters?
            continue
        if "-" in p:  # treat as pair
            first, second = p.split("-")
            result.spans.append(Span(int(first), int(second)))
        else:
            result.spans.append(Span(int(p)))
    return result


class MapABC(ABC):
    """base class for genomic map objects"""

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        init_sig = inspect.signature(cls.__init__)
        bargs = init_sig.bind_partial(cls, *args, **kwargs)
        bargs.apply_defaults()
        init_vals = bargs.arguments
        init_vals.pop("self", None)

        obj._serialisable = init_vals
        return obj

    @abstractmethod
    def __len__(self):
        ...

    @abstractmethod
    def __add__(self, other):
        ...

    @classmethod
    @abstractmethod
    def from_rich_dict(cls, data):
        ...

    @abstractmethod
    def gaps(self):
        ...

    @abstractmethod
    def nongap(self):
        ...

    @abstractmethod
    def get_coordinates(self):
        ...

    @abstractmethod
    def nucleic_reversed(self):
        ...

    @abstractmethod
    def to_rich_dict(self):
        ...

    @abstractmethod
    def inverse(self):
        ...


def _spans_from_locations(locations, parent_length) -> Tuple[Union[Span, _LostSpan]]:
    spans = []
    for start, end in locations:
        reverse = start > end
        if max(start, end) < 0 or min(start, end) > parent_length:
            raise RuntimeError(
                f"located outside sequence: {(start, end, parent_length)}"
            )
        if max(start, end) > parent_length and min(start, end) < 0:
            l_diff = min(start, end)
            r_diff = max(start, end) - parent_length
            start, end = (0, parent_length) if start < end else (parent_length, 0)
            spans += [
                LostSpan(abs(l_diff)),
                Span(start, end, reverse=reverse),
                LostSpan(abs(r_diff)),
            ]
        elif min(start, end) < 0:
            diff = min(start, end)
            start = max(start, 0)
            end = max(end, 0)
            spans += [
                LostSpan(abs(diff)),
                Span(start, end, reverse=reverse),
            ]
        elif max(start, end) > parent_length:
            diff = max(start, end) - parent_length
            start = min(start, parent_length)
            end = min(end, parent_length)
            spans += [
                Span(start, end, reverse=reverse),
                LostSpan(abs(diff)),
            ]
        else:
            spans += [Span(start, end, reverse=reverse)]
    return tuple(spans)


T = Union[List[int], Tuple[int]]


@dataclasses.dataclass
class IndelMap(MapABC):
    """store locations of deletions in a Aligned sequence"""

    # todo design notes
    #  I think this object should never try and store "history", i.e. it
    #  should directly relate to the current sequence only. Let that sequence,
    #  which is represented by a SeqView, store history.
    #  Following ths, storing reverse is also a bad idea for this object,
    #  also done by the SeqView only.
    # TODO reverse complement of Alignment -> Aligned -> SeqView, IndelMap
    #  should just do nucleic reverse. I think this is the next task.

    spans: dataclasses.InitVar[Optional[tuple]] = ()
    locations: dataclasses.InitVar[Optional[Sequence[T]]] = None
    termini_unknown: dataclasses.InitVar[bool] = dataclasses.field(default=False)
    parent_length: int = 0
    start: Optional[int] = dataclasses.field(init=False, default=0)
    end: Optional[int] = dataclasses.field(init=False, default=0)
    length: int = dataclasses.field(init=False, default=0)
    _serialisable: dict = dataclasses.field(init=False, repr=False)
    _spans: Tuple[Union[Span, _LostSpan, TerminalPadding]] = dataclasses.field(
        default=(), init=False
    )

    def __post_init__(self, spans, locations, termini_unknown):
        if locations:
            spans = _spans_from_locations(locations, parent_length=self.parent_length)
        elif isinstance(spans, property):
            # This clause is due a known issue with dataclasses.
            # As we have a spans property, the default spans value is
            # ignored, so we have to check for its value being a property
            # and then set the default value here
            spans = ()

        spans = tuple(spans)
        last_not_lost = None
        start, end = None, None
        for i, span in enumerate(spans):
            self.length += len(span)
            if span.lost:
                continue
            elif start is None:
                # this ugly logic because we're using spans!
                start = span.end if span.reverse else span.start

            last_not_lost = i

        if last_not_lost is not None:
            span = spans[last_not_lost]
            end = span.start if span.reverse else span.end

        if start is None:
            start = 0

        if end is None:
            end = self.parent_length

        self.start, self.end = start, end

        if termini_unknown:
            spans = list(spans)
            if spans[0].lost:
                spans[0] = TerminalPadding(spans[0].length)
            if spans[-1].lost:
                spans[-1] = TerminalPadding(spans[-1].length)

        self._spans = tuple(spans)

    def __getitem__(self, slice):
        # A possible shorter map at the same level
        new_map = as_map(slice, len(self), self.__class__)
        new_parts = []
        for span in new_map.spans:
            # we reset tidy start / end to false to avoid changes to
            # default behaviour of tidy on IndelMap (which ignores these
            # anyway)
            old_tidy = span.tidy_start, span.tidy_end
            span.tidy_start, span.tidy_end = False, False
            new_parts.extend(span.remap_with(self))
            span.tidy_start, span.tidy_end = old_tidy
        return self.__class__(spans=new_parts, parent_length=self.parent_length)

    def __len__(self):
        return self.length

    def __add__(self, other):
        if other.parent_length != self.parent_length:
            raise ValueError("Those maps belong to different sequences")
        return self.__class__(
            spans=self._spans + tuple(other.spans), parent_length=self.parent_length
        )

    def __mul__(self, scale):
        # For Protein -> DNA
        new_parts = []
        for span in self._spans:
            new_parts.append(span * scale)
        return self.__class__(spans=new_parts, parent_length=self.parent_length * scale)

    def __repr__(self):
        return repr(self._spans) + f"/{self.parent_length}"

    @property
    def offsets(self):
        return [0] + array([s.length for s in self._spans[:-1]]).cumsum().tolist()

    def gaps(self):
        """The gaps (lost spans) in this map"""
        locations = []
        offset = 0
        for s in self._spans:
            if s.lost:
                locations.append((offset, offset + s.length))
            offset += s.length
        return self.__class__(locations=locations, parent_length=len(self))

    def nongap(self):
        """ungappeed segments in this map"""
        locations = []
        offset = 0
        for s in self._spans:
            if not s.lost:
                locations.append((offset, offset + s.length))
            offset += s.length
        return self.__class__(locations=locations, parent_length=len(self))

    @property
    def spans(self) -> Iterator[Span]:
        """generator of spans"""
        yield from self._spans

    @property
    def complete(self):
        """whether any span represents a gap"""
        return not any(span.lost for span in self._spans)

    @property
    def useful(self):
        return not all(span.lost for span in self._spans)

    def get_coordinates(self):
        """returns span coordinates as [(start, end), ...]"""
        return [(s.start, s.end) for s in self._spans if not s.lost]

    def get_gap_coordinates(self):
        """returns [(gap pos, gap length), ...]"""
        gap_pos = []
        for i, span in enumerate(self._spans):
            if not span.lost:
                continue

            pos = self._spans[i - 1].end if i else 0
            gap_pos.append((pos, len(span)))

        return gap_pos

    def nucleic_reversed(self):
        """Same location on reversed parent"""
        spans = [s.reversed_relative_to(self.parent_length) for s in self._spans]
        return self.__class__(spans=spans, parent_length=self.parent_length)

    def strict_nucleic_reversed(self):
        """map for a sequence that has itself been reversed and complemented

        Notes
        -----
        discards reverse attribute on both spans and self
        """
        spans = []
        parent_length = self.parent_length
        for s in self._spans:
            if not s.lost:
                start = parent_length - s.end
                assert start >= 0
                end = start + s.length
                s = Span(start=start, end=end)
            spans.append(s)

        spans.reverse()
        return self.__class__(spans=spans, parent_length=self.parent_length)

    def to_rich_dict(self):
        """returns dicts for contained spans [dict(), ..]"""
        spans = [s.to_rich_dict() for s in self._spans]
        # exclude spans from deep copy since being overwritten
        data = copy.deepcopy(
            {k: v for k, v in self._serialisable.items() if k != "spans"}
        )
        data.pop("locations")
        data["spans"] = spans
        data["type"] = get_object_provenance(self)
        data["version"] = __version__
        return data

    @classmethod
    def from_rich_dict(cls, map_element):
        from cogent3.util.deserialise import _get_class

        map_element.pop("version", None)
        type_ = map_element.pop("type")
        assert _get_class(type_) == cls
        spans = []
        for element in map_element["spans"]:
            element.pop("version", None)
            klass = _get_class(element.pop("type"))
            instance = klass(**element)
            spans.append(instance)

        map_element["spans"] = spans
        return cls(**map_element)

    def with_termini_unknown(self):
        """returns new instance with terminal gaps indicated as unknown"""
        return self.__class__(
            spans=self._spans[:],
            parent_length=self.parent_length,
            termini_unknown=True,
        )

    def inverse(self):
        """returns instance with coordinates updated for aligned, unaligned"""
        # is this only required for parse_out_gaps?
        # NO also used in cogent3.align code

        # can't work if there are overlaps in the map
        # tidy ends don't survive inversion
        if self.parent_length is None:
            raise ValueError("Uninvertable. parent length not known")

        cum_posn = 0
        temp = []
        for span in self._spans:
            if not span.lost:
                if span.reverse:
                    temp.append(
                        (span.start, span.end, cum_posn + span.length, cum_posn)
                    )
                else:
                    temp.append(
                        (span.start, span.end, cum_posn, cum_posn + span.length)
                    )
            cum_posn += span.length

        temp.sort()
        new_spans = []
        last_start = 0
        for start, end, cum_start, cum_end in temp:
            if start > last_start:
                new_spans.append(LostSpan(start - last_start))
            elif start < last_start:
                raise ValueError(f"Uninvertable. Overlap: {start} < {last_start}")

            # we force tidy_<start/end> to be same as self, attribute has no meaning
            # for IndelMap, but retained for compatability for now
            new_spans.append(
                Span(
                    cum_start,
                    cum_end,
                    reverse=cum_start > cum_end,
                )
            )
            last_start = end

        if self.parent_length > last_start:
            new_spans.append(LostSpan(self.parent_length - last_start))

        return self.__class__(spans=new_spans, parent_length=len(self))

    T = Union[ndarray, int]

    def absolute_position(self, rel_pos: T) -> T:
        """converts rel_pos into an absolute position

        Raises
        ------
        raises ValueError if rel_pos < 0
        """
        check = array([rel_pos], dtype=int) if isinstance(rel_pos, int) else rel_pos
        if check.min() < 0:
            raise ValueError(f"must positive, not {rel_pos=}")

        if len(self) == self.parent_length:
            # handle case of reversed here?
            return rel_pos

        return self.start + rel_pos

    def relative_position(self, abs_pos: T) -> T:
        """converts abs_pos into an relative position

        Raises
        ------
        raises ValueError if abs_pos < 0
        """
        check = array([abs_pos], dtype=int) if isinstance(abs_pos, int) else abs_pos
        if check.min() < 0:
            raise ValueError(f"must positive, not {abs_pos=}")
        return abs_pos - self.start

    def get_covering_span(self):
        span = (self.end, self.start)
        return self.__class__(locations=[span], parent_length=self.parent_length)

    def zeroed(self):
        """returns a new instance with the first span starting at 0

        Note
        ----

        Useful when an annotated sequence is sliced, but the connection to
        the original parent is being deliberately broken as in the
        Sequence.deepcopy(sliced=True) case.
        """
        # todo is this really required, best ifn we can rely on SeqView to
        #  store all relationship to underlying sequence
        min_val = min(self.start, self.end)
        spans = []
        for span in self._spans:
            if span.lost:
                spans.append(span)
                continue
            kwargs = span.to_rich_dict()
            del kwargs["version"]
            del kwargs["type"]
            kwargs["start"] = kwargs["start"] - min_val
            kwargs["end"] = kwargs["end"] - min_val
            spans.append(Span(**kwargs))

        kwargs = self.to_rich_dict()
        del kwargs["version"]
        del kwargs["type"]
        kwargs["spans"] = spans
        kwargs["parent_length"] = abs(self.start - self.end)
        return self.__class__(**kwargs)

    def to_feature_map(self):
        """returns a Map type, suited to Features"""
        spans = []
        for span in self._spans:
            if span.lost:
                spans.append(span)
                continue
            kwargs = span.to_rich_dict()
            del kwargs["version"]
            del kwargs["type"]
            spans.append(Span(**kwargs))

        kwargs = self.to_rich_dict()
        del kwargs["version"]
        del kwargs["type"]
        del kwargs["termini_unknown"]  # not used for FeatureMap
        kwargs["spans"] = spans
        return FeatureMap(**kwargs)


class FeatureMap:
    """A map holds a list of spans."""

    def __init__(
        self,
        locations=None,
        spans=None,
        parent_length=None,
    ):
        assert parent_length is not None
        d = locals()
        exclude = ("self", "__class__", "__slots__")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}

        if spans is None:
            spans = _spans_from_locations(locations, parent_length=parent_length)

        self.offsets = []
        self.useful = False
        self.complete = True
        self.reverse = None
        posn = 0
        for span in spans:
            self.offsets.append(posn)
            posn += span.length
            if span.lost:
                self.complete = False
            elif not self.useful:
                self.useful = True
                self.start, self.end = span.start, span.end
                self.reverse = span.reverse
            else:
                self.start = min(self.start, span.start)
                self.end = max(self.end, span.end)
                if self.reverse is not None and (span.reverse != self.reverse):
                    self.reverse = None

        self.spans = tuple(spans)
        self.length = posn
        self.parent_length = parent_length

    def __len__(self):
        return self.length

    def __repr__(self):
        return repr(list(self.spans)) + f"/{self.parent_length}"

    def __getitem__(self, new_map):
        # A possible shorter map at the same level
        new_map = as_map(new_map, len(self), self.__class__)
        new_parts = []
        for span in new_map.spans:
            new_parts.extend(span.remap_with(self))
        return self.__class__(spans=new_parts, parent_length=self.parent_length)

    def __mul__(self, scale):
        new_parts = [span * scale for span in self.spans]
        return self.__class__(spans=new_parts, parent_length=self.parent_length * scale)

    def __div__(self, scale):
        new_parts = [span / scale for span in self.spans]
        return self.__class__(
            spans=new_parts, parent_length=self.parent_length // scale
        )

    def __add__(self, other):
        if other.parent_length != self.parent_length:
            raise ValueError("Those maps belong to different sequences")
        return self.__class__(
            spans=self.spans + other.spans, parent_length=self.parent_length
        )

    def _with_termini_unknown(self):
        return self.__class__(
            self,
            spans=self.spans[:],
            parent_length=self.parent_length,
            termini_unknown=True,
        )

    def get_covering_span(self):
        if self.reverse:
            span = (self.end, self.start)
        else:
            span = (self.start, self.end)
        return self.__class__([span], parent_length=self.parent_length)

    def covered(self):
        """>>> Map([(10,20), (15, 25), (80, 90)]).covered().spans
        [Span(10,25), Span(80, 90)]"""

        delta = {}
        for span in self.spans:
            if span.lost:
                continue
            delta[span.start] = delta.get(span.start, 0) + 1
            delta[span.end] = delta.get(span.end, 0) - 1
        positions = sorted(delta.keys())
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
        return self.__class__(locations=result, parent_length=self.parent_length)

    def reversed(self):
        """Reversed location on same parent"""
        spans = [s.reversed() for s in self.spans]
        spans.reverse()
        return self.__class__(spans=spans, parent_length=self.parent_length)

    def nucleic_reversed(self):
        """Same location on reversed parent"""
        spans = [s.reversed_relative_to(self.parent_length) for s in self.spans]
        return self.__class__(spans=spans, parent_length=self.parent_length)

    def get_gap_coordinates(self):
        """returns [(gap pos, gap length), ...]"""
        gap_pos = []
        for i, span in enumerate(self.spans):
            if not span.lost:
                continue

            pos = self.spans[i - 1].end if i else 0
            gap_pos.append((pos, len(span)))

        return gap_pos

    def gaps(self):
        """The gaps (lost spans) in this map"""
        locations = []
        offset = 0
        for s in self.spans:
            if s.lost:
                locations.append((offset, offset + s.length))
            offset += s.length
        return self.__class__(locations, parent_length=len(self))

    def shadow(self):
        """The 'negative' map of the spans not included in this map"""
        return self.inverse().gaps()

    def nongap(self):
        locations = []
        offset = 0
        for s in self.spans:
            if not s.lost:
                locations.append((offset, offset + s.length))
            offset += s.length
        return self.__class__(locations, parent_length=len(self))

    def without_gaps(self):
        return self.__class__(
            spans=[s for s in self.spans if not s.lost],
            parent_length=self.parent_length,
        )

    def inverse(self):
        """returns instance with coordinates updated for aligned, unaligned"""
        # is this only required for parse_out_gaps?
        # NO also used in cogent3.align code

        # can't work if there are overlaps in the map
        # tidy ends don't survive inversion
        if self.parent_length is None:
            raise ValueError("Uninvertable. parent length not known")

        cum_posn = 0
        temp = []
        for span in self.spans:
            if not span.lost:
                if span.reverse:
                    temp.append(
                        (span.start, span.end, cum_posn + span.length, cum_posn)
                    )
                else:
                    temp.append(
                        (span.start, span.end, cum_posn, cum_posn + span.length)
                    )
            cum_posn += span.length

        temp.sort()
        new_spans = []
        last_start = 0
        for start, end, cum_start, cum_end in temp:
            if start > last_start:
                new_spans.append(LostSpan(start - last_start))
            elif start < last_start:
                raise ValueError(f"Uninvertable. Overlap: {start} < {last_start}")

            # we force tidy_<start/end> to be same as self, attribute has no meaning
            # for IndelMap, but retained for compatability for now
            new_spans.append(
                Span(
                    cum_start,
                    cum_end,
                    reverse=cum_start > cum_end,
                )
            )
            last_start = end

        if self.parent_length > last_start:
            new_spans.append(LostSpan(self.parent_length - last_start))

        return self.__class__(spans=new_spans, parent_length=len(self))

    def _inverse(self):
        # can't work if there are overlaps in the map
        # tidy ends don't survive inversion
        if self.parent_length is None:
            raise ValueError("Uninvertable. parent length not known")
        posn = 0
        temp = []
        for span in self.spans:
            if not span.lost:
                if span.reverse:
                    temp.append((span.start, span.end, posn + span.length, posn))
                else:
                    temp.append((span.start, span.end, posn, posn + span.length))
            posn += span.length

        temp.sort()
        new_spans = []
        last_hi = 0
        for lo, hi, start, end in temp:
            if lo > last_hi:
                new_spans.append(LostSpan(lo - last_hi))
            elif lo < last_hi:
                raise ValueError(f"Uninvertable. Overlap: {lo} < {last_hi}")
            new_spans.append(Span(start, end, reverse=start > end))
            last_hi = hi
        if self.parent_length > last_hi:
            new_spans.append(LostSpan(self.parent_length - last_hi))
        return self.__class__(spans=new_spans, parent_length=len(self))

    def get_coordinates(self):
        """returns span coordinates as [(v1, v2), ...]

        v1/v2 are (start, end) unless the map is reversed, in which case it will
        be (end, start)"""

        if self.reverse:
            order_func = lambda x: (max(x), min(x))
        else:
            order_func = lambda x: x

        coords = list(
            map(order_func, [(s.start, s.end) for s in self.spans if not s.lost])
        )

        return coords

    def to_rich_dict(self):
        """returns dicts for contained spans [dict(), ..]"""
        spans = [s.to_rich_dict() for s in self.spans]
        data = copy.deepcopy(self._serialisable)
        data.pop("locations")
        data["spans"] = spans
        data["type"] = get_object_provenance(self)
        data["version"] = __version__
        return data

    def zeroed(self):
        """returns a new instance with the first span starting at 0

        Note
        ----

        Useful when an Annotatable object is sliced, but the connection to
        the original parent is being deliberately broken as in the
        Sequence.deepcopy(sliced=True) case.
        """
        # todo there's probably a more efficient way to do this
        # create the new instance
        from cogent3.util.deserialise import deserialise_map_spans

        data = self.to_rich_dict()
        zeroed = deserialise_map_spans(data)
        zeroed.parent_length = len(zeroed.get_covering_span())
        shift = min(zeroed.start, zeroed.end)
        new_end = 0
        for span in zeroed.spans:
            if span.lost:
                continue
            span.start -= shift
            span.end -= shift
            new_end = max(new_end, span.end)

        zeroed.start = 0
        zeroed.end = new_end

        return zeroed

    T = Union[ndarray, int]

    def absolute_position(self, rel_pos: T) -> T:
        """converts rel_pos into an absolute position

        Raises
        ------
        raises ValueError if rel_pos < 0
        """
        check = array([rel_pos], dtype=int) if isinstance(rel_pos, int) else rel_pos
        if check.min() < 0:
            raise ValueError(f"must positive, not {rel_pos=}")

        if len(self) == self.parent_length:
            # handle case of reversed here?
            return rel_pos

        return self.start + rel_pos

    def relative_position(self, abs_pos: T) -> T:
        """converts abs_pos into an relative position

        Raises
        ------
        raises ValueError if abs_pos < 0
        """
        check = array([abs_pos], dtype=int) if isinstance(abs_pos, int) else abs_pos
        if check.min() < 0:
            raise ValueError(f"must positive, not {abs_pos=}")
        return abs_pos - self.start


def gap_coords_to_map(gaps_lengths: dict, seq_length: int) -> FeatureMap:
    """
    Parameters
    ----------
    gaps_lengths
        {gap insertion pos: gap length, ...}
    seq_length : int
        length of unaligned sequence

    Returns
    -------
    Map
    """

    if not gaps_lengths:
        return FeatureMap([(0, seq_length)], parent_length=seq_length)

    spans = []
    last = pos = 0
    for pos in sorted(gaps_lengths):
        if pos > seq_length:
            raise ValueError(
                f"cannot have gap at position {pos} beyond seq_length= {seq_length}"
            )

        gap = LostSpan(length=gaps_lengths[pos])
        spans.extend([gap] if pos == 0 else [Span(last, pos), gap])
        last = pos

    if pos < seq_length:
        spans.append(Span(last, seq_length))

    return FeatureMap(spans=spans, parent_length=seq_length)
