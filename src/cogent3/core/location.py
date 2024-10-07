import copy
import dataclasses
import functools
import inspect
import json
from abc import ABC, abstractmethod
from bisect import bisect_left, bisect_right
from functools import total_ordering
from typing import Any, Iterator, Optional, Sequence, Union

import numpy
from numpy.typing import NDArray

from cogent3._version import __version__
from cogent3.util import warning as c3warn
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import get_object_provenance

strip = str.strip

_DEFAULT_GAP_DTYPE = numpy.int32

OptInt = Optional[int]


def _norm_index(i, length, default):
    """For converting s[:3] to s[0:3], s[-1] to s[len(s)-1] and s[0:lots] to s[0:len(s)]"""
    if i is None:
        i = default
    elif i < 0:
        i += length
    return min(max(i, 0), length)


def as_map(slice, length, cls):
    """Take anything that might be used as a subscript: Integer, Slice,
    or MapABC, and return cls."""

    if isinstance(slice, (list, tuple)):
        spans = []
        for i in slice:
            spans.extend(as_map(i, length, cls).spans)
        return cls(spans=spans, parent_length=length)
    elif isinstance(slice, (FeatureMap, IndelMap, Map)):
        return slice
    else:
        lo, hi, step = _norm_slice(slice, length)
        assert (step or 1) == 1
        # since we disallow step, a reverse slice means an empty series
        locations = [] if lo > hi else [(lo, hi)]
        # IndelMap has this class method
        cls = getattr(cls, "from_locations", cls)
        return cls(locations=locations, parent_length=length)


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

    def __truediv__(self, scale: int):
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

    def __truediv__(self, scale):
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


IntTypes = Union[int, numpy.int32, numpy.int64]
IntArrayTypes = NDArray[int]
SpanTypes = Union[Span, _LostSpan]
SeqSpanTypes = Sequence[SpanTypes]
SeqCoordTypes = Sequence[Sequence[IntTypes]]


class Map:  # pragma: no cover
    """A map holds a list of spans."""

    @c3warn.deprecated_callable(
        version="2024.6",
        reason="Replaced by IndelMap and FeatureMap",
        is_discontinued=True,
    )
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
        data.pop("locations", None)
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

    T = Union[numpy.ndarray, int]

    def absolute_position(self, rel_pos: T) -> T:
        """converts rel_pos into an absolute position

        Raises
        ------
        raises ValueError if rel_pos < 0
        """
        check = (
            numpy.array([rel_pos], dtype=int) if isinstance(rel_pos, int) else rel_pos
        )
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
        check = (
            numpy.array([abs_pos], dtype=int) if isinstance(abs_pos, int) else abs_pos
        )
        if check.min() < 0:
            raise ValueError(f"must positive, not {abs_pos=}")
        return abs_pos - self.start


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
    def __len__(self) -> int: ...

    @abstractmethod
    def __add__(self, other: "MapABC") -> "MapABC": ...

    @abstractmethod
    def nongap(self) -> Iterator[SpanTypes]: ...

    @abstractmethod
    def get_coordinates(self) -> SeqCoordTypes: ...

    @abstractmethod
    def nucleic_reversed(self) -> "MapABC": ...

    @abstractmethod
    def to_rich_dict(self) -> dict[str, Any]: ...

    @classmethod
    def from_locations(
        cls, locations: SeqCoordTypes, parent_length: int, **kwargs
    ) -> "MapABC":
        if len(locations):
            spans = _spans_from_locations(locations, parent_length=parent_length)
        else:
            spans = ()

        return cls.from_spans(spans=spans, parent_length=parent_length, **kwargs)

    @classmethod
    @abstractmethod
    def from_spans(
        cls, spans: SeqSpanTypes, parent_length: int, **kwargs
    ) -> "MapABC": ...


def _spans_from_locations(locations: SeqCoordTypes, parent_length: int) -> SeqSpanTypes:
    if not len(locations):
        # using len() because locations can be a numpy array
        return ()

    if locations[0][0] > locations[-1][1]:
        raise ValueError(f"locations must be ordered smallest-> largest {locations}")

    spans = []
    for start, end in locations:
        if start > end or min(start, end) < 0:
            raise ValueError("locations must be ordered smallest-> largest and >= 0")
        if start > parent_length:
            raise RuntimeError(
                f"located outside sequence: {(start, end, parent_length)}"
            )
        if end > parent_length:
            diff = end - parent_length
            end = min(end, parent_length)
            spans += [
                Span(start, end),
                LostSpan(abs(diff)),
            ]
        else:
            spans += [Span(start, end)]

    return tuple(spans)


def spans_to_gap_coords(
    indel_spans: SeqSpanTypes, dtype: IntTypes = _DEFAULT_GAP_DTYPE
) -> tuple[IntArrayTypes, IntArrayTypes]:
    """returns coordinates of sequence gaps

    Parameters
    ----------
    indel_spans
        sequence of Span/LostSpan instances
    dtype
        numpy type, default to _DEFAULT_GAP_DTYPE

    Returns
    -------
    numpy.array([gap pos,...]), numpy.array([cum gap length,...]),
    """
    gap_pos = []
    cum_lengths = []
    cum_length = 0
    for i, span in enumerate(indel_spans):
        if not span.lost:
            continue

        pos = indel_spans[i - 1].end if i else 0
        cum_length += len(span)
        gap_pos.append(pos)
        cum_lengths.append(cum_length)

    return numpy.array(gap_pos, dtype=dtype), numpy.array(cum_lengths, dtype=dtype)


def _gap_spans(
    gap_pos: IntArrayTypes, cum_gap_lengths: IntArrayTypes
) -> tuple[IntArrayTypes, IntArrayTypes]:
    """returns 1D arrays in alignment coordinates of
    gap start, gap stop"""
    if not len(gap_pos):
        r = numpy.array([], dtype=_DEFAULT_GAP_DTYPE)
        return r, r

    ends = gap_pos + cum_gap_lengths
    starts = gap_pos.copy()
    starts[1:] += cum_gap_lengths[:-1]

    return starts, ends


def _update_lengths(
    result_pos: IntArrayTypes,
    result_lengths: IntArrayTypes,
    gap_pos: IntArrayTypes,
    gap_lengths: IntArrayTypes,
) -> None:
    """modifies result_lengths in place with gap_lengths
    where elements in gap_pos occur in result_pos

    Notes
    -----
    result_pos is a superset of gap_pos
    """
    _, result_indices, other_indices = numpy.intersect1d(
        result_pos, gap_pos, assume_unique=True, return_indices=True
    )
    result_lengths[result_indices] += gap_lengths[other_indices]


def _step_adjustment(gap_start: int, start: int, step: int) -> int:
    """adjustment to a start, given a step size. The adjustment determines how
    many more steps we need to take from the beginning of a gap/ungapped segment
    to reach a position that is a multiple of the step"""
    overstep = (gap_start - start) % step
    return step - overstep if overstep else 0


def _step_adjusted_length(start: int, end: int, adj: int, step: int) -> int:
    """returns the adjusted length of a segment given a step size
    and its start and stop index. The adjustment is the number of steps
    needed to reach the next multiple of the step size"""
    adjusted_length = -((end - start - adj) // -step)
    return max(adjusted_length, 0)


def _input_vals_pos_step(
    seqlen: int, start: OptInt, stop: OptInt, step: int
) -> tuple[int, int, int]:
    """returns standardised start, stop, step values for positive step slicing.

    Notes
    -----
    The input values for start and stop can be +ve, -ve or None. The returned
    start and stop are strictly +ve and within the sequence length, and if not
    provided default to 0 and seqlen respectively.

    If the slice results in a zero-length slice, the returned start, stop, step
    will be 0, 0, 1 respectively.
    """
    start = start or 0
    if start > 0 and start >= seqlen:
        # start beyond seq is an empty slice
        return 0, 0, 1

    stop = seqlen if stop is None else stop
    if stop < 0 and abs(stop) >= seqlen:
        # finished slice before we started seq!
        return 0, 0, 1

    start = max(seqlen + start, 0) if start < 0 else start

    if stop > 0:
        stop = min(seqlen, stop)
    elif stop < 0:
        stop += seqlen

    if start >= stop:
        start = stop = 0
        step = 1

    return start, stop, step


def _input_vals_neg_step(
    seqlen: int, start: OptInt, stop: OptInt, step: int
) -> tuple[int, int, int]:
    """returns standardised start, stop, step values for negative step slicing.

    Notes
    -----
    The input values for start and stop can be positive, negative or None. The
    returned start and stop are strictly negative and correspond to indices
    within the sequence length. If not provided, start defaults to -1 and stop
    to -seqlen-1.

    If the slice results in a zero-length slice, the returned start, stop, step
    will be 0, 0, 1 respectively.
    """
    # Note how Python reverse slicing works
    # we need to make sure the start and stop are both
    # negative, for example "abcd"[-1:-5:-1] returns "dcba"
    # returning the complete string in reverse order is only
    # possible with negative indexing
    if start is None or start >= seqlen:  # set default
        start = -1  # Done
    elif start >= 0:  # convert to -ve index
        start = start - seqlen
    elif start < -seqlen:  # start is bounded by len(seq)
        return 0, 0, 1

    if stop is None:  # set default
        stop = -seqlen - 1
    elif stop >= 0:
        stop -= seqlen

    stop = max(stop, -seqlen - 1)  # stop should always be <= len(seq)

    # checking for zero-length slice
    return (0, 0, 1) if start < stop else (start, stop, step)


@dataclasses.dataclass
class IndelMap(MapABC):
    """store locations of deletions in a Aligned sequence

    Parameters
    ----------
    gap_pos
        start positions of gap segments in sequence (ungapped) coordinates
    cum_gap_lengths
        cumulative gap lengths per segment
    gap_lengths
        length of each gap segment
    termini_unknown
        if ``True``, returns new instance with terminal gaps indicated as
        unknown character '?'
    parent_length
        length of parent sequence (i.e. aligned sequence with gaps)
    """

    # gap data is gap positions, gap lengths on input, stored
    gap_pos: IntArrayTypes
    cum_gap_lengths: Optional[IntArrayTypes] = None
    gap_lengths: dataclasses.InitVar[Optional[IntArrayTypes]] = None
    termini_unknown: bool = False
    parent_length: int = 0
    _serialisable: dict = dataclasses.field(init=False, repr=False)
    num_gaps: int = dataclasses.field(init=False, repr=False, default=0)

    def __post_init__(self, gap_lengths: IntArrayTypes):
        assert gap_lengths is None or self.cum_gap_lengths is None
        if gap_lengths is not None:
            self.cum_gap_lengths = gap_lengths.cumsum()

        if len(self.gap_pos) != len(self.cum_gap_lengths):
            raise ValueError(
                f"length of gap_ pos {len(self.gap_pos)} != "
                f"length of gap lengths {len(self.cum_gap_lengths)}"
            )

        self.num_gaps = self.gap_pos.shape[0]
        if self.num_gaps and self.gap_pos[-1] > self.parent_length:
            raise ValueError(
                f"gap position {self.gap_pos[-1]} outside parent_length {self.parent_length}"
            )

        # make gap array immutable
        self.gap_pos.flags.writeable = False
        self.cum_gap_lengths.flags.writeable = False
        self._serialisable.pop("gap_lengths", None)

    @classmethod
    def from_spans(
        cls, spans: SeqSpanTypes, parent_length: int, termini_unknown: bool = False
    ) -> "IndelMap":
        gap_pos, cum_lengths = spans_to_gap_coords(spans)
        return cls(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_lengths,
            parent_length=parent_length,
            termini_unknown=termini_unknown,
        )

    @classmethod
    def from_aligned_segments(
        cls, locations: SeqCoordTypes, aligned_length: int
    ) -> "IndelMap":
        """
        converts coordinates from aligned segments into IndelMap for ungapped sequence

        Parameters
        ----------
        locations
            sequence of ungapped segment in alignment coordinates
        aligned_length
            length of the alignment
        """
        locations = list(locations)
        if not locations or (
            len(locations) == 1
            and locations[0][0] == 0
            and locations[0][-1] == aligned_length
        ):
            empty = numpy.array([], dtype=_DEFAULT_GAP_DTYPE)
            return cls(
                gap_pos=empty,
                cum_gap_lengths=empty.copy(),
                parent_length=aligned_length,
            )

        if locations[0][0] != 0:
            # starts with a gap
            locations = [(0, 0)] + locations
        if locations[-1][1] < aligned_length:
            # ends with a gap
            locations += [(aligned_length, aligned_length)]

        locations = numpy.array(locations, dtype=_DEFAULT_GAP_DTYPE).flatten()[1:-1]
        gap_coords = locations.reshape((locations.shape[0] // 2, 2))
        gap_ends, gap_starts = gap_coords[:, ::-1].T
        gap_lengths = gap_ends - gap_starts
        cum_lens = gap_lengths.cumsum()
        # convert to sequence coords
        gap_pos = gap_starts.copy()
        gap_pos[1:] -= cum_lens[:-1]
        seq_length = aligned_length - cum_lens[-1]

        return cls(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_lens,
            parent_length=seq_length,
        )

    # NOTE: cannot use string type hints with singledispatchmethod
    @functools.singledispatchmethod
    def __getitem__(self, item):
        raise NotImplementedError(f"cannot slice using {type(item)}")

    @__getitem__.register
    def _(self, item: int):
        return self[item : item + 1]

    @__getitem__.register
    def _(self, item: slice):
        # we're assuming that this gap object is associated with a sequence
        # that will also be sliced. Hence, we need to shift the gap insertion
        # positions relative to this newly sliced sequence.

        # refactor: simplify

        step = 1 if item.step is None else item.step
        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
        start, stop, step = func(len(self), item.start, item.stop, step)

        if step < 0:
            start = -start - 1
            stop = -stop - 1
            return self.nucleic_reversed()[start:stop:-step]

        if self._check_no_gap_slice(start, stop):
            return self._zero_imap(start, stop, step)

        gap_starts, gap_ends = _gap_spans(self.gap_pos, self.cum_gap_lengths)

        # determine the gap index where the slice starts
        l = numpy.searchsorted(gap_ends, start + 1, side="left")

        # check if entire slice is within a single gap
        if gap_starts[l] <= start < gap_ends[l] and stop <= gap_ends[l]:
            return self._single_imap(start, stop, step)

        cum_seq_length = 0
        adj_gaps = []
        lengths = self.get_gap_lengths()

        # determine the beginning of the slice
        if gap_starts[l] <= start < gap_ends[l]:
            # start is within a gap
            adj = start - gap_starts[l]
            # we want "ceiling division" to determine how many steps are
            # included which we can achieve with upside down "floor division"
            # e.g., a gap of 3 with a step of 2 should be 2 steps
            adj_gap_len = -((lengths[l] - adj) // -step)
            adj_gaps.append([0, adj_gap_len])
        else:
            # start is within a ungapped segment
            if stop <= gap_starts[l]:
                # the stop is within the same ungapped segment
                return self._zero_imap(start, stop, step)

            # determine the preceding seq length
            cum_seq_length += _step_adjusted_length(
                start=start, end=gap_starts[l], adj=0, step=step
            )
            # adj is how many more steps we need to take from the beginning of
            # the gap to reach a position which is a multiple of the step
            adj = _step_adjustment(gap_starts[l], start, step)
            # check gap is not stepped over entirely
            if adj < lengths[l]:
                adj_gap_len = -((lengths[l] - adj) // -step)
                adj_gaps.append([cum_seq_length, adj_gap_len])

        # start search for rhs index
        r = numpy.searchsorted(gap_ends[l:], stop, side="right") + l

        # iterate through the gaps that are fully within the slice
        for j in range(l + 1, r):
            # determine how long the preceding ungapped segment was
            adj = _step_adjustment(gap_ends[j - 1], start, step)
            adj_seq_len = _step_adjusted_length(
                (gap_ends[j - 1]), gap_starts[j], adj, step
            )
            cum_seq_length += adj_seq_len

            # calculate the adjustment to the gap
            adj = _step_adjustment(gap_starts[j], start, step)
            if adj < lengths[j]:
                # if adj > length, then the gap is "stepped over" by the stride
                adj_gap_len = -((lengths[j] - adj) // -step)

                if adj_gaps and adj_gaps[-1][0] == cum_seq_length:
                    # previous gap is contiguous with this one
                    adj_gaps[-1][1] += adj_gap_len
                else:
                    adj_gaps.append([cum_seq_length, adj_gap_len])

        # determine the end of the slice
        if l == r:
            # the stop was inside the first gap
            adj = _step_adjustment(gap_starts[r], start, step)

            if adj_gaps:
                # start was in the gap
                adj_gaps[-1][1] = _step_adjusted_length(gap_starts[r], stop, adj, step)
            else:
                # start was in the seq before the gap
                adj_gaps.append(
                    [
                        cum_seq_length,
                        _step_adjusted_length(gap_starts[r], stop, adj, step),
                    ]
                )

        elif stop >= gap_ends[-1] or stop < gap_starts[r]:
            # stop is within an ungapped segment
            adj = _step_adjustment(gap_ends[r - 1], start, step)
            cum_seq_length += _step_adjusted_length(gap_ends[r - 1], stop, adj, step)
        else:
            # stop is within a gap
            # determine length of previous ungapped segment
            adj = _step_adjustment(gap_ends[r - 1], start, step)
            cum_seq_length += _step_adjusted_length(
                gap_ends[r - 1], gap_starts[r], adj, step
            )

            # determine length of the gap when considering the stop
            adj = _step_adjustment(gap_starts[r], start, step)
            adj_gap_len = _step_adjusted_length(gap_starts[r], stop, adj, step)

            # check that we do not step over the gap entirely
            if adj_gap_len > 0:
                if adj_gaps and adj_gaps[-1][0] == cum_seq_length:
                    # the previous gap is contiguous with this one
                    adj_gaps[-1][1] += adj_gap_len
                else:
                    adj_gaps.append([cum_seq_length, adj_gap_len])

        pos_result = numpy.array([x for x, _ in adj_gaps])
        lengths = numpy.array([y for _, y in adj_gaps])
        parent_length = cum_seq_length

        return self.__class__(
            gap_pos=pos_result,
            gap_lengths=lengths,
            parent_length=parent_length,
        )

    def _zero_imap(self, start: int, stop: int, step: int) -> "IndelMap":
        """returns a new IndelMap with zero gaps"""
        zero_array = numpy.array([], dtype=_DEFAULT_GAP_DTYPE)
        return self.__class__(
            gap_pos=zero_array.copy(),
            cum_gap_lengths=zero_array.copy(),
            parent_length=_step_adjusted_length(start, stop, 0, step),
        )

    def _single_imap(self, start: int, stop: int, step: int) -> "IndelMap":
        """returns a new IndelMap with a single gap"""
        gap_pos = numpy.array([0], dtype=_DEFAULT_GAP_DTYPE)
        cum_gap_length = numpy.array([_step_adjusted_length(start, stop, 0, step)])
        return self.__class__(
            gap_pos=gap_pos, cum_gap_lengths=cum_gap_length, parent_length=0
        )

    def _check_no_gap_slice(self, start: int, stop: int) -> bool:
        """check for easy cases where slicing an indel maps results in no gaps.

        Returns
            True if no gaps are present in the slice

        Notes
        -----
        this method assumes positive indexing and a positive step
        (i.e. start, stop, step > 0)
        """
        # when slicing an indel maps, we can have four easy cases:
        # 1 - invalid slice
        if start == stop:
            return True

        # 2 - no gaps
        if not self.num_gaps:
            return True

        first_gap = self.gap_pos[0]  # start of first gap
        last_gap = self.gap_pos[-1] + self.cum_gap_lengths[-1]  # end of last gap

        # 3 - slice before first gap; 4 - after last gap (for positive step)
        if stop <= first_gap or start >= last_gap:
            return True

    def get_align_index(self, seq_index: int, slice_stop: bool = False) -> int:
        """convert a sequence index into an alignment index

        Parameters
        ----------
        seq_index
            coordinate on the sequence, must be < parent_length
        slice_stop
            set to True if the index is to be the end of an alignment slice.
            In that case, and if seq_index is in gap_pos then it returns
            the first alignment index of the gap run.
        """
        cum_lengths = self.cum_gap_lengths
        gap_pos = self.gap_pos
        # NOTE I explicitly cast all returned values to python int's due to
        # need for json serialisation, which does not support numpy int classes
        if seq_index < 0:
            seq_index += self.parent_length

        if seq_index < 0:
            raise IndexError(f"{seq_index} negative seq_index beyond limit ")

        if not self.num_gaps or seq_index < gap_pos[0]:
            return int(seq_index)

        # if stop_index, check if the seq_index corresponds to a gap position
        if slice_stop and (match := seq_index == gap_pos).any():
            # if so, we return the alignment coord for the first gap position
            (idx,) = numpy.where(match)[0]
            if idx:
                gap_len = cum_lengths[idx] - cum_lengths[idx - 1]
            else:
                gap_len = cum_lengths[idx]
            gap_end = gap_pos[idx] + cum_lengths[idx]
            return int(gap_end - gap_len)

        if seq_index >= gap_pos[-1]:
            return int(seq_index + cum_lengths[-1])

        # find gap position before seq_index
        index = numpy.searchsorted(gap_pos, seq_index, side="left")
        if seq_index < gap_pos[index]:
            gap_lengths = cum_lengths[index - 1] if index else 0
        else:
            gap_lengths = cum_lengths[index]

        return int(seq_index + gap_lengths)

    def get_seq_index(self, align_index: int) -> int:
        """converts alignment index to sequence index"""
        # NOTE I explicitly cast all returned values to python int's due to
        # need for json serialisation, which does not support numpy int classes
        if align_index < 0:
            align_index = len(self) + align_index
        if align_index < 0:
            raise IndexError(f"{align_index} align_index beyond limit")

        if not self.num_gaps or align_index < self.gap_pos[0]:
            return align_index

        # these are alignment indices for gaps
        cum_lengths = self.cum_gap_lengths
        gap_starts, gap_ends = _gap_spans(self.gap_pos, cum_lengths)
        if align_index >= gap_ends[-1]:
            return int(align_index - cum_lengths[-1])

        index = numpy.searchsorted(gap_ends, align_index, side="left")
        if align_index < gap_starts[index]:
            # before the gap at index
            return int(align_index - cum_lengths[index - 1])

        if align_index == gap_ends[index]:
            # after the gap at index
            return int(align_index - cum_lengths[index])

        if gap_starts[index] <= align_index < gap_ends[index]:
            # within the gap at index
            # so the gap insertion position is the sequence position
            return int(self.gap_pos[index])

    def __len__(self) -> int:
        length_gaps = self.cum_gap_lengths[-1] if self.num_gaps else 0
        return int(self.parent_length + length_gaps)

    def __add__(self, other: "IndelMap") -> "IndelMap":
        """designed to support concatenation of two aligned sequences"""
        gap_pos = self.gap_pos.tolist() + (self.parent_length + other.gap_pos).tolist()
        gap_pos = numpy.array(gap_pos, dtype=_DEFAULT_GAP_DTYPE)

        cum_length = self.cum_gap_lengths[-1] if self.num_gaps else 0
        cum_gap_lengths = (
            self.cum_gap_lengths.tolist()
            + (cum_length + other.cum_gap_lengths).tolist()
        )
        cum_gap_lengths = numpy.array(cum_gap_lengths, dtype=_DEFAULT_GAP_DTYPE)

        return self.__class__(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_gap_lengths,
            parent_length=self.parent_length + other.parent_length,
        )

    def __mul__(self, scale: int) -> "IndelMap":
        """used for going from amino-acid alignment to codon alignment"""
        gap_pos = self.gap_pos * scale
        cum_gap_lengths = self.cum_gap_lengths * scale
        return self.__class__(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_gap_lengths,
            parent_length=self.parent_length * scale,
        )

    def __repr__(self) -> str:
        gap_data = numpy.array([self.gap_pos, self.cum_gap_lengths]).T
        return f"{gap_data.tolist()!r}/{self.parent_length}"

    def get_gap_lengths(self) -> IntArrayTypes:
        lengths = self.cum_gap_lengths.copy()
        lengths[1:] = numpy.diff(lengths)
        return lengths

    def nongap(self) -> Iterator[SpanTypes]:
        """ungappeed segments in this map in aligned coordinates"""
        # we want to know the coordinates of the ungapped segments on
        # the aligned sequence. The gap_pos attribute is in sequence
        # coordinates
        prev_pos = 0
        for i, pos in enumerate(self.gap_pos):
            if pos == 0:
                # we start with a gap
                prev_pos = pos
                continue

            cum_length = 0 if i == 0 else self.cum_gap_lengths[i - 1]
            start = 0 if i == 0 else prev_pos
            start += cum_length
            end = self.gap_pos[i] + cum_length
            yield Span(start, end)
            prev_pos = pos

        if self.num_gaps and self.gap_pos[-1] + self.cum_gap_lengths[-1] < len(self):
            yield Span(self.gap_pos[-1] + self.cum_gap_lengths[-1], len(self))

    @property
    def spans(self) -> Iterator[SpanTypes]:
        """generator of spans"""
        if not self.num_gaps:
            yield Span(0, self.parent_length)
            return

        for i, pos in enumerate(self.gap_pos):
            cum_length = self.cum_gap_lengths[i]
            if pos == 0:
                cls = TerminalPadding if self.termini_unknown else LostSpan
                yield cls(cum_length)
                continue

            if i == 0:
                start = 0
                prev_length = 0
            else:
                start = self.gap_pos[i - 1]
                prev_length = self.cum_gap_lengths[i - 1]
            yield Span(start, pos)
            cls = (
                TerminalPadding
                if self.termini_unknown and i == self.num_gaps - 1
                else LostSpan
            )
            yield cls(cum_length - prev_length)

        if self.num_gaps and self.gap_pos[-1] < self.parent_length:
            yield Span(self.gap_pos[-1], self.parent_length)

    @property
    def complete(self) -> bool:
        """whether any span represents a gap"""
        return self.num_gaps != 0 and self.useful

    @property
    def useful(self) -> bool:
        return self.parent_length != 0

    def get_coordinates(self) -> SeqCoordTypes:
        """returns sequence coordinates of ungapped segments

        Returns
        -------
        [(start, end), ...]
        """
        if not self.num_gaps or (self.num_gaps == 1 and not self.gap_pos[0]):
            return [(0, int(self.parent_length))]
        elif self.num_gaps == 1:
            # does not start with a gap
            starts = [0, int(self.gap_pos[0])]
            ends = [int(self.gap_pos[0]), self.parent_length]
            return list(zip(starts, ends))

        starts = self.gap_pos[:-1].tolist()
        ends = self.gap_pos[1:].tolist()
        if self.gap_pos[0]:
            # does not start with a gap
            ends = starts[:1] + ends
            starts = [0] + starts

        if self.gap_pos[-1] + self.cum_gap_lengths[-1] < self.parent_length:
            # does end with a gap
            starts.append(ends[-1])
            ends.append(self.parent_length)

        return list(zip(starts, ends))

    def get_gap_coordinates(self) -> SeqCoordTypes:
        """returns [(gap pos, gap length), ...]"""
        lengths = self.get_gap_lengths()
        return numpy.array([self.gap_pos, lengths]).T.tolist()

    def get_gap_align_coordinates(self) -> SeqCoordTypes:
        """returns [(gap start, gap end), ...] in alignment indices

        Returns
        -------
        A 2D numpy array of integers. If the result is empty, it still
        has shape (0, 2).
        """
        starts, ends = _gap_spans(self.gap_pos, self.cum_gap_lengths)
        result = numpy.array([starts, ends]).T
        if not len(result):
            result = result.reshape((0, 2))
        return result

    def merge_maps(
        self, other: "IndelMap", parent_length: Optional[int] = None
    ) -> "IndelMap":
        """merge gaps of other with self

        Parameters
        ----------
        indel_map
            instance for same sequence
        parent_length
            overrides property
        """
        unique_pos = numpy.union1d(self.gap_pos, other.gap_pos)
        gap_lengths = numpy.zeros(unique_pos.shape, dtype=_DEFAULT_GAP_DTYPE)
        self_lengths = self.get_gap_lengths()
        other_lengths = other.get_gap_lengths()
        _update_lengths(unique_pos, gap_lengths, self.gap_pos, self_lengths)
        _update_lengths(unique_pos, gap_lengths, other.gap_pos, other_lengths)
        parent_length = parent_length or self.parent_length
        return self.__class__(
            gap_pos=unique_pos,
            gap_lengths=gap_lengths,
            parent_length=parent_length,
        )

    def joined_segments(self, coords: SeqCoordTypes) -> "IndelMap":
        """returns new map with disjoint gapped segments joined

        Parameters
        ----------
        coords
            sequence insert gap coordinates [(gap start, gap end), ...]
        """
        coords = sorted(coords)
        # using a dict here because joining can produce a gap merge
        gaps = {}
        cum_length = 0
        cum_parent_length = 0
        for start, end in coords:
            im = self[start:end]
            for i in range(im.num_gaps):
                pos = im.gap_pos[i] + cum_parent_length
                gaps[pos] = gaps.get(pos, cum_length) + im.cum_gap_lengths[i]
            cum_parent_length += im.parent_length
            if im.num_gaps:
                cum_length += im.cum_gap_lengths[-1]
        gap_pos = numpy.empty(len(gaps), dtype=_DEFAULT_GAP_DTYPE)
        cum_lengths = numpy.empty(len(gaps), dtype=_DEFAULT_GAP_DTYPE)
        for i, (pos, length) in enumerate(sorted(gaps.items())):
            gap_pos[i] = pos
            cum_lengths[i] = length
        return self.__class__(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_lengths,
            parent_length=cum_parent_length,
        )

    def nucleic_reversed(self) -> "IndelMap":
        """map for a sequence that has itself been reversed and complemented

        Notes
        -----
        discards reverse attribute on both spans and self
        """
        new_pos = self.gap_pos.copy()
        lengths = self.get_gap_lengths()

        if len(new_pos):
            new_pos = self.parent_length - new_pos
            new_pos = new_pos[::-1]
            lengths = lengths[::-1]

        return self.__class__(
            gap_pos=new_pos,
            gap_lengths=lengths,
            parent_length=self.parent_length,
        )

    def to_rich_dict(self) -> dict[str, Any]:
        """returns dicts for contained spans [dict(), ..]"""
        # exclude spans from deep copy since being overwritten
        data = copy.deepcopy(dict(self._serialisable.items()))
        data["type"] = get_object_provenance(self)
        data["version"] = __version__
        data["gap_pos"] = self.gap_pos.tolist()
        data["cum_gap_lengths"] = self.cum_gap_lengths.tolist()
        data["parent_length"] = int(self.parent_length)
        return data

    def to_json(self) -> str:
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, map_element) -> "IndelMap":
        from cogent3.util.deserialise import _get_class

        map_element.pop("version", None)
        type_ = map_element.pop("type")
        assert _get_class(type_) == cls
        map_element["gap_pos"] = numpy.array(map_element["gap_pos"])
        map_element["cum_gap_lengths"] = numpy.array(map_element["cum_gap_lengths"])

        return cls(**map_element)

    def with_termini_unknown(self) -> "IndelMap":
        """returns new instance with terminal gaps indicated as unknown"""
        return self.__class__(
            gap_pos=self.gap_pos.copy(),
            cum_gap_lengths=self.cum_gap_lengths.copy(),
            parent_length=self.parent_length,
            termini_unknown=True,
        )

    def to_feature_map(self) -> "FeatureMap":
        """returns a Map type, suited to Features"""
        return FeatureMap(spans=list(self.spans), parent_length=self.parent_length)

    def make_seq_feature_map(self, align_feature_map: "FeatureMap") -> "FeatureMap":
        """converts align_feature_map to a FeatureMap with sequence coordinates

        Parameters
        ----------
        align_feature_map
            with alignment coordinates

        Notes
        -----
        LostSpans in align_feature_map are skipped
        """
        spans = []
        for span in align_feature_map.spans:
            if span.lost:
                continue

            start = self.get_seq_index(span.start)
            end = self.get_seq_index(span.end)
            spans.append(Span(start, end))

        return FeatureMap(spans=spans, parent_length=self.parent_length)

    @functools.singledispatchmethod
    def shared_gaps(self, other: Union["IndelMap", numpy.ndarray]) -> numpy.ndarray:
        """returns a numpy array of the shared [(gap start, gap end), ...]

        Notes
        -----
        The result is in alignment coordinates
        """
        if len(self) != len(other):
            raise AssertionError(
                f"{len(other)=} != {len(self)=}, from a different alignment?"
            )
        if self.num_gaps == 0 or other.num_gaps == 0:
            return numpy.array([], dtype=_DEFAULT_GAP_DTYPE)

        return self.shared_gaps(other.get_gap_align_coordinates())

    @shared_gaps.register
    def _(self, other: numpy.ndarray) -> numpy.ndarray:
        """returns a numpy array of the shared [(gap start, gap end), ...]

        Notes
        -----
        The result is in alignment coordinates
        """
        other_gaps = other
        if not len(other_gaps):
            return numpy.array([], dtype=_DEFAULT_GAP_DTYPE)

        if other_gaps[-1][-1] > len(self):
            raise AssertionError(
                f"other_gaps {other_gaps!r} from a different alignment?"
            )
        self_gaps = self.get_gap_align_coordinates()

        result = coords_intersect(self_gaps, other_gaps)
        return numpy.array(result, dtype=_DEFAULT_GAP_DTYPE)

    @functools.singledispatchmethod
    def minus_gaps(self, other_gaps: Union["IndelMap", numpy.ndarray]):
        """returns new map with gaps in other_gaps removed from self

        Parameters
        ----------
        other_gaps
            IndelMap instance, or numpy array in alignment coordinates

        Returns
        -------
        New instance with gap spans represented by other_gaps removed from self
        """
        if len(self) != len(other_gaps):
            raise AssertionError(
                f"{len(other_gaps)=} != {len(self)=}, from a different alignment?"
            )

        return self.minus_gaps(other_gaps.get_gap_align_coordinates())

    @minus_gaps.register
    def _(self, other_gaps: numpy.ndarray):
        if not len(other_gaps):
            return self

        if other_gaps[-1][-1] > len(self):
            raise AssertionError(
                f"other_gaps {other_gaps!r} from a different alignment?"
            )
        self_gaps = self.get_gap_align_coordinates()
        unique = coords_minus_coords(self_gaps, other_gaps)
        new_gaps = numpy.empty((len(unique), 2), dtype=_DEFAULT_GAP_DTYPE)
        for i, (start, end) in enumerate(unique):
            new_gaps[i] = self.get_seq_index(start), end - start
        gap_pos, gap_lengths = new_gaps.T
        return self.__class__(
            gap_pos=gap_pos,
            gap_lengths=gap_lengths,
            parent_length=self.parent_length,
        )


_empty = None, None


def coords_minus_coords(
    coords1: numpy.ndarray, coords2: numpy.ndarray
) -> numpy.ndarray:
    """returns the coords1 minus any overlaps with coords2

    Parameters
    ----------
    coords1, coords2
        coordinates defining segments as [(gap start, gap end)..]

    Notes
    ------
    Assumes both coordinate series are sorted.
    """
    dtype = getattr(coords1, "dtype", _DEFAULT_GAP_DTYPE)
    unique_segments = []
    for a1, a2 in coords1:
        total_intersect = None
        for b1, b2 in coords2:
            if b2 < a1:
                continue

            if a2 <= b1:
                # no intersection
                break

            result = span_and_span((a1, a2), (b1, b2))
            if result[0] is not None:
                total_intersect = result[1] - result[0] + (total_intersect or 0)

        end = a2 - (total_intersect or 0)
        if end < 0:
            raise ValueError(
                f"new length negative segment length {a1=} {a2=} {total_intersect=}"
            )
        elif a2 - a1 != total_intersect:
            unique_segments.append((a1, end))

    return numpy.array(unique_segments, dtype=dtype)


def span_and_span(
    spans1: tuple[int, int], spans2: tuple[int, int]
) -> Union[tuple[int, int], tuple[None, None]]:
    """returns the intersection of two spans

    Parameters
    ----------
    spans1, span2
        each span is a tuple of start, end coordinates

    Returns
    -------
        If no intersection, returns (None, None)

    Raises
    ------
    ValueError
        if a span is not in span[0] < span[1] order
    """
    a1, a2 = spans1
    b1, b2 = spans2
    # return intersection of the spans
    if a1 >= a2 or b1 >= b2:
        raise ValueError("coordinates must be start < end")

    if a1 < b1 and a2 > b2:
        # span1 contains span2
        return b1, b2
    elif a1 >= b1 and a2 <= b2:
        # span1 equal to or within span2
        return a1, a2
    elif a1 == b1:
        # same start, different end
        return a1, min(a2, b2)
    elif a2 == b2:
        # different start, same end
        return max(a1, b1), a2
    elif a1 < b1 < a2:
        # span1 overlaps start of span2
        return b1, min(a2, b2)
    elif a1 < b2 < a2:
        # span1 overlaps end of span2
        return max(a1, b1), b2
    else:
        return _empty


def coords_intersect(
    coords1: numpy.ndarray, coords2: numpy.ndarray
) -> list[tuple[int, int]]:
    """returns the intersecting spans between two sets of coordinates

    Parameters
    ----------
    coords1, coords2
        coordinates defining segments as [(gap start, gap end)..]

    Notes
    ------
    Assumes both coordinate series are sorted.
    """
    intersects = []
    for a1, a2 in coords1:
        for b1, b2 in coords2:
            if a1 <= b2 and b1 <= a2:
                i1, i2 = span_and_span((a1, a2), (b1, b2))
                if i1 is not None:
                    intersects.append((i1, i2))
            elif a2 < b1:
                break
    return intersects


@dataclasses.dataclass
class FeatureMap(MapABC):
    """A map holds a list of spans."""

    spans: dataclasses.InitVar[Optional[SeqSpanTypes]] = ()
    parent_length: int = 0
    offsets: list[int] = dataclasses.field(init=False, repr=False)
    useful: bool = dataclasses.field(init=False, repr=False, default=False)
    complete: bool = dataclasses.field(init=False, repr=False, default=True)
    _serialisable: dict = dataclasses.field(init=False, repr=False)
    _spans: SeqSpanTypes = dataclasses.field(default=(), init=False)
    _start: Optional[int] = dataclasses.field(default=None, init=False)
    _end: Optional[int] = dataclasses.field(default=None, init=False)

    def __post_init__(self, spans: SeqSpanTypes):
        assert self.parent_length is not None
        self.parent_length = int(self.parent_length)
        if isinstance(spans, property):
            # This clause is due a known issue with dataclasses.
            # As we have a spans property, the default spans value is
            # ignored, so we have to check for its value being a property
            # and then set the default value here
            spans = ()

        spans = tuple(spans)

        self.offsets = []
        self.useful = False
        self.complete = True
        posn = 0
        for span in spans:
            self.offsets.append(posn)
            posn += span.length
            if span.lost:
                self.complete = False
            elif not self.useful:
                self.useful = True
                self._start, self._end = span.start, span.end
            else:
                self._start = min(self._start, span.start)
                self._end = max(self._end, span.end)

        self._spans = tuple(spans)
        self.length = posn

    @classmethod
    def from_spans(
        cls, spans: SeqSpanTypes, parent_length: int, **kwargs
    ) -> "FeatureMap":
        return cls(spans=spans, parent_length=parent_length)

    def __len__(self) -> int:
        return self.length

    def __repr__(self) -> str:
        return f"{list(self.spans)!r}/{self.parent_length}"

    def __getitem__(self, new_map) -> "FeatureMap":
        # A possible shorter map at the same level
        new_map = as_map(new_map, len(self), self.__class__)
        new_parts = []
        for span in new_map.spans:
            new_parts.extend(span.remap_with(self))
        return self.__class__(spans=new_parts, parent_length=self.parent_length)

    def __mul__(self, scale) -> "FeatureMap":
        new_parts = [span * scale for span in self.spans]
        return self.__class__(spans=new_parts, parent_length=self.parent_length * scale)

    def __truediv__(self, scale: int) -> "FeatureMap":
        new_parts = [span / scale for span in self.spans]
        return self.__class__(
            spans=new_parts, parent_length=self.parent_length // scale
        )

    def __add__(self, other) -> "FeatureMap":
        if other.parent_length != self.parent_length:
            raise ValueError("Those maps belong to different sequences")
        return self.__class__(
            spans=list(self.spans) + list(other.spans), parent_length=self.parent_length
        )

    @property
    def spans(self) -> Iterator[SeqSpanTypes]:
        yield from self._spans

    @property
    def num_spans(self):
        return len(self._spans)

    def get_covering_span(self) -> "FeatureMap":
        span = (self.start, self.end)
        return self.__class__.from_locations(
            locations=[span], parent_length=self.parent_length
        )

    def covered(self) -> "FeatureMap":
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
        return self.__class__.from_locations(
            locations=result, parent_length=self.parent_length
        )

    def nucleic_reversed(self) -> "FeatureMap":
        """map for a sequence that has itself been reversed and complemented

        Notes
        -----
        discards reverse attribute on both spans and self
        """
        spans = []
        parent_length = self.parent_length
        for s in self.spans:
            if not s.lost:
                start = parent_length - s.end
                assert start >= 0
                end = start + s.length
                s = Span(start=start, end=end)
            spans.append(s)

        spans.reverse()
        return self.__class__(spans=spans, parent_length=self.parent_length)

    def get_gap_coordinates(self) -> SeqCoordTypes:
        """returns [(gap pos, gap length), ...]"""
        gap_pos = []
        spans = list(self.spans)
        for i, span in enumerate(spans):
            if not span.lost:
                continue

            pos = spans[i - 1].end if i else 0
            gap_pos.append((pos, len(span)))

        return gap_pos

    def gaps(self) -> "FeatureMap":
        """The gaps (lost spans) in this map"""
        locations = []
        offset = 0
        for s in self.spans:
            if s.lost:
                locations.append((offset, offset + s.length))
            offset += s.length
        return self.__class__.from_locations(
            locations=locations, parent_length=len(self)
        )

    def shadow(self) -> "FeatureMap":
        """The 'negative' map of the spans not included in this map"""
        return self.inverse().gaps()

    def nongap(self) -> SeqSpanTypes:
        locations = []
        offset = 0
        for s in self.spans:
            if not s.lost:
                locations.append((offset, offset + s.length))
            offset += s.length
        return _spans_from_locations(locations=locations, parent_length=len(self))

    def without_gaps(self) -> "FeatureMap":
        return self.__class__(
            spans=[s for s in self.spans if not s.lost],
            parent_length=self.parent_length,
        )

    def inverse(self) -> "FeatureMap":
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

    def get_coordinates(self) -> SeqCoordTypes:
        """returns span coordinates as [(v1, v2), ...]

        v1/v2 are (start, end) unless the map is reversed, in which case it will
        be (end, start)"""

        return [(s.start, s.end) for s in self.spans if not s.lost]

    def to_rich_dict(self) -> dict[str, Any]:
        """returns dicts for contained spans [dict(), ..]"""
        spans = [s.to_rich_dict() for s in self.spans]
        data = copy.deepcopy(self._serialisable)
        data.pop("locations", None)
        data["spans"] = spans
        data["type"] = get_object_provenance(self)
        data["version"] = __version__
        data["parent_length"] = int(self.parent_length)
        return data

    def to_json(self) -> str:
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, map_element) -> "FeatureMap":
        from cogent3.util.deserialise import _get_class

        map_element.pop("version", None)
        type_ = map_element.pop("type")
        assert _get_class(type_) == cls
        spans = []
        for element in map_element.pop("spans"):
            element.pop("version", None)
            klass = _get_class(element.pop("type"))
            instance = klass(**element)
            spans.append(instance)

        map_element["spans"] = spans
        return cls(**map_element)

    def zeroed(self) -> "FeatureMap":
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

        zeroed._start = 0
        zeroed._end = new_end

        return zeroed

    def absolute_position(self, rel_pos: IntTypes) -> IntTypes:
        """converts rel_pos into an absolute position

        Raises
        ------
        raises ValueError if rel_pos < 0
        """
        check = (
            numpy.array([rel_pos], dtype=_DEFAULT_GAP_DTYPE)
            if isinstance(rel_pos, int)
            else rel_pos
        )
        if check.min() < 0:
            raise ValueError(f"must positive, not {rel_pos=}")

        return rel_pos if len(self) == self.parent_length else self.start + rel_pos

    def relative_position(self, abs_pos: IntTypes) -> IntTypes:
        """converts abs_pos into an relative position

        Raises
        ------
        raises ValueError if abs_pos < 0
        """
        check = (
            numpy.array([abs_pos], dtype=_DEFAULT_GAP_DTYPE)
            if isinstance(abs_pos, int)
            else abs_pos
        )
        if check.min() < 0:
            raise ValueError(f"must positive, not {abs_pos=}")
        return abs_pos - self.start

    @property
    def start(self):
        return self._start or 0

    @property
    def end(self):
        return self._end or 0


def gap_coords_to_map(
    gaps_lengths: dict[IntTypes, IntTypes], seq_length: int
) -> IndelMap:
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
        gap_pos = numpy.array([], dtype=_DEFAULT_GAP_DTYPE)
        lengths = gap_pos.copy()
    else:
        gap_pos, lengths = list(zip(*sorted(gaps_lengths.items())))
        gap_pos = numpy.array(gap_pos, dtype=_DEFAULT_GAP_DTYPE)
        lengths = numpy.array(lengths, dtype=_DEFAULT_GAP_DTYPE)

    return IndelMap(gap_pos=gap_pos, gap_lengths=lengths, parent_length=seq_length)


@register_deserialiser(get_object_provenance(IndelMap))
def deserialise_indelmap(data: dict) -> IndelMap:
    return IndelMap.from_rich_dict(data)


@register_deserialiser(get_object_provenance(FeatureMap))
def deserialise_featuremap(data: dict) -> FeatureMap:
    return FeatureMap.from_rich_dict(data)


@functools.singledispatch
def _norm_slice(index, length: int) -> tuple[int, int, Union[int, None]]:
    """_norm_slice(slice(1, -2, 3), 10) -> (1,8,3)"""
    start = index
    if start < 0:
        start += length
    if start >= length:
        raise IndexError(index)
    return start, start + 1, 1


@_norm_slice.register
def _(index: slice, length: int) -> tuple[int, int, Union[int, None]]:
    start = _norm_index(index.start, length, 0)
    end = _norm_index(index.stop, length, length)
    return start, end, index.step


@_norm_slice.register
def _(index: Span, length) -> tuple[int, int, Union[int, None]]:
    start = _norm_index(index.start, length, 0)
    end = _norm_index(index.end, length, length)
    return start, end, None


@_norm_slice.register
def _(index: FeatureMap, length) -> tuple[int, int, Union[int, None]]:
    start = _norm_index(index.start, length, 0)
    end = _norm_index(index.end, length, length)
    return start, end, None
