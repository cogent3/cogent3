from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Self, SupportsIndex

from cogent3._version import __version__
from cogent3.core.location import _input_vals_neg_step, _input_vals_pos_step
from cogent3.util.misc import get_object_provenance


class SliceRecordABC(ABC):
    """Abstract base class for recording the history of operations to be applied
    to some underlying data. Provides slicing functionality for the underlying data.

    Parameters
    ----------
    start
        start of the slice (inclusive indexing)
    stop
        stop of the slice (exclusive indexing)
    step
        step of the slice
    offset
        can be set with any additional offset that exists before the start of
        the underlying data
    parent_len
        length of the underlying data (not including offset)
    """

    __slots__ = ("_offset", "start", "step", "stop")

    @abstractmethod
    def __init__(
        self,
        *,
        parent_len: int,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        offset: int = 0,
    ) -> None:
        self.start: int
        self.stop: int
        self.step: int

    @abstractmethod
    def __eq__(self, other: object) -> bool: ...

    @abstractmethod
    def __ne__(self, other: object) -> bool: ...
    @property
    @abstractmethod
    def parent_len(self) -> int: ...

    @abstractmethod
    def _get_init_kwargs(self) -> dict[str, Any]:
        """return required arguments for construction that are unique to the
        subclass"""
        ...

    @abstractmethod
    def copy(self) -> Self: ...

    # refactor: design
    # can we remove the need for this method on the ABC and inheriting
    @property
    @abstractmethod
    def _zero_slice(self) -> Self: ...

    @property
    def offset(self) -> int:
        return self._offset

    @offset.setter
    def offset(self, value: int) -> None:
        value = value or 0
        self._offset = int(value)

    @property
    def plus_start(self) -> int:
        """start on plus strand"""
        if self.is_reversed:
            # all indices are negative, so the stop becomes the start. To convert
            # to positive indices, we add the length of the parent sequence to
            # the stop. Note that when abs(self.step) > 1, we instead calculate
            # the "true stop", i.e., the index that immediately follows (in the
            # direction of the step) the last selected index in the slice. We then
            # add the abs(self.step) to account for the fact that the stop uses
            # exclusive indexing, and the start is inclusive.

            assert self.stop < 0, "expected stop on reverse strand SeqView < 0"
            stop = self.stop if self.step == -1 else self.start + len(self) * self.step
            start = self.parent_len + stop + abs(self.step)
        else:
            start = self.start
        return start

    @property
    def plus_stop(self) -> int:
        """stop on plus strand"""
        if self.is_reversed:
            # self.start becomes the stop, self.start will be negative
            assert self.start < 0, "expected start on reverse strand SeqView < 0"
            stop = self.start + self.parent_len + 1
        else:
            # we ensure that the plus_stop is the index that immediately follows
            # the last selected index in the slice.
            stop = (
                self.stop
                if self.step == 1
                else (self.start + len(self) * self.step - self.step) + 1
            )
        return stop

    @property
    def plus_step(self) -> int:
        """step on plus strand"""
        return abs(self.step)

    @property
    def parent_start(self) -> int:
        """returns the start on the parent plus strand

        Returns
        -------
        offset + start, taking into account whether reversed. Result
        is positive.

        Notes
        -----
        This should NOT be used for slicing on the parent object.
        """
        return self.offset + self.plus_start

    @property
    def is_reversed(self) -> bool:
        return self.step < 0

    @property
    def parent_stop(self) -> int:
        """returns the stop on the parent plus strand

        Returns
        -------
        offset + stop, taking into account whether reversed. Result
        is positive.

        Notes
        -----
        This should NOT be used for slicing on the parent object.
        """
        return self.offset + self.plus_stop

    def absolute_position(self, rel_index: int, include_boundary: bool = False) -> int:
        """Converts an index relative to the current view to be with respect
        to the coordinates of the original "Python sequence".

        Parameters
        ----------
        rel_index
            relative position with respect to the current view
        include_boundary
            whether considering index as part of a range

        Returns
        -------
        the absolute index with respect to the coordinates of the original
        sequence (including offset if present).
        """
        if not self:
            return 0

        if rel_index < 0:
            msg = "only positive indexing supported!"
            raise IndexError(msg)

        # _get_index return the absolute position relative to the underlying sequence
        seq_index, _, _ = self._get_index(rel_index, include_boundary=include_boundary)

        offset = self.offset
        return (
            offset + self.parent_len + seq_index + 1
            if self.is_reversed
            else offset + seq_index
        )

    def relative_position(self, abs_index: int, stop: bool = False) -> int:
        """converts an index on the original "Python sequence" into an index
        on this "view"

        Notes
        -----
        The returned value DOES NOT reflect python indexing. Importantly,
        negative values represent positions that precede the current view.
        """
        if not self:
            return 0

        if abs_index < 0:
            msg = "Index must be +ve and relative to the + strand"
            raise IndexError(msg)

        if self.is_reversed:
            offset = self.offset

            if (
                tmp := (self.parent_len - abs_index + offset + self.start + 1)
            ) % self.step == 0 or stop:
                rel_pos = tmp // abs(self.step)
            else:
                rel_pos = (tmp // abs(self.step)) + 1

        else:
            offset = self.offset + self.start

            if (tmp := abs_index - offset) % self.step == 0 or stop:
                rel_pos = tmp // self.step
            else:
                rel_pos = (tmp // self.step) + 1

        return rel_pos

    def __len__(self) -> int:
        return abs((self.start - self.stop) // self.step)

    def __getitem__(self, segment: SupportsIndex | slice) -> Self:
        kwargs = self._get_init_kwargs()

        if isinstance(segment, SupportsIndex):
            start, stop, step = self._get_index(segment)
            return self.__class__(
                start=start,
                stop=stop,
                step=step,
                offset=self.offset,
                parent_len=self.parent_len,
                **kwargs,
            )

        if segment.start is segment.stop is segment.step is None:
            return self.copy()

        if len(self) == 0:
            return self

        if segment.start is not None and segment.start == segment.stop:
            return self._zero_slice

        slice_step = 1 if segment.step is None else segment.step

        if slice_step > 0:
            return self._get_slice(segment, slice_step, **kwargs)
        if slice_step < 0:
            return self._get_reverse_slice(segment, slice_step, **kwargs)
        msg = f"{self.__class__.__name__} cannot be sliced with a step of 0"
        raise ValueError(
            msg,
        )

    def _get_index(
        self,
        val: SupportsIndex,
        include_boundary: bool = False,
    ) -> tuple[int, int, int]:
        if len(self) == 0:
            raise IndexError(val)

        val = int(val)
        if (
            (val > 0 and include_boundary and val > len(self))
            or (val > 0 and not include_boundary and val >= len(self))
            or (val < 0 and include_boundary and abs(val) > (len(self) + 1))
            or (val < 0 and not include_boundary and abs(val) > len(self))
        ):
            raise IndexError(val)

        if self.step > 0:
            if val >= 0:
                val = self.start + val * self.step
            else:
                val = self.start + len(self) * self.step + val * abs(self.step)

            return val, val + 1, 1

        if self.step < 0:
            if val >= 0:
                val = self.start + val * self.step
            else:
                val = self.start + len(self) * self.step + val * self.step

            return val, val - 1, -1
        msg = f"Invalid step: {self.step}"
        raise ValueError(msg)

    def _get_slice(self, segment: slice, step: int, **kwargs: Any) -> Self:
        slice_start = segment.start if segment.start is not None else 0
        slice_stop = segment.stop if segment.stop is not None else len(self)

        if self.step > 0:
            return self._get_forward_slice_from_forward_seqview_(
                slice_start,
                slice_stop,
                step,
                **kwargs,
            )

        if self.step < 0:
            return self._get_forward_slice_from_reverse_seqview_(
                slice_start,
                slice_stop,
                step,
                **kwargs,
            )
        msg = f"Invalid step: {self.step}"
        raise ValueError(msg)

    def _get_forward_slice_from_forward_seqview_(
        self,
        slice_start: int,
        slice_stop: int,
        step: int,
        **kwargs: Any,
    ) -> Self:
        start = (
            self.start + slice_start * self.step
            if slice_start >= 0
            else max(
                self.start + len(self) * self.step + slice_start * self.step,
                self.start,
            )
        )
        if slice_stop > self.stop:
            stop = self.stop
        elif slice_stop >= 0:
            stop = self.start + slice_stop * self.step
        else:
            # "true stop" adjust for if abs(stop-start) % step != 0
            # "true stop" = self.start + len(self) * self.step
            stop = self.start + len(self) * self.step + slice_stop * self.step

        # if -ve, it's an invalid slice
        if start < 0 or stop < 0:
            return self._zero_slice

        # checking for zero-length slice
        if stop < start:
            return self._zero_slice
        if start > self.parent_len:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=min(self.stop, stop),
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )

    def _get_forward_slice_from_reverse_seqview_(
        self,
        slice_start: int,
        slice_stop: int,
        step: int,
        **kwargs: Any,
    ) -> Self:
        if slice_start >= 0:
            start = self.start + slice_start * self.step
        elif abs(slice_start) > len(self):
            start = self.start
        else:
            start = self.start + len(self) * self.step + slice_start * self.step

        if slice_stop >= 0:
            stop = self.start + slice_stop * self.step
        else:  # slice_stop < 0
            stop = self.start + len(self) * self.step + slice_stop * self.step

        # if +ve, it's an invalid slice
        if start >= 0 or stop >= 0:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=max(self.stop, stop),
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )

    def _get_reverse_slice(
        self,
        segment: slice,
        step: int,
        **kwargs: Any,
    ) -> Self:
        slice_start = segment.start if segment.start is not None else -1
        slice_stop = segment.stop if segment.stop is not None else -len(self) - 1

        if self.step < 0:
            return self._get_reverse_slice_from_reverse_seqview_(
                slice_start,
                slice_stop,
                step,
                **kwargs,
            )
        if self.step > 0:
            return self._get_reverse_slice_from_forward_seqview_(
                slice_start,
                slice_stop,
                step,
                **kwargs,
            )

        msg = f"Invalid step: {self.step}"
        raise ValueError(msg)

    def _get_reverse_slice_from_forward_seqview_(
        self,
        slice_start: int,
        slice_stop: int,
        step: int,
        **kwargs: Any,
    ) -> Self:
        # "true stop" adjust for if abs(stop-start) % step != 0
        # max possible start is "true stop" - step, because stop is not inclusive
        # "true stop" - step is converted to -ve index via subtracting len(self)
        if slice_start >= len(self):
            start = (self.start + len(self) * self.step - self.step) - self.parent_len
        elif slice_start >= 0:
            start = (self.start + slice_start * self.step) - self.parent_len
        else:
            start = (
                self.start
                + len(self) * self.step
                + slice_start * self.step
                - self.parent_len
            )

        if slice_stop >= self.parent_len:
            return self._zero_slice

        if slice_stop >= 0:
            stop = self.start + (slice_stop * self.step) - self.parent_len
        else:
            stop = (
                self.start
                + (len(self) * self.step)
                + (slice_stop * self.step)
                - self.parent_len
            )

        if start >= 0 or stop >= 0:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=max(stop, self.start - self.parent_len - 1),
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )

    def _get_reverse_slice_from_reverse_seqview_(
        self,
        slice_start: int,
        slice_stop: int,
        step: int,
        **kwargs: Any,
    ) -> Self:
        # refactor: simplify
        # are there other places in slice_record where we can use plus_start/stop/step
        if slice_start >= len(self):
            start = self.plus_start
        elif slice_start >= 0:
            start = self.parent_len + (self.start + slice_start * self.step)
        else:
            start = self.parent_len + (
                self.start + len(self) * self.step + slice_start * self.step
            )

        if slice_stop >= 0:
            stop = self.parent_len + (self.start + slice_stop * self.step)
            if stop <= self.parent_len + self.stop:
                return self._zero_slice
        else:
            stop = self.parent_len + (
                self.start + len(self) * self.step + slice_stop * self.step
            )
            if stop > self.parent_len + self.start:
                stop = self.parent_len + self.start + 1

        # if -ve, it's an invalid slice becomes zero
        # checking for zero-length slice
        if stop < start or start > self.parent_len or min(start, stop) < 0:
            return self._zero_slice

        return self.__class__(
            start=start,
            stop=stop,
            step=self.step * step,
            offset=self.offset,
            parent_len=self.parent_len,
            **kwargs,
        )


class SliceRecord(SliceRecordABC):
    """records cumulative slice operations on an object without modifying it

    Notes
    -----
    Basis for lazy evaluation of slicing operations on sequences and alignments.
    A reference to an instance of this class is used by different view objects.
    """

    __slots__ = ("_parent_len",)

    def __init__(
        self,
        *,
        parent_len: int,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        offset: int = 0,
    ) -> None:
        """
        Parameters
        ----------
        start
            start of the slice (inclusive indexing)
        stop
            stop of the slice (exclusive indexing)
        step
            step of the slice
        offset
            can be set with any additional offset that exists before the start of
            the underlying data
        parent_len
            length of the underlying data (not including offset)
        """
        if step == 0:
            msg = "step cannot be 0"
            raise ValueError(msg)
        step = step if step is not None else 1
        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
        start, stop, step = func(parent_len, start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        self._parent_len = parent_len
        self._offset = offset or 0

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SliceRecordABC):
            return False
        return (
            self.start == other.start
            and self.stop == other.stop
            and self.step == other.step
            and self._parent_len == other.parent_len
            and self._offset == other.offset
        )

    def __ne__(self, other: object) -> bool:
        return not self == other

    @property
    def parent_len(self) -> int:
        return self._parent_len

    def _get_init_kwargs(self) -> dict[str, Any]:
        return {}

    def copy(self, sliced: bool = False) -> Self:
        """returns self"""
        return self

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(start={self.start}, stop={self.stop}, step={self.step}, "
            f"parent_len={self.parent_len}, offset={self.offset})"
        )

    @property
    def _zero_slice(self) -> Self:
        return self.__class__(start=0, stop=0, step=1, parent_len=self.parent_len)

    def to_rich_dict(self) -> dict[str, str | dict[str, int]]:
        """returns dict suitable for serialisation"""
        data: dict[str, str | dict[str, int]] = {
            "type": get_object_provenance(self),
            "version": __version__,
        }
        data["init_args"] = {
            "parent_len": self.parent_len,
            "start": self.start,
            "stop": self.stop,
            "step": self.step,
            "offset": self.offset,
        }
        return data
