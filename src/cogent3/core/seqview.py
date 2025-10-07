from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, SupportsIndex, cast

import numpy
import numpy.typing as npt
from typing_extensions import Self

from cogent3._version import __version__
from cogent3.core.location import IndelMap
from cogent3.core.slice_record import SliceRecord, SliceRecordABC
from cogent3.util.misc import get_object_provenance

if TYPE_CHECKING:
    import cogent3.core.alphabet as c3_alphabet
    from cogent3.core.seq_storage import AlignedSeqsDataABC, SeqsDataABC

NumpyIntArrayType = npt.NDArray[numpy.integer]


class SeqViewABC(ABC):
    """
    An abstract base class for providing a view of a sequence.

    This class defines an interface for sequence views, by which operations
    performed on the sequence do not modify the original sequence data. Instead,
    modifications and operations are stored on the returned view of the sequence
    and can be realised by accessing the values of the view.
    """

    __slots__ = ("_parent_len", "_slice_record", "alphabet", "parent")

    def __init__(self) -> None:
        self.alphabet: c3_alphabet.CharAlphabet[Any]
        self.parent: str | bytes | NumpyIntArrayType | SeqsDataABC
        self._parent_len: int
        self._slice_record: SliceRecordABC

    @property
    @abstractmethod
    def seqid(self) -> str | None: ...

    @property
    @abstractmethod
    def parent_len(self) -> int: ...

    @property
    @abstractmethod
    def slice_record(self) -> SliceRecordABC: ...

    @property
    def parent_offset(self) -> int:
        """returns the offset from the true parent

        Notes
        -----
        If from storage with an offset attribute, returns
        that value or 0
        """
        return (
            cast("SeqsDataABC", self.parent).offset.get(cast("str", self.seqid), 0)
            if hasattr(self.parent, "offset")
            else 0
        )

    @property
    @abstractmethod
    def offset(self) -> int: ...

    @property
    def is_reversed(self) -> bool:
        """whether the sequence is reversed"""
        return self.slice_record.is_reversed

    @property
    @abstractmethod
    def str_value(self) -> str: ...

    @property
    @abstractmethod
    def array_value(self) -> NumpyIntArrayType: ...

    @property
    @abstractmethod
    def bytes_value(self) -> bytes: ...

    @abstractmethod
    def copy(self, sliced: bool = False) -> Self: ...

    @abstractmethod
    def __str__(self) -> str: ...

    @abstractmethod
    def __array__(
        self,
        dtype: numpy.dtype[numpy.integer] | None = None,
        copy: bool | None = None,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def __bytes__(self) -> bytes: ...

    @abstractmethod
    def __getitem__(self, segment: SupportsIndex | slice) -> Self: ...

    def __len__(self) -> int:
        return len(self.slice_record)

    def _get_init_kwargs(self) -> dict[str, Any]:
        return {
            "parent": self.parent,
            "parent_len": self._parent_len,
            "seqid": self.seqid,
            "alphabet": self.alphabet,
            "slice_record": self.slice_record,
        }

    def with_offset(self, offset: int) -> Self:
        """returns new instance with annotation offset set"""
        if self._slice_record.offset:
            msg = f"cannot set {offset=} on a SeqView with an offset {self._slice_record.offset=}"
            raise ValueError(
                msg,
            )

        init_kwargs = self._get_init_kwargs()
        init_kwargs["offset"] = offset
        return self.__class__(**init_kwargs)

    @abstractmethod
    def parent_coords(
        self, *, apply_offset: bool = False, **kwargs: Any
    ) -> tuple[str, int, int, int]: ...


class SeqView(SeqViewABC):
    """
    Provides a view of a sequence with support for slicing operations.

    This class represents a view of a sequence. It uses ``SliceRecord`` to
    enable efficient slicing without altering the original sequence data.

    Parameters
    ----------
    parent
        the original sequence data
    alphabet
        the alphabet object defining valid characters for the sequence
    seqid
        the name or identifier of the sequence
    parent_len
        the length of the sequence. Defaults to the length of the input sequence

    Notes
    -----
    It utilises the alphabet object to allow providing different types
    such as strings or numpy arrays, corresponding to its underlying data.
    """

    __slots__ = (
        "_parent_len",
        "_seqid",
        "_slice_record",
        "alphabet",
        "parent",
    )

    def __init__(
        self,
        *,
        parent: str | bytes | NumpyIntArrayType | SeqsDataABC,
        alphabet: c3_alphabet.CharAlphabet[Any],
        parent_len: int,
        seqid: str | None = None,
        slice_record: SliceRecordABC | None = None,
        offset: int = 0,
    ) -> None:
        self.alphabet = alphabet
        self.parent = parent
        self._seqid = seqid
        assert parent_len is not None
        self._parent_len = parent_len
        self._slice_record = (
            slice_record
            if slice_record is not None
            else SliceRecord(parent_len=self._parent_len)
        )
        if offset and self._slice_record.offset:
            msg = f"cannot set {offset=} on a SeqView with an offset {self._slice_record.offset=}"
            raise ValueError(
                msg,
            )
        if offset:
            self._slice_record.offset = offset

    @property
    def seqid(self) -> str | None:
        """name of the sequence"""
        return self._seqid

    @property
    def slice_record(self) -> SliceRecordABC:
        """the slice state of this view"""
        return self._slice_record

    @property
    def parent_len(self) -> int:
        """length of the parent sequence"""
        return self._parent_len

    @property
    def offset(self) -> int:
        """the annotation offset of this view"""
        return self.slice_record.offset

    def _get_init_kwargs(self) -> dict[str, Any]:
        return {
            "parent": self.parent,
            "parent_len": self._parent_len,
            "seqid": self.seqid,
            "alphabet": self.alphabet,
            "slice_record": self.slice_record,
        }

    @property
    def str_value(self) -> str:
        """returns the sequence as a string"""
        self.parent = cast("str | bytes | NumpyIntArrayType", self.parent)
        return self.alphabet.from_indices(
            self.parent[
                self.slice_record.start : self.slice_record.stop : self.slice_record.step
            ],
        )

    @property
    def array_value(self) -> NumpyIntArrayType:
        """returns the sequence as a array of indices"""
        self.parent = cast("str | bytes | NumpyIntArrayType", self.parent)
        return self.alphabet.to_indices(
            self.parent[
                self.slice_record.start : self.slice_record.stop : self.slice_record.step
            ],
            validate=False,
        )

    @property
    def bytes_value(self) -> bytes:
        """returns the sequence as a bytes string"""
        return self.str_value.encode("utf-8")

    def __str__(self) -> str:
        return self.str_value

    def __array__(
        self,
        dtype: numpy.dtype[numpy.integer] | None = None,
        copy: bool | None = None,
    ) -> NumpyIntArrayType:
        arr = self.array_value
        if dtype is not None:
            arr = arr.astype(dtype)
        return arr

    def __bytes__(self) -> bytes:
        return self.bytes_value

    def __getitem__(self, segment: SupportsIndex | slice) -> Self:
        return self.__class__(
            parent=self.parent,
            seqid=self.seqid,
            alphabet=self.alphabet,
            parent_len=self.parent_len,
            slice_record=self.slice_record[segment],
        )

    def __repr__(self) -> str:
        self.parent = cast("str | bytes | NumpyIntArrayType", self.parent)
        seq_preview = (
            f"{self.parent[:10]}...{self.parent[-5:]}"  # type: ignore[str-bytes-safe]
            if self.parent_len > 15
            else self.parent
        )
        return (
            f"{self.__class__.__name__}(seqid={self.seqid!r}, parent={seq_preview!r}, "
            f"slice_record={self.slice_record.__repr__()})"
        )

    def copy(self, sliced: bool = False) -> Self:
        """returns copy

        Parameters
        ----------
        sliced
            if True, the underlying sequence is truncated and the start/stop
            adjusted
        """
        if not sliced:
            return self.__class__(
                parent=self.parent,
                seqid=self.seqid,
                alphabet=self.alphabet,
                slice_record=self.slice_record.copy(),
                parent_len=self.parent_len,
            )
        # we slice parent with the plus attributes only, so, if reversed, we need to
        # retain that information in the slice_record
        sr = SliceRecord(parent_len=len(self), step=-1) if self.is_reversed else None
        self.parent = cast("str | bytes | NumpyIntArrayType", self.parent)
        return self.__class__(
            parent=self.parent[
                self.slice_record.plus_start : self.slice_record.plus_stop : self.slice_record.plus_step
            ],
            parent_len=len(self),
            seqid=self.seqid,
            alphabet=self.alphabet,
            slice_record=sr,
        )

    def parent_coords(
        self, *, apply_offset: bool = False, **kwargs: Any
    ) -> tuple[str, int, int, int]:
        """returns coordinates on parent

        Parameters
        ----------
        apply_offset
            if True adds annotation offset from parent

        Returns
        -------
        parent seqid, start, stop, strand
        """
        offset = self.parent_offset if apply_offset else 0
        return (
            cast("str", self.seqid),
            self.slice_record.parent_start + offset,
            self.slice_record.parent_stop + offset,
            self.slice_record.step,
        )


class SeqDataView(SeqView):
    """
    A view class for ``SeqsData``, providing methods for different
    representations of a single sequence.

    Notes
    -----
    ``str_value`` / ``array_value`` are not complemented, but can be reversed.
    The latter is done by the ``Sequence`` object which has a moltype.
    """

    __slots__ = ("_parent_len", "_seqid", "_slice_record", "alphabet", "parent")

    def __init__(
        self,
        *,
        parent: SeqsDataABC,
        alphabet: c3_alphabet.CharAlphabet[Any],
        parent_len: int,
        seqid: str | None = None,
        slice_record: SliceRecordABC | None = None,
        offset: int = 0,
    ) -> None:
        self.parent: SeqsDataABC
        super().__init__(
            parent=parent,
            alphabet=alphabet,
            parent_len=parent_len,
            seqid=seqid,
            slice_record=slice_record,
            offset=offset,
        )

    @property
    def offset(self) -> int:
        """the annotation offset of this view"""
        return self.slice_record.offset

    @property
    def str_value(self) -> str:
        """returns the sequence as a string"""
        return self.alphabet.from_indices(self.array_value)

    @property
    def array_value(self) -> NumpyIntArrayType:
        """returns the sequence as a numpy array"""
        # we select the data using plus strand coords
        raw = self.parent.get_seq_array(
            seqid=cast("str", self.seqid),
            start=self.slice_record.plus_start,
            stop=self.slice_record.plus_stop,
            step=self.slice_record.plus_step,
        )
        if self.slice_record.is_reversed:
            # and reverse result when a reversed slice
            raw = raw[::-1]
        return raw

    @property
    def bytes_value(self) -> bytes:
        """returns the sequence as bytes"""
        return self.str_value.encode("utf8")

    def __repr__(self) -> str:
        seq = f"{self[:10]!s}...{self[-5:]}" if len(self) > 15 else str(self)
        return (
            f"{self.__class__.__name__}(seqid={self.seqid!r}, parent={seq}, "
            f"slice_record={self.slice_record!r})"
        )

    # refactor: design, do we support copy? do we support copy with sliced?
    def copy(self, sliced: bool = False) -> Self:
        """returns copy"""
        return self

    def to_rich_dict(self) -> dict[str, Any]:
        """returns a json serialisable dict.

        Notes
        -----
        This method will slice the underlying sequence to the start and stop values

        Warnings
        --------
        This method is not intended to provide serialisation of this object,
        instead, it is intended for usage by an enclosing class.
        """

        data: dict[str, Any] = {
            "type": get_object_provenance(self),
            "version": __version__,
        }
        data["init_args"] = self._get_init_kwargs()

        if self.slice_record.is_reversed:
            adj = self.parent_len + 1
            start, stop = self.slice_record.stop + adj, self.slice_record.start + adj
        else:
            start, stop = self.slice_record.start, self.slice_record.stop

        data["init_args"]["parent"] = self.str_value[start:stop]
        new_sr = SliceRecord(
            parent_len=(stop - start),
            step=self.slice_record.step,
            offset=self.slice_record.parent_start,
        )
        data["init_args"]["slice_record"] = new_sr.to_rich_dict()
        data["init_args"]["alphabet"] = self.alphabet.to_rich_dict()
        return data

    @property
    def is_reversed(self) -> bool:
        if self.seqid in self.parent.reversed_seqs:
            # seqid is reversed relative to everything else
            # hence is_reversed is the opposite of the slice record
            return not self.slice_record.is_reversed
        return self.slice_record.is_reversed

    def parent_coords(
        self, *, apply_offset: bool = False, **kwargs: Any
    ) -> tuple[str, int, int, int]:
        """returns coordinates on parent

        Parameters
        ----------
        apply_offset
            if True adds annotation offset from parent

        Returns
        -------
        parent seqid, start, stop, strand
        """
        offset = self.parent_offset if apply_offset else 0
        start = self.slice_record.parent_start
        stop = self.slice_record.parent_stop
        step = self.slice_record.step
        return cast("str", self.seqid), start + offset, stop + offset, step


class AlignedDataViewABC(SeqViewABC):
    __slots__ = ()

    def __init__(self) -> None:
        self.parent: AlignedSeqsDataABC

    @abstractmethod
    def get_seq_view(self) -> SeqViewABC: ...

    @property
    @abstractmethod
    def map(self) -> IndelMap: ...

    @property
    @abstractmethod
    def slice_record(self) -> SliceRecordABC: ...

    @property
    @abstractmethod
    def gapped_str_value(self) -> str: ...

    @property
    @abstractmethod
    def gapped_array_value(self) -> NumpyIntArrayType: ...

    @property
    @abstractmethod
    def gapped_bytes_value(self) -> bytes: ...


class AlignedDataView(AlignedDataViewABC):
    """
    A view class for ``AlignedSeqsData``, providing methods for different representations
    of a single sequence.

    Notes
    -----
    ``str_value`` / ``array_value`` are not complemented, but can be reversed. The latter
    is done by the ``Aligned`` object which has a moltype. The ``slice_record`` attribute
    is shared with the containing ``Alignment``.
    """

    __slots__ = (
        "_offset",
        "_parent_len",
        "_seqid",
        "_slice_record",
        "alphabet",
        "parent",
    )

    def __init__(
        self,
        *,
        parent: AlignedSeqsDataABC,
        seqid: str,
        alphabet: c3_alphabet.CharAlphabet[Any],
        slice_record: SliceRecordABC | None = None,
    ) -> None:
        self.parent = parent
        self._seqid = seqid
        self.alphabet = alphabet
        self._parent_len = parent.align_len
        self._slice_record = (
            slice_record
            if slice_record is not None
            else SliceRecord(parent_len=self._parent_len)
        )

    @property
    def slice_record(self) -> SliceRecordABC:
        """the slice record for this view"""
        return self._slice_record

    @slice_record.setter
    def slice_record(self, value: SliceRecordABC) -> None:
        self._slice_record = value

    @property
    def offset(self) -> int:
        """the slice offset of this view"""
        return self.slice_record.offset

    @property
    def seqid(self) -> str:
        """the name of the sequence"""
        return self._seqid

    @property
    def parent_len(self) -> int:
        """length of the parent sequence"""
        return self._parent_len

    @property
    def map(self) -> IndelMap:
        """indel map (gaps) for the sequence"""
        imap = self._parent_map()
        start, stop, step = (
            self.slice_record.start,
            self.slice_record.stop,
            self.slice_record.step,
        )
        return imap[start:stop:step]

    def _parent_map(self) -> IndelMap:
        gap_pos_gap_length = self.parent.get_gaps(self.seqid)
        if gap_pos_gap_length.size > 0:
            gap_pos = numpy.array(gap_pos_gap_length[:, 0], dtype=int)
            cum_gap_lengths = numpy.array(gap_pos_gap_length[:, 1], dtype=int)
        else:
            gap_pos, cum_gap_lengths = (
                numpy.array([], dtype=int),
                numpy.array([], dtype=int),
            )
        return IndelMap(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_gap_lengths,
            parent_length=self.parent.get_seq_length(self.seqid),
        )

    @property
    def str_value(self) -> str:
        """returns the string value of the ungapped sequence"""
        return self.alphabet.from_indices(self.array_value)

    @property
    def gapped_str_value(self) -> str:
        """returns the string value of the gapped sequence"""
        return self.alphabet.from_indices(self.gapped_array_value)

    @property
    def array_value(self) -> NumpyIntArrayType:
        """returns the numpy array of indices for the ungapped sequence"""
        value = self.parent.get_seq_array(
            seqid=self.seqid,
            start=self.map.get_seq_index(self.slice_record.plus_start),
            stop=self.map.get_seq_index(self.slice_record.plus_stop),
            step=self.map.get_seq_index(self.slice_record.plus_step),
        )
        return value[::-1] if self.slice_record.is_reversed else value

    @property
    def gapped_array_value(self) -> NumpyIntArrayType:
        """returns the numpy array of indices for the gapped sequence"""
        value = self.parent.get_gapped_seq_array(
            seqid=self.seqid,
            start=self.slice_record.plus_start,
            stop=self.slice_record.plus_stop,
            step=self.slice_record.plus_step,
        )
        return value[::-1] if self.slice_record.is_reversed else value

    @property
    def bytes_value(self) -> bytes:
        """returns the bytes value of the ungapped sequence"""
        return self.str_value.encode("utf8")

    @property
    def gapped_bytes_value(self) -> bytes:
        """returns the bytes value of the gapped sequence"""
        return self.gapped_str_value.encode("utf8")

    def __str__(self) -> str:
        return self.gapped_str_value

    def __array__(
        self,
        dtype: numpy.dtype[numpy.integer] | None = None,
        copy: bool | None = None,
    ) -> NumpyIntArrayType:
        arr = self.gapped_array_value
        if dtype:
            arr = arr.astype(dtype)
        return arr

    def __bytes__(self) -> bytes:
        return self.gapped_bytes_value

    def __getitem__(self, segment: SupportsIndex | slice) -> Self:
        return self.__class__(
            parent=self.parent,
            seqid=self.seqid,
            alphabet=self.alphabet,
            slice_record=self.slice_record[segment],
        )

    def __repr__(self) -> str:
        seq_preview = (
            f"{self.parent.get_seq_array(seqid=self.seqid, start=0, stop=10)}..."
            f"{self.parent.get_seq_array(seqid=self.seqid, start=self.parent_len - 5)}"
            if self.parent_len > 15
            else self.parent.get_seq_array(seqid=self.seqid)
        )
        seq_preview = self.alphabet.from_indices(seq_preview)
        return (
            f"{self.__class__.__name__}(seqid={self.seqid!r}, map={self.map!r}, parent={seq_preview!r}, "
            f"slice_record={self.slice_record.__repr__()})"
        )

    def parent_coords(
        self, *, seq_coords: bool = False, apply_offset: bool = False, **kwargs: Any
    ) -> tuple[str, int, int, int]:
        """returns seqid, start, stop, strand on the parent

        Parameters
        ----------
        seq_coords
            if True, parent is the ungapped sequence
        apply_offset
            if True and seq_coords, adds annotation offset from parent
        """
        strand = -1 if self.is_reversed else 1
        if not seq_coords:
            return (
                self.seqid,
                self.slice_record.parent_start,
                self.slice_record.parent_stop,
                strand,
            )

        # AlignedDataView.parent_coords uses it's indelmap, etc..
        # to return the necessary coordinates

        # we want the coordinates on the parent sequence, which means we
        # need to use the parent's IndelMap for findings the correct indices.
        parent_map = self._parent_map()
        start = parent_map.get_seq_index(self.slice_record.parent_start)
        stop = parent_map.get_seq_index(self.slice_record.parent_stop)
        offset = self.parent_offset if apply_offset else 0

        return self.seqid, start + offset, stop + offset, strand

    def copy(self, sliced: bool = False) -> Self:
        """just returns self"""
        return self

    def _get_init_kwargs(self) -> dict[str, Any]:
        return {
            "parent": self.parent,
            "seqid": self.seqid,
            "alphabet": self.alphabet,
            "slice_record": self.slice_record,
        }

    def get_seq_view(self) -> SeqViewABC:
        """returns view of ungapped sequence data for seqid"""
        # we want the parent coordinates in sequence coordinates
        # parent_coords does not account for the stride
        seqid, start, stop, _ = self.parent_coords(seq_coords=True, apply_offset=False)
        parent_len = self.parent.get_seq_length(seqid)
        sr = SliceRecord(
            start=start,
            stop=stop,
            parent_len=parent_len,
        )[:: self.slice_record.step]

        return SeqDataView(
            parent=self.parent,
            seqid=seqid,
            alphabet=self.alphabet,
            parent_len=parent_len,
            slice_record=sr,
        )

    @property
    def is_reversed(self) -> bool:
        """whether the sliced view is reversed relative to the parent"""
        if self.seqid in self.parent.reversed_seqs:
            # seqid is reversed relative to everything else
            # hence is_reversed is the opposite of the slice record
            return not self.slice_record.is_reversed
        return self.slice_record.is_reversed
