from __future__ import annotations

import re
import typing

from collections import defaultdict
from dataclasses import InitVar, dataclass, field
from functools import singledispatch, singledispatchmethod
from pathlib import Path
from typing import Callable, Iterator, Mapping, Optional, Union

import numpy

from cogent3 import get_moltype
from cogent3.core.alphabet import CharAlphabet
from cogent3.core.location import IndelMap
from cogent3.core.moltype import MolType
from cogent3.core.sequence import (
    Sequence,
    SliceRecordABC,
    _input_vals_neg_step,
    _input_vals_pos_step,
)
from cogent3.util import warning as c3warn
from cogent3.util.misc import get_setting_from_environ


PrimitiveSeqTypes = Union[str, bytes, numpy.ndarray]


def assign_sequential_names(num_seqs, base_name="seq", start_at=0):
    """Returns list of num_seqs sequential, unique names."""
    return [f"{base_name}_{i}" for i in range(start_at, start_at + num_seqs)]


@singledispatch
def seq_index(seq: PrimitiveSeqTypes, alphabet: CharAlphabet) -> numpy.ndarray:
    raise NotImplementedError(
        f"{seq_index.__name__} not implemented for type {type(seq)}"
    )


@seq_index.register
def _(seq: str, alphabet: CharAlphabet) -> numpy.ndarray:
    return alphabet.to_indices(seq)


@seq_index.register
def _(seq: bytes, alphabet: CharAlphabet) -> numpy.ndarray:
    return alphabet.to_indices(seq)


@seq_index.register
def _(seq: numpy.ndarray, alphabet: CharAlphabet) -> numpy.ndarray:
    return seq.astype(alphabet.array_type)


class SeqDataView(SliceRecordABC):
    """
    A view class for SeqsData, providing properties for different
    representations.

    self.seq is a SeqsData() instance, but other properties are a reference to a
    single seqid only.

    Example
    -------
    data = {"seq1": "ACGT", "seq2": "GTTTGCA"}
    sd = SeqsData(data=data)
    sdv = sd.get_seq_view(seqid="seq1")
    """

    __slots__ = ("seq", "start", "stop", "step", "_offset", "_seqid", "_seq_len")

    def __init__(
        self,
        *,
        seq: SeqsData,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        offset: int = 0,
        seqid: Optional[str] = None,
        seq_len: Optional[int] = None,
    ):
        if step == 0:
            raise ValueError("step cannot be 0")
        step = 1 if step is None else step

        self._seq_len = self._checked_seq_len(seq_len)
        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
        start, stop, step = func(self._seq_len, start, stop, step)
        self.seq = seq
        self.start = start
        self.stop = stop
        self.step = step
        self._offset = offset
        self._seqid = seqid

    def _checked_seq_len(self, seq_len: int) -> int:
        assert seq_len is not None
        return seq_len

    @property
    def _zero_slice(self):
        return self.__class__(
            seq=self.seq, seqid=self.seqid, seq_len=self.seq_len, start=0, stop=0
        )

    @property
    def seqid(self) -> str:
        return self._seqid

    @property
    def seq_len(self) -> int:
        return self._seq_len

    def _get_init_kwargs(self):
        return {"seq": self.seq, "seqid": self.seqid}

    @property
    def str_value(self) -> str:
        raw = self.seq.get_seq_str(
            seqid=self.seqid, start=self.parent_start, stop=self.parent_stop
        )
        return raw if self.step == 1 else raw[:: self.step]

    @property
    def array_value(self) -> numpy.ndarray:
        raw = self.seq.get_seq_array(
            seqid=self.seqid, start=self.parent_start, stop=self.parent_stop
        )
        return raw if self.step == 1 else raw[:: self.step]

    @property
    def bytes_value(self) -> bytes:
        raw = self.seq.get_seq_bytes(
            seqid=self.seqid, start=self.parent_start, stop=self.parent_stop
        )
        return raw if self.step == 1 else raw[:: self.step]

    def replace(self, old, new):
        # todo: kath, placeholder for testing
        return self

    def __str__(self):
        return self.str_value

    def __array__(self):
        return self.array_value

    def __bytes__(self):
        return self.bytes_value

    def __repr__(self) -> str:
        seq = f"{self[:10]!s}...{self[-5:]}" if len(self) > 15 else str(self)
        return (
            f"{self.__class__.__name__}(seq={seq}, start={self.start}, "
            f"stop={self.stop}, step={self.step}, offset={self.offset}, "
            f"seqid={self.seqid!r}, seq_len={self.seq_len})"
        )

    @classmethod
    def to_rich_dict(cls):
        # todo: kath, placeholder until i delete to_rich_dict from SliceRecordABC
        ...

    @classmethod
    def from_rich_dict(cls, data: dict):
        init_args = data.pop("init_args")
        if "offset" in data:
            init_args["offset"] = data.pop("offset")
        return cls(**init_args)

    # todo: do we support copy? do we support copy with sliced?
    def copy(self, sliced: bool = False):
        """returns copy"""

        return self.__class__(
            seq=self.seq,
            start=self.start,
            stop=self.stop,
            step=self.step,
            offset=self.offset,
            seqid=self.seqid,
            seq_len=self.seq_len,
        )


class SeqsData:
    __slots__ = ("_data", "_alphabet", "_names", "_make_seq")

    def __init__(
        self,
        data: dict[str, PrimitiveSeqTypes],
        alphabet: CharAlphabet,
        make_seq: Optional[type] = None,
    ):
        self._alphabet = alphabet
        self._make_seq = make_seq
        # todo: kath, convert seq_index to using alphabet to do this directly
        self._data: dict[str, numpy.ndarray] = {
            k: seq_index(v, self._alphabet) for k, v in data.items()
        }

    @property
    def names(self) -> list:
        return list(self._data.keys())

    @property
    def make_seq(self) -> Union[SeqDataView, Sequence]:
        return self._make_seq

    @make_seq.setter
    def make_seq(self, make_seq: Callable) -> None:
        self._make_seq = make_seq

    @property
    def num_seqs(self) -> int:
        return len(self._data)

    @singledispatchmethod
    def __getitem__(self, index: Union[str, int]) -> SeqDataView:
        raise NotImplementedError(f"__getitem__ not implemented for {type(index)}")

    @__getitem__.register
    def _(self, index: str) -> SeqDataView:
        sdv = self.get_seq_view(seqid=index)
        return sdv if self._make_seq is None else self.make_seq(sdv, name=index)

    @__getitem__.register
    def _(self, index: int) -> SeqDataView:
        return self[self.names[index]]

    def get_seq_array(
        self, *, seqid: str, start: int = None, stop: int = None
    ) -> numpy.ndarray:
        return self._data[seqid][start:stop]

    def get_seq_str(self, *, seqid: str, start: int = None, stop: int = None) -> str:
        return self._alphabet.from_indices(
            self.get_seq_array(seqid=seqid, start=start, stop=stop)
        )

    def get_seq_bytes(
        self, *, seqid: str, start: int = None, stop: int = None
    ) -> bytes:
        return self.get_seq_str(seqid=seqid, start=start, stop=stop).encode("utf8")

    def get_seq_view(self, seqid: str) -> SeqDataView:
        seq_len = len(self._data[seqid])
        return SeqDataView(seq=self, seqid=seqid, seq_len=seq_len)


@singledispatch
def seq_to_gap_coords(
    seq: Union[str, numpy.ndarray], moltype: MolType
) -> tuple[str, IndelMap]:
    """
    Takes a sequence with (or without) gaps and returns an ungapped sequence
    and records the position and length of gaps in the original parent sequence
    """
    raise NotImplementedError(f"{seq} not implemented for type {type(seq)}")


@seq_to_gap_coords.register
def _(seq: str, moltype: MolType) -> tuple[str, IndelMap]:
    seq = moltype.make_seq(seq)
    indel_map, ungapped_seq = seq.parse_out_gaps()

    if indel_map.num_gaps == 0:
        return str(ungapped_seq), numpy.array([], dtype=int)

    return str(ungapped_seq), indel_map


@seq_to_gap_coords.register
def _(seq: numpy.ndarray, moltype: MolType) -> tuple[str, IndelMap]:
    gap_char = moltype.alphabets.degen_gapped.index(moltype.gap)
    # Create boolean array of gaps to get ungapped
    gaps_bool = seq == gap_char
    ungapped = seq[~gaps_bool]

    # no gaps in seq
    if numpy.array_equal(ungapped, seq):
        return ungapped, numpy.array([], dtype=int)

    parent_len = len(seq)
    in_gap = False
    parent_coords = []
    start = 0
    for i, gapped in enumerate(gaps_bool):
        if gapped and not in_gap:
            start = i
            in_gap = True
        if not gapped and in_gap:
            parent_coords.append([start, i])
            in_gap = False
        if gapped and i == parent_len - 1:
            # End of the sequence
            parent_coords.append([start, i + 1])

    # Format for IndelMap
    parent_coords = numpy.array(parent_coords)
    gap_lengths = numpy.array([end - start for start, end in parent_coords])
    cum_lengths = gap_lengths.cumsum()
    # Get gap start positions in sequence coords
    gap_pos = parent_coords.T[0] - numpy.append(0, cum_lengths[:-1, numpy.newaxis])

    indel_map = IndelMap(
        gap_pos=gap_pos, cum_gap_lengths=cum_lengths, parent_length=parent_len
    )

    return ungapped, indel_map


@singledispatch
def get_seqs_data(
    data: Union[dict, SeqsData], alphabet: CharAlphabet
) -> dict[str, numpy.ndarray]:
    raise NotImplementedError(f"SeqsData can not be constructed for type {type(data)}")


@get_seqs_data.register
def _(data: SeqsData, alphabet: CharAlphabet) -> dict[str, numpy.ndarray]:
    # todo: kath, perform integrity check between alphabet provided and SeqsData alphabet
    return data


@get_seqs_data.register
def _(data: dict, alphabet: CharAlphabet) -> dict[str, numpy.ndarray]:
    # apply name conversion here?
    seq_data = SeqsData(data=data, alphabet=alphabet)
    return get_seqs_data(seq_data, alphabet)


def make_unaligned_seqs(
    *,
    data: Union[dict[str, PrimitiveSeqTypes], SeqsData],
    moltype: Union[str, MolType],
    label_to_name: Callable = None,
    info: dict = None,
    source: Union[str, Path] = None,
    **kwargs,
):
    """Initialize an unaligned collection of sequences.

    Parameters
    ----------
    data
        sequences, assumes a dict {name: seq, ...} or a SeqsData
    moltype
        string representation of the moltype, e.g., 'dna', 'protein'.
    label_to_name
        function for converting original names into other names.
    info
        a dict from which to make an info object
    source
        origins of this data, defaults to 'unknown'. Converted to a string
        and added to info["source"].
    **kwargs
        other keyword arguments to be passed to SequenceCollection
    """
    moltype = get_moltype(moltype)
    alphabet = moltype.alphabet

    if len(data) == 0:
        raise ValueError("data must be at least one sequence.")
    seq_data = get_seqs_data(data, alphabet)

    return SequenceCollection(seq_data=seq_data, moltype=moltype)


class SequenceCollection:
    def __init__(self, *, seq_data: SeqsData, moltype: MolType):
        self.moltype = moltype
        self._seq_data = seq_data
        self._seq_data.make_seq = self.moltype.make_seq
        self._repr_policy = dict(num_seqs=10, num_pos=60, ref_name="longest", wrap=60)

    @property
    def names(self) -> Iterator[str]:
        yield from self._seq_data.names

    @property
    def seqs(self) -> SeqsData:
        return self._seq_data

    def iter_seqs(self, seq_order: list = None) -> Iterator[Union[Sequence, SeqsData]]:
        """Iterates over values (sequences) in the alignment, in order.

        Parameters
        ----------
        seq_order:
            list of seqids giving the order in which seqs will be returned.
            Defaults to self.names
        """

        seqs = self.seqs
        get = seqs.__getitem__
        for key in seq_order or self.names:
            yield get(key)

    @property
    def num_seqs(self) -> int:
        return self._seq_data.num_seqs

    @property
    @c3warn.deprecated_callable(
        version="2025.5", reason=".seqs can now be indexed by name", new=".seqs"
    )
    def named_seqs(self) -> SeqsData:  # pragma: no cover
        return self.seqs



class AlignedDataView(SeqDataView):
    """
    Example
    -------
    data = {"seq1": "ACGT", "seq2": "GTTTGCA"}
    ad = AlignedData.from_strings(data)
    adv = sd.get_aligned_view(seqid="seq1")
    """


@dataclass
class AlignedData:
    # Look out for any overlaps with SeqsData
    seqs: Optional[dict[str, numpy.ndarray]] = None
    gaps: Optional[dict[str, numpy.ndarray]] = None
    _moltype: MolType = field(init=False)
    _names: tuple[str] = field(init=False)
    _alpha: CharAlphabet = field(init=False)
    align_len: int = 0
    moltype: InitVar[Union[str, None]] = "dna"
    names: InitVar[Union[tuple[str], None]] = None

    def __post_init__(self, moltype, names):
        self._moltype = get_moltype(moltype)
        self._alpha = self._moltype.alphabets.degen_gapped
        self.seqs = {k: seq_index(v, self._alpha) for k, v in self.seqs.items()}

    @classmethod
    def from_gapped_seqs(
        cls,
        data: dict[str, Union[str, numpy.ndarray]],
        moltype: Union[str, MolType] = "dna",
        names: Optional[tuple[str]] = None,
    ):
        """
        Convert dict of {"seq_name": "seq"} to two dicts for seqs and gaps
        """
        seq_lengths = {len(v) for v in data.values()}
        if len(seq_lengths) != 1:
            raise ValueError("All sequence lengths must be the same.")

        align_len = seq_lengths.pop()
        moltype = get_moltype(moltype)

        seqs = {}
        gaps = {}
        for name, seq in data.items():
            seqs[name], gaps[name] = seq_to_gap_coords(seq, moltype)

        names = names or list(data.keys)

        return cls(
            seqs=seqs,
            gaps=gaps,
            moltype=moltype,
            names=names,
            align_len=align_len,
        )

    def get_aligned_view(self, seqid: str) -> AlignedDataView:
        # Need to revisit what the variable is called i.e. parent_length
        return AlignedDataView(seq=self, seqid=seqid, seq_len=self.align_len)

    def get_gaps(self, seqid: str) -> numpy.ndarray:
        return self.gaps[seqid]
