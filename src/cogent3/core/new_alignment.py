from __future__ import annotations

import re
import typing

from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import InitVar, dataclass, field
from functools import singledispatch, singledispatchmethod
from pathlib import Path
from typing import Callable, Iterator, Mapping, Optional, Union

import numpy

import cogent3.core.new_alphabet as new_alpha
import cogent3.core.new_moltype as new_moltype
import cogent3.core.new_sequence as new_seq

from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import (
    BasicAnnotationDb,
    GenbankAnnotationDb,
    SupportsFeatures,
    load_annotations,
)
from cogent3.core.info import Info as InfoClass
from cogent3.core.location import IndelMap
from cogent3.util import warning as c3warn
from cogent3.util.misc import get_setting_from_environ, negate_condition


DEFAULT_ANNOTATION_DB = BasicAnnotationDb

PrimitiveSeqTypes = Union[str, bytes, numpy.ndarray]


def assign_sequential_names(num_seqs: int, base_name: str = "seq", start_at: int = 0):
    """Returns list of num_seqs sequential, unique names."""
    return [f"{base_name}_{i}" for i in range(start_at, start_at + num_seqs)]


class SeqDataView(new_seq.SliceRecordABC):
    """
    A view class for SeqsData, providing properties for different
    representations.

    self.seqs is a SeqsData() instance, but other properties are a reference to a
    single seqid only.

    Example
    -------
    data = {"seq1": "ACGT", "seq2": "GTTTGCA"}
    sd = SeqsData(data=data)
    sdv = sd.get_seq_view(seqid="seq1")
    """

    __slots__ = ("seqs", "start", "stop", "step", "_offset", "_seqid", "_seq_len")

    def __init__(
        self,
        *,
        seqs: SeqsData,
        seq_len: int,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        offset: int = 0,
        seqid: Optional[str] = None,
    ):
        if step == 0:
            raise ValueError("step cannot be 0")
        step = 1 if step is None else step

        self._seq_len = self._checked_seq_len(seq_len)
        func = (
            new_seq._input_vals_pos_step if step > 0 else new_seq._input_vals_neg_step
        )
        start, stop, step = func(self._seq_len, start, stop, step)
        self.seqs = seqs
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
            seqs=self.seqs, seqid=self.seqid, seq_len=self.seq_len, start=0, stop=0
        )

    @property
    def seqid(self) -> str:
        return self._seqid

    @property
    def seq_len(self) -> int:
        return self._seq_len

    def _get_init_kwargs(self):
        return {"seqs": self.seqs, "seqid": self.seqid}

    @property
    def str_value(self) -> str:
        raw = self.seqs.get_seq_str(
            seqid=self.seqid, start=self.parent_start, stop=self.parent_stop
        )
        return raw if self.step == 1 else raw[:: self.step]

    @property
    def array_value(self) -> numpy.ndarray:
        raw = self.seqs.get_seq_array(
            seqid=self.seqid, start=self.parent_start, stop=self.parent_stop
        )
        return raw if self.step == 1 else raw[:: self.step]

    @property
    def bytes_value(self) -> bytes:
        raw = self.seqs.get_seq_bytes(
            seqid=self.seqid, start=self.parent_start, stop=self.parent_stop
        )
        return raw if self.step == 1 else raw[:: self.step]

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

    # todo: do we support copy? do we support copy with sliced?
    def copy(self, sliced: bool = False):
        """returns copy"""
        return self


class SeqsDataABC(ABC):
    """Abstract base class for respresenting the collection of sequences underlying
    a SequenceCollection
    """

    alphabet: new_alpha.CharAlphabet
    make_seq: Callable = None

    @abstractmethod
    def seq_lengths(self) -> dict[str, int]: ...

    @property
    @abstractmethod
    def names(self) -> list:  # refactor: design
        ...

    @property
    @abstractmethod
    def make_seq(self): ...

    @make_seq.setter
    @abstractmethod
    def make_seq(self, make_seq: Callable) -> None: ...

    @property
    @abstractmethod
    def alphabet(self) -> new_alpha.CharAlphabet: ...

    @abstractmethod
    def get_seq_array(
        self, *, seqid: str, start: int = None, stop: int = None
    ) -> numpy.ndarray: ...

    @abstractmethod
    def get_seq_str(
        self, *, seqid: str, start: int = None, stop: int = None
    ) -> str: ...

    @abstractmethod
    def get_seq_bytes(
        self, *, seqid: str, start: int = None, stop: int = None
    ) -> bytes: ...

    @abstractmethod
    def get_seq_view(self, seqid: str) -> new_seq.SliceRecordABC: ...

    @abstractmethod
    def subset(self, names: Union[str, typing.Sequence[str]]) -> SeqsDataABC: ...

    @abstractmethod
    def to_alphabet(self, alphabet: new_alpha.CharAlphabet) -> SeqsDataABC: ...

    @abstractmethod
    def __len__(self): ...

    @abstractmethod
    def __getitem__(
        self, index: Union[str, int]
    ) -> Union[new_seq.Sequence, new_seq.SliceRecordABC]: ...


class SeqsData(SeqsDataABC):
    __slots__ = ("_data", "_alphabet", "_names", "_make_seq")

    def __init__(
        self,
        *,
        data: dict[str, PrimitiveSeqTypes],
        alphabet: new_alpha.CharAlphabet,
        make_seq: Callable = None,
    ):
        self._alphabet = alphabet
        self._make_seq = make_seq
        self._data: dict[str, numpy.ndarray] = {
            str(name): self._alphabet.to_indices(seq) for name, seq in data.items()
        }

    @property
    def names(self) -> list:  # refactor: design
        return list(self._data.keys())

    @property
    def make_seq(self) -> Union[SeqDataView, new_seq.Sequence]:
        """if set, returns a function that takes 'seq' and 'name' as keyword
        arguments and returns a given Sequence from the collection.

        Notes
        -----
        Can be set with any callable function that takes 'seq' and 'name' as
        keyword arguments. Typically set with '<moltype-instance>.make_seq'.
        """
        return self._make_seq

    @make_seq.setter
    def make_seq(self, make_seq: Callable) -> None:
        self._make_seq = make_seq

    # todo: kath, do we want a num_seqs property for SeqsData?

    @property
    def alphabet(self) -> new_alpha.CharAlphabet:
        return self._alphabet

    def seq_lengths(self) -> dict[str, int]:
        return {name: seq.shape[0] for name, seq in self._data.items()}

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
        return SeqDataView(seqs=self, seqid=seqid, seq_len=seq_len)

    def subset(self, names: Union[str, typing.Sequence[str]]) -> SeqsData:
        """Returns a new SeqsData object with only the specified names."""
        names = [names] if isinstance(names, str) else names
        if data := {name: self._data.get(name) for name in names}:
            return self.__class__(
                data=data, alphabet=self.alphabet, make_seq=self.make_seq
            )
        else:
            raise ValueError(f"provided {names=} not found in collection")

    def to_alphabet(self, alphabet: new_alpha.CharAlphabet) -> SeqsData:
        # todo: kath
        ...

    def __len__(self):
        return len(self.names)

    @singledispatchmethod
    def __getitem__(
        self, index: Union[str, int]
    ) -> Union[SeqDataView, new_seq.Sequence]:
        raise NotImplementedError(f"__getitem__ not implemented for {type(index)}")

    @__getitem__.register
    def _(self, index: str) -> Union[SeqDataView, new_seq.Sequence]:
        sdv = self.get_seq_view(seqid=index)
        return (
            sdv
            if self.make_seq is None
            else self.make_seq(seq=sdv.str_value, name=index)
        )

    @__getitem__.register
    def _(self, index: int) -> Union[SeqDataView, new_seq.Sequence]:
        return self[self.names[index]]


@singledispatch

@singledispatch
def make_seqsdata_dict(
    data, label_to_name: Callable = None
) -> dict[str, Union[PrimitiveSeqTypes, new_seq.Sequence]]:
    if isinstance(data, SeqsData):
        return data
    raise NotImplementedError(f"make_data_dict not implemented for {type(data)}")


@make_seqsdata_dict.register
def _(
    data: dict, label_to_name: Callable = None
) -> dict[str, Union[PrimitiveSeqTypes, new_seq.Sequence]]:
    return {label_to_name(k): v for k, v in data.items()} if label_to_name else data


@make_seqsdata_dict.register
def _(
    data: list, label_to_name: Callable = None
) -> dict[str, Union[PrimitiveSeqTypes, new_seq.Sequence]]:
    names = assign_sequential_names(len(data))
    return make_seqsdata_dict(dict(zip(names, data)), label_to_name=label_to_name)


def make_unaligned_seqs(
    *,
    data: Union[dict[str, PrimitiveSeqTypes], SeqsData, list],
    moltype: Union[str, new_moltype.MolType],
    label_to_name: Callable = None,
    info: dict = None,
    source: Union[str, Path] = None,
    annotation_db: SupportsFeatures = None,
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
    moltype = new_moltype.get_moltype(moltype)

    if len(data) == 0:
        raise ValueError("data must be at least one sequence.")

    seqs_data = make_seqsdata_dict(data, label_to_name=label_to_name)

    if isinstance(seqs_data, dict):
        seqs_data = SeqsData(data=seqs_data, alphabet=moltype.degen_gapped_alphabet)
    elif not moltype.is_compatible_alphabet(data.alphabet):
        raise ValueError(
            f"Provided moltype: {moltype} is not compatible with SeqsData alphabet {data.alphabet}"
        )

    return SequenceCollection(
        seqs_data=seqs_data, moltype=moltype, info=info, annotation_db=annotation_db
    )


class SequenceCollection:
    def __init__(
        self,
        *,
        seqs_data: SeqsData,
        moltype: new_moltype.MolType,
        names: list[str] = None,
        info: Union[dict, InfoClass] = None,
        annotation_db: SupportsFeatures = None,
    ):
        self._seqs_data = seqs_data
        self.moltype = moltype
        self.names = names
        self._seqs_data.make_seq = self.moltype.make_seq
        if not isinstance(info, InfoClass):
            info = InfoClass(info) if info else InfoClass()
        self.info = info
        self._repr_policy = dict(num_seqs=10, num_pos=60, ref_name="longest", wrap=60)
        self._annotation_db = annotation_db or DEFAULT_ANNOTATION_DB()

    @property
    def names(self) -> Iterator[str]:
        return self._names

    @names.setter
    def names(self, names: list[str]):  # refactor: design
        if names is None:
            self._names = self._seqs_data.names
        elif set(names) <= set(self._seqs_data.names):
            self._names = names
        else:
            left_diff = set(names) - set(self._seqs_data.names)
            raise ValueError(f"Provided names not found in collection: {left_diff}")

    @property
    def seqs(self) -> SeqsData:
        return self._seqs_data

    @property
    def num_seqs(self) -> int:
        return len(self.names)

    @property
    def annotation_db(self):
        return self._annotation_db

    @annotation_db.setter
    def annotation_db(self, value):
        if value == self._annotation_db:
            return

        self._annotation_db = value

        # todo: kath
        for seq in self.seqs:
            seq.replace_annotation_db(value, check=False)

    @property
    @c3warn.deprecated_callable(
        version="2025.5", reason=".seqs can now be indexed by name", new=".seqs"
    )
    def named_seqs(self) -> SeqsData:  # pragma: no cover
        return self.seqs

    def iter_seqs(
        self, seq_order: list = None
    ) -> Iterator[Union[new_seq.Sequence, SeqsData]]:
        """Iterates over sequences in the collection, in order.

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

    def take_seqs(
        self,
        names: Union[str, typing.Sequence[str]],
        negate: bool = False,
        copy_annotations: bool = False,
        **kwargs,
    ):
        """Returns new collection containing only specified seqs.

        Parameters
        ----------
        names
            sequences to select (or exclude if negate=True)
        negate
            select all sequences EXCEPT names
        kwargs
            keyword arguments to be passed to the constructor of the new collection
        copy_annotations
            if True, only annotations from selected seqs are copied to the annotation_db
            of the new collection

        Notes
        -----
        The seqs in the new collection will be references to the same objects as
        the seqs in the old collection.
        """
        # todo: kath, note gavins comment that this method could operate on only the names
        # list and not the seqs_data object

        if negate:
            names = [name for name in self.names if name not in names]
        else:
            names = [name for name in names if name in self.names]

        if not names:
            raise ValueError(f"{names=} and {negate=} resulted in no names")

        seqs_data = self.seqs.subset(names)

        result = self.__class__(
            seqs_data=seqs_data,
            moltype=self.moltype,
            info=self.info,
            **kwargs,
        )
        if self.annotation_db:
            if copy_annotations:
                result.annotation_db = type(self.annotation_db)()
                result.annotation_db.update(
                    annot_db=self.annotation_db, seqids=result.names
                )
            else:
                result.annotation_db = self.annotation_db
        return result

    def get_seq_names_if(
        self, f: Callable[[new_seq.Sequence], bool], negate: bool = False
    ):
        """Returns list of names of seqs where f(seq) is True."""
        get = self.seqs.__getitem__

        new_f = negate_condition(f) if negate else f

        return [name for name in self.names if new_f(get(name))]

    def take_seqs_if(
        self, f: Callable[[new_seq.Sequence], bool], negate: bool = False, **kwargs
    ):
        """Returns new collection containing seqs where f(seq) is True.

        Notes
        -----
        The seqs in the new collection are the same objects as the
        seqs in the old collection, not copies.
        """
        # pass take_seqs the result of get_seq_indices
        return self.take_seqs(self.get_seq_names_if(f, negate), **kwargs)

    def __eq__(self, other: Union[new_seq.SequenceCollection, dict]) -> bool:
        return id(self) == id(other)

    def __ne__(self, other: new_seq.SequenceCollection) -> bool:
        return not self.__eq__(other)

    def __repr__(self):
        seqs = []
        limit = 10
        delimiter = ""

        repr_seq_names = [min(self.names, key=lambda name: len(self.seqs[name]))]
        if len(self.names) > 1:
            # In case of a tie, min and max return first.
            # reversed ensures if all seqs are of same length, different seqs are returned
            repr_seq_names.append(
                max(reversed(self.names), key=lambda name: len(self.seqs[name]))
            )

        for name in repr_seq_names:
            elts = list(str(self.seqs[name])[: limit + 1])
            if len(elts) > limit:
                elts[-1] = "..."
            seqs.append(f"{name}[{delimiter.join(elts)}]")

        if len(self.names) > 2:
            seqs.insert(1, "...")

        seqs = ", ".join(seqs)

        return f"{len(self.names)}x ({seqs}) {self.moltype.label} seqcollection"

    def _repr_html_(self) -> str:
        settings = self._repr_policy.copy()
        env_vals = get_setting_from_environ(
            "COGENT3_ALIGNMENT_REPR_POLICY",
            dict(num_seqs=int, num_pos=int, wrap=int),
        )
        settings.update(env_vals)
        return self.to_html(
            name_order=self.names[: settings["num_seqs"]],
            limit=settings["num_pos"],
            wrap=settings["wrap"],
        )

    def to_html(
        self,
        name_order: Optional[typing.Sequence[str]] = None,
        wrap: int = 60,
        limit: Optional[int] = None,
        colors: Optional[Mapping[str, str]] = None,
        font_size: int = 12,
        font_family: str = "Lucida Console",
    ) -> str:
        """returns html with embedded styles for sequence colouring

        Parameters
        ----------
        name_order
            order of names for display.
        wrap
            number of alignment columns per row
        limit
            truncate alignment to this length
        colors
            {character
            moltype.
        font_size
            in points. Affects labels and sequence and line spacing
            (proportional to value)
        font_family
            string denoting font family

        Examples
        ---------

        In a jupyter notebook, this code is used to provide the representation.

        .. code-block:: python

            seq_col # is rendered by jupyter

        You can directly use the result for display in a notebook as

        .. code-block:: python

            from IPython.core.display import HTML
            HTML(seq_col.to_html())
        """
        css, styles = self.moltype.get_css_style(
            colors=colors, font_size=font_size, font_family=font_family
        )

        seq_lengths = numpy.array([len(seq) for seq in self.seqs])
        min_val = seq_lengths.min()
        max_val = seq_lengths.max()
        med_val = numpy.median(seq_lengths)

        if name_order:
            selected = self.take_seqs(name_order)
        else:
            name_order = self.names
            selected = self

        # Stylise each character in each sequence
        gaps = "".join(frozenset([selected.moltype.gap, selected.moltype.missing]))
        template = '<span class="%s">%%s</span>'
        styled_seqs = defaultdict(list)
        max_truncated_len = 0
        for name in name_order:
            sequence = str(self.seqs[name])[:limit]
            seq_len = len(sequence)
            max_truncated_len = max(seq_len, max_truncated_len)
            start_gap = re.search(f"^[{gaps}]+", sequence)
            end_gap = re.search(f"[{gaps}]+$", sequence)
            start = 0 if start_gap is None else start_gap.end()
            end = seq_len if end_gap is None else end_gap.start()

            seq = []
            for i, char in enumerate(sequence):
                if i < start or i >= end:
                    style = f"terminal_ambig_{self.moltype.label}"
                else:
                    style = styles[char]
                s = template % style
                s = s % char
                seq.append(s)

            styled_seqs[name] = seq

        # Ensure all sublists are of same length
        for name in styled_seqs:
            if len(styled_seqs[name]) < max_truncated_len:
                styled_seqs[name].extend(
                    [""] * (max_truncated_len - len(styled_seqs[name]))
                )

        # Make html table
        seqs = numpy.array([styled_seqs[n] for n in name_order], dtype="O")
        table = ["<table>"]
        seq_ = "<td>%s</td>"
        label_ = '<td class="label">%s</td>'
        num_row_ = '<tr class="num_row"><td></td><td><b>{:,d}</b></td></tr>'
        for i in range(0, max_truncated_len, wrap):
            table.append(num_row_.format(i))
            seqblock = seqs[:, i : i + wrap].tolist()
            for n, s in zip(name_order, seqblock):
                s = "".join(s)
                # Filter out rows that are empty (due to combination of shorter sequences + wrapping)
                if s != "":
                    row = "".join([label_ % n, seq_ % s])
                    table.append(f"<tr>{row}</tr>")
        table.append("</table>")
        if (
            limit
            and limit < len(selected.names)
            or name_order
            and len(name_order) < len(selected.names)
        ):
            summary = (
                "%s x {min=%s, median=%s, max=%s} (truncated to %s x %s) %s sequence collection"
            ) % (
                self.num_seqs,
                min_val,
                med_val,
                max_val,
                len(name_order) if name_order else len(selected.names),
                limit or len(selected),
                selected.moltype.label,
            )
        else:
            summary = ("%s x {min=%s, median=%s, max=%s} %s sequence collection") % (
                self.num_seqs,
                min_val,
                med_val,
                max_val,
                selected.moltype.label,
            )

        text = [
            "<style>",
            ".c3align table {margin: 10px 0;}",
            ".c3align td { border: none !important; text-align: left !important; }",
            ".c3align tr:not(.num_row) td span {margin: 0 2px;}",
            ".c3align tr:nth-child(even) {background: #f7f7f7;}",
            ".c3align .num_row {background-color:rgba(161, 195, 209, 0.5) !important; border-top: solid 1px black; }",
            ".c3align .label { font-size: %dpt ; text-align: right !important; "
            "color: black !important; padding: 0 4px; display: table-cell !important; "
            "font-weight: normal !important; }" % font_size,
            "\n".join([".c3align " + style for style in css]),
            "</style>",
            '<div class="c3align">',
            "\n".join(table),
            f"<p><i>{summary}</i></p>",
            "</div>",
        ]
        return "\n".join(text)

    def set_repr_policy(
        self,
        num_seqs: Optional[int] = None,
        num_pos: Optional[int] = None,
        ref_name: Optional[str] = None,
        wrap: Optional[int] = None,
    ):
        """specify policy for repr(self)

        Parameters
        ----------
        num_seqs
            number of sequences to include in represented display.
        num_pos
            length of sequences to include in represented display.
        ref_name
            name of sequence to be placed first, or "longest" (default).
            If latter, indicates longest sequence will be chosen.
        wrap
            number of printed bases per row
        """
        if num_seqs:
            if not isinstance(num_seqs, int):
                raise TypeError("num_seqs is not an integer")
            self._repr_policy["num_seqs"] = num_seqs

        if num_pos:
            if not isinstance(num_pos, int):
                raise TypeError("num_pos is not an integer")
            self._repr_policy["num_pos"] = num_pos

        if ref_name:
            if not isinstance(ref_name, str):
                raise TypeError("ref_name is not a string")

            if ref_name != "longest" and ref_name not in self.names:
                raise ValueError(f"no sequence name matching {ref_name}")

            self._repr_policy["ref_name"] = ref_name

        if wrap:
            if not isinstance(wrap, int):
                raise TypeError("wrap is not an integer")
            self._repr_policy["wrap"] = wrap


@singledispatch
def seq_to_gap_coords(
    seq: Union[str, numpy.ndarray], moltype: new_moltype.MolType
) -> tuple[str, IndelMap]:
    """
    Takes a sequence with (or without) gaps and returns an ungapped sequence
    and records the position and length of gaps in the original parent sequence
    """
    raise NotImplementedError(f"{seq} not implemented for type {type(seq)}")


@seq_to_gap_coords.register
def _(seq: str, moltype: new_moltype.MolType) -> tuple[str, IndelMap]:
    seq = moltype.make_seq(seq=seq)
    indel_map, ungapped_seq = seq.parse_out_gaps()

    if indel_map.num_gaps == 0:
        return str(ungapped_seq), numpy.array([], dtype=int)

    return str(ungapped_seq), indel_map


@seq_to_gap_coords.register
def _(seq: numpy.ndarray, moltype: new_moltype.MolType) -> tuple[str, IndelMap]:
    gap_char = moltype.degen_gapped_alphabet.index(moltype.gap)
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
    _moltype: new_moltype.MolType = field(init=False)
    _names: tuple[str] = field(init=False)
    _alpha: new_alpha.CharAlphabet = field(init=False)
    align_len: int = 0
    moltype: InitVar[Union[str, None]] = "dna"
    names: InitVar[Union[tuple[str], None]] = None

    def __post_init__(self, moltype, names):
        self._moltype = new_moltype.get_moltype(moltype)
        self._alpha = self._moltype.degen_gapped_alphabet
        self.seqs = {k: self._alpha.to_indices(v) for k, v in self.seqs.items()}

    @classmethod
    def from_gapped_seqs(
        cls,
        data: dict[str, Union[str, numpy.ndarray]],
        moltype: Union[str, new_moltype.MolType] = "dna",
        names: Optional[tuple[str]] = None,
    ):
        """
        Convert dict of {"seq_name": "seq"} to two dicts for seqs and gaps
        """
        seq_lengths = {len(v) for v in data.values()}
        if len(seq_lengths) != 1:
            raise ValueError("All sequence lengths must be the same.")

        align_len = seq_lengths.pop()
        moltype = new_moltype.get_moltype(moltype)

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
        return AlignedDataView(seqs=self, seqid=seqid, seq_len=self.align_len)

    def get_gaps(self, seqid: str) -> numpy.ndarray:
        return self.gaps[seqid]
