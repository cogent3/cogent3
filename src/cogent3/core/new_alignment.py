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
    FeatureDataType,
    SupportsFeatures,
)
from cogent3.core.info import Info as InfoClass
from cogent3.core.location import IndelMap
from cogent3.format.fasta import alignment_to_fasta
from cogent3.format.phylip import alignment_to_phylip
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
            seqs=self.seqs, seqid=self.seqid, seq_len=self._seq_len, start=0, stop=0
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
    make_seq: Callable = (
        None  # refactor: type hint for callable should specificy input/return type
    )

    @abstractmethod
    def seq_lengths(self) -> dict[str, int]: ...

    @property
    @abstractmethod
    def names(self) -> list: ...

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
    def add_seqs(self, seqs) -> SeqsDataABC: ...

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
        # todo: kath, make underlying data unchangable with `arrl.flags.writeable = False`
        self._data: dict[str, numpy.ndarray] = {
            str(name): self._alphabet.to_indices(seq) for name, seq in data.items()
        }

    @property
    def names(self) -> list:
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
        if data := {name: self._data.get(name) for name in names if name in self.names}:
            return self.__class__(
                data=data, alphabet=self.alphabet, make_seq=self.make_seq
            )
        else:
            raise ValueError(f"provided {names=} not found in collection")

    def add_seqs(
        self, seqs: dict[str, PrimitiveSeqTypes], force_unique_keys=True
    ) -> SeqsData:
        # todo: kath
        ...

    def to_alphabet(
        self, alphabet: new_alpha.CharAlphabet, check_valid=True
    ) -> SeqsData:
        # refactor: design -- map directly between arrays?
        # refactor: design -- better way to check if alphabets are DNA and RNA?
        if len(alphabet) == len(self.alphabet):
            return self.__class__(data=self._data, alphabet=alphabet)

        new_data = {}
        for seqid in self.names:
            seq_data = self.get_seq_array(seqid=seqid)
            old = self.alphabet.as_bytes()
            new = alphabet.as_bytes()

            convert_old_to_bytes = new_alpha.array_to_bytes(old)
            convert_bytes_to_new = new_alpha.bytes_to_array(
                new, dtype=new_alpha.get_array_type(len(new))
            )

            as_new_alpha = convert_bytes_to_new(convert_old_to_bytes(seq_data))

            if check_valid and not alphabet.is_valid(as_new_alpha):
                raise ValueError(
                    f"Changing from old alphabet={self.alphabet} to new {alphabet=} is not valid for this data"
                )

            new_data[seqid] = as_new_alpha

        return self.__class__(data=new_data, alphabet=alphabet)

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

    def get_seq(self, seqname: str, copy_annotations: bool = False):
        """Return a sequence object for the specified seqname."""
        seq = self.seqs[seqname]

        if copy_annotations:
            seq.annotation_db = type(self.annotation_db)()
            seq.annotation_db.update(annot_db=self.annotation_db, seqids=seqname)

        return seq

    def add_seqs(self, seqs: dict[str, PrimitiveSeqTypes], **kwargs):
        # todo: kath, add this method once its functional on SeqsData
        ...

    def to_dict(self, as_array: bool = False) -> dict[str, Union[str, numpy.ndarray]]:
        """Return a dictionary of sequences.

        Parameters
        ----------
        as_array
            if True, sequences are returned as numpy arrays, otherwise as strings
        """

        get = self.seqs.get_seq_array if as_array else self.seqs.get_seq_str
        return {name: get(seqid=name) for name in self.names}

    def degap(self, **kwargs):
        """Returns copy in which sequences have no gaps.

        Parameters
        ----------
        kwargs
            passed to class constructor
        """
        # refactor: array - chars > gap char
        seqs = {name: self.seqs[name].degap() for name in self.names}
        return make_unaligned_seqs(seqs, moltype=self.moltype, info=self.info, **kwargs)

    def get_lengths(
        self, include_ambiguity: bool = False, allow_gap: bool = False
    ) -> dict[str, int]:
        """returns sequence lengths as a dict of {seqid: length}

        Parameters
        ----------
        include_ambiguity
            if True, motifs containing ambiguous characters
            from the seq moltype are included. No expansion of those is attempted.
        allow_gap
            if True, motifs containing a gap character are included.

        """
        counts = self.counts_per_seq(
            motif_length=1, include_ambiguity=include_ambiguity, allow_gap=allow_gap
        )
        return counts.row_sum()

    def copy_annotations(self, seq_db: SupportsFeatures) -> None:
        """copy annotations into attached annotation db

        Parameters
        ----------
        seq_db
            compatible annotation db

        Notes
        -----
        Only copies annotations for records with seqid in self.names
        """
        if not isinstance(seq_db, SupportsFeatures):
            raise TypeError(
                f"type {type(seq_db)} does not match SupportsFeatures interface"
            )

        num = 0
        for seqid in self.names:
            num += seq_db.num_matches(seqid=seqid)
            if num > 0:
                break
        else:
            # no matching ID's, nothing to do
            return

        if self.annotation_db is None:
            self.annotation_db = type(seq_db)()

        if self.annotation_db.compatible(seq_db, symmetric=False):
            # our db contains the tables in other, so we update in place
            self.annotation_db.update(annot_db=seq_db, seqids=self.names)
        else:
            # we use the union method to define a new one
            # the setter handles propagation of the new instance to bound
            # sequences
            self.annotation_db = self.annotation_db.union(seq_db)

    def make_feature(
        self,
        *,
        feature: FeatureDataType,
    ) -> Feature:
        """
        create a feature on named sequence, or on the collection itself

        Parameters
        ----------
        feature
            a dict with all the necessary data to construct a feature

        Returns
        -------
        Feature

        Notes
        -----
        To get a feature AND add it to annotation_db, use add_feature().
        """
        return self.seqs[feature["seqid"]].make_feature(feature)

    def add_feature(
        self,
        *,
        seqid: str,
        biotype: str,
        name: str,
        spans: list[tuple[int, int]],
        parent_id: Optional[str] = None,
        strand: str = "+",
    ) -> Feature:
        """
        add feature on named sequence

        Parameters
        ----------
        seqid
            seq name to associate with
        parent_id
            name of the parent feature
        biotype
            biological type
        name
            feature name
        spans
            plus strand coordinates
        strand
            either '+' or '-'

        Returns
        -------
        Feature
        """
        if not self.annotation_db:
            self.annotation_db = DEFAULT_ANNOTATION_DB()

        if seqid and seqid not in self.names:
            raise ValueError(f"unknown {seqid=}")

        feature = {k: v for k, v in locals().items() if k != "self"}

        self.annotation_db.add_feature(**feature)
        feature.pop("parent_id", None)
        return self.make_feature(feature=feature)

    def get_features(
        self,
        *,
        seqid: Union[str, Iterator[str]] = None,
        biotype: Optional[str] = None,
        name: Optional[str] = None,
        allow_partial: bool = False,
    ) -> Iterator[Feature]:
        """yields Feature instances

        Parameters
        ----------
        seqid
            limit search to features on this named sequence, defaults to search all
        biotype
            biotype of the feature, e.g. CDS, gene
        name
            name of the feature
        allow_partial
            allow features partially overlaping self

        Notes
        -----
        When dealing with a nucleic acid moltype, the returned features will
        yield a sequence segment that is consistently oriented irrespective
        of strand of the current instance.
        """
        if not self.annotation_db:
            return None

        seqid_to_seqname = {seq.parent_coordinates()[0]: seq.name for seq in self.seqs}
        if seqid and (
            seqid not in seqid_to_seqname and seqid not in seqid_to_seqname.values()
        ):
            raise ValueError(f"unknown {seqid=}")

        for feature in self.annotation_db.get_features_matching(
            seqid=seqid,
            biotype=biotype,
            name=name,
            on_alignment=False,
            allow_partial=allow_partial,
        ):
            seqname = seqid_to_seqname[feature["seqid"]]
            seq = self.seqs[seqname]
            if offset := seq.annotation_offset:
                feature["spans"] = (numpy.array(feature["spans"]) - offset).tolist()
            yield seq.make_feature(feature, self)

    def to_fasta(self, block_size: int = 60) -> str:
        """Return collection in Fasta format.

        Parameters
        ----------
        block_size
            the sequence length to write to each line,
            by default 60

        Returns
        -------
        The collection in Fasta format.
        """
        return alignment_to_fasta(self.to_dict(), block_size=block_size)

    def to_phylip(self):
        """
        Return collection in PHYLIP format and mapping to sequence ids

        Notes
        -----
        raises exception if invalid sequences do not all have the same length
        """
        if len(set(self.seqs.seq_lengths().values())) > 1:
            raise ValueError("not all seqs same length, cannot convert to phylip")

        return alignment_to_phylip(self.to_dict())

    def has_terminal_stop(self, gc: typing.Any = None, strict: bool = False) -> bool:
        """Returns True if any sequence has a terminal stop codon.

        Parameters
        ----------
        gc
            valid input to cogent3.get_code(), a genetic code object, number
            or name
        strict
            If True, raises an exception if a seq length not divisible by 3
        """

        for seq_name in self.names:
            seq = self.seqs[seq_name]
            if seq.has_terminal_stop(gc=gc, strict=strict):
                return True
        return False

    def get_identical_sets(self, mask_degen: bool = False):  # refactor: array
        """returns sets of names for sequences that are identical

        Parameters
        ----------
        mask_degen
            if True, degenerate characters are ignored

        """

        def reduced(seq, indices):
            return "".join(seq[i] for i in range(len(seq)) if i not in indices)

        identical_sets = []
        mask_posns = {}
        seen = []
        # if strict, we do a sort and one pass through the list
        seqs = self.to_dict()
        if not mask_degen:
            seqs_names = [(s, n) for n, s in seqs.items()]
            seqs_names.sort()
            matched = None
            dupes = defaultdict(set)
            for i in range(len(seqs_names) - 1):
                if seqs_names[i][0] == seqs_names[i + 1][0]:
                    matched = seqs_names[i][1] if matched is None else matched
                    dupes[matched].update([seqs_names[i + 1][1], matched])
                else:
                    matched = None
            identical_sets = list(dupes.values())
            return identical_sets

        for i in range(len(self.names) - 1):
            n1 = self.names[i]
            if n1 in seen:
                continue

            seq1 = seqs[n1]
            if n1 not in mask_posns:
                pos = self.moltype.get_degenerate_positions(seq1, include_gap=True)
                mask_posns[n1] = pos

            group = set()
            for j in range(i + 1, len(self.names)):
                n2 = self.names[j]
                if n2 in seen:
                    continue

                seq2 = seqs[n2]
                if n2 not in mask_posns:
                    pos = self.moltype.get_degenerate_positions(seq2, include_gap=True)
                    mask_posns[n2] = pos

                seq2 = seqs[n2]

                s1, s2 = seq1, seq2

                pos = mask_posns[n1] + mask_posns[n2]
                if pos:
                    s1 = reduced(seq1, pos)
                    s2 = reduced(seq2, pos)

                if s1 == s2:
                    seen.append(n2)
                    group.update([n1, n2])

            if group:
                identical_sets.append(group)

        return identical_sets

    def get_similar(
        self,
        target: new_seq.Sequence,
        min_similarity: float = 0.0,
        max_similarity: float = 1.0,
        metric: float = new_seq.frac_same,
        transform: bool = None,
    ) -> new_seq.SequenceCollection:
        """Returns new SequenceCollection containing sequences similar to target.

        Parameters
        ----------
        target
            sequence object to compare to. Can be in the collection.
        min_similarity
            minimum similarity that will be kept. Default 0.0.
        max_similarity
            maximum similarity that will be kept. Default 1.0.
            (Note that both min_similarity and max_similarity are inclusive.)
        metric
            similarity function to use. Must be f(first_seq, second_seq).
            The default metric is fraction similarity, ranging from 0.0 (0%
            identical) to 1.0 (100% identical). The Sequence classes have lots
            of methods that can be passed in as unbound methods to act as the
            metric, e.g. frac_same_gaps.
        transform
            transformation function to use on the sequences before
            the metric is calculated. If None, uses the whole sequences in each
            case. A frequent transformation is a function that returns a specified
            range of a sequence, e.g. eliminating the ends. Note that the
            transform applies to both the real sequence and the target sequence.

        WARNING: if the transformation changes the type of the sequence (e.g.
        extracting a string from an RnaSequence object), distance metrics that
        depend on instance data of the original class may fail.
        """
        if transform:
            target = transform(target)

        def m(x):
            return metric(target, x)

        if transform:

            def f(x):
                result = m(transform(x))
                return min_similarity <= result <= max_similarity

        else:

            def f(x):
                result = m(x)
                return min_similarity <= result <= max_similarity

        return self.take_seqs_if(f)

    def __str__(self):
        """Returns self in FASTA-format, respecting name order."""
        from cogent3.format.alignment import FORMATTERS

        return FORMATTERS["fasta"](self.to_dict())

    def __eq__(
        self, other: Union[new_seq.SequenceCollection, dict]
    ) -> bool:  # refactor: design
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
            number of columns per row
        limit
            truncate view of collection to this length
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

        seq_lengths = numpy.array(list(self.seqs.seq_lengths().values()))
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
def merged_db_collection(seqs) -> SupportsFeatures:
    """return one AnnotationDb from a collection of sequences

    Parameters
    ----------
    seqs
        iterable list of data

    Returns
    -------
    list of all annotation db's

    Raises
    ------
    TypeError if different classes of AnnotationDb
    """
    first = None
    merged = None
    for seq in seqs:
        if not isinstance(seq, new_seq.Sequence):
            continue

        db = seq.annotation_db

        if first is None and db:
            # todo gah should this be a copy so immutable?
            first = db
            merged = first
            continue

        if first is None or db is None or first is db:
            continue
        first.update(db)

    return merged


@merged_db_collection.register
def _(seqs: dict) -> SupportsFeatures:
    return merged_db_collection(seqs.values())


@singledispatch
def coerce_to_seqs_data_dict(
    data, label_to_name: Callable = None
) -> dict[
    str, Union[PrimitiveSeqTypes, new_seq.Sequence]
]:  # refactor: type hint for callable should specificy input/return type
    raise NotImplementedError(
        f"coerce_to_seqs_data_dict not implemented for {type(data)}"
    )


@coerce_to_seqs_data_dict.register
def _(data: dict, label_to_name: Callable = None) -> dict[str, PrimitiveSeqTypes]:
    is_sequence = isinstance(next(iter(data.values()), None), new_seq.Sequence)
    return {
        (label_to_name(k) if label_to_name else k): (str(v) if is_sequence else v)
        for k, v in data.items()
    }


@coerce_to_seqs_data_dict.register
def _(data: list, label_to_name: Callable = None) -> dict[str, PrimitiveSeqTypes]:
    first = data[0]
    labelled_seqs = assign_names(first, data=data)
    return coerce_to_seqs_data_dict(labelled_seqs, label_to_name=label_to_name)


@coerce_to_seqs_data_dict.register
def _(data: set, label_to_name: Callable = None) -> dict[str, PrimitiveSeqTypes]:
    first = next(iter(data))
    labelled_seqs = assign_names(first, data=data)
    return coerce_to_seqs_data_dict(labelled_seqs, label_to_name=label_to_name)


@singledispatch
def assign_names(
    first, data: Union[list, set]
) -> dict[str, Union[PrimitiveSeqTypes, new_seq.Sequence]]:
    if isinstance(first, (str, bytes, numpy.ndarray)):
        names = assign_sequential_names(len(data))
        return dict(zip(names, data))
    raise NotImplementedError(f"assign_names not implemented for {type(data)}")


@assign_names.register
def _(first: new_seq.Sequence, data: Union[list, set]) -> dict[str, new_seq.Sequence]:
    return {seq.name: seq for seq in data}


@assign_names.register
def _(first: list, data: Union[list, set]) -> dict[str, PrimitiveSeqTypes]:
    return {name_seq[0]: name_seq[1] for name_seq in data}


@singledispatch
def make_unaligned_seqs(
    data: Union[dict[str, PrimitiveSeqTypes], SeqsData, list],
    *,
    moltype: Union[str, new_moltype.MolType],
    label_to_name: Callable = None,
    info: dict = None,
    source: Union[str, Path] = None,
    annotation_db: SupportsFeatures = None,
) -> SequenceCollection:  # refactor: design
    """Initialise an unaligned collection of sequences.

    Parameters
    ----------
    data
        sequence data, a SeqsData, a dict {name: seq, ...}, an iterable of sequences
    moltype
        string representation of the moltype, e.g., 'dna', 'protein'.
    label_to_name
        function for converting original names into other names.
    info
        a dict from which to make an info object
    source
        origins of this data, defaults to 'unknown'. Converted to a string
        and added to info["source"].
    """
    if len(data) == 0:
        raise ValueError("data must be at least one sequence.")

    # create a db from sequences if they're annotated
    seqs_anno_db = merged_db_collection(data)
    # if db from seqs and provided db, merge
    # todo: kath, should we default to not merging unless specified?
    if annotation_db and seqs_anno_db:
        annotation_db.update(seqs_anno_db)
    else:
        annotation_db = annotation_db or seqs_anno_db

    seqs_data = coerce_to_seqs_data_dict(data, label_to_name=label_to_name)

    moltype = new_moltype.get_moltype(moltype)
    seqs_data = SeqsData(data=seqs_data, alphabet=moltype.degen_gapped_alphabet)
    return make_unaligned_seqs(
        seqs_data,
        moltype=moltype,
        label_to_name=label_to_name,
        info=info,
        source=source,
        annotation_db=annotation_db,
    )


@make_unaligned_seqs.register
def _(
    data: SeqsData,
    *,
    moltype: Union[str, new_moltype.MolType],
    label_to_name: Callable = None,
    info: dict = None,
    source: Union[str, Path] = None,
    annotation_db: SupportsFeatures = None,
) -> SequenceCollection:

    moltype = new_moltype.get_moltype(moltype)
    if not moltype.is_compatible_alphabet(data.alphabet):
        raise ValueError(
            f"Provided moltype: {moltype} is not compatible with SeqsData alphabet {data.alphabet}"
        )

    info = info if isinstance(info, dict) else {}
    info["source"] = str(info.get("source", "unknown"))

    return SequenceCollection(
        seqs_data=data, moltype=moltype, info=info, annotation_db=annotation_db
    )


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
