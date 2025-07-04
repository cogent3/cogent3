from __future__ import annotations

import collections
import contextlib
import copy
import dataclasses
import hashlib
import json
import pathlib
import re
import types
import typing
import warnings
from abc import ABC, abstractmethod
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Iterator, Mapping
from collections.abc import Sequence as PySeq
from functools import singledispatch, singledispatchmethod
from pathlib import Path
from typing import Optional

import numba
import numpy
import typing_extensions

import cogent3
from cogent3._version import __version__
from cogent3.core import (
    alphabet as c3_alphabet,
)
from cogent3.core import (
    genetic_code as c3_genetic_code,
)
from cogent3.core import (
    moltype as c3_moltype,
)
from cogent3.core import (
    sequence as c3_sequence,
)
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import (
    AnnotatableMixin,
    FeatureDataType,
    SupportsFeatures,
)
from cogent3.core.info import Info as InfoClass
from cogent3.core.location import (
    FeatureMap,
    IndelMap,
    Strand,
)
from cogent3.core.profile import PSSM, MotifCountsArray, MotifFreqsArray, load_pssm
from cogent3.maths.stats.number import CategoryCounter
from cogent3.util import progress_display as UI
from cogent3.util import warning as c3warn
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.dict_array import DictArray, DictArrayTemplate
from cogent3.util.io import atomic_write, get_format_suffixes
from cogent3.util.misc import (
    extend_docstring_from,
    get_object_provenance,
    get_setting_from_environ,
    negate_condition,
)
from cogent3.util.union_dict import UnionDict

# DESIGN NOTES
# the sequence data collections (SeqsDataABC and AlignedSeqsDataABC)
# have no concept of strand. All transformations with respect to strand
# are applied by the sequence record objects that have a .moltype
# attribute, i.e. Sequence and Aligned.
# both collections can indicate sequences that are reverse complemented
# by providing their names to the reversed_seqs argument.

if typing.TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.tree import PhyloNode
    from cogent3.evolve.fast_distance import DistanceMatrix
    from cogent3.maths.stats.contingency import TestResult


OptInt = int | None
OptFloat = float | None
OptStr = str | None
OptList = list | None
OptIterStr = Iterable[str] | None
PySeqStr = PySeq[str]
OptPySeqStr = PySeqStr | None
OptDict = dict | None
OptBool = bool | None
OptSliceRecord = c3_sequence.SliceRecord | None
DictStrStr = dict[str, str]
DictStrInt = dict[str, int]
OptCallable = Callable | None
OptRenamerCallable = Callable[[str], str] | None
OptPathType = str | Path | None
StrORArray = str | numpy.ndarray[int]
StrORBytesORArray = str | bytes | numpy.ndarray[int]
StrORBytesORArrayOrSeq = str | bytes | numpy.ndarray[int] | c3_sequence.Sequence
MolTypes = str | c3_moltype.MolType

# small number: 1-EPS is almost 1, and is used for things like the
# default number of gaps to allow in a column.
EPS = 1e-6


def array_hash64(data: numpy.ndarray) -> int:
    """returns 64-bit hash of numpy array.

    Notes
    -----
    This function does not introduce randomisation and so
    is reproducible between processes.
    """
    h = hashlib.md5(data.tobytes(), usedforsecurity=False)
    return int.from_bytes(h.digest()[:8], byteorder="little")


class _SeqNamer:
    def __init__(
        self,
        name_func: OptRenamerCallable = None,
        base_name: str = "seq",
        start_at: int = 0,
    ) -> None:
        self._base_name = base_name
        self._num = start_at
        self._name_func = name_func

    def __call__(self, seq: StrORBytesORArray, name: OptStr = None) -> str:
        name = name or getattr(seq, "name", name)

        if not name:
            name = f"{self._base_name}_{self._num}"
            self._num += 1
        elif self._name_func:
            name = self._name_func(name)

        return name


@numba.njit(cache=True)
def _gap_ok_vector_single(
    data: numpy.ndarray[numpy.uint8],
    gap_index: int,
    missing_index: int,
    num_allowed: int,
) -> bool:  # pragma: no cover
    """returns indicies for which the number of gaps & missing data is less than or equal to num_allowed"""
    num = 0
    for i in range(len(data)):
        if data[i] == gap_index or (
            missing_index is not None and data[i] == missing_index
        ):
            num += 1

        if num > num_allowed:
            break

    return num <= num_allowed


@numba.njit(cache=True)
def _gap_ok_vector_multi(
    motifs: numpy.ndarray,
    gap_index: int,
    missing_index: int,
    motif_length: int,
    num_allowed: int,
) -> numpy.ndarray[bool]:  # pragma: no cover
    """returns indicies for which the number of gaps & missing data in a vector of motifs is less than or equal to num_allowed"""
    num = 0
    for motif in motifs:
        for j in range(motif_length):
            if motif[j] == gap_index or (
                missing_index is not None and motif[j] == missing_index
            ):
                num += 1
                break

        if num > num_allowed:
            break

    return num <= num_allowed


class SeqDataView(c3_sequence.SeqView):
    """
    A view class for ``SeqsData``, providing methods for different
    representations of a single sequence.

    Notes
    -----
    ``str_value`` / ``array_value`` are not complemented, but can be reversed.
    The latter is done by the ``Sequence`` object which has a moltype.
    """

    __slots__ = ("_parent_len", "_seqid", "_slice_record", "alphabet", "parent")

    @property
    def str_value(self) -> str:
        """returns the sequence as a string"""
        return self.alphabet.from_indices(self.array_value)

    @property
    def array_value(self) -> numpy.ndarray:
        """returns the sequence as a numpy array"""
        # we select the data using plus strand coords
        raw = self.parent.get_seq_array(
            seqid=self.seqid,
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
    def copy(self, sliced: bool = False) -> typing_extensions.Self:
        """returns copy"""
        return self

    def to_rich_dict(self) -> dict[str, str | dict[str, str]]:
        """returns a json serialisable dict.

        Notes
        -----
        This method will slice the underlying sequence to the start and stop values

        Warnings
        --------
        This method is not intended to provide serialisation of this object,
        instead, it is intended for usage by an enclosing class.
        """

        data = {"type": get_object_provenance(self), "version": __version__}
        data["init_args"] = self._get_init_kwargs()

        if self.slice_record.is_reversed:
            adj = self.parent_len + 1
            start, stop = self.slice_record.stop + adj, self.slice_record.start + adj
        else:
            start, stop = self.slice_record.start, self.slice_record.stop

        data["init_args"]["parent"] = self.str_value[start:stop]
        new_sr = c3_sequence.SliceRecord(
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


class SeqsDataABC(ABC):
    """Abstract base class for respresenting the storage object for sequences underlying
    a SequenceCollection.
    """

    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        *,
        data: dict[str, StrORBytesORArray],
        alphabet: c3_alphabet.AlphabetABC,
        offset: dict[str, int] | None = None,
        check: bool = True,
        reversed_seqs: set[str] | None = None,
    ) -> None: ...

    @classmethod
    @abstractmethod
    def from_seqs(
        cls,
        *,
        data: dict[str, StrORBytesORArray],
        alphabet: c3_alphabet.AlphabetABC,
        **kwargs,
    ): ...

    @abstractmethod
    def __eq__(self, value: object) -> bool: ...

    @abstractmethod
    def __ne__(self, value: object) -> bool: ...

    @abstractmethod
    def get_seq_length(self, seqid: str) -> int: ...

    @property
    @abstractmethod
    def reversed_seqs(self) -> frozenset[str]: ...

    @property
    @abstractmethod
    def names(self) -> tuple[str, ...]: ...

    @property
    @abstractmethod
    def alphabet(self) -> c3_alphabet.CharAlphabet: ...

    @property
    @abstractmethod
    def offset(self) -> dict[str, int]: ...

    @abstractmethod
    def get_seq_array(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
    ) -> numpy.ndarray: ...

    @abstractmethod
    def get_seq_str(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
    ) -> str: ...

    @abstractmethod
    def get_seq_bytes(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
    ) -> bytes: ...

    @abstractmethod
    def get_view(self, seqid: str) -> c3_sequence.SeqViewABC: ...

    @abstractmethod
    def to_alphabet(self, alphabet: c3_alphabet.AlphabetABC) -> SeqsDataABC: ...

    @abstractmethod
    def add_seqs(self, seqs, **kwargs) -> SeqsDataABC: ...

    @abstractmethod
    def __len__(self) -> int: ...

    @abstractmethod
    def __getitem__(
        self,
        index: str | int,
    ) -> c3_sequence.Sequence | c3_sequence.SeqViewABC: ...

    @abstractmethod
    def copy(self, **kwargs) -> SeqsDataABC: ...


class SeqsData(SeqsDataABC):
    """The builtin ``cogent3`` implementation of sequence storage underlying
    a ``SequenceCollection``. The sequence data is stored as numpy arrays. Indexing
    this object (using an int or seq name) returns a ``SeqDataView``, which can realise
    the corresponding slice as a string, bytes, or numpy array via the alphabet.

    Notes
    -----
    Methods on this object only accepts plust strand start, stop and step
    indices for selecting segments of data. It can return the gap coordinates
    for a sequence as used by IndelMap.
    """

    __slots__ = ("_alphabet", "_data", "_offset", "_reversed")

    def __init__(
        self,
        *,
        data: dict[str, StrORBytesORArray],
        alphabet: c3_alphabet.CharAlphabet,
        offset: dict[str, int] | None = None,
        check: bool = True,
        reversed_seqs: set[str] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        data
            raw data as {seq name: sequence, ...} where the sequence can be converted
            to a numpy array using the provided alphabet.
        alphabet
            a cogent3 CharAlphabet instance, typically defined as
            <moltype>.most_degen_alphabet()
        offset
            dict indicating annotation offsets for each sequence
        check
            use the alphabet to check the sequences are valid
        reversed_seqs
            names of seqs that are reverse complemented

        Raises
        ------
        c3_alphabet.AlphabetError if the check fails
        """
        self._alphabet = alphabet
        self._offset = offset or {}
        self._reversed = frozenset(reversed_seqs or set())
        if check:
            assert self._offset.keys() <= data.keys(), (
                "sequence name provided in offset not found in data"
            )
            if any(not alphabet.is_valid(seq) for seq in data.values()):
                msg = f"One or more sequences are invalid for alphabet {alphabet}"
                raise c3_alphabet.AlphabetError(
                    msg,
                )
        self._data: dict[str, numpy.ndarray] = {}
        for name, seq in data.items():
            arr = self._alphabet.to_indices(seq)
            arr.flags.writeable = False
            self._data[str(name)] = arr

    def __eq__(self, other: SeqsDataABC) -> bool:
        if not isinstance(other, self.__class__):
            return False
        for attr_name in ("_alphabet", "_offset"):
            self_attr = getattr(self, attr_name)
            other_attr = getattr(other, attr_name)
            if self_attr != other_attr:
                return False

        # compare individuals sequences
        if self._data.keys() != other._data.keys():
            return False
        return all(
            numpy.array_equal(self._data[name], other._data[name])
            for name in self._data
        )

    def __ne__(self, other: object) -> bool:
        return not self == other

    @classmethod
    def from_seqs(
        cls,
        *,
        data: dict[str, StrORBytesORArray],
        alphabet: c3_alphabet.AlphabetABC,
        **kwargs,
    ):
        return cls(data=data, alphabet=alphabet, **kwargs)

    @property
    def names(self) -> tuple[str, ...]:
        """returns the names of the sequences in the storage"""
        return tuple(self._data.keys())

    @property
    def reversed_seqs(self) -> frozenset[str]:
        """names of sequences that are reverse complemented"""
        return self._reversed

    @property
    def alphabet(self) -> c3_alphabet.CharAlphabet:
        """the character alphabet for validating, encoding, decoding sequences"""
        return self._alphabet

    @property
    def offset(self) -> dict[str, int]:
        """annotation offsets for each sequence"""
        return {name: self._offset.get(name, 0) for name in self.names}

    def get_seq_length(self, seqid: str) -> int:
        """return length for seqid"""
        return self._data[seqid].shape[0]

    def get_seq_array(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray:
        start = start or 0
        stop = stop if stop is not None else self.get_seq_length(seqid)
        step = step or 1

        if start < 0 or stop < 0 or step < 1:
            msg = f"{start=}, {stop=}, {step=} not >= 1"
            raise ValueError(msg)

        out_len = (stop - start + step - 1) // step
        out = numpy.empty(out_len, dtype=self.alphabet.dtype)

        out[:] = self._data[seqid][start:stop:step]
        return out

    def get_seq_str(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> str:
        return self._alphabet.from_indices(
            self.get_seq_array(seqid=seqid, start=start, stop=stop, step=step),
        )

    def get_seq_bytes(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> bytes:
        return self.get_seq_str(seqid=seqid, start=start, stop=stop, step=step).encode(
            "utf8",
        )

    def get_view(self, seqid: str) -> SeqDataView:
        """reurns view of sequence data for seqid"""
        seq_len = len(self._data[seqid])
        return SeqDataView(
            parent=self,
            seqid=seqid,
            parent_len=seq_len,
            alphabet=self.alphabet,
            offset=self._offset.get(seqid, 0),
        )

    def add_seqs(
        self,
        seqs: dict[str, StrORBytesORArray],
        force_unique_keys: bool = True,
        offset: dict[str, int] | None = None,
    ) -> SeqsData:
        """Returns a new SeqsData object with added sequences. If force_unique_keys
        is True, raises ValueError if any names already exist in the collection."""
        if force_unique_keys and any(name in self.names for name in seqs):
            msg = "One or more sequence names already exist in collection"
            raise ValueError(msg)
        new_data = {
            **self._data,
            **{name: self.alphabet.to_indices(seq) for name, seq in seqs.items()},
        }
        return self.copy(
            data=new_data,
            alphabet=self.alphabet,
            offset={**self._offset, **(offset or {})},
        )

    def to_alphabet(
        self,
        alphabet: c3_alphabet.AlphabetABC,
        check_valid: bool = True,
    ) -> SeqsData:
        if (
            len(self.alphabet) == len(alphabet)
            and len(
                {
                    (a, b)
                    for a, b in zip(self.alphabet, alphabet, strict=False)
                    if a != b
                },
            )
            == 1
        ):
            # rna <-> dna swap just replace alphabet
            return self.copy(alphabet=alphabet)

        new_data = {}
        for seqid in self.names:
            seq_data = self.get_seq_array(seqid=seqid)
            as_new_alpha = self.alphabet.convert_seq_array_to(
                seq=seq_data,
                alphabet=alphabet,
                check_valid=False,
            )
            if check_valid and not alphabet.is_valid(as_new_alpha):
                msg = (
                    f"Changing from old alphabet={self.alphabet} to new "
                    f"{alphabet=} is not valid for this data"
                )
                raise c3_alphabet.AlphabetError(msg)
            new_data[seqid] = as_new_alpha

        return self.copy(
            data=new_data,
            alphabet=alphabet,
            check=False,
        )

    def __len__(self) -> int:
        return len(self.names)

    @singledispatchmethod
    def __getitem__(self, index: str | int) -> c3_sequence.SeqViewABC:
        msg = f"__getitem__ not implemented for {type(index)}"
        raise NotImplementedError(msg)

    @__getitem__.register
    def _(self, index: str) -> c3_sequence.SeqViewABC:
        return self.get_view(seqid=index)

    @__getitem__.register
    def _(self, index: int) -> c3_sequence.SeqViewABC:
        return self[self.names[index]]

    def copy(self, **kwargs) -> typing_extensions.Self:
        """shallow copy of self

        Notes
        -----
        kwargs are passed to constructor and will over-ride existing values
        """
        init_args = {
            "data": self._data,
            "alphabet": self._alphabet,
            "offset": self._offset,
            "reversed_seqs": self._reversed,
            **kwargs,
        }
        return self.__class__(**init_args)


class SequenceCollection(AnnotatableMixin):
    """A container of unaligned sequences.

    Notes
    -----
    Should be constructed using ``make_unaligned_seqs()``.
    """

    def __init__(
        self,
        *,
        seqs_data: SeqsDataABC,
        moltype: c3_moltype.MolType,
        info: dict | InfoClass | None = None,
        source: OptPathType = None,
        annotation_db: SupportsFeatures | None = None,
        name_map: OptDict = None,
        is_reversed: bool = False,
    ) -> None:
        """
        Parameters
        ----------
        seqs_data
            a SeqsDataABC instance containg the sequence data
        moltype
            the molecular type of the sequences
        info
            additional information about the collection
        source
            the source of the sequence data
        annotation_db
            a database of annotations for the sequences
        name_map
            map between the names specified in the collection and names used in
            the underlying seqs_data. Used for when the names have been changed,
            but we want to query for annotations using the original names.
        is_reversed
            entire collection is reverse complemented
        """
        self._seqs_data = seqs_data
        self.moltype = moltype
        name_map = name_map or {name: name for name in seqs_data.names}
        self._name_map = types.MappingProxyType(name_map)
        self._is_reversed = is_reversed
        if not isinstance(info, InfoClass):
            info = InfoClass(info) if info else InfoClass()
        self.info = info
        self.source = source
        self._repr_policy = {
            "num_seqs": 10,
            "num_pos": 60,
            "ref_name": "longest",
            "wrap": 60,
        }
        self._annotation_db: list[SupportsFeatures] = self._init_annot_db_value(
            annotation_db
        )
        self._seqs = None
        self._post_init()

    def _post_init(self) -> None:
        # override in subclasses
        self._seqs = _IndexableSeqs(self, make_seq=self._make_seq)

    def __getstate__(self) -> dict:
        return self._get_init_kwargs()

    def __setstate__(self, state: dict) -> None:
        state["name_map"] = types.MappingProxyType(state.pop("name_map"))
        obj = self.__class__(**state)

        self.__dict__.update(obj.__dict__)

    def _make_seq(self, name: str) -> c3_sequence.Sequence:
        # seqview is given the name of the parent (if different from the current name)
        # the sequence is given the current name
        seqid = self._name_map.get(name, name)
        sv = self._seqs_data.get_view(seqid)
        if self._is_reversed:
            sv = sv[::-1]
        return self.moltype.make_seq(
            seq=sv,
            name=name,
            annotation_db=self._annotation_db,
        )

    def _get_init_kwargs(self) -> dict:
        """dict of all the arguments needed to initialise a new instance"""
        # both SequenceCollection and Alignment implement _get_init_kwargs,
        # ensuring methods in SequenceCollection that are inherited by Alignment
        # capture initialisation arguments unique to the subclass.
        # mutable arguments are copied
        return {
            "seqs_data": self._seqs_data,
            "moltype": self.moltype,
            "name_map": dict(self._name_map),
            "info": self.info.copy(),
            "annotation_db": self._annotation_db,
            "source": self.source,
            "is_reversed": self._is_reversed,
        }

    @property
    def storage(self) -> SeqsDataABC:
        """the unaligned sequence storage instance of the collection"""
        return self._seqs_data

    @storage.setter
    def storage(self, value: object) -> None:
        # storage cannot be set after initialisation
        msg = "storage cannot be set after initialisation"
        raise TypeError(msg)

    @property
    def seqs(self) -> _IndexableSeqs:
        """iterable of sequences in the collection

        Notes
        -----
        Can be indexed by a sequence name or integer index.
        Cannot be sliced.

        Returns
        -------
        Instance of ``MolType`` sequence or ``Aligned`` sequence if
        self is an ``Alignment``.
        """
        return self._seqs

    @property
    def names(self) -> tuple[str, ...]:
        """returns the names of the sequences in the collection"""
        return tuple(self._name_map.keys())

    @property
    def name_map(self) -> types.MappingProxyType:
        """returns mapping of seq names to parent seq names

        Notes
        -----
        The underlying SeqsData may have different names for the same
        sequences. This object maps the names of sequences in self to
        the names of sequences in SeqsData.
        MappingProxyType is an immutable mapping, so it cannot be
        changed. Use self.rename_seqs() to do that.
        """
        return self._name_map

    @name_map.setter
    def name_map(self, value: OptDict) -> None:  # noqa: ARG002
        """name_map can only be set at initialisation"""
        msg = "name_map cannot be set after initialisation"
        raise TypeError(msg)

    @property
    def num_seqs(self) -> int:
        """the number of sequences in the collection"""
        return len(self.names)

    def iter_seqs(
        self,
        seq_order: OptList = None,
    ) -> Iterator[c3_sequence.Sequence | c3_sequence.SeqViewABC]:
        """Iterates over sequences in the collection, in order.

        Parameters
        ----------
        seq_order:
            list of seqids giving the order in which seqs will be returned.
            Defaults to self.names
        """
        if seq_order is None:
            yield from self.seqs
        else:
            for name in seq_order:
                yield self.seqs[name]

    def take_seqs(
        self,
        names: str | typing.Sequence[str],
        negate: bool = False,
        copy_annotations: bool = False,
        **kwargs,
    ) -> typing_extensions.Self:
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
        """

        # to return a new collection with a subset of the sequences we dont
        # want to modify the underlying data, instead we create a new collection
        # with a subset of the names, recorded in the name_map dict.

        # refactor: design, reimplement on Alignment. on which, if self.array_seqs
        # defined, assign result of self._array_seqs.take(subset_name_indices) to
        # resulting alignments _array_seqs attribute

        if isinstance(names, str):
            names = [names]

        if negate:
            names = [name for name in self.names if name not in names]

        if not names:
            msg = f"{names=} and {negate=} resulted in no names"
            raise ValueError(msg)

        if diff := set(names) - set(self.names):
            msg = f"The following provided names not found in collection: {diff}"
            raise ValueError(msg)

        selected_name_map = {name: self._name_map[name] for name in names}

        init_kwargs = self._get_init_kwargs()
        init_kwargs["name_map"] = selected_name_map
        if self._annotation_db:
            if copy_annotations:
                ann_db = type(self.annotation_db)()
                ann_db.update(
                    annot_db=self.annotation_db,
                    seqids=list(selected_name_map),
                )
            else:
                ann_db = self.annotation_db
            init_kwargs["annotation_db"] = ann_db

        return self.__class__(**init_kwargs)

    def get_seq_names_if(
        self,
        f: Callable[[c3_sequence.Sequence], bool],
        negate: bool = False,
    ) -> list[str]:
        """Returns list of names of seqs where f(seq) is True.

        Parameters
        ----------
        f
            function that takes a sequence object and returns True or False
        negate
            select all sequences EXCEPT those where f(seq) is True

        Notes
        -----
        Sequence objects can be converted into strings or numpy arrays using
        str() and numpy.array() respectively.
        """
        get = self.seqs

        new_f = negate_condition(f) if negate else f

        return [name for name in self.names if new_f(get[name])]

    def take_seqs_if(
        self,
        f: Callable[[c3_sequence.Sequence], bool],
        negate: bool = False,
    ) -> typing_extensions.Self:
        """Returns new collection containing seqs where f(seq) is True.

        Parameters
        ----------
        f
            function that takes a sequence object and returns True or False
        negate
            select all sequences EXCEPT those where f(seq) is True

        Notes
        -----
        Sequence objects can be converted into strings or numpy arrays using
        str() and numpy.array() respectively.
        """
        return self.take_seqs(self.get_seq_names_if(f, negate))

    def get_seq(
        self,
        seqname: str,
        copy_annotations: bool = False,
    ) -> c3_sequence.Sequence:
        """Return a Sequence object for the specified seqname.

        Parameters
        ----------
        seqname
            name of the sequence to return
        copy_annotations
            if True, only the annotations for the specified sequence are copied
            to the annotation database of the Sequence object which is decoupled
            from this collection. If False, the connection to this collections db
            is retained.
        """
        seq = self.seqs[seqname]
        if copy_annotations and self._annotation_db:
            # we need to copy the sequence too to break the link to self.annotation_db
            seq = seq.copy(exclude_annotations=True)
            seq.annotation_db = type(self.annotation_db)()
            seq.annotation_db.update(annot_db=self.annotation_db, seqids=seqname)
            return seq

        seq.annotation_db = self._annotation_db

        return seq

    def add_seqs(
        self,
        seqs: dict[str, StrORBytesORArray] | SeqsData | list,
        **kwargs,
    ) -> typing_extensions.Self:
        """Returns new collection with additional sequences.

        Parameters
        ----------
        seqs
            sequences to add
        """
        assign_names = _SeqNamer()
        seqs = _make_name_seq_mapping(seqs, assign_names)
        name_map = make_name_map(seqs)
        data, offsets, _ = prep_for_seqs_data(
            seqs,
            self.moltype,
            assign_names,
        )

        if not name_map:
            name_map = dict(zip(data, data, strict=False))

        kwargs["offset"] = offsets
        seqs_data = self._seqs_data.add_seqs(data, **kwargs)
        return self.__class__(
            seqs_data=seqs_data,
            moltype=self.moltype,
            name_map={**self._name_map, **name_map},
            info=self.info,
            source=self.source,
            annotation_db=self._annotation_db,
        )

    def rename_seqs(self, renamer: Callable[[str], str]) -> typing_extensions.Self:
        """Returns new collection with renamed sequences.

        Parameters
        ----------
        renamer
            callable that takes a name string and returns a string

        Raises
        ------
        ValueError if renamer produces duplicate names.

        Notes
        -----
        The resulting object stores the mapping of new to old names in
        self.name_map.
        """
        new_name_map = {}
        for name, parent_name in self._name_map.items():
            new_name = renamer(name)
            # we retain the parent_name when it differs from the name,
            # this can happen after multiple renames on the same collection
            parent_name = parent_name if name != parent_name else name  # noqa: PLW2901
            new_name_map[new_name] = parent_name

        if len(new_name_map) != len(self._name_map):
            msg = f"non-unique names produced by {renamer=}"
            raise ValueError(msg)

        init_args = self._get_init_kwargs()
        init_args["name_map"] = new_name_map

        return self.__class__(**init_args)

    def to_dict(self, as_array: bool = False) -> dict[str, str | numpy.ndarray]:
        """Return a dictionary of sequences.

        Parameters
        ----------
        as_array
            if True, sequences are returned as numpy arrays, otherwise as strings
        """
        return {s.name: (numpy.array(s) if as_array else str(s)) for s in self.seqs}

    def to_rich_dict(self) -> dict[str, str | dict[str, str]]:
        """returns a json serialisable dict

        Notes
        -----
        Deserialisation the object produced by this method will not include the
        annotation_db if present.
        """
        kwargs = self._get_init_kwargs()
        kwargs.pop("is_reversed", None)  # reversal is realised
        kwargs["moltype"] = self.moltype.label
        kwargs.pop("annotation_db", None)
        kwargs.pop(
            "offset",
            None,
        )  # no need for offset if we dont have an annotation_db
        kwargs.pop("seqs_data", None)  # we serialise the seqs_data directly

        data = {
            "init_args": kwargs,
            "type": get_object_provenance(self),
            "version": __version__,
        }
        data["seqs"] = {self._name_map[s.name]: str(s) for s in self.seqs}

        return data

    @classmethod
    def from_rich_dict(
        cls,
        data: dict[str, str | dict[str, str]],
    ) -> SequenceCollection:
        """returns a new instance from a rich dict"""
        return make_unaligned_seqs(data["seqs"], **data["init_args"])

    def to_json(self) -> str:
        """returns json formatted string"""
        return json.dumps(self.to_rich_dict())

    def degap(self, storage_backend: str | None = None, **kwargs) -> SequenceCollection:
        """returns collection sequences without gaps or missing characters.

        Parameters
        ----------
        storage_backend
            name of the storage backend to use for the SeqsData object, defaults to
            cogent3 builtin.
        kwargs
            keyword arguments for the storage driver

        Notes
        -----
        The returned collection will not retain an annotation_db if present.
        """
        if storage_backend:
            make_storage = cogent3._plugin.get_unaligned_storage_driver(  # noqa: SLF001
                storage_backend,
            ).from_seqs
        else:
            make_storage = self._seqs_data.from_seqs
        data = {}
        for name in self.names:
            # because we are in a SequenceCollection, which cannot be sliced, so
            # we can just interrogate the bound _seqs_data directly.
            seq = self._seqs_data.get_seq_array(seqid=self._name_map.get(name, name))
            data[self.name_map[name]] = self.moltype.degap(seq)

        init_kwargs = self._get_init_kwargs()
        init_kwargs["seqs_data"] = make_storage(
            data=data,
            alphabet=self._seqs_data.alphabet,
            check=False,
            **kwargs,
        )

        return self.__class__(**init_kwargs)

    def to_moltype(self, moltype: MolTypes) -> typing_extensions.Self:
        """returns copy of self with changed moltype

        Parameters
        ----------
        moltype
            name of the new moltype, e.g, 'dna', 'rna'.

        Notes
        -----
        Cannot convert from nucleic acids to proteins. Use get_translation() for that.

        """
        mtype = c3_moltype.get_moltype(moltype)
        if mtype is self.moltype:
            return self  # nothing to be done

        alpha = mtype.most_degen_alphabet()
        try:
            new_seqs_data = self._seqs_data.to_alphabet(alpha)
        except c3_moltype.MolTypeError as e:
            msg = f"Failed to convert moltype from {self.moltype.label} to {moltype}"
            raise c3_moltype.MolTypeError(
                msg,
            ) from e

        init_kwargs = self._get_init_kwargs()
        init_kwargs["seqs_data"] = new_seqs_data
        init_kwargs["moltype"] = mtype

        return self.__class__(**init_kwargs)

    def to_dna(self) -> typing_extensions.Self:
        """returns copy of self as a collection of DNA moltype seqs"""
        return self.to_moltype("dna")

    def to_rna(self) -> typing_extensions.Self:
        """returns copy of self as a collection of RNA moltype seqs"""
        return self.to_moltype("rna")

    def get_translation(
        self,
        gc: c3_genetic_code.GeneticCodeChoiceType = 1,
        incomplete_ok: bool = False,
        include_stop: bool = False,
        trim_stop: bool = True,
        **kwargs,
    ) -> typing_extensions.Self:
        """translate sequences from nucleic acid to protein

        Parameters
        ----------
        gc
            genetic code, either the number or name
            (use cogent3.core.genetic_code.available_codes)
        incomplete_ok
            codons that are mixes of nucleotide and gaps converted to '?'.
            raises a ValueError if False
        include_stop
            whether to allow a stops in the translated sequence
        trim_stop
            exclude terminal stop codons if they exist
        kwargs
            related to construction of the resulting object

        Returns
        -------
        A new instance of self translated into protein

        Notes
        -----
        Translating will break the relationship to an annotation_db if present.
        """
        if not self.moltype.is_nucleic:
            msg = f"moltype must be a DNA/RNA, not {self.moltype.name!r}"
            raise c3_moltype.MolTypeError(
                msg,
            )

        translated = {}
        for seq in self.seqs:
            pep = seq.get_translation(
                gc,
                incomplete_ok=incomplete_ok,
                include_stop=include_stop,
                trim_stop=trim_stop,
            )
            translated[self.name_map[seq.name]] = numpy.array(pep)

        pep_moltype = c3_moltype.get_moltype(
            "protein_with_stop" if include_stop else "protein",
        )
        seqs_data = self._seqs_data.from_seqs(
            data=translated,
            alphabet=pep_moltype.most_degen_alphabet(),
        )
        return self.__class__(
            seqs_data=seqs_data,
            moltype=pep_moltype,
            name_map=self._name_map,
            info=self.info,
            source=self.source,
            **kwargs,
        )

    def rc(self) -> typing_extensions.Self:
        """Returns the reverse complement of all sequences in the collection.
        A synonym for reverse_complement.
        """
        init_kwargs = self._get_init_kwargs()
        init_kwargs["is_reversed"] = not self._is_reversed
        return self.__class__(**init_kwargs)

    def reverse_complement(self) -> typing_extensions.Self:
        """Returns the reverse complement of all sequences in the collection.
        A synonym for rc.
        """
        return self.rc()

    def distance_matrix(self, calc: str = "pdist") -> DistanceMatrix:
        """Estimated pairwise distance between sequences

        Parameters
        ----------
        calc
            The distance calculation method to use, either "pdist" or "jc69".
            - "pdist" is an approximation of the proportion sites different.
            - "jc69" is an approximation of the Jukes Cantor distance.

        Returns
        -------
        DistanceMatrix
            Estimated pairwise distances between sequences in the collection

        Notes
        -----
        pdist approximates the proportion sites different from the Jaccard
        distance. Coefficients for the approximation were derived from a
        polynomial fit between Jaccard distance of kmers with k=10 and the
        proportion of sites different using mammalian 106 protein coding
        gene DNA sequence alignments.

        jc69 approximates the Jukes Cantor distance using the approximated
        proportion sites different, i.e., a transformation of the above.
        """
        from cogent3.app.dist import get_approx_dist_calc

        # check moltype
        if len(self.moltype.alphabet) != 4:
            msg = "only defined for DNA/RNA molecular types"
            raise NotImplementedError(msg)

        # assert we have more than one sequence in the SequenceCollection
        if self.num_seqs == 1:
            msg = (
                "Pairwise distance cannot be computed for a single sequence. "
                "Please provide at least two sequences."
            )
            raise ValueError(
                msg,
            )

        dist_calc_app = get_approx_dist_calc(
            dist=calc,
            num_states=len(self.moltype.alphabet),
        )

        return dist_calc_app(self)

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
            msg = f"type {type(seq_db)} does not match SupportsFeatures interface"
            raise TypeError(
                msg,
            )

        num = 0
        for seqid in self.names:
            num += seq_db.num_matches(seqid=seqid)
            if num > 0:
                break
        else:
            # no matching ID's, nothing to do
            return

        if not self._annotation_db:
            self.replace_annotation_db(type(seq_db)())

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
        parent_id: OptStr = None,
        strand: str | int = "+",
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
        if seqid and seqid not in self.names:
            msg = f"unknown {seqid=}"
            raise ValueError(msg)

        feature = {k: v for k, v in locals().items() if k != "self"}
        feature["strand"] = Strand.from_value(strand).value
        # property ensures db is created
        self.annotation_db.add_feature(**feature)
        feature.pop("parent_id", None)
        return self.make_feature(feature=feature)

    def get_features(
        self,
        *,
        seqid: str | Iterator[str] | None = None,
        biotype: OptStr = None,
        name: OptStr = None,
        start: OptInt = None,
        stop: OptInt = None,
        allow_partial: bool = False,
        **kwargs,
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
        start
            start position of the feature (not inclusive)
        stop
            stop position of the feature (inclusive)
        allow_partial
            allow features partially overlaping self
        kwargs
            additional keyword arguments to query the annotation db

        Notes
        -----
        - When dealing with a nucleic acid moltype, the returned features will
        yield a sequence segment that is consistently oriented irrespective
        of strand of the current instance.
        - start is non-inclusive, so if allow_partial is False, only features
        strictly starting after start will be returned.

        """

        if not self._annotation_db:
            return None

        if seqid and (seqid not in self.names):
            msg = f"unknown {seqid=}"
            raise ValueError(msg)

        for feature in self.annotation_db.get_features_matching(
            seqid=seqid,
            biotype=biotype,
            name=name,
            on_alignment=False,
            start=start,
            stop=stop,
            allow_partial=allow_partial,
            **kwargs,
        ):
            seqname = feature["seqid"]
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
        fasta = cogent3._plugin.get_seq_format_writer_plugin(format_name="fasta")  # noqa: SLF001
        return fasta.formatted(self, block_size=block_size)

    def write(self, filename: str, file_format: OptStr = None, **kwargs) -> None:
        """Write the sequences to a file, preserving order of sequences.

        Parameters
        ----------
        filename
            name of the sequence file
        file_format
            format of the sequence file

        Notes
        -----

        If file_format is None, will attempt to infer format from the filename
        suffix.
        """

        suffix, _ = get_format_suffixes(filename)
        if file_format is None and suffix:
            file_format = suffix

        if file_format == "json":
            with atomic_write(filename, mode="wt") as f:
                f.write(self.to_json())
            return

        if "order" not in kwargs:
            kwargs["order"] = self.names

        writer = cogent3._plugin.get_seq_format_writer_plugin(  # noqa: SLF001
            format_name=file_format,
            file_suffix=suffix,
            unaligned_seqs=type(self) == SequenceCollection,
        )
        _ = writer.write(seqcoll=self, path=filename, **kwargs)

    def dotplot(
        self,
        name1: OptStr = None,
        name2: OptStr = None,
        window: int = 20,
        threshold: OptInt = None,
        k: OptInt = None,
        min_gap: int = 0,
        width: int = 500,
        title: OptStr = None,
        rc: bool = False,
        biotype: str | tuple[str] | None = None,
        show_progress: bool = False,
    ):
        """make a dotplot between specified sequences. Random sequences
        chosen if names not provided.

        Parameters
        ----------
        name1, name2
            names of sequences -- if not provided, a random choice is made
        window
            segment size for comparison between sequences
        threshold
            windows where the sequences are identical >= threshold are a match
        k
            size of k-mer to break sequences into. Larger values increase
            speed but reduce resolution. If not specified, and
            window == threshold, then k is set to window. Otherwise, it is
            computed as the maximum of {threshold // (window - threshold), 5}.
        min_gap
            permitted gap for joining adjacent line segments, default is no gap
            joining
        width
            figure width. Figure height is computed based on the ratio of
            len(seq1) / len(seq2)
        title
            title for the plot
        rc
            include dotplot of reverse compliment also. Only applies to Nucleic
            acids moltypes
        biotype
            if selected sequences are annotated, display only these biotypes

        Returns
        -------
        a Drawable or AnnotatedDrawable
        """
        from cogent3.draw.dotplot import Dotplot
        from cogent3.draw.drawable import AnnotatedDrawable

        randgen = numpy.random.default_rng()
        if k is not None and not (0 < k < window):
            msg = f"{k=} must be positive, < {window=} and < {threshold=}"
            raise ValueError(msg)

        if len(self.names) == 1:
            name1 = name2 = self.names[0]
        elif name1 is None and name2 is None:
            name1, name2 = randgen.choice(self.names, size=2, replace=False).tolist()
        elif not (name1 and name2):
            names = list({*self.names, None} ^ {name1, name2})
            name = next(iter(randgen.choice(names, size=1))).item()
            name1 = name1 or name
            name2 = name2 or name

        if not {name1, name2} <= set(self.names):
            msg = f"{name1}, {name2} missing"
            raise ValueError(msg)

        seq1 = self.seqs[name1]
        seq2 = self.seqs[name2]
        if not self._annotation_db:
            annotated = False
        else:
            annotated = any(
                self.annotation_db.num_matches(seqid=self._name_map[n], biotype=biotype)
                for n in [name1, name2]
            )
        dotplot = Dotplot(
            seq1,
            seq2,
            isinstance(self, Alignment),
            window=window,
            threshold=threshold,
            k=k,
            min_gap=min_gap,
            xtitle=None if annotated else seq1.name,
            ytitle=None if annotated else seq2.name,
            title=title,
            moltype=self.moltype,
            rc=rc,
            show_progress=show_progress,
            width=width,
        )

        if annotated:
            data = getattr(seq1, "seq", seq1)
            bottom = data.get_drawable(biotype=biotype)
            data = getattr(seq2, "seq", seq2)
            left = data.get_drawable(biotype=biotype, vertical=True)
            dotplot = AnnotatedDrawable(
                dotplot,
                left_track=left,
                bottom_track=bottom,
                xtitle=seq1.name,
                ytitle=seq2.name,
                title=title,
                xrange=[0, len(seq1)],
                yrange=[0, len(seq2)],
            )

        return dotplot

    @UI.display_wrap
    def apply_pssm(
        self,
        pssm: PSSM = None,
        path: OptStr = None,
        background: numpy.ndarray = None,
        pseudocount: int = 0,
        names: OptList = None,
        ui=None,
    ) -> numpy.array:  # refactor: design: move to rich for progress bars?
        """scores sequences using the specified pssm

        Parameters
        ----------
        pssm :
            A profile.PSSM instance, if not provided, will be loaded from path
        path
            path to either a jaspar or cisbp matrix (path must end have a suffix
            matching the format).
        background
            background frequencies distribution
        pseudocount
            adjustment for zero in matrix
        names
            returns only scores for these sequences and in the name order

        Returns
        -------
        numpy array of log2 based scores at every position
        """
        assert not self.is_ragged(), "all sequences must have same length"
        assert pssm or path, "Must specify a PSSM or a path"
        assert not (pssm and path), "Can only specify one of pssm, path"

        if isinstance(names, str):
            names = [names]

        if path:
            pssm = load_pssm(path, background=background, pseudocount=pseudocount)

        assert set(pssm.motifs) == set(self.moltype)

        seqs = [self.seqs[n] for n in names] if names else self.seqs
        result = [pssm.score_seq(seq) for seq in ui.series(seqs)]

        return numpy.array(result)

    def get_ambiguous_positions(self):
        """Returns dict of seq:{position:char} for ambiguous chars.

        Used in likelihood calculations.
        """
        result = {}
        alpha = self.moltype.most_degen_alphabet()
        for name in self.names:
            result[name] = ambig = {}
            array = numpy.array(self.seqs[name])
            for i in numpy.where(array > alpha.gap_index)[0]:
                ambig[i] = alpha[array[i]]
        return result

    def trim_stop_codons(
        self, gc: c3_genetic_code.GeneticCodeChoiceType = 1, strict: bool = False
    ) -> typing_extensions.Self:
        """Removes any terminal stop codons from the sequences

        Parameters
        ----------
        gc
            valid input to cogent3.get_code(), a genetic code object, number
            or name, defaults to standard code
        strict
            If True, raises an exception if a seq length not divisible by 3
        """
        if not self.has_terminal_stop(gc=gc, strict=strict):
            return self

        new_seqs = {
            self.name_map[s.name]: s.trim_stop_codon(gc=gc, strict=strict).to_array(
                apply_transforms=False,
            )
            for s in self.seqs
        }

        init_kwargs = self._get_init_kwargs()
        init_kwargs["seqs_data"] = self._seqs_data.from_seqs(
            data=new_seqs,
            alphabet=self._seqs_data.alphabet,
            offset=self._seqs_data.offset,
            check=False,
        )
        init_kwargs["annotation_db"] = self._annotation_db
        return self.__class__(**init_kwargs)

    def counts_per_seq(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = False,
        warn: bool = False,
    ) -> MotifCountsArray:  # refactor: using array
        """counts of motifs per sequence

        Parameters
        ----------
        motif_length
            number of characters per tuple.
        include_ambiguity
            if True, motifs containing ambiguous characters from the seq moltype
            are included. No expansion of those is attempted.
        allow_gap
            if True, motifs containing a gap character are included.
        warn
            warns if motif_length > 1 and collection trimmed to produce motif
            columns.

        Notes
        -----

        only non-overlapping motifs are counted
        """
        counts = []
        motifs = set()
        for name in self.names:
            seq = self.get_seq(name)
            c = seq.counts(
                motif_length=motif_length,
                include_ambiguity=include_ambiguity,
                allow_gap=allow_gap,
                warn=warn,
            )
            motifs.update(c.keys())
            counts.append(c)
        motifs = sorted(motifs)
        for i, c in enumerate(counts):
            counts[i] = c.tolist(motifs)
        return MotifCountsArray(counts, motifs, row_indices=self.names)

    def counts(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = False,
    ) -> MotifCountsArray:
        """counts of motifs

        Parameters
        ----------
        motif_length
            number of elements per character.
        include_ambiguity
            if True, motifs containing ambiguous characters from the seq moltype
            are included. No expansion of those is attempted.
        allow_gap
            if True, motifs containing a gap character are included.
        exclude_unobserved
            if True, unobserved motif combinations are excluded.

        Notes
        -----

        only non-overlapping motifs are counted
        """
        per_seq = self.counts_per_seq(
            motif_length=motif_length,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
            exclude_unobserved=exclude_unobserved,
        )
        return per_seq.motif_totals()

    def count_ambiguous_per_seq(self) -> DictArray:
        """Counts of ambiguous characters per sequence."""

        darr = DictArrayTemplate(self.names)
        counts = numpy.array(
            [self.seqs[name].count_ambiguous() for name in self.names],
            dtype=numpy.uint32,
        )

        return darr.wrap(counts)

    def get_motif_probs(
        self,
        alphabet: c3_alphabet.AlphabetABC = None,
        include_ambiguity: bool = False,
        exclude_unobserved: bool = False,
        allow_gap: bool = False,
        pseudocount: int = 0,
    ) -> dict:  # refactor: using array
        """Return a dictionary of motif probs, calculated as the averaged
        frequency across sequences.

        Parameters
        ----------
        alphabet
            alphabet to use for motifs
        include_ambiguity
            if True resolved ambiguous codes are included in estimation of
            frequencies.
        exclude_unobserved
            if True, motifs that are not present in the alignment are excluded
            from the returned dictionary.
        allow_gap
            allow gap motif
        pseudocount
            value to add to each count

        Notes
        -----

        only non-overlapping motifs are counted
        """
        moltype = self.moltype
        if alphabet is None:
            alphabet = moltype.alphabet
            if allow_gap:
                alphabet = moltype.gapped_alphabet

        motif_len = alphabet.motif_len
        counts = Counter()
        for seq_name in self.names:
            sequence = self.seqs[seq_name]
            if motif_len > 1:
                sequence = [
                    str(sequence[i : i + motif_len])
                    for i in range(0, len(sequence) + 1 - motif_len, motif_len)
                ]
            for motif in sequence:
                if not allow_gap and self.moltype.gap in motif:
                    continue

                counts[motif] += 1

        probs = {}
        if not exclude_unobserved:
            for motif in alphabet:
                probs[motif] = pseudocount

        for motif, count in list(counts.items()):
            motif_set = moltype.resolve_ambiguity(motif, alphabet=alphabet)
            if len(motif_set) > 1:
                if include_ambiguity:
                    count = float(count) / len(motif_set)
                else:
                    continue
            for motif in motif_set:
                probs[motif] = probs.get(motif, pseudocount) + count

        total = float(sum(probs.values()))
        for motif in probs:
            probs[motif] /= total

        return probs

    def probs_per_seq(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = False,
        warn: bool = False,
    ) -> MotifFreqsArray:
        """return frequency array of motifs per sequence

        Parameters
        ----------
        motif_length
            number of characters per motif
        include_ambiguity
            if True, include motifs containing ambiguous characters
        allow_gap
            if True, include motifs containing a gap character
        exclude_unobserved
            if True, exclude motifs not present in the sequences in
            the resulting array
        warn
            warns if motif_length > 1 and collection trimmed to produce motif
            columns.
        """

        counts = self.counts_per_seq(
            motif_length=motif_length,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
            exclude_unobserved=exclude_unobserved,
            warn=warn,
        )
        return None if counts is None else counts.to_freq_array()

    def entropy_per_seq(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = True,
        warn: bool = False,
    ) -> numpy.ndarray:
        """Returns the Shannon entropy per sequence.

        Parameters
        ----------
        motif_length: int
            number of characters per tuple.
        include_ambiguity: bool
            if True, motifs containing ambiguous characters
            from the seq moltype are included. No expansion of those is attempted.
        allow_gap: bool
            if True, motifs containing a gap character are included.
        exclude_unobserved: bool
            if True, unobserved motif combinations are excluded.
        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns

        Notes
        -----
        For motif_length > 1, it's advisable to specify exclude_unobserved=True,
        this avoids unnecessary calculations.
        """
        probs = self.probs_per_seq(
            motif_length=motif_length,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
            exclude_unobserved=exclude_unobserved,
            warn=warn,
        )

        return None if probs is None else probs.entropy()

    def get_lengths(
        self,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
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
            motif_length=1,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
        )
        return counts.row_sum()

    def pad_seqs(self, pad_length: OptInt = None):
        """Returns copy in which sequences are padded with the gap character to same length.

        Parameters
        ----------
        pad_length
            Length all sequences are to be padded to. Will pad to max sequence
            length if pad_length is None or less than max length.
        """

        max_len = max(
            self._seqs_data.get_seq_length(n) for n in self._name_map.values()
        )

        if pad_length is None or pad_length < max_len:
            pad_length = max_len

        padded_seqs = {}
        for seq in self.seqs:
            padded_seq = numpy.full(
                shape=pad_length,
                fill_value=self.moltype.gapped_alphabet.gap_index,
                dtype=self.moltype.most_degen_alphabet().dtype,
            )
            padded_seq[: len(seq)] = numpy.array(seq)
            # the padded_seqs dict will be used to create the seqs_data, so the
            # keys should be the seqids from the original seqs_data, if this differs
            # from the seq name, this will be recorded in the name_map
            padded_seqs[self.name_map[seq.name]] = padded_seq

        init_kwargs = self._get_init_kwargs()
        # when we access .seqs, if the collection has been reversed, this will have
        # been applied to the returned Sequence, so we reset it to False for the
        # returned collection
        init_kwargs["seqs_data"] = self._seqs_data.from_seqs(
            data=padded_seqs,
            alphabet=self._seqs_data.alphabet,
            offset=self._seqs_data.offset,
        )
        init_kwargs["is_reversed"] = False
        return self.__class__(**init_kwargs)

    def strand_symmetry(self, motif_length: int = 1):
        """returns dict of strand symmetry test results per seq"""
        return {s.name: s.strand_symmetry(motif_length=motif_length) for s in self.seqs}

    def is_ragged(self) -> bool:
        """rerturns True if sequences are of different lengths"""
        return (
            len({self._seqs_data.get_seq_length(n) for n in self._name_map.values()})
            > 1
        )

    def has_terminal_stop(
        self, gc: c3_genetic_code.GeneticCodeChoiceType = 1, strict: bool = False
    ) -> bool:
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

    def get_identical_sets(
        self,
        mask_degen: bool = False,
    ) -> list[set]:  # refactor: array/simplify
        """returns sets of names for sequences that are identical

        Parameters
        ----------
        mask_degen
            if True, degenerate characters are ignored

        """

        if self.is_ragged():
            msg = "not all seqs same length, cannot get identical sets"
            raise ValueError(msg)

        if mask_degen and not self.moltype.degen_alphabet:
            warnings.warn(
                "in get_identical_sets, mask_degen has no effect as moltype "
                f"{self.moltype.label!r} has no degenerate characters",
                UserWarning,
                stacklevel=2,
            )
            mask_degen = False

        def reduced(seq, indices):
            return "".join(seq[i] for i in range(len(seq)) if i not in indices)

        identical_sets = []
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
            return list(dupes.values())

        mask_posns = {
            name: self.moltype.get_degenerate_positions(seq, include_gap=True)
            for name, seq in seqs.items()
        }

        for i in range(len(self.names) - 1):
            n1 = self.names[i]
            if n1 in seen:
                continue

            seq1 = seqs[n1]
            group = set()
            for j in range(i + 1, len(self.names)):
                n2 = self.names[j]
                if n2 in seen:
                    continue

                seq2 = seqs[n2]
                pos = mask_posns[n1] + mask_posns[n2]

                if pos:
                    seq1 = reduced(seq1, pos)
                    seq2 = reduced(seq2, pos)

                if seq1 == seq2:
                    seen.append(n2)
                    group.update([n1, n2])

            if group:
                identical_sets.append(group)

        return identical_sets

    def get_similar(
        self,
        target: c3_sequence.Sequence,
        min_similarity: float = 0.0,
        max_similarity: float = 1.0,
        metric: Callable[
            [c3_sequence.Sequence, c3_sequence.Sequence],
            float,
        ] = c3_sequence.frac_same,
        transform: bool | None = None,
    ) -> SequenceCollection:
        """Returns new SequenceCollection containing sequences similar to target.

        Parameters
        ----------
        target
            sequence object to compare to. Can be in the collection.
        min_similarity
            minimum similarity that will be kept. Default 0.0.
        max_similarity
            maximum similarity that will be kept. Default 1.0.
        metric
            a similarity function to use. Must be f(first_seq, second_seq).
            The default metric is fraction similarity, ranging from 0.0 (0%
            identical) to 1.0 (100% identical). The Sequence class have lots
            of methods that can be passed in as unbound methods to act as the
            metric, e.g. frac_same_gaps.
        transform
            transformation function to use on the sequences before the metric
            is calculated. If None, uses the whole sequences in each case. A
            frequent transformation is a function that returns a specified range
            of a sequence, e.g. eliminating the ends. Note that the transform
            applies to both the real sequence and the target sequence.

        Notes
        -----
        both min_similarity and max_similarity are inclusive.

        Warnings
        --------
        if the transformation changes the type of the sequence (e.g. extracting
        a string from an RnaSequence object), distance metrics that depend on
        instance data of the original class may fail.
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

    def __str__(self) -> str:
        """Returns self in FASTA-format, respecting name order."""
        from cogent3.format.sequence import FORMATTERS

        return FORMATTERS["fasta"](self.to_dict())

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            return False
        self_init = self._get_init_kwargs()
        other_init = other._get_init_kwargs()
        for key, self_val in self_init.items():
            if key in ("annotation_db", "slice_record"):
                continue
            other_val = other_init.get(key)
            if self_val != other_val:
                return False
        return True

    def __ne__(self, other: object) -> bool:
        return not self == other

    def __repr__(self) -> str:
        seqs = []
        limit = 10
        delimiter = ""

        repr_seq_names = [min(self.names, key=lambda name: len(self.seqs[name]))]
        if len(self.names) > 1:
            # In case of a tie, min and max return first.
            # reversed ensures if all seqs are of same length, different seqs are returned
            repr_seq_names.append(
                max(reversed(self.names), key=lambda name: len(self.seqs[name])),
            )

        for name in repr_seq_names:
            elts = list(str(self.seqs[name])[: limit + 1])
            if len(elts) > limit:
                elts[-1] = "..."
            seqs.append(f"{name}[{delimiter.join(elts)}]")

        if len(self.names) > 2:
            seqs.insert(1, "...")

        seqs = ", ".join(seqs)

        return f"{len(self.names)}x {self.moltype.label} seqcollection: ({seqs})"

    def _repr_html_(self) -> str:
        settings = self._repr_policy.copy()
        env_vals = get_setting_from_environ(
            "COGENT3_ALIGNMENT_REPR_POLICY",
            {"num_seqs": int, "num_pos": int, "wrap": int},
        )
        settings.update(env_vals)
        return self.to_html(
            name_order=self.names[: settings["num_seqs"]],
            limit=settings["num_pos"],
            wrap=settings["wrap"],
        )

    def to_html(
        self,
        name_order: typing.Sequence[str] | None = None,
        wrap: int = 60,
        limit: OptInt = None,
        colors: Mapping[str, str] | None = None,
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
        --------

        In a jupyter notebook, this code is used to provide the representation.

        .. code-block:: python

            seq_col  # is rendered by jupyter

        You can directly use the result for display in a notebook as

        .. code-block:: python

            from IPython.core.display import HTML

            HTML(seq_col.to_html())
        """
        css, styles = self.moltype.get_css_style(
            colors=colors,
            font_size=font_size,
            font_family=font_family,
        )

        seq_lengths = numpy.array(
            [self._seqs_data.get_seq_length(n) for n in self._name_map.values()],
        )
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
                    [""] * (max_truncated_len - len(styled_seqs[name])),
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
            for n, s in zip(name_order, seqblock, strict=False):
                s = "".join(s)
                # Filter out rows that are empty (due to combination of shorter sequences + wrapping)
                if s != "":
                    row = "".join([label_ % n, seq_ % s])
                    table.append(f"<tr>{row}</tr>")
        table.append("</table>")
        if (limit and limit < len(selected.names)) or (
            name_order and len(name_order) < len(selected.names)
        ):
            summary = f"{self.num_seqs} x {{min={min_val}, median={med_val}, max={max_val}}} (truncated to {len(name_order) if name_order else len(selected.names)} x {limit or len(selected)}) {selected.moltype.label} sequence collection"
        else:
            summary = f"{self.num_seqs} x {{min={min_val}, median={med_val}, max={max_val}}} {selected.moltype.label} sequence collection"

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
        num_seqs: OptInt = None,
        num_pos: OptInt = None,
        ref_name: OptInt = None,
        wrap: OptInt = None,
    ) -> None:
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
                msg = "num_seqs is not an integer"
                raise TypeError(msg)
            self._repr_policy["num_seqs"] = num_seqs

        if num_pos:
            if not isinstance(num_pos, int):
                msg = "num_pos is not an integer"
                raise TypeError(msg)
            self._repr_policy["num_pos"] = num_pos

        if ref_name:
            if not isinstance(ref_name, str):
                msg = "ref_name is not a string"
                raise TypeError(msg)

            if ref_name != "longest" and ref_name not in self.names:
                msg = f"no sequence name matching {ref_name}"
                raise ValueError(msg)

            self._repr_policy["ref_name"] = ref_name

        if wrap:
            if not isinstance(wrap, int):
                msg = "wrap is not an integer"
                raise TypeError(msg)
            self._repr_policy["wrap"] = wrap

    def duplicated_seqs(self) -> list[list[str]]:
        """returns the names of duplicated sequences"""
        seq_hashes = collections.defaultdict(list)
        for s in self.seqs:
            seq_hashes[array_hash64(numpy.array(s))].append(s.name)
        return [v for v in seq_hashes.values() if len(v) > 1]

    def drop_duplicated_seqs(self) -> typing_extensions.Self:
        """returns self without duplicated sequences

        Notes
        -----
        Retains the first sequence of each duplicte group.
        """
        dupes = self.duplicated_seqs()
        if not dupes:
            return self

        omit = []
        for group in dupes:
            omit.extend(group[1:])
        return self.take_seqs(omit, negate=True)


@register_deserialiser(
    get_object_provenance(SequenceCollection),
    "cogent3.core.alignment.SequenceCollection",
)
def deserialise_sequence_collection(
    data: dict[str, str | dict[str, str]],
) -> SequenceCollection:
    """deserialise SequenceCollection"""
    if "init_args" not in data:
        return deserialise_old_to_new_type_seqcoll(data)

    return SequenceCollection.from_rich_dict(data)


def deserialise_old_to_new_type_seqcoll(
    data: dict[str, str | dict[str, str]],
) -> SequenceCollection:
    """deserialise old type SequenceCollection as a new type collection"""
    moltype_name = data["moltype"]
    moltype_name = "text" if moltype_name == "bytes" else moltype_name
    info_data = data["info"]
    source = info_data.pop("source", None)
    seq_data = {
        seqid: record["seq"]["init_args"]["seq"]
        for seqid, record in data["seqs"].items()
    }
    return make_unaligned_seqs(
        seq_data, moltype=moltype_name, info=info_data, source=source
    )


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
        if not isinstance(seq, c3_sequence.Sequence):
            continue

        db = seq.annotation_db

        if first is None and db:
            # TODO gah should this be a copy so immutable?
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


@dataclasses.dataclass
class raw_seq_data:
    seq: bytes | numpy.ndarray[int]
    name: OptStr = None
    parent_name: OptStr = None
    offset: int = 0
    is_reversed: bool = False


@singledispatch
def coerce_to_raw_seq_data(
    seq: StrORBytesORArrayOrSeq,
    moltype: c3_moltype.MolType,
    name: OptStr = None,
) -> raw_seq_data:
    """aggregates sequence data into a single object

    Parameters
    ----------
    seq
        sequence data, can be a string, bytes, numpy array or Sequence
        instance. The latter is converted to a numpy array.
    moltype
        name of a cogent3 molecular type, or a cogent3 MolType instance
    name
        name of the sequence

    Returns
    -------
        raw_seq_data
    """
    if isinstance(seq, Aligned):
        # convert the Aligned instance
        # into a Sequence instance that includes the gaps
        return coerce_to_raw_seq_data(seq.gapped_seq, moltype, name)
    msg = f"coerce_to_seq_data not implemented for {type(seq)}"
    raise TypeError(msg)


@coerce_to_raw_seq_data.register
def _(
    seq: c3_sequence.Sequence,
    moltype: c3_moltype.MolType,
    name: str,
) -> raw_seq_data:
    # converts the sequence to a numpy array
    seq = seq.to_moltype(moltype)
    parent_name, start, _, step = seq.parent_coordinates()
    raw_seq = numpy.array(seq)
    return raw_seq_data(
        seq=raw_seq,
        name=name or seq.name,
        parent_name=parent_name,
        offset=start,
        is_reversed=step < 0,
    )


@coerce_to_raw_seq_data.register
def _(
    seq: str,
    moltype: c3_moltype.MolType,
    name: OptStr = None,
) -> raw_seq_data:
    return coerce_to_raw_seq_data(seq.encode("utf8"), moltype, name)


@coerce_to_raw_seq_data.register
def _(
    seq: numpy.ndarray,
    moltype: c3_moltype.MolType,  # noqa: ARG001
    name: OptStr = None,
) -> raw_seq_data:
    return raw_seq_data(seq=seq, name=name)


@coerce_to_raw_seq_data.register
def _(
    seq: bytes,
    moltype: c3_moltype.MolType,
    name: OptStr = None,
) -> raw_seq_data:
    # converts the sequence to a upper case bytes, and applies
    # moltype coercion if needed (e.g. RNA to DNA replaces U with T)
    seq = seq.upper()
    seq = moltype.coerce_to(seq) if moltype.coerce_to else seq
    return raw_seq_data(seq=seq, name=name)


CT = tuple[dict[str, StrORBytesORArray], dict[str, int], set[str]]


def prep_for_seqs_data(
    data: dict[str, StrORBytesORArrayOrSeq],
    moltype: c3_moltype.MolType,
    seq_namer: _SeqNamer,
) -> CT:
    """normalises input data for constructing a SeqsData object

    Parameters
    ----------
    data
        a dict[str, StrORBytesORArray] where the key is the sequence name and
        the value the sequence, or a series of Sequences instances
    moltype
        name of a cogent3 molecular type, or a cogent3 MolType instance
    seq_namer
        callback that takes the sequence name and transforms it to a new name

    Returns
    -------
    seq data as dict[str, bytes | numpy.ndarray], offsets as dict[str, int],
    reversed sequences as set[str], name_map as dict[str, str]
    """
    seqs = {}  # for the (Aligned)SeqsDataABC
    offsets = {}  # for the (Aligned)SeqsDataABC
    rvd = set()
    for name, seq in data.items():
        name = seq_namer(seq=seq, name=name)  # noqa: PLW2901
        seq_data = coerce_to_raw_seq_data(seq, moltype, name=name)
        if seq_data.offset:
            offsets[seq_data.parent_name or name] = seq_data.offset
        seqs[seq_data.parent_name or seq_data.name] = seq_data.seq
        if seq_data.is_reversed:
            rvd.add(name)

    return seqs, offsets, rvd


@singledispatch
def _make_name_seq_mapping(
    data: typing.Sequence[StrORBytesORArrayOrSeq],
    seq_namer: _SeqNamer,
) -> dict[str, StrORBytesORArrayOrSeq]:
    """returns a dict mapping names to sequences

    Parameters
    ----------
    data
        a dict of {name: seq, ...}, or python sequence of StrORBytesORArrayOrSeq
    seq_namer
        callback that takes the sequence and optionally a name and
        returns a new name
    """
    if not hasattr(data, "seqs"):
        msg = f"_make_name_seq_mapping not implemented for {type(data)}"
        raise NotImplementedError(msg)

    return {seq_namer(seq=record): record for record in data.seqs}


@_make_name_seq_mapping.register
def _(
    data: dict,
    seq_namer: _SeqNamer,
) -> dict[str, StrORBytesORArray]:
    return {seq_namer(seq=seq, name=name): seq for name, seq in data.items()}


@_make_name_seq_mapping.register
def _(
    data: list,
    seq_namer: _SeqNamer,
) -> dict[str, StrORBytesORArray]:
    if isinstance(data[0], (list, tuple)):
        # handle case where we've been given data like
        # for example [[name, seq], ...]
        with contextlib.suppress(ValueError):
            return _make_name_seq_mapping(dict(data), seq_namer)
    return {seq_namer(seq=record): record for record in data}


@_make_name_seq_mapping.register
def _(
    data: tuple,
    seq_namer: _SeqNamer,
) -> dict[str, StrORBytesORArray]:
    return _make_name_seq_mapping(list(data), seq_namer)


@_make_name_seq_mapping.register
def _(
    data: set,
    seq_namer: _SeqNamer,
) -> dict[str, StrORBytesORArray]:
    return _make_name_seq_mapping(list(data), seq_namer)


def _seqname_parent_name(
    record: StrORBytesORArray,
    name: str | None = None,
) -> tuple[str, str]:
    if hasattr(record, "parent_coordinates"):
        parent_name, *_ = record.parent_coordinates()
        return name or record.name, parent_name or name
    if hasattr(record, "name"):
        return name or record.name, name or record.name

    return name, name


def make_name_map(data: dict[str, StrORBytesORArray]) -> dict[str, str]:
    """returns a dict mapping names to parent names

    Parameters
    ----------
    data
        a dict of {name: seq, ...}

    Returns
    -------
        empty dict if names and parent names are always equal
    """
    name_map = {}
    for name, record in data.items():
        new_name, parent_name = _seqname_parent_name(record, name=name)
        if new_name == parent_name:
            continue
        name_map[new_name] = parent_name

    return name_map


def make_unaligned_storage(
    data: dict[str, bytes | numpy.ndarray[int]],
    *,
    moltype: str | c3_moltype.MolType,
    label_to_name: OptRenamerCallable = None,
    offset: DictStrInt | None = None,
    reversed_seqs: set[str] | None = None,
    storage_backend: str | None = None,
    **kwargs,
) -> SeqsData:
    """makes the unaligned storage instance for a SequenceCollection

    Parameters
    ----------
    data
        {name: seq}
    moltype
        label or instance of a cogent3 MolType
    label_to_name
        callable to convert the original name to a new name
    offset
        {name: offset} where the offset is the start position of the
        sequence in the parent sequence
    reversed_seqs
        set of names that are on the reverse strand of the parent sequence
    storage_backend
        name of a third-party storage driver to provide storage functionality
    kwargs
        additional keyword arguments for the storage driver

    Notes
    -----
    This function is intended for use primarly by make_unaligned_seqs function.
    """
    moltype = c3_moltype.get_moltype(moltype)
    alphabet = moltype.most_degen_alphabet()
    # if we have Sequences, we need to construct the name map before we construct
    # the SeqsData object - however, if a name_map is provided, we assume that it
    # corrects for any naming differences in data and skip this step
    assign_names = _SeqNamer(name_func=label_to_name)
    seqs_data, offs, rvd = prep_for_seqs_data(data, moltype, assign_names)
    offset = offset or {}
    offset = {**offs, **offset}
    # seqs_data keys should be the same as the value of name_map
    # name_map keys correspond to names in the sequence collection
    # name_map values correspond to names in seqs_data
    sd_kwargs = {
        "data": seqs_data,
        "alphabet": alphabet,
        "offset": offset,
        "reversed_seqs": reversed_seqs or rvd,
        **kwargs,
    }
    klass = cogent3._plugin.get_unaligned_storage_driver(storage_backend)  # noqa: SLF001
    return klass.from_seqs(**sd_kwargs)


@c3warn.deprecated_args("2025.9", "no longer has an effect", discontinued="new_type")
@singledispatch
def make_unaligned_seqs(
    data: dict[str, StrORBytesORArray] | list | SeqsDataABC,
    *,
    moltype: str | c3_moltype.MolType,
    label_to_name: OptRenamerCallable = None,
    info: OptDict = None,
    source: OptPathType = None,
    annotation_db: SupportsFeatures = None,
    offset: DictStrInt | None = None,
    name_map: DictStrStr | None = None,
    is_reversed: bool = False,
    reversed_seqs: set[str] | None = None,
    storage_backend: str | None = None,
    **kwargs,
) -> SequenceCollection:
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
    annotation_db
        annotation database to attach to the collection
    offset
        a dict mapping names to annotation offsets
    name_map
        a dict mapping sequence names to "parent" sequence names. The parent
        name will be used for querying a annotation_db.
    is_reversed
        entire collection has been reverse complemented
    reversed_seqs
        set of names that are on the reverse strand of the parent sequence
    storage_backend
        name of the storage backend to use for the SeqsData object, defaults to
        cogent3 builtin.
    kwargs
        keyword arguments for the storage driver

    Notes
    -----
    If no annotation_db is provided, but the sequences are annotated, an
    annotation_db is created by merging any annotation db's found in the sequences.
    If the sequences are annotated AND an annotation_db is provided, only the
    annotation_db is used.
    """
    # refactor: design
    # rename offset to offsets as it could track potentially multiple offsets

    if len(data) == 0:
        msg = "data must be at least one sequence."
        raise ValueError(msg)

    annotation_db = annotation_db or merged_db_collection(data)
    assign_names = _SeqNamer(name_func=label_to_name)
    data = _make_name_seq_mapping(data, assign_names)
    if name_map is None:
        name_map = make_name_map(data) or None

    seqs_data = make_unaligned_storage(
        data,
        label_to_name=label_to_name,
        moltype=moltype,
        offset=offset,
        reversed_seqs=reversed_seqs,
        storage_backend=storage_backend,
        **kwargs,
    )
    # as they were handled in this function, we do not pass on:
    # offset
    # label_to_name
    # reversed_seqs
    # storage_backend

    return make_unaligned_seqs(
        seqs_data,
        moltype=moltype,
        info=info,
        source=source,
        annotation_db=annotation_db,
        name_map=name_map,
        is_reversed=is_reversed,
    )


@make_unaligned_seqs.register
def _(
    data: SeqsDataABC,
    *,
    moltype: str | c3_moltype.MolType,
    label_to_name: OptRenamerCallable = None,
    info: dict | None = None,
    source: OptPathType = None,
    annotation_db: SupportsFeatures = None,
    offset: dict[str, int] | None = None,
    name_map: dict[str, str] | None = None,
    is_reversed: bool = False,
) -> SequenceCollection:
    moltype = c3_moltype.get_moltype(moltype)
    if not moltype.is_compatible_alphabet(data.alphabet):
        msg = f"Provided moltype: {moltype} is not compatible with SeqsData alphabet {data.alphabet}"
        raise ValueError(
            msg,
        )

    # we cannot set offset when creating from an SeqsData
    if offset:
        msg = f"Setting offset is not supported for {data=}"
        raise ValueError(msg)

    info = info if isinstance(info, dict) else {}
    source = str(source) if source else str(info.pop("source", "unknown"))
    seqs = SequenceCollection(
        seqs_data=data,
        moltype=moltype,
        info=info,
        annotation_db=annotation_db,
        source=source,
        name_map=name_map,
        is_reversed=is_reversed,
    )
    if label_to_name:
        seqs = seqs.rename_seqs(label_to_name)
    return seqs


@singledispatch
def decompose_gapped_seq(
    seq: typing.union[StrORBytesORArray, c3_sequence.Sequence],
    *,
    alphabet: c3_alphabet.AlphabetABC,
    missing_as_gap: bool = True,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Takes a sequence with (or without) gaps and returns an ungapped sequence
    and a map of the position and length of gaps in the original parent sequence
    """
    msg = f"decompose_gapped_seq not implemented for type {type(seq)}"
    raise NotImplementedError(
        msg,
    )


@decompose_gapped_seq.register
def _(
    seq: numpy.ndarray,
    *,
    alphabet: c3_alphabet.AlphabetABC,
    missing_as_gap: bool = True,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    if missing_as_gap and alphabet.missing_index:
        missing_index = numpy.uint8(alphabet.missing_index)
    else:
        missing_index = -1
    return decompose_gapped_seq_array(
        seq.astype(alphabet.dtype),
        alphabet.gap_index,
        missing_index=missing_index,
    )


@decompose_gapped_seq.register
def _(
    seq: str,
    *,
    alphabet: c3_alphabet.AlphabetABC,
    missing_as_gap: bool = True,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    if not alphabet.is_valid(seq):
        msg = f"Sequence is invalid for alphabet {alphabet}"
        raise c3_alphabet.AlphabetError(msg)

    return decompose_gapped_seq(
        alphabet.to_indices(seq),
        alphabet=alphabet,
        missing_as_gap=missing_as_gap,
    )


@decompose_gapped_seq.register
def _(
    seq: bytes,
    *,
    alphabet: c3_alphabet.AlphabetABC,
    missing_as_gap: bool = True,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    return decompose_gapped_seq(
        seq.decode("utf-8"),
        alphabet=alphabet,
        missing_as_gap=missing_as_gap,
    )


@decompose_gapped_seq.register
def _(
    seq: c3_sequence.Sequence,
    *,
    alphabet: c3_alphabet.AlphabetABC,
    missing_as_gap: bool = True,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    return decompose_gapped_seq(
        numpy.array(seq),
        alphabet=alphabet,
        missing_as_gap=missing_as_gap,
    )


@numba.jit(cache=True)
def decompose_gapped_seq_array(
    seq: numpy.ndarray,
    gap_index: int,
    missing_index: int = -1,
) -> tuple[numpy.ndarray, numpy.ndarray]:  # pragma: no cover
    """
    extracts the ungapped sequence and gap data from a gapped sequence

    Parameters
    ----------
    seq
        numpy array representing a gapped sequence
    gap_index
        from an alphabet
    missing_index
        from an alphabet, represents index for missing character

    Returns
    -------
    ungapped, [[gap_pos in sequence coords, cumulative gap length]]

    Notes
    -----
    being called by decompose_gapped_seq
    A missing_index is an ambiguity code that includes the gap character.
    Be careful in providing this value when dealing with sequences that
    may have had a feature masking applied.
    """
    seqlen = len(seq)
    working = numpy.empty((seqlen, numpy.int64(2)), dtype=numpy.int64)

    in_gap = False
    num_gaps = 0
    start = 0
    for i, base in enumerate(seq):
        gapped = base == gap_index or base == missing_index  # noqa
        if gapped and not in_gap:
            start = i
            in_gap = True
        elif not gapped and in_gap:
            working[num_gaps][:] = start, i - start
            num_gaps += 1
            in_gap = False

        if gapped and i == seqlen - 1:
            # end of sequence
            working[num_gaps][:] = start, i - start + 1
            num_gaps += 1

    if num_gaps == 0:
        return seq, numpy.empty((0, 2), dtype=numpy.int64)

    gap_coords = working[:num_gaps]
    gap_coords.T[1] = gap_coords.T[1].cumsum()
    # get gap start positions in sequence coords
    for index, cum_length in enumerate(numpy.append(0, gap_coords.T[1][:-1])):
        gap_coords[index][0] -= cum_length

    ungapped = numpy.empty(seqlen - gap_coords.T[1][-1], dtype=seq.dtype)
    seqpos = 0
    for base in seq:
        gapped = base == gap_index or base == missing_index  # noqa
        if not gapped:
            ungapped[seqpos] = base
            seqpos += 1

    return ungapped, gap_coords


@numba.jit(cache=True)
def compose_gapped_seq(
    ungapped_seq: numpy.ndarray,
    gaps: numpy.ndarray,
    gap_index: int,
) -> numpy.ndarray:  # pragma: no cover
    """reconstruct a gapped sequence from an ungapped sequence and gap data"""
    if not len(gaps):
        return ungapped_seq

    gapped_len = len(ungapped_seq) + gaps[-1, 1]

    gapped_seq = numpy.empty(gapped_len, dtype=ungapped_seq.dtype)

    pos = 0
    ungapped_pos = 0
    prev_gap_len = 0
    for gap_pos, cum_gap_len in gaps:
        gap_len = cum_gap_len - prev_gap_len
        prev_gap_len = cum_gap_len

        gapped_seq[pos : pos + gap_pos - ungapped_pos] = ungapped_seq[
            ungapped_pos:gap_pos
        ]
        pos += gap_pos - ungapped_pos
        ungapped_pos = gap_pos

        gapped_seq[pos : pos + gap_len] = gap_index
        pos += gap_len

    gapped_seq[pos:] = ungapped_seq[ungapped_pos:]

    return gapped_seq


class Aligned(AnnotatableMixin):
    """A single sequence in an alignment.

    Notes
    -----
    This is a wrapper around a ``AlignedDataView``. This class performs any
    complementing needed. It can be cast directly to a string or numpy array,
    e.g. ``numpy.array(<aligned instance>)`` returns a numpy unsigned 8-bit
    integer array.
    """

    __slots__ = ("_annotation_db", "_data", "_moltype", "_name")

    def __init__(
        self,
        data: AlignedDataView,
        moltype: c3_moltype.MolType,
        name: OptStr = None,
        annotation_db: SupportsFeatures | list[SupportsFeatures] | None = None,
    ) -> None:
        self._data = data
        self._moltype = moltype
        self._name = name or data.seqid
        self._annotation_db: list[SupportsFeatures] = self._init_annot_db_value(
            annotation_db
        )

    @classmethod
    def from_map_and_seq(
        cls, indel_map: IndelMap, seq: c3_sequence.Sequence
    ) -> Aligned:
        """Creates an Aligned instance from an indel map and a Sequence."""
        moltype = seq.moltype
        # refactor: design
        # this is a temporary approach during migration to new_types
        # to support the sequence alignment algorithms
        # a better solution is to create a AlignedDataView instance the
        # map and seq directly without requiring a parent AlignedSeqsData
        asd = AlignedSeqsData.from_seqs_and_gaps(
            seqs={seq.name: numpy.array(seq)},
            gaps={seq.name: indel_map.array},
            alphabet=moltype.most_degen_alphabet(),
        )

        return cls(asd.get_view(seq.name), moltype)

    @classmethod
    def from_map_and_aligned_data_view(
        cls,
        indel_map: IndelMap,
        seq: AlignedDataViewABC,
    ) -> Aligned:
        """Creates an Aligned instance from an indel map and AlignedDataView."""
        moltype = seq.alphabet.moltype
        seqid = seq.seqid
        seq = seq.array_value
        # refactor: design
        # see above comment in from_map_and_seq
        asd = AlignedSeqsData.from_seqs_and_gaps(
            seqs={seqid: seq},
            gaps={seqid: indel_map.array},
            alphabet=moltype.most_degen_alphabet(),
        )

        return cls(asd.get_view(seqid), moltype)

    def __len__(self) -> int:
        return len(self.map)

    @property
    def data(self) -> AlignedDataView:
        return self._data

    @property
    def map(self) -> IndelMap:
        return self.data.map

    @property
    def seq(self) -> c3_sequence.Sequence:
        """the ungapped sequence."""
        # if the slice record has abs(step) > 1, we cannot retain a connection
        # to the underlying aligned seq data container because the gaps are
        # not going to be modulo the step.
        if self.data.slice_record.plus_step == 1:
            # to complement or not handled by the view
            seq = self.data.get_seq_view()
        elif self.data.slice_record.step > 1:
            # we have a step, but no complementing will be required
            seq = self.moltype.degap(self.data.gapped_array_value)
        else:
            # gapped_array_value gives the reverse of the plus strand
            # so we need to complement it. We do that here because with a
            # step != 1, we cannot retain a connection to the underlying
            # annotations
            seq = self.moltype.degap(self.data.gapped_array_value)
            seq = self.moltype.complement(seq)

        mt_seq = self.moltype.make_seq(seq=seq, name=self.data.seqid)
        ann_db = self._annotation_db if self.data.slice_record.plus_step == 1 else None
        mt_seq.replace_annotation_db(ann_db)
        return mt_seq

    @property
    def gapped_seq(self) -> c3_sequence.Sequence:
        """Returns Sequence object, including gaps."""
        seq = self.data.gapped_array_value
        if self.data.slice_record.step < 0:
            seq = self.moltype.complement(seq)
        return self.moltype.make_seq(seq=seq, name=self.data.seqid)

    @property
    def moltype(self) -> c3_moltype.MolType:
        return self._moltype

    @property
    def name(self) -> str:
        return self._name

    def gap_vector(self) -> list[bool]:
        """Returns gap_vector of positions."""
        return self.gapped_seq.gap_vector()

    def make_feature(self, feature: FeatureDataType, alignment: Alignment) -> Feature:
        """returns a feature, not written into annotation_db"""
        annot = self.seq.make_feature(feature)
        inverted = self.map.to_feature_map().inverse()
        return annot.remapped_to(alignment, inverted)

    def __repr__(self) -> str:
        seq = f"{str(self)[:7]}... {len(self):,}" if len(self) > 10 else str(self)
        return (
            f"Aligned(name={self.name!r}, seq={seq!r}, moltype={self.moltype.name!r})"
        )

    def __str__(self) -> str:
        return str(self.gapped_seq)

    def __array__(
        self,
        dtype: numpy.dtype | None = None,
        copy: bool | None = None,
    ) -> numpy.ndarray[int]:
        return numpy.array(self.gapped_seq, dtype=dtype)

    def __bytes__(self) -> bytes:
        return bytes(self.gapped_seq)

    def __iter__(self):
        """Iterates over sequence one motif (e.g. char) at a time, incl. gaps"""
        yield from self.gapped_seq

    @singledispatchmethod
    def __getitem__(self, span: int | slice):
        msg = f"__getitem__ not implemented for {type(span)}"
        raise NotImplementedError(msg)

    @__getitem__.register
    def _(self, span: int):
        return self.__getitem__(slice(span, span + 1))

    @__getitem__.register
    def _(self, span: slice):
        return self.__class__(data=self.data[span], moltype=self.moltype)

    @__getitem__.register
    def _(self, span: FeatureMap):
        # we assume the feature map is in align coordinates
        data, gaps = self.slice_with_map(span)
        seqid = self.data.seqid
        seqs_data = self.data.parent.from_seqs_and_gaps(
            seqs={seqid: data},
            gaps={seqid: gaps},
            alphabet=self.moltype.most_degen_alphabet(),
        )
        view = seqs_data.get_view(seqid)

        return Aligned(view, self.moltype)

    def slice_with_map(self, span: FeatureMap) -> tuple[numpy.ndarray, numpy.ndarray]:
        start, end = span.start, span.end
        if span.useful and len(list(span.spans)) == 1:
            im = self.map[start:end]
            seq_start = self.map.get_seq_index(start)
            seq_end = self.map.get_seq_index(end)
            data = self.data.array_value[seq_start:seq_end]
            # .array_value will return the data in the correct orientation
            # which means we need to complement it if the data is reversed
            data = (
                self.moltype.complement(data)
                if self.data.slice_record.is_reversed
                else data
            )
        elif not span.useful:
            im = self.map[start:end]
            data = self.data.array_value[:0]
        else:
            # multiple spans
            align_coords = span.get_coordinates()
            im = self.map.joined_segments(align_coords)
            seq_map = self.map.make_seq_feature_map(span)
            # self.seq will return the data in the correct orientation
            # and will complement it if the data is reversed
            data = numpy.array(self.seq.gapped_by_map(seq_map))

        gaps = numpy.array([im.gap_pos, im.cum_gap_lengths]).T
        return data, gaps

    def parent_coordinates(self, seq_coords=False):
        """returns seqid, start, stop, strand on the parent sequence

        Parameters
        ----------
        seq_coords
            if True, the coordinates for the unaligned sequence
        """
        strand = -1 if self.data.is_reversed else 1
        seqid = self.data.seqid
        if not seq_coords:
            start = self.data.slice_record.parent_start
            stop = self.data.slice_record.parent_stop
        else:
            # AlignedDataView.parent_seq_coords uses it's indelmap, etc..
            # to return the necessary coordinates
            seqid, start, stop, strand = self.data.parent_seq_coords()

        return (
            seqid,
            start,
            stop,
            strand,
        )

    @extend_docstring_from(c3_sequence.Sequence.annotate_matches_to)
    def annotate_matches_to(
        self,
        pattern: str,
        biotype: str,
        name: str,
        allow_multiple: bool = False,
    ):
        return self.seq.annotate_matches_to(
            pattern=pattern,
            biotype=biotype,
            name=name,
            allow_multiple=allow_multiple,
        )


class AlignedSeqsDataABC(SeqsDataABC):
    """Abstract base class for respresenting the storage object for sequences underlying
    a Alignment.
    """

    # all methods that are from SeqsDataABC should work in sequence coordinates,
    # all methods unique to AlignedSeqsDataABC should work in aligned coordinates.
    # all indices provided to AlignedSeqsDataABC should be on the plus strand.
    __slots__ = ()

    @classmethod
    @abstractmethod
    def from_seqs_and_gaps(
        cls,
        *,
        seqs: dict[str, StrORBytesORArray],
        gaps: dict[str, numpy.ndarray],
        alphabet: c3_alphabet.AlphabetABC,
    ) -> typing_extensions.Self: ...

    @abstractmethod
    def __init__(
        self,
        *,
        gapped_seqs: numpy.ndarray,
        names: tuple[str],
        alphabet: c3_alphabet.AlphabetABC,
        ungapped_seqs: dict[str, numpy.ndarray] | None = None,
        gaps: dict[str, numpy.ndarray] | None = None,
        offset: DictStrInt | None = None,
        align_len: OptInt = None,
        check: bool = True,
        reversed_seqs: set[str] | None = None,
    ) -> None: ...

    @abstractmethod
    def get_view(
        self,
        seqid: str,
        slice_record: c3_sequence.SliceRecord | None = None,
    ) -> AlignedDataViewABC:
        # overriding the SeqsDataABC method as we support directly
        # providing the slice_record instance
        ...

    @classmethod
    @abstractmethod
    def from_names_and_array(
        cls,
        *,
        names: list[str],
        data: numpy.ndarray,
        alphabet: c3_alphabet.AlphabetABC,
    ) -> typing_extensions.Self: ...

    @property
    @abstractmethod
    def align_len(self) -> int: ...

    @abstractmethod
    def get_seq_array(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray: ...

    @abstractmethod
    def get_seq_str(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> str: ...

    @abstractmethod
    def get_seq_bytes(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> bytes: ...

    @abstractmethod
    def get_gapped_seq_array(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray: ...

    @abstractmethod
    def get_gapped_seq_str(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> str: ...

    @abstractmethod
    def get_gapped_seq_bytes(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> bytes: ...

    @abstractmethod
    def get_ungapped(
        self,
        name_map: dict[str, str],
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> tuple[dict, dict]:
        """
        Returns a dictionary of sequence data with no gaps or missing characters and
        a dictionary with information to construct a new SequenceCollection via
        make_unaligned_seqs.

        Parameters
        ----------
        name_map
            A dict of {aln_name: data_name, ...} indicating the mapping between
            names in the encompassing Alignment (aln_name) and the names in self
            (data_name).
        start
            The alignment starting position.
        stop
            The alignment stopping position.
        step
            The step size.

        Returns
        -------
        tuple
            A tuple containing the following:
            - seqs (dict): A dictionary of {name: seq, ...} where the sequences have no gaps
              or missing characters.
            - kwargs (dict): A dictionary of keyword arguments for make_unaligned_seqs, e.g.,
              {"offset": self.offset, "name_map": name_map}.
        """
        ...

    @abstractmethod
    def get_positions(
        self,
        names: PySeqStr,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray: ...

    @abstractmethod
    def copy(self, **kwargs) -> typing_extensions.Self: ...


def _gapped_seq_len(seq: numpy.ndarray, gap_map: numpy.ndarray) -> int:
    """calculate the gapped sequence length from a ungapped sequence and gap map

    Parameters
    ----------
    seq
        numpy array of sequence indices
    gap_map
        numpy array of [gap index, cumulative gap length] pairs
    """
    if isinstance(gap_map, IndelMap):
        gap_map = gap_map.array
    try:
        gap_len = gap_map[-1][1]
    except IndexError:  # no gaps
        return len(seq)

    return len(seq) + gap_len


class AlignedSeqsData(AlignedSeqsDataABC):
    """The builtin ``cogent3`` implementation of aligned sequences storage
    underlying an ``Alignment``. Indexing this object returns an ``AlignedDataView``
    which can realise the corresponding slice as a string, bytes, or numpy array,
    gapped or ungapped.

    Notes
    -----
    Methods on this object only accepts plust strand start, stop and step
    indices for selecting segments of data. It can return the gap coordinates
    for a sequence as used by ``IndelMap``.
    """

    __slots__ = (
        "_align_len",
        "_alphabet",
        "_gapped",
        "_gaps",
        "_name_to_index",
        "_names",
        "_offset",
        "_reversed",
        "_ungapped",
    )

    def __init__(
        self,
        *,
        gapped_seqs: numpy.ndarray,
        names: tuple[str],
        alphabet: c3_alphabet.CharAlphabet,
        ungapped_seqs: dict[str, numpy.ndarray] | None = None,
        gaps: dict[str, numpy.ndarray] | None = None,
        offset: DictStrInt | None = None,
        align_len: OptInt = None,
        check: bool = True,
        reversed_seqs: set[str] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        gapped_seqs
            2D numpy.uint8 array of aligned sequences. axis 0 are sequences,
            axis 1 are alignment positions
        names
            sequence names in order matching the axis 0 of gapped_seqs
        alphabet
            caharacter alphabet for the sequences
        ungapped_seqs
            a dictionary mapping names to 1D numpy.uint8 arrays of individual
            sequences without gaps. If not provided, computed on demand.
        gaps
            a dictionary mapping names to 1D numpy.int32 arrays of gap data,
            axis 0 is a gap axis 1 is [gap position in sequence coordinates,
            cumulative gap length].  If not provided, computed on demand.
        offset
            a dictionary of annotation offsets
        align_len
            length of the alignment, which must equal the gapped_seqs.shape[1]
        check
            validate any keys in offset, ungapped_seqs, gaps are a subset of names
        reversed_seqs
            names of seqs that are reverse complemented
        """
        self._alphabet = alphabet
        self._names = tuple(names)
        self._name_to_index = {name: i for i, name in enumerate(names)}
        self._gapped = gapped_seqs
        self._ungapped = ungapped_seqs or {}
        self._gaps = gaps or {}
        align_len = align_len or gapped_seqs.shape[1]
        self._reversed = frozenset(reversed_seqs or set())
        if align_len:
            assert align_len == gapped_seqs.shape[1], "mismatch in alignment length"

        self._align_len = align_len
        self._offset = offset or {}

        if check:
            if not set(names) >= set(self._gaps.keys()) or not set(names) >= set(
                self._ungapped.keys(),
            ):
                msg = "Keys in ungapped seqs and gaps must be subsets of names."
                raise ValueError(
                    msg,
                )
            if not set(names) >= set(self._offset):
                msg = "Keys in offset must be a subset of names."
                raise ValueError(msg)

            if len(names) != gapped_seqs.shape[0]:
                msg = f"{len(names)=} != {gapped_seqs.shape[0]=}"
                raise ValueError(msg)

    def __eq__(self, other: AlignedSeqsDataABC) -> bool:
        if not isinstance(other, self.__class__):
            return False
        attrs = (
            "_names",
            "_name_to_index",
            "_alphabet",
            "_align_len",
            "_offset",
        )
        for attr_name in attrs:
            self_attr = getattr(self, attr_name)
            other_attr = getattr(other, attr_name)
            if self_attr != other_attr:
                return False

        return numpy.all(self._gapped == other._gapped)

    def __ne__(self, other: object) -> bool:
        return not self == other

    @classmethod
    def from_seqs(
        cls,
        *,
        data: dict[str, StrORArray],
        alphabet: c3_alphabet.AlphabetABC,
        **kwargs,
    ) -> typing_extensions.Self:
        """Construct an AlignedSeqsData object from a dict of aligned sequences

        Parameters
        ----------
        data
            dict of gapped sequences {name: seq, ...}. sequences must all be
            the same length
        alphabet
            alphabet object for the sequences
        """
        seq_lengths = {len(v) for v in data.values()}
        if len(seq_lengths) != 1:
            msg = "All sequence lengths must be the same."
            raise ValueError(msg)

        align_len = seq_lengths.pop()
        names = tuple(data.keys())
        array_seqs = numpy.empty((len(names), align_len), dtype=alphabet.dtype)
        for i, name in enumerate(names):
            array_seqs[i] = alphabet.to_indices(data[name])
            if not alphabet.is_valid(data[name]):
                msg = f"Sequence {name} contains invalid characters."
                raise c3_alphabet.AlphabetError(
                    msg,
                )

        array_seqs.flags.writeable = False
        return cls(
            gapped_seqs=array_seqs,
            alphabet=alphabet,
            align_len=align_len,
            check=False,
            names=names,
            **kwargs,
        )

    @classmethod
    def from_seqs_and_gaps(
        cls,
        *,
        seqs: dict[str, StrORBytesORArray],
        gaps: dict[str, numpy.ndarray],
        alphabet: c3_alphabet.AlphabetABC,
        **kwargs,
    ) -> typing_extensions.Self:
        """Construct an AlignedSeqsData object from a dict of ungapped sequences
        and a corresponding dict of gap data.

        Parameters
        ----------
        seqs
            dict of ungapped sequences {name: seq, ...}
        gaps
            gap data {name: [[seq gap position, cumulative gap length], ...], ...}
        alphabet
            alphabet object for the sequences
        """
        names = tuple(kwargs.pop("names", seqs.keys()))
        if not names:
            msg = "seqs cannot be empty"
            raise ValueError(msg)

        align_len = kwargs.pop("align_len", None)
        if align_len is None:
            align_len = _gapped_seq_len(seqs[names[0]], gaps[names[0]])

        gapped_seqs = numpy.empty((len(names), align_len), dtype=alphabet.dtype)
        for i, name in enumerate(names):
            seq = alphabet.to_indices(seqs[name])
            seqs[name] = seq
            if name not in gaps:
                msg = f"Missing gap data for sequence {name!r}"
                raise ValueError(msg)
            gapped_seqs[i] = compose_gapped_seq(seq, gaps[name], alphabet.gap_index)
            assert len(gapped_seqs[i]) == align_len, "aligned lengths do not match"

        gapped_seqs.flags.writeable = False
        return cls(
            ungapped_seqs=seqs,
            gaps=gaps,
            gapped_seqs=gapped_seqs,
            alphabet=alphabet,
            names=names,
            align_len=align_len,
            **kwargs,
        )

    @classmethod
    def from_names_and_array(
        cls,
        *,
        names: PySeq[str],
        data: numpy.ndarray,
        alphabet: c3_alphabet.AlphabetABC,
    ) -> typing_extensions.Self:
        """Construct an AlignedSeqsData object from a list of names and a numpy
        array of aligned sequence data.

        Parameters
        ----------
        names
            list of sequence names
        data
            numpy array of aligned sequence data
        alphabet
            alphabet object for the sequences
        """
        if len(names) != data.shape[0] or not len(names):
            msg = "Number of names must match number of rows in data."
            raise ValueError(msg)

        gapped_seqs = data.astype(alphabet.dtype)
        gapped_seqs.flags.writeable = False
        return cls(
            gapped_seqs=gapped_seqs,
            names=names,
            alphabet=alphabet,
        )

    @property
    def names(self) -> tuple[str, ...]:
        """returns the names of the sequences in the storage"""
        return self._names

    @property
    def reversed_seqs(self) -> frozenset[str]:
        """names of sequences that are reverse complemented"""
        return self._reversed

    @property
    def alphabet(self) -> c3_alphabet.CharAlphabet:
        """the character alphabet for validating, encoding, decoding sequences"""
        return self._alphabet

    @property
    def align_len(self) -> int:
        """Return the length of the alignment."""
        return self._align_len

    @property
    def offset(self) -> dict[str, int]:
        """returns the offset of each sequence in the Alignment"""
        return {name: self._offset.get(name, 0) for name in self.names}

    def __len__(self) -> int:
        return self.align_len

    @singledispatchmethod
    def __getitem__(self, index: str | int) -> AlignedDataViewABC:
        return self.get_view(index)

    def get_seq_length(self, seqid: str) -> int:
        """return length of the unaligned seq for seqid"""
        return len(self._get_ungapped(seqid))

    @singledispatchmethod
    def get_view(
        self,
        seqid: str,
        slice_record: c3_sequence.SliceRecord | None = None,
    ) -> AlignedDataView:
        """reurns view of aligned sequence data for seqid

        Parameters
        ----------
        seqid
            sequence name
        slice_record
            slice record to use for slicing the data. If None, uses the
            default slice record for the entire sequence.
        """
        return AlignedDataView(
            parent=self,
            seqid=seqid,
            alphabet=self.alphabet,
            slice_record=slice_record,
        )

    @get_view.register
    def _(self, seqid: int):
        return self.get_view(self.names[seqid])

    def _get_gaps(self, seqid: str) -> numpy.ndarray:
        if seqid not in self._gaps:
            self._make_gaps_and_ungapped(seqid)
        return self._gaps[seqid]

    def _get_ungapped(self, seqid: str) -> numpy.ndarray:
        if seqid not in self._ungapped:
            self._make_gaps_and_ungapped(seqid)
        return self._ungapped[seqid]

    def get_gaps(self, seqid: str) -> numpy.ndarray:
        """returns the gap data for seqid"""
        return self._get_gaps(seqid)

    def _make_gaps_and_ungapped(self, seqid: str) -> None:
        if seqid in self._gaps and seqid in self._ungapped:
            # job already done
            return

        index = self._name_to_index[seqid]
        ungapped, gaps = decompose_gapped_seq(
            self._gapped[index],
            alphabet=self.alphabet,
        )
        self._gaps[seqid] = gaps
        self._ungapped[seqid] = ungapped

    def get_seq_array(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray:
        """Return ungapped sequence corresponding to seqid as an array of indices.

        Notes
        -----
        Assumes start/stop are in sequence coordinates. If seqid is in
        reversed_seqs, that sequence will be in plus strand orientation.
        It is client codes responsibility to ensure the coordinates are
        consistent with that.
        """
        start = start or 0
        stop = stop if stop is not None else self.get_seq_length(seqid)
        step = step or 1

        if start < 0 or stop < 0 or step < 1:
            msg = f"{start=}, {stop=}, {step=} not >= 1"
            raise ValueError(msg)

        out_len = (stop - start + step - 1) // step

        seq = numpy.empty(out_len, dtype=self.alphabet.dtype)
        seq[:] = self._get_ungapped(seqid)[start:stop:step]

        return seq

    def get_gapped_seq_array(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray:
        """Return sequence data corresponding to seqid as an array of indices.
        start/stop are in alignment coordinates. Includes gaps.
        """
        start = start or 0
        stop = stop if stop is not None else self.align_len
        step = step or 1
        if start < 0 or stop < 0 or step < 1:
            msg = f"{start=}, {stop=}, {step=} not >= 1"
            raise ValueError(msg)

        index = self._name_to_index[seqid]
        out_len = (stop - start + step - 1) // step
        gapped = numpy.empty(out_len, dtype=self.alphabet.dtype)
        gapped[:] = self._gapped[index][start:stop:step]

        return gapped

    def get_seq_str(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> str:
        """Return ungapped sequence corresponding to seqid as a string.
        start/stop are in sequence coordinates. Excludes gaps."""
        return self.alphabet.from_indices(
            self.get_seq_array(seqid=seqid, start=start, stop=stop, step=step),
        )

    def get_gapped_seq_str(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> str:
        """Return sequence corresponding to seqid as a string.
        start/stop are in alignment coordinates. Includes gaps."""
        return self.alphabet.from_indices(
            self.get_gapped_seq_array(seqid=seqid, start=start, stop=stop, step=step),
        )

    def get_seq_bytes(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> bytes:
        """Return ungapped sequence corresponding to seqid as a bytes string.
        start/stop are in sequence coordinates. Excludes gaps."""
        return self.get_seq_str(seqid=seqid, start=start, stop=stop, step=step).encode(
            "utf8",
        )

    def get_gapped_seq_bytes(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> bytes:
        """Return sequence corresponding to seqid as a bytes string.
        start/stop are in alignment coordinates. Includes gaps."""
        return self.get_gapped_seq_str(
            seqid=seqid,
            start=start,
            stop=stop,
            step=step,
        ).encode("utf8")

    @extend_docstring_from(AlignedSeqsDataABC.get_ungapped)
    def get_ungapped(
        self,
        name_map: dict[str, str],
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> tuple[dict, dict]:
        # redesign
        # if gaps exist, don't go via gapped seq
        # convert alignment coords into sequence coords using the location.align_to_seq_index function
        # this means we will need to convert coordinates to a plus strand slice
        if (start or 0) < 0 or (stop or 0) < 0 or (step or 1) <= 0:
            msg = f"{start=}, {stop=}, {step=} not >= 0"
            raise ValueError(msg)

        seq_array = numpy.empty(
            (len(name_map), self.align_len),
            dtype=self.alphabet.dtype,
        )
        names = tuple(name_map.values())
        for i, name in enumerate(names):
            index = self._name_to_index[name]
            seq_array[i] = self._gapped[index]
        seq_array = seq_array[:, start:stop:step]
        # now exclude gaps and missing
        seqs = {}
        for i, name in enumerate(names):
            seq = seq_array[i]
            indices = seq != self.alphabet.gap_index
            if self.alphabet.missing_index is not None:
                indices &= seq != self.alphabet.missing_index
            seqs[name] = seq[indices]

        offset = {n: v for n, v in self._offset.items() if n in names}
        return seqs, {
            "offset": offset,
            "name_map": name_map,
            "reversed_seqs": self._reversed,
        }

    def add_seqs(
        self,
        seqs: dict[str, StrORArray],
        force_unique_keys: bool = True,
        offset: dict[str, int] | None = None,
    ) -> AlignedSeqsData:
        """Returns a new AlignedSeqsData object with added sequences.

        Parameters
        ----------
        seqs
            dict of sequences to add {name: seq, ...}
        force_unique_keys
            if True, raises ValueError if any sequence names already exist in the collection
        offset
            dict of offsets relative to for the new sequences.
        """
        if force_unique_keys and any(name in self.names for name in seqs):
            msg = "One or more sequence names already exist in collection"
            raise ValueError(msg)

        new_seq_lens = {len(seq) for seq in seqs.values()}
        if len(new_seq_lens) != 1 or new_seq_lens.pop() != self.align_len:
            msg = "All sequences must be the same length as existing sequences"
            raise ValueError(
                msg,
            )

        new_seqs = dict(zip(self.names, self._gapped, strict=False))
        for name, seq in seqs.items():
            seq = self.alphabet.to_indices(seq)
            if not self.alphabet.is_valid(seq):
                msg = f"Sequence {name!r} contains invalid characters."
                raise c3_alphabet.AlphabetError(
                    msg,
                )
            seq.flags.writeable = False
            new_seqs[name] = seq

        names = tuple(new_seqs.keys())
        gapped = numpy.empty((len(names), self.align_len), dtype=self.alphabet.dtype)
        for i, name in enumerate(names):
            gapped[i] = new_seqs[name]

        return self.__class__(
            gapped_seqs=gapped,
            names=names,
            alphabet=self.alphabet,
            offset={**self._offset, **(offset or {})},
            align_len=self.align_len,
        )

    def to_alphabet(
        self,
        alphabet: c3_alphabet.AlphabetABC,
        check_valid: bool = True,
    ) -> typing_extensions.Self:
        """Returns a new AlignedSeqsData object with the same underlying data
        with a new alphabet."""
        if (
            len(alphabet) == len(self.alphabet)
            and len(
                {
                    (a, b)
                    for a, b in zip(self.alphabet, alphabet, strict=False)
                    if a != b
                },
            )
            == 1
        ):
            # special case where mapping between dna and rna
            return self.__class__(
                gapped_seqs=self._gapped,
                alphabet=alphabet,
                offset=self._offset,
                align_len=self.align_len,
                names=self.names,
            )

        gapped = numpy.empty(
            (len(self.names), self.align_len),
            dtype=self.alphabet.dtype,
        )

        for i in range(len(self.names)):
            seq_data = self._gapped[i]
            as_new_alpha = self.alphabet.convert_seq_array_to(
                seq=seq_data,
                alphabet=alphabet,
                check_valid=False,
            )
            if check_valid and not alphabet.is_valid(as_new_alpha):
                msg = (
                    f"Changing from old alphabet={self.alphabet} to new "
                    f"{alphabet=} is not valid for this data"
                )
                raise c3_alphabet.AlphabetError(msg)
            gapped[i] = as_new_alpha

        return self.__class__(
            gapped_seqs=gapped,
            alphabet=alphabet,
            offset=self._offset,
            names=self.names,
        )

    def get_positions(
        self,
        names: PySeqStr,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray:
        """returns an array of the selected positions for names."""
        indices = tuple(self._name_to_index[name] for name in names)
        start = start or 0
        stop = stop or self.align_len
        step = step or 1
        if start < 0 or stop < 0 or step < 1:
            msg = f"{start=}, {stop=}, {step=} not >= 1"
            raise ValueError(msg)

        if abs((start - stop) // step) == self.align_len:
            array_seqs = self._gapped[indices, :]
        else:
            array_seqs = self._gapped[indices, start:stop:step]

        return array_seqs

    def copy(self, **kwargs) -> typing_extensions.Self:
        """shallow copy of self

        Notes
        -----
        kwargs are passed to constructor and will over-ride existing values
        """
        init_args = {
            "gapped_seqs": self._gapped,
            "names": self._names,
            "alphabet": self._alphabet,
            "ungapped_seqs": self._ungapped,
            "gaps": self._gaps,
            "offset": self._offset,
            "align_len": self._align_len,
            "check": False,
            "reversed_seqs": self._reversed,
            **kwargs,
        }

        return self.__class__(**init_args)


class AlignedDataViewABC(c3_sequence.SeqViewABC):
    __slots__ = ()

    @abstractmethod
    def get_seq_view(self) -> c3_sequence.SeqViewABC: ...

    @property
    @abstractmethod
    def map(self) -> IndelMap: ...

    @property
    @abstractmethod
    def slice_record(self) -> c3_sequence.SliceRecordABC: ...

    @property
    @abstractmethod
    def gapped_str_value(self) -> str: ...

    @property
    @abstractmethod
    def gapped_array_value(self) -> numpy.ndarray: ...

    @property
    @abstractmethod
    def gapped_bytes_value(self) -> bytes: ...


class AlignedDataView(c3_sequence.SeqViewABC):
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
        alphabet: c3_alphabet.AlphabetABC,
        slice_record: OptSliceRecord = None,
    ) -> None:
        self.parent = parent
        self._seqid = seqid
        self.alphabet = alphabet
        self._parent_len = parent.align_len
        self._slice_record = (
            slice_record
            if slice_record is not None
            else c3_sequence.SliceRecord(parent_len=self._parent_len)
        )

    @property
    def slice_record(self) -> c3_sequence.SliceRecordABC:
        """the slice record for this view"""
        return self._slice_record

    @slice_record.setter
    def slice_record(self, value: c3_sequence.SliceRecordABC) -> None:
        self._slice_record = value

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
    def array_value(self) -> numpy.ndarray:
        """returns the numpy array of indices for the ungapped sequence"""
        value = self.parent.get_seq_array(
            seqid=self.seqid,
            start=self.map.get_seq_index(self.slice_record.plus_start),
            stop=self.map.get_seq_index(self.slice_record.plus_stop),
            step=self.map.get_seq_index(self.slice_record.plus_step),
        )
        return value[::-1] if self.slice_record.is_reversed else value

    @property
    def gapped_array_value(self) -> numpy.ndarray:
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
        dtype: numpy.dtype | None = None,
        copy: bool | None = None,
    ) -> numpy.ndarray[int]:
        arr = self.gapped_array_value
        if dtype:
            arr = arr.astype(dtype)
        return arr

    def __bytes__(self) -> bytes:
        return self.gapped_bytes_value

    def __getitem__(self, segment) -> typing_extensions.Self:
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

    def parent_seq_coords(self) -> tuple[str, int, int, int]:
        """returns seqid, start, stop, strand on the parent sequence"""
        # we want the coordinates on the parent sequence, which means we
        # need to use the parent's IndelMap for findings the correct indices.
        parent_map = self._parent_map()
        start = parent_map.get_seq_index(self.slice_record.parent_start)
        stop = parent_map.get_seq_index(self.slice_record.parent_stop)
        strand = -1 if self.slice_record.step < 0 else 1

        return self.seqid, start, stop, strand

    def copy(self, sliced: bool = False) -> typing_extensions.Self:
        """just returns self"""
        return self

    def _get_init_kwargs(self) -> dict:
        return {
            "parent": self.parent,
            "seqid": self.seqid,
            "alphabet": self.alphabet,
            "slice_record": self.slice_record,
        }

    def get_seq_view(self) -> c3_sequence.SeqViewABC:
        """returns view of ungapped sequence data for seqid"""
        # we want the parent coordinates in sequence coordinates
        # parent_seq_coords does not account for the stride
        seqid, start, stop, _ = self.parent_seq_coords()
        parent_len = self.parent.get_seq_length(seqid)
        offset = self.parent.offset.get(seqid)
        sr = c3_sequence.SliceRecord(
            start=start,
            stop=stop,
            parent_len=parent_len,
            offset=offset,
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


def make_gap_filter(template, gap_fraction, gap_run):
    """Returns f(seq) -> True if no gap runs and acceptable gap fraction.

    Calculations relative to template.
    gap_run = number of consecutive gaps allowed in either the template or seq
    gap_fraction = fraction of positions that either have a gap in the template
        but not in the seq or in the seq but not in the template
    NOTE: template and seq must both be ArraySequence objects.
    """
    template_gaps = numpy.array(template.gap_vector())

    def result(seq) -> bool:
        """Returns True if seq adhers to the gap threshold and gap fraction."""
        seq_gaps = numpy.array(seq.gap_vector())
        # check if gap amount bad
        if sum(seq_gaps != template_gaps) / float(len(seq)) > gap_fraction:
            return False
        # check if gap runs bad
        return not (
            b"\x01" * gap_run
            in numpy.logical_and(seq_gaps, numpy.logical_not(template_gaps))
            .astype(numpy.uint8)
            .tobytes()
            or b"\x01" * gap_run
            in numpy.logical_and(template_gaps, numpy.logical_not(seq_gaps))
            .astype(numpy.uint8)
            .tobytes()
        )

    return result


class _IndexableSeqs:
    """container that is created by SequenceCollection and Alignment instances"""

    def __init__(
        self,
        parent: SequenceCollection | Alignment,
        make_seq: typing.Callable[[str], c3_sequence.Sequence | Aligned],
    ) -> None:
        """
        Parameters
        ----------
        parent
            either a SequenceCollection or Alignment instance
        make_seq
            method on the parent that creates the correct object type when given a seqid
        """
        self.parent = parent
        self._make_seq = make_seq

    @singledispatchmethod
    def __getitem__(
        self,
        key: str | int | slice,
    ) -> c3_sequence.Sequence | Aligned:
        msg = f"indexing not supported for {type(key)}, try .take_seqs()"
        raise TypeError(msg)

    @__getitem__.register
    def _(self, key: int) -> c3_sequence.Sequence | Aligned:
        return self[self.parent.names[key]]

    @__getitem__.register
    def _(self, key: str) -> c3_sequence.Sequence | Aligned:
        return self._make_seq(key)

    def __repr__(self) -> str:
        one_seq = self[self.parent.names[0]]
        return f"({one_seq!r}, + {self.parent.num_seqs - 1} seqs)"

    def __len__(self) -> int:
        return self.parent.num_seqs

    def __iter__(self):
        for name in self.parent.names:
            yield self._make_seq(name)


class Alignment(SequenceCollection):
    """A collection of aligned sequences.

    Notes
    -----
    Should be constructed using ``make_aligned_seqs()``.
    """

    def __init__(
        self,
        seqs_data: AlignedSeqsDataABC,  # seqs_data
        slice_record: OptSliceRecord = None,
        **kwargs,
    ) -> None:
        super().__init__(seqs_data=seqs_data, **kwargs)
        self._slice_record = (
            slice_record
            if slice_record is not None
            else c3_sequence.SliceRecord(parent_len=self._seqs_data.align_len)
        )
        self._array_seqs = None

    def _post_init(self) -> None:
        self._seqs = _IndexableSeqs(self, make_seq=self._make_aligned)

    def __eq__(self, other: object) -> bool:
        return super().__eq__(other) and self._slice_record == other._slice_record

    def __ne__(self, other: object) -> bool:
        return not self == other

    @property
    def storage(self) -> AlignedSeqsDataABC:
        """the aligned sequence storage instance of the collection"""
        return self._seqs_data

    @storage.setter
    def storage(self, value: object) -> None:
        # storage cannot be set after initialisation
        msg = "storage cannot be set after initialisation"
        raise TypeError(msg)

    def _get_init_kwargs(self) -> dict:
        """returns the kwargs needed to re-instantiate the object"""
        return {
            "seqs_data": self._seqs_data,
            "moltype": self.moltype,
            "name_map": dict(self._name_map),
            "info": self.info.copy(),
            "annotation_db": self._annotation_db,
            "slice_record": self._slice_record,
            "source": self.source,
        }

    @singledispatchmethod
    def __getitem__(self, index) -> typing_extensions.Self:
        msg = f"__getitem__ not implemented for {type(index)}"
        raise NotImplementedError(msg)

    @__getitem__.register
    def _(self, index: str) -> typing_extensions.Self:
        return self.seqs[index]

    @__getitem__.register
    def _(self, index: int) -> typing_extensions.Self:
        new_slice = self._slice_record[index]
        kwargs = self._get_init_kwargs()
        kwargs["slice_record"] = new_slice
        return self.__class__(**kwargs)

    @__getitem__.register
    def _(self, index: slice) -> typing_extensions.Self:
        new_slice = self._slice_record[index]
        kwargs = self._get_init_kwargs()
        kwargs["slice_record"] = new_slice
        if new_slice.plus_step > 1:
            # we retain the annotation database only for "simple" slices
            kwargs.pop("annotation_db", None)

        return self.__class__(**kwargs)

    @__getitem__.register
    def _(self, index: FeatureMap) -> typing_extensions.Self:
        return self._mapped(index)

    @__getitem__.register
    def _(self, index: Feature):
        if index.parent is not self:
            msg = "This feature applied to the wrong sequence / alignment"
            raise ValueError(msg)
        return index.get_slice()

    def __repr__(self) -> str:
        seqs = []
        limit = 10
        delimiter = ""
        for count, name in enumerate(self.names):
            if count == 3:
                seqs.append("...")
                break
            elts = list(str(self.seqs[name])[: limit + 1])
            if len(elts) > limit:
                elts[-1] = "..."
            seqs.append(f"{name}[{delimiter.join(elts)}]")
        seqs = ", ".join(seqs)

        return f"{len(self.names)} x {len(self)} {self.moltype.label} alignment: {seqs}"

    def __len__(self) -> int:
        return len(self._slice_record)

    def __array__(
        self,
        dtype: numpy.dtype | None = None,
        copy: bool | None = None,
    ) -> numpy.ndarray[int]:
        return self.array_seqs

    def _make_aligned(self, seqid: str) -> Aligned:
        adv = self._seqs_data.get_view(
            self._name_map.get(seqid, seqid),
            slice_record=self._slice_record,
        )
        aligned = Aligned(data=adv, moltype=self.moltype, name=seqid)
        aligned.annotation_db = self._annotation_db
        return aligned

    @property
    def positions(self) -> list[list[str]]:
        # refactor: design
        # possibly rename to str_positions since we have array_positions
        from_indices = self.moltype.most_degen_alphabet().from_indices
        return [list(from_indices(pos)) for pos in self.array_positions]

    @property
    def array_seqs(self) -> numpy.ndarray:
        """Returns a numpy array of sequences, axis 0 is seqs in order
        corresponding to names"""
        if self._array_seqs is None:
            names = [self._name_map[n] for n in self.names]
            # create the dest array dim
            arr_seqs = self._seqs_data.get_positions(
                names=names,
                start=self._slice_record.plus_start,
                stop=self._slice_record.plus_stop,
                step=self._slice_record.plus_step,
            )
            if self.moltype.is_nucleic and self._slice_record.is_reversed:
                rev_complement = self.moltype.rc
                arr_seqs = arr_seqs.copy()
                arr_seqs.flags.writeable = True
                for i in range(arr_seqs.shape[0]):
                    arr_seqs[i] = rev_complement(arr_seqs[i])

            arr_seqs.flags.writeable = False  # make sure data is immutable
            self._array_seqs = arr_seqs

        return self._array_seqs

    @property
    def array_positions(self) -> numpy.ndarray:
        """Returns a numpy array of positions, axis 0 is alignment positions
        columns in order corresponding to names."""
        return self.array_seqs.T

    def get_seq(
        self,
        seqname: str,
        copy_annotations: bool = False,
    ) -> c3_sequence.Sequence:
        """Return a Sequence object for the specified seqname.

        Parameters
        ----------
        seqname
            name of the sequence to return
        copy_annotations
            if True, only the annotations for the specified sequence are copied
            to the annotation database of the Sequence object which is decoupled
            from this collection. If False, the connection to this collections db
            is retained.
        """
        seq = self.seqs[seqname].seq
        if copy_annotations and self._annotation_db:
            # we need to copy the sequence too to break the link to self.annotation_db
            seq = seq.copy(exclude_annotations=True)
            seq.annotation_db = type(self.annotation_db)()
            seq.annotation_db.update(annot_db=self.annotation_db, seqids=seqname)
            return seq

        seq.annotation_db = self._annotation_db
        return seq

    def get_gapped_seq(
        self,
        seqname: str,
        recode_gaps: bool = False,
    ) -> c3_sequence.Sequence:
        """Return a gapped Sequence object for the specified seqname.


        Parameters
        ----------
        seqname
            sequence name
        recode_gaps
            if True, gap characters are replaced by the most general
            ambiguity code, e.g. N for DNA and RNA

        Notes
        -----
        This method breaks the connection to the annotation database.
        """
        s = self.seqs[seqname].gapped_seq
        if recode_gaps:
            s = str(s)
            non_ambig = list(self.moltype)
            ambig = self.moltype.degenerate_from_seq(non_ambig)
            for gapchar in self.moltype.gaps:
                s = s.replace(gapchar, ambig)

        return self.moltype.make_seq(seq=s, name=seqname)

    def rc(self):
        """Returns the reverse complement of all sequences in the alignment.
        A synonym for reverse_complement.
        """
        init_kwargs = self._get_init_kwargs()
        init_kwargs["slice_record"] = self._slice_record[::-1]
        return self.__class__(**init_kwargs)

    def alignment_quality(self, app_name: str = "ic_score", **kwargs):
        """
        Computes the alignment quality using the indicated app

        Parameters
        ----------
        app_name
            name of an alignment score calculating app, e.g. 'ic_score',
            'cogent3_score', 'sp_score'

        kwargs
            keyword arguments to be passed to the app. Use
            ``cogent3.app_help(app_name)`` to see the available options.

        Returns
        -------
        float or a NotCompleted instance if the score could not be computed
        """
        app = cogent3.get_app(app_name, **kwargs)
        return app(self)

    def rename_seqs(self, renamer: Callable[[str], str]):
        """Returns new alignment with renamed sequences."""
        new = super().rename_seqs(renamer)

        if self._array_seqs is not None:
            new._array_seqs = self._array_seqs

        return new

    def iter_positions(
        self,
        pos_order: list | None = None,
    ) -> typing.Iterator[list, list, list]:
        """Iterates over positions in the alignment, in order.

        Parameters
        ----------
        pos_order
            list of indices specifying the column order. If None, the
            positions are iterated in order.

        Returns
        -------
        yields lists of elemenets for each position (column) in the alignment
        """
        # refactor: array
        # this could also iter columns of indices as a numpy array - could be an optional arg
        # refactor: add motif_length argument
        pos_order = pos_order or range(len(self))
        for pos in pos_order:
            yield [str(self[seq][pos]) for seq in self.names]

    def get_position_indices(
        self,
        f: Callable[[str], bool],
        negate: bool = False,
    ) -> list[int]:
        """Returns list of column indices for which f(col) is True.

        Parameters
        ----------
        f
          function that returns true/false given an alignment position
        negate
          if True, not f() is used
        """
        # refactor:
        # type hint for f
        # implement native
        new_f = negate_condition(f) if negate else f

        # refactor: design
        # use array_positions here
        return [i for i, col in enumerate(self.positions) if new_f(col)]

    def take_positions(
        self,
        cols: list[int] | numpy.ndarray[int],
        negate: bool = False,
    ) -> typing_extensions.Self:
        """Returns new Alignment containing only specified positions.

        Parameters
        ----------
        cols
            list of column indices to keep
        negate
            if True, all columns except those in cols are kept
        """
        # refactor: array - use array operations throughout method
        if negate:
            col_lookup = dict.fromkeys(cols)
            cols = [i for i in range(len(self)) if i not in col_lookup]

        new_data = {
            self.name_map[aligned.name]: numpy.array(aligned).take(cols)
            for aligned in self.seqs
        }
        seqs_data = self._seqs_data.from_seqs(
            data=new_data,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        kwargs = self._get_init_kwargs()
        kwargs["seqs_data"] = seqs_data
        kwargs.pop("annotation_db", None)
        kwargs.pop("slice_record", None)
        return self.__class__(**kwargs)

    def take_positions_if(
        self, f: Callable[[str], bool], negate: bool = False
    ) -> typing_extensions.Self:
        """Returns new Alignment containing cols where f(col) is True."""
        return self.take_positions(self.get_position_indices(f, negate=negate))

    def get_gap_array(self, include_ambiguity: bool = True) -> numpy.ndarray:
        """returns bool array with gap state True, False otherwise

        Parameters
        ----------
        include_ambiguity
            if True, ambiguity characters that include the gap state are
            included
        """
        alpha = self.moltype.most_degen_alphabet()
        gapped = self.array_seqs == alpha.gap_index
        if include_ambiguity:
            gapped = gapped | (self.array_seqs == alpha.missing_index)
        return gapped

    def iupac_consensus(self, allow_gap: bool = True) -> str:
        """Returns string containing IUPAC consensus sequence of the alignment."""
        exclude = set() if allow_gap else set(self.moltype.gaps)
        consensus = []
        degen = self.moltype.degenerate_from_seq
        for col in self.iter_positions():
            col = set(col) - exclude
            consensus.append(degen("".join(col)))
        return "".join(consensus)

    def majority_consensus(self) -> c3_sequence.Sequence:
        """Returns consensus sequence containing most frequent item at each
        position."""
        states = []
        data = zip(*map(str, self.seqs), strict=False)
        for pos in data:
            pos = CategoryCounter(pos)
            states.append(pos.mode)

        return self.moltype.make_seq(seq="".join(states))

    def counts_per_pos(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        warn: bool = False,
    ) -> DictArray:
        """return DictArray of counts per position

        Parameters
        ----------
        motif_length
            number of elements per character.
        include_ambiguity
            if True, motifs containing ambiguous characters from the seq moltype
            are included. No expansion of those is attempted.
        allow_gap
            if True, motifs containing a gap character are included.
        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns
        """
        # refactor: performance
        # use self.variable_positions and a numba decorated
        # function for counting k-mers the latter should allow returning the
        # first allowed state for when position is not variable

        align_len = len(self._slice_record)
        length = (align_len // motif_length) * motif_length
        if warn and align_len != length:
            warnings.warn(f"trimmed {align_len - length}", UserWarning, stacklevel=2)

        data = list(self.to_dict().values())
        alpha = self.moltype.alphabet.get_kmer_alphabet(motif_length)
        all_motifs = set()
        exclude_chars = set()
        if not allow_gap:
            exclude_chars.update(self.moltype.gap)

        if not include_ambiguity and self.moltype.degen_alphabet:
            ambigs = [c for c, v in self.moltype.ambiguities.items() if len(v) > 1]
            exclude_chars.update(ambigs)

        result = []
        for i in range(0, align_len - motif_length + 1, motif_length):
            counts = CategoryCounter([s[i : i + motif_length] for s in data])
            all_motifs.update(list(counts))
            result.append(counts)

        if all_motifs:
            alpha += tuple(sorted(set(alpha) ^ all_motifs))

        if exclude_chars:
            # this additional clause is required for the bytes moltype
            # That moltype includes '-' as a character
            alpha = [m for m in alpha if not (set(m) & exclude_chars)]

        for i, counts in enumerate(result):
            result[i] = counts.tolist(alpha)

        return MotifCountsArray(result, alpha)

    def probs_per_pos(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        warn: bool = False,
    ) -> MotifFreqsArray:
        """returns MotifFreqsArray per position"""
        counts = self.counts_per_pos(
            motif_length=motif_length,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
            warn=warn,
        )
        return counts.to_freq_array()

    def entropy_per_pos(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        warn: bool = False,
    ) -> numpy.ndarray:
        """returns shannon entropy per position"""
        probs = self.probs_per_pos(
            motif_length=motif_length,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
            warn=warn,
        )
        return probs.entropy()

    def counts_per_seq(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = False,
        warn: bool = False,
    ) -> MotifCountsArray:
        """counts of non-overlapping motifs per sequence

        Parameters
        ----------
        motif_length
            number of elements per character.
        include_ambiguity
            if True, motifs containing ambiguous characters
            from the seq moltype are included. No expansion of those is attempted.
        allow_gap
            if True, motifs containing a gap character are included.
        exclude_unobserved
            if False, all canonical states included
        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns
        """
        length = (len(self) // motif_length) * motif_length
        if warn and len(self) != length:
            warnings.warn(f"trimmed {len(self) - length}", UserWarning, stacklevel=2)

        counts = []
        motifs = set()
        for name in self.names:
            seq = self.get_gapped_seq(name)
            c = seq.counts(
                motif_length=motif_length,
                include_ambiguity=include_ambiguity,
                allow_gap=allow_gap,
            )
            motifs.update(c.keys())
            counts.append(c)

        # if type motifs not same as type element in moltype
        if not exclude_unobserved:
            motifs.update(self.moltype.alphabet.get_kmer_alphabet(motif_length))

        motifs = sorted(motifs)
        if not motifs:
            return None

        for i, c in enumerate(counts):
            counts[i] = c.tolist(motifs)
        return MotifCountsArray(counts, motifs, row_indices=self.names)

    def probs_per_seq(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = False,
        warn: bool = False,
    ) -> MotifFreqsArray:
        """return MotifFreqsArray per sequence

        Parameters
        ----------
        motif_length
            number of characters per tuple.
        include_ambiguity
            if True, motifs containing ambiguous characters
            from the seq moltype are included. No expansion of those is attempted.
        allow_gap
            if True, motifs containing a gap character are included.
        exclude_unobserved
            if True, unobserved motif combinations are excluded.
        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns
        """

        counts = self.counts_per_seq(
            motif_length=motif_length,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
            exclude_unobserved=exclude_unobserved,
            warn=warn,
        )
        return None if counts is None else counts.to_freq_array()

    def entropy_per_seq(
        self,
        motif_length: int = 1,
        include_ambiguity: bool = False,
        allow_gap: bool = False,
        exclude_unobserved: bool = True,
        warn: bool = False,
    ) -> numpy.ndarray:
        """returns the Shannon entropy per sequence

        Parameters
        ----------
        motif_length
            number of characters per tuple.
        include_ambiguity
            if True, motifs containing ambiguous characters
            from the seq moltype are included. No expansion of those is attempted.
        allow_gap
            if True, motifs containing a gap character are included.
        exclude_unobserved
            if True, unobserved motif combinations are excluded.
        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns

        Notes
        -----
        For motif_length > 1, it's advisable to specify exclude_unobserved=True,
        this avoids unnecessary calculations.
        """

        probs = self.probs_per_seq(
            motif_length=motif_length,
            include_ambiguity=include_ambiguity,
            allow_gap=allow_gap,
            exclude_unobserved=exclude_unobserved,
            warn=warn,
        )
        return None if probs is None else probs.entropy()

    def count_gaps_per_pos(self, include_ambiguity: bool = True) -> DictArray:
        """return counts of gaps per position as a DictArray

        Parameters
        ----------
        include_ambiguity
            if True, ambiguity characters that include the gap state are
            included
        """
        gap_array = self.get_gap_array(include_ambiguity=include_ambiguity)
        darr = DictArrayTemplate(range(len(self)))

        result = gap_array.sum(axis=0)
        return darr.wrap(result)

    def count_gaps_per_seq(
        self,
        induced_by: bool = False,
        unique: bool = False,
        include_ambiguity: bool = True,
        drawable: bool = False,
    ) -> DictArray:
        """return counts of gaps per sequence as a DictArray

        Parameters
        ----------
        induced_by
            a gapped column is considered to be induced by a seq if the seq
            has a non-gap character in that column.
        unique
            count is limited to gaps uniquely induced by each sequence
        include_ambiguity
            if True, ambiguity characters that include the gap state are
            included
        drawable
            if True, resulting object is capable of plotting data via specified
            plot type 'bar', 'box' or 'violin'
        """
        from cogent3.draw.drawable import Drawable

        gap_array = self.get_gap_array(include_ambiguity=include_ambiguity)
        darr = DictArrayTemplate(self.names)

        if unique:
            # we identify cols with a single non-gap character
            gap_cols = gap_array.sum(axis=0) == self.num_seqs - 1
            gap_array = gap_array[:, gap_cols] == False  # noqa
        elif induced_by:
            # identify all columns with gap opposite
            gap_cols = gap_array.sum(axis=0) > 0
            gap_array = gap_array[:, gap_cols] == False  # noqa
        else:
            gap_cols = gap_array.sum(axis=0) > 0
            gap_array = gap_array[:, gap_cols]

        result = gap_array.sum(axis=1)
        result = darr.wrap(result)
        if drawable:
            drawable = drawable.lower()
            trace_name = pathlib.Path(self.source).name if self.source else None
            draw = Drawable("Gaps Per Sequence", showlegend=False)
            draw.layout |= {"yaxis": {"title": "Gap counts"}}
            if drawable == "bar":
                trace = UnionDict(type="bar", y=result.array, x=self.names)
            else:
                trace = UnionDict(
                    type=drawable,
                    y=result.array,
                    text=self.names,
                    name=trace_name,
                )

            draw.add_trace(trace)
            result = draw.bound_to(result)

        return result

    def count_ambiguous_per_seq(self) -> DictArray:
        """Return the counts of ambiguous characters per sequence as a DictArray."""

        gap_index = self.moltype.most_degen_alphabet().gap_index
        ambigs_pos = self.array_seqs > gap_index
        ambigs = ambigs_pos.sum(axis=1)

        return DictArray.from_array_names(ambigs, self.names)

    def variable_positions(
        self,
        include_gap_motif: bool = True,
        include_ambiguity: bool = False,
        motif_length: int = 1,
    ) -> tuple[int]:
        """Return a list of variable position indexes.

        Parameters
        ----------
        include_gap_motif
            if False, sequences with a gap motif in a column are ignored.
        include_ambiguity
            if True, all states are considered.
        motif_length
            if any position within a motif is variable, the entire motif is
            considered variable.

        Returns
        -------
        tuple of integers, if motif_length > 1, the returned positions are
        motif_length long sequential indices.

        Notes
        -----
        Truncates alignment to be modulo motif_length.
        """
        align_len = len(self) // motif_length * motif_length
        array_seqs = self.array_seqs[:, :align_len]
        # both gap and missing data are treated as gaps
        alpha = self.moltype.most_degen_alphabet()
        gap_index = alpha.gap_index
        missing_index = alpha.missing_index
        # at the individual position level
        if include_gap_motif and include_ambiguity:
            indices = _var_pos_allow_all(array_seqs)
        elif include_gap_motif:
            indices = _var_pos_canonical_or_gap(array_seqs, gap_index, missing_index)
        elif include_ambiguity:
            indices = _var_pos_not_gap(array_seqs, gap_index)
        else:
            indices = _var_pos_canonical(array_seqs, gap_index)

        if motif_length > 1:
            indices = indices.reshape(-1, motif_length).any(axis=1).repeat(motif_length)

        return tuple(numpy.where(indices)[0].tolist())

    def omit_bad_seqs(self, quantile: OptFloat = None):
        """Returns new alignment without sequences with a number of uniquely
        introduced gaps exceeding quantile

        Uses count_gaps_per_seq(unique=True) to obtain the counts of gaps
        uniquely introduced by a sequence. The cutoff is the quantile of
        this distribution.

        Parameters
        ----------
        quantile
            sequences whose unique gap count is in a quantile larger than this
            cutoff are excluded. The default quantile is (num_seqs - 1) / num_seqs
        """
        gap_counts = self.count_gaps_per_seq(unique=True)
        quantile = quantile or (self.num_seqs - 1) / self.num_seqs
        cutoff = numpy.quantile(gap_counts.array, quantile)
        names = [name for name, count in gap_counts.items() if count <= cutoff]
        return self.take_seqs(names)

    def degap(self, storage_backend: str | None = None, **kwargs) -> SequenceCollection:
        """returns collection sequences without gaps or missing characters.

        Parameters
        ----------
        storage_backend
            name of the storage backend to use for the SeqsData object, defaults to
            cogent3 builtin.
        kwargs
            keyword arguments for the storage driver

        Notes
        -----
        The returned collection will not retain an annotation_db if present.
        """
        # because SequenceCollection does not track slice operations, we need
        # to apply any slice record to the underlying data
        sr = self._slice_record
        data, kw = self._seqs_data.get_ungapped(
            name_map=self._name_map,
            start=sr.plus_start,
            stop=sr.plus_stop,
            step=sr.plus_step,
        )
        kwargs = kw | kwargs
        # the SeqsData classes will return the data corresponding to the slice,
        # however, will not complement the data if the step is negative. We do
        # this here.
        rev_complement = self.moltype.rc
        data = (
            {name: rev_complement(seq) for name, seq in data.items()}
            if sr.step < 0
            else data
        )
        kwargs["annotation_db"] = self._annotation_db
        kwargs["storage_backend"] = storage_backend
        return make_unaligned_seqs(data, moltype=self.moltype, info=self.info, **kwargs)

    def get_degapped_relative_to(self, name: str) -> typing_extensions.Self:
        """Remove all columns with gaps in sequence with given name.

        Parameters
        ----------
        name
            sequence name

        Notes
        -----
        The returned alignment will not retain an annotation_db if present.
        """

        if name not in self.names:
            msg = f"Alignment missing sequence named {name!r}"
            raise ValueError(msg)

        gapindex = self.moltype.most_degen_alphabet().gap_index
        seqindex = self.names.index(name)
        indices = self.array_seqs[seqindex] != gapindex
        new = self.array_seqs[:, indices]

        new_seq_data = self._seqs_data.from_names_and_array(
            names=self._name_map.values(),
            data=new,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        kwargs = self._get_init_kwargs()
        kwargs["seqs_data"] = new_seq_data
        kwargs.pop("annotation_db", None)
        kwargs.pop("slice_record", None)
        return self.__class__(**kwargs)

    def matching_ref(
        self, ref_name: str, gap_fraction: float, gap_run: int
    ) -> typing_extensions.Self:
        """Returns new alignment with seqs well aligned with a reference.

        Parameters
        ----------
        ref_name
            name of the sequence to use as the reference
        gap_fraction
            fraction of positions that either have a gap in the
            template but not in the seq or in the seq but not in the template
        gap_run
            number of consecutive gaps tolerated in query relative to
            sequence or sequence relative to query
        """
        template = self.seqs[ref_name]
        gap_filter = make_gap_filter(template, gap_fraction, gap_run)
        return self.take_seqs_if(gap_filter)

    def sliding_windows(
        self,
        window: int,
        step: int,
        start: OptInt = None,
        end: OptInt = None,
    ) -> typing.Generator[typing_extensions.Self, None, None]:
        """Generator yielding new alignments of given length and interval.

        Parameters
        ----------
        window
            The length of each returned alignment.
        step
            The interval between the start of the successive windows.
        start
            first window start position
        end
            last window start position
        """
        start = start or 0
        end = [end, len(self) - window + 1][end is None]
        end = min(len(self) - window + 1, end)
        if start < end and len(self) - end >= window - 1:
            for pos in range(start, end, step):
                yield self[pos : pos + window]

    def gapped_by_map(self, keep: FeatureMap, **kwargs) -> typing_extensions.Self:
        # refactor: docstring
        # TODO: kath, not explicitly tested
        seqs = {}
        for seq in self.seqs:
            selected = seq[keep]
            seqs[self.name_map[seq.name]] = numpy.array(selected.gapped_seq)

        seqs_data = self._seqs_data.from_seqs(
            data=seqs,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        init_kwargs = self._get_init_kwargs()
        init_kwargs.pop("annotation_db", None)
        init_kwargs |= kwargs
        init_kwargs["seqs_data"] = seqs_data
        init_kwargs.pop("slice_record", None)
        return self.__class__(**init_kwargs)

    def filtered(
        self,
        predicate: typing.Callable[[typing_extensions.Self], bool],
        motif_length: int = 1,
        drop_remainder: bool = True,
        **kwargs,
    ) -> typing_extensions.Self:
        """The alignment positions where predicate(column) is true.

        Parameters
        ----------
        predicate
            a callback function that takes an tuple of motifs and returns
            True/False
        motif_length
            length of the motifs the sequences should be split  into, eg. 3 for
            filtering aligned codons.
        drop_remainder
            If length is not modulo motif_length, allow dropping the terminal
            remaining columns
        """
        # refactor: type hint for predicate
        length = len(self)
        drop = length % motif_length
        if drop != 0 and not drop_remainder:
            msg = f"aligned length not divisible by motif_length={motif_length}"
            raise ValueError(
                msg,
            )
        length -= drop
        kept = numpy.zeros(length, dtype=bool)
        for pos in range(0, length, motif_length):
            seqs = [seq[pos : pos + motif_length] for seq in self.seqs]
            if predicate(seqs):
                kept[pos : pos + motif_length] = True

        indices = numpy.where(kept)[0]

        return self.take_positions(indices.tolist())

    def no_degenerates(
        self, motif_length: int = 1, allow_gap: bool = False
    ) -> typing_extensions.Self:
        """returns new alignment without degenerate characters

        Parameters
        ----------
        motif_length
            sequences are segmented into units of this size and the segments are
            excluded if they contain degenerate characters.
        allow_gap
            whether gaps are allowed or whether they are treated as a degenerate
            character (latter is default, as most evolutionary modelling treats
            gaps as N).
        """
        if self.moltype.degen_alphabet is None:
            msg = (
                f"Invalid MolType={self.moltype.label} (no degenerate characters), "
                "create the alignment using DNA, RNA or PROTEIN"
            )
            raise c3_moltype.MolTypeError(
                msg,
            )

        chars = len(self.moltype)

        array_pos = self.array_positions
        # by design, char alphabets are organised such that the canonical
        # characters always occur first, followed by gap, then ambiguity
        # characters. so we can define a cutoff as follows:
        cutoff = chars + 1 if allow_gap else chars
        indices = (array_pos < cutoff).all(axis=1)

        if motif_length > 1:
            num_motif = len(self) // motif_length

            if remainder := len(self) % motif_length:
                indices = indices[:-remainder]
                array_pos = array_pos[:-remainder]

            motif_valid = indices.reshape(num_motif, motif_length).all(axis=1).flatten()
            indices = numpy.repeat(motif_valid, motif_length)

        selected = array_pos[indices].T

        aligned_seqs_data = self._seqs_data.from_names_and_array(
            names=self._name_map.values(),
            data=selected,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        kwargs = self._get_init_kwargs()
        kwargs["seqs_data"] = aligned_seqs_data
        kwargs.pop("annotation_db", None)
        kwargs.pop("slice_record", None)
        return self.__class__(**kwargs)

    def _omit_gap_pos_single(
        self,
        gap_index: int,
        missing_index: OptInt,
        allowed_num: int,
    ) -> numpy.ndarray[bool]:
        # for motif_length == 1
        indices = numpy.empty(len(self), dtype=bool)
        for i, col in enumerate(self.array_seqs.T):
            indices[i] = _gap_ok_vector_single(
                col,
                gap_index,
                missing_index,
                allowed_num,
            )
        return indices

    def _omit_gap_pos_multi(
        self,
        gap_index: int,
        missing_index: OptInt,
        allowed_num: int,
        motif_length: int,
    ) -> numpy.ndarray[bool]:
        # for motif_length > 1
        num_motifs = len(self) // motif_length
        indices = numpy.empty(num_motifs * motif_length, dtype=bool)
        for i in range(0, len(self), motif_length):
            col = self.array_seqs.T[i : i + motif_length].T
            ok = _gap_ok_vector_multi(
                col,
                gap_index,
                missing_index,
                motif_length,
                allowed_num,
            )
            indices[i : i + motif_length] = ok
        return indices

    def omit_gap_pos(
        self,
        allowed_gap_frac: float | None = None,
        motif_length: int = 1,
    ):
        """Returns new alignment where all cols (motifs) have <= allowed_gap_frac gaps.

        Parameters
        ----------
        allowed_gap_frac
            specifies proportion of gaps is allowed in each column. Set to 0 to
            exclude columns with any gaps, 1 to include all columns. Default is
            None which is equivalent to (num_seqs-1)/num_seqs and leads to
            elimination of columns that are only gaps.
        motif_length
            set's the "column" width, e.g. setting to 3 corresponds to codons.
            A motif that includes a gap at any position is included in the
            counting.
        """
        alpha = self.moltype.most_degen_alphabet()
        gap_index = alpha.gap_index
        missing_index = alpha.missing_index
        if not gap_index and not missing_index:
            return self

        allowed_num = (
            self.num_seqs - 1
            if allowed_gap_frac is None
            else int(numpy.floor(self.num_seqs * allowed_gap_frac))
        )

        # we are assuming gap and missing data have same value on other strand
        positions = self.array_positions
        if motif_length == 1:
            indices = self._omit_gap_pos_single(gap_index, missing_index, allowed_num)
        else:
            positions = positions[: motif_length * (len(positions) // motif_length)]
            indices = self._omit_gap_pos_multi(
                gap_index,
                missing_index,
                allowed_num,
                motif_length,
            )
        selected = positions[indices, :].T

        aligned_seqs_data = self._seqs_data.from_names_and_array(
            names=self._name_map.values(),
            data=selected,
            alphabet=alpha,
        )
        kwargs = self._get_init_kwargs()
        kwargs["seqs_data"] = aligned_seqs_data
        # slice record now needs to be reset since we likely have disjoint positions
        kwargs.pop("slice_record", None)
        return self.__class__(**kwargs)

    def has_terminal_stop(
        self, gc: c3_genetic_code.GeneticCodeChoiceType = 1, strict: bool = False
    ) -> bool:
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
            seq = self.seqs[seq_name].seq
            if seq.has_terminal_stop(gc=gc, strict=strict):
                return True
        return False

    def sample(
        self,
        *,
        n: int | None = None,
        with_replacement: bool = False,
        motif_length: int = 1,
        randint: typing.Callable[
            [int, int | None, int | None], numpy.ndarray
        ] = numpy.random.randint,
        permutation: typing.Callable[
            [numpy.ndarray], numpy.ndarray
        ] = numpy.random.permutation,
    ) -> typing_extensions.Self:
        """Returns random sample of positions from self, e.g. to bootstrap.

        Parameters
        ----------
        n
            number of positions to sample. If None, all positions are sampled.
        with_replacement
            if True, samples with replacement.
        motif_length
            number of positions to sample as a single motif.
        randint
            random number generator, default is numpy.randint
        permutation
            function to generate a random permutation of positions, default is
            numpy.permutation

        Notes
        -----
        By default (resampling all positions without replacement), generates
        a permutation of the positions of the alignment.

        Setting with_replacement to True and otherwise leaving parameters as
        defaults generates a standard bootstrap resampling of the alignment.
        """
        # refactor: type hint for randint, permutation
        # refactor: array
        #  Given array_pos property, it will be efficient to generate random
        # indices and use numpy.take() using that. In the case of motif_length
        # != 1, the total number of positions is just len(self) // motif_length.
        # Having produced that, those indices can be scaled back up, or the
        # numpy array reshaped.

        population_size = len(self) // motif_length
        if not with_replacement and n and n > population_size:
            msg = f"cannot sample without replacement when {n=} > {population_size=}"
            raise ValueError(msg)

        n = n or population_size

        if with_replacement:
            locations = randint(0, population_size, n)
        else:
            locations = permutation(population_size)[:n]

        if motif_length == 1:
            positions = locations
        else:
            positions = numpy.empty(n * motif_length, dtype=int)
            for i, loc in enumerate(locations):
                positions[i * motif_length : (i + 1) * motif_length] = range(
                    loc * motif_length,
                    (loc + 1) * motif_length,
                )

        return self.take_positions(positions)

    def distance_matrix(
        self,
        calc: str = "pdist",
        drop_invalid: bool = False,
        parallel: bool = False,
    ):
        """Returns pairwise distances between sequences.

        Parameters
        ----------
        calc
            a pairwise distance calculator name. Presently only
            'pdist', 'jc69', 'tn93', 'hamming', 'paralinear' are supported.
        drop_invalid
            If True, sequences for which a pairwise distance could not be
            calculated are excluded. If False, an ArithmeticError is raised if
            a distance could not be computed on observed data.
        """
        from cogent3.evolve.pairwise_distance_numba import get_distance_calculator

        calc = get_distance_calculator(calc)
        try:
            result = calc(self, invalid_raises=not drop_invalid, parallel=parallel)
        except ArithmeticError as e:
            msg = "not all pairwise distances could be computed, try drop_invalid=True"
            raise ArithmeticError(msg) from e

        if drop_invalid:
            result = result.drop_invalid()

        return result

    def quick_tree(
        self,
        calc: str = "pdist",
        drop_invalid: bool = False,
        parallel: bool = False,
        use_hook: str | None = None,
    ) -> PhyloNode:
        """Returns a phylogenetic tree.

        Parameters
        ----------
        calc
            a pairwise distance calculator or name of one. For options see
            cogent3.evolve.fast_distance.available_distances
        drop_invalid
            If True, sequences for which a pairwise distance could not be
            calculated are excluded. If False, an ArithmeticError is raised if
            a distance could not be computed on observed data.
        parallel
            parallel execution of distance calculations
        use_hook
            name of a third-party package that implements the quick_tree
            hook. If not specified, defaults to the first available hook or
            the cogent3 quick_tree() app. To force default, set
            use_hook="cogent3".

        Returns
        -------
        a phylogenetic tree
        """
        dm = self.distance_matrix(
            calc=calc,
            drop_invalid=drop_invalid,
            parallel=parallel,
        )
        return dm.quick_tree(use_hook=use_hook)

    def trim_stop_codons(
        self,
        gc: c3_genetic_code.GeneticCodeChoiceType = 1,
        strict: bool = False,
        **kwargs,
    ) -> typing_extensions.Self:
        # refactor: array
        if not self.has_terminal_stop(gc=gc, strict=strict):
            return self

        # define a regex for finding stop codons followed by terminal gaps
        gc = c3_genetic_code.get_code(gc)
        gaps = "".join(self.moltype.gaps)
        pattern = f"({'|'.join(gc['*'])})[{gaps}]*$"
        terminal_stop = re.compile(pattern)

        data = self.to_dict()
        result = {}
        for name, seq in data.items():
            if match := terminal_stop.search(seq):
                diff = len(seq) - match.start()
                seq = terminal_stop.sub("-" * diff, seq)

            result[self.name_map[name]] = seq

        seqs_data = self._seqs_data.from_seqs(
            data=result,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        init_kwargs = self._get_init_kwargs()
        init_kwargs["seqs_data"] = seqs_data
        init_kwargs.pop("slice_record", None)
        init_kwargs |= kwargs
        return self.__class__(**init_kwargs)

    @extend_docstring_from(SequenceCollection.get_translation)
    def get_translation(
        self,
        gc: c3_genetic_code.GeneticCodeChoiceType = 1,
        incomplete_ok: bool = False,
        include_stop: bool = False,
        trim_stop: bool = True,
        **kwargs,
    ) -> typing_extensions.Self:
        if not self.moltype.is_nucleic:
            msg = f"moltype must be a DNA/RNA, not {self.moltype.name!r}"
            raise c3_moltype.MolTypeError(msg)

        if not trim_stop or include_stop:
            seqs = self
        else:
            seqs = self.trim_stop_codons(gc=gc, strict=not incomplete_ok)

        translated = {}
        for seqname in seqs.names:
            seq = seqs.get_gapped_seq(seqname)
            pep = seq.get_translation(
                gc,
                incomplete_ok=incomplete_ok,
                include_stop=include_stop,
                trim_stop=trim_stop,
            )
            translated[self.name_map[seqname]] = numpy.array(pep)

        pep_moltype = c3_moltype.get_moltype(
            "protein_with_stop" if include_stop else "protein",
        )
        seqs_data = self._seqs_data.from_seqs(
            data=translated,
            alphabet=pep_moltype.most_degen_alphabet(),
            offset=None,
        )
        return self.__class__(
            seqs_data=seqs_data,
            moltype=pep_moltype,
            name_map=self._name_map,
            info=self.info,
            source=self.source,
            **kwargs,
        )

    def make_feature(
        self,
        *,
        feature: FeatureDataType,
        on_alignment: OptBool = None,
    ) -> Feature:
        """
        create a feature on named sequence, or on the alignment itself

        Parameters
        ----------
        feature
            a dict with all the necessary data rto construct a feature
        on_alignment
            the feature is in alignment coordinates, incompatible with setting
            'seqid'. Set to True if 'seqid' not provided.

        Returns
        -------
        Feature

        Raises
        ------
        ValueError if define a 'seqid' not on alignment or use 'seqid' and
        on_alignment.

        Notes
        -----
        To get a feature AND add it to annotation_db, use add_feature().
        """
        if on_alignment is None:
            on_alignment = feature.pop("on_alignment", None)

        if not on_alignment and feature["seqid"]:
            return self.seqs[feature["seqid"]].make_feature(feature, self)

        feature["seqid"] = feature.get("seqid", None)
        # there's no sequence to bind to, the feature is directly on self
        revd = Strand.from_value(feature.pop("strand", None)) is Strand.MINUS
        feature["strand"] = Strand.MINUS.value if revd else Strand.PLUS.value
        fmap = FeatureMap.from_locations(
            locations=feature.pop("spans"),
            parent_length=len(self),
        )
        if revd:
            fmap = fmap.nucleic_reversed()
        return Feature(parent=self, map=fmap, **feature)

    def add_feature(
        self,
        *,
        biotype: str,
        name: str,
        spans: list[tuple[int, int]],
        seqid: OptStr = None,
        parent_id: OptStr = None,
        strand: str = "+",
        on_alignment: OptBool = None,
    ) -> Feature:
        """
        add feature on named sequence, or on the alignment itself

        Parameters
        ----------
        seqid
            sequence name, incompatible with on_alignment
        parent_id
            name of the parent feature
        biotype
            biological type, e.g. CDS
        name
            name of the feature
        spans
            plus strand coordinates of feature
        strand
            '+' (default) or '-'
        on_alignment
            the feature is in alignment coordinates, incompatible with setting
            seqid. Set to True if seqid not provided.

        Returns
        -------
        Feature

        Raises
        ------
        ValueError if define a seqid not on alignment or use seqid and
        on_alignment.
        """
        if seqid and on_alignment is None:
            on_alignment = False
        elif not on_alignment:
            on_alignment = on_alignment is None

        if seqid and on_alignment:
            msg = "seqid and on_alignment are incomatible"
            raise ValueError(msg)

        if seqid and seqid not in self.names:
            msg = f"unknown {seqid=}"
            raise ValueError(msg)

        feature = {k: v for k, v in locals().items() if k != "self"}
        feature["strand"] = Strand.from_value(strand).value
        # property ensures db is created
        self.annotation_db.add_feature(**feature)
        for discard in ("on_alignment", "parent_id"):
            feature.pop(discard, None)
        return self.make_feature(feature=feature, on_alignment=on_alignment)

    def _get_seq_features(
        self,
        *,
        seqid: str | None = None,
        biotype: str | None = None,
        name: str | None = None,
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
        if not self._annotation_db:
            return None

        seqid_to_seqname = {v: k for k, v in self._name_map.items()}

        seqids = [seqid] if isinstance(seqid, str) else seqid
        if seqids is None:
            seqids = tuple(seqid_to_seqname)
        elif set(seqids) & set(self.names):
            # we've been given seq names, convert to parent names
            seqids = [self.seqs[seqid].parent_coordinates()[0] for seqid in seqids]
        elif seqids and set(seqids) <= seqid_to_seqname.keys():
            # already correct
            pass
        else:
            msg = f"unknown {seqid=}"
            raise ValueError(msg)

        for seqid in seqids:
            seqname = seqid_to_seqname[seqid]
            seq = self.seqs[seqname]
            # we use parent seqid
            parent_id, start, stop, _ = seq.parent_coordinates()
            offset = seq.data.offset

            for feature in self.annotation_db.get_features_matching(
                seqid=parent_id,
                biotype=biotype,
                name=name,
                on_alignment=False,
                allow_partial=allow_partial,
                start=start,
                stop=stop,
            ):
                if offset:
                    feature["spans"] = (numpy.array(feature["spans"]) - offset).tolist()
                # passing self only used when self is an Alignment
                yield seq.make_feature(feature, self)

    def get_features(
        self,
        *,
        seqid: str | None = None,
        biotype: str | None = None,
        name: str | None = None,
        on_alignment: bool | None = None,
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
        on_alignment
            limit query to features on Alignment, ignores sequences. Ignored on
            SequenceCollection instances.
        allow_partial
            allow features partially overlaping self

        Notes
        -----
        When dealing with a nucleic acid moltype, the returned features will
        yield a sequence segment that is consistently oriented irrespective
        of strand of the current instance.
        """
        if not self._annotation_db or not len(self._annotation_db):
            return None

        # we only do on-alignment in here
        if not on_alignment:
            local_vars = locals()
            kwargs = {k: v for k, v in local_vars.items() if k != "self"}
            kwargs.pop("on_alignment")
            yield from self._get_seq_features(**kwargs)

        if on_alignment == False:  # noqa
            return

        seq_map = None
        for feature in self.annotation_db.get_features_matching(
            biotype=biotype,
            name=name,
            on_alignment=on_alignment,
            allow_partial=allow_partial,
        ):
            if feature["seqid"]:
                continue
            on_al = feature.pop("on_alignment", on_alignment)
            if feature["seqid"]:
                msg = f"{on_alignment=} {feature=}"
                raise RuntimeError(msg)
            if seq_map is None:
                seq_map = self.seqs[0].map.to_feature_map()
                *_, strand = self.seqs[0].seq.parent_coordinates()
            else:
                strand = feature.pop("strand", None)

            spans = numpy.array(feature["spans"])
            spans = seq_map.relative_position(spans)
            feature["spans"] = spans.tolist()
            # and if i've been reversed...?
            feature["strand"] = Strand.from_value(strand).value
            yield self.make_feature(feature=feature, on_alignment=on_al)

    def get_projected_feature(self, *, seqid: str, feature: Feature) -> Feature:
        """returns an alignment feature projected onto the seqid sequence

        Parameters
        ----------
        seqid
            name of the sequence to project the feature onto
        feature
            a Feature, bound to self, that will be projected

        Returns
        -------
        a new Feature bound to seqid

        Notes
        -----
        The alignment coordinates of feature are converted into the seqid
        sequence coordinates and the object is bound to that sequence.

        The feature is added to the annotation_db.
        """
        target_aligned = self.seqs[seqid]
        if feature.parent is not self:
            msg = "Feature does not belong to this alignment"
            raise ValueError(msg)
        result = feature.remapped_to(target_aligned.seq, target_aligned.map)
        # property ensures db is created
        self.annotation_db.add_feature(**feature.to_dict())
        return result

    def get_projected_features(self, *, seqid: str, **kwargs) -> list[Feature]:
        """projects all features from other sequences onto seqid"""
        annots = []
        for name in self.names:
            if name == seqid:
                continue
            annots.extend(list(self.get_features(seqid=name, **kwargs)))
        return [self.get_projected_feature(seqid=seqid, feature=a) for a in annots]

    def get_drawables(
        self,
        *,
        biotype: Optional[str, typing.Iterable[str]] = None,
    ) -> dict:
        """returns a dict of drawables, keyed by type

        Parameters
        ----------
        biotype
            passed to get_features(biotype). Can be a single biotype or
            series. Only features matching this will be included.
        """
        result = defaultdict(list)
        for f in self.get_features(biotype=biotype, allow_partial=True):
            result[f.biotype].append(f.get_drawable())
        return result

    def get_drawable(
        self,
        *,
        biotype: Optional[str, typing.Iterable[str]] = None,
        width: int = 600,
        vertical: int = False,
        title: OptStr = None,
    ):
        """make a figure from sequence features

        Parameters
        ----------
        biotype
            passed to get_features(biotype). Can be a single biotype or
            series. Only features matching this will be included.
        width
            width in pixels
        vertical
            rotates the drawable
        title
            title for the plot

        Returns
        -------
        a Drawable instance
        """
        # TODO gah I think this needs to be modified to make row-blocks
        # for each sequence in the alignment, or are we displaying the
        # sequence name in the feature label?
        from cogent3.draw.drawable import Drawable

        drawables = self.get_drawables(biotype=biotype)
        if not drawables:
            return None
        # we order by tracks
        top = 0
        space = 0.25
        annotes = []
        for feature_type in drawables:
            new_bottom = top + space
            for i, annott in enumerate(drawables[feature_type]):
                annott.shift(y=new_bottom - annott.bottom)
                if i > 0:
                    # refactor: design
                    # modify the api on annott, we should not be using
                    # a private attribute!
                    annott._showlegend = False
                annotes.append(annott)

            top = annott.top

        top += space
        height = max((top / len(self)) * width, 300)
        xaxis = {"range": [0, len(self)], "zeroline": False, "showline": True}
        yaxis = {
            "range": [0, top],
            "visible": False,
            "zeroline": True,
            "showline": True,
        }

        if vertical:
            all_traces = [t.T.as_trace() for t in annotes]
            width, height = height, width
            xaxis, yaxis = yaxis, xaxis
        else:
            all_traces = [t.as_trace() for t in annotes]

        drawer = Drawable(title=title, traces=all_traces, width=width, height=height)
        drawer.layout.update(xaxis=xaxis, yaxis=yaxis)
        return drawer

    def seqlogo(
        self,
        width: float = 700,
        height: float = 100,
        wrap: OptInt = None,
        vspace: float = 0.005,
        colours: dict | None = None,
    ):
        """returns Drawable sequence logo using mutual information

        Parameters
        ----------
        width, height
            plot dimensions in pixels
        wrap
            number of alignment columns per row
        vspace
            vertical separation between rows, as a proportion of total plot
        colours
            mapping of characters to colours. If note provided, defaults to
            custom for everything ecept protein, which uses protein moltype
            colours.

        Notes
        -----
        Computes MI based on log2 and includes the gap state, so the maximum
        possible value is -log2(1/num_states)
        """
        assert 0 <= vspace <= 1, "vertical space must be in range 0, 1"
        freqs = self.counts_per_pos(allow_gap=True).to_freq_array()
        if colours is None and "protein" in self.moltype.label:
            colours = self.moltype._colors

        return freqs.logo(
            width=width,
            height=height,
            wrap=wrap,
            vspace=vspace,
            colours=colours,
        )

    def coevolution(
        self,
        stat: str = "nmi",
        segments: list[tuple[int, int]] | None = None,
        drawable: OptStr = None,
        show_progress: bool = False,
        parallel: bool = False,
        par_kw: OptDict = None,
    ):
        """performs pairwise coevolution measurement

        Parameters
        ----------
        stat
            coevolution metric, defaults to 'nmi' (Normalized Mutual
            Information). Valid choices are 'rmi' (Resampled Mutual Information)
            and 'mi', mutual information.
        segments
            coordinates of the form [(start, end), ...] where all possible
            pairs of alignment positions within and between segments are
            examined.
        drawable
            Result object is capable of plotting data specified type. str value
            must be one of plot type 'box', 'heatmap', 'violin'.
        show_progress
            shows a progress bar.
        parallel
            run in parallel, according to arguments in par_kwargs.
        par_kw
            dict of values for configuring parallel execution.

        Returns
        -------
        DictArray of results with lower-triangular values. Upper triangular
        elements and estimates that could not be computed for numerical reasons
        are set as nan
        """
        from cogent3.draw.drawable import AnnotatedDrawable, Drawable
        from cogent3.evolve import coevolution as coevo
        from cogent3.util.union_dict import UnionDict

        # refactor: design
        # These graphical representations of matrices should be separate functions
        # in the drawing submodule somewhere
        stat = stat.lower()
        if segments:
            segments = [range(*segment) for segment in segments]

        result = coevo.coevolution_matrix(
            alignment=self,
            stat=stat,
            positions=segments,
            show_progress=show_progress,
            parallel=parallel,
            par_kw=par_kw,
        )
        if drawable is None:
            return result

        trace_name = pathlib.Path(self.source).name if self.source else None
        if drawable in ("box", "violin"):
            trace = UnionDict(
                type=drawable,
                y=result.array.flatten(),
                showlegend=False,
                name="",
            )
            draw = Drawable(
                width=500,
                height=500,
                title=trace_name,
                ytitle=stat.upper(),
            )
            draw.add_trace(trace)
            result = draw.bound_to(result)
        elif drawable:
            axis_title = "Alignment Position"
            axis_args = {
                "showticklabels": True,
                "mirror": True,
                "showgrid": False,
                "showline": True,
                "zeroline": False,
            }
            height = 500
            width = height
            draw = Drawable(
                width=width,
                height=height,
                xtitle=axis_title,
                ytitle=axis_title,
                title=trace_name,
            )

            trace = UnionDict(
                type="heatmap",
                z=result.array,
                colorbar={"title": {"text": stat.upper(), "font": {"size": 16}}},
            )
            draw.add_trace(trace)
            draw.layout.xaxis.update(axis_args)
            draw.layout.yaxis.update(axis_args)

            try:
                bottom = self.get_drawable()
                left = self.get_drawable(vertical=True)
            except AttributeError:
                bottom = False

            if bottom and drawable != "box":
                xlim = 1.2
                draw.layout.width = height * xlim
                layout = {"legend": {"x": xlim, "y": 1}}
                draw = AnnotatedDrawable(
                    draw,
                    left_track=left,
                    bottom_track=bottom,
                    xtitle=axis_title,
                    ytitle=axis_title,
                    xrange=[0, len(self)],
                    yrange=[0, len(self)],
                    layout=layout,
                )

            result = draw.bound_to(result)

        return result

    def information_plot(
        self,
        width: OptInt = None,
        height: OptInt = None,
        window: OptInt = None,
        stat: str = "median",
        include_gap: bool = True,
    ):
        """plot information per position

        Parameters
        ----------
        width
            figure width in pixels
        height
            figure height in pixels
        window
            used for smoothing line, defaults to sqrt(length)
        stat
            'mean' or 'median, used as the summary statistic for each window
        include_gap
            whether to include gap counts, shown on right y-axis
        """
        from cogent3.draw.drawable import AnnotatedDrawable, Drawable

        # refactor: design
        # These graphical representations of matrices should be separate functions
        # in the drawing submodule somewhere
        window = window or numpy.sqrt(len(self))
        window = int(window)
        y = self.entropy_per_pos()
        nan_indices = numpy.isnan(y)
        if nan_indices.sum() == y.shape[0]:  # assuming 1D array
            y.fill(0.0)
        max_entropy = y[nan_indices == False].max()  # noqa
        y = max_entropy - y  # convert to information
        # now make all nan's 0
        y[nan_indices] = 0
        stats = {"mean": numpy.mean, "median": numpy.median}
        if stat not in stats:
            msg = 'stat must be either "mean" or "median"'
            raise ValueError(msg)
        calc_stat = stats[stat]
        num = len(y) - window
        v = [calc_stat(y[i : i + window]) for i in range(num)]
        x = numpy.arange(num)
        x += window // 2  # shift x-coordinates to middle of window
        trace_line = UnionDict(
            type="scatter",
            x=x,
            y=v,
            mode="lines",
            name=f"smoothed {stat}",
            line={"shape": "spline", "smoothing": 1.3},
        )
        trace_marks = UnionDict(
            type="scatter",
            x=numpy.arange(y.shape[0]),
            y=y,
            mode="markers",
            opacity=0.5,
            name="per position",
        )
        layout = UnionDict(
            title="Information per position",
            width=width,
            height=height,
            showlegend=True,
            yaxis={"range": [0, max(y) * 1.2], "showgrid": False},
            xaxis={
                "showgrid": False,
                "range": [0, len(self)],
                "mirror": True,
                "showline": True,
            },
        )

        traces = [trace_marks, trace_line]
        if include_gap:
            gap_counts = self.count_gaps_per_pos()
            y = [calc_stat(gap_counts[i : i + window]) for i in range(num)]
            trace_g = UnionDict(
                type="scatter",
                x=x,
                y=y,
                yaxis="y2",
                name="Gaps",
                mode="lines",
                line={"shape": "spline", "smoothing": 1.3},
            )
            traces += [trace_g]
            layout.yaxis2 = {
                "title": "Count",
                "side": "right",
                "overlaying": "y",
                "range": [0, max(gap_counts) * 1.2],
                "showgrid": False,
                "showline": True,
            }

        draw = Drawable(
            title="Information per position",
            xtitle="Position",
            ytitle=f"Information (window={window})",
        )
        draw.traces.extend(traces)
        draw.layout |= layout
        draw.layout.legend = {"x": 1.1, "y": 1}

        try:
            drawable = self.get_drawable()
        except AttributeError:
            drawable = False

        if drawable:
            draw = AnnotatedDrawable(
                draw,
                bottom_track=drawable,
                xtitle="position",
                ytitle=f"Information (window={window})",
                layout=layout,
            )

        return draw

    def to_phylip(self) -> str:
        """
        Return collection in PHYLIP format and mapping to sequence ids

        Notes
        -----
        raises exception if sequences do not all have the same length
        """
        phylip = cogent3._plugin.get_seq_format_writer_plugin(format_name="phylip")  # noqa: SLF001
        return phylip.formatted(self)

    def to_pretty(
        self,
        name_order: list[str] | None = None,
        wrap: int | None = None,
    ) -> str:
        """returns a string representation of the alignment in pretty print format

        Parameters
        ----------
        name_order
            order of names for display.
        wrap
            maximum number of printed bases
        """
        names, output = self._get_raw_pretty(name_order=name_order)
        label_width = max(list(map(len, names)))
        name_template = "{:>%d}" % label_width
        display_names = {n: name_template.format(n) for n in names}

        def make_line(label, seq) -> str:
            return f"{label}    {seq}"

        if wrap is None:
            result = [make_line(display_names[n], "".join(output[n])) for n in names]
            return "\n".join(result)

        align_length = len(self)
        result = []
        for start in range(0, align_length, wrap):
            for n in names:
                result.append(
                    make_line(
                        display_names[n],
                        "".join(output[n][start : start + wrap]),
                    ),
                )

            result.append("")

        if not result[-1]:
            del result[-1]

        return "\n".join(result)

    def to_html(
        self,
        name_order: typing.Sequence[str] | None = None,
        wrap: int = 60,
        limit: int | None = None,
        ref_name: str = "longest",
        colors: Mapping[str, str] | None = None,
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
        ref_name
            Name of an existing sequence or 'longest'. If the latter, the
            longest sequence (excluding gaps and ambiguities) is selected as the
            reference.
        colors
            {character
            moltype.
        font_size
            in points. Affects labels and sequence and line spacing
            (proportional to value)
        font_family
            string denoting font family

        Examples
        --------

        In a jupyter notebook, this code is used to provide the representation.

        .. code-block:: python

            aln  # is rendered by jupyter

        You can directly use the result for display in a notebook as

        .. code-block:: python

            from IPython.core.display import HTML

            HTML(aln.to_html())
        """
        css, styles = self.moltype.get_css_style(
            colors=colors,
            font_size=font_size,
            font_family=font_family,
        )
        if name_order:
            selected = self.take_seqs(name_order)
            name_order = list(name_order)
        else:
            name_order = list(self.names)
            ref_name = ref_name or "longest"
            selected = self

        if ref_name == "longest":
            lengths = selected.get_lengths(include_ambiguity=False, allow_gap=False)

            length_names = defaultdict(list)
            for n, l in lengths.items():
                length_names[l].append(n)

            longest = max(length_names)
            ref = sorted(length_names[longest])[0]

        elif ref_name:
            if ref_name not in selected.names:
                msg = f"Unknown sequence name {ref_name}"
                raise ValueError(msg)
            ref = ref_name

        name_order.remove(ref)
        name_order.insert(0, ref)

        if limit is None:
            names, output = selected._get_raw_pretty(name_order)
        else:
            names, output = selected[:limit]._get_raw_pretty(name_order)

        refname = names[0]
        refseq = output[refname]
        seqlen = len(refseq)

        if selected.moltype.gaps:
            gaps = "".join(selected.moltype.gaps)
            start_gap = re.search(f"^[{gaps}]+", "".join(refseq))
            end_gap = re.search(f"[{gaps}]+$", "".join(refseq))
            start = 0 if start_gap is None else start_gap.end()
            end = seqlen if end_gap is None else end_gap.start()
        else:
            start = 0
            end = seqlen

        seq_style = []
        template = '<span class="%s">%%s</span>'
        styled_seqs = defaultdict(list)
        for i in range(seqlen):
            char = refseq[i]
            if i < start or i >= end:
                style = f"terminal_ambig_{selected.moltype.label}"
            else:
                style = styles[char]

            seq_style.append(template % style)
            styled_seqs[refname].append(seq_style[-1] % char)

        for name in names:
            if name == refname:
                continue

            seq = []
            for i, c in enumerate(output[name]):
                if c == ".":
                    s = seq_style[i] % c
                else:
                    s = template % (styles[c])
                    s = s % c
                seq.append(s)

            styled_seqs[name] = seq

        # make a html table
        seqs = numpy.array([styled_seqs[n] for n in names], dtype="O")
        table = ["<table>"]
        seq_ = "<td>%s</td>"
        label_ = '<td class="label">%s</td>'
        num_row_ = '<tr class="num_row"><td></td><td><b>{:,d}</b></td></tr>'
        for i in range(0, seqlen, wrap):
            table.append(num_row_.format(i))
            seqblock = seqs[:, i : i + wrap].tolist()
            for n, s in zip(names, seqblock, strict=False):
                s = "".join(s)
                row = "".join([label_ % n, seq_ % s])
                table.append(f"<tr>{row}</tr>")
        table.append("</table>")
        if (limit and limit < len(selected)) or (
            name_order and len(name_order) < len(selected.names)
        ):
            summary = (
                f"{self.num_seqs} x {len(self)} (truncated to "
                f"{len(name_order) if name_order else len(selected.names)} x "
                f"{limit or len(selected)}) {selected.moltype.label} alignment"
            )
        else:
            summary = (
                f"{self.num_seqs} x {len(self)} {selected.moltype.label} alignment"
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
            "\n".join([f".c3align {style}" for style in css]),
            "</style>",
            '<div class="c3align">',
            "\n".join(table),
            f"<p><i>{summary}</i></p>",
            "</div>",
        ]
        return "\n".join(text)

    def _get_raw_pretty(self, name_order):
        """returns dict {name: seq, ...} for pretty print"""
        if name_order is not None:
            assert set(name_order) <= set(self.names), "names don't match"

        output = defaultdict(list)
        names = name_order or self.names
        num_seqs = len(names)

        seqs = [str(self.seqs[name]) for name in names]
        positions = list(zip(*seqs, strict=False))

        for position in positions:
            ref = position[0]
            output[names[0]].append(ref)
            for seq_num in range(1, num_seqs):
                val = "." if position[seq_num] == ref else position[seq_num]
                output[names[seq_num]].append(val)

        return names, output

    def _repr_html_(self) -> str:
        settings = self._repr_policy.copy()
        env_vals = get_setting_from_environ(
            "COGENT3_ALIGNMENT_REPR_POLICY",
            {"num_seqs": int, "num_pos": int, "wrap": int, "ref_name": str},
        )
        settings.update(env_vals)
        return self.to_html(
            name_order=self.names[: settings["num_seqs"]],
            ref_name=settings["ref_name"],
            limit=settings["num_pos"],
            wrap=settings["wrap"],
        )

    def _mapped(self, slicemap) -> typing_extensions.Self:
        seqs = {}
        maps = {}
        for aligned in self.seqs:
            seq, im = aligned.slice_with_map(slicemap)
            name = self.name_map[aligned.name]
            seqs[name] = seq
            maps[name] = im

        data = self._seqs_data.from_seqs_and_gaps(
            seqs=seqs,
            gaps=maps,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        kwargs = self._get_init_kwargs()
        kwargs["seqs_data"] = data
        kwargs.pop("slice_record", None)
        kwargs.pop("annotation_db", None)
        return self.__class__(**kwargs)

    def apply_scaled_gaps(
        self,
        other: SequenceCollection,
        aa_to_codon: bool = True,
    ) -> typing_extensions.Self:
        """applies gaps in self to ungapped sequences"""
        assert set(other.names) == set(self.names), "Names must match"
        if aa_to_codon and not all(
            (not self.moltype.is_nucleic, other.moltype.is_nucleic),
        ):
            msg = "aa_to_codon True requires nucleic moltype destination not {self.moltype.name!r}"
            raise ValueError(
                msg,
            )
        if not aa_to_codon and not all(
            (self.moltype.is_nucleic, not other.moltype.is_nucleic),
        ):
            msg = f"aa_to_codon False requires protein moltype destination not {other.moltype.name!r}"
            raise ValueError(
                msg,
            )

        assert aa_to_codon is not None
        gaps = {}
        seqs = {}
        scale = 3 if aa_to_codon else 1 / 3

        for seq in self.seqs:
            parent_name = self._name_map[seq.name]
            gaps[parent_name] = seq.map.scaled(scale).array
            ungapped = numpy.array(other.seqs[seq.name])
            seqs[parent_name] = ungapped
            seq_len = self._seqs_data.get_seq_length(parent_name)
            if len(ungapped) != int(seq_len * scale):
                msg = f"Destination sequence for {seq.name!r} != {scale:.2f} x {seq_len} sequence"
                raise ValueError(
                    msg,
                )

        seq_data = self._seqs_data.from_seqs_and_gaps(
            seqs=seqs,
            gaps=gaps,
            alphabet=other.moltype.most_degen_alphabet(),
        )
        init_args = self._get_init_kwargs()
        init_args.pop("slice_record")  # slice is realised
        init_args["moltype"] = other.moltype
        init_args["seqs_data"] = seq_data
        init_args["annotation_db"] = other.annotation_db
        return self.__class__(**init_args)

    def copy(self, copy_annotations: bool = False) -> typing_extensions.Self:
        """creates new instance, only mutable attributes are copied"""
        kwargs = self._get_init_kwargs()
        if copy_annotations:
            kwargs["annotation_db"] = copy.deepcopy(self._annotation_db)
        return self.__class__(**kwargs)

    def deepcopy(self, **kwargs) -> typing_extensions.Self:
        """returns deep copy of self

        Notes
        -----
        Reduced to sliced sequences in self, kwargs are ignored.
        Annotation db is not copied if the alignment has been sliced.
        """
        import copy

        kwargs = self._get_init_kwargs()
        kwargs.pop("seqs_data")
        kwargs["annotation_db"] = (
            None
            if len(self) != self._seqs_data.align_len
            else copy.deepcopy(self._annotation_db)
        )
        new_seqs_data = {
            n: self._seqs_data.get_gapped_seq_array(
                seqid=n,
                start=self._slice_record.plus_start,
                stop=self._slice_record.plus_stop,
                step=self._slice_record.plus_step,
            )
            for n in self._name_map.values()
        }
        new_seqs_data = self._seqs_data.from_seqs(
            data=new_seqs_data,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        kwargs["seqs_data"] = new_seqs_data
        kwargs["slice_record"] = (
            c3_sequence.SliceRecord(parent_len=new_seqs_data.align_len, step=-1)
            if self._slice_record.is_reversed
            else None
        )
        return self.__class__(**kwargs)

    def with_masked_annotations(
        self,
        biotypes: PySeqStr,
        mask_char: str = "?",
        shadow: bool = False,
        seqid: str | None = None,
    ) -> typing_extensions.Self:
        """returns an alignment with regions replaced by mask_char

        Parameters
        ----------
        biotypes
            annotation type(s)
        mask_char
            must be a character valid for the moltype. The default value is
            the most ambiguous character, eg. '?' for DNA
        shadow
            If True, masks everything but the biotypes
        seqid
            name of sequence to mask, defaults to all
        """
        # by doing numpy.array(seq), this method applies the slice_record
        # and modifies the underlying sequence data. We therefore split
        # gapped sequences into their gaps and seq components and discard
        # the current slice_record.
        gaps = {}
        ungapped = {}
        for seq in self.seqs:
            parent_name = self._name_map[seq.name]
            gaps[parent_name] = seq.map.array
            if seqid is None or seq.name == seqid:
                ungapped[parent_name] = numpy.array(
                    seq.seq.with_masked_annotations(biotypes, mask_char, shadow),
                )
            else:
                ungapped[parent_name] = numpy.array(seq.seq)

        seq_data = self._seqs_data.from_seqs_and_gaps(
            seqs=ungapped,
            gaps=gaps,
            alphabet=self.moltype.most_degen_alphabet(),
        )
        init = self._get_init_kwargs()
        init["seqs_data"] = seq_data
        # we have to drop the slice record since we have changed the data
        init["slice_record"] = None
        return self.__class__(**init)

    def to_rich_dict(self) -> dict[str, str | dict[str, str]]:
        """returns a json serialisable dict"""
        kwargs = self._get_init_kwargs()
        kwargs.pop("slice_record")  # slice is realised
        kwargs["moltype"] = self.moltype.label
        kwargs.pop("annotation_db", None)  # we dont serialise the annotation db
        kwargs.pop(
            "offset",
            None,
        )  # no need for offset since annotation db is not serialised
        kwargs.pop("seqs_data")

        seqs = {self._name_map[s.name]: str(s) for s in self.seqs}
        return {
            "init_args": kwargs,
            "type": get_object_provenance(self),
            "version": __version__,
            "seqs": seqs,
        }

    @classmethod
    def from_rich_dict(cls, data: dict[str, str | dict[str, str]]) -> Alignment:
        data["init_args"].pop("annotation_db", None)
        return make_aligned_seqs(data["seqs"], **data["init_args"])

    def is_ragged(self) -> bool:
        """by definition False for an Alignment"""
        return False

    def strand_symmetry(self, motif_length: int = 1) -> dict[str, TestResult]:
        """returns dict of strand symmetry test results per ungapped seq"""
        return {
            s.name: s.seq.strand_symmetry(motif_length=motif_length) for s in self.seqs
        }

    def duplicated_seqs(self) -> list[list[str]]:
        """returns the names of duplicated sequences

        Notes
        -----
        The gapped sequence is used.
        """
        if not len(self):
            # all have zero lengths
            return [] if self.num_seqs < 2 else [list(self.names)]

        return super().duplicated_seqs()


@register_deserialiser(
    get_object_provenance(Alignment),
    "cogent3.core.alignment.Alignment",
    "cogent3.core.c3_alignment.Alignment",
)
def deserialise_alignment(data: dict[str, str | dict[str, str]]) -> Alignment:
    if "init_args" not in data:
        # old type Alignment class
        return deserialise_alignment_to_new_type_alignment(data)
    return Alignment.from_rich_dict(data)


@register_deserialiser("cogent3.core.alignment.ArrayAlignment")
def deserialise_array_align_to_new_type_alignment(
    data: dict[str, str | dict[str, str]],
) -> Alignment:
    """deserialise old type ArrayAlignment as a new type alignment."""
    moltype_name = data["moltype"]
    moltype_name = "text" if moltype_name == "bytes" else moltype_name
    info_data = data["info"]
    source = info_data.pop("source", None) if info_data else None

    seq_data = {}
    for seqid, record in data["seqs"].items():
        if isinstance(record["seq"], str):
            seq = record["seq"]
        else:
            seq = record["seq"]["init_args"]["seq"]

        seq_data[seqid] = seq

    return make_aligned_seqs(
        seq_data, moltype=moltype_name, info=info_data, source=source
    )


def deserialise_alignment_to_new_type_alignment(
    data: dict[str, str | dict[str, str]],
) -> Alignment:
    """deserialise old type Alignment as a new type alignment"""
    from cogent3 import get_moltype
    from cogent3.core.location import deserialise_indelmap

    moltype_name = data["moltype"]
    mt = get_moltype(moltype_name)
    alpha = mt.most_degen_alphabet()
    info_data = data["info"]
    source = info_data.pop("source", None)

    seqs = {}
    gaps = {}
    for name, record in data["seqs"].items():
        seq_data = record["seq_init"]
        minit = record["map_init"]
        imap = deserialise_indelmap(minit)
        raw_seq = seq_data["seq"]["init_args"]["seq"]
        seqs[name] = raw_seq
        gaps[name] = imap.array

    asd = AlignedSeqsData.from_seqs_and_gaps(seqs=seqs, gaps=gaps, alphabet=alpha)
    return make_aligned_seqs(asd, moltype=moltype_name, info=info_data, source=source)


def make_aligned_storage(
    data: dict[str, bytes | numpy.ndarray[int]],
    *,
    moltype: str | c3_moltype.MolType,
    label_to_name: OptRenamerCallable = None,
    offset: DictStrInt | None = None,
    reversed_seqs: set[str] | None = None,
    storage_backend: str | None = None,
    **kwargs,
) -> AlignedSeqsDataABC:
    """makes the aligned storage instance for Alignment class

    Parameters
    ----------
    data
        {seq name: sequence, ...}
    moltype
        label for looking up the molecular type
    label_to_name
        renamer function
    offset
        mapping of {name: annotation offset}
    reversed_seqs
        name of sequences that are reverse complemented relative to their
        parent
    storage_backend
        name of a third-party storage driver to provide storage functionality
    kwargs
        additional keyword arguments for the storage driver

    Notes
    -----
    This function is intended for use primarly by make_aligned_seqs function.
    """
    assign_names = _SeqNamer(name_func=label_to_name)
    moltype = c3_moltype.get_moltype(moltype)
    alphabet = moltype.most_degen_alphabet()
    seqs_data, offs, rvd = prep_for_seqs_data(data, moltype, assign_names)
    asd_kwargs = {
        "alphabet": alphabet,
        "offset": offset or offs,
        "reversed_seqs": reversed_seqs or rvd,
        "data": seqs_data,
        **kwargs,
    }
    # plugin module is private only to exclude users, not developers
    klass = cogent3._plugin.get_aligned_storage_driver(storage_backend)  # noqa: SLF001
    return klass.from_seqs(**asd_kwargs)


@c3warn.deprecated_args(
    "2025.9", "no longer has an effect", discontinued=["new_type", "array_align"]
)
@singledispatch
def make_aligned_seqs(
    data: dict[str, StrORBytesORArray] | list | AlignedSeqsDataABC,
    *,
    moltype: str | c3_moltype.MolType,
    label_to_name: OptRenamerCallable = None,
    info: OptDict = None,
    source: OptPathType = None,
    annotation_db: SupportsFeatures | None = None,
    offset: DictStrInt | None = None,
    name_map: DictStrStr | None = None,
    is_reversed: OptBool = None,
    reversed_seqs: set[str] | None = None,
    storage_backend: str | None = None,
    **kwargs,
) -> Alignment:
    """Initialise an aligned collection of sequences.

    Parameters
    ----------
    data
        sequence data, a AlignedSeqsData, a dict {name: seq, ...}, an iterable of sequences
    moltype
        string representation of the moltype, e.g., 'dna', 'protein'.
    label_to_name
        function for converting original names into other names.
    info
        a dict from which to make an info object
    source
        origins of this data, defaults to 'unknown'. Converted to a string
        and added to info["source"].
    annotation_db
        annotation database to attach to the collection
    offset
        a dict mapping names to annotation offsets
    name_map
        a dict mapping sequence names to "parent" sequence names. The parent
        name will be used for querying a annotation_db.
    is_reversed
        entire collection has been reverse complemented
    reversed_seqs
        set of names that are on the reverse strand of the parent sequence
    storage_backend
        name of the storage backend to use for the SeqsData object, defaults to
        cogent3 builtin.
    kwargs
        keyword arguments for the AlignedSeqsData constructor

    Notes
    -----
    If no annotation_db is provided, but the sequences are annotated, an
    annotation_db is created by merging any annotation db's found in the sequences.
    If the sequences are annotated AND an annotation_db is provided, only the
    annotation_db is used.
    """
    if len(data) == 0:
        msg = "data must be at least one sequence."
        raise ValueError(msg)

    moltype = c3_moltype.get_moltype(moltype)
    annotation_db = kwargs.pop("annotation_db", annotation_db) or merged_db_collection(
        data,
    )

    # if we have Sequences, we need to construct the name map before we construct
    # the AlignedSeqsData object - however, if a name_map is provided, we assume that it
    # corrects for any naming differences in data and skip this step
    assign_names = _SeqNamer(name_func=label_to_name)
    data = _make_name_seq_mapping(data, assign_names)
    if name_map is None:
        name_map = make_name_map(data) or None

    seqs_data = make_aligned_storage(
        data,
        moltype=moltype,
        label_to_name=label_to_name,
        offset=offset,
        reversed_seqs=reversed_seqs,
        storage_backend=storage_backend,
        **kwargs,
    )
    # as they were handled in this function, we do not pass on:
    # offset
    # label_to_name
    # reversed_seqs
    # storage_backend
    return make_aligned_seqs(
        seqs_data,
        moltype=moltype,
        info=info,
        annotation_db=annotation_db,
        source=source,
        name_map=name_map,
        is_reversed=is_reversed,
    )


# I'm explicitly encoding the alternate variants of the following as numba
# can then cache the byte-compiled functions.
@make_aligned_seqs.register
def _(
    data: AlignedSeqsDataABC,
    *,
    moltype: str | c3_moltype.MolType,
    label_to_name: OptRenamerCallable = None,
    info: OptDict = None,
    source: OptPathType = None,
    annotation_db: SupportsFeatures | None = None,
    offset: DictStrInt | None = None,
    name_map: DictStrStr | None = None,
    is_reversed: OptBool = None,
    **kwargs,
) -> Alignment:
    moltype = c3_moltype.get_moltype(moltype)
    if not moltype.is_compatible_alphabet(data.alphabet):
        msg = (
            f"Provided moltype: {moltype.label} is not compatible with AlignedSeqsData"
        )
        raise ValueError(
            msg,
            f" alphabet: {data.alphabet}",
        )

    # we cannot set offset when creating from an AlignedSeqsData
    if offset:
        msg = f"Setting offset is not supported for {data=}"
        raise ValueError(msg)

    info = info if isinstance(info, dict) else {}
    source = str(source) if source else str(info.pop("source", "unknown"))
    sr = (
        c3_sequence.SliceRecord(parent_len=data.align_len, step=-1)
        if is_reversed
        else None
    )
    aln = Alignment(
        seqs_data=data,
        moltype=moltype,
        info=info,
        source=source,
        annotation_db=annotation_db,
        name_map=name_map,
        slice_record=sr,
    )
    if label_to_name:
        aln = aln.rename_seqs(label_to_name)
    return aln


@numba.jit(cache=True)
def _var_pos_canonical_or_gap(
    arr: numpy.ndarray,
    gap_index: int,
    missing_index: int,
) -> numpy.ndarray:  # pragma: no cover
    """return boolean array indicating columns with more than one value below threshold

    Parameters
    ----------
    arr
        a 2D array with rows being sequences and columns positions
    gap_index
        value of gap state
    missing_index
        value of missing state

    Returns
    ------
    a boolean array
    """
    # relying on consistent ordering of gap as num canonical + 1
    m, n = arr.shape
    result = numpy.zeros(n, dtype=numpy.bool_)
    for pos in numpy.arange(n):
        last = -1
        for seq in numpy.arange(m):
            state = arr[seq, pos]
            if state <= gap_index or state == missing_index:
                if last == -1:
                    last = state
                elif state != last:
                    result[pos] = True
                    break

    return result


@numba.jit(cache=True)
def _var_pos_canonical(
    arr: numpy.ndarray,
    gap_index: int,
) -> numpy.ndarray:  # pragma: no cover
    """return boolean array indicating columns with more than one value below threshold

    Parameters
    ----------
    arr
        a 2D array with rows being sequences and columns positions
    gap_index
        value of gap state

    Returns
    ------
    a boolean array
    """
    # relying on consistent ordering of gap as num canonical + 1
    m, n = arr.shape
    result = numpy.zeros(n, dtype=numpy.bool_)
    for pos in numpy.arange(n):
        last = -1
        for seq in numpy.arange(m):
            state = arr[seq, pos]
            if state < gap_index:
                if last == -1:
                    last = state
                elif state != last:
                    result[pos] = True
                    break

    return result


@numba.jit(cache=True)
def _var_pos_allow_all(
    arr: numpy.ndarray,
) -> numpy.ndarray:  # pragma: no cover
    """return boolean array indicating columns with more than one value below threshold

    Parameters
    ----------
    arr
        a 2D array with rows being sequences and columns positions

    Returns
    ------
    a boolean array
    """
    m, n = arr.shape
    result = numpy.zeros(n, dtype=numpy.bool_)
    for pos in numpy.arange(n):
        last = -1
        for seq in numpy.arange(m):
            state = arr[seq, pos]
            if last == -1:
                last = state
            elif state != last:
                result[pos] = True
                break
    return result


@numba.jit(cache=True)
def _var_pos_not_gap(
    arr: numpy.ndarray,
    gap_index: int,
) -> numpy.ndarray:  # pragma: no cover
    """return boolean array indicating columns with more than one value below threshold

    Parameters
    ----------
    arr
        a 2D array with rows being sequences and columns positions
    gap_index
        value of gap state

    Returns
    ------
    a boolean array
    """
    m, n = arr.shape
    result = numpy.zeros(n, dtype=numpy.bool_)
    for pos in numpy.arange(n):
        last = -1
        for seq in numpy.arange(m):
            state = arr[seq, pos]
            if state != gap_index:
                if last == -1:
                    last = state
                elif state != last:
                    result[pos] = True
                    break

    return result
