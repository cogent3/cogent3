from __future__ import annotations

import json
import os
import re
import typing
import warnings
from abc import ABC, abstractmethod
from collections import defaultdict
from functools import singledispatch, singledispatchmethod
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping, Optional, Union

import numpy

from cogent3 import get_app
from cogent3._version import __version__
from cogent3.core import new_alphabet, new_moltype, new_sequence
from cogent3.core.annotation import Feature
from cogent3.core.annotation_db import (
    BasicAnnotationDb,
    FeatureDataType,
    SupportsFeatures,
)
from cogent3.core.info import Info as InfoClass
from cogent3.core.location import (
    FeatureMap,
    IndelMap,
    _input_vals_neg_step,
    _input_vals_pos_step,
)
from cogent3.core.profile import PSSM, MotifCountsArray, MotifFreqsArray, load_pssm
from cogent3.format.alignment import save_to_filename
from cogent3.format.fasta import seqs_to_fasta
from cogent3.format.phylip import alignment_to_phylip
from cogent3.maths.stats.number import CategoryCounter
from cogent3.util import progress_display as UI
from cogent3.util import warning as c3warn
from cogent3.util.deserialise import deserialise_object, register_deserialiser
from cogent3.util.dict_array import DictArray, DictArrayTemplate
from cogent3.util.io import atomic_write, get_format_suffixes
from cogent3.util.misc import (
    get_object_provenance,
    get_setting_from_environ,
    negate_condition,
)
from cogent3.util.union_dict import UnionDict

DEFAULT_ANNOTATION_DB = BasicAnnotationDb

OptInt = Optional[int]
OptFloat = Optional[float]
OptStr = Optional[str]
OptList = Optional[list]
OptDict = Optional[dict]
OptCallable = Optional[Callable]
OptRenamerCallable = Optional[Callable[[str], str]]
OptPathType = Union[str, Path, None]
StrORArray = Union[str, numpy.ndarray]
StrORBytesORArray = Union[str, bytes, numpy.ndarray]
MolTypes = Union[str, new_moltype.MolType]


def assign_sequential_names(num_seqs: int, base_name: str = "seq", start_at: int = 0):
    """Returns list of sequential, unique names, e.g., ['seq_0' ... 'seq_n'] for
     num_seqs = n+1

    Parameters
    ----------
    num_seqs
        number of names to generate
    base_name
        the sequence name prefix
    start_at
        the number to start the sequence names at
    """
    return [f"{base_name}_{i}" for i in range(start_at, start_at + num_seqs)]


class SeqDataView(new_sequence.SeqView):
    """
    A view class for SeqsData, providing methods for different representations
    of a single sequence.

    self.seq is a SeqsData() instance, but other properties are a reference to a
    single seqid only.

    Example
    -------
    data = {"seq1": "ACGT", "seq2": "GTTTGCA"}
    sd = SeqsData(data=data)
    sdv = sd.get_seq_view(seqid="seq1")
    sdv.array_value
    # array([3, 1, 2, 0], dtype=int8)
    """

    __slots__ = ("parent", "alphabet", "_seqid", "_parent_len", "_slice_record")

    @property
    def str_value(self) -> str:
        """returns the sequence as a string"""
        # todo: kath, in ADV, the .get_seq_str method gets passed the step and
        # the returned sequence is sliced by the step. In SDV, the step is not
        # applied in the get_seq_str method, but is applied in this method here.
        # can we make this consistent?

        # also, keep using parent_start and parent_stop or to use start and stop?
        # this is another inconsistency between SDV and ADV
        raw = self.parent.get_seq_str(
            seqid=self.seqid,
            start=self.slice_record.parent_start,
            stop=self.slice_record.parent_stop,
        )
        return raw if self.slice_record.step == 1 else raw[:: self.slice_record.step]

    @property
    def array_value(self) -> numpy.ndarray:
        """returns the sequence as a numpy array"""
        raw = self.parent.get_seq_array(
            seqid=self.seqid,
            start=self.slice_record.parent_start,
            stop=self.slice_record.parent_stop,
        )
        return raw if self.slice_record.step == 1 else raw[:: self.slice_record.step]

    @property
    def bytes_value(self) -> bytes:
        """returns the sequence as bytes"""
        raw = self.parent.get_seq_bytes(
            seqid=self.seqid,
            start=self.slice_record.parent_start,
            stop=self.slice_record.parent_stop,
        )
        return raw if self.slice_record.step == 1 else raw[:: self.slice_record.step]

    def __repr__(self) -> str:
        seq = f"{self[:10]!s}...{self[-5:]}" if len(self) > 15 else str(self)
        return (
            f"{self.__class__.__name__}(seqid={self.seqid!r}, parent={seq}, "
            f"slice_record={self.slice_record!r})"
        )

    # refactor: design, do we support copy? do we support copy with sliced?
    def copy(self, sliced: bool = False):
        """returns copy"""
        return self

    def to_rich_dict(self) -> dict[str, str | dict[str, str]]:
        """returns a json serialisable dict.

        Notes
        -----
        This method will slice the underlying sequence to the start and stop values

        Warning
        -------
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
        new_sr = new_sequence.SliceRecord(
            parent_len=(stop - start),
            step=self.slice_record.step,
            offset=self.slice_record.parent_start,
        )
        data["init_args"]["slice_record"] = new_sr.to_rich_dict()
        data["init_args"]["alphabet"] = self.alphabet.to_rich_dict()
        return data


class SeqsDataABC(ABC):
    """Abstract base class for respresenting the collection of sequences underlying
    a SequenceCollection
    """

    __slots__ = ()

    # refactor: design
    # SeqsData needs a class method for creating instances from a dict, which will need to be
    # on AlignedSeqsData as well. If we rename AlignedSeqsData.from_aligned_seqs to from_seqs
    # (or something similar), we can use the same method name for both classes.

    @abstractmethod
    def seq_lengths(self) -> dict[str, int]: ...

    @property
    @abstractmethod
    def names(self) -> list: ...

    @property
    @abstractmethod
    def alphabet(self) -> new_alphabet.CharAlphabet: ...

    @property
    @abstractmethod
    def offset(self) -> dict[str, int]: ...

    @abstractmethod
    def get_seq_array(
        self, *, seqid: str, start: OptInt = None, stop: OptInt = None
    ) -> numpy.ndarray: ...

    @abstractmethod
    def get_seq_str(
        self, *, seqid: str, start: OptInt = None, stop: OptInt = None
    ) -> str: ...

    @abstractmethod
    def get_seq_bytes(
        self, *, seqid: str, start: OptInt = None, stop: OptInt = None
    ) -> bytes: ...

    @abstractmethod
    def get_view(self, seqid: str) -> new_sequence.SeqViewABC: ...

    @abstractmethod
    def subset(self, names: Union[str, typing.Sequence[str]]) -> SeqsDataABC: ...

    @abstractmethod
    def to_alphabet(self, alphabet: new_alphabet.AlphabetABC) -> SeqsDataABC: ...

    @abstractmethod
    def add_seqs(self, seqs) -> SeqsDataABC: ...

    @abstractmethod
    def to_rich_dict(self) -> dict: ...

    @abstractmethod
    def __len__(self) -> int: ...

    @abstractmethod
    def __getitem__(
        self, index: Union[str, int]
    ) -> Union[new_sequence.Sequence, new_sequence.SeqViewABC]: ...


class SeqsData(SeqsDataABC):
    """A collection of sequences underlying a SequenceCollection. The sequence
    data is stored as numpy arrays, however the underlying data can be accessed
    as strings, bytes, or numpy arrays.

    Attributes
    ----------
    data
        a dictionary of {name: sequence} pairs
    alphabet
        an instance of CharAlphabet valid for the sequences
    offset
        a dictionary of {name: offset} pairs indicating the offset of the sequence
    reversed
        a boolean indicating if the sequences are reversed
    check
        a boolean indicating if the data should be checked for naming consistency
        between arguments
    """

    __slots__ = ("_data", "_alphabet", "_offset", "_reversed")

    def __init__(
        self,
        *,
        data: dict[str, StrORBytesORArray],
        alphabet: new_alphabet.AlphabetABC,
        offset: dict[str, int] = None,
        reversed: bool = False,
        check: bool = True,
    ):
        self._alphabet = alphabet
        self._data: dict[str, numpy.ndarray] = {}
        for name, seq in data.items():
            arr = self._alphabet.to_indices(seq)
            arr.flags.writeable = False
            self._data[str(name)] = arr
        self._offset = offset or {}
        if check:
            assert set(self._offset.keys()) <= set(
                self._data.keys()
            ), "sequence name provided in offset not found in data"
        self._reversed = reversed
        # refactor: rename reversed argument to not clash with python built-in

    @property
    def names(self) -> list:
        return list(self._data.keys())

    @property
    def alphabet(self) -> new_alphabet.CharAlphabet:
        return self._alphabet

    @property
    def offset(self) -> dict[str, int]:
        return {name: self._offset.get(name, 0) for name in self.names}

    @property
    def is_reversed(self) -> bool:
        return self._reversed

    def reverse(self) -> SeqsData:
        """Reverse the orientation of all sequences in the collection."""
        # refactor: design, this can't be put on the SeqsDataABC unless we add it to the
        # AlignedSeqsData class as well -- however, ASD does not have the responsibility
        # of managing the orientation of the sequences, the lies with the
        # Alignment.slice_record

        # we dont touch the data, just toggle the reversed attribute
        return self.__class__(
            data=self._data,
            alphabet=self.alphabet,
            offset=self._offset,
            reversed=not self._reversed,
            check=False,
        )

    def seq_lengths(self) -> dict[str, int]:
        """Returns lengths of sequences as dict of {name: length, ... }."""
        return {name: seq.shape[0] for name, seq in self._data.items()}

    def get_seq_array(
        self, *, seqid: str, start: OptInt = None, stop: OptInt = None
    ) -> numpy.ndarray:
        return self._data[seqid][start:stop]

    def get_seq_str(
        self, *, seqid: str, start: OptInt = None, stop: OptInt = None
    ) -> str:
        return self._alphabet.from_indices(
            self.get_seq_array(seqid=seqid, start=start, stop=stop)
        )

    def get_seq_bytes(
        self, *, seqid: str, start: OptInt = None, stop: OptInt = None
    ) -> bytes:
        return self.get_seq_str(seqid=seqid, start=start, stop=stop).encode("utf8")

    def get_view(self, seqid: str) -> SeqDataView:
        seq_len = len(self._data[seqid])
        offset = self._offset.get(seqid, 0)
        slice_record = new_sequence.SliceRecord(
            step=-1 if self.is_reversed else 1, parent_len=seq_len, offset=offset
        )
        return SeqDataView(
            parent=self,
            seqid=seqid,
            parent_len=seq_len,
            alphabet=self.alphabet,
            slice_record=slice_record,
        )

    def subset(self, names: Union[str, typing.Sequence[str]]) -> SeqsData:
        """Returns a new SeqsData object with only the specified names."""
        names = [names] if isinstance(names, str) else names
        if data := {name: self._data.get(name) for name in names if name in self.names}:
            return self.__class__(
                data=data,
                alphabet=self.alphabet,
                offset={
                    name: offset
                    for name, offset in self._offset.items()
                    if name in names
                },
                reversed=self._reversed,
                check=False,
            )
        else:
            raise ValueError(f"provided {names=} not found in collection")

    def rename_seqs(self, renamer: Callable[[str], str]) -> SeqsData:
        # todo: kath, create a map between original seqid and current seqname for feature querying?
        renamed_data = {}
        renamed_offset = {}

        for name in self._data:
            new_name = renamer(name)
            renamed_data[new_name] = self._data[name]
            if name in self._offset:
                renamed_offset[new_name] = self._offset[name]

        return self.__class__(
            data=renamed_data,
            alphabet=self.alphabet,
            offset=renamed_offset,
            reversed=self._reversed,
            check=False,
        )

    def add_seqs(
        self,
        seqs: dict[str, StrORBytesORArray],
        force_unique_keys=True,
        offset=None,
    ) -> SeqsData:
        """Returns a new SeqsData object with added sequences. If force_unique_keys
        is True, raises ValueError if any names already exist in the collection."""
        # refactor: design
        # if the collection has been reversed, we assume the same for
        # the new data. If this is not a desirable behaviour, we could realise
        # the reversal of the current collection, add the new data, and set
        # reversed to False
        if force_unique_keys and any(name in self.names for name in seqs):
            raise ValueError("One or more sequence names already exist in collection")
        new_data = {
            **self._data,
            **{name: self.alphabet.to_indices(seq) for name, seq in seqs.items()},
        }
        return self.__class__(
            data=new_data,
            alphabet=self.alphabet,
            offset={**self._offset, **(offset or {})},
            reversed=self._reversed,
        )

    def to_alphabet(
        self, alphabet: new_alphabet.AlphabetABC, check_valid=True
    ) -> SeqsData:
        # refactor: design -- map directly between arrays?
        # refactor: design -- better way to check if alphabets are DNA and RNA?
        if len(alphabet) == len(self.alphabet):
            return self.__class__(data=self._data, alphabet=alphabet)

        new_data = {}
        old = self.alphabet.as_bytes()
        new = alphabet.as_bytes()
        convert_old_to_bytes = new_alphabet.array_to_bytes(old)
        convert_bytes_to_new = new_alphabet.bytes_to_array(
            new, dtype=new_alphabet.get_array_type(len(new))
        )

        for seqid in self.names:
            seq_data = self.get_seq_array(seqid=seqid)
            as_new_alpha = convert_bytes_to_new(convert_old_to_bytes(seq_data))

            if check_valid and not alphabet.is_valid(as_new_alpha):
                raise ValueError(
                    f"Changing from old alphabet={self.alphabet} to new "
                    f"{alphabet=} is not valid for this data"
                )
            new_data[seqid] = as_new_alpha

        return self.__class__(
            data=new_data,
            alphabet=alphabet,
            offset=self._offset,
            reversed=self._reversed,
            check=False,
        )

    def __len__(self):
        return len(self.names)

    @singledispatchmethod
    def __getitem__(self, index: Union[str, int]) -> new_sequence.SeqViewABC:
        raise NotImplementedError(f"__getitem__ not implemented for {type(index)}")

    @__getitem__.register
    def _(self, index: str) -> new_sequence.SeqViewABC:
        return self.get_view(seqid=index)

    @__getitem__.register
    def _(self, index: int) -> new_sequence.SeqViewABC:
        return self[self.names[index]]

    def to_rich_dict(self) -> dict[str, str | dict[str, str]]:
        """returns a json serialisable dict"""
        return {
            "init_args": {
                "data": {name: self.get_seq_str(seqid=name) for name in self.names},
                "alphabet": self.alphabet.to_rich_dict(),
                "offset": self._offset,
                "reversed": self.is_reversed,
            },
            "type": get_object_provenance(self),
            "version": __version__,
        }

    @classmethod
    def from_rich_dict(cls, data: dict[str, str | dict[str, str]]) -> SeqsData:
        """returns a new instance from a rich dict"""
        alphabet = deserialise_object(data["init_args"]["alphabet"])
        return cls(
            data=data["init_args"]["data"],
            alphabet=alphabet,
            offset=data["init_args"]["offset"],
            reversed=data["init_args"]["reversed"],
        )


@register_deserialiser(get_object_provenance(SeqsData))
def deserialise_seqs_data(data: dict[str, str | dict[str, str]]) -> SeqsData:
    return SeqsData.from_rich_dict(data)


class SequenceCollection:
    def __init__(
        self,
        *,
        seqs_data: SeqsDataABC,
        moltype: new_moltype.MolType,
        names: OptList = None,
        info: Optional[Union[dict, InfoClass]] = None,
        source: OptPathType = None,
        annotation_db: Optional[SupportsFeatures] = None,
    ):
        self._seqs_data = seqs_data
        self.moltype = moltype
        self.names = names
        if not isinstance(info, InfoClass):
            info = InfoClass(info) if info else InfoClass()
        self.info = info
        self.source = source
        self._repr_policy = dict(num_seqs=10, num_pos=60, ref_name="longest", wrap=60)
        self._annotation_db = annotation_db or DEFAULT_ANNOTATION_DB()
        self._seqs = None
        self._post_init()

    def _post_init(self):
        # override in subclasses
        self._seqs = _IndexableSeqs(self, make_seq=self._make_seq)

    def _make_seq(self, name: str) -> new_sequence.Sequence:
        sv = self._seqs_data.get_view(name)
        return self.moltype.make_seq(seq=sv, name=name)

    @property
    def seqs(self) -> _IndexableSeqs:
        return self._seqs

    @property
    def names(self) -> list:
        return self._names

    @names.setter
    def names(self, names: list[str]):
        # refactor: design, do we accept the names as an arg, or should we rely on seqs_data keys?
        # if we are allowing the set of names to be a subset of the seqs_data keys, we should
        # ensure that we are using the names and not the keys in the rest of the code!
        if names is None:
            self._names = self._seqs_data.names
        elif set(names) <= set(self._seqs_data.names):
            self._names = names
        else:
            left_diff = set(names) - set(self._seqs_data.names)
            raise ValueError(f"Provided names not found in collection: {left_diff}")

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

    @property
    @c3warn.deprecated_callable(
        version="2025.5", reason=".seqs can now be indexed by name", new=".seqs"
    )
    def named_seqs(self) -> SeqsDataABC:  # pragma: no cover
        return self.seqs

    def iter_seqs(
        self, seq_order: OptList = None
    ) -> Iterator[Union[new_sequence.Sequence, new_sequence.SeqViewABC]]:
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
        # refactor: design
        # consider gavins comment that this method could operate on only the names
        # list and not the seqs_data object

        if isinstance(names, str):
            names = [names]

        if negate:
            names = [name for name in self.names if name not in names]
        else:
            names = [name for name in names if name in self.names]

        if not names:
            raise ValueError(f"{names=} and {negate=} resulted in no names")

        seqs_data = self._seqs_data.subset(names)

        result = self.__class__(
            seqs_data=seqs_data,
            moltype=self.moltype,
            names=names,
            info=self.info,
            source=self.source,
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
        self, f: Callable[[new_sequence.Sequence], bool], negate: bool = False
    ):
        """Returns list of names of seqs where f(seq) is True."""
        get = self.seqs

        new_f = negate_condition(f) if negate else f

        return [name for name in self.names if new_f(get[name])]

    def take_seqs_if(
        self, f: Callable[[new_sequence.Sequence], bool], negate: bool = False, **kwargs
    ):
        """Returns new collection containing seqs where f(seq) is True.

        Parameters
        ----------
        f
            function that takes a sequence object and returns True or False
        negate
            select all sequences EXCEPT those where f(seq) is True

        Notes
        -----
        The seqs in the new collection are the same objects as the
        seqs in the old collection, not copies.
        """
        return self.take_seqs(self.get_seq_names_if(f, negate), **kwargs)

    def get_seq(
        self, seqname: str, copy_annotations: bool = False
    ) -> new_sequence.Sequence:
        """Return a sequence object for the specified seqname.

        Parameters
        ----------
        seqname
            name of the sequence to return
        copy_annotations
            if True, only annotations from the selected seq are copied to the
            annotation_db of the new collection, the same annotation db is used.

        """
        # refactor: design
        # This method provides a mechanism for binding the annotation db to the sequence instance,
        # which self.seqs[seq_name] does not. This is a difference to the original implementation,
        # so it needs more thought.
        seq = self.seqs[seqname]
        if copy_annotations:
            seq.annotation_db = type(self.annotation_db)()
            seq.annotation_db.update(annot_db=self.annotation_db, seqids=seqname)
        else:
            seq.annotation_db = self.annotation_db

        return seq

    def add_seqs(
        self, seqs: Union[dict[str, StrORBytesORArray], SeqsData, list], **kwargs
    ) -> SequenceCollection:
        """Returns new collection with additional sequences.

        Parameters
        ----------
        seqs
            sequences to add
        """
        data = coerce_to_seqs_data_dict(seqs, label_to_name=None)
        seqs_data = self._seqs_data.add_seqs(data, **kwargs)
        return self.__class__(
            seqs_data=seqs_data,
            moltype=self.moltype,
            info=self.info,
            annotation_db=self.annotation_db,
        )

    def rename_seqs(self, renamer: Callable[[str], str]):
        """Returns new collection with renamed sequences."""
        result = self.__class__(
            seqs_data=self._seqs_data.rename_seqs(renamer),
            moltype=self.moltype,
            info=self.info,
            source=self.source,
        )
        if self.annotation_db:
            # refactor: design
            # how to manage to keep track of aliases? possibly stored
            # in the annotation db
            name_map = {renamer(name): name for name in self.names}
            result.info.name_map = name_map
            result.annotation_db = self.annotation_db

        return result

    def to_dict(self, as_array: bool = False) -> dict[str, Union[str, numpy.ndarray]]:
        """Return a dictionary of sequences.

        Parameters
        ----------
        as_array
            if True, sequences are returned as numpy arrays, otherwise as strings
        """
        return {
            name: (numpy.array(self.seqs[name]) if as_array else str(self.seqs[name]))
            for name in self.names
        }

    def to_rich_dict(self) -> dict[str, str | dict[str, str]]:
        """returns a json serialisable dict

        Notes
        -----
        Serialisation will not include the annotation_db if present.
        """
        moltype = self.moltype.label
        info = {} if self.info is None else self.info

        if info.get("Refs", None) is not None and "Refs" in info:
            info.pop("Refs")

        info = info or None

        init_args = dict(
            moltype=moltype,
            names=self.names,
            info=info,
        )

        if hasattr(self, "annotation_db") and self.annotation_db:
            init_args["annotation_db"] = self.annotation_db.to_rich_dict()

        data = {
            "init_args": init_args,
            "type": get_object_provenance(self),
            "version": __version__,
        }
        data["seqs_data"] = self._seqs_data.to_rich_dict()

        return data

    @classmethod
    def from_rich_dict(
        cls, data: dict[str, str | dict[str, str]]
    ) -> SequenceCollection:
        """returns a new instance from a rich dict"""
        seqs_data = deserialise_object(data["seqs_data"])
        data["init_args"].pop("annotation_db", None)
        moltype = data["init_args"].pop("moltype")
        moltype = new_moltype.get_moltype(moltype)

        return cls(seqs_data=seqs_data, **data["init_args"], moltype=moltype)

    def to_json(self):
        """returns json formatted string"""
        return json.dumps(self.to_rich_dict())

    def degap(self) -> SequenceCollection:
        """Returns new collection in which sequences have no gaps.

        Notes
        -----
        The returned collection will not retain an annotation_db if present.
        """
        # todo: kath, constructing the SeqsData object here violates loose coupling
        # can we put this in the SeqsData class?
        seqs = {
            name: self.moltype.degap(self._seqs_data.get_seq_array(seqid=name))
            for name in self.names
        }
        seqs_data = self._seqs_data.__class__(
            data=coerce_to_seqs_data_dict(seqs),
            alphabet=self._seqs_data.alphabet,
            offset=self._seqs_data.offset,
            reversed=self._seqs_data.is_reversed,
        )
        return self.__class__(
            seqs_data=seqs_data,
            moltype=self.moltype,
            info=self.info,
            source=self.source,
        )

    def to_moltype(self, moltype: str) -> SequenceCollection:
        """returns copy of self with changed moltype

        Parameters
        ----------
        moltype
            name of the new moltype, e.g, 'dna', 'rna'.

        Notes
        -----
        Cannot convert from nucleic acids to proteins. Use get_translation() for that.

        """
        mtype = new_moltype.get_moltype(moltype)
        if mtype is self.moltype:
            return self  # nothing to be done

        alpha = mtype.most_degen_alphabet()
        try:
            new_seqs_data = self._seqs_data.to_alphabet(alpha)
        except ValueError as e:
            raise ValueError(
                f"Failed to convert moltype from {self.moltype.label} to {moltype}"
            ) from e

        return self.__class__(
            seqs_data=new_seqs_data,
            moltype=mtype,
            info=self.info,
            source=self.source,
            annotation_db=self.annotation_db,
        )

    def to_dna(self):
        """returns copy of self as a collection of DNA moltype seqs"""
        return self.to_moltype("dna")

    def to_rna(self):
        """returns copy of self as a collection of RNA moltype seqs"""
        return self.to_moltype("rna")

    def get_translation(
        self,
        gc: int = 1,
        incomplete_ok: bool = False,
        include_stop: bool = False,
        trim_stop: bool = True,
        **kwargs,
    ):
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

        if len(self.moltype.alphabet) != 4:
            raise new_alphabet.AlphabetError("Sequences must be a DNA/RNA")

        translated = {}
        # do the translation
        for seqname in self.names:
            seq = self.seqs[seqname]
            pep = seq.get_translation(
                gc,
                incomplete_ok=incomplete_ok,
                include_stop=include_stop,
                trim_stop=trim_stop,
            )
            translated[seqname] = pep

        seqs_data = coerce_to_seqs_data_dict(translated, label_to_name=None)
        pep_moltype = pep.moltype

        seqs_data = self._seqs_data.__class__(
            data=seqs_data,
            alphabet=pep_moltype.most_degen_alphabet(),
            offset=self._seqs_data.offset,
            reversed=self._seqs_data.is_reversed,
        )
        return self.__class__(
            seqs_data=seqs_data,
            moltype=pep_moltype,
            info=self.info,
            source=self.source,
            **kwargs,
        )

    def rc(self):
        """Returns the reverse complement of all sequences in the collection.
        A synonym for reverse_complement.

        Notes
        -----
        Reverse complementing the collection will break the relationship to an
        annotation_db if present.
        """
        return self.__class__(
            seqs_data=self._seqs_data.reverse(),
            names=self.names,
            info=self.info,
            moltype=self.moltype,
            source=self.source,
        )

    def reverse_complement(self):
        """Returns the reverse complement of all sequences in the collection.
        A synonym for rc.

        Notes
        -----
        Reverse complementing the collection will break the relationship to an
        annotation_db if present.

        """
        return self.rc()

    def distance_matrix(self, calc: str = "pdist"):
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
            raise NotImplementedError("only defined for DNA/RNA molecular types")

        # assert we have more than one sequence in the SequenceCollection
        if self.num_seqs == 1:
            raise ValueError(
                "Pairwise distance cannot be computed for a single sequence. "
                "Please provide at least two sequences."
            )

        dist_calc_app = get_approx_dist_calc(
            dist=calc, num_states=len(self.moltype.alphabet)
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

        # refactor
        # if we update a db with the same db, it just hangs... need to fix
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

        if not self.annotation_db:
            return None

        if seqid and (seqid not in self.names):
            raise ValueError(f"unknown {seqid=}")

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
        return seqs_to_fasta(self.to_dict(), block_size=block_size)

    def to_phylip(self):
        """
        Return collection in PHYLIP format and mapping to sequence ids

        Notes
        -----
        raises exception if sequences do not all have the same length
        """
        if self.is_ragged():
            raise ValueError("not all seqs same length, cannot convert to phylip")

        return alignment_to_phylip(self.to_dict())

    def write(self, filename: str, file_format: OptStr = None, **kwargs):
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

        save_to_filename(self.to_dict(), filename, file_format, **kwargs)

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

        Returns
        -------
        a Drawable or AnnotatedDrawable
        """
        from cogent3.draw.dotplot import Dotplot
        from cogent3.draw.drawable import AnnotatedDrawable

        if k is not None:
            assert 0 < k < window, "k must be smaller than window size"

        if len(self.names) == 1:
            name1 = name2 = self.names[0]
        elif name1 is None and name2 is None:
            name1, name2 = list(numpy.random.choice(self.names, size=2, replace=False))
        elif not (name1 and name2):
            names = list(set(self.names + [None]) ^ {name1, name2})
            name = list(numpy.random.choice(names, size=1))[0]
            name1 = name1 or name
            name2 = name2 or name

        if not {name1, name2} <= set(self.names):
            msg = f"{name1}, {name2} missing"
            raise ValueError(msg)

        seq1 = self.get_seq(seqname=name1, copy_annotations=False)
        seq2 = self.get_seq(seqname=name2, copy_annotations=False)

        if seq1.is_annotated() or seq2.is_annotated():
            annotated = True
            data = getattr(seq1, "data", seq1)
            bottom = data.get_drawable()
            data = getattr(seq2, "data", seq2)
            left = data.get_drawable(vertical=True)
        else:
            annotated = False

        dotplot = Dotplot(
            seq1,
            seq2,
            False,
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
        # refactor: performance
        result = {}
        for name in self.names:
            result[name] = ambig = {}
            for i, motif in enumerate(self.seqs[name]):
                if self.moltype.is_ambiguity(motif):
                    ambig[i] = motif
        return result

    def trim_stop_codons(self, gc: Any = 1, strict: bool = False):
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

        new_seqs = {s.name: s.trim_stop_codon(gc=gc, strict=strict) for s in self.seqs}
        # todo: kath
        # should we be using the seqs_data object to do this?
        # creating the new SeqsData object here violates the loose coupling
        # principle which we have been trying to adhere to.
        # or we warn of the potential change in the seqs_data object

        seqs_data = self._seqs_data.__class__(
            data=coerce_to_seqs_data_dict(new_seqs),
            alphabet=self._seqs_data.alphabet,
            offset=self._seqs_data.offset,
            reversed=self._seqs_data.is_reversed,
        )
        result = self.__class__(
            seqs_data=seqs_data,
            moltype=self.moltype,
            info=self.info,
            source=self.source,
        )
        if self.annotation_db:
            result.annotation_db = self.annotation_db
        return result

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
                exclude_unobserved=exclude_unobserved,
                warn=warn,
            )
            motifs.update(c.keys())
            counts.append(c)
        motifs = list(sorted(motifs))
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

    def get_motif_probs(
        self,
        alphabet: new_alphabet.AlphabetABC = None,
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

        counts = {}
        for seq_name in self.names:
            sequence = self.seqs[seq_name]
            motif_len = alphabet.motif_len
            if motif_len > 1:
                posns = list(range(0, len(sequence) + 1 - motif_len, motif_len))
                sequence = [sequence[i : i + motif_len] for i in posns]
            for motif in sequence:
                if not allow_gap and self.moltype.gap in motif:
                    continue

                if motif in counts:
                    counts[motif] += 1
                else:
                    counts[motif] = 1

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

    def pad_seqs(self, pad_length: OptInt = None):
        """Returns copy in which sequences are padded to same length.

        Parameters
        ----------
        pad_length
            Length all sequences are to be padded to. Will pad to max sequence
            length if pad_length is None or less than max length.
        """

        max_len = max(self._seqs_data.seq_lengths().values())

        if pad_length is None:
            pad_length = max_len
        elif pad_length < max_len:
            raise ValueError(
                f"pad_length must be at greater or equal to maximum sequence length: {str(max_len)}"
            )

        new_seqs = {}
        for seq_name in self.names:
            seq = self._seqs_data.get_seq_str(seqid=seq_name)
            padded_seq = seq + "-" * (pad_length - len(seq))
            new_seqs[seq_name] = padded_seq

        seqs_data = self._seqs_data.__class__(
            data=new_seqs,
            alphabet=self._seqs_data.alphabet,
            offset=self._seqs_data.offset,
            reversed=self._seqs_data.is_reversed,
        )
        return self.__class__(
            seqs_data=seqs_data,
            moltype=self.moltype,
            info=self.info,
            source=self.source,
            annotation_db=self.annotation_db,
        )

    def strand_symmetry(self, motif_length: int = 1):
        """returns dict of strand symmetry test results per seq"""
        return {s.name: s.strand_symmetry(motif_length=motif_length) for s in self.seqs}

    def is_ragged(self) -> bool:
        return len(set(self._seqs_data.seq_lengths().values())) > 1

    def has_terminal_stop(self, gc: Any = None, strict: bool = False) -> bool:
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
        self, mask_degen: bool = False
    ) -> list[set]:  # refactor: array/simplify
        """returns sets of names for sequences that are identical

        Parameters
        ----------
        mask_degen
            if True, degenerate characters are ignored

        """

        if self.is_ragged():
            raise ValueError("not all seqs same length, cannot get identical sets")

        if mask_degen and not self.moltype.degen_alphabet:
            warnings.warn(
                "in get_identical_sets, mask_degen has no effect as moltype "
                f"{self.moltype.label!r} has no degenerate characters",
                UserWarning,
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
            identical_sets = list(dupes.values())
            return identical_sets

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
        target: new_sequence.Sequence,
        min_similarity: float = 0.0,
        max_similarity: float = 1.0,
        metric: Callable[
            [new_sequence.Sequence, new_sequence.Sequence], float
        ] = new_sequence.frac_same,
        transform: bool = None,
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

        Warning
        -------
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

    def __str__(self):
        """Returns self in FASTA-format, respecting name order."""
        from cogent3.format.alignment import FORMATTERS

        return FORMATTERS["fasta"](self.to_dict())

    def __eq__(
        self, other: Union[SequenceCollection, dict]
    ) -> bool:  # refactor: design
        return id(self) == id(other)

    def __ne__(self, other: SequenceCollection) -> bool:
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
        limit: OptInt = None,
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
            colors=colors, font_size=font_size, font_family=font_family
        )

        seq_lengths = numpy.array(list(self._seqs_data.seq_lengths().values()))
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
        num_seqs: OptInt = None,
        num_pos: OptInt = None,
        ref_name: OptInt = None,
        wrap: OptInt = None,
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


@register_deserialiser(get_object_provenance(SequenceCollection))
def deserialise_sequence_collection(data) -> SequenceCollection:
    return SequenceCollection.from_rich_dict(data)


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
        if not isinstance(seq, new_sequence.Sequence):
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
    data, label_to_name: OptRenamerCallable = None
) -> dict[str, StrORBytesORArray]:
    # refactor: handle conversion of SeqView to SeqDataView.
    raise NotImplementedError(
        f"coerce_to_seqs_data_dict not implemented for {type(data)}"
    )


@coerce_to_seqs_data_dict.register
def _(
    data: dict, label_to_name: OptRenamerCallable = None
) -> dict[str, StrORBytesORArray]:
    is_sequence = isinstance(next(iter(data.values()), None), new_sequence.Sequence)
    return {
        (label_to_name(k) if label_to_name else k): (
            numpy.array(v) if is_sequence else v
        )
        for k, v in data.items()
    }


@coerce_to_seqs_data_dict.register
def _(
    data: list, label_to_name: OptRenamerCallable = None
) -> dict[str, StrORBytesORArray]:
    first = data[0]
    labelled_seqs = assign_names(first, data=data)
    return coerce_to_seqs_data_dict(labelled_seqs, label_to_name=label_to_name)


@coerce_to_seqs_data_dict.register
def _(
    data: set, label_to_name: OptRenamerCallable = None
) -> dict[str, StrORBytesORArray]:
    first = next(iter(data))
    labelled_seqs = assign_names(first, data=data)
    return coerce_to_seqs_data_dict(labelled_seqs, label_to_name=label_to_name)


@coerce_to_seqs_data_dict.register
def _(
    data: SequenceCollection, label_to_name: OptRenamerCallable = None
) -> dict[str, StrORBytesORArray]:
    return coerce_to_seqs_data_dict(
        data.to_dict(as_array=True), label_to_name=label_to_name
    )


@singledispatch
def assign_names(
    first, data: Union[list, set]
) -> dict[str, Union[StrORBytesORArray, new_sequence.Sequence]]:
    if isinstance(first, (str, bytes, numpy.ndarray)):
        names = assign_sequential_names(len(data))
        return dict(zip(names, data))
    raise NotImplementedError(f"assign_names not implemented for {type(first)}")


@assign_names.register
def _(
    first: new_sequence.Sequence, data: Union[list, set]
) -> dict[str, new_sequence.Sequence]:
    return {seq.name: seq for seq in data}


@assign_names.register
def _(first: list, data: Union[list, set]) -> dict[str, StrORBytesORArray]:
    return {name_seq[0]: name_seq[1] for name_seq in data}


@assign_names.register
def _(first: tuple, data: Union[list, set]) -> dict[str, StrORBytesORArray]:
    return {name_seq[0]: name_seq[1] for name_seq in data}


@singledispatch
def make_unaligned_seqs(
    data: Union[dict[str, StrORBytesORArray], SeqsData, list],
    *,
    moltype: Union[str, new_moltype.MolType],
    label_to_name: OptRenamerCallable = None,
    info: dict = None,
    source: OptPathType = None,
    annotation_db: SupportsFeatures = None,
    offset: dict[str, int] = None,
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

    Notes
    -----
    - If no annotation_db is provided, but the sequences are annotated, an
    annotation_db is created by merging any annotation db's found in the sequences.
    - If the sequences are annotated AND an annotation_db is provided, only the
    annotation_db is used.

    """
    # refactor: design/simplify
    # the dispatches above handle the different type of data that the previous
    # make_unaligned_seqs handled. This includes dicts, lists (where items are
    # either pairs of [name, seq], or just seqs), tuples, Sequences, etc.
    # We could simplify it greatly by only supporting dicts, SeqsDataABC,
    # or SequenceCollections.

    # refactor: design
    # rename offset to offsets as it could track potentially multiple offsets
    # refactor: add offset/s to docstring

    if len(data) == 0:
        raise ValueError("data must be at least one sequence.")

    annotation_db = annotation_db or merged_db_collection(data)

    seqs_data = coerce_to_seqs_data_dict(data, label_to_name=label_to_name)

    moltype = new_moltype.get_moltype(moltype)
    alphabet = moltype.most_degen_alphabet()

    seqs_data = SeqsData(data=seqs_data, alphabet=alphabet, offset=offset)
    return make_unaligned_seqs(
        seqs_data,
        moltype=moltype,
        label_to_name=label_to_name,
        info=info,
        source=source,
        annotation_db=annotation_db,
        offset=offset,
    )


@make_unaligned_seqs.register
def _(
    data: SeqsDataABC,
    *,
    moltype: Union[str, new_moltype.MolType],
    label_to_name: OptRenamerCallable = None,
    info: dict = None,
    source: OptPathType = None,
    annotation_db: SupportsFeatures = None,
    offset: dict[str, int] = None,
) -> SequenceCollection:
    moltype = new_moltype.get_moltype(moltype)
    if not moltype.is_compatible_alphabet(data.alphabet):
        raise ValueError(
            f"Provided moltype: {moltype} is not compatible with SeqsData alphabet {data.alphabet}"
        )

    info = info if isinstance(info, dict) else {}
    if source:
        info["source"] = str(source)
    else:
        info["source"] = str(info.get("source", "unknown"))

    return SequenceCollection(
        seqs_data=data,
        moltype=moltype,
        info=info,
        annotation_db=annotation_db,
        source=source,
    )


@singledispatch
def seq_to_gap_coords(
    seq: StrORBytesORArray,
    *,
    alphabet: new_alphabet.AlphabetABC,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Takes a sequence with (or without) gaps and returns an ungapped sequence
    and a map of the position and length of gaps in the original parent sequence
    """
    raise NotImplementedError(f"seq_to_gap_coords not implemented for type {type(seq)}")


@seq_to_gap_coords.register
def _(
    seq: numpy.ndarray,
    *,
    alphabet: new_alphabet.AlphabetABC,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    gaps_bool = seq == alphabet.gap_index
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
            # end of sequence
            parent_coords.append([start, i + 1])

    # format for IndelMap
    parent_coords = numpy.array(parent_coords)
    gap_lengths = numpy.array([end - start for start, end in parent_coords])
    cum_lengths = gap_lengths.cumsum()
    # get gap start positions in sequence coords
    gap_pos = parent_coords.T[0] - numpy.append(0, cum_lengths[:-1, numpy.newaxis])

    return ungapped, numpy.array([gap_pos, cum_lengths]).T


@seq_to_gap_coords.register
def _(seq: str, *, alphabet: new_alphabet.AlphabetABC) -> tuple[str, numpy.ndarray]:
    if not alphabet.is_valid(seq):
        raise ValueError(f"Sequence is invalid for alphabet {alphabet}")

    return seq_to_gap_coords(alphabet.to_indices(seq), alphabet=alphabet)


@seq_to_gap_coords.register
def _(
    seq: bytes, *, alphabet: new_alphabet.AlphabetABC
) -> tuple[numpy.ndarray, numpy.ndarray]:
    return seq_to_gap_coords(seq.decode("utf-8"), alphabet=alphabet)


class Aligned:
    """A single sequence in an alignment."""

    def __init__(self, data: AlignedDataView, moltype: new_moltype.MolType):
        self._data = data
        self._moltype = moltype

    def __len__(self) -> int:
        return len(self.map)

    @property
    def data(self) -> AlignedDataView:
        return self._data

    @property
    def map(self) -> IndelMap:
        return self.data.map

    @property
    def seq(self) -> new_sequence.Sequence:
        """Returns Sequence object, excluding gaps."""
        # if the slice record has abs(step) > 1, we cannot retain a connection to the underlying aligned
        # seq data container because the gaps are not going to be modulo the step.
        rev = False
        if abs(self.data.slice_record.step) == 1:
            seq = self.data.get_seq_view()
        elif self.data.slice_record.step < -1:
            # refactor: design
            # gapped_str_value will apply the step to the underlying data, so
            # we need to re-reverse the underlying data AND reverse the Sequence
            # so that the seq knows to complement the data on output
            # we should revisit this design
            seq = self.moltype.degap(self.data.gapped_str_value)[::-1]
            rev = True
        else:
            seq = self.moltype.degap(self.data.gapped_str_value)
        return (
            self.moltype.make_seq(seq=seq, name=self.data.seqid)[::-1]
            if rev
            else self.moltype.make_seq(seq=seq, name=self.data.seqid)
        )

    @property
    def gapped_seq(self) -> new_sequence.Sequence:
        """Returns Sequence object, including gaps."""
        seq = self.data.gapped_str_value
        if self.data.slice_record.step < 0:
            seq = self.moltype.complement(seq)
        return self.moltype.make_seq(seq=seq, name=self.data.seqid)

    @property
    def moltype(self) -> new_moltype.MolType:
        return self._moltype

    @property
    def name(self) -> str:
        return self.data.seqid

    def gap_vector(self) -> list[bool]:
        """Returns gap_vector of GappedSeq, for omit_gap_pos."""
        return self.gapped_seq.gap_vector()

    def __str__(self) -> str:
        return str(self.gapped_seq)

    def __array__(self, dtype=None, copy=None) -> numpy.ndarray:
        return numpy.array(self.gapped_seq, dtype=dtype)

    def __bytes__(self) -> bytes:
        return bytes(self.gapped_seq)

    def __iter__(self):
        """Iterates over sequence one motif (e.g. char) at a time, incl. gaps"""
        yield from self.gapped_seq

    @singledispatchmethod
    def __getitem__(self, span: Union[int, slice]):
        raise NotImplementedError(f"__getitem__ not implemented for {type(span)}")

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
        elif not span.useful:
            im = self.map[start:end]
            data = self.data.array_value[:0]
        else:
            # multiple spans
            align_coords = span.get_coordinates()
            im = self.map.joined_segments(align_coords)
            seq_map = self.map.make_seq_feature_map(span)
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

    def __repr__(self) -> str:
        # refactor: design
        # todo: when design is finalised, add tests for this
        return f"Aligned(map={self.data.map}, data={self.seq})"


class AlignedSeqsDataABC(SeqsDataABC):
    # all methods that are from SeqsDataABC should work in sequence coordinates
    # all methods unique to AlignedSeqsDataABC should work in aligned coordinates
    __slots__ = ()

    @classmethod
    @abstractmethod
    def from_aligned_seqs(
        cls,
        *,
        data: dict[str, StrORArray],
        alphabet: new_alphabet.AlphabetABC,
    ): ...

    @classmethod
    @abstractmethod
    def from_seqs_and_maps(
        cls,
        *,
        seqs: dict[str, StrORBytesORArray],
        gaps: dict[str, numpy.ndarray],
        alphabet: new_alphabet.AlphabetABC,
    ): ...

    @property
    @abstractmethod
    def align_len(self) -> int: ...

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


def _gapped_seq_len(seq: numpy.ndarray, gap_map: numpy.ndarray) -> int:
    """calculate the gapped sequence length from a ungapped sequence and gap map

    Parameters
    ----------
    seq
        numpy array of sequence indices
    gap_map
        numpy array of [gap index, cumulative gap length] pairs
    """
    try:
        gap_len = gap_map[-1][1]
    except IndexError:  # no gaps
        return len(seq)

    return len(seq) + gap_len


class AlignedSeqsData(AlignedSeqsDataABC):
    # refactor: docstring
    # refactor: design

    __slots__ = (
        "_seqs",
        "_gaps",
        "_alphabet",
        "_align_len",
        "_offset",
    )

    def __init__(
        self,
        *,
        seqs: Optional[dict[str, StrORBytesORArray]],
        gaps: Optional[dict[str, numpy.ndarray]],
        alphabet: new_alphabet.AlphabetABC,
        offset: dict[str, int] = None,
        align_len: OptInt = None,
        check: bool = True,
    ):
        self._alphabet = alphabet
        if check or not align_len:
            if not seqs or not gaps:
                raise ValueError("Both seqs and gaps must be provided.")
            if set(seqs.keys()) != set(gaps.keys()):
                raise ValueError("Keys in seqs and gaps must be identical.")
            if offset:
                assert set(offset.keys()).issubset(seqs.keys())
            gapped_seq_lengths = {
                _gapped_seq_len(v, l) for v, l in zip(seqs.values(), gaps.values())
            }
            if len(gapped_seq_lengths) != 1:
                raise ValueError("All sequence lengths must be the same.")
        self._seqs = {}
        for k, v in seqs.items():
            seq = self._alphabet.to_indices(v)
            self._seqs[k] = seq
            seq.flags.writeable = False
        self._gaps = {}
        for k, v in gaps.items():
            self._gaps[k] = v
            v.flags.writeable = False
        self._align_len = align_len or gapped_seq_lengths.pop()
        self._offset = offset or {}

    @classmethod
    def from_aligned_seqs(
        cls,
        *,
        data: dict[str, StrORArray],
        alphabet: new_alphabet.AlphabetABC,
        **kwargs,
    ):
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
            raise ValueError("All sequence lengths must be the same.")

        align_len = seq_lengths.pop()
        seqs = {}
        gaps = {}
        for name, seq in data.items():
            seq, gap_map = seq_to_gap_coords(seq, alphabet=alphabet)
            seq = alphabet.to_indices(seq)
            seq.flags.writeable = False
            gap_map.flags.writeable = False
            seqs[name], gaps[name] = seq, gap_map

        return cls(
            seqs=seqs,
            gaps=gaps,
            alphabet=alphabet,
            align_len=align_len,
            check=False,
            **kwargs,
        )

    @classmethod
    def from_seqs_and_maps(
        cls,
        *,
        seqs: dict[str, StrORBytesORArray],
        gaps: dict[str, numpy.ndarray],
        alphabet: new_alphabet.AlphabetABC,
        **kwargs,
    ):
        """Construct an AlignedSeqsData object from a dict of ungapped sequences
        and a corresponding dict of gap maps.

        Parameters
        ----------
        seqs
            dict of ungapped sequences {name: seq, ...}
        gap_maps
            dict of gap maps {name: gap_map, ...} where a gap_map is a 2D numpy array
            of [gap index, cumulative gap length] pairs
        alphabet
            alphabet object for the sequences
        """
        return cls(
            seqs=seqs,
            gaps=gaps,
            alphabet=alphabet,
            **kwargs,
        )

    @property
    def names(self) -> tuple[str]:
        return list(self._seqs.keys())

    @property
    def alphabet(self) -> new_alphabet.CharAlphabet:
        return self._alphabet

    @property
    def align_len(self) -> int:
        return self._align_len

    @property
    def offset(self) -> dict[str, int]:
        """returns the offset of each sequence in the Alignment"""
        return {name: self._offset.get(name, 0) for name in self.names}

    def __len__(self) -> int:
        return self.align_len

    @singledispatchmethod
    def __getitem__(self, index: Union[str, int]):
        raise NotImplementedError(f"__getitem__ not implemented for {type(index)}")

    @__getitem__.register
    def _(self, index: str):
        return self.get_view(index)

    @__getitem__.register
    def _(self, index: int):
        return self[self.names[index]]

    def seq_lengths(self) -> dict[str, int]:
        """Returns lengths of ungapped sequences as dict of {name: length, ... }."""
        # refactor: design
        # change to method get_seq_length(name) so API is not too closely tied to implementation
        return {name: len(seq) for name, seq in self._seqs.items()}

    @singledispatchmethod
    def get_view(self, seqid: str, slice_record: new_sequence.SliceRecord = None):
        return AlignedDataView(
            parent=self,
            seqid=seqid,
            alphabet=self.alphabet,
            slice_record=slice_record,
        )

    @get_view.register
    def _(self, seqid: int, slice_record: new_sequence.SliceRecord = None):
        return self.get_view(self.names[seqid], slice_record)

    def get_gaps(self, seqid: str) -> numpy.ndarray:
        return self._gaps[seqid]

    def get_seq_array(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> numpy.ndarray:
        """Return ungapped sequence corresponding to seqid as an array of indices.
        assumes start/stop are in sequence coordinates. Excludes gaps.
        """
        # refactor: design
        seq = self._seqs[seqid]
        return seq[start:stop:step]

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

        # refactor: design
        # should this pass the data sliced data to the ADV which then applies the step???

        step = 1 if step is None else step
        func = _input_vals_pos_step if step > 0 else _input_vals_neg_step
        start, stop, step = func(self.align_len, start, stop, step)
        seq = self._seqs[seqid]

        # if there's no gaps, just slice the sequence
        if len(seq) == self.align_len:
            return seq[start:stop:step]

        unknown = self.alphabet.gap_index

        # there is only gaps
        if len(seq) == 0:
            return numpy.full(
                int(-((stop - start) // -step)), unknown, dtype=numpy.uint8
            )

        gaps = self._gaps[seqid]
        indel_map = IndelMap(
            gap_pos=gaps[:, 0], cum_gap_lengths=gaps[:, 1], parent_length=len(seq)
        )

        if step < 0:
            seq = seq[::-1]
            indel_map = indel_map[::-1]
            start = -start - 1
            stop = -stop - 1
            step = -step

        # convert from alignment coordinates to sequences coordinates
        # selecting only the indices that are in the frame of the stride
        indices = []
        for span in indel_map.nongap():
            indices += [
                indel_map.get_seq_index(x)
                for x in range(span.start, span.end)
                if start <= x < stop and (x - start) % step == 0
            ]

        seq_sliced = seq[indices]
        map_sliced = indel_map[start:stop:step]
        concat = numpy.array([], dtype=numpy.uint8)

        # iterate through spans
        for span in map_sliced.spans:
            if span.lost:
                seg = numpy.full(span.length, unknown, dtype=numpy.uint8)
            else:
                seg = seq_sliced[span.start : span.end]
            concat = numpy.concatenate((concat, seg))

        return concat

    def get_seq_str(
        self,
        *,
        seqid: str,
        start: OptInt = None,
        stop: OptInt = None,
        step: OptInt = None,
    ) -> str:
        """Return ungapped sequence corresponding to seqid as a string.
        start/stop are in alignment coordinates and the alignment coordinates
        are converted to sequence coordinates. Excludes gaps."""
        return self.alphabet.from_indices(
            self.get_seq_array(seqid=seqid, start=start, stop=stop, step=step)
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
            self.get_gapped_seq_array(seqid=seqid, start=start, stop=stop, step=step)
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
        start/stop are in alignment coordinates and the alignment coordinates
        are converted to sequence coordinates. Excludes gaps."""
        return self.get_seq_str(seqid=seqid, start=start, stop=stop, step=step).encode(
            "utf8"
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
            seqid=seqid, start=start, stop=stop, step=step
        ).encode("utf8")

    def add_seqs(
        self,
        seqs: dict[str, StrORArray],
        force_unique_keys=True,
        offset: dict[str, int] = None,
    ) -> AlignedSeqsData:
        """Returns a new AlignedSeqsData object with added sequences.

        Parameters
        ----------
        seqs
            dict of sequences to add {name: seq, ...}
        force_unique_keys
            if True, raises ValueError if any sequence names already exist in the collection
        """
        if force_unique_keys and any(name in self.names for name in seqs):
            raise ValueError("One or more sequence names already exist in collection")

        new_seq_lens = {len(seq) for seq in seqs.values()}
        if len(new_seq_lens) != 1 or new_seq_lens.pop() != self.align_len:
            raise ValueError(
                "All sequences must be the same length as existing sequences"
            )

        new_seqs = {}
        new_gaps = {}
        for name, seq in seqs.items():
            seq, gap_map = seq_to_gap_coords(seq, alphabet=self.alphabet)
            seq = self.alphabet.to_indices(seq)
            seq.flags.writeable = False
            gap_map.flags.writeable = False
            new_seqs[name], new_gaps[name] = seq, gap_map

        return self.__class__(
            seqs={**self._seqs, **new_seqs},
            gaps={**self._gaps, **new_gaps},
            alphabet=self.alphabet,
            offset={**self._offset, **(offset or {})},
            align_len=self.align_len,
        )

    def subset(self, names: Union[str, typing.Sequence[str]]):
        """Returns a new AlignedSeqsData object with only the specified names."""
        names = [names] if isinstance(names, str) else names
        if seq_data := {
            name: self._seqs.get(name) for name in names if name in self.names
        }:
            gap_data = {name: self._gaps.get(name) for name in seq_data}
            return self.__class__(
                seqs=seq_data,
                gaps=gap_data,
                alphabet=self.alphabet,
                offset={
                    name: offset
                    for name, offset in self._offset.items()
                    if name in names
                },
                align_len=self.align_len,
                check=False,
            )
        else:
            raise ValueError(f"provided {names=} not found in collection")

    def to_alphabet(self, alphabet: new_alphabet.AlphabetABC):
        # todo: kath
        ...

    def to_rich_dict(self):
        # todo: kath
        ...


class AlignedDataViewABC(new_sequence.SeqViewABC):
    __slots__ = ()

    @abstractmethod
    def get_seq_view(self) -> new_sequence.SeqViewABC: ...

    @property
    @abstractmethod
    def gapped_str_value(self) -> str: ...

    @property
    @abstractmethod
    def gapped_array_value(self) -> numpy.ndarray: ...

    @property
    @abstractmethod
    def gapped_bytes_value(self) -> bytes: ...


class AlignedDataView(new_sequence.SeqViewABC):
    # refactor: docstring

    __slots__ = (
        "parent",
        "_seqid",
        "alphabet",
        "_offset",
        "_parent_len",
        "_slice_record",
    )

    def __init__(
        self,
        *,
        parent: AlignedSeqsDataABC,
        seqid: str,
        alphabet: new_alphabet.AlphabetABC,
        slice_record: new_sequence.SliceRecord = None,
    ):
        self.parent = parent
        self._seqid = seqid
        self.alphabet = alphabet
        self._parent_len = parent.align_len
        self._slice_record = (
            slice_record
            if slice_record is not None
            else new_sequence.SliceRecord(parent_len=self._parent_len)
        )

    @property
    def slice_record(self):
        return self._slice_record

    @slice_record.setter
    def slice_record(self, value: new_sequence.SliceRecordABC):
        self._slice_record = value

    @property
    def seqid(self) -> str:
        return self._seqid

    @property
    def parent_len(self) -> int:
        return self._parent_len

    @property
    def map(self) -> IndelMap:
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
            parent_length=self.parent.seq_lengths()[self.seqid],
        )

    @property
    def str_value(self) -> str:
        return self.parent.get_seq_str(
            seqid=self.seqid,
            start=self.slice_record.start,
            stop=self.slice_record.stop,
            step=self.slice_record.step,
        )

    @property
    def gapped_str_value(self) -> str:
        return self.parent.get_gapped_seq_str(
            seqid=self.seqid,
            start=self.slice_record.start,
            stop=self.slice_record.stop,
            step=self.slice_record.step,
        )

    @property
    def array_value(self) -> numpy.ndarray:
        return self.parent.get_seq_array(
            seqid=self.seqid,
            start=self.slice_record.start,
            stop=self.slice_record.stop,
            step=self.slice_record.step,
        )

    @property
    def gapped_array_value(self) -> numpy.ndarray:
        return self.parent.get_gapped_seq_array(
            seqid=self.seqid,
            start=self.slice_record.start,
            stop=self.slice_record.stop,
            step=self.slice_record.step,
        )

    @property
    def bytes_value(self) -> bytes:
        return self.parent.get_seq_bytes(
            seqid=self.seqid,
            start=self.slice_record.start,
            stop=self.slice_record.stop,
            step=self.slice_record.step,
        )

    @property
    def gapped_bytes_value(self) -> bytes:
        return self.parent.get_gapped_seq_bytes(
            seqid=self.seqid,
            start=self.slice_record.start,
            stop=self.slice_record.stop,
            step=self.slice_record.step,
        )

    def __str__(self) -> str:
        return self.str_value

    def __array__(self, dtype=None, copy=None) -> numpy.ndarray:
        arr = self.array_value
        if dtype:
            arr = arr.astype(dtype)
        return arr

    def __bytes__(self) -> bytes:
        return self.bytes_value

    def __getitem__(self, segment) -> AlignedDataViewABC:
        return self.__class__(
            parent=self.parent,
            seqid=self.seqid,
            alphabet=self.alphabet,
            slice_record=self.slice_record[segment],
        )

    def __repr__(self) -> str:
        # refactor: design
        # todo: once design is settled, need to be tested
        seq_preview = (
            f"{self.parent.get_seq_array(seqid=self.seqid, start=0, stop=10)}..."
            f"{self.parent.get_seq_array(seqid=self.seqid, start=self.parent_len-5)}"
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
        start = self.map.get_seq_index(self.slice_record.parent_start)
        stop = self.map.get_seq_index(self.slice_record.parent_stop)
        strand = -1 if self.slice_record.step < 0 else 1

        return self.seqid, start, stop, strand

    def copy(self, sliced: bool = False):
        return self

    def _get_init_kwargs(self) -> dict:
        return {
            "parent": self.parent,
            "seqid": self.seqid,
            "alphabet": self.alphabet,
            "slice_record": self.slice_record,
        }

    def to_rich_dict(self) -> dict:
        ...
        # todo: kath...

    def get_seq_view(self) -> new_sequence.SeqViewABC:
        # we want the parent coordinates in sequence coordinates
        # parent_seq_coords does not account for the stride
        seqid, start, stop, _ = self.parent_seq_coords()
        parent_len = self.parent.seq_lengths()[seqid]
        offset = self.parent.offset.get(seqid)
        sr = new_sequence.SliceRecord(
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


def make_gap_filter(template, gap_fraction, gap_run):
    """Returns f(seq) -> True if no gap runs and acceptable gap fraction.

    Calculations relative to template.
    gap_run = number of consecutive gaps allowed in either the template or seq
    gap_fraction = fraction of positions that either have a gap in the template
        but not in the seq or in the seq but not in the template
    NOTE: template and seq must both be ArraySequence objects.
    """
    template_gaps = numpy.array(template.gap_vector())

    def result(seq):
        """Returns True if seq adhers to the gap threshold and gap fraction."""
        seq_gaps = numpy.array(seq.gap_vector())
        # check if gap amount bad
        if sum(seq_gaps != template_gaps) / float(len(seq)) > gap_fraction:
            return False
        # check if gap runs bad
        if (
            b"\x01" * gap_run
            in numpy.logical_and(seq_gaps, numpy.logical_not(template_gaps))
            .astype(numpy.uint8)
            .tobytes()
        ):
            return False
        # check if insertion runs bad
        elif (
            b"\x01" * gap_run
            in numpy.logical_and(template_gaps, numpy.logical_not(seq_gaps))
            .astype(numpy.uint8)
            .tobytes()
        ):
            return False
        return True

    return result


class _IndexableSeqs:
    """container that is created by SequenceCollection and Alignment instances"""

    def __init__(
        self,
        parent: Union[SequenceCollection, Alignment],
        make_seq: typing.Callable[[str], typing.Union[new_sequence.Sequence, Aligned]],
    ):
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
        self, key: Union[str, int, slice]
    ) -> Union[new_sequence.Sequence, Aligned]:
        raise TypeError(f"indexing not supported for {type(key)}, try .take_seqs()")

    @__getitem__.register
    def _(self, key: int) -> Union[new_sequence.Sequence, Aligned]:
        return self[self.parent.names[key]]

    @__getitem__.register
    def _(self, key: str) -> Union[new_sequence.Sequence, Aligned]:
        return self._make_seq(key)

    def __repr__(self) -> str:
        one_seq = self[self.parent.names[0]]
        return f"({one_seq!r}, + {self.parent.num_seqs-1} seqs)"

    def __len__(self) -> int:
        return self.parent.num_seqs

    def __iter__(self):
        for name in self.parent.names:
            yield self._make_seq(name)


class Alignment(SequenceCollection):
    def __init__(
        self,
        slice_record: new_sequence.SliceRecord = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._slice_record = (
            slice_record
            if slice_record is not None
            else new_sequence.SliceRecord(parent_len=self._seqs_data.align_len)
        )

    def _post_init(self):
        self._seqs = _IndexableSeqs(self, make_seq=self._make_aligned)

    @singledispatchmethod
    def __getitem__(self, index):
        raise NotImplementedError(f"__getitem__ not implemented for {type(index)}")

    @__getitem__.register
    def _(self, index: str):
        return self._make_aligned(seqid=index)

    @__getitem__.register
    def _(self, index: int):
        new_slice = self._slice_record[index]
        return self.__class__(
            seqs_data=self._seqs_data,
            slice_record=new_slice,
            moltype=self.moltype,
            info=self.info,
        )

    @__getitem__.register
    def _(self, index: slice):
        new_slice = self._slice_record[index]
        return self.__class__(
            seqs_data=self._seqs_data,
            slice_record=new_slice,
            moltype=self.moltype,
            info=self.info,
        )

    def __repr__(self):
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

    def __len__(self):
        return len(self._seqs_data)

    def _make_aligned(self, seqid: str) -> Aligned:
        # refactor: design
        # add slice_record as an argument to the get_view() method
        data = self._seqs_data.get_view(seqid, self._slice_record)
        return Aligned(data=data, moltype=self.moltype)

    def get_seq(
        self, seqname: str, copy_annotations: bool = False
    ) -> new_sequence.Sequence:
        """Return a Sequence object for the specified seqname.

        Parameters
        ----------
        seqname
            name of the sequence to return
        copy_annotations
            if True, only the annotations for the specified sequence are copied
            to the annotation database of the Sequence object. If False, all
            annotations are copied.
        """
        # todo
        # more efficient to constuct the sequence object directly?
        seq = self.seqs[seqname].seq
        if copy_annotations:
            seq.annotation_db = type(self.annotation_db)()
            seq.annotation_db.update(annot_db=self.annotation_db, seqids=seqname)
        else:
            seq.annotation_db = self.annotation_db

        return seq

    @c3warn.deprecated_args(
        version="2025.6", reason="naming consistency", old_new=[("seq_name", "seqname")]
    )
    def get_gapped_seq(self, seqname: str) -> new_sequence.Sequence:
        """Return a gapped Sequence object for the specified seqname.

        Notes
        -----
        This method breaks the connection to the annotation database.
        """
        return self.seqs[seqname].gapped_seq

    def rc(self):
        """Returns the reverse complement of all sequences in the alignment.
        A synonym for reverse_complement.

        Notes
        -----
        Reverse complementing the collection will break the relationship to an
        annotation_db if present.
        """
        return self.__class__(
            seqs_data=self._seqs_data,
            slice_record=self._slice_record[::-1],
            names=self.names,
            info=self.info,
            moltype=self.moltype,
            source=self.source,
        )

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
        app = get_app(app_name, **kwargs)
        return app(self)

    def iter_positions(
        self, pos_order: list = None
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
        pos_order = pos_order or range(self._seqs_data.align_len)
        for pos in pos_order:
            yield [str(self[seq][pos]) for seq in self.names]

    positions = property(iter_positions)

    def get_position_indices(
        self, f: Callable, native: bool = False, negate: bool = False
    ) -> list[int]:
        """Returns list of column indices for which f(col) is True.

        f : callable
          function that returns true/false given an alignment position
        native : boolean
          if True, and ArrayAlignment, f is provided with slice of array
          otherwise the string is used
        negate : boolean
          if True, not f() is used
        """
        # refactor:
        # type hint for f
        # implement native
        if negate:

            def new_f(x):
                return not f(x)

        else:
            new_f = f

        result = [i for i, col in enumerate(self.positions) if new_f(col)]

        return result

    def get_gap_array(self, include_ambiguity: bool = True) -> numpy.ndarray:
        """returns bool array with gap state True, False otherwise

        Parameters
        ----------
        include_ambiguity
            if True, ambiguity characters that include the gap state are
            included
        """
        result = numpy.full((self.num_seqs, len(self)), False, dtype=bool)
        for i, seqid in enumerate(self.names):
            if include_ambiguity:
                seq_array = self._seqs_data.get_gapped_seq_array(seqid=seqid)
                result[i][seq_array >= self.moltype.most_degen_alphabet().gap_index] = (
                    True
                )
            else:
                # refactor: design
                # this should be on IndelMap
                gaps = self._seqs_data.get_gaps(seqid)
                for gap_pos, cum_gap_length in gaps:
                    result[i][list(range(gap_pos, gap_pos + cum_gap_length))] = True
        return result

    def iupac_consensus(self, allow_gap: bool = True) -> str:
        """Returns string containing IUPAC consensus sequence of the alignment."""
        exclude = set() if allow_gap else set(self.moltype.gaps)
        consensus = []
        degen = self.moltype.degenerate_from_seq
        for col in self.iter_positions():
            col = set(col) - exclude
            consensus.append(degen("".join(col)))
        return "".join(consensus)

    def majority_consensus(self) -> new_sequence.Sequence:
        """Returns consensus sequence containing most frequent item at each
        position."""
        states = []
        data = zip(*map(str, self.seqs))
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

        warn
            warns if motif_length > 1 and alignment trimmed to produce
            motif columns
        """
        length = (len(self) // motif_length) * motif_length
        if warn and len(self) != length:
            warnings.warn(f"trimmed {len(self) - length}", UserWarning)

        data = list(self.to_dict().values())
        alpha = self.moltype.alphabet.get_kmer_alphabet(motif_length)
        all_motifs = set()
        exclude_chars = set()
        if not allow_gap:
            exclude_chars.update(self.moltype.gap)

        if not include_ambiguity:
            ambigs = [c for c, v in self.moltype.ambiguities.items() if len(v) > 1]
            exclude_chars.update(ambigs)

        result = []
        for i in range(0, len(self) - motif_length + 1, motif_length):
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

        result = MotifCountsArray(result, alpha)
        return result

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
            warnings.warn(f"trimmed {len(self) - length}", UserWarning)

        counts = []
        motifs = set()
        for name in self.names:
            seq = self.get_gapped_seq(name)
            c = seq.counts(
                motif_length=motif_length,
                include_ambiguity=include_ambiguity,
                allow_gap=allow_gap,
                exclude_unobserved=exclude_unobserved,
            )
            motifs.update(c.keys())
            counts.append(c)

        if not exclude_unobserved:
            motifs.update(self.moltype.alphabet.get_kmer_alphabet(motif_length))

        motifs = list(sorted(motifs))
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
        result = darr.wrap(result)
        return result

    def count_gaps_per_seq(
        self,
        induced_by: bool = False,
        unique: bool = False,
        include_ambiguity: bool = True,
        drawable: bool = False,
    ):
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
            gap_array = gap_array[:, gap_cols] == False
        elif induced_by:
            # identify all columns with gap opposite
            gap_cols = gap_array.sum(axis=0) > 0
            gap_array = gap_array[:, gap_cols] == False
        else:
            gap_cols = gap_array.sum(axis=0) > 0
            gap_array = gap_array[:, gap_cols]

        result = gap_array.sum(axis=1)
        result = darr.wrap(result)
        if drawable:
            drawable = drawable.lower()
            trace_name = (
                os.path.basename(self.info.source) if self.info.source else None
            )
            draw = Drawable("Gaps Per Sequence", showlegend=False)
            draw.layout |= dict(yaxis=dict(title="Gap counts"))
            if drawable == "bar":
                trace = UnionDict(type="bar", y=result.array, x=self.names)
            else:
                trace = UnionDict(
                    type=drawable, y=result.array, text=self.names, name=trace_name
                )

            draw.add_trace(trace)
            result = draw.bound_to(result)

        return result

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

    def matching_ref(self, ref_name, gap_fraction, gap_run):
        """Returns new alignment with seqs well aligned with a reference.

        gap_fraction = fraction of positions that either have a gap in the
            template but not in the seq or in the seq but not in the template
        gap_run = number of consecutive gaps tolerated in query relative to
            sequence or sequence relative to query
        """
        template = self.seqs[ref_name]
        gap_filter = make_gap_filter(template, gap_fraction, gap_run)
        return self.take_seqs_if(gap_filter)

    def sliding_windows(
        self, window: int, step: int, start: OptInt = None, end: OptInt = None
    ):
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

    def _get_seq_features(
        self,
        *,
        seqid: Optional[str] = None,
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
        if self.annotation_db is None:
            return None

        seqid_to_seqname = {seq.parent_coordinates()[0]: seq.name for seq in self.seqs}

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
            raise ValueError(f"unknown {seqid=}")

        for seqid in seqids:
            seqname = seqid_to_seqname[seqid]
            seq = self.seqs[seqname]
            # we use parent seqid
            parent_id, start, stop, _ = seq.parent_coordinates()
            offset = seq.annotation_offset

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
        seqid: Optional[str] = None,
        biotype: Optional[str] = None,
        name: Optional[str] = None,
        on_alignment: Optional[bool] = None,
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
        # we only do on-alignment in here
        if not on_alignment:
            kwargs = {
                k: v for k, v in locals().items() if k not in ("self", "__class__")
            }
            kwargs.pop("on_alignment")
            yield from self._get_seq_features(**kwargs)

        if on_alignment == False:
            return

        if self.annotation_db is None:
            anno_db = merged_db_collection(self._seqs_data)
            self.annotation_db = anno_db

        if self.annotation_db is None:
            return None

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
                raise RuntimeError(f"{on_alignment=} {feature=}")
            if seq_map is None:
                seq_map = self.seqs[0].map.to_feature_map()
                *_, strand = self.seqs[0].data.parent_coordinates()

            spans = numpy.array(feature["spans"])
            spans = seq_map.relative_position(spans)
            feature["spans"] = spans.tolist()
            # and if i've been reversed...?
            feature["strand"] = "-" if strand == -1 else "+"
            yield self.make_feature(feature=feature, on_alignment=on_al)

    def to_pretty(self, name_order=None, wrap=None):
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
        display_names = dict([(n, name_template.format(n)) for n in names])

        def make_line(label, seq):
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
                    )
                )

            result.append("")

        if not result[-1]:
            del result[-1]

        return "\n".join(result)

    def to_html(
        self,
        name_order: Optional[typing.Sequence[str]] = None,
        wrap: int = 60,
        limit: Optional[int] = None,
        ref_name: str = "longest",
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
        ---------

        In a jupyter notebook, this code is used to provide the representation.

        .. code-block:: python

            aln  # is rendered by jupyter

        You can directly use the result for display in a notebook as

        .. code-block:: python

            from IPython.core.display import HTML

            HTML(aln.to_html())
        """
        css, styles = self.moltype.get_css_style(
            colors=colors, font_size=font_size, font_family=font_family
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
                raise ValueError(f"Unknown sequence name {ref_name}")
            ref = ref_name

        name_order.remove(ref)
        name_order.insert(0, ref)

        if limit is None:
            names, output = selected._get_raw_pretty(name_order)
        else:
            names, output = selected[:limit]._get_raw_pretty(name_order)

        gaps = "".join(selected.moltype.gaps)
        refname = names[0]
        refseq = output[refname]
        seqlen = len(refseq)
        start_gap = re.search(f"^[{gaps}]+", "".join(refseq))
        end_gap = re.search(f"[{gaps}]+$", "".join(refseq))
        start = 0 if start_gap is None else start_gap.end()
        end = len(refseq) if end_gap is None else end_gap.start()
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
            for n, s in zip(names, seqblock):
                s = "".join(s)
                row = "".join([label_ % n, seq_ % s])
                table.append(f"<tr>{row}</tr>")
        table.append("</table>")
        if (
            limit
            and limit < len(selected)
            or name_order
            and len(name_order) < len(selected.names)
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
        positions = list(zip(*seqs))

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
            dict(num_seqs=int, num_pos=int, wrap=int, ref_name=str),
        )
        settings.update(env_vals)
        return self.to_html(
            name_order=self.names[: settings["num_seqs"]],
            ref_name=settings["ref_name"],
            limit=settings["num_pos"],
            wrap=settings["wrap"],
        )


@singledispatch
def make_aligned_seqs(
    data: Union[dict[str, StrORBytesORArray], AlignedSeqsData],
    *,
    moltype: str,
    info: dict = None,
    offset: dict[str, int] = None,
    annotation_db: SupportsFeatures = None,
) -> Alignment:
    raise NotImplementedError(f"make_aligned_seqs not implemented for {type(data)}")


@make_aligned_seqs.register
def _(
    data: dict,
    *,
    moltype: str,
    info: dict = None,
    offset: dict[str, int] = None,
    annotation_db: SupportsFeatures = None,
) -> Alignment:
    # todo: implement a coerce to AlignedSeqsData dict function
    moltype = new_moltype.get_moltype(moltype)
    alphabet = moltype.most_degen_alphabet()
    aligned_data = AlignedSeqsData.from_aligned_seqs(
        data=data,
        alphabet=alphabet,
        offset=offset,
    )
    annotation_db = annotation_db or merged_db_collection(data)
    return make_aligned_seqs(
        aligned_data, moltype=moltype, info=info, annotation_db=annotation_db
    )


@make_aligned_seqs.register
def _(
    data: AlignedSeqsDataABC,
    *,
    moltype: str,
    info: dict = None,
    offset: dict[str, int] = None,
    annotation_db: SupportsFeatures = None,
) -> Alignment:
    moltype = new_moltype.get_moltype(moltype)
    if not moltype.is_compatible_alphabet(data.alphabet):
        raise ValueError(
            f"Provided moltype: {moltype.label} is not compatible with AlignedSeqsData",
            f" alphabet: {data.alphabet}",
        )

    return Alignment(
        seqs_data=data, moltype=moltype, info=info, annotation_db=annotation_db
    )
