from __future__ import annotations

import hashlib
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, cast

import numba
import numpy
import numpy.typing as npt
from typing_extensions import Self

import cogent3.core.alphabet as c3_alphabet
import cogent3.core.moltype as c3_moltype
import cogent3.core.sequence as c3_sequence
from cogent3.core.location import IndelMap
from cogent3.core.seqview import AlignedDataView, AlignedDataViewABC, SeqDataView
from cogent3.util.misc import extend_docstring_from

if TYPE_CHECKING:
    from collections.abc import Collection, Mapping
    from collections.abc import Sequence as PySeq

    from cogent3.core.seqview import SeqViewABC
    from cogent3.core.slice_record import SliceRecord

NumpyIntArrayType = npt.NDArray[numpy.integer]


class SeqsDataABC(ABC):
    """Abstract base class for respresenting the storage object for sequences underlying
    a SequenceCollection.
    """

    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        *,
        data: dict[str, str | bytes | NumpyIntArrayType],
        alphabet: c3_alphabet.AlphabetABC[Any],
        offset: dict[str, int] | None = None,
        check: bool = True,
        reversed_seqs: set[str] | None = None,
    ) -> None: ...
    @classmethod
    @abstractmethod
    def from_seqs(
        cls,
        *,
        data: Mapping[str, str | bytes | NumpyIntArrayType],
        alphabet: c3_alphabet.CharAlphabet[Any],
        **kwargs: Any,
    ) -> Self: ...

    @abstractmethod
    def __eq__(self, value: object) -> bool: ...

    def __ne__(self, other: object) -> bool:
        return not self == other

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
    def alphabet(self) -> c3_alphabet.CharAlphabet[Any]: ...

    @property
    @abstractmethod
    def offset(self) -> dict[str, int]: ...

    @abstractmethod
    def get_seq_array(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType: ...

    def get_seq_str(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> str:
        """Return ungapped sequence corresponding to seqid as a string.
        start/stop are in sequence coordinates. Excludes gaps."""
        return self.alphabet.from_indices(
            self.get_seq_array(seqid=seqid, start=start, stop=stop, step=step),
        )

    def get_seq_bytes(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> bytes:
        """Return ungapped sequence corresponding to seqid as a bytes string.
        start/stop are in sequence coordinates. Excludes gaps."""
        return self.get_seq_str(seqid=seqid, start=start, stop=stop, step=step).encode(
            "utf8",
        )

    @abstractmethod
    def get_view(self, seqid: str) -> SeqViewABC: ...

    @abstractmethod
    def to_alphabet(
        self,
        alphabet: c3_alphabet.CharAlphabet[Any],
        check_valid: bool = True,
    ) -> Self: ...

    @abstractmethod
    def add_seqs(
        self,
        seqs: Mapping[str, str | bytes | NumpyIntArrayType],
        force_unique_keys: bool = True,
        offset: dict[str, int] | None = None,
    ) -> SeqsDataABC: ...

    @abstractmethod
    def __len__(self) -> int: ...

    @abstractmethod
    def __getitem__(
        self,
        index: str | int,
    ) -> c3_sequence.Sequence | SeqViewABC: ...

    @abstractmethod
    def copy(self, **kwargs: Any) -> Self: ...

    @abstractmethod
    def get_hash(self, seqid: str) -> str | None: ...


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

    __slots__ = ("_alphabet", "_data", "_hashes", "_offset", "_reversed")

    def __init__(
        self,
        *,
        data: Mapping[str, str | bytes | NumpyIntArrayType],
        alphabet: c3_alphabet.CharAlphabet[Any],
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
        AlphabetError if the check fails
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
        self._data: dict[str, NumpyIntArrayType] = {}
        self._hashes: dict[str, str] = {}
        for name, seq in data.items():
            arr = self._alphabet.to_indices(seq)
            self._hashes[name] = array_hash64(arr)
            arr.flags.writeable = False
            self._data[str(name)] = arr

    def __eq__(self, other: object) -> bool:
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

    def __len__(self) -> int:
        return len(self.names)

    def __getitem__(self, index: str | int) -> SeqViewABC:
        if isinstance(index, int):
            return self[self.names[index]]
        if isinstance(index, str):
            return self.get_view(seqid=index)

        msg = f"__getitem__ not implemented for {type(index)}"
        raise NotImplementedError(msg)

    @classmethod
    def from_seqs(
        cls,
        *,
        data: Mapping[str, str | bytes | NumpyIntArrayType],
        alphabet: c3_alphabet.CharAlphabet[Any],
        **kwargs: Any,
    ) -> Self:
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
    def alphabet(self) -> c3_alphabet.CharAlphabet[Any]:
        """the character alphabet for validating, encoding, decoding sequences"""
        return self._alphabet

    @property
    def offset(self) -> dict[str, int]:
        """annotation offsets for each sequence"""
        return {name: self._offset.get(name, 0) for name in self.names}

    def get_seq_length(self, seqid: str) -> int:
        """return length for seqid"""
        return self._data[seqid].shape[0]

    def get_hash(self, seqid: str) -> str | None:
        """returns hash of seqid"""
        return self._hashes.get(seqid)

    def get_seq_array(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType:
        start = start or 0
        stop = stop if stop is not None else self.get_seq_length(seqid)
        step = step or 1

        if start < 0 or stop < 0 or step < 1:
            msg = f"{start=}, {stop=}, {step=} not >= 1"
            raise ValueError(msg)

        out_len = (stop - start + step - 1) // step

        seq = numpy.empty(out_len, dtype=self.alphabet.dtype)
        seq[:] = self._data[seqid][start:stop:step]

        return seq

    def get_view(self, seqid: str) -> SeqDataView:
        """reurns view of sequence data for seqid"""
        seq_len = len(self._data[seqid])
        return SeqDataView(
            parent=self,
            seqid=seqid,
            parent_len=seq_len,
            alphabet=self.alphabet,
        )

    def add_seqs(
        self,
        seqs: Mapping[str, str | bytes | NumpyIntArrayType],
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

    def copy(self, **kwargs: Any) -> Self:
        """shallow copy of self

        Notes
        -----
        kwargs are passed to constructor and will over-ride existing values
        """
        init_args: dict[str, Any] = {
            "data": self._data,
            "alphabet": self._alphabet,
            "offset": self._offset,
            "reversed_seqs": self._reversed,
            **kwargs,
        }
        return self.__class__(**init_args)

    def to_alphabet(
        self,
        alphabet: c3_alphabet.CharAlphabet[Any],
        check_valid: bool = True,
    ) -> Self:
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
                check_valid=check_valid,
            )
            new_data[seqid] = as_new_alpha

        return self.copy(
            data=new_data,
            alphabet=alphabet,
            check=False,
        )


class AlignedSeqsDataABC(SeqsDataABC):
    """Abstract base class for respresenting the storage object for sequences underlying
    a Alignment.
    """

    # all methods that are from SeqsDataABC should work in sequence coordinates,
    # all methods unique to AlignedSeqsDataABC should work in aligned coordinates.
    # all indices provided to AlignedSeqsDataABC should be on the plus strand.
    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        *,
        gapped_seqs: NumpyIntArrayType,
        names: Collection[str],
        alphabet: c3_alphabet.CharAlphabet[Any],
        ungapped_seqs: dict[str, NumpyIntArrayType] | None = None,
        gaps: Mapping[str, NumpyIntArrayType] | None = None,
        offset: dict[str, int] | None = None,
        align_len: int | None = None,
        check: bool = True,
        reversed_seqs: set[str] | None = None,
    ) -> None: ...

    @classmethod
    @abstractmethod
    def from_seqs_and_gaps(
        cls,
        *,
        seqs: Mapping[str, str | bytes | NumpyIntArrayType],
        gaps: Mapping[str, NumpyIntArrayType],
        alphabet: c3_alphabet.CharAlphabet[Any],
    ) -> Self: ...

    @classmethod
    @abstractmethod
    def from_names_and_array(
        cls,
        *,
        names: Collection[str],
        data: NumpyIntArrayType,
        alphabet: c3_alphabet.CharAlphabet[Any],
    ) -> Self: ...

    @property
    @abstractmethod
    def align_len(self) -> int: ...

    @abstractmethod
    def get_view(
        self,
        seqid: str | int,
        slice_record: SliceRecord | None = None,
    ) -> AlignedDataViewABC:
        # overriding the SeqsDataABC method as we support directly
        # providing the slice_record instance
        ...

    @abstractmethod
    def get_gapped_seq_array(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def get_gapped_seq_str(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> str: ...

    @abstractmethod
    def get_gapped_seq_bytes(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> bytes: ...

    @abstractmethod
    def get_ungapped(
        self,
        name_map: Mapping[str, str],
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> tuple[dict[str, NumpyIntArrayType], dict[str, Any]]:
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
    def get_pos_range(
        self,
        names: PySeq[str],
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def get_positions(
        self,
        names: PySeq[str],
        positions: PySeq[int] | NumpyIntArrayType,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def variable_positions(
        self,
        names: PySeq[str],
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def get_gaps(self, seqid: str) -> NumpyIntArrayType: ...


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
        "_hashes",
        "_name_to_index",
        "_names",
        "_offset",
        "_reversed",
        "_ungapped",
    )

    def __init__(
        self,
        *,
        gapped_seqs: NumpyIntArrayType,
        names: Collection[str],
        alphabet: c3_alphabet.CharAlphabet[Any],
        ungapped_seqs: dict[str, NumpyIntArrayType] | None = None,
        gaps: Mapping[str, NumpyIntArrayType] | None = None,
        offset: dict[str, int] | None = None,
        align_len: int | None = None,
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
        self._gaps = dict(gaps) if gaps else {}
        self._hashes: dict[str, str] = {}
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

    def __eq__(self, other: object) -> bool:
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

        return bool(numpy.all(self._gapped == other._gapped))

    def __len__(self) -> int:
        return self.align_len

    def __getitem__(self, index: str | int) -> AlignedDataViewABC:
        return self.get_view(index)

    @property
    def names(self) -> tuple[str, ...]:
        """returns the names of the sequences in the storage"""
        return self._names

    @property
    def reversed_seqs(self) -> frozenset[str]:
        """names of sequences that are reverse complemented"""
        return self._reversed

    @property
    def alphabet(self) -> c3_alphabet.CharAlphabet[Any]:
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

    @classmethod
    def from_seqs(
        cls,
        *,
        data: Mapping[str, str | bytes | NumpyIntArrayType],
        alphabet: c3_alphabet.CharAlphabet[Any],
        **kwargs: Any,
    ) -> Self:
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
            array_seqs[i] = alphabet.to_indices(data[name], validate=True)

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
        seqs: Mapping[str, str | bytes | NumpyIntArrayType],
        gaps: Mapping[str, NumpyIntArrayType],
        alphabet: c3_alphabet.CharAlphabet[Any],
        **kwargs: Any,
    ) -> Self:
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
        names: tuple[str, ...] = tuple(kwargs.pop("names", seqs.keys()))
        if not names:
            msg = "seqs cannot be empty"
            raise ValueError(msg)

        align_len: int | None = kwargs.pop("align_len", None)
        if align_len is None:
            align_len = _gapped_seq_len(seqs[names[0]], gaps[names[0]])

        gapped_seqs = numpy.empty((len(names), align_len), dtype=alphabet.dtype)
        seq_arrays: dict[str, NumpyIntArrayType] = {}
        for i, name in enumerate(names):
            seq = alphabet.to_indices(seqs[name])
            seq_arrays[name] = seq
            if name not in gaps:
                msg = f"Missing gap data for sequence {name!r}"
                raise ValueError(msg)
            gapped_seqs[i] = compose_gapped_seq(
                seq, gaps[name], cast("int", alphabet.gap_index)
            )
            assert len(gapped_seqs[i]) == align_len, "aligned lengths do not match"

        gapped_seqs.flags.writeable = False
        return cls(
            ungapped_seqs=seq_arrays,
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
        names: Collection[str],
        data: NumpyIntArrayType,
        alphabet: c3_alphabet.CharAlphabet[Any],
    ) -> Self:
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

    def get_hash(self, seqid: str) -> str:
        """returns hash of seqid"""
        if seqid not in self._hashes:
            arr = self.get_gapped_seq_array(seqid=seqid)
            self._hashes[seqid] = array_hash64(arr)
        return self._hashes[seqid]

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

    def _get_ungapped(self, seqid: str) -> NumpyIntArrayType:
        if seqid not in self._ungapped:
            self._make_gaps_and_ungapped(seqid)
        return self._ungapped[seqid]

    def get_seq_length(self, seqid: str) -> int:
        """return length of the unaligned seq for seqid"""
        return len(self._get_ungapped(seqid))

    def get_view(
        self,
        seqid: str | int,
        slice_record: SliceRecord | None = None,
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
        if isinstance(seqid, int):
            seqid = self.names[seqid]

        return AlignedDataView(
            parent=self,
            seqid=seqid,
            alphabet=self.alphabet,
            slice_record=slice_record,
        )

    def _get_gaps(self, seqid: str) -> NumpyIntArrayType:
        if seqid not in self._gaps:
            self._make_gaps_and_ungapped(seqid)
        return self._gaps[seqid]

    def get_gaps(self, seqid: str) -> NumpyIntArrayType:
        """returns the gap data for seqid"""
        return self._get_gaps(seqid)

    def get_seq_array(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType:
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
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType:
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

    def get_gapped_seq_str(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> str:
        """Return sequence corresponding to seqid as a string.
        start/stop are in alignment coordinates. Includes gaps."""
        return self.alphabet.from_indices(
            self.get_gapped_seq_array(seqid=seqid, start=start, stop=stop, step=step),
        )

    def get_gapped_seq_bytes(
        self,
        *,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
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
        name_map: Mapping[str, str],
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> tuple[dict[str, NumpyIntArrayType], dict[str, Any]]:
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
        seqs: dict[str, NumpyIntArrayType] = {}
        for i, name in enumerate(names):
            seq: NumpyIntArrayType = seq_array[i]
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

    def get_pos_range(
        self,
        names: PySeq[str],
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType:
        """returns an array of the selected positions for names."""
        start = start or 0
        stop = stop or self.align_len
        step = step or 1
        if start < 0 or stop < 0 or step < 1:
            msg = f"{start=}, {stop=}, {step=} not >= 1"
            raise ValueError(msg)

        indices = tuple(self._name_to_index[name] for name in names)
        if abs((start - stop) // step) == self.align_len:
            array_seqs = self._gapped[indices, :]
        else:
            array_seqs = self._gapped[indices, start:stop:step]

        return array_seqs

    def get_positions(
        self,
        names: PySeq[str],
        positions: PySeq[int] | NumpyIntArrayType,
    ) -> NumpyIntArrayType:
        """returns alignment positions for names

        Parameters
        ----------
        names
            series of sequence names
        positions
            indices lying within self

        Returns
        -------
            2D numpy.array, oriented by sequence

        Raises
        ------
        IndexError if a provided position is negative or
        greater then alignment length.
        """
        if diff := set(names) - set(self.names):
            msg = f"these names not present {diff}"
            raise ValueError(msg)

        min_index, max_index = numpy.min(positions), numpy.max(positions)
        if min_index < 0 or max_index > self.align_len:
            msg = f"Out of range: {min_index=} and / or {max_index=}"
            raise IndexError(msg)

        seq_indices = numpy.array([self._name_to_index[n] for n in names])
        return self._gapped[numpy.ix_(seq_indices, numpy.asarray(positions))]

    def copy(self, **kwargs: Any) -> Self:
        """shallow copy of self

        Notes
        -----
        kwargs are passed to constructor and will over-ride existing values
        """
        init_args: dict[str, Any] = {
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

    def add_seqs(
        self,
        seqs: Mapping[str, str | bytes | NumpyIntArrayType],
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
            seq = self.alphabet.to_indices(seq, validate=True)
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
        alphabet: c3_alphabet.CharAlphabet[Any],
        check_valid: bool = True,
    ) -> Self:
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
                check_valid=check_valid,
            )
            gapped[i] = as_new_alpha

        return self.__class__(
            gapped_seqs=gapped,
            alphabet=alphabet,
            offset=self._offset,
            names=self.names,
        )

    def variable_positions(
        self,
        names: PySeq[str],
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> NumpyIntArrayType:
        """returns absolute indices of positions that have more than one state

        Parameters
        ----------
        names
            selected seqids
        start
            absolute start
        stop
            absolute stop
        step
            step

        Returns
        -------
        Absolute indices (as distinct from an index relative to start) of
        variable positions.
        """
        start = start or 0
        if len(names) < 2:
            return numpy.array([])

        array_seqs = self.get_pos_range(names, start=start, stop=stop, step=step)
        if array_seqs.size == 0:
            return numpy.array([])

        indices = (array_seqs != array_seqs[0]).any(axis=0)
        return numpy.where(indices)[0] + start


def _gapped_seq_len(
    seq: str | bytes | NumpyIntArrayType, gap_map: IndelMap | NumpyIntArrayType
) -> int:
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


@numba.jit(cache=True)
def compose_gapped_seq(
    ungapped_seq: NumpyIntArrayType,
    gaps: NumpyIntArrayType,
    gap_index: int,
) -> NumpyIntArrayType:  # pragma: no cover
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


def decompose_gapped_seq(
    seq: str | bytes | NumpyIntArrayType | c3_sequence.Sequence,
    *,
    alphabet: c3_alphabet.AlphabetABC[Any],
    missing_as_gap: bool = True,
) -> tuple[NumpyIntArrayType, NumpyIntArrayType]:
    """
    Takes a sequence with (or without) gaps and returns an ungapped sequence
    and a map of the position and length of gaps in the original parent sequence
    """

    if isinstance(seq, bytes):
        seq = seq.decode("utf-8")

    if isinstance(seq, str):
        if not alphabet.is_valid(seq):
            msg = f"Sequence is invalid for alphabet {alphabet}"
            raise c3_alphabet.AlphabetError(msg)
        seq = alphabet.to_indices(seq)

    if isinstance(seq, c3_sequence.Sequence):
        seq = numpy.array(seq)

    if isinstance(seq, numpy.ndarray):
        if missing_as_gap and alphabet.missing_index:
            missing_index = int(numpy.uint8(alphabet.missing_index))
        else:
            missing_index = -1
        return decompose_gapped_seq_array(
            seq.astype(alphabet.dtype),
            cast("int", alphabet.gap_index),
            missing_index=missing_index,
        )

    msg = f"decompose_gapped_seq not implemented for type {type(seq)}"
    raise NotImplementedError(
        msg,
    )


@numba.jit(cache=True)
def decompose_gapped_seq_array(
    seq: NumpyIntArrayType,
    gap_index: int,
    missing_index: int = -1,
) -> tuple[NumpyIntArrayType, NumpyIntArrayType]:  # pragma: no cover
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
        gapped = base in (gap_index, missing_index)
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
        gapped = base in (gap_index, missing_index)
        if not gapped:
            ungapped[seqpos] = base
            seqpos += 1

    return ungapped, gap_coords


def array_hash64(data: npt.NDArray[numpy.number]) -> str:
    """returns 64-bit hash of numpy array.

    Notes
    -----
    This function does not introduce randomisation and so
    is reproducible between processes.
    """
    h = hashlib.md5(data.tobytes(), usedforsecurity=False)
    return h.hexdigest()


if __name__ == "__main__":
    seqs_data = SeqsData(data={}, alphabet=c3_moltype.DNA.alphabet)

    aln_seqs_data = AlignedSeqsData(
        gapped_seqs=numpy.array([1]), names=("a",), alphabet=c3_moltype.DNA.alphabet
    )
