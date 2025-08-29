import functools
import itertools
import json
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable, Sized
from collections.abc import Sequence as PySeq
from typing import TYPE_CHECKING, Any, Generic, TypeVar, cast

import numba
import numpy
import numpy.typing as npt
from typing_extensions import Self

from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import get_object_provenance

if TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.moltype import MolType

NumpyIntArrayType = npt.NDArray[numpy.integer]


StrOrBytes = TypeVar("StrOrBytes", str, bytes)


def _coerce_to_type(
    orig: StrOrBytes | PySeq[StrOrBytes],
    text: str | bytes,
) -> StrOrBytes:
    if isinstance(orig, str):
        return cast("StrOrBytes", str(text))
    if isinstance(orig, (bytes, bytearray)):
        return cast(
            "StrOrBytes", text.encode("utf8") if isinstance(text, str) else text
        )

    first = orig[0]
    if isinstance(first, str):
        return cast("StrOrBytes", str(text))
    if isinstance(first, (bytes, bytearray)):
        return cast(
            "StrOrBytes", text.encode("utf8") if isinstance(text, str) else text
        )

    msg = f"{type(orig)} is invalid"
    raise TypeError(msg)


class AlphabetError(TypeError): ...


class convert_alphabet:
    """convert one character set into another"""

    def __init__(
        self,
        src: bytes,
        dest: bytes,
        delete: bytes | None = None,
        allow_duplicates: bool = False,
    ) -> None:
        """
        Parameters
        ----------
        src
            unique series of characters
        new
            characters in src will be mapped to these in order
        delete
            characters to be deleted from the result
        allow_duplicates
            whether to allow duplicates in src and dest

        Notes
        -----
        orig and new must be the same lengths.

        Examples
        --------
        >>> dna_to_ry = convert_alphabet(src=b"TCAG", dest=b"YYRR")
        >>> dna_to_ry(b"GCTA")
        b'RYYR'
        >>> dna_to_ry = convert_alphabet(src=b"TCAG", dest=b"YYRR", delete=b"-")
        >>> dna_to_ry(b"GC-TA")
        b'RYYR'
        """
        if len(src) != len(dest):
            msg = f"length of src={len(src)} != length of new {len(dest)}"
            raise ValueError(msg)

        if not allow_duplicates:
            consistent_words(src, length=1)

        self._table = b"".maketrans(src, dest)
        self._delete = delete or b""

    def __call__(self, seq: bytes) -> bytes:
        return seq.translate(self._table, delete=self._delete)


class AlphabetABC(ABC, Generic[StrOrBytes]):
    def __init__(self) -> None:
        self.dtype: type[numpy.unsignedinteger]

    @property
    @abstractmethod
    def num_canonical(self) -> int: ...

    @property
    @abstractmethod
    def motif_len(self) -> int: ...

    @abstractmethod
    def is_valid(self, seq: str | bytes | NumpyIntArrayType) -> bool: ...

    @abstractmethod
    def to_indices(
        self,
        seq: str
        | bytes
        | list[str | bytes]
        | tuple[str | bytes, ...]
        | NumpyIntArrayType,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def from_indices(
        self, seq: str | bytes | NumpyIntArrayType
    ) -> str | list[str] | NumpyIntArrayType: ...

    @property
    @abstractmethod
    def gap_char(self) -> str | bytes | None: ...

    @property
    @abstractmethod
    def gap_index(self) -> int | None: ...

    @property
    @abstractmethod
    def missing_char(self) -> str | bytes | None: ...

    @property
    @abstractmethod
    def missing_index(self) -> int | None: ...

    @abstractmethod
    def to_rich_dict(self, for_pickle: bool = False) -> dict[str, Any]: ...

    @abstractmethod
    def to_json(self) -> str: ...

    @classmethod
    @abstractmethod
    def from_rich_dict(cls, data: dict[str, Any]) -> Self: ...

    def __deepcopy__(self, memo: dict[int, Self]) -> Self:
        return type(self).from_rich_dict(self.to_rich_dict())

    @property
    def moltype(self) -> "MolType[StrOrBytes] | None":
        return _alphabet_moltype_map.get(self)


class MonomerAlphabetABC(ABC, Generic[StrOrBytes]):
    @abstractmethod
    def get_kmer_alphabet(
        self,
        k: int,
        include_gap: bool = True,
    ) -> "KmerAlphabet[StrOrBytes] | CharAlphabet[StrOrBytes]": ...

    @abstractmethod
    def as_bytes(self) -> bytes: ...

    @abstractmethod
    def convert_seq_array_to(
        self,
        *,
        alphabet: "CharAlphabet[Any]",
        seq: NumpyIntArrayType,
        check_valid: bool = True,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def with_gap_motif(
        self,
        gap_char: str
        | bytes = "-",  # Strictly speaking this should use StrOrBytes but generics don't play well with these default arguments
        missing_char: str | bytes = "?",
        include_missing: bool = False,
        gap_as_state: bool = False,
    ) -> "MonomerAlphabetABC[StrOrBytes]": ...


def get_array_type(
    num_elements: int,
) -> type[numpy.unsignedinteger]:
    """Returns the smallest numpy integer dtype that can contain elements
    within num_elements.
    """
    if num_elements <= 2**8:
        return numpy.uint8
    if num_elements <= 2**16:
        return numpy.uint16
    if num_elements <= 2**32:
        return numpy.uint32
    if num_elements <= 2**64:
        return numpy.uint64

    msg = f"{num_elements} is too big for 64-bit integer"
    raise NotImplementedError(msg)


def consistent_words(
    words: PySeq[Sized] | PySeq[int],
    length: int | None = None,
) -> None:
    """make sure all alphabet elements are unique and have the same length"""
    if not words:
        msg = "no words provided"
        raise ValueError(msg)

    if len(set(words)) != len(words):
        msg = f"duplicate elements in {words}"
        raise ValueError(msg)

    lengths = (
        {1} if isinstance(words[0], int) else {len(cast("Sized", w)) for w in words}
    )
    if len(lengths) != 1:
        msg = f"mixed lengths {lengths} in {words}"
        raise ValueError(msg)
    actual_length = next(iter(lengths))
    if length and actual_length != length:
        msg = f"word length {actual_length} does not match expected {length}"
        raise ValueError(msg)


class bytes_to_array:
    """wrapper around convert_alphabet. It defines a linear mapping from provided
    characters to uint8. The resulting object is callable, taking a bytes object
    and returning a numpy array."""

    def __init__(
        self,
        chars: bytes,
        dtype: type[numpy.unsignedinteger],
        delete: bytes | None = None,
    ) -> None:
        """
        Parameters
        ----------
        chars
            unique series of characters
        delete
            characters to be deleted from the result
        """
        # we want a bytes translation map
        self._converter = convert_alphabet(
            chars,
            bytes(bytearray(range(len(chars)))),
            delete=delete,
        )
        self.dtype = dtype

    def __call__(self, seq: bytes) -> NumpyIntArrayType:
        b = self._converter(seq)
        return numpy.array(memoryview(b), dtype=self.dtype)


class array_to_bytes:
    """wrapper around convert_alphabet. It defines a linear mapping from uint8
    integers to the provided characters. The resulting object is callable,
    taking a numpy array and returning a bytes object.
    """

    def __init__(self, chars: bytes) -> None:
        # we want a bytes translation map
        self._converter = convert_alphabet(
            bytes(bytearray(range(len(chars)))),
            chars,
        )

    def __call__(self, seq: NumpyIntArrayType) -> bytes:
        b = self._converter(seq.tobytes())
        return bytes(bytearray(b))


class CharAlphabet(
    tuple[StrOrBytes, ...],
    Generic[StrOrBytes],
    AlphabetABC[StrOrBytes],
    MonomerAlphabetABC[StrOrBytes],
):
    """representing fundamental monomer character sets.

    Notes
    -----
    Provides methods for efficient conversion between characters and integers
    from fundamental types of strings, bytes and numpy arrays.
    """

    __slots__ = ()

    def __new__(
        cls,
        chars: StrOrBytes | PySeq[StrOrBytes],
        gap: StrOrBytes | None = None,
        missing: StrOrBytes | None = None,
    ) -> "CharAlphabet[StrOrBytes]":
        """
        Parameters
        ----------
        chars
            the characters in the alphabet
        gap
            character representing the gap character (universally '-' in cogent3)
        missing
            character representing the missing data, typically '?'.
        """
        if not chars:
            msg = f"cannot create empty {cls.__name__!r}"
            raise ValueError(msg)

        if gap is not None and _coerce_to_type(chars, gap) not in chars:
            msg = f"gap {gap!r} not in chars {chars!r}"
            raise ValueError(msg)

        if missing is not None and _coerce_to_type(chars, missing) not in chars:
            msg = f"char missing={_coerce_to_type(chars, missing)!r} not in {chars!r}"
            raise ValueError(msg)

        consistent_words(chars, length=1)
        if isinstance(chars, bytes):
            # we need to convert to tuple in a way that preserves elements
            # as bytes
            chars = cast("PySeq[StrOrBytes]", tuple(bytes([c]) for c in chars))

        return tuple.__new__(cls, chars)

    def __init__(
        self,
        chars: PySeq[StrOrBytes] | StrOrBytes,  # noqa: ARG002
        gap: StrOrBytes | None = None,
        missing: StrOrBytes | None = None,
    ) -> None:
        self.dtype = get_array_type(len(self))
        self._gap_char: StrOrBytes | None = gap
        self._gap_index = self.dtype(self.index(gap)) if gap else None
        self._missing_char: StrOrBytes | None = missing
        self._missing_index = self.dtype(self.index(missing)) if missing else None

        # the number of canonical states are non-gap, non-missing states
        adj = (1 if gap else 0) + (1 if missing else 0)

        self._num_canonical = len(self) - adj

        self._chars: set[StrOrBytes] = set(self)  # for quick lookup
        byte_chars = self.as_bytes()
        self._bytes2arr = bytes_to_array(byte_chars, dtype=self.dtype)
        self._arr2bytes = array_to_bytes(byte_chars)

    @property
    def num_canonical(self) -> int:
        return self._num_canonical

    @property
    def gap_char(self) -> StrOrBytes | None:
        return self._gap_char

    @property
    def gap_index(self) -> int | None:
        return int(self._gap_index) if self._gap_index is not None else None

    @property
    def missing_char(self) -> StrOrBytes | None:
        return self._missing_char

    @property
    def missing_index(self) -> int | None:
        return int(self._missing_index) if self._missing_index is not None else None

    @property
    def motif_len(self) -> int:
        # this is enforced by consistent_words(length=1) in the constructor
        return 1

    def to_indices(
        self,
        seq: str
        | bytes
        | tuple[str | bytes, ...]
        | list[str | bytes]
        | NumpyIntArrayType,
    ) -> NumpyIntArrayType:
        """returns a sequence of indices for the characters in seq"""
        if isinstance(seq, tuple):
            indices: list[int] = []
            for c in seq:
                indices.extend(self.to_indices(c).tolist())
            return numpy.array(indices, dtype=self.dtype)

        if isinstance(seq, str):
            seq = seq.encode("utf8")

        if isinstance(seq, bytes):
            # any non-canonical characters should lie outside the range
            # we replace these with a single value
            return self._bytes2arr(seq)

        if isinstance(seq, numpy.ndarray):
            return seq.astype(self.dtype)

        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    def from_indices(self, seq: str | bytes | NumpyIntArrayType) -> str:
        """returns a string from a sequence of indices"""
        if isinstance(seq, str):
            return seq
        if isinstance(seq, bytes):
            return seq.decode("utf8")
        if isinstance(seq, numpy.ndarray):
            return self._arr2bytes(seq).decode("utf8")

        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    def with_gap_motif(
        self,
        gap_char: str
        | bytes = "-",  # Strictly speaking this should use StrOrBytes but generics don't play well with these default arguments
        missing_char: str | bytes | None = "?",
        include_missing: bool = False,
        gap_as_state: bool = False,
    ) -> Self:
        """returns new monomer alphabet with gap and missing characters added

        Parameters
        ----------
        gap_char
            the IUPAC gap character "-"
        missing_char
            the IUPAC missing character "?"
        include_missing
            if True, and self.missing_char, it is included in the new alphabet
        gap_as_state
            include the gap character as a state in the alphabet, drops gap_char
            attribute in resulting KmerAlphabet
        """
        if self.gap_char and self.missing_char:
            return self

        missing_char = self._missing_char or missing_char if include_missing else None
        chars = cast(
            "tuple[StrOrBytes, ...]",
            (
                (*tuple(self[: self._num_canonical]), gap_char, missing_char)
                if missing_char
                else (*tuple(self[: self._num_canonical]), gap_char)
            ),
        )
        return self.__class__(
            chars,
            gap=None if gap_as_state else cast("StrOrBytes", gap_char),
            missing=cast("StrOrBytes", missing_char),
        )

    def get_kmer_alphabet(
        self, k: int, include_gap: bool = True
    ) -> "KmerAlphabet[StrOrBytes] | CharAlphabet[StrOrBytes]":
        """returns kmer alphabet with words of size k

        Parameters
        ----------
        k
            word size
        include_gap
            if True, and self.gap_char, we set
            KmerAlphabet.gap_char = self.gap_char * k

        Notes
        -----
        If self.missing_char is present, it is included in the new alphabet as
        missing_char * k
        """
        # refactor: revisit whether the include_gap argument makes sense
        if k == 1:
            return self

        return _get_kmer_alpha(
            self,
            self.num_canonical,
            self.gap_char,
            self.missing_char,
            include_gap,
            k,
            self.moltype,
        )

    def is_valid(self, seq: str | bytes | NumpyIntArrayType) -> bool:
        """seq is valid for alphabet"""
        if isinstance(seq, (str, bytes)):
            seq = self.to_indices(seq)

        if isinstance(seq, numpy.ndarray):
            return bool(seq.min() >= 0 and seq.max() < len(self) if len(seq) else True)

        # refactor: design
        if hasattr(seq, "alphabet"):
            # assume a SeqView instance
            return seq.alphabet == self
        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    def as_bytes(self) -> bytes:
        """returns self as a byte string"""
        if isinstance(self[0], bytes):
            return b"".join(cast("tuple[bytes, ...]", self))

        return "".join(cast("tuple[str, ...]", self)).encode("utf8")

    def array_to_bytes(self, seq: NumpyIntArrayType) -> bytes:
        """returns seq as a byte string"""
        return self._arr2bytes(seq)

    def to_rich_dict(self, for_pickle: bool = False) -> dict[str, Any]:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        data: dict[str, Any] = {
            "chars": list(self),
            "gap": self.gap_char,
            "missing": self.missing_char,
        }
        if not for_pickle:
            data["type"] = get_object_provenance(self)
            data["version"] = __version__

        return data

    def to_json(self) -> str:
        """returns a serialisable string"""
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, data: dict[str, Any]) -> "CharAlphabet[StrOrBytes]":
        """returns an instance from a serialised dictionary"""
        return cls(data["chars"], gap=data["gap"], missing=data["missing"])

    def get_subset(
        self,
        motif_subset: PySeq[StrOrBytes],
        excluded: bool = False,
    ) -> "CharAlphabet[StrOrBytes]":
        """Returns a new Alphabet object containing a subset of motifs in self.

        Raises an exception if any of the items in the subset are not already
        in self.
        """
        self_set = set(self)
        self_set.update(c for c in (self.gap_char, self.missing_char) if c is not None)

        if diff := set(motif_subset) - self_set:
            msg = f"{diff!r} not members in self"
            raise AlphabetError(msg)

        if excluded:
            motif_subset = tuple(m for m in self if m not in motif_subset)

        motif_subset = tuple(motif_subset)

        gap = self.gap_char if self.gap_char in motif_subset else None
        missing = self.missing_char if self.missing_char in motif_subset else None
        return make_alphabet(
            chars=motif_subset,
            gap=gap,
            missing=missing,
            moltype=cast("MolType[StrOrBytes]", self.moltype),
        )

    def convert_seq_array_to(
        self,
        *,
        alphabet: "CharAlphabet[Any]",
        seq: NumpyIntArrayType,
        check_valid: bool = True,
    ) -> NumpyIntArrayType:
        """converts a numpy array with indices from self to other

        Parameters
        ----------
        alphabet
            alphabet to convert to
        seq
            ndarray of uint8 integers
        check_valid
            validates both input and out sequences are valid for self and
            other respectively. Validation failure raises an AlphabetError.

        Returns
        -------
        the indices of characters in common between self and other
        are swapped
        """
        if check_valid and not self.is_valid(seq):
            msg = f"input sequence not valid for {self}"
            raise AlphabetError(msg)

        converter = make_converter(self, alphabet)
        output = numpy.frombuffer(converter(seq.tobytes()), dtype=alphabet.dtype)
        if check_valid and not alphabet.is_valid(output):
            msg = f"output sequence not valid for {alphabet}"
            raise AlphabetError(msg)

        return output


@functools.cache
def _get_kmer_alpha(
    monomers: CharAlphabet[StrOrBytes],
    num_canonical: int,
    gap_char: StrOrBytes | None,
    missing_char: StrOrBytes | None,
    include_gap: bool,
    k: int,
    moltype: "MolType[StrOrBytes]",
) -> "KmerAlphabet[StrOrBytes]":
    chars = tuple(monomers[:num_canonical])

    gap = gap_char * k if include_gap and gap_char is not None else None
    missing = missing_char * k if missing_char is not None and include_gap else None
    joiner: Callable[[Iterable[str | bytes]], StrOrBytes] = cast(
        "Callable[[Iterable[str | bytes]], StrOrBytes]",
        (b"".join if isinstance(monomers[0], bytes) else "".join),
    )
    words = tuple(joiner(e) for e in itertools.product(chars, repeat=k))
    words += (gap,) if include_gap and gap else ()
    words += (missing,) if include_gap and missing else ()

    kalpha = KmerAlphabet(words=words, monomers=monomers, gap=gap, k=k, missing=missing)
    if moltype:
        _alphabet_moltype_map[kalpha] = moltype
    return kalpha


@functools.cache
def make_converter(
    src_alpha: CharAlphabet[Any],
    dest_alpha: CharAlphabet[Any],
) -> convert_alphabet:
    """makes a convert_alphabet instance between two CharAlphabet alphabets

    Parameters
    ----------
    src_alpha, dest_alpha
        the two CharAlphabet instances

    Notes
    -----
    The common characters between the two alphabets are used to create the
    convert_alphabet() instance. The result of applying the convertor to a
    sequence should be validated by the caller to only shared characters
    were present.
    """
    # we need to make a converter that maps the indices of common characters
    # to new their values in dest_alpha and maps chars unique to src_alpha
    # to a single value that exceeds the length of dest_alpha so it's
    # easy to validate
    src_bytes = {bytes([c]) for c in src_alpha.as_bytes()}
    other_bytes = {bytes([c]) for c in dest_alpha.as_bytes()}
    common = b"".join(src_bytes & other_bytes)

    diff = b"".join(src_bytes - other_bytes)
    src_indices: list[NumpyIntArrayType] = [
        src_alpha.to_indices(common),
        src_alpha.to_indices(diff),
    ]
    dest_indices: list[NumpyIntArrayType] = [
        dest_alpha.to_indices(common),
        numpy.array([len(dest_alpha)] * len(diff), dtype=dest_alpha.dtype),
    ]

    if dest_alpha.moltype is None:
        msg = "Expected dest_alpha to have a moltype but is None."
        raise ValueError(msg)
    # if we are going to a nucleic acid alphabet, we need to add an alias
    # mapping for the character NOT present in that alphabet, i.e. if dest is
    # RNA, we need to map T to U, if dest is DNA we need to map U to T
    # we do this by adding these on at the end of the common characters
    if dest_alpha.moltype.is_nucleic and "T" in dest_alpha:
        src_indices.append(src_alpha.to_indices("U"))
        dest_indices.append(dest_alpha.to_indices("T"))
    elif dest_alpha.moltype.is_nucleic:
        src_indices.append(src_alpha.to_indices("T"))
        dest_indices.append(dest_alpha.to_indices("U"))

    return convert_alphabet(
        numpy.concatenate(src_indices).tobytes(),
        numpy.concatenate(dest_indices).tobytes(),
        allow_duplicates=True,
    )


def _get_closest_char_alphabet(
    moltype_name: str, motifset: list[str] | set[str]
) -> CharAlphabet[Any]:
    from cogent3.core.moltype import get_moltype

    mtype = get_moltype(moltype_name)
    motifset = set(motifset)
    num_motifs = len(motifset)
    alphas = sorted([(abs(len(a) - num_motifs), a) for a in mtype.iter_alphabets()])
    return alphas[0][1]


@register_deserialiser(
    get_object_provenance(CharAlphabet), "cogent3.core.alphabet.CharAlphabet"
)
def deserialise_char_alphabet(data: dict[str, Any]) -> CharAlphabet[Any]:
    if "motifset" not in data:
        return CharAlphabet.from_rich_dict(data)

    from cogent3.core.moltype import get_moltype

    # this is a legacy format, we assume this was derived from a
    # moltype so we get the moltype and look for an alphabet with
    # matching characters and gap
    mtype = get_moltype(data["moltype"])
    motifset = set(data.pop("motifset"))
    return _get_closest_char_alphabet(mtype, motifset)


@numba.jit(nopython=True)
def coord_conversion_coeffs(
    num_states: int,
    k: int,
    dtype: type[numpy.integer] | None = None,
) -> NumpyIntArrayType:  # pragma: no cover
    """coefficients for multi-dimensional coordinate conversion into 1D index"""
    return numpy.array([num_states ** (i - 1) for i in range(k, 0, -1)], dtype=dtype)


@numba.jit(nopython=True)
def coord_to_index(
    coord: NumpyIntArrayType,
    coeffs: NumpyIntArrayType,
) -> numpy.integer:  # pragma: no cover
    """converts a multi-dimensional coordinate into a 1D index"""
    return (coord * coeffs).sum()


@numba.jit(nopython=True)
def index_to_coord(
    index: int, coeffs: NumpyIntArrayType
) -> NumpyIntArrayType:  # pragma: no cover
    """converts a 1D index into a multi-dimensional coordinate"""
    ndim = len(coeffs)
    coord = numpy.zeros(ndim, dtype=coeffs.dtype)
    remainder = index
    for i in range(ndim):
        n, remainder = numpy.divmod(remainder, coeffs[i])
        coord[i] = n
    return coord


@numba.jit(nopython=True)
def seq_to_kmer_indices(
    seq: NumpyIntArrayType,
    result: NumpyIntArrayType,
    coeffs: NumpyIntArrayType,
    num_states: int,
    k: int,
    gap_char_index: int | numpy.integer = -1,
    gap_index: int = -1,
    independent_kmer: bool = True,
) -> NumpyIntArrayType:  # pragma: no cover
    """return 1D indices for valid k-mers

    Parameters
    ----------
    seq
        numpy array of uint, assumed that canonical characters have
        sequential indexes which are all < num_states
    result
        array to be written into
    coeffs
        result of calling coord_conversion_coeffs() with num_states and k
    num_states
        the number of canonical states, which defines range of allowed
        ints at a position
    k
        k-mer size
    gap_char_index
        value of a gap character in seq. If > 0, any k-mer containing
        a gap character is assigned an index of num_states**k. Any
        k-mer containing an integer greater than gap_index is assigned an
        index of (num_states**k) + 1. If gap_char_index < 0, then any
        k-mer containing an integer greater than num_states is assigned an
        index of num_states**k.
    gap_index
        value to be assigned when gap_char_index is != -1, must be positive
    independent_kmer
        if True, sets returns indices for non-overlapping k-mers, otherwise
        returns indices for all k-mers
    """
    # check the result length is consistent with the settings
    step: int = k if independent_kmer else 1
    size: int = int(numpy.ceil((len(seq) - k + 1) / step))
    if len(result) < size:
        msg = f"size of result {len(result)} <= {size}"
        raise ValueError(msg)

    missing_index: int = gap_index + 1 if gap_char_index > 0 else num_states**k
    if gap_char_index > 0 and gap_index <= 0:
        msg = f"gap_index={gap_index} but gap_char_index={gap_char_index}"
        raise ValueError(msg)

    for result_idx, i in enumerate(range(0, len(seq) - k + 1, step)):
        segment = seq[i : i + k]
        if (num_states > segment).all():
            result[result_idx] = coord_to_index(seq[i : i + k], coeffs)
        elif gap_char_index > 0 and segment.max() == gap_char_index:
            result[result_idx] = gap_index
        else:
            result[result_idx] = missing_index
    return result


# TODO: profile this against pure python
@numba.jit(nopython=True)
def kmer_indices_to_seq(
    kmer_indices: NumpyIntArrayType,
    result: NumpyIntArrayType,
    coeffs: NumpyIntArrayType,
    k: int,
    gap_index: int | numpy.integer = -1,
    gap_char_index: int | numpy.integer = -1,
    independent_kmer: bool = True,
) -> NumpyIntArrayType:  # pragma: no cover
    if gap_index > 0:
        if gap_char_index <= 0:
            msg = f"gap_index={gap_index} but gap_char_index={gap_char_index}"
            raise ValueError(msg)
        gap_state = numpy.array([gap_char_index] * k, dtype=result.dtype)
    else:
        gap_state = numpy.empty(0, dtype=result.dtype)
    for index, kmer_index in enumerate(kmer_indices):
        if kmer_index == gap_index:
            coord = gap_state
        else:
            coord = index_to_coord(kmer_index, coeffs)

        if index == 0:
            result[:k] = coord
        elif independent_kmer:
            seq_index = index * k
            result[seq_index : seq_index + k] = coord
        else:
            # we are just adding the last monomer
            result[index + k - 1] = coord[-1]

    return result


class KmerAlphabetABC(ABC, Generic[StrOrBytes]):
    @abstractmethod
    def to_index(self, seq: str | bytes | NumpyIntArrayType) -> int: ...

    @abstractmethod
    def from_index(self, kmer_index: int) -> str | NumpyIntArrayType: ...

    @abstractmethod
    def with_gap_motif(
        self,
        include_missing: bool = False,
        **kwargs: Any,
    ) -> "KmerAlphabetABC[StrOrBytes]": ...


class KmerAlphabet(
    tuple[StrOrBytes, ...],
    Generic[StrOrBytes],
    AlphabetABC[StrOrBytes],
    KmerAlphabetABC[StrOrBytes],
):
    """k-mer alphabet represents complete non-monomer alphabets

    Notes
    -----
    Differs from SenseCodonAlphabet case by representing all possible permutations of
    k-length of the provided monomer alphabet. More efficient mapping between integers
    and k-length strings
    """

    __slots__ = ()

    def __new__(
        cls,
        words: tuple[StrOrBytes, ...],
        monomers: CharAlphabet[StrOrBytes],
        k: int,
        gap: StrOrBytes | None = None,
        missing: StrOrBytes | None = None,
    ) -> "KmerAlphabet[StrOrBytes]":
        """
        Parameters
        ----------
        words
            len(monomers)**k k-length strings
        monomers
            the underlying character alphabet
        k
            the number of monomers per word
        gap
            k monomer gap characters
        missing, optional
            k monomer missing characters
        """
        if not words:
            msg = f"cannot create empty {cls.__name__!r}"
            raise ValueError(msg)

        if gap is not None and _coerce_to_type(words[0], gap) not in words:
            msg = f"gap {gap!r} not in words {words!r}"
            raise ValueError(msg)

        consistent_words(words, length=k)
        return tuple.__new__(cls, words)

    def __init__(
        self,
        words: tuple[StrOrBytes, ...],
        monomers: CharAlphabet[StrOrBytes],
        k: int,
        gap: StrOrBytes | None = None,
        missing: StrOrBytes | None = None,
    ) -> None:
        """
        Parameters
        ----------
        words
            series of k length strings
        monomers
            base CharAlphabet instance
        k
            word size
        gap
            the gap state ("-" * k) if present
        """
        self.monomers: CharAlphabet[StrOrBytes] = monomers
        self.dtype = get_array_type(len(self))
        self._words: set[StrOrBytes] = set(self)  # for quick lookup
        self.k = k

        self._gap_char: StrOrBytes | None = gap
        self._gap_index = self.index(gap) if gap else None
        self._missing_char: StrOrBytes | None = missing
        self._missing_index = self.index(missing) if missing else None
        if gap and (
            self.monomers.gap_char is None
            or gap != self.monomers.gap_char * k
            or gap not in self
        ):
            msg = f"Unexpected gap {gap!r}"
            raise ValueError(msg)

        if missing and (
            self.monomers.missing_char is None
            or missing != self.monomers.missing_char * k
            or missing not in self
        ):
            msg = f"Unexpected missing char {missing!r}"
            raise ValueError(msg)
        self._coeffs = coord_conversion_coeffs(
            self.monomers.num_canonical,
            k,
            dtype=self.dtype,
        )
        self._num_canonical = self.monomers.num_canonical**k

    def __reduce__(
        self,
    ) -> tuple[
        type,
        tuple[
            tuple[StrOrBytes, ...],
            CharAlphabet[StrOrBytes],
            int,
            StrOrBytes | None,
            StrOrBytes | None,
        ],
        None,
    ]:
        """support for pickling"""
        return (
            self.__class__,
            (
                tuple(self),
                self.monomers,
                self.k,
                self.gap_char,
                self.missing_char,
            ),
            None,
        )

    @property
    def gap_char(self) -> StrOrBytes | None:
        return self._gap_char

    @property
    def gap_index(self) -> int | None:
        return self._gap_index

    @property
    def missing_char(self) -> StrOrBytes | None:
        return self._missing_char

    @property
    def missing_index(self) -> int | None:
        return self._missing_index

    @property
    def num_canonical(self) -> int:
        return self._num_canonical

    def to_indices(
        self,
        seq: str
        | bytes
        | tuple[str | bytes, ...]
        | list[str | bytes]
        | NumpyIntArrayType,
        independent_kmer: bool = True,
    ) -> NumpyIntArrayType:
        """returns a sequence of k-mer indices

        Parameters
        ----------
        seq
            a sequence of monomers
        independent_kmer
            if True, returns non-overlapping k-mers

        Notes
        -----
        If self.gap_char is not None, then the following rules apply:
        If a sequence k-mer contains a gap character it is
        assigned an index of (num. monomer states**k). If
        a k-mer contains a non-canonical and non-gap character,
        it is assigned an index of (num. monomer states**k) + 1.
        If self.gap_char is None, then both of the above cases
        are defined as (num. monomer states**k)."""
        # TODO: handle case of non-modulo sequences
        if isinstance(seq, (tuple, list)):
            return numpy.array([self.to_index(c) for c in seq], dtype=self.dtype)

        if isinstance(seq, (str, bytes)):
            seq = self.monomers.to_indices(seq)

        if isinstance(seq, numpy.ndarray):
            size = len(seq) - self.k + 1
            if independent_kmer:
                size = int(numpy.ceil(size / self.k))
            result = numpy.zeros(size, dtype=self.dtype)
            gap_char_index = self.monomers.gap_index or -1
            gap_index = self.gap_index or -1
            return seq_to_kmer_indices(
                seq,
                result,
                self._coeffs,
                self.monomers.num_canonical,
                self.k,
                gap_char_index=gap_char_index,
                gap_index=gap_index,
                independent_kmer=independent_kmer,
            )

        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    def from_indices(  # type: ignore[override]
        self,
        kmer_indices: NumpyIntArrayType,
        independent_kmer: bool = True,
    ) -> NumpyIntArrayType:
        """converts array of k-mer indices into an array of monomer indices

        Parameters
        ----------
        kmer_indices
            a sequence of k-mer indices
        independent_kmer
            whether the k-mers are overlapping or not
            _description_
        """
        if independent_kmer:
            size = len(kmer_indices) * self.k
        else:
            size = len(kmer_indices) + self.k - 1

        result = numpy.zeros(size, dtype=self.dtype)
        gap_char_index = self.monomers.gap_index or -1
        gap_index = self.gap_index or -1
        return kmer_indices_to_seq(
            kmer_indices,
            result,
            self._coeffs,
            self.k,
            gap_char_index=gap_char_index,
            gap_index=gap_index,
            independent_kmer=independent_kmer,
        )

    def with_gap_motif(  # type: ignore[override]
        self,
        include_missing: bool = False,
        **kwargs: Any,
    ) -> "KmerAlphabet[StrOrBytes] | CharAlphabet[StrOrBytes]":
        """returns a new KmerAlphabet with the gap motif added

        Notes
        -----
        Adds gap state to monomers and recreates k-mer alphabet
        for self

        kwargs is for compatibility with the CharAlphabet method
        """
        if self.gap_char:
            return self
        missing = "?" if include_missing else None
        monomers = self.monomers.with_gap_motif(missing_char=missing)
        return monomers.get_kmer_alphabet(self.k, include_gap=True)

    def to_index(self, seq: str | bytes | NumpyIntArrayType) -> int:
        """encodes a k-mer as a single integer

        Parameters
        ----------
        seq
            sequence to be encoded, can be either a string or numpy array
        overlapping
            if False, performs operation on sequential k-mers, e.g. codons


        Notes
        -----
        If self.gap_char is defined, then the following rules apply:
        returns num_states**k if a k-mer contains a gap character,
        otherwise returns num_states**k + 1 if a k-mer contains a
        non-canonical character. If self.gap_char is not defined,
        returns num_states**k for both cases."""
        if isinstance(seq, (str, bytes)):
            seq = self.monomers.to_indices(seq)

        if isinstance(seq, numpy.ndarray):
            if len(seq) != self.k:
                msg = f"Length of seq ({len(seq)}) is not k={self.k}"
                raise ValueError(msg)
            if (self.monomers.num_canonical > seq).all():
                return int(coord_to_index(seq, self._coeffs))

            gap_char_index = self.monomers.gap_index or self.monomers.num_canonical
            gap_index = self.gap_index or self.num_canonical
            missing_index = gap_index + 1 if self.gap_char else gap_index
            return gap_index if seq.max() == gap_char_index else missing_index

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def from_index(self, kmer_index: int) -> NumpyIntArrayType:
        """decodes an integer into a k-mer"""
        if self.gap_index is not None and kmer_index == self.gap_index:
            return numpy.array([self.monomers.gap_index] * self.k, dtype=self.dtype)
        if self.missing_index is not None and kmer_index == self.missing_index:
            return numpy.array([self.monomers.missing_index] * self.k, dtype=self.dtype)
        if kmer_index >= len(self):
            msg = f"{kmer_index} is out of range"
            raise ValueError(msg)

        return index_to_coord(kmer_index, self._coeffs)

    def is_valid(self, seq: str | bytes | NumpyIntArrayType) -> bool:
        """seq is valid for alphabet

        Parameters
        ----------
        seq
            a numpy array of integers

        Notes
        -----
        This will raise a TypeError for string or bytes. Using to_indices() to
        convert those ensures a valid result.
        """
        if isinstance(seq, numpy.ndarray):
            max_val = max(
                self.missing_index or 0, self.gap_index or 0, self.num_canonical
            )
            return bool(seq.min() >= 0 and seq.max() <= max_val if len(seq) else True)

        msg = f"{type(seq)} is invalid, must be numpy.ndarray"
        raise TypeError(msg)

    def to_rich_dict(self, for_pickle: bool = False) -> dict[str, Any]:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        data: dict[str, Any] = {
            "words": list(self),
            "monomers": self.monomers.to_rich_dict(),
            "k": self.k,
            "gap": self.gap_char,
            "missing": self.missing_char,
        }
        if not for_pickle:
            data["type"] = get_object_provenance(self)
            data["version"] = __version__
        return data

    def to_json(self) -> str:
        """returns a serialisable string"""
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, data: dict[str, Any]) -> "KmerAlphabet[StrOrBytes]":
        """returns an instance from a serialised dictionary"""
        monomers = CharAlphabet[StrOrBytes].from_rich_dict(data["monomers"])
        return cls(
            words=data["words"],
            monomers=monomers,
            k=data["k"],
            gap=data["gap"],
            missing=data["missing"],
        )

    @property
    def motif_len(self) -> int:
        return self.k


@register_deserialiser(
    get_object_provenance(KmerAlphabet), "cogent3.core.alphabet.JointEnumeration"
)
def deserialise_kmer_alphabet(
    data: dict[str, Any],
) -> KmerAlphabet[Any] | CharAlphabet[Any]:
    if data["type"] != "cogent3.core.alphabet.JointEnumeration":
        return KmerAlphabet.from_rich_dict(data)

    monomers = data["data"][0]
    motif_length = len(data["data"])
    # which is the closest char alphabet
    mtype = data["moltype"]
    alpha = _get_closest_char_alphabet(mtype, monomers)
    return alpha.get_kmer_alphabet(motif_length, include_gap=bool(data["gap"]))


class SenseCodonAlphabet(
    tuple[str, ...],
    AlphabetABC[str],
    KmerAlphabetABC[str],
):
    """represents the sense-codons of a GeneticCode"""

    __slots__ = ()

    def __new__(
        cls,
        words: tuple[str, ...],
        monomers: CharAlphabet[str],
        gap: str | None = None,
    ) -> "SenseCodonAlphabet":
        if not words:
            msg = f"cannot create empty {cls.__name__!r}"
            raise ValueError(msg)

        if gap is not None:
            assert _coerce_to_type(words[0], gap) in words

        consistent_words(words, length=3)
        return tuple.__new__(cls, words)

    def __init__(
        self,
        words: tuple[str, ...],  # noqa: ARG002
        monomers: CharAlphabet[str],
        gap: str | None = None,
    ) -> None:
        """
        Parameters
        ----------
        words
            series of 3 character strings representing sense codons
        monomers
            CharAlphabet instance of DNA monomers, including degenerate
            characters
        gap
            the gap state "---" if present
        """
        self._gap_char = gap
        self.monomers = monomers
        self.dtype = get_array_type(len(self))
        self._words = set(self)  # for quick lookup
        self._to_indices = {codon: i for i, codon in enumerate(self)}
        self._from_indices = {i: codon for codon, i in self._to_indices.items()}
        self._motif_len = 3
        if monomers.moltype:
            _alphabet_moltype_map[self] = monomers.moltype

    def __reduce__(
        self,
    ) -> tuple[
        type,
        tuple[tuple[str, ...], CharAlphabet[str], str | None],
        None,
    ]:
        """support for pickling"""
        return (
            self.__class__,
            (tuple(self), self.monomers, self.gap_char),
            None,
        )

    @property
    def gap_char(self) -> str | None:
        return self._gap_char

    @property
    def gap_index(self) -> int | None:
        return self._to_indices[self.gap_char] if self.gap_char else None

    @property
    def missing_char(self) -> None:
        """not supported on CodonAlphabet"""
        return None

    @property
    def missing_index(self) -> None:
        """not supported on CodonAlphabet"""
        return None

    def to_indices(  # type: ignore[override]
        self, seq: str | list[str] | tuple[str, ...] | NumpyIntArrayType
    ) -> NumpyIntArrayType:
        """returns a sequence of codon indices"""

        if isinstance(seq, str):
            seq = [seq[i : i + 3] for i in range(0, len(seq), 3)]

        if isinstance(seq, numpy.ndarray):
            size = len(seq) // 3
            seq = [self.monomers.from_indices(c) for c in seq.reshape(size, 3)]

        if isinstance(seq, (list, tuple)):
            return numpy.array([self.to_index(c) for c in seq], dtype=self.dtype)

        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    def to_index(self, codon: str) -> int:  # type: ignore
        """encodes a codon as a single integer"""
        if len(codon) != 3:
            msg = f"{codon=!r} is not of length 3"
            raise ValueError(msg)
        if not self.monomers.is_valid(codon):
            msg = f"{codon=!r} elements not nucleotides"
            raise AlphabetError(msg)

        if self.moltype and self.moltype.has_ambiguity(codon):
            return len(self)

        try:
            return self._to_indices[codon]
        except KeyError as e:
            msg = f"{codon=!r} not in alphabet"
            raise AlphabetError(msg) from e

    def from_index(self, index: int) -> str:  # type: ignore[override]
        """returns a single codon from an index"""
        if index > len(self) or index < 0:
            msg = f"{index=!r} is not within range"
            raise ValueError(msg)

        try:
            return self._from_indices[index]
        except KeyError as e:
            msg = f"invalid {index=}"
            raise ValueError(msg) from e

    def from_indices(self, indices: NumpyIntArrayType) -> list[str]:  # type: ignore[override]
        """returns a list of codons from a numpy array of indices"""
        return [self.from_index(index) for index in indices]

    @property
    def num_canonical(self) -> int:
        """returns the number of canonical states"""
        return len(self._words)

    def is_valid(self, seq: str | bytes | NumpyIntArrayType) -> bool:
        """seq is valid for alphabet"""

        if isinstance(seq, numpy.ndarray):
            return bool(seq.min() >= 0 and seq.max() < len(self) if len(seq) else True)

        if isinstance(seq, str):
            try:
                _ = self.to_indices(seq)
            except (ValueError, AlphabetError):
                return False
            else:
                return True

        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    def with_gap_motif(self, **kwargs: Any) -> Self:  # type: ignore[override]
        """returns a new SenseCodonAlphabet with the gap motif '---' added

        Notes
        -----
        kwargs is for compatibility with the ABC method
        """
        if self.gap_char:
            return self
        monomers = self.monomers.with_gap_motif()
        if not monomers.gap_char:
            msg = f"{monomers.gap_char!r} not a valid gap character"
            raise ValueError(msg)

        gap_char = monomers.gap_char * 3
        words = (*tuple(self), gap_char)
        return self.__class__(words=words, monomers=monomers, gap=gap_char)

    def to_rich_dict(self, for_pickle: bool = False) -> dict[str, Any]:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        data: dict[str, Any] = {
            "monomers": self.monomers.to_rich_dict(),
            "words": list(self),
        }
        if not for_pickle:
            data["type"] = get_object_provenance(self)
            data["version"] = __version__
        return data

    def to_json(self) -> str:
        """returns a serialisable string"""
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, data: dict[str, Any]) -> "SenseCodonAlphabet":
        """returns an instance from a serialised dictionary"""
        data["monomers"] = deserialise_char_alphabet(data["monomers"])
        data.pop("type", None)
        data.pop("version", None)
        return cls(**data)

    @property
    def motif_len(self) -> int:
        """always 3 for codon alphabets"""
        return self._motif_len


@register_deserialiser("cogent3.core.alphabet.Alphabet")
def deserialise_old_alphabet(
    data: dict[str, Any],
) -> KmerAlphabet[Any] | CharAlphabet[Any] | SenseCodonAlphabet:
    motif_length = len(data["motifset"][0])
    if motif_length == 1:
        # this is a char alphabet
        return deserialise_char_alphabet(data)

    motifs = set(data["motifset"])
    num_motifs = len(motifs)
    num_monomers = round(num_motifs ** (1 / motif_length))
    if num_monomers**motif_length != num_motifs:
        # if the number of motifs is not the motif_length exponent of the number
        # of monomers, we have a codon alphabet
        return deserialise_codon_alphabet(data)

    # we get the closest alphabet
    mtype = data["moltype"]
    monomers = list(set("".join(motifs)))
    alpha = _get_closest_char_alphabet(mtype, monomers)
    return alpha.get_kmer_alphabet(motif_length, include_gap=bool(data["gap"]))


@register_deserialiser(get_object_provenance(SenseCodonAlphabet))
def deserialise_codon_alphabet(data: dict[str, Any]) -> SenseCodonAlphabet:
    if data["type"] != "cogent3.core.alphabet.Alphabet":
        return SenseCodonAlphabet.from_rich_dict(data)

    # we search through the genetic codes for the one with a matching
    # motifset
    from cogent3.core.genetic_code import GeneticCode, available_codes, get_code

    include_gap = bool(data["gap"])
    motifs = set(data["motifset"])
    alphabets: list[tuple[int, GeneticCode]] = []
    code_ids = available_codes().columns["Code ID"]
    for code_id in sorted(code_ids):
        code = get_code(code_id)
        code_alpha = code.get_alphabet(include_gap=include_gap)
        diff = len(set(code_alpha) ^ motifs)
        alphabets.append((diff, code))

    alphabets = sorted(alphabets, key=lambda x: x[0])
    code = alphabets[0][1]
    return code.get_alphabet(include_gap=include_gap)


_alphabet_moltype_map: dict[AlphabetABC[Any], "MolType[Any]"] = {}


def make_alphabet(
    *,
    chars: PySeq[StrOrBytes] | StrOrBytes,
    gap: StrOrBytes | None,
    missing: StrOrBytes | None,
    moltype: "MolType[StrOrBytes]",
) -> CharAlphabet[StrOrBytes]:
    """constructs a character alphabet and registers the associated moltype

    Notes
    -----
    The moltype is associated with the alphabet, available as the
    alphabet.moltype property.
    If an alphabet has already been constructed by a moltype no new
    entry is made. An alphabet shared by multiple moltype instances
    will therefore only return the first moltype associated with it.
    """
    alpha = CharAlphabet(chars, gap=gap, missing=missing)
    if alpha not in _alphabet_moltype_map:
        _alphabet_moltype_map[alpha] = moltype

    return alpha
