import functools
import itertools
import json
import typing

from abc import ABC, abstractmethod

import numba
import numpy

from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import get_object_provenance


StrORBytes = typing.Union[str, bytes]
StrORBytesORArray = typing.Union[str, bytes, numpy.ndarray]
OptInt = typing.Optional[int]
OptStr = typing.Optional[str]
OptBytes = typing.Optional[bytes]


def _element_type(val: typing.Sequence) -> typing.Any:
    return type(val[0])


@functools.singledispatch
def _coerce_to_type(orig: StrORBytes, text: str) -> StrORBytes:
    raise TypeError(f"{type(orig)} is invalid")


@_coerce_to_type.register
def _(orig: str, text: str) -> StrORBytes:
    return text


@_coerce_to_type.register
def _(orig: bytes, text: str) -> StrORBytes:
    return text.encode("utf8")


@_coerce_to_type.register
def _(orig: bytearray, text: str) -> StrORBytes:
    return text.encode("utf8")


@_coerce_to_type.register
def _(orig: tuple, text: str) -> StrORBytes:
    return (text.encode("utf8"),)


@_coerce_to_type.register
def _(orig: list, text: str) -> StrORBytes:
    return (text.encode("utf8"),)


class AlphabetError(TypeError): ...


class convert_alphabet:
    """convert one character set into another"""

    def __init__(self, src: bytes, dest: bytes, delete: OptBytes = None):
        """
        Parameters
        ----------
        src
            unique series of characters
        new
            characters in src will be mapped to these in order
        delete
            characters to be deleted from the result

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
            raise ValueError(f"length of src={len(src)} != length of new {len(dest)}")

        consistent_words(src, length=1)

        self._table = b"".maketrans(src, dest)
        self._delete = delete or b""

    def __call__(self, seq: bytes) -> bytes:
        return seq.translate(self._table, delete=self._delete)


class AlphabetABC(ABC):
    @property
    @abstractmethod
    def num_canonical(self) -> int: ...

    @property
    @abstractmethod
    def motif_len(self) -> int: ...

    @abstractmethod
    def is_valid(self, seq) -> bool: ...

    @abstractmethod
    def with_gap_motif(self): ...

    @abstractmethod
    def to_indices(self, seq: str) -> numpy.ndarray: ...

    @abstractmethod
    def from_indices(self, seq: numpy.ndarray) -> str: ...

    @property
    @abstractmethod
    def gap_char(self) -> OptStr: ...

    @property
    @abstractmethod
    def gap_index(self) -> OptInt: ...

    @property
    @abstractmethod
    def missing_char(self) -> OptStr: ...

    @property
    @abstractmethod
    def missing_index(self) -> OptInt: ...

    @abstractmethod
    def to_rich_dict(self) -> dict: ...

    @abstractmethod
    def to_json(self) -> str: ...

    @classmethod
    @abstractmethod
    def from_rich_dict(cls, data: dict) -> None: ...


class MonomerAlphabetABC(ABC):
    @abstractmethod
    def get_kmer_alphabet(self, size: int): ...

    @abstractmethod
    def as_bytes(self) -> bytes: ...


def get_array_type(num_elements: int):
    """Returns the smallest numpy integer dtype that can contain elements
    within num_elements.
    """
    if num_elements < 2**8:
        dtype = numpy.uint8
    elif num_elements < 2**16:
        dtype = numpy.uint16
    elif num_elements < 2**32:
        dtype = numpy.uint32
    elif num_elements < 2**64:
        dtype = numpy.uint64
    else:
        raise NotImplementedError(f"{num_elements} is too big for 64-bit integer")

    return dtype


def consistent_words(
    words: typing.Sequence[typing.Union[str, int]], length: OptInt = None
) -> None:
    """make sure all alphabet elements are unique and have the same length"""
    if not words:
        raise ValueError("no words provided")

    if len(set(words)) != len(words):
        raise ValueError(f"duplicate elements in {words}")

    lengths = {1} if isinstance(words[0], int) else {len(w) for w in words}
    if len(lengths) != 1:
        raise ValueError(f"mixed lengths {lengths} in {words}")
    l = list(lengths)[0]
    if length and l != length:
        raise ValueError(f"word length {l} does not match expected {length}")


class bytes_to_array:
    """wrapper around convert_alphabet. It defines a linear mapping from provided
    characters to uint8. The resulting object is callable, taking a bytes object
    and returning a numpy array."""

    def __init__(self, chars: bytes, dtype, delete: OptBytes = None):
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
            chars, bytes(bytearray(range(len(chars)))), delete=delete
        )
        self.dtype = dtype

    def __call__(self, seq: bytes) -> numpy.ndarray:
        b = self._converter(seq)
        return numpy.array(memoryview(b), dtype=self.dtype)


class array_to_bytes:
    """wrapper around convert_alphabet. It defines a linear mapping from uint8
    integers to the provided characters. The resulting object is callable,
    taking a numpy array and returning a bytes object.
    """

    def __init__(self, chars: bytes):
        # we want a bytes translation map
        self._converter = convert_alphabet(
            bytes(bytearray(range(len(chars)))),
            chars,
        )

    def __call__(self, seq: numpy.ndarray) -> bytes:
        b = self._converter(seq.tobytes())
        return bytes(bytearray(b))


class CharAlphabet(tuple, AlphabetABC, MonomerAlphabetABC):
    def __new__(
        cls,
        chars: typing.Sequence[StrORBytes],
        gap: OptStr = None,
        missing: OptStr = None,
    ):
        if not chars:
            raise ValueError(f"cannot create empty {cls.__name__!r}")

        if gap is not None:
            assert _coerce_to_type(chars[0], gap) in chars

        if missing is not None:
            # if a bytes alphabet and missing provided this will fail
            # since individual elements of a bytes object are integers
            # leaving as is because we're considering a full bytes
            # alphabet only, in which case the missing character is already
            # present
            assert _coerce_to_type(chars[0], missing) in chars

        consistent_words(chars, length=1)
        return tuple.__new__(cls, chars, gap=gap, missing=missing)

    def __init__(
        self,
        chars: typing.Sequence[StrORBytes],
        gap: OptStr = None,
        missing: OptStr = None,
    ):
        self._gap_char = gap
        self._gap_index = self.index(gap) if gap else None
        self._missing_char = missing
        self._missing_index = self.index(missing) if missing else None

        # the number of canonical states are non-gap, non-missing states
        adj = (1 if gap else 0) + (1 if missing else 0)

        self._num_canonical = len(self) - adj

        self.dtype = get_array_type(len(self))
        self._chars = set(self)  # for quick lookup
        byte_chars = self.as_bytes()
        self._bytes2arr = bytes_to_array(byte_chars, dtype=self.dtype)
        self._arr2bytes = array_to_bytes(byte_chars)

    @property
    def num_canonical(self) -> int:
        return self._num_canonical

    @property
    def gap_char(self) -> OptStr:
        return self._gap_char

    @property
    def gap_index(self) -> OptInt:
        return self._gap_index

    @property
    def missing_char(self) -> OptStr:
        return self._missing_char

    @property
    def missing_index(self) -> OptInt:
        return self._missing_index

    @property
    def motif_len(self) -> int:
        # this is enforced by consistent_words(length=1) in the constructor
        return 1

    @functools.singledispatchmethod
    def to_indices(self, seq: StrORBytesORArray) -> numpy.ndarray:
        raise TypeError(f"{type(seq)} is invalid")

    @to_indices.register
    def _(self, seq: bytes) -> numpy.ndarray:
        # any non-canonical characters should lie outside the range
        # we replace these with a single value
        return self._bytes2arr(seq)

    @to_indices.register
    def _(self, seq: str) -> numpy.ndarray:
        return self.to_indices(seq.encode("utf8"))

    @to_indices.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        return seq

    @functools.singledispatchmethod
    def from_indices(self, seq: numpy.ndarray) -> str:
        raise TypeError(f"{type(seq)} is invalid")

    @from_indices.register
    def _(self, seq: numpy.ndarray) -> str:
        return self._arr2bytes(seq).decode("utf8")

    def with_gap_motif(self, gap_char="-", missing_char="?"):
        """returns new monomer alphabet with gap and missing characters added

        Parameters
        ----------
        gap_char
            the IUPAC gap character "-"
        missing_char
            the IUPAC missing character "?"
        """
        if self.gap_char and self.missing_char:
            return self
        chars = (
            tuple(self[: self._num_canonical]) + (gap_char, missing_char)
            if missing_char
            else (gap_char,)
        )
        return self.__class__(chars, gap=gap_char, missing=missing_char)

    @functools.cache
    def get_kmer_alphabet(self, k: int, include_gap: bool = True) -> "KmerAlphabet":
        """returns kmer alphabet with words of size k

        Parameters
        ----------
        k
            word size
        include_gap
            if True, and self.gap_char, we set
            KmerAlphabet.gap_char = self.gap_char * k
        """
        # refactor: revisit whether the include_gap argument makes sense
        if k == 1:
            return self

        chars = tuple(self[: self.num_canonical])

        gap = self._gap_char * k if include_gap and self._gap_char is not None else None
        missing = (
            self._missing_char * k
            if self._missing_char is not None and include_gap
            else None
        )
        words = tuple("".join(e) for e in itertools.product(chars, repeat=k))
        words += (gap,) if include_gap and gap else ()
        words += (missing,) if include_gap and missing else ()

        return KmerAlphabet(words=words, monomers=self, gap=gap, k=k, missing=missing)

    @functools.singledispatchmethod
    def is_valid(self, seq: StrORBytesORArray) -> bool:
        # refactor: design
        if hasattr(seq, "alphabet"):
            # assume a SeqView instance
            return seq.alphabet == self
        raise TypeError(f"{type(seq)} is invalid")

    @is_valid.register
    def _(self, seq: str) -> bool:
        return self.is_valid(self.to_indices(seq))

    @is_valid.register
    def _(self, seq: bytes) -> bool:
        return self.is_valid(self.to_indices(seq))

    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        return seq.min() >= 0 and seq.max() < len(self) if len(seq) else True

    def as_bytes(self) -> bytes:
        """returns self as a byte string"""
        if not isinstance(self[0], str):
            return bytes(bytearray(self))

        return "".join(self).encode("utf8")

    def array_to_bytes(self, seq: numpy.ndarray) -> bytes:
        """returns seq as a byte string"""
        return self._arr2bytes(seq)

    def to_rich_dict(self) -> dict:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        return {
            "chars": list(self),
            "gap": self.gap_char,
            "missing": self.missing_char,
            "version": __version__,
            "type": get_object_provenance(self),
        }

    def to_json(self) -> str:
        """returns a serialisable string"""
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, data: dict) -> "CharAlphabet":
        """returns an instance from a serialised dictionary"""
        return cls(data["chars"], gap=data["gap"], missing=data["missing"])


@register_deserialiser(get_object_provenance(CharAlphabet))
def deserialise_char_alphabet(data: dict) -> CharAlphabet:
    return CharAlphabet.from_rich_dict(data)


@numba.jit(nopython=True)
def coord_conversion_coeffs(num_states, k, dtype=None):  # pragma: no cover
    """coefficients for multi-dimensional coordinate conversion into 1D index"""
    return numpy.array([num_states ** (i - 1) for i in range(k, 0, -1)], dtype=dtype)


@numba.jit(nopython=True)
def coord_to_index(coord, coeffs):  # pragma: no cover
    """converts a multi-dimensional coordinate into a 1D index"""
    return (coord * coeffs).sum()


@numba.jit(nopython=True)
def index_to_coord(index: int, coeffs: numpy.ndarray):  # pragma: no cover
    """converts a 1D index into a multi-dimensional coordinate"""
    ndim = len(coeffs)
    coord = numpy.zeros(ndim, dtype=coeffs.dtype)
    remainder = index
    for i in range(ndim):
        n, remainder = numpy.divmod(remainder, coeffs[i])
        coord[i] = n
    return coord


@numba.jit()
def seq_to_kmer_indices(
    seq: numpy.ndarray,
    result: numpy.ndarray,
    coeffs: numpy.ndarray,
    num_states: int,
    k: int,
    gap_char_index: int = -1,
    gap_index: int = -1,
    independent_kmer: bool = True,
) -> numpy.ndarray:  # pragma: no cover
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
        raise ValueError(f"size of result {len(result)} <= {size}")

    missing_index: int = gap_index + 1 if gap_char_index > 0 else num_states**k
    if gap_char_index > 0:
        assert (
            gap_index > 0
        ), f"gap_index={gap_index} but gap_char_index={gap_char_index}"

    for result_idx, i in enumerate(range(0, len(seq) - k + 1, step)):
        segment = seq[i : i + k]
        if (num_states > segment).all():
            result[result_idx] = coord_to_index(seq[i : i + k], coeffs)
        elif gap_char_index > 0 and segment.max() == gap_char_index:
            result[result_idx] = gap_index
        else:
            result[result_idx] = missing_index
    return result


# todo: profile this against pure python
@numba.jit()
def kmer_indices_to_seq(
    kmer_indices: numpy.ndarray,
    result: numpy.ndarray,
    coeffs: numpy.ndarray,
    k: int,
    gap_index: int = -1,
    gap_char_index: int = -1,
    independent_kmer: bool = True,
) -> numpy.ndarray:  # pragma: no cover
    if gap_index > 0:
        assert (
            gap_char_index > 0
        ), f"gap_index={gap_index} but gap_char_index={gap_char_index}"
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


class KmerAlphabetABC(ABC):
    @abstractmethod
    def to_index(self, seq) -> int: ...

    @abstractmethod
    def from_index(self, kmer_index: int) -> str: ...


class KmerAlphabet(tuple, AlphabetABC, KmerAlphabetABC):
    def __new__(
        cls,
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        k: int,
        gap: OptStr = None,
        missing: OptStr = None,
    ):
        if not words:
            raise ValueError(f"cannot create empty {cls.__name__!r}")

        if gap is not None:
            assert _coerce_to_type(words[0], gap) in words

        consistent_words(words, length=k)
        return tuple.__new__(cls, words, monomers=monomers, gap=gap, missing=missing)

    def __init__(
        self,
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        k: int,
        gap: OptStr = None,
        missing: OptStr = None,
    ):
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
        self.monomers = monomers
        self.dtype = get_array_type(len(self))
        self._words = set(self)  # for quick lookup
        self.k = k

        self._gap_char = gap
        self._gap_index = self.index(gap) if gap else None
        self._missing_char = missing
        self._missing_index = self.index(missing) if missing else None
        if gap:
            assert self.monomers.gap_char is not None
            assert gap == self.monomers.gap_char * k
            assert gap in self

        if missing:
            assert self.monomers.missing_char is not None
            assert missing == self.monomers.missing_char * k
            assert missing in self

        self._coeffs = coord_conversion_coeffs(
            self.monomers.num_canonical, k, dtype=self.dtype
        )
        self._num_canonical = self.monomers.num_canonical**k

    @property
    def gap_char(self) -> OptStr:
        return self._gap_char

    @property
    def gap_index(self) -> OptInt:
        return self._gap_index

    @property
    def missing_char(self) -> OptStr:
        return self._missing_char

    @property
    def missing_index(self) -> OptInt:
        return self._missing_index

    @property
    def num_canonical(self) -> int:
        return self._num_canonical

    @functools.singledispatchmethod
    def to_indices(self, seq, independent_kmer: bool = True) -> numpy.ndarray:
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
        # todo: handle case of non-modulo sequences
        raise TypeError(f"{type(seq)} is invalid")

    @to_indices.register
    def _(self, seq: str, independent_kmer: bool = True) -> numpy.ndarray:
        seq = self.monomers.to_indices(seq)
        return self.to_indices(seq, independent_kmer=independent_kmer)

    @to_indices.register
    def _(self, seq: numpy.ndarray, independent_kmer: bool = True) -> numpy.ndarray:
        size = len(seq) - self.k + 1
        if independent_kmer:
            size = int(numpy.ceil(size / self.k))
        result = numpy.zeros(size, dtype=get_array_type(size))
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

    def from_indices(
        self, kmer_indices: numpy.ndarray, independent_kmer: bool = True
    ) -> numpy.ndarray:
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

    def with_gap_motif(self, include_missing: bool = False):
        """returns a new KmerAlphabet with the gap motif added

        Notes
        -----
        Adds gap state to monomers and recreates k-mer alphabet
        for self
        """
        if self.gap_char:
            return self
        missing = "?" if include_missing else None
        monomers = self.monomers.with_gap_motif(missing_char=missing)
        return monomers.get_kmer_alphabet(self.k, include_gap=True)

    @functools.singledispatchmethod
    def to_index(self, seq) -> numpy.ndarray:
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
        raise TypeError(f"{type(seq)} not supported")

    @to_index.register
    def _(self, seq: str) -> numpy.ndarray:
        seq = self.monomers.to_indices(seq)
        return self.to_index(seq)

    @to_index.register
    def _(self, seq: bytes) -> numpy.ndarray:
        seq = self.monomers.to_indices(seq)
        return self.to_index(seq)

    @to_index.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        assert len(seq) == self.k
        if (self.monomers.num_canonical > seq).all():
            return coord_to_index(seq, self._coeffs)

        gap_char_index = self.monomers.gap_index or self.monomers.num_canonical
        gap_index = self.gap_index or self.num_canonical
        missing_index = gap_index + 1 if self.gap_char else gap_index
        return gap_index if seq.max() == gap_char_index else missing_index

    def from_index(self, kmer_index: int) -> numpy.ndarray:
        """decodes an integer into a k-mer"""
        if self.gap_index is not None and kmer_index == self.gap_index:
            return numpy.array([self.monomers.gap_index] * self.k, dtype=self.dtype)
        elif self.missing_index is not None and kmer_index == self.missing_index:
            return numpy.array([self.monomers.missing_index] * self.k, dtype=self.dtype)
        elif kmer_index >= len(self):
            raise ValueError(f"{kmer_index} is out of range")

        return index_to_coord(kmer_index, self._coeffs)

    @functools.singledispatchmethod
    def is_valid(self, seq: numpy.ndarray) -> bool:
        """whether integers are within the valid range

        Parameters
        ----------
        seq
            a numpy array of integers

        Notes
        -------
        This will raise a TypeError for string or bytes. Using to_indices() to
        convert those ensures a valid result.
        """
        raise TypeError(f"{type(seq)} is invalid, must be numpy.ndarray")

    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        max_val = max(self.missing_index or 0, self.gap_index or 0, self.num_canonical)
        return seq.min() >= 0 and seq.max() <= max_val if len(seq) else True

    def to_rich_dict(self) -> dict:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        return {
            "words": list(self),
            "monomers": self.monomers.to_rich_dict(),
            "k": self.k,
            "gap": self.gap_char,
            "missing": self.missing_char,
            "version": __version__,
            "type": get_object_provenance(self),
        }

    def to_json(self) -> str:
        """returns a serialisable string"""
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, data: dict) -> "KmerAlphabet":
        """returns an instance from a serialised dictionary"""
        monomers = CharAlphabet.from_rich_dict(data["monomers"])
        return cls(
            words=data["words"],
            monomers=monomers,
            k=data["k"],
            gap=data["gap"],
            missing=data["missing"],
        )


@register_deserialiser(get_object_provenance(KmerAlphabet))
def deserialise_kmer_alphabet(data: dict) -> KmerAlphabet:
    return KmerAlphabet.from_rich_dict(data)


class CodonAlphabet(tuple):
    """represents the sense-codons of a GeneticCode"""

    def __new__(
        cls,
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        gap: OptStr = None,
    ):
        if not words:
            raise ValueError(f"cannot create empty {cls.__name__!r}")

        if gap is not None:
            assert _coerce_to_type(words[0], gap) in words

        consistent_words(words, length=3)
        return tuple.__new__(cls, words, monomers=monomers, gap=gap)

    def __init__(
        self,
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        gap: OptStr = None,
    ):
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
        self.motif_length = 3

    @property
    def gap_char(self):
        return self._gap_char

    @property
    def gap_index(self):
        return self._to_indices(self.gap_char) if self._gap_char else None

    @property
    def missing_char(self):
        """not supported on CodonAlphabet"""
        return None

    @property
    def missing_index(self):
        """not supported on CodonAlphabet"""
        return None

    @functools.singledispatchmethod
    def to_indices(self, seq) -> numpy.ndarray:
        """returns a sequence of codon indices"""
        raise TypeError(f"{type(seq)} is not supported")

    @to_indices.register
    def _(self, seq: str) -> numpy.ndarray:
        return self.to_indices([seq[i : i + 3] for i in range(0, len(seq), 3)])

    @to_indices.register
    def _(self, seq: list) -> numpy.ndarray:
        return [self.to_index(c) for c in seq]

    def to_index(self, codon: str) -> int:
        if len(codon) != 3:
            raise ValueError(f"{codon=!r} is not of length 3")
        try:
            return self._to_indices[codon]
        except KeyError as e:
            raise ValueError(f"{codon=!r} not in alphabet") from e

    def from_index(self, index: int) -> str:
        if index > len(self) or index < 0:
            raise ValueError(f"{index=!r} is not within range")

        try:
            return self._from_indices[index]
        except KeyError as e:
            raise ValueError(f"invalid {index=}") from e

    def from_indices(self, indices: numpy.ndarray) -> list[str]:
        return [self.from_index(index) for index in indices]

    @property
    def num_canonical(self):
        return len(self._words)

    @functools.singledispatchmethod
    def is_valid(self, seq) -> bool:
        """seq is valid for alphabet"""
        raise TypeError(f"{type(seq)} not supported")

    @is_valid.register
    def _(self, seq: str) -> bool:
        try:
            _ = self.to_indices(seq)
            return True
        except ValueError:
            return False

    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        return seq.min() >= 0 and seq.max() < len(self)

    def with_gap_motif(self):
        if self.gap_char:
            return self
        monomers = self.monomers.with_gap_motif()
        gap_char = monomers.gap_char * 3
        words = tuple(self) + (gap_char,)
        return self.__class__(words=words, monomers=monomers, gap=gap_char)

    def to_rich_dict(self):
        from cogent3._version import __version__

        return {
            "monomers": self.monomers.to_rich_dict(),
            "type": get_object_provenance(self),
            "version": __version__,
            "words": list(self),
        }

    def to_json(self):
        return json.dumps(self.to_rich_dict())

    @classmethod
    def from_rich_dict(cls, data):
        data["monomers"] = deserialise_char_alphabet(data["monomers"])
        data.pop("type", None)
        data.pop("version", None)
        return cls(**data)


@register_deserialiser(get_object_provenance(CodonAlphabet))
def deserialise_codon_alphabet(data: dict) -> CodonAlphabet:
    return CodonAlphabet.from_rich_dict(data)


_alphabet_moltype_map = {}


def make_alphabet(*, chars, gap, missing, moltype):
    """constructs a character alphabet and registers the associated moltype

    Notes
    -----
    The moltype is associated with the alphabet, available as an
    alphabet property.
    """
    alpha = CharAlphabet(chars, gap=gap, missing=missing)
    _alphabet_moltype_map[alpha] = moltype
    return alpha
