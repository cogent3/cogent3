import functools
import itertools
import typing

from abc import ABC, abstractmethod

import numba
import numpy


StrORBytes = typing.Union[str, bytes]
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
    def motif_len(self) -> int: ...

    @abstractmethod
    def is_valid(self, seq) -> bool: ...

    @abstractmethod
    def with_gap_motif(self): ...


class MonomerAlphabetABC(ABC):
    @abstractmethod
    def get_kmer_alphabet(self, size: int): ...

    @abstractmethod
    def to_indices(self, seq: str) -> numpy.ndarray: ...

    @abstractmethod
    def from_indices(self, seq: numpy.ndarray) -> str: ...

    @abstractmethod
    def to_bytes(self) -> bytes: ...


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
        return bytearray(b)


class CharAlphabet(tuple, AlphabetABC, MonomerAlphabetABC):
    def __new__(cls, chars: typing.Sequence[StrORBytes], gap: OptStr = None):
        if not chars:
            raise ValueError(f"cannot create empty {cls.__name__!r}")

        if gap is not None:
            assert _coerce_to_type(chars[0], gap) in chars

        consistent_words(chars, length=1)
        return tuple.__new__(cls, chars, gap=gap)

    def __init__(self, chars: typing.Sequence[StrORBytes], gap: str = None):
        self._gap_char = gap
        self._gap_index = self.index(gap) if gap else None
        self.dtype = get_array_type(len(self))
        self._chars = set(self)  # for quick lookup
        try:
            byte_chars = b"".join([c.encode("utf8") for c in chars])
        except AttributeError:
            byte_chars = bytes(bytearray(chars))
        self._bytes2arr = bytes_to_array(byte_chars, dtype=self.dtype)
        self._arr2bytes = array_to_bytes(byte_chars)

    @property
    def gap_char(self) -> OptStr:
        return self._gap_char

    @property
    def gap_index(self) -> OptInt:
        return self._gap_index

    @property
    def motif_len(self) -> int:
        # this is enforced by consistent_words(length=1) in the constructor
        return 1

    @functools.singledispatchmethod
    def to_indices(self, seq: str) -> numpy.ndarray:
        raise TypeError(f"{type(seq)} is invalid")

    @to_indices.register
    def _(self, seq: bytes) -> numpy.ndarray:
        return self._bytes2arr(seq)

    @to_indices.register
    def _(self, seq: str) -> numpy.ndarray:
        return self._bytes2arr(seq.encode("utf8"))

    @to_indices.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        return seq

    @functools.singledispatchmethod
    def from_indices(self, seq: numpy.ndarray) -> str:
        raise TypeError(f"{type(seq)} is invalid")

    @from_indices.register
    def _(self, seq: numpy.ndarray) -> str:
        return self._arr2bytes(seq).decode("utf8")

    def with_gap_motif(self): ...

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

        if self.gap_char:
            chars = tuple(c for c in self if c != self.gap_char)
        else:
            chars = tuple(self)

        gap = self._gap_char * k if include_gap and self._gap_char is not None else None
        words = tuple("".join(e) for e in itertools.product(chars, repeat=k))
        if include_gap and gap:
            words += (gap,)
        return KmerAlphabet(words=words, monomers=self, gap=gap, k=k)

    @functools.singledispatchmethod
    def is_valid(self, seq) -> bool:
        raise TypeError(f"{type(seq)} is invalid")

    @is_valid.register
    def _(self, seq: str) -> bool:
        return set(seq) <= self._chars

    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        return seq.min() >= 0 and seq.max() < len(self)

    def to_bytes(self) -> bytes:
        if not isinstance(self[0], str):
            return bytes(bytearray(self))

        return "".join(self).encode("utf8")


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


# todo: profile this against pure python,
#  appears slower as a numba function!
@numba.jit()
def seq_to_kmer_indices(
    seq: numpy.ndarray,
    result: numpy.ndarray,
    coeffs: numpy.ndarray,
    num_states: int,
    k: int,
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
        defines range of possible ints at a position
    k
        k-mer size
    independent_kmer
        if True, sets returns indices for non-overlapping k-mers, otherwise
        returns indices for all k-mers

    Notes
    -----
    If a k-mer has an element > num_states, that k-mer index is assigned
    to 1+num_states**k
    """
    # check the result length is consistent with the settings
    step: int = k if independent_kmer else 1
    size: int = int(numpy.ceil((len(seq) - k + 1) / step))
    if len(result) < size:
        raise ValueError(f"size of result {len(result)} <= {size}")

    for result_idx, i in enumerate(range(0, len(seq) - k + 1, step)):
        for j in range(i, i + k):
            if seq[j] < num_states:
                result[result_idx] = coord_to_index(seq[i : i + k], coeffs)
            else:
                continue
    return result


# todo: profile this against pure python
@numba.jit()
def kmer_indices_to_seq(
    kmer_indices: numpy.ndarray,
    result: numpy.ndarray,
    coeffs: numpy.ndarray,
    k: int,
    independent_kmer: bool = True,
) -> numpy.ndarray:  # pragma: no cover
    for index, kmer_index in enumerate(kmer_indices):
        coord = index_to_coord(kmer_index, coeffs)
        if index == 0:
            result[0:3] = coord
        elif independent_kmer:
            seq_index = index * k
            result[seq_index : seq_index + k] = coord
        else:
            # we are just adding the last monomer
            result[index + k - 1] = coord[-1]

    return result


class KmerAlphabet(tuple, AlphabetABC):
    def __new__(
        cls,
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        k: int,
        gap: OptStr = None,
    ):
        if not words:
            raise ValueError(f"cannot create empty {cls.__name__!r}")

        if gap is not None:
            assert _coerce_to_type(words[0], gap) in words

        consistent_words(words, length=k)
        return tuple.__new__(cls, words, monomers=monomers, gap=gap)

    def __init__(
        self,
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        k: int,
        gap: OptStr = None,
    ):
        self._gap_char = gap
        self.monomers = monomers
        self.dtype = get_array_type(len(self))
        self._words = set(self)  # for quick lookup
        self.k = k

        size = len(self.monomers) - 1 if self.monomers.gap_char else len(self.monomers)
        self._coeffs = coord_conversion_coeffs(size, k, dtype=self.dtype)

    @functools.singledispatchmethod
    def to_indices(self, seq, independent_kmer: bool = True) -> numpy.ndarray:
        """returns a sequence of k-mer indices

        Parameters
        ----------
        seq
            a sequence of monomers
        independent_kmer
            if True, returns non-overlapping k-mers
        """
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
        return seq_to_kmer_indices(
            seq, result, self._coeffs, len(self.monomers), self.k, independent_kmer
        )

    def from_indices(
        self, kmer_indices: numpy.ndarray, independent_kmer: bool = True
    ) -> numpy.ndarray:

        if independent_kmer:
            size = len(kmer_indices) * self.k
        else:
            size = len(kmer_indices) + self.k - 1

        result = numpy.zeros(size, dtype=self.dtype)
        return kmer_indices_to_seq(
            kmer_indices,
            result,
            self._coeffs,
            self.k,
            independent_kmer=independent_kmer,
        )

    def with_gap_motif(self): ...

    @functools.singledispatchmethod
    def to_index(self, seq) -> numpy.ndarray:
        """encodes a k-mer as a single integer

        Parameters
        ----------
        seq
            sequence to be encoded, can be either a string or numpy array
        overlapping
            if False, performs operation on sequential k-mers, e.g. codons
        """
        raise TypeError(f"{type(seq)} not supported")

    @to_index.register
    def _(self, seq: str) -> numpy.ndarray:
        """encodes a k-mer as a single integer

        Parameters
        ----------
        seq
            sequence to be encoded
        overlapping
            if False, performs operation on sequential k-mers, e.g. codons
        """
        seq = self.monomers.to_indices(seq)
        return self.to_index(seq)

    @to_index.register
    def _(self, seq: bytes) -> numpy.ndarray:
        """encodes a k-mer as a single integer

        Parameters
        ----------
        seq
            sequence to be encoded
        overlapping
            if False, performs operation on sequential k-mers, e.g. codons
        """
        seq = self.monomers.to_indices(seq)
        return self.to_index(seq)

    @to_index.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        """encodes a k-mer as a single integer

        Parameters
        ----------
        seq
            sequence to be encoded
        """
        return coord_to_index(seq, self._coeffs)

    def from_index(self, kmer_index: int) -> numpy.ndarray:
        """decodes an integer into a k-mer"""
        return index_to_coord(kmer_index, self._coeffs)


_alphabet_moltype_map = {}


def make_alphabet(*, chars, gap, moltype):
    """constructs an alphabet and registers the associated moltype

    Notes
    -----
    The moltype is associated with the alphabet, available as an
    alphabet property.
    """
    alpha = CharAlphabet(chars, gap=gap)
    _alphabet_moltype_map[alpha] = moltype
    return alpha
