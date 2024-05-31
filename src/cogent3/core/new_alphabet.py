import functools
import typing

from abc import ABC, abstractmethod

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
        b"RYYR"
        >>> dna_to_ry = convert_alphabet(src=b"TCAG", dest=b"YYRR", delete=b"-")
        >>> dna_to_ry(b"GC-TA")
        b"RYYR"
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
    def to_indices(self, seq: str) -> numpy.ndarray: ...

    @abstractmethod
    def from_indices(self, seq: numpy.ndarray) -> str: ...

    @abstractmethod
    def is_valid(self, seq) -> bool: ...

    @abstractmethod
    def with_gap_motif(self): ...


class MonomerAlphabetABC(ABC):
    @abstractmethod
    def get_word_alphabet(self, size: int): ...


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
    """converts utf8 to numpy array"""

    def __init__(self, chars: bytes, dtype):
        # we want a bytes translation map
        self._table = b"".maketrans(
            chars,
            bytes(bytearray(range(len(chars)))),
        )
        self.dtype = dtype

    def __call__(self, seq: bytes) -> numpy.ndarray:
        b = seq.translate(self._table)
        return numpy.array(memoryview(b), dtype=self.dtype)


class array_to_bytes:
    """converts utf8 to numpy array"""

    def __init__(self, chars, dtype):
        # we want a bytes translation map
        self._table = b"".maketrans(
            bytes(bytearray(range(len(chars)))),
            chars,
        )
        self.dtype = dtype

    def __call__(self, seq: numpy.ndarray) -> bytes:
        b = seq.tobytes().translate(self._table)
        return bytearray(b)


class CharAlphabet(tuple, AlphabetABC, MonomerAlphabetABC):
    def __new__(cls, chars: typing.Sequence[StrORBytes], gap: OptStr = None):
        if not chars:
            raise ValueError(f"cannot create empty {cls.__name__!r}")

        if gap is not None:
            assert _coerce_to_type(chars[0], gap) in chars

        consistent_words(chars, length=1)
        return tuple.__new__(cls, chars, gap=gap)

    def __init__(self, chars: typing.Sequence[StrORBytes], gap: str = "-"):
        self._gap = gap
        self.dtype = get_array_type(len(self))
        self._chars = set(self)  # for quick lookup
        try:
            byte_chars = b"".join([c.encode("utf8") for c in chars])
        except AttributeError:
            byte_chars = bytes(bytearray(chars))
        self._bytes2arr = bytes_to_array(byte_chars, dtype=self.dtype)
        self._arr2bytes = array_to_bytes(byte_chars, dtype=self.dtype)

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

    def get_word_alphabet(self, size: int): ...

    @functools.singledispatchmethod
    def is_valid(self, seq) -> bool:
        raise TypeError(f"{type(seq)} is invalid")

    @is_valid.register
    def _(self, seq: str) -> bool:
        return set(seq) <= self._chars

    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        return seq.min() >= 0 and seq.max() < len(self)


class KmerAlphabet(AlphabetABC):
    def to_indices(self) -> numpy.ndarray: ...

    def from_indices(self) -> numpy.ndarray: ...

    def with_gap_motif(self): ...


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
