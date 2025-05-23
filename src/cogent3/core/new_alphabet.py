import functools
import itertools
import json
import typing
from abc import ABC, abstractmethod

import numba
import numpy
import typing_extensions

from cogent3.util import warning as c3warn
from cogent3.util.deserialise import register_deserialiser
from cogent3.util.misc import get_object_provenance

if typing.TYPE_CHECKING:
    from cogent3.core.new_moltype import MolType

NumpyIntType = numpy.dtype[numpy.integer]
NumpyIntArrayType = numpy.typing.NDArray[numpy.integer]
StrORBytes = str | bytes
StrORArray = str | NumpyIntArrayType
StrORBytesORArray = str | bytes | NumpyIntArrayType
OptInt = int | None
OptStr = str | None
OptBytes = bytes | None
PySeqStrOrBytes = typing.Sequence[str | bytes]


@functools.singledispatch
def _coerce_to_type(orig: StrORBytes, text: str) -> StrORBytes:
    msg = f"{type(orig)} is invalid"
    raise TypeError(msg)


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
    return _coerce_to_type(orig[0], text)


@_coerce_to_type.register
def _(orig: list, text: str) -> StrORBytes:
    return _coerce_to_type(orig[0], text)


class AlphabetError(TypeError): ...


class convert_alphabet:
    """convert one character set into another"""

    def __init__(
        self,
        src: bytes,
        dest: bytes,
        delete: OptBytes = None,
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
    def to_indices(self, seq: StrORBytesORArray) -> NumpyIntArrayType: ...

    @abstractmethod
    def from_indices(self, seq: StrORBytesORArray) -> str: ...

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
    def to_rich_dict(self, for_pickle: bool = False) -> dict[str, typing.Any]: ...

    @abstractmethod
    def to_json(self) -> str: ...

    @classmethod
    @abstractmethod
    def from_rich_dict(cls, data: dict) -> typing_extensions.Self: ...

    def __deepcopy__(self, memo) -> typing_extensions.Self:
        return type(self).from_rich_dict(self.to_rich_dict())

    @property
    def moltype(self) -> typing.Union["MolType", None]:
        return _alphabet_moltype_map.get(self)


class MonomerAlphabetABC(ABC):
    @abstractmethod
    def get_kmer_alphabet(
        self,
        k: int,
        include_gap: bool = True,
    ) -> "KmerAlphabet": ...

    @abstractmethod
    def as_bytes(self) -> bytes: ...

    @abstractmethod
    def convert_seq_array_to(
        self,
        *,
        alphabet: typing_extensions.Self,
        seq: NumpyIntArrayType,
        check_valid: bool = True,
    ) -> NumpyIntArrayType: ...

    @abstractmethod
    def with_gap_motif(
        self,
        gap_char: str = "-",
        missing_char: str = "?",
        include_missing: bool = False,
        gap_as_state: bool = False,
    ) -> typing_extensions.Self: ...


def get_array_type(num_elements: int) -> NumpyIntType:
    """Returns the smallest numpy integer dtype that can contain elements
    within num_elements.
    """
    if num_elements <= 2**8:
        dtype = numpy.uint8
    elif num_elements <= 2**16:
        dtype = numpy.uint16
    elif num_elements <= 2**32:
        dtype = numpy.uint32
    elif num_elements <= 2**64:
        dtype = numpy.uint64
    else:
        msg = f"{num_elements} is too big for 64-bit integer"
        raise NotImplementedError(msg)

    return dtype


def consistent_words(
    words: PySeqStrOrBytes,
    length: OptInt = None,
) -> None:
    """make sure all alphabet elements are unique and have the same length"""
    if not words:
        msg = "no words provided"
        raise ValueError(msg)

    if len(set(words)) != len(words):
        msg = f"duplicate elements in {words}"
        raise ValueError(msg)

    lengths = {1} if isinstance(words[0], int) else {len(w) for w in words}
    if len(lengths) != 1:
        msg = f"mixed lengths {lengths} in {words}"
        raise ValueError(msg)
    l = next(iter(lengths))
    if length and l != length:
        msg = f"word length {l} does not match expected {length}"
        raise ValueError(msg)


class bytes_to_array:
    """wrapper around convert_alphabet. It defines a linear mapping from provided
    characters to uint8. The resulting object is callable, taking a bytes object
    and returning a numpy array."""

    def __init__(
        self, chars: bytes, dtype: NumpyIntType, delete: OptBytes = None
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


class CharAlphabet(tuple, AlphabetABC, MonomerAlphabetABC):
    """representing fundamental monomer character sets.

    Notes
    -----
    Provides methods for efficient conversion between characters and integers
    from fundamental types of strings, bytes and numpy arrays.
    """

    __slots__ = ()

    def __new__(
        cls,
        chars: typing.Sequence[StrORBytes],
        gap: OptStr = None,
        missing: OptStr = None,
    ) -> "CharAlphabet":
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

        if gap is not None:
            assert _coerce_to_type(chars[0], gap) in chars

        if missing is not None and _coerce_to_type(chars[0], missing) not in chars:
            msg = (
                f"char missing={_coerce_to_type(chars[0], missing)!r} not in {chars!r}"
            )
            raise ValueError(msg)

        consistent_words(chars, length=1)
        if isinstance(chars, bytes):
            # we need to convert to tuple in a way that preserves elements
            # as bytes
            chars = tuple(bytes([c]) for c in chars)
        return tuple.__new__(cls, chars, gap=gap, missing=missing)

    def __init__(
        self,
        chars: typing.Sequence[StrORBytes],
        gap: OptStr = None,
        missing: OptStr = None,
    ) -> None:
        self.dtype = get_array_type(len(self))
        self._gap_char = gap
        self._gap_index = self.dtype(self.index(gap)) if gap else None
        self._missing_char = missing
        self._missing_index = self.dtype(self.index(missing)) if missing else None

        # the number of canonical states are non-gap, non-missing states
        adj = (1 if gap else 0) + (1 if missing else 0)

        self._num_canonical = len(self) - adj

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

    @c3warn.deprecated_callable(
        version="2025.6",
        reason="Replaced by .motif_len",
        is_discontinued=True,
    )
    def get_motif_len(self) -> int:
        """the size of each member of the alphabet"""
        # added to maintain compatibility with the old API
        return self.motif_len

    @functools.singledispatchmethod
    def to_indices(self, seq: StrORBytesORArray | tuple) -> NumpyIntArrayType:
        """returns a sequence of indices for the characters in seq"""
        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    @to_indices.register
    def _(self, seq: tuple) -> NumpyIntArrayType:
        indices = []
        for c in seq:
            indices.extend(self.to_indices(c).tolist())
        return numpy.array(indices, dtype=self.dtype)

    @to_indices.register
    def _(self, seq: bytes) -> NumpyIntArrayType:
        # any non-canonical characters should lie outside the range
        # we replace these with a single value
        return self._bytes2arr(seq)

    @to_indices.register
    def _(self, seq: str) -> NumpyIntArrayType:
        return self.to_indices(seq.encode("utf8"))

    @to_indices.register
    def _(self, seq: numpy.ndarray) -> NumpyIntArrayType:
        return seq.astype(self.dtype)

    @functools.singledispatchmethod
    def from_indices(self, seq: StrORBytesORArray) -> str:
        """returns a string from a sequence of indices"""
        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    @from_indices.register
    def _(self, seq: str) -> str:
        return seq

    @from_indices.register
    def _(self, seq: bytes) -> str:
        return seq.decode("utf8")

    @from_indices.register
    def _(self, seq: numpy.ndarray) -> str:
        return self._arr2bytes(seq).decode("utf8")

    def with_gap_motif(
        self,
        gap_char: str = "-",
        missing_char: str = "?",
        include_missing: bool = False,
        gap_as_state: bool = False,
    ) -> typing_extensions.Self:
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
        chars = (
            (*tuple(self[: self._num_canonical]), gap_char, missing_char)
            if missing_char
            else (*tuple(self[: self._num_canonical]), gap_char)
        )
        if gap_as_state:
            gap_char = None
        return self.__class__(chars, gap=gap_char, missing=missing_char)

    def get_kmer_alphabet(self, k: int, include_gap: bool = True) -> "KmerAlphabet":
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

    @functools.singledispatchmethod
    def is_valid(self, seq: StrORBytesORArray) -> bool:
        """seq is valid for alphabet"""
        # refactor: design
        if hasattr(seq, "alphabet"):
            # assume a SeqView instance
            return seq.alphabet == self
        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    @is_valid.register
    def _(self, seq: str) -> bool:
        return self.is_valid(self.to_indices(seq))

    @is_valid.register
    def _(self, seq: bytes) -> bool:
        return self.is_valid(self.to_indices(seq))

    # cannot use NumpyIntArrayType as a type hint here
    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        return bool(seq.min() >= 0 and seq.max() < len(self) if len(seq) else True)

    def as_bytes(self) -> bytes:
        """returns self as a byte string"""
        if isinstance(self[0], bytes):
            return b"".join(self)

        return "".join(self).encode("utf8")

    def array_to_bytes(self, seq: NumpyIntArrayType) -> bytes:
        """returns seq as a byte string"""
        return self._arr2bytes(seq)

    def to_rich_dict(self, for_pickle: bool = False) -> dict[str, typing.Any]:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        data = {
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
    def from_rich_dict(cls, data: dict) -> typing_extensions.Self:
        """returns an instance from a serialised dictionary"""
        return cls(data["chars"], gap=data["gap"], missing=data["missing"])

    @c3warn.deprecated_callable(
        version="2025.6",
        reason="Replaced by get_kmer_alphabet",
        is_discontinued=True,
    )
    def get_word_alphabet(self, k: int, include_gap: bool = True) -> "KmerAlphabet":
        return self.get_kmer_alphabet(k, include_gap=include_gap)

    def get_subset(
        self,
        motif_subset: PySeqStrOrBytes,
        excluded: bool = False,
    ) -> typing_extensions.Self:
        """Returns a new Alphabet object containing a subset of motifs in self.

        Raises an exception if any of the items in the subset are not already
        in self.
        """
        self_set = set(self)
        self_set |= {self.gap_char, self.missing_char}
        self_set.discard(None)

        if diff := set(motif_subset) - self_set:
            msg = f"{diff!r} not members in self"
            raise AlphabetError(msg)

        if excluded:
            motif_subset = [m for m in self if m not in motif_subset]
        gap = self.gap_char if self.gap_char in motif_subset else None
        missing = self.missing_char if self.missing_char in motif_subset else None
        return make_alphabet(
            chars=tuple(motif_subset),
            gap=gap,
            missing=missing,
            moltype=self.moltype,
        )

    def convert_seq_array_to(
        self,
        *,
        alphabet: typing_extensions.Self,
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
    monomers: typing.Sequence[str],
    num_canonical: int,
    gap_char: str | None,
    missing_char: str | None,
    include_gap: bool,
    k: int,
    moltype: "MolType",
) -> "KmerAlphabet":
    chars = tuple(monomers[:num_canonical])

    gap = gap_char * k if include_gap and gap_char is not None else None
    missing = missing_char * k if missing_char is not None and include_gap else None
    words = tuple("".join(e) for e in itertools.product(chars, repeat=k))
    words += (gap,) if include_gap and gap else ()
    words += (missing,) if include_gap and missing else ()

    kalpha = KmerAlphabet(words=words, monomers=monomers, gap=gap, k=k, missing=missing)
    if moltype:
        _alphabet_moltype_map[kalpha] = moltype
    return kalpha


@functools.cache
def make_converter(
    src_alpha: CharAlphabet,
    dest_alpha: CharAlphabet,
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
    src_indices = [src_alpha.to_indices(common), src_alpha.to_indices(diff)]
    dest_indices = [
        dest_alpha.to_indices(common),
        numpy.array([len(dest_alpha)] * len(diff), dtype=dest_alpha.dtype),
    ]

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

    src_indices = numpy.concatenate(src_indices)
    dest_indices = numpy.concatenate(dest_indices)

    return convert_alphabet(
        src_indices.tobytes(),
        dest_indices.tobytes(),
        allow_duplicates=True,
    )


@register_deserialiser(get_object_provenance(CharAlphabet))
def deserialise_char_alphabet(data: dict) -> CharAlphabet:
    return CharAlphabet.from_rich_dict(data)


@numba.jit(nopython=True)
def coord_conversion_coeffs(
    num_states: int,
    k: int,
    dtype: numpy.dtype | None = None,
) -> numpy.ndarray:  # pragma: no cover
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


@numba.jit(nopython=True)
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
        msg = f"size of result {len(result)} <= {size}"
        raise ValueError(msg)

    missing_index: int = gap_index + 1 if gap_char_index > 0 else num_states**k
    if gap_char_index > 0:
        assert gap_index > 0, (
            f"gap_index={gap_index} but gap_char_index={gap_char_index}"
        )

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
    kmer_indices: numpy.ndarray,
    result: numpy.ndarray,
    coeffs: numpy.ndarray,
    k: int,
    gap_index: int = -1,
    gap_char_index: int = -1,
    independent_kmer: bool = True,
) -> numpy.ndarray:  # pragma: no cover
    if gap_index > 0:
        assert gap_char_index > 0, (
            f"gap_index={gap_index} but gap_char_index={gap_char_index}"
        )
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

    @abstractmethod
    def with_gap_motif(
        self,
        include_missing: bool = False,
        **kwargs,
    ) -> typing_extensions.Self: ...


class KmerAlphabet(tuple, AlphabetABC, KmerAlphabetABC):
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
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        k: int,
        gap: OptStr = None,
        missing: OptStr = None,
    ) -> "KmerAlphabet":
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
            self.monomers.num_canonical,
            k,
            dtype=self.dtype,
        )
        self._num_canonical = self.monomers.num_canonical**k

    def __reduce__(
        self,
    ) -> tuple[type, tuple[tuple[str, ...], CharAlphabet, int, OptStr, OptStr]]:
        """support for pickling"""
        return (
            self.__class__,
            (tuple(self), self.monomers, self.k, self.gap_char, self.missing_char),
            None,
        )

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
        # TODO: handle case of non-modulo sequences
        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    @to_indices.register
    def _(self, seq: tuple, **kwargs) -> numpy.ndarray[int]:
        return numpy.array([self.to_index(c) for c in seq], dtype=self.dtype)

    @to_indices.register
    def _(self, seq: list, **kwargs) -> numpy.ndarray:
        return numpy.array([self.to_index(c) for c in seq], dtype=self.dtype)

    @to_indices.register
    def _(self, seq: str, independent_kmer: bool = True) -> numpy.ndarray:
        seq = self.monomers.to_indices(seq)
        return self.to_indices(seq, independent_kmer=independent_kmer)

    @to_indices.register
    def _(self, seq: numpy.ndarray, independent_kmer: bool = True) -> numpy.ndarray:
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

    def from_indices(
        self,
        kmer_indices: numpy.ndarray,
        independent_kmer: bool = True,
    ) -> numpy.ndarray:
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

    def with_gap_motif(
        self,
        include_missing: bool = False,
        **kwargs,
    ) -> typing_extensions.Self:
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

    @functools.singledispatchmethod
    def to_index(self, seq) -> int:
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
        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    @to_index.register
    def _(self, seq: str) -> int:
        seq = self.monomers.to_indices(seq)
        return self.to_index(seq)

    @to_index.register
    def _(self, seq: bytes) -> int:
        seq = self.monomers.to_indices(seq)
        return self.to_index(seq)

    @to_index.register
    def _(self, seq: numpy.ndarray) -> int:
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
        if self.missing_index is not None and kmer_index == self.missing_index:
            return numpy.array([self.monomers.missing_index] * self.k, dtype=self.dtype)
        if kmer_index >= len(self):
            msg = f"{kmer_index} is out of range"
            raise ValueError(msg)

        return index_to_coord(kmer_index, self._coeffs)

    @functools.singledispatchmethod
    def is_valid(self, seq) -> bool:
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
        msg = f"{type(seq)} is invalid, must be numpy.ndarray"
        raise TypeError(msg)

    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        max_val = max(self.missing_index or 0, self.gap_index or 0, self.num_canonical)
        return bool(seq.min() >= 0 and seq.max() <= max_val if len(seq) else True)

    def to_rich_dict(self, for_pickle: bool = False) -> dict:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        data = {
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

    @property
    def motif_len(self) -> int:
        return self.k


@register_deserialiser(get_object_provenance(KmerAlphabet))
def deserialise_kmer_alphabet(data: dict) -> KmerAlphabet:
    return KmerAlphabet.from_rich_dict(data)


class SenseCodonAlphabet(tuple, AlphabetABC, KmerAlphabetABC):
    """represents the sense-codons of a GeneticCode"""

    __slots__ = ()

    def __new__(
        cls,
        words: tuple[StrORBytes, ...],
        monomers: CharAlphabet,
        gap: OptStr = None,
    ) -> "SenseCodonAlphabet":
        if not words:
            msg = f"cannot create empty {cls.__name__!r}"
            raise ValueError(msg)

        if gap is not None:
            assert _coerce_to_type(words[0], gap) in words

        consistent_words(words, length=3)
        return tuple.__new__(cls, words, monomers=monomers, gap=gap)

    def __init__(
        self,
        words: tuple[StrORBytes, ...],  # noqa: ARG002
        monomers: CharAlphabet,
        gap: OptStr = None,
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
    ) -> tuple[type, tuple[tuple[str, ...], CharAlphabet, OptStr], None]:
        """support for pickling"""
        return (
            self.__class__,
            (tuple(self), self.monomers, self.gap_char),
            None,
        )

    @property
    def gap_char(self) -> str:
        return self._gap_char

    @property
    def gap_index(self) -> OptInt:
        return self._to_indices[self.gap_char] if self._gap_char else None

    @property
    def missing_char(self) -> None:
        """not supported on CodonAlphabet"""
        return None

    @property
    def missing_index(self) -> None:
        """not supported on CodonAlphabet"""
        return None

    @functools.singledispatchmethod
    def to_indices(self, seq) -> numpy.ndarray:
        """returns a sequence of codon indices"""
        msg = f"{type(seq)} is invalid"
        raise TypeError(msg)

    @to_indices.register
    def _(self, seq: str) -> numpy.ndarray:
        return self.to_indices([seq[i : i + 3] for i in range(0, len(seq), 3)])

    @to_indices.register
    def _(self, seq: list) -> numpy.ndarray:
        return numpy.array([self.to_index(c) for c in seq], dtype=self.dtype)

    @to_indices.register
    def _(self, seq: tuple) -> numpy.ndarray[int]:
        return self.to_indices(list(seq))

    @to_indices.register
    def _(self, seq: numpy.ndarray) -> numpy.ndarray:
        # we assume that this is a dna sequence encoded as a numpy array
        size = len(seq) // 3
        return self.to_indices(
            [self.monomers.from_indices(c) for c in seq.reshape(size, 3)],
        )

    def to_index(self, codon: str) -> int:
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

    def from_index(self, index: int) -> str:
        """returns a single codon from an index"""
        if index > len(self) or index < 0:
            msg = f"{index=!r} is not within range"
            raise ValueError(msg)

        try:
            return self._from_indices[index]
        except KeyError as e:
            msg = f"invalid {index=}"
            raise ValueError(msg) from e

    def from_indices(self, indices: NumpyIntArrayType) -> list[str]:
        """returns a list of codons from a numpy array of indices"""
        return [self.from_index(index) for index in indices]

    @property
    def num_canonical(self) -> int:
        """returns the number of canonical states"""
        return len(self._words)

    @functools.singledispatchmethod
    def is_valid(self, seq) -> bool:
        """seq is valid for alphabet"""
        msg = f"{type(seq)} not supported"
        raise TypeError(msg)

    @is_valid.register
    def _(self, seq: str) -> bool:
        try:
            _ = self.to_indices(seq)
        except (ValueError, AlphabetError):
            return False
        else:
            return True

    @is_valid.register
    def _(self, seq: numpy.ndarray) -> bool:
        return bool(seq.min() >= 0 and seq.max() < len(self) if len(seq) else True)

    def with_gap_motif(self, **kwargs) -> typing_extensions.Self:
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

    def to_rich_dict(self, for_pickle: bool = False) -> dict:
        """returns a serialisable dictionary"""
        from cogent3._version import __version__

        data = {
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
    def from_rich_dict(cls, data: dict[str, typing.Any]) -> "SenseCodonAlphabet":
        """returns an instance from a serialised dictionary"""
        data["monomers"] = deserialise_char_alphabet(data["monomers"])
        data.pop("type", None)
        data.pop("version", None)
        return cls(**data)

    @property
    def motif_len(self) -> int:
        """always 3 for codon alphabets"""
        return self._motif_len


@register_deserialiser(get_object_provenance(SenseCodonAlphabet))
def deserialise_codon_alphabet(data: dict) -> SenseCodonAlphabet:
    return SenseCodonAlphabet.from_rich_dict(data)


_alphabet_moltype_map = {}


def make_alphabet(
    *,
    chars: typing.Sequence[str],
    gap: str | None,
    missing: str | None,
    moltype: "MolType",
) -> CharAlphabet:
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
