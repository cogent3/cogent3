import numpy
import pytest

from cogent3.core.new_alphabet import (
    CharAlphabet,
    array_to_bytes,
    bytes_to_array,
    consistent_words,
    convert_alphabet,
    get_array_type,
)


@pytest.mark.parametrize(
    "num,expect",
    (
        (2**8 - 1, numpy.uint8),
        (2**16 - 1, numpy.uint16),
        (2**32 - 1, numpy.uint32),
        (2**64 - 1, numpy.uint64),
    ),
)
def test_get_array_type(num, expect):
    assert get_array_type(num) is expect


def test_get_array_type_invalid():
    with pytest.raises(NotImplementedError):
        get_array_type(2**64)


def test_make_char_alphabet():
    alpha = CharAlphabet(list("TCAG"))
    assert len(alpha) == 4


def test_convert_alphabet_mixed_length():
    with pytest.raises(ValueError):
        convert_alphabet(orig="ACG", new="TGCA")


@pytest.mark.parametrize("words", ("AACG", ["A", "CG"], ["CC", "AA"], "", (), []))
def test_invalid_words(words):
    with pytest.raises(ValueError):
        consistent_words(words, length=1)


def test_bytes2arr():
    b2a = bytes_to_array(b"TCAG", dtype=numpy.uint8)
    got = b2a(b"AAGTA")
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


def test_arr2bytes():
    a2b = array_to_bytes(b"TCAG", dtype=numpy.uint8)
    got = a2b(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == b"AAGTA"


@pytest.mark.parametrize("seq", (b"AAGTA", "AAGTA"))
def test_charalphabet_to_indices(seq):
    alpha = CharAlphabet(list("TCAG"))
    got = alpha.to_indices(seq)
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


def test_charalphabet_to_indices_invalid():
    alpha = CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.to_indices(list("TCAG"))


def test_charalphabet_from_indices():
    alpha = CharAlphabet(list("TCAG"))
    got = alpha.from_indices(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == "AAGTA"


def test_charalphabet_from_indices_invalid():
    alpha = CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.from_indices([0, 1])


def test_convert_alphabet():
    complement = convert_alphabet(b"ACGT", b"TGCA")
    seq = b"AAGG"
    assert complement(seq) == b"TTCC"


def test_charalphabet_with_gap():
    alpha = CharAlphabet(list("TCAG-"), gap="-")
    got = alpha.from_indices(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == "AAGTA"


def test_charalphabet_with_gap_missing():
    with pytest.raises(AssertionError):
        _ = CharAlphabet(list("TCAG"), gap="-")


def test_charalphabet_motif_len():
    alpha = CharAlphabet(list("TCAG"))
    assert alpha.motif_len == 1


@pytest.mark.parametrize("seq", ("ACGT", numpy.array([2, 1, 3, 0], dtype=numpy.uint8)))
def test_is_valid(seq):
    alpha = CharAlphabet(list("TCAG"))
    assert alpha.is_valid(seq)


@pytest.mark.parametrize(
    "seq", ("-ACGT", numpy.array([4, 2, 1, 3, 0], dtype=numpy.uint8))
)
def test_not_is_valid(seq):
    alpha = CharAlphabet(list("TCAG"))
    assert not alpha.is_valid(seq)


def test_is_valid_typerror():
    alpha = CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.is_valid(list("ACGT"))
