import itertools

import numpy
import pytest

from numpy.testing import assert_allclose

from cogent3.core import new_alphabet, new_moltype


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
    assert new_alphabet.get_array_type(num) is expect


def test_get_array_type_invalid():
    with pytest.raises(NotImplementedError):
        new_alphabet.get_array_type(2**64)


def test_make_char_alphabet():
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    assert len(alpha) == 4


def test_convert_alphabet_mixed_length():
    with pytest.raises(ValueError):
        new_alphabet.convert_alphabet(src="ACG", dest="TGCA")


@pytest.mark.parametrize("words", ("AACG", ["A", "CG"], ["CC", "AA"], "", (), []))
def test_invalid_words(words):
    with pytest.raises(ValueError):
        new_alphabet.consistent_words(words, length=1)


def test_bytes2arr():
    b2a = new_alphabet.bytes_to_array(b"TCAG", dtype=numpy.uint8)
    got = b2a(b"AAGTA")
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


def test_arr2bytes():
    a2b = new_alphabet.array_to_bytes(b"TCAG", dtype=numpy.uint8)
    got = a2b(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == b"AAGTA"


@pytest.mark.parametrize(
    "seq", (b"AAGTA", "AAGTA", numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
)
def test_charalphabet_to_indices(seq):
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    got = alpha.to_indices(seq)
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


def test_charalphabet_to_indices_invalid():
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.to_indices(list("TCAG"))


def test_charalphabet_from_indices():
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    got = alpha.from_indices(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == "AAGTA"


def test_charalphabet_from_indices_invalid():
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.from_indices([0, 1])


def test_convert_alphabet():
    complement = new_alphabet.convert_alphabet(b"ACGT", b"TGCA")
    seq = b"AAGG"
    assert complement(seq) == b"TTCC"


def test_convert_alphabet_with_delete():
    complement = new_alphabet.convert_alphabet(b"ACGT", b"TGCA", delete=b"-")
    seq = b"AAG--G"
    assert complement(seq) == b"TTCC"


def test_charalphabet_with_gap():
    alpha = new_alphabet.CharAlphabet(list("TCAG-"), gap="-")
    got = alpha.from_indices(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == "AAGTA"


def test_charalphabet_with_gap_missing():
    with pytest.raises(AssertionError):
        _ = new_alphabet.CharAlphabet(list("TCAG"), gap="-")


def test_charalphabet_motif_len():
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    assert alpha.motif_len == 1


@pytest.mark.parametrize("seq", ("ACGT", numpy.array([2, 1, 3, 0], dtype=numpy.uint8)))
def test_is_valid(seq):
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    assert alpha.is_valid(seq)


@pytest.mark.parametrize(
    "seq", ("-ACGT", numpy.array([4, 2, 1, 3, 0], dtype=numpy.uint8))
)
def test_not_is_valid(seq):
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    assert not alpha.is_valid(seq)


def test_is_valid_typerror():
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.is_valid(list("ACGT"))


@pytest.mark.parametrize("num_states,ndim", ((2, 1), (4, 2), (4, 4)))
def test_interconversion_of_coords_indices(num_states, ndim):
    # make sure match numpy functions
    coeffs = numpy.array(new_alphabet.coord_conversion_coeffs(num_states, ndim))
    for coord in itertools.product(*(list(range(num_states)),) * ndim):
        coord = numpy.array(coord, dtype=numpy.uint64)
        nidx = numpy.ravel_multi_index(coord, dims=(num_states,) * ndim)
        idx = new_alphabet.coord_to_index(coord, coeffs)
        assert idx == nidx
        ncoord = numpy.unravel_index(nidx, shape=(num_states,) * ndim)
        got = new_alphabet.index_to_coord(idx, coeffs)
        assert (got == ncoord).all()
        assert (got == coord).all()


def test_coord2index_fail():
    coord = numpy.array((0, 1), dtype=numpy.uint64)
    coeffs = numpy.array(new_alphabet.coord_conversion_coeffs(2, 4))
    # dimension of coords inconsistent with coeffs
    with pytest.raises(ValueError):
        new_alphabet.coord_to_index(coord, coeffs)


@pytest.mark.parametrize("k", (2, 3))
def test_kmer_alphabet_construction(k):
    dna = new_moltype.get_moltype("dna")
    monomers = dna.alphabet
    kmers = monomers.get_word_alphabet(k)
    assert len(kmers) == 4**k


@pytest.mark.parametrize("k", (2, 3))
def test_kmer_gapped_alphabet_construction(k):
    dna = new_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_word_alphabet(k, include_gap=True)
    assert len(kmers) == 1 + 4**k


def get_rand_coord_seq(k, cast):
    dna = new_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    coord = numpy.random.randint(0, 4, size=k, dtype=numpy.uint8)
    seq = monomers.from_indices(coord)
    if cast == bytes:
        seq = seq.encode("utf8")
    elif cast is None:
        seq = coord
    return coord, seq


@pytest.mark.parametrize("cast", (None, bytes, str))
@pytest.mark.parametrize("k", (2, 3))
def test_kmer_alphabet_pack(k, cast):
    # should match numpy ravel
    dna = new_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_word_alphabet(k, include_gap=True)
    coord, seq = get_rand_coord_seq(k, cast)
    expect = numpy.ravel_multi_index(coord, dims=(4,) * k)
    assert kmers.pack(seq) == expect


def test_kmer_alphabet_pack_invalid():
    # should match numpy ravel
    dna = new_moltype.get_moltype("dna")
    kmers = dna.gapped_alphabet.get_word_alphabet(3, include_gap=True)
    with pytest.raises(TypeError):
        kmers.pack([0, 1, 2])


@pytest.mark.parametrize("k", (2, 3))
def test_kmer_alphabet_unpack(k):
    # should match numpy ravel
    dna = new_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_word_alphabet(k, include_gap=True)
    coord = numpy.random.randint(0, 4, size=k, dtype=numpy.uint8)
    index = numpy.ravel_multi_index(coord, dims=(4,) * k)
    assert (kmers.unpack(index) == coord).all()


def test_seq_to_kmer_indices():
    numpy_func = numpy.ravel_multi_index
    dims = 4, 4, 4
    dna = new_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    coeffs = new_alphabet.coord_conversion_coeffs(4, 3)
    num = len(seq) // 3
    result = numpy.zeros(num, dtype=numpy.uint8)
    got = new_alphabet.seq_to_kmer_indices(arr, result, coeffs, 4, 3)
    expect = numpy.array(
        [numpy_func(arr[i : i + 3], dims=dims) for i in range(0, 9, 3)]
    )
    assert_allclose(got, expect)


def test_kmer_alphabet_to_indices_independent():
    numpy_func = numpy.ravel_multi_index
    dims = 4, 4, 4
    dna = new_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_word_alphabet(k=3, include_gap=False)
    got = trinuc_alpha.to_indices(arr, independent_kmer=True)
    expect = numpy.array(
        [numpy_func(arr[i : i + 3], dims=dims) for i in range(0, 9, 3)]
    )
    assert_allclose(got, expect)


def test_kmer_alphabet_to_indices_non_independent():
    numpy_func = numpy.ravel_multi_index
    dims = 4, 4, 4
    dna = new_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_word_alphabet(k=3, include_gap=False)
    got = trinuc_alpha.to_indices(arr, independent_kmer=False)
    assert len(got) == len(seq) - 3 + 1
    expect = numpy.array(
        [numpy_func(arr[i : i + 3], dims=dims) for i in range(0, 9 - 3 + 1)]
    )
    assert_allclose(got, expect)


def test_kmer_alphabet_to_indices_invalid():
    dna = new_moltype.get_moltype("dna")
    trinuc_alpha = dna.alphabet.get_word_alphabet(k=3, include_gap=False)
    with pytest.raises(TypeError):
        trinuc_alpha.to_indices(["ACG", "TGG"])


def test_kmer_alphabet_from_indices_independent():
    dna = new_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_word_alphabet(k=3, include_gap=False)
    kmer_seq = trinuc_alpha.to_indices(arr, independent_kmer=True)
    got = trinuc_alpha.from_indices(kmer_seq)
    assert_allclose(got, arr)


def test_kmer_alphabet_from_indices_not_independent():
    dna = new_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_word_alphabet(k=3, include_gap=False)
    kmer_seq = trinuc_alpha.to_indices(arr, independent_kmer=False)
    got = trinuc_alpha.from_indices(kmer_seq, independent_kmer=False)
    assert_allclose(got, arr)
