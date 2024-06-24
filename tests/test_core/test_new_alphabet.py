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


@pytest.mark.parametrize(
    "moltype", (new_moltype.DNA, new_moltype.DNA, new_moltype.ASCII)
)
def test_to_bytes(moltype):
    alpha = moltype.alphabet
    got = alpha.as_bytes()
    expect = "".join(alpha).encode("utf8")
    assert len(got) == len(alpha)
    assert got == expect


def test_to_bytes_bytes():
    alpha = new_moltype.BYTES.alphabet
    got = alpha.as_bytes()
    assert len(got) == len(alpha)
    expect = bytes(bytearray(range(256)))
    assert got == expect
    assert chr(got[65]) == chr(65)


@pytest.mark.parametrize("words", ("AACG", ["A", "CG"], ["CC", "AA"], "", (), []))
def test_invalid_words(words):
    with pytest.raises(ValueError):
        new_alphabet.consistent_words(words, length=1)


def test_bytes2arr():
    b2a = new_alphabet.bytes_to_array(b"TCAG", dtype=numpy.uint8)
    got = b2a(b"AAGTA")
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


def test_arr2bytes():
    a2b = new_alphabet.array_to_bytes(b"TCAG")
    got = a2b(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == b"AAGTA"


@pytest.mark.parametrize(
    "seq", (b"AAGTA", "AAGTA", numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
)
def test_charalphabet_to_indices(seq):
    alpha = new_alphabet.CharAlphabet(list("TCAG"))
    got = alpha.to_indices(seq)
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


@pytest.mark.parametrize("gap", ("", "-"))
def test_charalphabet_to_indices_non_canonical(gap):
    # non-canonical characters should be beyond limit of
    # canonical
    canonical = list(f"TCAG{gap}")
    alpha = new_alphabet.CharAlphabet(canonical, gap=gap or None)
    got = alpha.to_indices("ACNGT")
    assert got.max() >= len(canonical)


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


@pytest.mark.parametrize("gap", ("", "-"))
def test_char_alphabet_with_gap(gap):
    alpha = new_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    gapped = alpha.with_gap_motif()
    assert "-" in gapped


@pytest.mark.parametrize("k", (2, 3))
def test_kmer_alphabet_construction(k):
    dna = new_moltype.get_moltype("dna")
    monomers = dna.alphabet
    kmers = monomers.get_kmer_alphabet(k)
    assert len(kmers) == 4**k


@pytest.mark.parametrize("gap_index", (-1, 4))
@pytest.mark.parametrize("indices", (2, 3, (2, 3)))
def test_seq_to_kmer_indices_handle_gap(indices, gap_index):
    k = 2
    coeffs = new_alphabet.coord_conversion_coeffs(4, k)
    seq = numpy.array([0, 1, 1, 3], dtype=numpy.uint8)
    result = numpy.zeros(len(seq) // k, dtype=numpy.uint8)
    # the gap index is 4
    seq.put(indices, 4)
    got = new_alphabet.seq_to_kmer_indices(
        seq,
        result,
        coeffs,
        4,
        k,
        gap_char_index=gap_index,
        gap_index=-1 if gap_index == -1 else 16,
        independent_kmer=True,
    )
    expect = numpy.array([1, 16], dtype=numpy.uint8)
    assert (got == expect).all()


@pytest.mark.parametrize("gap_index", (-1, 4))
def test_seq_to_kmer_indices_handle_missing(gap_index):
    k = 2
    coeffs = new_alphabet.coord_conversion_coeffs(4, k)
    # first dinuc has missing, last dinuc has a gap
    seq = numpy.array([5, 1, 0, 1, 0, 4], dtype=numpy.uint8)
    result = numpy.zeros(len(seq) // k, dtype=numpy.uint8)
    expect = numpy.array([17 if gap_index > 0 else 16, 1, 16], dtype=numpy.uint8)
    got = new_alphabet.seq_to_kmer_indices(
        seq,
        result,
        coeffs,
        4,
        k,
        gap_char_index=gap_index,
        gap_index=-1 if gap_index == -1 else 16,
        independent_kmer=True,
    )
    assert (got == expect).all()


@pytest.mark.parametrize("gap", ("", "-"))
def test_kmer_alpha_to_indices_to_index_with_gap(gap):
    alpha = new_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    dinucs = alpha.get_kmer_alphabet(k=2, include_gap=bool(gap))
    got1 = int(dinucs.to_indices("--")[0])
    got2 = dinucs.to_index("--")
    assert got1 == got2 == 16


@pytest.mark.parametrize("seq", ("AC", numpy.array([2, 1], dtype=numpy.uint8)))
def test_kmer_alphabet_to_index_mixed(seq):
    dna = new_moltype.get_moltype("dna")
    kmers = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    kmer_index = kmers.to_index(seq)
    # we compare result to tuple result
    assert kmer_index == 9


@pytest.mark.parametrize("gap", ("", "-"))
def test_kmer_index_to_seq(gap):
    alpha = new_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    dinuc_alpha = alpha.get_kmer_alphabet(k=2, include_gap=bool(gap))
    seq = alpha.to_indices("ACTG")
    expect = numpy.array(
        [dinuc_alpha.index(d) for d in ("AC", "TG")], dtype=numpy.uint8
    )
    assert (dinuc_alpha.to_indices(seq) == expect).all()


def test_kmer_index_to_gapped_seq():
    alpha = new_alphabet.CharAlphabet(list("TCAG-"), gap="-")
    dinuc_alpha = alpha.get_kmer_alphabet(k=2, include_gap=True)
    seq = alpha.to_indices("AC--TG")
    expect = numpy.array(
        [dinuc_alpha.index(d) for d in ("AC", "--", "TG")], dtype=numpy.uint8
    )
    assert (dinuc_alpha.to_indices(seq) == expect).all()


@pytest.mark.parametrize("gap", ("", "-"))
@pytest.mark.parametrize("k", (2, 3))
def test_kmer_alphabet_with_gap_motif(gap, k):
    alpha = new_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    kmers = alpha.get_kmer_alphabet(k=k)
    gap_state = "-" * k
    gapped = kmers.with_gap_motif()
    assert gap_state in gapped
    assert gapped.gap_char == gap_state
    assert gapped.gap_index == len(gapped) - 1


@pytest.mark.parametrize("gap", ("", "-"))
def test_kmer_alpha_to_indices_with_gap(gap):
    alpha = new_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    dinucs = alpha.get_kmer_alphabet(k=2, include_gap=bool(gap))
    seq = numpy.array([5, 1, 0, 1, 0, 4], dtype=numpy.uint8)
    got = dinucs.to_indices(seq)
    expect = numpy.array([17 if gap else 16, 1, 16], dtype=numpy.uint8)
    assert (got == expect).all()


@pytest.mark.parametrize("k", (2, 3))
def test_kmer_gapped_alphabet_construction(k):
    dna = new_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_kmer_alphabet(k, include_gap=True)
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
def test_kmer_alphabet_to_index(k, cast):
    # should match numpy ravel
    dna = new_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_kmer_alphabet(k, include_gap=True)
    coord, seq = get_rand_coord_seq(k, cast)
    expect = numpy.ravel_multi_index(coord, dims=(4,) * k)
    assert kmers.to_index(seq) == expect


def test_kmer_alphabet_to_index_invalid():
    # should match numpy ravel
    dna = new_moltype.get_moltype("dna")
    kmers = dna.gapped_alphabet.get_kmer_alphabet(3, include_gap=True)
    with pytest.raises(TypeError):
        kmers.to_index([0, 1, 2])


@pytest.mark.parametrize("k", (2, 3))
def test_kmer_alphabet_from_index(k):
    # should match numpy ravel
    dna = new_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_kmer_alphabet(k, include_gap=True)
    coord = numpy.random.randint(0, 4, size=k, dtype=numpy.uint8)
    index = numpy.ravel_multi_index(coord, dims=(4,) * k)
    assert (kmers.from_index(index) == coord).all()


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
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
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
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    got = trinuc_alpha.to_indices(arr, independent_kmer=False)
    assert len(got) == len(seq) - 3 + 1
    expect = numpy.array(
        [numpy_func(arr[i : i + 3], dims=dims) for i in range(0, 9 - 3 + 1)]
    )
    assert_allclose(got, expect)


def test_kmer_alphabet_to_indices_invalid():
    dna = new_moltype.get_moltype("dna")
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    with pytest.raises(TypeError):
        trinuc_alpha.to_indices(["ACG", "TGG"])


@pytest.mark.parametrize("seq", ("ACGT", "ACYG", "NN"))
def test_kmer_alphabet_is_valid(seq):
    dna = new_moltype.get_moltype("dna")
    dinuc_alpha = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    dinucs = dinuc_alpha.to_indices(seq)
    assert dinuc_alpha.is_valid(dinucs)


def test_kmer_alphabet_not_is_valid():
    dna = new_moltype.get_moltype("dna")
    seq = numpy.array([2, 5, 3, 17], dtype=numpy.uint8)
    dinuc_alpha = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    assert not dinuc_alpha.is_valid(seq)


def test_kmer_alphabet_invalid_is_valid():
    dna = new_moltype.get_moltype("dna")
    dinuc_alpha = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    with pytest.raises(TypeError):
        dinuc_alpha.is_valid("ACGT")


def test_kmer_alphabet_from_indices_independent():
    dna = new_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    kmer_seq = trinuc_alpha.to_indices(arr, independent_kmer=True)
    got = trinuc_alpha.from_indices(kmer_seq)
    assert_allclose(got, arr)


def test_kmer_alphabet_from_indices_independent_with_gap():
    dna = new_moltype.get_moltype("dna")
    seq = "AT--"
    array_seq = dna.gapped_alphabet.to_indices(seq)
    dinuc_alpha = dna.gapped_alphabet.get_kmer_alphabet(k=2, include_gap=True)
    kmer_seq = dinuc_alpha.to_indices(array_seq, independent_kmer=True)
    got = dinuc_alpha.from_indices(kmer_seq)
    assert_allclose(got, array_seq)


def test_kmer_alphabet_from_indices_not_independent():
    dna = new_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    kmer_seq = trinuc_alpha.to_indices(arr, independent_kmer=False)
    got = trinuc_alpha.from_indices(kmer_seq, independent_kmer=False)
    assert_allclose(got, arr)


@pytest.mark.parametrize("gap", ("", "-"))
@pytest.mark.parametrize("missing", ("", "?"))
def test_char_alphabet_num_canonical(gap, missing):
    alpha = new_alphabet.CharAlphabet(
        list(f"TCAG{gap}{missing}"), gap=gap or None, missing=missing or None
    )
    assert alpha.num_canonical == 4


@pytest.mark.parametrize("include_gap", (True, False))
@pytest.mark.parametrize("gap", ("", "-"))
@pytest.mark.parametrize("missing", ("", "?"))
def test_char_kmer_alphabet_construct_gap_state(gap, include_gap, missing):
    gap_state = "--"
    missing_state = "??"
    alpha = new_alphabet.CharAlphabet(
        list(f"TCAG{gap}{missing}"), gap=gap or None, missing=missing or None
    )
    dinucs = alpha.get_kmer_alphabet(k=2, include_gap=include_gap)
    assert gap_state in dinucs if include_gap and gap else gap_state not in dinucs
    assert (
        missing_state in dinucs
        if include_gap and missing
        else missing_state not in dinucs
    )


def test_kmeralpha_from_index_gap_or_missing():
    alpha = new_alphabet.CharAlphabet(list("TCAG-?"), gap="-", missing="?")
    dinuc = alpha.get_kmer_alphabet(k=2, include_gap=True)
    # should be a gap motif
    expect = alpha.to_indices("--")
    got = dinuc.from_index(16)
    assert (got == expect).all()
    expect = alpha.to_indices("??")
    got = dinuc.from_index(17)
    assert (got == expect).all()


def test_kmeralpha_from_index_invalid():
    alpha = new_alphabet.CharAlphabet(list("TCAG-?"), gap="-", missing="?")
    dinuc = alpha.get_kmer_alphabet(k=2, include_gap=False)
    with pytest.raises(ValueError):
        dinuc.from_index(16)
