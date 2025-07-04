import itertools

import numpy
import pytest
from numpy.testing import assert_allclose

from cogent3.core import alphabet as c3_alphabet
from cogent3.core import genetic_code as c3_genetic_code
from cogent3.core import moltype as c3_moltype


@pytest.mark.parametrize(
    ("num", "expect"),
    [
        (2**8 - 1, numpy.uint8),
        (2**16 - 1, numpy.uint16),
        (2**32 - 1, numpy.uint32),
        (2**64 - 1, numpy.uint64),
    ],
)
def test_get_array_type(num, expect):
    assert c3_alphabet.get_array_type(num) is expect


def test_get_array_type_invalid():
    with pytest.raises(NotImplementedError):
        c3_alphabet.get_array_type(2**64 + 1)


def test_make_char_alphabet():
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    assert len(alpha) == 4


def test_convert_alphabet_mixed_length():
    with pytest.raises(ValueError):
        c3_alphabet.convert_alphabet(src="ACG", dest="TGCA")


@pytest.mark.parametrize(
    "moltype",
    [c3_moltype.DNA, c3_moltype.ASCII],
)
def test_to_bytes(moltype):
    alpha = moltype.alphabet
    got = alpha.as_bytes()
    expect = "".join(alpha).encode("utf8")
    assert len(got) == len(alpha)
    assert got == expect


def test_to_bytes_bytes():
    alpha = c3_moltype.BYTES.alphabet
    got = alpha.as_bytes()
    assert len(got) == len(alpha)
    expect = bytes(bytearray(range(256)))
    assert got == expect
    assert chr(got[65]) == chr(65)


@pytest.mark.parametrize("words", ["AACG", ["A", "CG"], ["CC", "AA"], "", (), []])
def test_invalid_words(words):
    with pytest.raises(ValueError):
        c3_alphabet.consistent_words(words, length=1)


def test_bytes2arr():
    b2a = c3_alphabet.bytes_to_array(b"TCAG", dtype=numpy.uint8)
    got = b2a(b"AAGTA")
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


def test_arr2bytes():
    a2b = c3_alphabet.array_to_bytes(b"TCAG")
    got = a2b(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == b"AAGTA"


def test_charalphabet_missing_error():
    with pytest.raises(ValueError):
        c3_alphabet.CharAlphabet(list("TCAG"), missing=".")


@pytest.mark.parametrize(
    "seq",
    [b"AAGTA", "AAGTA", numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)],
)
def test_charalphabet_to_indices(seq):
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    got = alpha.to_indices(seq)
    assert (got == numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8)).all()


@pytest.mark.parametrize("gap", ["", "-"])
def test_charalphabet_to_indices_non_canonical(gap):
    # non-canonical characters should be beyond limit of
    # canonical
    canonical = list(f"TCAG{gap}")
    alpha = c3_alphabet.CharAlphabet(canonical, gap=gap or None)
    got = alpha.to_indices("ACNGT")
    assert got.max() >= len(canonical)


def test_charalphabet_to_indices_invalid():
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.to_indices(list("TCAG"))


def test_charalphabet_from_indices():
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    got = alpha.from_indices(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == "AAGTA"


def test_charalphabet_from_indices_invalid():
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.from_indices([0, 1])


def test_convert_alphabet():
    complement = c3_alphabet.convert_alphabet(b"ACGT", b"TGCA")
    seq = b"AAGG"
    assert complement(seq) == b"TTCC"


def test_convert_alphabet_with_delete():
    complement = c3_alphabet.convert_alphabet(b"ACGT", b"TGCA", delete=b"-")
    seq = b"AAG--G"
    assert complement(seq) == b"TTCC"


def test_charalphabet_with_gap():
    alpha = c3_alphabet.CharAlphabet(list("TCAG-"), gap="-")
    got = alpha.from_indices(numpy.array([2, 2, 3, 0, 2], dtype=numpy.uint8))
    assert got == "AAGTA"


def test_charalphabet_with_gap_missing():
    with pytest.raises(AssertionError):
        _ = c3_alphabet.CharAlphabet(list("TCAG"), gap="-")


def test_charalphabet_motif_len():
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    assert alpha.motif_len == 1


@pytest.mark.parametrize("seq", ["ACGT", numpy.array([2, 1, 3, 0], dtype=numpy.uint8)])
def test_is_valid(seq):
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    assert alpha.is_valid(seq)


@pytest.mark.parametrize(
    "seq",
    ["-ACGT", numpy.array([4, 2, 1, 3, 0], dtype=numpy.uint8)],
)
def test_not_is_valid(seq):
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    assert not alpha.is_valid(seq)


def test_is_valid_typerror():
    alpha = c3_alphabet.CharAlphabet(list("TCAG"))
    with pytest.raises(TypeError):
        _ = alpha.is_valid(list("ACGT"))


@pytest.mark.parametrize(("num_states", "ndim"), [(2, 1), (4, 2), (4, 4)])
def test_interconversion_of_coords_indices(num_states, ndim):
    # make sure match numpy functions
    coeffs = numpy.array(c3_alphabet.coord_conversion_coeffs(num_states, ndim))
    for coord in itertools.product(*(list(range(num_states)),) * ndim):
        coord = numpy.array(coord, dtype=numpy.uint64)
        nidx = numpy.ravel_multi_index(coord, dims=(num_states,) * ndim)
        idx = c3_alphabet.coord_to_index(coord, coeffs)
        assert idx == nidx
        ncoord = numpy.unravel_index(nidx, shape=(num_states,) * ndim)
        got = c3_alphabet.index_to_coord(idx, coeffs)
        assert (got == ncoord).all()
        assert (got == coord).all()


def test_coord2index_fail():
    coord = numpy.array((0, 1), dtype=numpy.uint64)
    coeffs = numpy.array(c3_alphabet.coord_conversion_coeffs(2, 4))
    # dimension of coords inconsistent with coeffs
    with pytest.raises(ValueError):
        c3_alphabet.coord_to_index(coord, coeffs)


@pytest.mark.parametrize("gap", ["", "-"])
def test_char_alphabet_with_gap(gap):
    alpha = c3_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    gapped = alpha.with_gap_motif()
    assert "-" in gapped


@pytest.mark.parametrize("k", [2, 3])
def test_kmer_alphabet_construction(k):
    dna = c3_moltype.get_moltype("dna")
    monomers = dna.alphabet
    kmers = monomers.get_kmer_alphabet(k)
    assert len(kmers) == 4**k
    assert kmers.motif_len == k


@pytest.mark.parametrize("gap_index", [-1, 4])
@pytest.mark.parametrize("indices", [2, 3, (2, 3)])
def test_seq_to_kmer_indices_handle_gap(indices, gap_index):
    k = 2
    coeffs = c3_alphabet.coord_conversion_coeffs(4, k)
    seq = numpy.array([0, 1, 1, 3], dtype=numpy.uint8)
    result = numpy.zeros(len(seq) // k, dtype=numpy.uint8)
    # the gap index is 4
    seq.put(indices, 4)
    got = c3_alphabet.seq_to_kmer_indices(
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


@pytest.mark.parametrize("gap_index", [-1, 4])
def test_seq_to_kmer_indices_handle_missing(gap_index):
    k = 2
    coeffs = c3_alphabet.coord_conversion_coeffs(4, k)
    # first dinuc has missing, last dinuc has a gap
    seq = numpy.array([5, 1, 0, 1, 0, 4], dtype=numpy.uint8)
    result = numpy.zeros(len(seq) // k, dtype=numpy.uint8)
    expect = numpy.array([17 if gap_index > 0 else 16, 1, 16], dtype=numpy.uint8)
    got = c3_alphabet.seq_to_kmer_indices(
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


@pytest.mark.parametrize("gap", ["", "-"])
def test_kmer_alpha_to_indices_to_index_with_gap(gap):
    alpha = c3_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    dinucs = alpha.get_kmer_alphabet(k=2, include_gap=bool(gap))
    got1 = int(dinucs.to_indices("--")[0])
    got2 = dinucs.to_index("--")
    assert got1 == got2 == 16


@pytest.mark.parametrize("seq", ["AC", numpy.array([2, 1], dtype=numpy.uint8)])
def test_kmer_alphabet_to_index_mixed(seq):
    dna = c3_moltype.get_moltype("dna")
    kmers = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    kmer_index = kmers.to_index(seq)
    # we compare result to tuple result
    assert kmer_index == 9


@pytest.mark.parametrize("gap", ["", "-"])
def test_kmer_index_to_seq(gap):
    alpha = c3_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    dinuc_alpha = alpha.get_kmer_alphabet(k=2, include_gap=bool(gap))
    seq = alpha.to_indices("ACTG")
    expect = numpy.array(
        [dinuc_alpha.index(d) for d in ("AC", "TG")],
        dtype=numpy.uint8,
    )
    assert (dinuc_alpha.to_indices(seq) == expect).all()


def test_kmer_index_to_gapped_seq():
    alpha = c3_alphabet.CharAlphabet(list("TCAG-"), gap="-")
    dinuc_alpha = alpha.get_kmer_alphabet(k=2, include_gap=True)
    seq = alpha.to_indices("AC--TG")
    expect = numpy.array(
        [dinuc_alpha.index(d) for d in ("AC", "--", "TG")],
        dtype=numpy.uint8,
    )
    assert (dinuc_alpha.to_indices(seq) == expect).all()


@pytest.mark.parametrize("gap", ["", "-"])
@pytest.mark.parametrize("k", [2, 3])
def test_kmer_alphabet_with_gap_motif(gap, k):
    alpha = c3_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    kmers = alpha.get_kmer_alphabet(k=k)
    gap_state = "-" * k
    gapped = kmers.with_gap_motif()
    assert gap_state in gapped
    assert gapped.gap_char == gap_state
    assert gapped.gap_index == len(gapped) - 1


@pytest.mark.parametrize("gap", ["", "-"])
def test_kmer_alpha_to_indices_with_gap(gap):
    alpha = c3_alphabet.CharAlphabet(list(f"TCAG{gap}"), gap=gap or None)
    dinucs = alpha.get_kmer_alphabet(k=2, include_gap=bool(gap))
    seq = numpy.array([5, 1, 0, 1, 0, 4], dtype=numpy.uint8)
    got = dinucs.to_indices(seq)
    expect = numpy.array([17 if gap else 16, 1, 16], dtype=numpy.uint8)
    assert (got == expect).all()


@pytest.mark.parametrize("k", [2, 3])
def test_kmer_gapped_alphabet_construction(k):
    dna = c3_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_kmer_alphabet(k, include_gap=True)
    assert len(kmers) == 1 + 4**k


def get_rand_coord_seq(k, cast):
    dna = c3_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    coord = numpy.random.randint(0, 4, size=k, dtype=numpy.uint8)
    seq = monomers.from_indices(coord)
    if cast == bytes:
        seq = seq.encode("utf8")
    elif cast is None:
        seq = coord
    return coord, seq


@pytest.mark.parametrize("cast", [None, bytes, str])
@pytest.mark.parametrize("k", [2, 3])
def test_kmer_alphabet_to_index(k, cast):
    # should match numpy ravel
    dna = c3_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_kmer_alphabet(k, include_gap=True)
    coord, seq = get_rand_coord_seq(k, cast)
    expect = numpy.ravel_multi_index(coord, dims=(4,) * k)
    assert kmers.to_index(seq) == expect


def test_kmer_alphabet_to_index_invalid():
    # should match numpy ravel
    dna = c3_moltype.get_moltype("dna")
    kmers = dna.gapped_alphabet.get_kmer_alphabet(3, include_gap=True)
    with pytest.raises(TypeError):
        kmers.to_index([0, 1, 2])


@pytest.mark.parametrize("k", [2, 3])
def test_kmer_alphabet_from_index(k):
    # should match numpy ravel
    dna = c3_moltype.get_moltype("dna")
    monomers = dna.gapped_alphabet
    kmers = monomers.get_kmer_alphabet(k, include_gap=True)
    coord = numpy.random.randint(0, 4, size=k, dtype=numpy.uint8)
    index = numpy.ravel_multi_index(coord, dims=(4,) * k)
    assert (kmers.from_index(index) == coord).all()


def test_seq_to_kmer_indices():
    numpy_func = numpy.ravel_multi_index
    dims = 4, 4, 4
    dna = c3_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    coeffs = c3_alphabet.coord_conversion_coeffs(4, 3)
    num = len(seq) // 3
    result = numpy.zeros(num, dtype=numpy.uint8)
    got = c3_alphabet.seq_to_kmer_indices(arr, result, coeffs, 4, 3)
    expect = numpy.array(
        [numpy_func(arr[i : i + 3], dims=dims) for i in range(0, 9, 3)],
    )
    assert_allclose(got, expect)


def test_kmer_alphabet_to_indices_independent():
    numpy_func = numpy.ravel_multi_index
    dims = 4, 4, 4
    dna = c3_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    got = trinuc_alpha.to_indices(arr, independent_kmer=True)
    expect = numpy.array(
        [numpy_func(arr[i : i + 3], dims=dims) for i in range(0, 9, 3)],
    )
    assert_allclose(got, expect)


def test_kmer_alphabet_to_indices_non_independent():
    numpy_func = numpy.ravel_multi_index
    dims = 4, 4, 4
    dna = c3_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    got = trinuc_alpha.to_indices(arr, independent_kmer=False)
    assert len(got) == len(seq) - 3 + 1
    expect = numpy.array(
        [numpy_func(arr[i : i + 3], dims=dims) for i in range(9 - 3 + 1)],
    )
    assert_allclose(got, expect)


@pytest.mark.parametrize("cast", [list, tuple])
def test_kmer_alphabet_to_indices_list_tuple(cast):
    seq = cast(["ACG", "TCC"])
    dna = c3_moltype.get_moltype("dna")
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    got = trinuc_alpha.to_indices(seq)
    expect = numpy.array(
        [trinuc_alpha.to_index(kmer) for kmer in seq],
        dtype=trinuc_alpha.dtype,
    )
    assert_allclose(got, expect)


def test_kmer_alphabet_to_indices_invalid():
    dna = c3_moltype.get_moltype("dna")
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    with pytest.raises(TypeError):
        trinuc_alpha.to_indices({"ACG", "TGG"})


@pytest.mark.parametrize("seq", ["ACGT", "ACYG", "NN"])
def test_kmer_alphabet_is_valid(seq):
    dna = c3_moltype.get_moltype("dna")
    dinuc_alpha = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    dinucs = dinuc_alpha.to_indices(seq)
    assert dinuc_alpha.is_valid(dinucs)


def test_kmer_alphabet_not_is_valid():
    dna = c3_moltype.get_moltype("dna")
    seq = numpy.array([2, 5, 3, 17], dtype=numpy.uint8)
    dinuc_alpha = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    assert not dinuc_alpha.is_valid(seq)


def test_kmer_alphabet_invalid_is_valid():
    dna = c3_moltype.get_moltype("dna")
    dinuc_alpha = dna.alphabet.get_kmer_alphabet(k=2, include_gap=False)
    with pytest.raises(TypeError):
        dinuc_alpha.is_valid("ACGT")


def test_kmer_alphabet_from_indices_independent():
    dna = c3_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    kmer_seq = trinuc_alpha.to_indices(arr, independent_kmer=True)
    got = trinuc_alpha.from_indices(kmer_seq)
    assert_allclose(got, arr)


def test_kmer_alphabet_from_indices_independent_with_gap():
    dna = c3_moltype.get_moltype("dna")
    seq = "AT--"
    array_seq = dna.gapped_alphabet.to_indices(seq)
    dinuc_alpha = dna.gapped_alphabet.get_kmer_alphabet(k=2, include_gap=True)
    kmer_seq = dinuc_alpha.to_indices(array_seq, independent_kmer=True)
    got = dinuc_alpha.from_indices(kmer_seq)
    assert_allclose(got, array_seq)


def test_kmer_alphabet_from_indices_not_independent():
    dna = c3_moltype.get_moltype("dna")
    seq = "ATGGGCAGA"
    arr = dna.alphabet.to_indices(seq)
    trinuc_alpha = dna.alphabet.get_kmer_alphabet(k=3, include_gap=False)
    kmer_seq = trinuc_alpha.to_indices(arr, independent_kmer=False)
    got = trinuc_alpha.from_indices(kmer_seq, independent_kmer=False)
    assert_allclose(got, arr)


@pytest.mark.parametrize("gap", ["", "-"])
@pytest.mark.parametrize("missing", ["", "?"])
def test_char_alphabet_num_canonical(gap, missing):
    alpha = c3_alphabet.CharAlphabet(
        list(f"TCAG{gap}{missing}"),
        gap=gap or None,
        missing=missing or None,
    )
    assert alpha.num_canonical == 4


@pytest.mark.parametrize("include_gap", [True, False])
@pytest.mark.parametrize("gap", ["", "-"])
@pytest.mark.parametrize("missing", ["", "?"])
def test_char_kmer_alphabet_construct_gap_state(gap, include_gap, missing):
    gap_state = "--"
    missing_state = "??"
    alpha = c3_alphabet.CharAlphabet(
        list(f"TCAG{gap}{missing}"),
        gap=gap or None,
        missing=missing or None,
    )
    dinucs = alpha.get_kmer_alphabet(k=2, include_gap=include_gap)
    assert gap_state in dinucs if include_gap and gap else gap_state not in dinucs
    assert (
        missing_state in dinucs
        if include_gap and missing
        else missing_state not in dinucs
    )


def test_kmeralpha_from_index_gap_or_missing():
    alpha = c3_alphabet.CharAlphabet(list("TCAG-?"), gap="-", missing="?")
    dinuc = alpha.get_kmer_alphabet(k=2, include_gap=True)
    # should be a gap motif
    expect = alpha.to_indices("--")
    got = dinuc.from_index(16)
    assert (got == expect).all()
    expect = alpha.to_indices("??")
    got = dinuc.from_index(17)
    assert (got == expect).all()


def test_kmeralpha_from_index_invalid():
    alpha = c3_alphabet.CharAlphabet(list("TCAG-?"), gap="-", missing="?")
    dinuc = alpha.get_kmer_alphabet(k=2, include_gap=False)
    with pytest.raises(ValueError):
        dinuc.from_index(16)


@pytest.mark.parametrize("gap", ["", "-"])
@pytest.mark.parametrize("missing", ["", "?"])
def test_serialise_charalphabet(gap, missing):
    alpha = c3_alphabet.CharAlphabet(
        list(f"TCAG{gap}{missing}"),
        gap=gap or None,
        missing=missing or None,
    )
    ser = alpha.to_rich_dict()
    deser = c3_alphabet.CharAlphabet.from_rich_dict(ser)
    assert deser == alpha


@pytest.mark.parametrize("gap", ["", "-"])
@pytest.mark.parametrize("missing", ["", "?"])
@pytest.mark.parametrize("include_gap", [True, False])
def test_serialise_kmeralphabet(gap, missing, include_gap):
    alpha = c3_alphabet.CharAlphabet(
        list(f"TCAG{gap}{missing}"),
        gap=gap or None,
        missing=missing or None,
    )
    kmers = alpha.get_kmer_alphabet(k=2, include_gap=include_gap)
    ser = kmers.to_rich_dict()
    deser = c3_alphabet.KmerAlphabet.from_rich_dict(ser)
    assert deser == kmers


@pytest.mark.parametrize(
    "alpha",
    [
        c3_moltype.DNA.alphabet,
        c3_moltype.DNA.gapped_alphabet,
        c3_moltype.DNA.gapped_alphabet.get_kmer_alphabet(k=2),
    ],
)
def test_deserialise_alphas(alpha):
    from cogent3.util.deserialise import deserialise_object

    ser = alpha.to_rich_dict()
    got = deserialise_object(ser)
    assert got == alpha

    # now from json
    ser = alpha.to_json()
    got = deserialise_object(ser)
    assert got == alpha


@pytest.fixture
def calpha():
    gc = c3_genetic_code.get_code(1)
    return c3_alphabet.SenseCodonAlphabet(
        words=gc.sense_codons,
        monomers=gc.moltype.alphabet,
    )


def test_codon_alphabet_init(calpha):
    assert len(calpha) == 61
    assert calpha.motif_len == 3


def test_codon_alphabet_index(calpha):
    assert calpha.to_index("TTT") == 0


@pytest.mark.parametrize("seq", [["TTT", "TTC"], "TTTTTC"])
def test_to_indices(calpha, seq):
    got = calpha.to_indices(seq)
    assert_allclose(got, numpy.array([0, 1], dtype=numpy.uint8))


def test_codon_alphabet_invalid_codon(calpha):
    with pytest.raises(c3_alphabet.AlphabetError):
        calpha.to_index("TGA")


def test_codon_alphabet_from_index(calpha):
    assert calpha.from_index(0) == "TTT"


def test_codon_alphabet_from_indices(calpha):
    assert calpha.from_indices(numpy.array([1, 0], dtype=numpy.uint8)) == ["TTC", "TTT"]


@pytest.mark.parametrize("seq", ["GGGAGA", numpy.array([1, 2]), numpy.array([2, 60])])
def test_codon_alphabet_is_valid(seq, calpha):
    assert calpha.is_valid(seq)


@pytest.mark.parametrize("seq", ["GGGTGA", numpy.array([-1, 2]), numpy.array([2, 66])])
def test_codon_alphabet_not_is_valid(seq, calpha):
    assert not calpha.is_valid(seq)


@pytest.mark.parametrize("as_array", [False, True])
def test_codon_alphabet_to_indices(calpha, as_array):
    seq = calpha.monomers.to_indices("GGGAAG") if as_array else "GGGAAG"
    got = calpha.to_indices(seq)
    expect = numpy.array(
        [calpha.to_index("GGG"), calpha.to_index("AAG")],
        dtype=numpy.uint8,
    )
    assert_allclose(got, expect)
    assert isinstance(got, numpy.ndarray)
    assert got.dtype == numpy.uint8


def test_codon_alphabet_with_gap_motif(calpha):
    with_gap = calpha.with_gap_motif()
    assert len(with_gap) == len(calpha) + 1
    assert with_gap.gap_char == "---"
    # on the standard genetic code, gap index is 61


def test_codon_alphabet_with_gap_motif_translation(calpha):
    with_gap = calpha.with_gap_motif()
    seq = "ATG---TAC"
    translated = with_gap.to_indices(seq)
    expected = [with_gap.to_index("ATG"), with_gap.gap_index, with_gap.to_index("TAC")]
    assert_allclose(translated, expected)


def test_codon_alphabet_missing(calpha):
    # not provided on codon alphabet
    assert calpha.missing_char is None
    assert calpha.missing_index is None


def test_codon_alphabet_serlialise_round_trip(calpha):
    from cogent3.util.deserialise import deserialise_object

    got = deserialise_object(calpha.to_json())
    assert isinstance(got, c3_alphabet.SenseCodonAlphabet)
    assert list(got) == list(calpha)
    assert calpha.to_index("TTC") == 1


@pytest.fixture(params=c3_moltype.DNA.iter_alphabets())
def dna_alpha(request):
    return request.param


def test_alphabet_moltype(dna_alpha):
    assert dna_alpha.moltype is c3_moltype.DNA


def test_alpha_no_moltype():
    alpha = c3_alphabet.CharAlphabet(list("ACG"))
    assert alpha.moltype is None


def test_codon_alphabet_moltype(calpha):
    assert calpha.moltype is c3_moltype.DNA


@pytest.mark.parametrize(
    "seq",
    [
        "GGTAC",
        b"GGTAC",
        numpy.array([3, 3, 0, 2, 1], dtype=numpy.uint8),
        tuple("GGTAC"),
    ],
)
def test_char_alphabet_to_indices_types(seq):
    dna = c3_moltype.get_moltype("dna")
    alpha = dna.alphabet
    got = alpha.to_indices(seq)
    expect = numpy.array([3, 3, 0, 2, 1], dtype=numpy.uint8)
    assert_allclose(got, expect)


@pytest.mark.parametrize(
    "seq",
    ["GGTAC", b"GGTAC", numpy.array([3, 3, 0, 2, 1], dtype=numpy.uint8)],
)
def test_char_alphabet_from_indices_types(seq):
    dna = c3_moltype.get_moltype("dna")
    alpha = dna.alphabet
    got = alpha.from_indices(seq)
    expect = "GGTAC"
    assert got == expect


@pytest.mark.parametrize("cast", [list, tuple])
def test_sensecodon_alphabet_to_indices_list_tuple(cast, calpha):
    seq = cast(["ACG", "TCC"])
    got = calpha.to_indices(seq)
    expect = numpy.array(
        [calpha.to_index(kmer) for kmer in seq],
        dtype=calpha.dtype,
    )
    assert_allclose(got, expect)


def test_get_subset():
    dna_alpha = c3_moltype.get_moltype("dna").gapped_alphabet
    got = dna_alpha.get_subset(("A", "G"))
    assert len(got) == 2
    assert got.gap_char is None
    assert got.missing_char is None
    got = dna_alpha.get_subset(("A", "G", "-"))
    assert got.gap_char == "-"
    assert len(got) == 3
    got = dna_alpha.get_subset("-", excluded=True)
    assert got.gap_char is None
    assert len(got) == 4
    assert got.moltype is c3_moltype.DNA


def test_get_subset_invalid():
    dna_alpha = c3_moltype.get_moltype("dna").alphabet
    with pytest.raises(c3_alphabet.AlphabetError):
        dna_alpha.get_subset(("A", "G", "-"))


def test_convert_seq_array_to_check_valid():
    bytes_alpha = c3_moltype.BYTES.most_degen_alphabet()
    dna_alpha = c3_moltype.DNA.most_degen_alphabet()
    # seq is invalid for DNA
    seq = numpy.array([65, 67, 71, 84], dtype=numpy.uint8)
    with pytest.raises(c3_alphabet.AlphabetError):
        dna_alpha.convert_seq_array_to(alphabet=bytes_alpha, seq=seq, check_valid=True)

    # bytes alphabet seq with these low ordinals invalid for DNA output
    seq = numpy.array([0, 1, 2, 3], dtype=numpy.uint8)
    with pytest.raises(c3_alphabet.AlphabetError):
        bytes_alpha.convert_seq_array_to(alphabet=dna_alpha, seq=seq, check_valid=True)


@pytest.mark.parametrize(
    "alpha",
    [
        c3_moltype.DNA.alphabet,
        c3_moltype.DNA.alphabet.get_kmer_alphabet(k=2),
        c3_genetic_code.DEFAULT.get_alphabet(include_stop=False),
    ],
)
def test_pickling_alphabet(alpha):
    import pickle

    ser = pickle.dumps(alpha)
    got = pickle.loads(ser)
    assert got == alpha


@pytest.mark.parametrize(
    "alpha",
    [
        c3_moltype.DNA.alphabet.get_kmer_alphabet(k=3),
        c3_genetic_code.DEFAULT.get_alphabet(include_stop=False),
    ],
)
def test_trinucs_to_indices_dtype(alpha):
    s = numpy.array([2, 0, 3] * 1000, dtype=numpy.uint8)
    small = alpha.to_indices(s[:24])
    assert len(small) == len(small.tobytes())
    full = alpha.to_indices(s)
    assert len(full) == len(full.tobytes())


@pytest.mark.parametrize(
    "alpha",
    [
        c3_moltype.DNA.alphabet,
        c3_moltype.DNA.alphabet.get_kmer_alphabet(k=2),
        c3_moltype.DNA.alphabet.get_kmer_alphabet(k=3),
        c3_moltype.DNA.alphabet.get_kmer_alphabet(k=4),
        c3_genetic_code.DEFAULT.get_alphabet(include_stop=False),
    ],
)
@pytest.mark.parametrize("cast", [numpy.array, str])
def test_consistent_dtype(alpha, cast):
    seq = c3_moltype.DNA.make_seq(seq="ATG" * 600, name="s1")
    seq = cast(seq)
    indices = alpha.to_indices(seq)
    assert indices.dtype == alpha.dtype
