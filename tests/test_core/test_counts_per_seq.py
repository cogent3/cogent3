"""Regression tests for non-overlapping, per-sequence motif counts."""

import itertools
import warnings
from collections import Counter

import numpy
import pytest

from cogent3.core import alignment as c3_alignment
from cogent3.core.profile import MotifCountsArray


@pytest.fixture(params=["aligned", "unaligned"])
def make_collection(request):
    if request.param == "aligned":
        return c3_alignment.make_aligned_seqs
    return c3_alignment.make_unaligned_seqs


def _reference_counts(data, moltype, motif_length, include_ambiguity, allow_gap):
    """Use Python blocks and explicit character sets, independent of NumPy counting."""
    gaps = set() if moltype == "bytes" else {"-", "?"}
    ambiguities = {
        "dna": set("RYSWKMBDHVN?"),
        "rna": set("RYSWKMBDHVN?"),
        "protein": set("BZX?"),
        "protein_with_stop": set("BZX?"),
        "text": set(),
        "bytes": set(),
    }[moltype]
    excluded = (set() if allow_gap else gaps) | (
        set() if include_ambiguity else ambiguities
    )
    result = {}
    for name, seq in data.items():
        # The final incomplete block is discarded; overlapping windows are not used.
        blocks = (
            seq[start : start + motif_length]
            for start in range(0, len(seq) - motif_length + 1, motif_length)
        )
        result[name] = Counter(
            block for block in blocks if not excluded.intersection(block)
        )
    return result


def _assert_counts(collection, data, motif_length=1, **kwargs):
    """Check values, public labels, integer dtype, and stable row/column ordering."""
    counts = _reference_counts(
        data,
        collection.moltype.label,
        motif_length,
        kwargs.get("include_ambiguity", False),
        kwargs.get("allow_gap", False),
    )
    motifs = set().union(*(set(row) for row in counts.values()))
    if not kwargs.get("exclude_unobserved", True):
        # Construct the canonical Cartesian product independently of the kmer API.
        alphabet = collection.moltype.alphabet
        empty = type(alphabet[0])()
        motifs.update(
            empty.join(chars)
            for chars in itertools.product(alphabet, repeat=motif_length)
        )
    motifs = sorted(motifs)
    expect = {name: {motif: counts[name][motif] for motif in motifs} for name in data}
    kwargs.setdefault("exclude_unobserved", True)
    got = collection.counts_per_seq(motif_length=motif_length, **kwargs)
    assert isinstance(got, MotifCountsArray)
    assert got.motifs == tuple(motifs)
    assert got.template.names[0] == list(data)
    assert numpy.issubdtype(got.array.dtype, numpy.integer)
    assert got.to_dict() == expect


@pytest.mark.parametrize(
    ("moltype", "data"),
    [
        ("dna", {"z": "ACGTRYN?-ACG", "a": "TTTGGGAAACCC"}),
        ("rna", {"z": "ACGURYN?-ACG", "a": "UUUGGGAAACCC"}),
        ("protein", {"z": "ACDBZX?-ACDE", "a": "WWWGGGAAACCC"}),
        ("protein_with_stop", {"z": "AC*BZX?-ACDE", "a": "WWWGGGAAACCC"}),
        ("text", {"z": "AAZZ?-BCDEFG", "a": "ZZZAAAAAABBB"}),
        (
            "bytes",
            {
                "z": b"\x00\xff\x80?-ABC\xfe\x00\xffA",
                "a": b"\xff\xff\xffAAABBB\x00\x00\x00",
            },
        ),
    ],
)
@pytest.mark.parametrize("motif_length", [1, 2, 3])
@pytest.mark.parametrize(
    ("include_ambiguity", "allow_gap"), list(itertools.product([False, True], repeat=2))
)
def test_counts_match_independent_reference(
    make_collection, moltype, data, motif_length, include_ambiguity, allow_gap
):
    collection = make_collection(data, moltype=moltype)
    _assert_counts(
        collection,
        data,
        motif_length,
        include_ambiguity=include_ambiguity,
        allow_gap=allow_gap,
    )


@pytest.mark.parametrize(
    ("moltype", "motif_length", "sequence"),
    [
        ("dna", 3, "ACGACG"),
        ("rna", 2, "ACUACU"),
        ("protein", 2, "ACDACE"),
        ("protein_with_stop", 2, "AC*AC*"),
        ("text", 1, "AABBAA"),
        ("bytes", 1, b"\xff\x00AA\x80\xff"),
    ],
)
def test_unobserved_canonical_columns(make_collection, moltype, motif_length, sequence):
    data = {"z": sequence, "a": sequence}
    collection = make_collection(data, moltype=moltype)
    _assert_counts(collection, data, motif_length, exclude_unobserved=False)


@pytest.mark.parametrize("motif_length", [1, 2, 3])
def test_counts_after_reverse_complement_rename_and_reorder(
    make_collection, motif_length
):
    data = {"first": "ACGTN?-ACGTA", "middle": "AAAAAAAAAAAA", "last": "TTGCA-RYN?GG"}
    collection = make_collection(data, moltype="dna")
    # Materialize counts first to expose stale cached arrays after transformations.
    collection.counts_per_seq()
    transformed = collection.rc().renamed_seqs(str.upper).take_seqs(["LAST", "FIRST"])
    complements = str.maketrans("ACGTRY-N?", "TGCAYR-N?")
    expected = {
        name.upper(): data[name].translate(complements)[::-1]
        for name in ("last", "first")
    }
    _assert_counts(
        transformed, expected, motif_length, include_ambiguity=True, allow_gap=True
    )


@pytest.mark.parametrize(
    "selection", [slice(1, 12, 2), slice(None, None, -1), slice(10, 1, -2)]
)
@pytest.mark.parametrize("motif_length", [1, 2, 3])
def test_counts_on_alignment_slice(selection, motif_length):
    data = {"z": "ACGTN?-ACGTA", "a": "TTGCATGCAACG"}
    alignment = c3_alignment.make_aligned_seqs(data, moltype="dna")
    expected = {name: seq[selection] for name, seq in data.items()}
    if selection.step < 0:
        # Negative-step DNA views also complement the selected characters.
        complements = str.maketrans("ACGTN?-", "TGCAN?-")
        expected = {name: seq.translate(complements) for name, seq in expected.items()}
    _assert_counts(
        alignment[selection],
        expected,
        motif_length,
        include_ambiguity=True,
        allow_gap=True,
    )


@pytest.mark.parametrize("motif_length", [1, 2, 3])
def test_ragged_collection_counts(motif_length):
    data = {"long": "ACGTN?-ACGTA", "medium": "TTTCCGG", "short": "ACG", "empty": ""}
    collection = c3_alignment.make_unaligned_seqs(data, moltype="dna")
    _assert_counts(
        collection, data, motif_length, include_ambiguity=True, allow_gap=True
    )


@pytest.mark.parametrize("exclude_unobserved", [False, True])
@pytest.mark.parametrize("motif_length", [1, 3])
def test_empty_sequences_preserve_columns(
    make_collection, exclude_unobserved, motif_length
):
    collection = make_collection({"z": "", "a": ""}, moltype="dna")
    got = collection.counts_per_seq(
        motif_length=motif_length, exclude_unobserved=exclude_unobserved
    )
    # Alignment has a historical monomer fallback in moltype order; collections
    # include canonical kmers unless exclude_unobserved is requested.
    if make_collection is c3_alignment.make_aligned_seqs:
        motifs = tuple(collection.moltype)
    elif exclude_unobserved:
        motifs = tuple(sorted(collection.moltype))
    else:
        motifs = tuple(
            sorted(
                "".join(chars)
                for chars in itertools.product("ACGT", repeat=motif_length)
            )
        )
    assert isinstance(got, MotifCountsArray)
    assert got.motifs == motifs
    assert got.template.names[0] == ["z", "a"]
    assert numpy.issubdtype(got.array.dtype, numpy.integer)
    assert not got.array.any()


def test_alignment_shorter_than_motif_retains_monomer_fallback():
    alignment = c3_alignment.make_aligned_seqs({"z": "AC", "a": "TT"}, moltype="dna")
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        got = alignment.counts_per_seq(motif_length=3, warn=True)
    assert not recorded
    assert got.motifs == tuple(alignment.moltype)
    assert not got.array.any()


def test_ragged_collection_shorter_than_motif_has_zero_counts():
    # A short row must not prevent valid rows in a ragged collection being counted.
    data = {"long": "ACGTAC", "short": "A", "empty": ""}
    collection = c3_alignment.make_unaligned_seqs(data, moltype="dna")
    _assert_counts(collection, data, motif_length=3)


def test_all_filtered_motifs_preserve_container_specific_result(make_collection):
    collection = make_collection({"z": "NN--??", "a": "??NN--"}, moltype="dna")
    got = collection.counts_per_seq(motif_length=2, exclude_unobserved=True)
    if make_collection is c3_alignment.make_aligned_seqs:
        assert got is None
    else:
        assert isinstance(got, MotifCountsArray)
        assert got.motifs == tuple("ACGT")
        assert not got.array.any()


def test_alignment_truncation_warning():
    collection = c3_alignment.make_aligned_seqs(
        {"z": "ACGTA", "a": "TTGCC"}, moltype="dna"
    )
    with pytest.warns(UserWarning, match="trimmed 2") as recorded:
        _assert_counts(collection, collection.to_dict(), motif_length=3, warn=True)
    assert len(recorded) == 1


def test_collection_truncation_warnings_identify_sequences():
    data = {"z": "ACGTA", "a": "TTGC", "complete": "ACGTAC", "empty": ""}
    collection = c3_alignment.make_unaligned_seqs(data, moltype="dna")
    with pytest.warns(
        UserWarning, match="length not divisible by 3, truncating"
    ) as recorded:
        _assert_counts(collection, data, motif_length=3, warn=True)
    assert [str(w.message) for w in recorded] == [
        "z length not divisible by 3, truncating",
        "a length not divisible by 3, truncating",
    ]


def test_truncation_is_silent_by_default(make_collection):
    data = {"z": "ACGTA", "a": "TTGCC"}
    collection = make_collection(data, moltype="dna")
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        _assert_counts(collection, data, motif_length=3)
    assert not recorded


def test_observed_long_motifs_do_not_require_cartesian_alphabet(make_collection):
    # 4**40 canonical columns cannot be materialized, but only three are observed.
    data = {"z": "ACGT" * 20, "a": "A" * 40 + "T" * 40}
    collection = make_collection(data, moltype="dna")
    _assert_counts(collection, data, motif_length=40)
