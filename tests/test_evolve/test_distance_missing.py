"""Pairwise deletion must not turn zero distance into sequence identity."""

import itertools
from fractions import Fraction

import numpy
import pytest
from numpy.testing import assert_allclose

import cogent3
from cogent3.evolve.fast_distance import get_distance_calculator


def _pair_counts(seq1, seq2, canonical):
    pairs = [
        (a, b)
        for a, b in zip(seq1, seq2, strict=True)
        if a in canonical and b in canonical
    ]
    return len(pairs), sum(a != b for a, b in pairs)


@pytest.mark.parametrize(
    ("moltype", "seqs", "canonical"),
    [
        ("dna", {"a": "A-", "b": "AA", "c": "AC"}, "ACGT"),
        ("dna", {"a": "AN", "b": "AA", "c": "AC"}, "ACGT"),
        ("rna", {"a": "U?", "b": "UU", "c": "UC"}, "ACGU"),
        ("protein", {"a": "M-", "b": "MM", "c": "MP"}, "ACDEFGHIKLMNPQRSTVWY"),
        ("protein", {"a": "MX", "b": "MM", "c": "MP"}, "ACDEFGHIKLMNPQRSTVWY"),
    ],
)
@pytest.mark.parametrize("method", ["pdist", "hamming"])
def test_pairwise_deletion_order_invariance(moltype, seqs, canonical, method):
    for names in itertools.permutations(seqs):
        aln = cogent3.make_aligned_seqs(
            {name: seqs[name] for name in names}, moltype=moltype
        )
        calc = get_distance_calculator(method, moltype=moltype)
        calc(aln)
        assert calc.duplicated is None
        for include_duplicates in (False, True):
            dists = calc.get_pairwise_distances(include_duplicates=include_duplicates)
            for a, b in itertools.combinations(seqs, 2):
                length, differences = _pair_counts(seqs[a], seqs[b], canonical)
                proportion = Fraction(differences, length)
                expected = float(proportion) if method == "pdist" else differences
                assert dists[a, b] == expected
                assert calc.lengths[a, b] == length
                assert calc.proportions[a, b] == float(proportion)


@pytest.mark.parametrize("method", ["pdist", "hamming", "jc69"])
def test_fast_slow_dist_pairwise_deletion(method):
    seqs = {"a": "A-", "b": "AA", "c": "AC"}
    for names in itertools.permutations(seqs):
        aln = cogent3.make_aligned_seqs(
            {name: seqs[name] for name in names}, moltype="dna"
        )
        key = "fast_calc" if method == "jc69" else "distance"
        app = cogent3.get_app("fast_slow_dist", moltype="dna", **{key: method})
        got = app(aln)
        for a, b in itertools.combinations(seqs, 2):
            length, differences = _pair_counts(seqs[a], seqs[b], "ACGT")
            p = float(Fraction(differences, length))
            expected = {
                "pdist": p,
                "hamming": differences,
                "jc69": -0.75 * numpy.log(1 - 4 * p / 3),
            }[method]
            assert_allclose(got[a, b], expected)
        # Positive control: the separate newer implementation already computes
        # these distances correctly and is not part of this repair.
        modern = aln.distance_matrix(calc="pdist", parallel=False)
        assert modern["b", "c"] == 0.5


@pytest.mark.parametrize(
    "method", ["pdist", "hamming", "jc69", "tn93", "paralinear", "logdet"]
)
@pytest.mark.parametrize("moltype", ["dna", "rna", "protein"])
def test_zero_shared_site_differences_are_not_duplicates(method, moltype):
    if moltype == "protein" and method in {"jc69", "tn93"}:
        pytest.skip("nucleotide-only calculator")
    aln = cogent3.make_aligned_seqs({"a": "A-", "b": "AA"}, moltype=moltype)
    calc = get_distance_calculator(method, moltype=moltype)
    calc(aln)
    assert calc.duplicated is None
    assert calc.dists["a", "b"] == 0.0
    assert calc.lengths["a", "b"] == 1
    assert calc.proportions["a", "b"] == 0.0
    assert calc.variances["a", "b"] == 0.0
    assert calc.stderr["a", "b"] == 0.0


@pytest.mark.parametrize(
    "method", ["pdist", "hamming", "jc69", "tn93", "paralinear", "logdet"]
)
def test_no_shared_sites_are_invalid_not_duplicates(method):
    aln = cogent3.make_aligned_seqs({"a": "--", "b": "NN", "c": "AC"}, moltype="dna")
    calc = get_distance_calculator(method, moltype="dna")
    with numpy.errstate(invalid="ignore", divide="ignore"):
        calc(aln)
    assert calc.duplicated is None
    for a, b in itertools.combinations(aln.names, 2):
        assert numpy.isnan(calc.dists[a, b])
    calc = get_distance_calculator(method, moltype="dna", invalid_raises=True)
    with (
        numpy.errstate(invalid="ignore", divide="ignore"),
        pytest.raises(ArithmeticError),
    ):
        calc(aln)


@pytest.mark.parametrize(
    "method", ["pdist", "hamming", "jc69", "tn93", "paralinear", "logdet"]
)
@pytest.mark.parametrize("sequence", ["ACGTACGT", "ACGT-CGT", "ACGTNCGT"])
def test_genuine_duplicate_reuse_is_preserved(method, sequence):
    seqs = {"a": sequence, "b": sequence, "c": "ACGAACGT"}
    aln = cogent3.make_aligned_seqs(seqs, moltype="dna")
    calc = get_distance_calculator(method, moltype="dna")
    calc(aln)
    assert calc.duplicated == {"a": ["b"]}
    all_dists = calc.get_pairwise_distances(include_duplicates=True)
    assert all_dists["a", "b"] == 0.0
    assert all_dists["a", "c"] == all_dists["b", "c"]
    assert set(calc.get_pairwise_distances(include_duplicates=False).names) == {
        "a",
        "c",
    }
    for name in ("lengths", "proportions", "variances", "stderr"):
        if method in {"hamming", "pdist"} and name == "stderr":
            # Their nonzero-distance variance is not implemented.
            continue
        if method == "logdet" and name == "stderr":
            continue
        stats = getattr(calc, name)
        assert stats["a", "c"] == stats["b", "c"]
