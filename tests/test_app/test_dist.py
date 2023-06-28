import itertools
import pathlib

from tempfile import TemporaryDirectory
from unittest import TestCase, main

import pytest

from numpy import log, polyval
from numpy.testing import assert_allclose

from cogent3 import (
    DNA,
    PROTEIN,
    get_app,
    load_aligned_seqs,
    make_aligned_seqs,
    make_unaligned_seqs,
    open_data_store,
)
from cogent3.app.composable import WRITER
from cogent3.app.dist import (
    JACCARD_PDIST_POLY_COEFFS,
    approx_jc69,
    approx_pdist,
    jaccard_dist,
)
from cogent3.evolve.fast_distance import DistanceMatrix, HammingPair, TN93Pair
from cogent3.maths.distance_transform import jaccard


DATADIR = pathlib.Path(__file__).parent.parent / "data"

_seqs1 = {
    "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
    "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
    "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
    "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
}

_seqs2 = {
    "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
    "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
    "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
}

_seqs3 = {
    "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
    "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
}

_seqs4 = {
    "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
    "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
}

_seqs5 = {"Human": "ASSLQHENSSLLLT", "Bandicoot": "XSLMLETSSLLSN"}


@pytest.fixture(scope="function")
def _seqs1_collection():
    return make_unaligned_seqs(data=_seqs1, moltype="dna")


@pytest.fixture(scope="function")
def _seqs2_collection():
    return make_unaligned_seqs(data=_seqs2, moltype="dna")


def _get_all_composable_apps():
    out_dstore = open_data_store(":memory:", mode="w")
    return [
        get_app("align_to_ref"),
        get_app("progressive_align", model="GY94"),
        get_app("fixed_length", 100),
        get_app("sample.min_length", 100),
        get_app("write_seqs", out_dstore),
        get_app(
            "omit_bad_seqs",
        ),
        get_app(
            "omit_degenerates",
        ),
        get_app("take_codon_positions", 1),
        get_app(
            "take_named_seqs",
        ),
        get_app("trim_stop_codons", gc=1),
    ]


class FastSlowDistTests(TestCase):
    seqs1 = make_unaligned_seqs(_seqs1, moltype=DNA)
    seqs2 = make_unaligned_seqs(_seqs2, moltype=DNA)
    seqs3 = make_unaligned_seqs(_seqs3, moltype=DNA)
    seqs4 = make_unaligned_seqs(_seqs4, moltype=DNA)
    seqs5 = make_unaligned_seqs(_seqs5, moltype=PROTEIN)

    def test_init(self):
        """tests if fast_slow_dist can be initialised correctly"""

        fast_slow_dist = get_app("fast_slow_dist", fast_calc="hamming", moltype="dna")
        self.assertIsInstance(fast_slow_dist.fast_calc, HammingPair)
        self.assertIsNone(fast_slow_dist._sm)

        fast_slow_dist = get_app("fast_slow_dist", distance="TN93")
        self.assertIsInstance(fast_slow_dist.fast_calc, TN93Pair)
        self.assertEqual(fast_slow_dist._sm.name, "TN93")
        fast_slow_dist = get_app("fast_slow_dist", distance="GTR")
        self.assertEqual(fast_slow_dist._sm.name, "GTR")

        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        self.assertEqual(fast_slow_dist._sm.name, "TN93")
        self.assertIsNone(fast_slow_dist.fast_calc)

        with self.assertRaises(ValueError):
            fast_slow_dist = get_app(
                "fast_slow_dist", distance="TN93", fast_calc="TN93", slow_calc="TN93"
            )

        with self.assertRaises(ValueError):
            fast_slow_dist = get_app("fast_slow_dist", fast_calc="GTR")

        with self.assertRaises(ValueError):
            fast_slow_dist = get_app("fast_slow_dist", slow_calc="hamming")

    def test_compatible_parameters(self):
        """tests if the input parameters are compatible with fast_slow_dist initialisation"""
        for kwargs in (
            dict(fast_calc="hamming", moltype="dna"),
            dict(fast_calc="TN93"),
            dict(slow_calc="GTR"),
            dict(fast_calc="TN93"),
        ):
            _ = get_app("fast_slow_dist", **kwargs)

    def test_incompatible_parameters(self):
        """tests incompatible input parameters with fast_slow_dist initialisation"""
        for kwargs in (
            dict(fast_calc="hamming"),
            dict(slow_calc="paralinear"),
            dict(fast_calc="GTR"),
            dict(slow_calc="hamming", moltype="dna"),
        ):
            with self.assertRaises(ValueError):
                _ = get_app("fast_slow_dist", **kwargs)

    def test_composable_apps(self):
        """tests two composable apps"""
        composable_apps = _get_all_composable_apps()
        calc_dist = get_app("fast_slow_dist", fast_calc="hamming", moltype="dna")
        for app in composable_apps:
            if app.app_type is WRITER:
                # cannot have a WRITER before a GENERIC
                continue
            # Compose two composable applications, there should not be exceptions.
            got = app + calc_dist
            self.assertIsInstance(got, type(calc_dist))
            self.assertIs(got.input, app)
            self.assertIsInstance(got._data_types, frozenset)
            self.assertIsInstance(got._return_types, frozenset)
            self.assertIs(got.input, app)
            app.disconnect()
            calc_dist.disconnect()

    def test_est_dist_pair_slow(self):
        """tests the distance between seq pairs in aln"""
        aligner = get_app(
            "align_to_ref",
        )
        aln3 = aligner(self.seqs3)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Human")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Mouse")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Mouse", "Human")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = get_app(
            "align_to_ref",
        )
        aln3 = aligner(self.seqs4)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Human")
        aln3 = aligner(self.seqs4)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Opossum")
        aln3 = aligner(self.seqs4)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        # now as a process
        proc = get_app(
            "align_to_ref",
        ) + get_app("fast_slow_dist", fast_calc="hamming", moltype="dna")
        got = proc(self.seqs1)
        self.assertEqual(got[("Human", "Rhesus")], 1)

        treestring = "(Human:0.2,Bandicoot:0.2)"
        aligner = get_app("progressive_align", model="WG01", guide_tree=treestring)
        _ = aligner(self.seqs5)

    def test_composes_with_write_tabular(self):
        """correctly links to tabular"""
        with TemporaryDirectory(dir=".") as dirname:
            out_dstore = open_data_store(dirname, suffix="tsv", mode="w")
            writer = get_app("write_tabular", out_dstore)
            dist_calc = get_app("fast_slow_dist", distance="hamming", moltype="protein")
            _ = dist_calc + writer

    def test_functions_as_composable(self):
        """works as a composable app"""
        from pathlib import Path

        loader = get_app("load_aligned", moltype="dna", format="paml")
        dist = get_app("fast_slow_dist", "hamming", moltype="dna")
        with TemporaryDirectory(dir=".") as dirname:
            dirname = Path(dirname)
            out_dstore = open_data_store(dirname, suffix="tsv", mode="w")
            writer = get_app("write_tabular", out_dstore)
            proc = loader + dist + writer
            _ = proc("data/brca1_5.paml")
            output = dirname / "brca1_5.tsv"
            self.assertTrue(output.exists())


if __name__ == "__main__":
    main()


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_jaccard_dist(moltype):
    """jaccard_dist app should work for the simple case

    ("s1", "ACGTA"),
    ("s2", "----C"),

    with k=2
    s1 kmers = "AC", "CG", "GT", "TA"
    s2 kmers = "AC", "CG", "GT", "TC"

    J(A,B) = 1 - |A ∩ B| / |A ∪ B|

    J(s1, s2) = 1 - |{"AC", "CG", "GT"}| / |{"AC", "CG", "GT", "TA", "TC"}|
    J(s1, s2) = 1 - 3 / 5
    J(s1, s2) = 0.4
    """
    data = dict([("s1", "ACGTA"), ("s2", "ACGTC")])
    collection = make_unaligned_seqs(data=data, moltype=moltype)

    jdist_k2 = jaccard_dist(k=2)
    dists = jdist_k2(collection)

    assert dists[("s1", "s2")] == 0.4
    assert dists[("s2", "s1")] == 0.4
    assert dists[("s1", "s1")] == 0.0
    assert dists[("s2", "s2")] == 0.0


def test_approx_pdist():
    """approx_pdist should work for the simple case

    y = polyval(JACCARD_PDIST_POLY_COEFFS, x)
    """

    data = dict(
        [
            (("s1", "s1"), 0.0),
            (("s1", "s2"), 0.4),
            (("s2", "s1"), 0.4),
            (("s2", "s2"), 0.0),
        ]
    )
    dm = DistanceMatrix(data)

    pdist_app = approx_pdist()
    pdists = pdist_app(dm)

    expect_diff = polyval(JACCARD_PDIST_POLY_COEFFS, 0.4)
    expect_same = polyval(JACCARD_PDIST_POLY_COEFFS, 0.0)

    assert pdists[("s1", "s2")] == expect_diff
    assert pdists[("s2", "s1")] == expect_diff
    assert pdists[("s1", "s1")] == expect_same
    assert pdists[("s2", "s2")] == expect_same


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_approx_jc69(moltype):
    """approx_jc69 should work the same as exact jc69 when given exact pdist"""
    seq_data = dict([("s1", "ACGAA"), ("s2", "ACGAC")])
    aln = make_aligned_seqs(data=seq_data, moltype=moltype)
    expected = aln.distance_matrix(calc="jc69")

    data = dict(
        [
            (("s1", "s1"), 0.0),
            (("s1", "s2"), 1 / 5),
            (("s2", "s1"), 1 / 5),
            (("s2", "s2"), 0.0),
        ]
    )

    dm = DistanceMatrix(data)
    jc_dist_app = approx_jc69()
    got = jc_dist_app(dm)

    assert got[("s1", "s2")] == expected[("s1", "s2")]
    assert got[("s2", "s1")] == expected[("s2", "s1")]
    assert got[("s1", "s1")] == expected[("s1", "s1")]
    assert got[("s2", "s2")] == expected[("s2", "s2")]


@pytest.mark.parametrize("moltype", ("dna", "rna"))
def test_approx_pdist_same_diff(moltype):
    """comparisons between seqs with the same position different should be equal.
    comparison between seqs with more positions different should yield a higher
    measure than comparisons between seqs with fewer positions different.

    NOTE: the coefficients used in Jaccard to Pdist fit
        were generated using k=10, here I used k=3

    ("s1", "ACGTA"),
    ("s2", "----C"),
    ("s3", "----T"),
    ("s4", "---AT"),
    """

    data = dict(
        [
            ("s1", "ACGTA"),
            ("s2", "ACGTC"),
            ("s3", "ACGTT"),
            ("s4", "ACGAT"),
        ]
    )
    pdist_app = jaccard_dist(k=3) + approx_pdist()
    collection = make_unaligned_seqs(data=data, moltype=moltype)
    dists = pdist_app(collection)

    # comparisons with one position different should be smaller than those with two
    assert dists[("s1", "s2")] < dists[("s1", "s4")]
    assert dists[("s1", "s3")] < dists[("s1", "s4")]
    assert dists[("s1", "s2")] < dists[("s1", "s4")]
    assert dists[("s2", "s3")] < dists[("s2", "s4")]

    # both (s1 and s2) and (s2 and s3) have the same position different
    assert dists[("s1", "s2")] == dists[("s1", "s3")]


def test_jaccard_dist_vals(_seqs1_collection):
    """values in the DistanceMatrix should match individually calculating the jaccard
    distance for pairs of sequence.
    """
    seqs = _seqs1_collection
    jaccard_dist_app = jaccard_dist(k=10)
    jdists = jaccard_dist_app(seqs)
    names = jdists.names

    for i, j in itertools.combinations(range(len(names)), 2):
        seq1, seq2 = names[i], names[j]
        got = jdists[(seq1, seq2)]
        expect = jaccard(
            set(seqs.get_seq(seq1).get_kmers(k=10, strict=True)),
            set(seqs.get_seq(seq2).get_kmers(k=10, strict=True)),
        )
        assert_allclose(got, expect)


def test_approx_pdist_vals(_seqs1_collection):
    """values in the DistanceMatrix should match individually calculating the pdist
    for pairs of sequence.

    testing integration of jaccard_dist() + approx_pdist() is identical to
    step-by-step calculation
    """

    seqs = _seqs1_collection

    jaccard_dist_app = jaccard_dist(k=10)
    jdists = jaccard_dist_app(seqs)

    pdist_app = jaccard_dist(k=10) + approx_pdist()
    pdists = pdist_app(seqs)
    names = pdists.names

    for i, j in itertools.combinations(range(len(names)), 2):
        seq1, seq2 = names[i], names[j]
        got = pdists[(seq1, seq2)]
        expect = polyval(JACCARD_PDIST_POLY_COEFFS, jdists[(seq1, seq2)])
        assert got == expect


def test_approx_jc69_vals(_seqs1_collection):
    """values in the DistanceMatrix should match individually calculating the jc distance
    for pairs of sequence.

    testing integration of jaccard_dist() + approx_pdist() + approx_jc69() is identical to
    step-by-step calculation
    """

    seqs = _seqs1_collection
    jaccard_dist_app = jaccard_dist(k=10)
    jdists = jaccard_dist_app(seqs)
    names = jdists.names

    pdist_app = approx_pdist()
    pdists = pdist_app(jdists)

    jc_app = jaccard_dist(k=10) + approx_pdist() + approx_jc69()
    jc_dists = jc_app(seqs)

    for i, j in itertools.combinations(range(len(names)), 2):
        seq1, seq2 = names[i], names[j]
        got = jc_dists[(seq1, seq2)]
        expect = -3 / 4 * log(1 - 4 / 3 * pdists[(seq1, seq2)])
        assert got == expect


def test_symmetry_of_dists():
    """distances are symmetric"""
    seqs = load_aligned_seqs(DATADIR / "primate_brca1.fasta", moltype="dna")
    dists = seqs.distance_matrix(calc="pdist")
    app = approx_jc69()
    got = app(dists)
    assert_allclose(got.array, got.array.T)


def test_gap_dist():
    app = get_app("gap_dist", gap_insert=10, gap_extend=1)
    # two sequences share a gap
    data = {
        "a": "TG----AATATGT------GAAAGAG",
        "b": "TTGAAGAATATGT------GAAAGAG",
        "c": "CTGAAGAACCTGTGAAAGTGAAAGAG",
    }
    aln = make_aligned_seqs(data, moltype="dna", array_align=True)
    expect = {
        ("a", "b"): 14.0,  # one gap diff of size 4
        ("a", "c"): 30.0,
        ("b", "c"): 16.0,
    }
    expect = DistanceMatrix(expect)
    dmat = app.main(aln)
    assert dmat.to_dict() == expect.to_dict()

    # shared gap actually not shared, 3 events
    data = {
        "a": "TG----AATATGTA-----GAAAGAG",
        "b": "TTGAAGAATATGTA------AAAGAG",
        "c": "CTGAAGAACCTGTGAAAGTGAAAGAG",
    }
    aln = make_aligned_seqs(data, moltype="dna", array_align=True)
    expect = {
        ("a", "b"): 45.0,  # 3 gaps diff of size 15
        ("a", "c"): 29.0,
        ("b", "c"): 16.0,
    }
    expect = DistanceMatrix(expect)
    dmat = app.main(aln)
    assert dmat.to_dict() == expect.to_dict()

    # additional gaps on either side of shared gap is two events
    data = {
        "a": "G--AG----A",
        "b": "TGGAGT--GA",
        "c": "TGGAGTGTGA",
    }
    aln = make_aligned_seqs(data, moltype="dna", array_align=True)
    expect = {
        ("a", "b"): 38,  # 3 gaps diff of size 8
        ("a", "c"): 26.0,
        ("b", "c"): 12.0,
    }
    expect = DistanceMatrix(expect)
    dmat = app.main(aln)
    assert dmat.to_dict() == expect.to_dict()
    data = {"a": "AAGAA-A", "b": "-ATAATG", "c": "C-TGG-G"}
    aln = make_aligned_seqs(data, moltype="dna", array_align=True)
    expect = {
        ("a", "b"): 22.0,  # 2 gaps diff of size 2
        ("a", "c"): 11.0,
        ("b", "c"): 33.0,  # 3 gaps diff of size 2
    }
    expect = DistanceMatrix(expect)
    dmat = app.main(aln)
    assert dmat.to_dict() == expect.to_dict()
