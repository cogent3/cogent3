import itertools

import numpy
import pytest

import cogent3.align.progressive
import cogent3.evolve.substitution_model
from cogent3 import (
    get_moltype,
    load_aligned_seqs,
    load_tree,
    make_unaligned_seqs,
)
from cogent3.align.align import (
    classic_align_pairwise,
    global_pairwise,
    local_pairwise,
    make_dna_scoring_dict,
    make_generic_scoring_dict,
)
from cogent3.evolve.models import HKY85, get_model

dna_model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
    model_gaps=False,
    equal_motif_probs=True,
)

DNA = get_moltype("dna")
PROTEIN = get_moltype("protein")
seq1 = DNA.make_seq(seq="AAACCGGACATTACGTGCGTA", name="FAKE01")
seq2 = DNA.make_seq(seq="CCGGTCAGGTTACGTACGTT", name="FAKE02")


def matched_columns(align):
    """Count the matched columns in an alignment"""

    def all_same(column):
        column = set(numpy.array(column).flatten())
        return len(column) == 1

    return len(align.filtered(all_same))


def _aligned_both_ways(seq1, seq2, **kw):
    S = make_dna_scoring_dict(10, -1, -8)
    a1 = classic_align_pairwise(seq1, seq2, S, 10, 2, **kw)
    a2 = classic_align_pairwise(seq2, seq1, S, 10, 2, **kw)
    return [a1, a2]


def test_local():
    for a in _aligned_both_ways(seq1, seq2, local=True):
        assert matched_columns(a) == 15
        assert len(a) == 19


def test_pwise_protein():
    """works for pairwise protein alignment"""
    S = make_generic_scoring_dict(1, PROTEIN)
    seq1 = PROTEIN.make_seq(seq="MAYPFQLGLQD", name="seq1")
    seq2 = PROTEIN.make_seq(seq="MAYPFGLQD", name="seq2")
    a1 = classic_align_pairwise(seq1, seq2, S, 10, 2, local=False)
    assert a1.to_dict() == {"seq1": "MAYPFQLGLQD", "seq2": "MAYPF--GLQD"}


def test_gap_at_one_end():
    for a in _aligned_both_ways(seq1, seq2, local=False):
        assert matched_columns(a) == 15
        assert len(a) == 23


def test_gaps_at_both_ends():
    s = "AAACCGGTTT"
    s1 = DNA.make_seq(seq=s[:-2], name="A")
    s2 = DNA.make_seq(seq=s[2:], name="B")
    for a in _aligned_both_ways(s1, s2, local=False):
        assert matched_columns(a) == 6
        assert len(a) == 10


def test_short():
    s1 = DNA.make_seq(seq="TACAGTA", name="A")
    s2 = DNA.make_seq(seq="TACGTC", name="B")
    for a in _aligned_both_ways(s1, s2, local=False):
        assert matched_columns(a) == 5
        assert len(a) == 7


def test_pairwise_returns_score():
    """exercise pairwise local/global returns alignment score"""
    S = make_dna_scoring_dict(10, -1, -8)
    _, score = local_pairwise(seq1, seq2, S, 10, 2, return_score=True)
    assert score > 100
    _, score = global_pairwise(seq1, seq2, S, 10, 2, return_score=True)
    assert score > 100


def test_codon():
    s1 = DNA.make_seq(seq="TACGCCGTA", name="A")
    s2 = DNA.make_seq(seq="TACGTA", name="B")
    codon_model = cogent3.evolve.substitution_model.TimeReversibleCodon(
        model_gaps=False,
        equal_motif_probs=True,
        mprob_model="conditional",
    )
    tree = cogent3.make_tree(tip_names=["A", "B"])
    lf = codon_model.make_likelihood_function(tree, aligned=False)
    lf.set_sequences({"A": s1, "B": s2})
    a = lf.get_log_likelihood().edge.get_viterbi_path().get_alignment()
    assert matched_columns(a) == 6
    assert len(a) == 9


def test_local_tiebreak():
    """Should pick the first best-equal hit rather than the last one"""
    # so that the Pyrex and Python versions give the same result.
    score_matrix = make_dna_scoring_dict(match=1, transition=-1, transversion=-1)
    pattern = DNA.make_seq(seq="CWC", name="pattern")
    two_hit = DNA.make_seq(seq="CACTC", name="target")
    aln = local_pairwise(pattern, two_hit, score_matrix, 5, 2)
    hit = aln.get_seq("target")
    assert str(hit).lower() == "cac"


def test_forward():
    tree = cogent3.make_tree(tip_names="AB")
    pc = dna_model.make_likelihood_function(tree, aligned=False)
    pc.set_sequences({"A": seq1, "B": seq2})
    LnL = pc.get_log_likelihood()
    assert isinstance(LnL, float)


def make_aln(
    orig,
    model=dna_model,
    param_vals=None,
    indel_rate=0.1,
    indel_length=0.5,
    **kw,
):
    """Helper function to create an alignment"""
    kw["indel_rate"] = indel_rate
    kw["indel_length"] = indel_length
    seqs = {key: DNA.make_seq(seq=value) for (key, value) in list(orig.items())}
    if len(seqs) == 2:
        tree = cogent3.make_tree(treestring="(A:.1,B:.1)")
    else:
        tree = cogent3.make_tree(treestring="(((A:.1,B:.1):.1,C:.1):.1,D:.1)")
    aln, tree = cogent3.align.progressive.tree_align(
        model,
        seqs,
        tree=tree,
        param_vals=param_vals,
        show_progress=False,
        **kw,
    )
    return aln


def _test_aln_standard(seqs, model=dna_model, param_vals=None, **kw):
    orig = {n: s.replace("-", "") for (n, s) in list(seqs.items())}
    aln = make_aln(orig, model=model, param_vals=param_vals, **kw)
    result = aln.to_dict()
    # assert the alignment result is correct
    assert seqs == result
    # and the moltype matches the model
    model = get_model(model)
    assert aln.moltype is model.moltype

    # assert the returned alignment has the correct parameter values in the
    # align.info object.
    if param_vals:
        for param, val in param_vals:
            assert aln.info.align_params[param] == val


def _test_aln_hirschberg(seqs, **kw):
    tmp = cogent3.align.pairwise.HIRSCHBERG_LIMIT
    try:
        cogent3.align.pairwise.HIRSCHBERG_LIMIT = 100
        _test_aln_standard(seqs, **kw)
    finally:
        cogent3.align.pairwise.HIRSCHBERG_LIMIT = tmp


@pytest.fixture(params=[_test_aln_standard, _test_aln_hirschberg])
def _test_aln(request):
    return request.param


def test_progressive1(_test_aln):
    """test progressive alignment, gaps in middle"""
    _test_aln(
        {"A": "TACAGTA", "B": "TAC-GTC", "C": "TA---TA", "D": "TAC-GTC"},
        model="F81",
    )


def test_progessive_model_name(_test_aln):
    """tree_align handles models specified by name"""
    _test_aln({"A": "TACAGTA", "B": "TAC-GTC", "C": "TA---TA", "D": "TAC-GTC"})


def test_progressive_est_tree():
    """exercise progressive alignment without a guide tree"""
    seqs = make_unaligned_seqs(
        data={
            "A": "TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA",
            "B": "TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT",
            "C": "TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT",
        },
        moltype="dna",
    )
    aln, _ = cogent3.align.progressive.tree_align(
        HKY85(),
        seqs,
        show_progress=False,
        param_vals={"kappa": 4.0},
    )

    expect = {
        "A": "TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA-------",
        "C": "TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT",
        "B": "TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT",
    }
    assert aln.to_dict() == expect

    aln, _ = cogent3.align.progressive.tree_align(
        HKY85(),
        seqs,
        show_progress=False,
        params_from_pairwise=True,
    )

    expect = {
        "A": "TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAACA-------",
        "C": "TGTGGCACAAGTACTCATGCCAGCTCAGTACAGCATGAGAACAGCAGTTT",
        "B": "TGTGGCACAGATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTT",
    }
    assert aln.to_dict() == expect


def test_align_info():
    """alignment info object has parameter values"""
    aln = make_aln(
        {"A": "GCCTCGG", "B": "GCCTCGG", "C": "GCCTCGGAAACGT", "D": "AAACGT"},
    )
    assert aln.info["align_params"]["lnL"] < 0


def test_progressive_params(_test_aln):
    """excercise progressive alignment providing model params"""
    _test_aln(
        {"A": "TACAGTA", "B": "TAC-GTC", "C": "TA---TA", "D": "CAC-CTA"},
        model=HKY85(),
        param_vals=[("kappa", 2.0)],
    )


def test_tree_align_does_pairs(_test_aln):
    """test tree_align handles pairs of sequences"""
    _test_aln({"A": "ACTTGTAC", "B": "AC--GTAC"})


def test_gap_at_start(_test_aln):
    """test progressive alignment, gaps at start"""
    _test_aln({"A": "-AC", "B": "-AC", "C": "-AC", "D": "GAC"})


def test_gap_at_end(_test_aln):
    """test progressive alignment, gaps at end"""
    _test_aln({"A": "GT-", "B": "GT-", "C": "GT-", "D": "GTA"})


def test_gaps2(_test_aln):
    """gaps have real costs, even end gaps"""
    _test_aln({"A": "G-", "B": "G-", "C": "GA", "D": "A-"})

    _test_aln({"A": "-G", "B": "-G", "C": "AG", "D": "-A"})


def test_difficult_end_gaps(_test_aln):
    _test_aln({"A": "--CCTC", "B": "--CCTC", "C": "GACCTC", "D": "GA----"})
    _test_aln(
        {
            "A": "GCCTCGG------",
            "B": "GCCTCGG------",
            "C": "GCCTCGGAAACGT",
            "D": "-------AAACGT",
        },
    )


def test_tree_align_handles_zero_lengths():
    seqs = make_unaligned_seqs(
        data={
            "A": "TTAATTTTAGTAGTGCTATCCCC",
            "B": "TTAATTTTAGTAGTGCTATCCCA",
            "C": "TTAATTTTAGTAGTGCTATCC",
        },
        moltype="dna",
    )

    expected = {
        "A": "TTAATTTTAGTAGTGCTATCCCC",
        "B": "TTAATTTTAGTAGTGCTATCCCA",
        "C": "TTAATTTTAGTAGTGCTATCC--",
    }

    tree_mapping = [
        "A: 0.0225400070648391",
        "B: 0.0225400070648391",
        "C: 0.0",
    ]
    tree_variants = itertools.permutations(tree_mapping, r=3)

    for tree_encoding in tree_variants:
        aln, _ = cogent3.align.progressive.tree_align(
            model="F81",
            seqs=seqs,
            tree=cogent3.make_tree("({},{},{})".format(*tree_encoding)),
            show_progress=False,
        )
        assert aln.to_dict() == expected


@pytest.fixture(scope="session")
def seqs(DATA_DIR):
    tree = load_tree(DATA_DIR / "brca1_5.tree")
    aln = load_aligned_seqs(DATA_DIR / "brca1.fasta", moltype="dna")
    return aln[200:1200].take_seqs(tree.get_tip_names()).degap()


@pytest.mark.xfail(reason="fails on linux due to no effect of iters")
def test_tree_align_pwise_iter(seqs):
    kwargs = {
        "model": "F81",
        "seqs": seqs,
        "show_progress": False,
        "indel_rate": 1e-3,
        "indel_length": 1e-1,
    }
    aln, _ = cogent3.align.progressive.tree_align(iters=None, **kwargs)
    one = aln.alignment_quality(app_name="sp_score", calc="pdist")
    for _ in range(10):
        aln, _ = cogent3.align.progressive.tree_align(
            iters=1,
            approx_dists=True,
            **kwargs,
        )
        two = aln.alignment_quality(app_name="sp_score", calc="pdist")
        # the quality scores will differ, but they're not deterministic
        # because the alignments are not deterministic
        if not numpy.allclose(two, one):
            break
    else:
        msg = "all attempts produced alignments with identical quality"
        raise AssertionError(msg)


def test_tree_align_dists_from_pairwise_align(seqs):
    # difficult to test reliably so only exercising use of option
    aln, tree = cogent3.align.progressive.tree_align(
        model="F81",
        seqs=seqs,
        show_progress=False,
        approx_dists=False,
    )
    assert aln


def test_tree_align_two(seqs):
    seqs = seqs.take_seqs(["Human", "NineBande"])
    aln, tree = cogent3.align.progressive.tree_align(
        model="F81",
        seqs=seqs,
        show_progress=False,
        iters=1,
        approx_dists=True,
    )
    # the tree should have equal branch lengths
    dist = set(tree.get_distances().values())
    assert len(dist) == 1
    assert isinstance(next(iter(dist)), float)
    assert len(aln) >= seqs.get_lengths().array.min()


def test_make_dna_scoring_dict():
    scoring_matrix = make_dna_scoring_dict(10, -1, -8)

    # all transitions equal
    assert (
        scoring_matrix[("A", "G")]
        == scoring_matrix[("G", "A")]
        == scoring_matrix[("C", "T")]
        == scoring_matrix[("T", "C")]
        == -1
    )

    # all transversions equal
    assert (
        scoring_matrix[("A", "T")]
        == scoring_matrix[("A", "C")]
        == scoring_matrix[("G", "T")]
        == scoring_matrix[("G", "C")]
        == scoring_matrix[("T", "A")]
        == scoring_matrix[("T", "G")]
        == scoring_matrix[("C", "A")]
        == scoring_matrix[("C", "G")]
        == -8
    )

    # all matches equal
    assert (
        scoring_matrix[("A", "A")]
        == scoring_matrix[("G", "G")]
        == scoring_matrix[("C", "C")]
        == scoring_matrix[("T", "T")]
        == 10
    )
