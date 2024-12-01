import itertools

import numpy
import pytest
from numpy.testing import assert_allclose

from cogent3 import get_moltype, make_aligned_seqs
from cogent3.evolve import coevolution as c3_coevo
from cogent3.maths.stats.number import CategoryCounter


def calc_mi_pair(alignment, pos1, pos2, normalised=False):
    """Calculate mutual information between two positions in an alignment"""
    num_states = len(alignment.moltype.alphabet)
    data = alignment.array_seqs[:, [pos1, pos2]]
    valid = data[(data < num_states).all(axis=1), :]
    col1, col2 = valid.T
    counts_1 = CategoryCounter(col1)
    counts_2 = CategoryCounter(col2)
    counts_12 = CategoryCounter([tuple(pr) for pr in valid])
    stat = counts_1.entropy + counts_2.entropy - counts_12.entropy
    if normalised:
        stat /= counts_12.entropy
    return stat


@pytest.fixture(params=["AAAAA", "ACGT", "AACC", "ACGNNGG"])
def column(request):
    return get_moltype("dna", new_type=True).alphabet.to_indices(request.param)


# revised code tests
def test_vector_entropy(column):
    orig = CategoryCounter(column[column < 4])
    counts = c3_coevo._count_states(column, 4)
    got = c3_coevo._vector_entropy(counts)
    assert_allclose(got, orig.entropy)


@pytest.fixture
def alignment():
    data = """
    AAGG-T
    NCGN-C
    NGTG-C
    NTTA-G
    """
    # 0,2,1.3,nan,1.5
    # entropies
    data = {f"s{i}": s.strip() for i, s in enumerate(data.splitlines()) if s.strip()}
    return make_aligned_seqs(data=data, moltype="dna")


def test_column_entropies(alignment):
    expect = alignment.entropy_per_pos()
    got = c3_coevo._calc_column_entropies(alignment.array_seqs, 4)
    assert_allclose(got, expect)


@pytest.mark.parametrize("pos", itertools.combinations(range(5), 2))
def test_joint_entropies(alignment, pos):
    data = alignment.array_seqs[:, pos]
    valid = data[(data < 4).all(axis=1), :]
    valid = [tuple(p) for p in valid]
    expect = CategoryCounter(valid).entropy
    counts = c3_coevo._count_joint_states(data, 4)
    got = c3_coevo._calc_joint_entropy(counts)
    assert_allclose(got, expect)


@pytest.mark.parametrize("normalised", [False, True])
def test_mi_pair_comp(alignment, normalised):
    pair = calc_mi_pair(alignment, pos1=1, pos2=5, normalised=normalised)
    got = c3_coevo.coevolution_matrix(
        alignment=alignment,
        stat="nmi" if normalised else "mi",
        show_progress=False,
    )[5, 1]
    assert_allclose(got, pair)


def test_aa_seq():
    aln = make_aligned_seqs(
        data={
            "FlyingFox": "SQ",
            "DogFaced": "SQ",
            "FreeTaile": "SQ",
            "LittleBro": "SQ",
            "TombBat": "SQ",
            "RoundEare": "SQ",
            "FalseVamp": "SQ",
            "LeafNose": "SQ",
            "Horse": "SQ",
            "Rhino": "SQ",
            "Pangolin": "SL",
            "Cat": "SQ",
            "Dog": "SQ",
            "Llama": "SQ",
            "Pig": "SQ",
            "Cow": "SQ",
            "Hippo": "SQ",
            "SpermWhale": "SQ",
            "HumpbackW": "SQ",
            "Mole": "SQ",
            "Hedgehog": "SQ",
            "TreeShrew": "SQ",
            "FlyingLem": "SQ",
            "Galago": "SQ",
            "HowlerMon": "SQ",
            "Rhesus": "SQ",
            "Orangutan": "SQ",
            "Gorilla": "SQ",
            "Human": "SQ",
            "Chimpanzee": "SQ",
            "Jackrabbit": "SQ",
            "FlyingSqu": "SQ",
            "OldWorld": "SQ",
            "Mouse": "SQ",
            "Rat": "SQ",
            "NineBande": "RQ",
            "HairyArma": "RQ",
            "Anteater": "SQ",
            "Sloth": "SQ",
            "Dugong": "SQ",
            "Manatee": "SQ",
            "AfricanEl": "SQ",
            "AsianElep": "SQ",
            "RockHyrax": "SQ",
            "TreeHyrax": "SQ",
            "Aardvark": "SQ",
            "GoldenMol": "SQ",
            "Madagascar": "SQ",
            "Tenrec": "SQ",
            "LesserEle": "SQ",
            "GiantElep": "SQ",
            "Caenolest": "NQ",
            "Phascogale": "NQ",
            "Wombat": "NQ",
            "Bandicoot": "NQ",
        },
        moltype="protein",
    )
    expect = calc_mi_pair(aln, pos1=0, pos2=1, normalised=True)

    got = c3_coevo.coevolution_matrix(
        alignment=aln,
        stat="nmi",
        show_progress=False,
    )
    assert_allclose(got[1, 0], expect)
    got = c3_coevo.coevolution_matrix(
        alignment=aln,
        stat="rmi",
        show_progress=False,
    )
    expect = 0.06083812598795946  # from original rmi calculation code
    got = got.array[~numpy.isnan(got.array)]
    assert_allclose(got[0], expect)


@pytest.fixture(params=["dna", "protein"])
def pos_pair(request):
    data = [
        "TC",
        "TC",
        "TA",
        "TC",
        "AC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "AC",
        "TC",
        "AC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TA",
        "TC",
        "TC",
        "TC",
        "TC",
        "TA",
        "TC",
        "TC",
        "TC",
        "TC",
        "AC",
        "AC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
        "TC",
    ]
    return make_aligned_seqs(
        data={f"s{i}": p for i, p in enumerate(data)},
        moltype=request.param,
    )


def test_scaled_mi(pos_pair):
    # should work with different molecular types
    # and give the same answer (for this rather simple case)
    got = c3_coevo.coevolution_matrix(
        alignment=pos_pair,
        stat="rmi",
        show_progress=False,
    )[1, 0]
    assert_allclose(got, 8 / 42)


@pytest.mark.parametrize("met1, met2", [("mi", "nmi"), ("mi", "rmi"), ("nmi", "rmi")])
def test_diff_mi_metrics(pos_pair, met1, met2):
    met1 = c3_coevo.coevolution_matrix(
        alignment=pos_pair,
        stat=met1,
        show_progress=False,
    )[1, 0]
    met2 = c3_coevo.coevolution_matrix(
        alignment=pos_pair,
        stat=met2,
        show_progress=False,
    )[1, 0]
    assert not numpy.allclose(met1, met2)


@pytest.fixture
def small_dna():
    return make_aligned_seqs(
        data={f"s{i}": s for i, s in enumerate(["AA", "AA", "GG", "GG", "GC"])},
        moltype="dna",
    )


def test_resampled_mi_interface(small_dna):
    """resampled_mi_alignment should correctly compute statistic from
    alignment"""
    dmat = c3_coevo.coevolution_matrix(
        alignment=small_dna,
        stat="rmi",
        show_progress=False,
    )
    # expected value from hand calculation
    assert_allclose(dmat[1, 0], 0.78333333)


@pytest.fixture
def cols():
    dna = get_moltype("dna", new_type=True)
    alpha = dna.alphabet
    counts1 = numpy.zeros(len(alpha), dtype=numpy.int64)
    counts2 = numpy.zeros(len(alpha), dtype=numpy.int64)
    for a, b in ["AA", "AA", "GG", "GG", "GC"]:
        counts1[alpha.to_indices(a)] += 1
        counts2[alpha.to_indices(b)] += 1
    return counts1, counts2


def _make_expected_array(alpha, data):
    expected = numpy.zeros((4, 4), dtype=numpy.float64)
    for f, weights in data:
        i = alpha.to_indices(f)
        for k, v in weights.items():
            j = alpha.to_indices(k)
            expected[i, j] = v
    return expected


def test_calc_weights_new(cols):
    """resampled mi weights should be correctly computed"""
    dna = get_moltype("dna", new_type=True)
    alpha = dna.alphabet

    e1 = [
        ("A", {"G": 0.1}),
        ("G", {"A": 0.1}),
    ]
    expect_w1 = _make_expected_array(alpha, e1)
    e2 = [
        ("A", {"C": 0.033333333333333333, "G": 0.066666666666666666}),
        ("G", {"A": 0.066666666666666666, "C": 0.033333333333333333}),
        ("C", {"A": 0.050000000000000003, "G": 0.050000000000000003}),
    ]
    expect_w2 = _make_expected_array(alpha, e2)

    c1, c2 = cols
    w1 = c3_coevo._make_weights(c1)
    w2 = c3_coevo._make_weights(c2)
    assert_allclose(w1, expect_w1)
    assert_allclose(w2, expect_w2)


@pytest.fixture
def protein_aln4():
    return make_aligned_seqs(
        dict([("A1", "AACF"), ("A12", "AADF"), ("A123", "ADCF"), ("A111", "AAD-")]),
        moltype="protein",
    )


@pytest.mark.parametrize("stat", ["mi", "nmi"])
def test_alignment_analyses_moltype_protein(protein_aln4, stat):
    """alignment methods work with moltype = PROTEIN"""
    r = c3_coevo.coevolution_matrix(
        alignment=protein_aln4,
        show_progress=False,
        stat=stat,
    )
    assert r.shape == (4, 4)


@pytest.fixture(params=[("AA", "AA"), ("AC",), ("NA", "AC")])
def zero_case(request):
    return make_aligned_seqs(
        {f"s{i}": s for i, s in enumerate(request.param)},
        moltype="dna",
    )


@pytest.mark.parametrize("stat", ["mi", "nmi", "rmi"])
def test_mi_pair_0(zero_case, stat):
    """all MI 0.0"""
    got = c3_coevo.coevolution_matrix(
        alignment=zero_case,
        stat=stat,
        show_progress=False,
    )[1, 0]
    assert_allclose(got, 0.0)


@pytest.fixture(params=[("CG", "AC"), ("CG", "AC", "NN")])
def one_case(request):
    return make_aligned_seqs(
        {f"s{i}": s for i, s in enumerate(request.param)},
        moltype="dna",
    )


@pytest.mark.parametrize("stat", ["mi", "nmi", "rmi"])
def test_mi_pair_1(one_case, stat):
    """all MI 1"""
    got = c3_coevo.coevolution_matrix(
        alignment=one_case,
        stat=stat,
        show_progress=False,
    )[1, 0]
    assert_allclose(got, 1.0)
