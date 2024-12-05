import warnings

import numpy
import pytest
from numpy.testing import assert_allclose, assert_equal

import cogent3
from cogent3.evolve import pairwise_distance_numba as pdist_numba
from cogent3.evolve.distance import EstimateDistances
from cogent3.evolve.fast_distance import (
    DistanceMatrix,
    HammingPair,
    JC69Pair,
    LogDetPair,
    ParalinearPair,
    ProportionIdenticalPair,
    TN93Pair,
    _calculators,
    available_distances,
    fill_diversity_matrix,
    get_distance_calculator,
    get_moltype_index_array,
    seq_to_indices,
)
from cogent3.evolve.models import F81, HKY85, JC69

warnings.filterwarnings("ignore", "Not using MPI as mpi4py not found")


# hides the warning from taking log of -ve determinant
numpy.seterr(invalid="ignore")

DNA = cogent3.get_moltype("dna")
PROTEIN = cogent3.get_moltype("protein")
RNA = cogent3.get_moltype("rna")


@pytest.fixture
def dna_char_indices():
    """Fixture providing DNA character indices for tests"""
    return get_moltype_index_array(DNA)


@pytest.fixture
def rna_char_indices():
    """Fixture providing RNA character indices for tests"""
    return get_moltype_index_array(RNA)


@pytest.fixture
def basic_alignment():
    """Fixture providing a basic alignment for tests"""
    return cogent3.make_aligned_seqs(
        data=[("s1", "ACGTACGTAC"), ("s2", "GTGTACGTAC")],
        moltype=DNA,
    )


@pytest.fixture
def ambig_alignment():
    """Fixture providing an ambiguous alignment for tests"""
    return cogent3.make_aligned_seqs(
        data=[("s1", "RACGTACGTACN"), ("s2", "AGTGTACGTACA")],
        moltype=DNA,
    )


@pytest.fixture
def diff_alignment():
    """Fixture providing a different alignment for tests"""
    return cogent3.make_aligned_seqs(
        data=[("s1", "ACGTACGTTT"), ("s2", "GTGTACGTAC")],
        moltype=DNA,
    )


@pytest.fixture
def alignment():
    return cogent3.make_aligned_seqs(
        data=[("s1", "ACGTACGTAC"), ("s2", "GTGTACGTAC")],
        moltype=DNA,
    )


@pytest.fixture(scope="session")
def brca1_5(DATA_DIR):
    return cogent3.load_aligned_seqs(DATA_DIR / "brca1_5.paml", moltype=DNA)


def test_char_to_index(dna_char_indices, rna_char_indices):
    """should correctly recode a DNA & RNA seqs into indices"""
    seq = "TCAGRNY?-"
    expected = [0, 1, 2, 3, -9, -9, -9, -9, -9]
    indices = seq_to_indices(seq, dna_char_indices)
    assert_equal(indices, expected)
    seq = "UCAGRNY?-"
    indices = seq_to_indices(seq, rna_char_indices)
    assert_equal(indices, expected)


@pytest.fixture(scope="session")
def al():
    return cogent3.make_aligned_seqs(
        data={
            "a": "GTACGTACGATC",
            "b": "GTACGTACGTAC",
            "c": "GTACGTACGTTC",
            "e": "GTACGTACTGGT",
        },
        moltype="dna",
    )


@pytest.fixture(scope="session")
def collection():
    return cogent3.make_unaligned_seqs(
        data={
            "a": "GTACGTACGATC",
            "b": "GTACGTACGTAC",
            "c": "GTACGTACGTTC",
            "e": "GTACGTACTGGT",
        },
        moltype="dna",
    )


def assert_dists_approx_equal(expected, observed, atol=1e-6):
    observed = {frozenset(k): v for (k, v) in list(observed.items())}
    expected = {frozenset(k): v for (k, v) in list(expected.items())}
    for key in expected:
        assert_allclose(expected[key], observed[key], rtol=1e-6, atol=atol)


@pytest.fixture
def min_working_example_dmat():
    return DistanceMatrix(
        {
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("B", "C"): 3,
        },
    )


def test_max_pair_mwe(min_working_example_dmat):
    assert min_working_example_dmat.max_pair() == ("B", "C")


def test_min_pair_mwe(min_working_example_dmat):
    assert min_working_example_dmat.min_pair() == ("A", "B")


def test_max_pair_has_max_val():
    aln = cogent3.load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    dmat = aln.distance_matrix()
    got = dmat[dmat.max_pair()]
    expect = dmat.array.max()
    assert got == expect


def test_min_pair_has_min_val():
    aln = cogent3.load_aligned_seqs("data/primate_brca1.fasta", moltype="dna")
    dmat = aln.distance_matrix()
    got = dmat[dmat.min_pair()]
    numpy.fill_diagonal(dmat.array, numpy.inf)
    expect = dmat.array.min()
    assert got == expect


def test_min_max_pair_single_pair():
    dmat = DistanceMatrix({("A", "B"): 2})
    assert dmat.max_pair() == ("A", "B")
    assert dmat.min_pair() == ("A", "B")


def test_max_pair_tied():
    dmat = DistanceMatrix(
        {
            ("A", "B"): 1,
            ("A", "C"): 1,
            ("A", "D"): 3,
            ("B", "C"): 3,
            ("B", "D"): 2,
            ("C", "D"): 2,
        },
    )

    got = set(dmat.max_pair())
    expect = {frozenset(("B", "C")), frozenset(("A", "D"))}
    assert got in expect


def test_min_pair_tied():
    dmat = DistanceMatrix(
        {
            ("A", "B"): 1,
            ("A", "C"): 1,
            ("A", "D"): 3,
            ("B", "C"): 3,
            ("B", "D"): 2,
            ("C", "D"): 2,
        },
    )

    got = set(dmat.min_pair())
    expect = {frozenset(("A", "B")), frozenset(("A", "C"))}
    assert got in expect


def test_dropping_from_matrix():
    """pairwise distances should have method for dropping invalid data"""
    data = {
        ("ABAYE2984", "Atu3667"): None,
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): None,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): None,
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): None,
        ("BAA10469", "Atu3667"): None,
        ("BAA10469", "Avin_42730"): 1.85,
    }

    darr = DistanceMatrix(data)
    new = darr.drop_invalid()
    assert new is None

    data = {
        ("ABAYE2984", "Atu3667"): 0.25,
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): 0.25,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): 0.25,
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): 0.25,
        ("BAA10469", "Atu3667"): 0.25,
        ("BAA10469", "Avin_42730"): 1.85,
    }
    darr = DistanceMatrix(data)
    new = darr.drop_invalid()
    assert new.shape == (2, 2)


def test_distance_matrix_from_array_names():
    data = numpy.array(
        [
            [0.7, 0.1, 0.2, 0.3],
            [0.1, 0.7, 0.1, 0.3],
            [0.3, 0.2, 0.6, 0.3],
            [0.4, 0.1, 0.1, 0.7],
        ],
    )
    names = list("abcd")
    got = DistanceMatrix.from_array_names(data, names)
    assert got["b", "c"] == 0.1


def test_fill_diversity_matrix_all(dna_char_indices):
    """make correct diversity matrix when all chars valid"""
    s1 = seq_to_indices("ACGTACGTAC", dna_char_indices)
    s2 = seq_to_indices("GTGTACGTAC", dna_char_indices)
    matrix = numpy.zeros((4, 4), float)
    # self-self should just be an identity matrix
    fill_diversity_matrix(matrix, s1, s1)
    assert_equal(matrix.sum(), len(s1))
    assert_equal(
        matrix,
        numpy.array([[2, 0, 0, 0], [0, 3, 0, 0], [0, 0, 3, 0], [0, 0, 0, 2]], float),
    )

    # small diffs
    matrix.fill(0)
    fill_diversity_matrix(matrix, s1, s2)
    assert_equal(
        matrix,
        numpy.array([[2, 0, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 0, 2]], float),
    )


def test_fill_diversity_matrix_some(dna_char_indices):
    """make correct diversity matrix when not all chars valid"""
    s1 = seq_to_indices("RACGTACGTACN", dna_char_indices)
    s2 = seq_to_indices("AGTGTACGTACA", dna_char_indices)
    matrix = numpy.zeros((4, 4), float)
    # self-self should just be an identity matrix
    fill_diversity_matrix(matrix, s1, s1)
    assert_equal(matrix.sum(), 10)  # 2 invalid chars
    assert_equal(
        matrix,
        numpy.array([[2, 0, 0, 0], [0, 3, 0, 0], [0, 0, 3, 0], [0, 0, 0, 2]], float),
    )

    # small diffs
    matrix.fill(0)
    fill_diversity_matrix(matrix, s1, s2)
    assert_equal(
        matrix,
        numpy.array([[2, 0, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 0, 2]], float),
    )


def test_hamming_from_matrix(dna_char_indices):
    """compute hamming from diversity matrix"""
    s1 = seq_to_indices("ACGTACGTAC", dna_char_indices)
    s2 = seq_to_indices("GTGTACGTAC", dna_char_indices)
    matrix = numpy.zeros((4, 4), float)
    fill_diversity_matrix(matrix, s1, s2)
    total = 0
    for i in range(4):
        total += matrix[i, i]
    assert_equal(total / len(s1), 0.8)


def test_hamming_pair(basic_alignment):
    """get distances dict"""
    calc = HammingPair(DNA, alignment=basic_alignment)
    calc.run(show_progress=False)
    dists = calc.get_pairwise_distances()
    assert_equal(dists.to_dict()["s1", "s2"], 2)


def test_prop_pair(basic_alignment):
    """get distances dict"""
    calc = ProportionIdenticalPair(DNA, alignment=basic_alignment)
    calc.run(show_progress=False)
    dists = calc.get_pairwise_distances()
    assert_equal(dists.to_dict()["s1", "s2"], 0.2)


def test_jc69_from_matrix(dna_char_indices):
    """compute JC69 from diversity matrix"""
    s1 = seq_to_indices("ACGTACGTAC", dna_char_indices)
    s2 = seq_to_indices("GTGTACGTAC", dna_char_indices)
    matrix = numpy.zeros((4, 4), float)
    fill_diversity_matrix(matrix, s1, s2)
    total = 0
    for i in range(4):
        total += matrix[i, i]
    p = 1 - (total / len(s1))
    assert_allclose(p, 0.2)


def test_wrong_moltype(basic_alignment):
    """specifying wrong moltype raises ValueError"""
    with pytest.raises(ValueError):
        _ = JC69Pair(PROTEIN, alignment=basic_alignment)


def test_jc69_from_alignment(basic_alignment, ambig_alignment, diff_alignment):
    """compute JC69 dists from an alignment"""
    calc = JC69Pair(DNA, alignment=basic_alignment)
    calc.run(show_progress=False)
    assert_equal(calc.lengths["s1", "s2"], 10)
    assert_equal(calc.proportions["s1", "s2"], 0.2)
    # value from OSX MEGA 5
    assert_allclose(calc.dists["s1", "s2"], 0.2326161962)
    # value**2 from OSX MEGA 5
    assert_allclose(calc.variances["s1", "s2"], 0.029752066125078681)
    # value from OSX MEGA 5
    assert_allclose(calc.stderr["s1", "s2"], 0.1724878724)

    # same answer when using ambiguous alignment
    calc.run(ambig_alignment, show_progress=False)
    assert_allclose(calc.dists["s1", "s2"], 0.2326161962)

    # but different answer if subsequent alignment is different
    calc.run(diff_alignment, show_progress=False)
    assert calc.dists["s1", "s2"] != 0.2326161962


def test_tn93_from_matrix(basic_alignment, ambig_alignment, diff_alignment):
    """compute TN93 distances"""
    calc = TN93Pair(DNA, alignment=basic_alignment)
    calc.run(show_progress=False)
    assert_equal(calc.lengths["s1", "s2"], 10)
    assert_equal(calc.proportions["s1", "s2"], 0.2)
    # value from OSX MEGA 5
    assert_allclose(calc.dists["s1", "s2"], 0.2554128119)
    # value**2 from OSX MEGA 5
    assert_allclose(calc.variances["s1", "s2"], 0.04444444445376601)
    # value from OSX MEGA 5
    assert_allclose(calc.stderr["s1", "s2"], 0.2108185107)

    # same answer when using ambiguous alignment
    calc.run(ambig_alignment, show_progress=False)
    assert_allclose(calc.dists["s1", "s2"], 0.2554128119)

    # but different answer if subsequent alignment is different
    calc.run(diff_alignment, show_progress=False)
    assert calc.dists["s1", "s2"] != 0.2554128119


def test_distance_pair(alignment):
    """get distances dict"""
    calc = TN93Pair(DNA, alignment=alignment)
    calc.run(show_progress=False)
    dists = calc.get_pairwise_distances()
    dists = dists.to_dict()
    dist = 0.2554128119
    expect = {("s1", "s2"): dist, ("s2", "s1"): dist}
    assert list(dists.keys()) == list(expect.keys())
    assert_allclose(list(dists.values()), list(expect.values()))


def test_logdet_pair_dna(brca1_5):
    """logdet should produce distances that match MEGA"""
    aln = brca1_5
    logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
    logdet_calc.run(use_tk_adjustment=True, show_progress=False)
    dists = logdet_calc.get_pairwise_distances().to_dict()
    all_expected = {
        ("Human", "NineBande"): 0.075336929999999996,
        ("NineBande", "DogFaced"): 0.0898575452,
        ("DogFaced", "Human"): 0.1061747919,
        ("HowlerMon", "DogFaced"): 0.0934480008,
        ("Mouse", "HowlerMon"): 0.26422862920000001,
        ("NineBande", "Human"): 0.075336929999999996,
        ("HowlerMon", "NineBande"): 0.062202897899999998,
        ("DogFaced", "NineBande"): 0.0898575452,
        ("DogFaced", "HowlerMon"): 0.0934480008,
        ("Human", "DogFaced"): 0.1061747919,
        ("Mouse", "Human"): 0.26539976700000001,
        ("NineBande", "HowlerMon"): 0.062202897899999998,
        ("HowlerMon", "Human"): 0.036571181899999999,
        ("DogFaced", "Mouse"): 0.2652555144,
        ("HowlerMon", "Mouse"): 0.26422862920000001,
        ("Mouse", "DogFaced"): 0.2652555144,
        ("NineBande", "Mouse"): 0.22754789210000001,
        ("Mouse", "NineBande"): 0.22754789210000001,
        ("Human", "Mouse"): 0.26539976700000001,
        ("Human", "HowlerMon"): 0.036571181899999999,
    }
    for pair in dists:
        got = dists[pair]
        expected = all_expected[pair]
        assert_allclose(got, expected)


def test_slice_dmatrix():
    data = {
        ("ABAYE2984", "Atu3667"): 0.25,
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): 0.25,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): 0.25,
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): 0.25,
        ("BAA10469", "Atu3667"): 0.25,
        ("BAA10469", "Avin_42730"): 1.85,
    }
    darr = DistanceMatrix(data)
    names = darr.template.names[0][:3]
    got = darr[:3, :3]
    assert list(got.template.names[0]) == names


def test_logdet_tk_adjustment(brca1_5):
    """logdet using tamura kumar differs from classic"""
    aln = brca1_5
    logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
    logdet_calc.run(use_tk_adjustment=True, show_progress=False)
    tk = logdet_calc.get_pairwise_distances()
    logdet_calc.run(use_tk_adjustment=False, show_progress=False)
    not_tk = logdet_calc.get_pairwise_distances()
    assert tk != not_tk


def test_logdet_pair_aa(brca1_5):
    """logdet shouldn't fail to produce distances for aa seqs"""
    aln = brca1_5
    aln = aln.get_translation()
    logdet_calc = LogDetPair(moltype=PROTEIN, alignment=aln)
    logdet_calc.run(use_tk_adjustment=True, show_progress=False)
    dists = logdet_calc.get_pairwise_distances()
    assert isinstance(dists, DistanceMatrix)


def test_logdet_missing_states():
    """should calculate logdet measurement with missing states"""
    data = [
        (
            "seq1",
            "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        ),
        (
            "seq2",
            "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTNTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        ),
    ]
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)
    logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
    logdet_calc.run(use_tk_adjustment=True, show_progress=False)

    dists = logdet_calc.get_pairwise_distances().to_dict()
    assert next(iter(dists.values())) is not None

    logdet_calc.run(use_tk_adjustment=False, show_progress=False)
    dists = logdet_calc.get_pairwise_distances().to_dict()
    assert next(iter(dists.values())) is not None


def test_logdet_variance():
    """calculate logdet variance consistent with hand calculation"""
    data = [
        (
            "seq1",
            "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        ),
        (
            "seq2",
            "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        ),
    ]
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)
    logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
    logdet_calc.run(use_tk_adjustment=True, show_progress=False)
    assert logdet_calc.variances[1, 1] is None

    index = dict(list(zip("ACGT", list(range(4)), strict=False)))
    J = numpy.zeros((4, 4))
    for p in zip(data[0][1], data[1][1], strict=False):
        J[index[p[0]], index[p[1]]] += 1
    for i in range(4):
        if J[i, i] == 0:
            J[i, i] += 0.5
    J /= J.sum()
    M = numpy.linalg.inv(J)
    var = 0.0
    for i in range(4):
        for j in range(4):
            var += M[j, i] ** 2 * J[i, j] - 1
    var /= 16 * len(data[0][1])

    logdet_calc.run(use_tk_adjustment=False, show_progress=False)
    logdet_calc.get_pairwise_distances()
    assert_allclose(logdet_calc.variances[1, 1], var, atol=1e-3)


def test_logdet_for_determinant_lte_zero():
    """returns distance of None if the determinant is <= 0"""
    data = {
        "seq1": "AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        "seq2": "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
    }
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)

    logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
    logdet_calc.run(use_tk_adjustment=True, show_progress=False)
    dists = logdet_calc.get_pairwise_distances().to_dict()
    assert numpy.isnan(next(iter(dists.values())))
    logdet_calc.run(use_tk_adjustment=False, show_progress=False)
    dists = logdet_calc.get_pairwise_distances().to_dict()
    assert numpy.isnan(next(iter(dists.values())))

    # but raises ArithmeticError if told to
    logdet_calc = LogDetPair(moltype=DNA, alignment=aln, invalid_raises=True)
    with pytest.raises(ArithmeticError):
        logdet_calc.run(use_tk_adjustment=True, show_progress=False)


def test_paralinear_pair_aa():
    """paralinear shouldn't fail to produce distances for aa seqs"""
    aln = cogent3.load_aligned_seqs("data/brca1_5.paml", moltype=DNA)
    aln = aln.get_translation()
    paralinear_calc = ParalinearPair(moltype=PROTEIN, alignment=aln)
    paralinear_calc.run(show_progress=False)
    dists = paralinear_calc.get_pairwise_distances()
    assert isinstance(dists, DistanceMatrix)


def test_paralinear_distance():
    """calculate paralinear variance consistent with hand calculation"""
    data = [
        (
            "seq1",
            "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        ),
        (
            "seq2",
            "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        ),
    ]
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)
    paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
    paralinear_calc.run(show_progress=False)

    index = dict(list(zip("ACGT", list(range(4)), strict=False)))
    J = numpy.zeros((4, 4))
    for p in zip(data[0][1], data[1][1], strict=False):
        J[index[p[0]], index[p[1]]] += 1
    for i in range(4):
        if J[i, i] == 0:
            J[i, i] += 0.5
    J /= J.sum()
    numpy.linalg.inv(J)
    f = J.sum(1), J.sum(0)
    dist = -0.25 * numpy.log(
        numpy.linalg.det(J) / numpy.sqrt(f[0].prod() * f[1].prod()),
    )

    assert_allclose(paralinear_calc.dists["seq1", "seq2"], dist)


def test_paralinear_variance():
    """calculate paralinear variance consistent with hand calculation"""
    data = [
        (
            "seq1",
            "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        ),
        (
            "seq2",
            "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        ),
    ]
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)
    paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
    paralinear_calc.run(show_progress=False)

    index = dict(list(zip("ACGT", list(range(4)), strict=False)))
    J = numpy.zeros((4, 4))
    for p in zip(data[0][1], data[1][1], strict=False):
        J[index[p[0]], index[p[1]]] += 1
    for i in range(4):
        if J[i, i] == 0:
            J[i, i] += 0.5
    J /= J.sum()
    M = numpy.linalg.inv(J)
    f = J.sum(1), J.sum(0)
    var = 0.0
    for i in range(4):
        for j in range(4):
            var += M[j, i] ** 2 * J[i, j]
        var -= 1 / numpy.sqrt(f[0][i] * f[1][i])
    var /= 16 * len(data[0][1])

    assert_allclose(paralinear_calc.variances[1, 1], var, atol=1e-3)


def test_paralinear_for_determinant_lte_zero():
    """returns distance of None if the determinant is <= 0"""
    data = {
        "seq1": "AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        "seq2": "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
    }
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)

    paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
    paralinear_calc.run(show_progress=False)
    dists = paralinear_calc.get_pairwise_distances().to_dict()
    assert numpy.isnan(next(iter(dists.values())))
    paralinear_calc.run(show_progress=False)
    dists = paralinear_calc.get_pairwise_distances().to_dict()
    assert numpy.isnan(next(iter(dists.values())))


def test_paralinear_pair_dna():
    """calculate paralinear distance consistent with logdet distance"""
    data = [
        (
            "seq1",
            "TAATTCATTGGGACGTCGAATCCGGCAGTCCTGCCGCAAAAGCTTCCGGAATCGAATTTTGGCA",
        ),
        (
            "seq2",
            "AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGG",
        ),
    ]
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)
    paralinear_calc = ParalinearPair(moltype=DNA, alignment=aln)
    paralinear_calc.run(show_progress=False)
    logdet_calc = LogDetPair(moltype=DNA, alignment=aln)
    logdet_calc.run(show_progress=False)
    assert logdet_calc.dists[1, 1] == paralinear_calc.dists[1, 1]
    assert paralinear_calc.variances[1, 1] == logdet_calc.variances[1, 1]


def get_calc(data):
    aln = cogent3.make_aligned_seqs(data=data, moltype=DNA)
    calc = ParalinearPair(moltype=DNA, alignment=aln)
    calc(show_progress=False)
    return calc


def test_duplicated():
    """correctly identifies duplicates"""

    # no duplicates
    data = [
        (
            "seq1",
            "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        ),
        (
            "seq2",
            "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        ),
    ]
    calc = get_calc(data)
    assert calc.duplicated is None
    data = [
        (
            "seq1",
            "GGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        ),
        (
            "seq2",
            "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        ),
        (
            "seq3",
            "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        ),
    ]
    calc = get_calc(data)
    assert calc.duplicated in ({"seq2": ["seq3"]}, {"seq3": ["seq2"]})
    # default to get all pairwise distances
    pwds = calc.get_pairwise_distances().to_dict()
    assert pwds[("seq2", "seq3")] == 0.0
    assert pwds[("seq2", "seq1")] == pwds[("seq3", "seq1")]

    # only unique seqs when using include_duplicates=False

    pwds = calc.get_pairwise_distances(include_duplicates=False).to_dict()
    present = next(iter(calc.duplicated.keys()))
    missing = calc.duplicated[present][0]
    assert {(present, missing)} == {("seq2", "seq3")}
    assert (present, "seq1") in pwds
    assert (missing, "seq1") not in pwds


def test_get_calculator():
    """exercising getting specified calculator"""
    for key in _calculators:
        get_distance_calculator(key)
        get_distance_calculator(key.upper())

    with pytest.raises(ValueError):
        get_distance_calculator("blahblah")


def test_available_distances():
    """available_distances has correct content"""
    content = available_distances()
    assert content.shape == (6, 2)
    assert content["tn93", 1] == "dna, rna"


def test_to_dict():
    """distance matrix correctly produces a 1D dict"""
    data = {("s1", "s2"): 0.25, ("s2", "s1"): 0.25}
    dmat = DistanceMatrix(data)
    got = dmat.to_dict()
    assert got == data


def test_matrix_dtype():
    """tests DistanceMatrix correctly accepts the data with proper dtype"""
    data = {
        ("ABAYE2984", "Atu3667"): None,
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): None,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): None,
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): None,
        ("BAA10469", "Atu3667"): None,
        ("BAA10469", "Avin_42730"): 1.85,
    }
    names = set()
    for p in data:
        names.update(p)

    # tests when data has None values and DistanceMatrix using dtype('float')
    darr = DistanceMatrix(data)
    assert darr.shape == (4, 4)
    assert set(darr.names) == names
    for (a, b), dist in data.items():
        if dist is None:
            assert numpy.isnan(darr[a, b])
        else:
            assert_allclose(dist, darr[a, b])

    data = {
        ("ABAYE2984", "Atu3667"): "None",
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): None,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): "None",
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): None,
        ("BAA10469", "Atu3667"): None,
        ("BAA10469", "Avin_42730"): 1.85,
    }

    # tests when data has str values and DistanceMatrix using dtype('float')
    with pytest.raises(ValueError):
        _ = DistanceMatrix(data)


def test_EstimateDistances(tmp_path, al):
    """testing (well, exercising at least), EstimateDistances"""
    d = EstimateDistances(al, JC69())
    d.run(show_progress=False)
    canned_result = {
        ("b", "e"): 0.440840,
        ("c", "e"): 0.440840,
        ("a", "c"): 0.088337,
        ("a", "b"): 0.188486,
        ("a", "e"): 0.440840,
        ("b", "c"): 0.0883373,
    }
    result = d.get_pairwise_distances().to_dict()
    assert_dists_approx_equal(canned_result, result)

    # excercise writing to file
    d.write(tmp_path / "junk.txt")


def test_take_dists():
    """subsets the distance matrix"""
    data = {
        ("ABAYE2984", "Atu3667"): 0.25,
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): 0.25,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): 0.25,
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): 0.25,
        ("BAA10469", "Atu3667"): 0.25,
        ("BAA10469", "Avin_42730"): 1.85,
    }
    darr = DistanceMatrix(data)
    got1 = darr.take_dists(["ABAYE2984", "Atu3667", "Avin_42730"])
    got2 = darr.take_dists("BAA10469", negate=True)
    assert_allclose(got1.array.astype(float), got2.array.astype(float))


def test_build_phylogeny():
    """build a NJ tree"""
    from cogent3 import make_tree

    dists = {
        ("DogFaced", "FlyingFox"): 0.05,
        ("DogFaced", "FreeTaile"): 0.14,
        ("DogFaced", "LittleBro"): 0.16,
        ("DogFaced", "TombBat"): 0.15,
        ("FlyingFox", "DogFaced"): 0.05,
        ("FlyingFox", "FreeTaile"): 0.12,
        ("FlyingFox", "LittleBro"): 0.13,
        ("FlyingFox", "TombBat"): 0.14,
        ("FreeTaile", "DogFaced"): 0.14,
        ("FreeTaile", "FlyingFox"): 0.12,
        ("FreeTaile", "LittleBro"): 0.09,
        ("FreeTaile", "TombBat"): 0.1,
        ("LittleBro", "DogFaced"): 0.16,
        ("LittleBro", "FlyingFox"): 0.13,
        ("LittleBro", "FreeTaile"): 0.09,
        ("LittleBro", "TombBat"): 0.12,
        ("TombBat", "DogFaced"): 0.15,
        ("TombBat", "FlyingFox"): 0.14,
        ("TombBat", "FreeTaile"): 0.1,
        ("TombBat", "LittleBro"): 0.12,
    }
    dists = DistanceMatrix(dists)
    got = dists.quick_tree(show_progress=False)
    expect = make_tree(
        treestring="((TombBat,(DogFaced,FlyingFox)),LittleBro,FreeTaile)",
    )
    assert expect.same_topology(got)


def test_names():
    """names property works"""
    data = {
        ("ABAYE2984", "Atu3667"): 0.25,
        ("ABAYE2984", "Avin_42730"): 0.638,
        ("ABAYE2984", "BAA10469"): None,
        ("Atu3667", "ABAYE2984"): 0.25,
        ("Atu3667", "Avin_42730"): 2.368,
        ("Atu3667", "BAA10469"): 0.25,
        ("Avin_42730", "ABAYE2984"): 0.638,
        ("Avin_42730", "Atu3667"): 2.368,
        ("Avin_42730", "BAA10469"): 1.85,
        ("BAA10469", "ABAYE2984"): 0.25,
        ("BAA10469", "Atu3667"): 0.25,
        ("BAA10469", "Avin_42730"): 1.85,
    }
    names = set()
    for p in data:
        names.update(p)
    darr = DistanceMatrix(data)
    assert set(darr.names) == names
    darr = darr.drop_invalid()
    for n in ("ABAYE2984", "BAA10469"):
        names.remove(n)
    assert set(darr.names) == names


def test_EstimateDistancesWithMotifProbs(al):
    """EstimateDistances with supplied motif probs"""
    motif_probs = {"A": 0.1, "C": 0.2, "G": 0.2, "T": 0.5}
    d = EstimateDistances(al, HKY85(), motif_probs=motif_probs)
    d.run(show_progress=False)
    canned_result = {
        ("a", "c"): 0.07537,
        ("b", "c"): 0.07537,
        ("a", "e"): 0.39921,
        ("a", "b"): 0.15096,
        ("b", "e"): 0.39921,
        ("c", "e"): 0.37243,
    }
    result = d.get_pairwise_distances().to_dict()
    assert_dists_approx_equal(canned_result, result, atol=1e-5)


def test_EstimateDistances_fromThreeway(al):
    """testing (well, exercising at least), EsimateDistances fromThreeway"""
    d = EstimateDistances(al, JC69(), threeway=True)
    d.run(show_progress=False)
    canned_result = {
        ("b", "e"): 0.495312,
        ("c", "e"): 0.479380,
        ("a", "c"): 0.089934,
        ("a", "b"): 0.190021,
        ("a", "e"): 0.495305,
        ("b", "c"): 0.0899339,
    }
    result = d.get_pairwise_distances(summary_function="mean").to_dict()
    assert_dists_approx_equal(canned_result, result, atol=1e-4)


def test_EstimateDistances_fromUnaligned(collection):
    """Excercising estimate distances from unaligned sequences"""
    d = EstimateDistances(collection, JC69(), do_pair_align=True, rigorous_align=True)
    d.run(show_progress=False)
    canned_result = {
        ("b", "e"): 0.440840,
        ("c", "e"): 0.440840,
        ("a", "c"): 0.088337,
        ("a", "b"): 0.188486,
        ("a", "e"): 0.440840,
        ("b", "c"): 0.0883373,
    }
    result = d.get_pairwise_distances().to_dict()
    assert_dists_approx_equal(canned_result, result)

    d = EstimateDistances(collection, JC69(), do_pair_align=True, rigorous_align=False)
    d.run(show_progress=False)
    canned_result = {
        ("b", "e"): 0.440840,
        ("c", "e"): 0.440840,
        ("a", "c"): 0.088337,
        ("a", "b"): 0.188486,
        ("a", "e"): 0.440840,
        ("b", "c"): 0.0883373,
    }
    result = d.get_pairwise_distances().to_dict()
    assert_dists_approx_equal(canned_result, result)


def test_EstimateDistances_other_model_params(al):
    """test getting other model params from EstimateDistances"""
    d = EstimateDistances(al, HKY85(), est_params=["kappa"])
    d.run(show_progress=False)
    # this will be a Number object with Mean, Median etc ..
    kappa = d.get_param_values("kappa")
    assert_allclose(kappa.mean, 0.8939, atol=1e-4)
    # this will be a dict with pairwise instances, it's called by the above
    # method, so the correctness of it's values is already checked
    _ = d.get_pairwise_param("kappa")


def test_EstimateDistances_modify_lf(al):
    """tests modifying the lf"""

    def constrain_fit(lf):
        lf.set_param_rule("kappa", is_constant=True)
        lf.optimise(local=True, show_progress=False)
        return lf

    d = EstimateDistances(al, HKY85(), modify_lf=constrain_fit)
    d.run(show_progress=False)
    result = d.get_pairwise_distances().to_dict()
    d = EstimateDistances(al, F81())
    d.run(show_progress=False)
    expect = d.get_pairwise_distances().to_dict()
    assert_dists_approx_equal(expect, result)


def test_get_raw_estimates(al):
    """correctly return raw result object"""
    d = EstimateDistances(al, HKY85(), est_params=["kappa"])
    d.run(show_progress=False)
    expect = {
        ("a", "b"): {
            "kappa": 1.0000226766004808e-06,
            "length": 0.18232155856115662,
        },
        ("a", "c"): {
            "kappa": 1.0010380037049357e-06,
            "length": 0.087070406623635604,
        },
        ("a", "e"): {"kappa": 2.3965871843412687, "length": 0.4389176272584539},
        ("b", "e"): {"kappa": 2.3965871854366592, "length": 0.43891762729173389},
        ("b", "c"): {
            "kappa": 1.0010380037049357e-06,
            "length": 0.087070406623635604,
        },
        ("c", "e"): {"kappa": 0.57046787478038707, "length": 0.43260232210282784},
    }
    got = d.get_all_param_values()
    for pair in expect:
        for param in expect[pair]:
            assert_allclose(got[pair][param], expect[pair][param], atol=1e-6)


def test_no_calc():
    """returns None if no calculation done"""
    al = cogent3.load_aligned_seqs("data/brca1_5.paml")
    d = EstimateDistances(al, submodel=HKY85())
    assert d.get_pairwise_distances() is None


def test_to_table():
    """converts a distance matrix to a Table"""
    data = {
        ("A", "B"): 2,
        ("A", "C"): 3,
        ("B", "C"): 1,
        ("B", "A"): 2,
        ("C", "A"): 3,
        ("C", "B"): 1,
    }
    darr = DistanceMatrix(data)
    table = darr.to_table()
    assert table.shape == (3, 4)
    assert table.columns["names"].tolist() == list(darr.names)
    assert table["A", "B"] == 2
    assert table["A", "A"] == 0


@pytest.mark.parametrize(
    "aln",
    ["basic_alignment", "ambig_alignment", "diff_alignment"],
)
def test_jc69_dists(request, aln):
    """full numba implementation matches original"""
    aln = request.getfixturevalue(aln)
    calc = JC69Pair(DNA, alignment=aln)
    calc.run(show_progress=False)
    expect = calc.dists.array
    got = pdist_numba.jc69(aln)
    assert_allclose(expect, got)


@pytest.mark.parametrize("aln", ["basic_alignment", "diff_alignment"])
def test_tn93_dists(request, aln):
    """full numba implementation matches original"""
    # note: I excluded the ambig_alignment case because the
    # old and new approaches differ in how they compute the
    # nucleotide frequencies. The old approach would effectively
    # exclude data in aligned columns that contained a
    # non-canonical state, the approach computes the nucleotide
    # frequencies from the whole alignment and so does not exclude
    # valid states just because another sequence has an ambiguity
    aln = request.getfixturevalue(aln)
    # convert to new type alignment
    aln = cogent3.make_aligned_seqs(data=aln.to_dict(), moltype="dna", new_type=True)
    calc = TN93Pair(aln.moltype, alignment=aln)
    calc.run(show_progress=False)
    expect = calc.dists.array
    got = pdist_numba.tn93(aln, parallel=False)
    assert_allclose(got.array, expect)


@pytest.mark.parametrize("aln", ["basic_alignment", "diff_alignment"])
def test_paralinear_dists(request, aln):
    aln = request.getfixturevalue(aln)
    # convert to new type alignment
    aln = cogent3.make_aligned_seqs(data=aln.to_dict(), moltype="dna", new_type=True)
    calc = ParalinearPair(aln.moltype, alignment=aln)
    calc.run(show_progress=False)
    expect = calc.dists.array
    got = pdist_numba.paralinear(aln, parallel=False)
    assert_allclose(got.array, expect, atol=1e-6)


@pytest.mark.parametrize(
    "aln",
    ["basic_alignment", "ambig_alignment", "diff_alignment"],
)
def test_hamming_dists(request, aln):
    aln = request.getfixturevalue(aln)
    # convert to new type alignment
    aln = cogent3.make_aligned_seqs(data=aln.to_dict(), moltype="dna", new_type=True)
    calc = HammingPair(aln.moltype, alignment=aln)
    calc.run(show_progress=False)
    expect = calc.dists.array
    got = pdist_numba.hamming(aln, parallel=False)
    assert_allclose(got.array, expect)


@pytest.mark.parametrize(
    "aln",
    ["basic_alignment", "ambig_alignment", "diff_alignment"],
)
def test_prop_dists(request, aln):
    aln = request.getfixturevalue(aln)
    # convert to new type alignment
    aln = cogent3.make_aligned_seqs(data=aln.to_dict(), moltype="dna", new_type=True)
    calc = ProportionIdenticalPair(aln.moltype, alignment=aln)
    calc.run(show_progress=False)
    expect = calc.dists.array
    got = pdist_numba.pdist(aln, parallel=False)
    assert_allclose(got.array, expect)


@pytest.mark.parametrize("calc", ["jc69", "tn93", "paralinear", "hamming", "pdist"])
def test_numba_get_dist(calc, basic_alignment) -> None:
    aln = cogent3.make_aligned_seqs(
        data=basic_alignment.to_dict(),
        moltype="dna",
        new_type=True,
    )
    calc = pdist_numba.get_distance_calculator(calc)
    got = calc(aln)
    assert isinstance(got, DistanceMatrix)


@pytest.mark.parametrize("calc", ["tn93", "paralinear"])
def test_invalid_moltype_fast_distances(calc, basic_alignment):
    from cogent3.core import new_moltype

    aln = cogent3.make_aligned_seqs(
        data=basic_alignment.to_dict(),
        moltype="protein",
        new_type=True,
    )
    calc = pdist_numba.get_distance_calculator(calc)
    with pytest.raises(new_moltype.MolTypeError):
        calc(aln, invalid_raises=True)


def test_new_paralinear_for_determinant_lte_zero():
    """returns distance of None if the determinant is <= 0"""
    data = {
        "seq1": "AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        "seq2": "TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
    }
    aln = cogent3.make_aligned_seqs(data=data, moltype="dna", new_type=True)

    dists = pdist_numba.paralinear(aln, parallel=False)
    assert numpy.isnan(dists.array).any()

    # raises ArithmeticError if told to
    with pytest.raises(ArithmeticError):
        dists = pdist_numba.paralinear(aln, parallel=False, invalid_raises=True)

    # also raised if alignment method used
    with pytest.raises(ArithmeticError):
        aln.distance_matrix(calc="paralinear", drop_invalid=False)


def test_unknown_dist():
    with pytest.raises(ValueError):
        pdist_numba.get_distance_calculator("blah")


@pytest.mark.parametrize("calc", ["jc69", "tn93", "paralinear", "hamming", "pdist"])
def test_compare_parallel_serial(DATA_DIR, calc):
    aln = cogent3.load_aligned_seqs(
        DATA_DIR / "brca1.fasta",
        moltype="dna",
        new_type=True,
    )
    serial = aln.distance_matrix(calc=calc, parallel=False)
    parallel = aln.distance_matrix(calc=calc, parallel=True)
    assert_allclose(serial.array, parallel.array)
