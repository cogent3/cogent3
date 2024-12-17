import os
import warnings

import pytest
from numpy.testing import assert_allclose, assert_almost_equal

import cogent3
import cogent3.evolve.parameter_controller
import cogent3.evolve.substitution_model

_NEW_TYPE = "COGENT3_NEW_TYPE" in os.environ
good_rule_sets = [
    [{"par_name": "length", "is_independent": True}],
    [{"par_name": "length", "is_independent": True}],
    [
        {
            "par_name": "length",
            "clade": True,
            "is_independent": True,
            "edges": ["a", "b"],
        },
    ],
    [{"par_name": "length", "is_independent": True, "edges": ["a", "c", "e"]}],
    [{"par_name": "length", "is_independent": True, "edge": "a"}],
]
bad_rule_sets = [[{"par_name": "length", "clade": True, "edges": ["b", "f"]}]]


@pytest.fixture
def setup_data():
    """Fixture to set up data for tests"""
    al = cogent3.make_aligned_seqs(
        data={"a": "tata", "b": "tgtc", "c": "gcga", "d": "gaac", "e": "gagc"},
        moltype="dna",
    )
    tree = cogent3.make_tree(treestring="((a,b),(c,d),e);")
    model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
        equal_motif_probs=True,
        model_gaps=True,
    )
    return al, tree, model


def test_scoped_local(setup_data):
    al, tree, _ = setup_data
    model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
        equal_motif_probs=True,
        model_gaps=True,
        predicates={"kappa": "transition"},
    )
    lf = model.make_likelihood_function(tree)
    lf.set_constant_lengths()
    lf.set_alignment(al)
    null = lf.get_num_free_params()
    lf.set_param_rule(par_name="kappa", is_independent=True, edges=["b", "d"])
    assert null + 2 == lf.get_num_free_params()


def test_set_get_motif_probs_nstat():
    from cogent3 import get_model

    aln = cogent3.make_aligned_seqs(
        data={
            "a": "AACGAAGCAGAGTCACGGCA",
            "b": "ACGGAAGTTGAGTCACCCCA",
            "c": "TGCATCGAAAAGTCACGCTG",
        },
        moltype="dna",
    )
    bases = "ACGT"
    expect = aln.get_motif_probs()
    expect = [expect[b] for b in bases]
    tree = cogent3.make_tree("(a,b,c)")
    gn = get_model("GN")
    lf = gn.make_likelihood_function(tree)
    lf.set_alignment(aln)
    got = lf.get_motif_probs().to_dict()
    got = [got[b] for b in bases]
    assert_allclose(got, expect)


def test_set_motif_probs(setup_data):
    """Mprobs supplied to the parameter controller"""
    al, tree, _ = setup_data

    def compare_mprobs(got, exp):
        # handle min val
        motifs = list(got)
        assert_almost_equal(
            [got[m] for m in motifs],
            [exp[m] for m in motifs],
            decimal=5,
        )

    model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
        model_gaps=True,
        motif_probs=None,
    )
    lf = model.make_likelihood_function(tree, motif_probs_from_align=False)

    mprobs = {"A": 0.1, "C": 0.2, "G": 0.2, "T": 0.5, "-": 0.0}
    lf.set_motif_probs(mprobs)
    # node the LF adjust motif probs so they are all >= 1e-6
    got = lf.get_motif_probs().to_dict()
    compare_mprobs(got, mprobs)

    lf.set_motif_probs_from_data(al[:1], is_constant=True)
    assert_almost_equal(lf.get_motif_probs()["G"], 0.6, decimal=4)

    lf.set_motif_probs_from_data(al[:1], pseudocount=1)
    assert lf.get_motif_probs()["G"] != 0.6

    # test with consideration of ambiguous states
    al = cogent3.make_aligned_seqs(
        data={"seq1": "ACGTAAGNA", "seq2": "ACGTANGTC", "seq3": "ACGTACGTG"},
        moltype="dna",
    )
    lf.set_motif_probs_from_data(al, include_ambiguity=True, is_constant=True)
    motif_probs = dict(lf.get_motif_probs())
    correct_probs = {
        "A": 8.5 / 27,
        "C": 5.5 / 27,
        "-": 0.0,
        "T": 5.5 / 27,
        "G": 7.5 / 27,
    }
    compare_mprobs(motif_probs, correct_probs)
    assert_allclose(sum(motif_probs.values()), 1.0)


def test_set_multilocus(setup_data):
    """2 loci each with own mprobs"""
    _, tree, _ = setup_data
    model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
        motif_probs=None,
    )
    lf = model.make_likelihood_function(
        tree,
        motif_probs_from_align=False,
        loci=["a", "b"],
    )

    mprobs_a = {"A": 0.2, "T": 0.2, "C": 0.3, "G": 0.3}
    mprobs_b = {"A": 0.1, "T": 0.2, "C": 0.3, "G": 0.4}

    for is_constant in [False, True]:
        lf.set_motif_probs(mprobs_a, is_constant=is_constant)
        lf.set_motif_probs(mprobs_b, locus="b")
        assert lf.get_motif_probs(locus="a") == mprobs_a
        assert lf.get_motif_probs(locus="b") == mprobs_b


def test_set_param_rules(setup_data):
    _, tree, model = setup_data
    lf = model.make_likelihood_function(tree)

    def do_rules(rule_set):
        for rule in rule_set:
            lf.set_param_rule(**rule)

    for rule_set in good_rule_sets:
        lf.set_default_param_rules()
        do_rules(rule_set)
    for rule_set in bad_rule_sets:
        lf.set_default_param_rules()
        with pytest.raises((KeyError, TypeError, AssertionError, ValueError)):
            do_rules(rule_set)


def test_set_constant_lengths(setup_data):
    al, _, model = setup_data
    t = cogent3.make_tree(treestring="((a:1,b:2):3,(c:4,d:5):6,e:7);")
    lf = model.make_likelihood_function(t)
    lf.set_param_rule("length", is_constant=True)
    lf.set_alignment(al)
    assert lf.get_param_value("length", "b") == 2
    assert lf.get_param_value("length", "d") == 5


@pytest.mark.skipif(
    _NEW_TYPE,
    reason="test env artifact for new_type, it passes if run alone",
)
def test_pairwise_clock(DATA_DIR):
    al = cogent3.load_aligned_seqs(DATA_DIR / "brca1.fasta", moltype="dna")
    a, b = "Human", "Mouse"
    al = al.take_seqs([a, b]).omit_gap_pos()
    tree = cogent3.make_tree(tip_names=[a, b])
    model = cogent3.evolve.substitution_model.TimeReversibleDinucleotide(
        equal_motif_probs=True,
        model_gaps=True,
        mprob_model="tuple",
    )
    lf = model.make_likelihood_function(tree)
    lf.set_local_clock(a, b)
    lf.set_alignment(al)
    lf.optimise(local=True, show_progress=False)
    rd = lf.get_param_value_dict(["edge"], params=["length"])
    assert rd["length"][a] == rd["length"][b]
    # the lnL is from an old_type fit
    # the following fails with new_type in the test suite (returns a
    # lnL -7435.289193790548), but if run alone it in a jupyter
    # notebook it passes
    assert_almost_equal(lf.lnL, -7376.697602025973)


def test_local_clock(setup_data):
    al, tree, model = setup_data
    lf = model.make_likelihood_function(tree)
    lf.set_local_clock("c", "d")
    lf.set_alignment(al)
    lf.optimise(local=True, tolerance=1e-8, max_restarts=2, show_progress=False)
    rd = lf.get_param_value_dict(["edge"], params=["length"])
    assert_almost_equal(lf.get_log_likelihood(), -27.84254174)
    assert rd["length"]["c"] == rd["length"]["d"]
    assert rd["length"]["a"] != rd["length"]["e"]


def test_complex_parameter_rules(setup_data):
    al, tree, _ = setup_data
    model = cogent3.evolve.substitution_model.TimeReversibleNucleotide(
        equal_motif_probs=True,
        model_gaps=True,
        predicates={"kappa": "transition"},
    )
    lf = model.make_likelihood_function(tree)
    lf.set_param_rule(par_name="kappa", is_independent=True)
    lf.set_param_rule(par_name="kappa", is_independent=False, edges=["b", "d"])
    lf.set_constant_lengths(
        cogent3.make_tree(treestring="((a:1,b:1):1,(c:2,d:1):1,e:1);"),
    )
    lf.set_alignment(al)
    lf.optimise(local=True, show_progress=False)
    rd = lf.get_param_value_dict(["edge"], params=["kappa"])
    assert_almost_equal(lf.get_log_likelihood(), -27.3252, 3)
    assert rd["kappa"]["b"] == rd["kappa"]["d"]
    assert rd["kappa"]["a"] != rd["kappa"]["b"]


def test_bounds(setup_data):
    _, tree, model = setup_data
    lf = model.make_likelihood_function(tree)
    lf.set_param_rule("length", value=3, lower=0, upper=5)

    # Out of bounds value should warn and keep bounded
    with warnings.catch_warnings(record=True) as w:
        lf.set_param_rule("length", lower=0, upper=2, warn=True)
        assert len(w), "No warning issued"
    assert lf.get_param_value("length", edge="a") == 2

    # upper < lower bounds should fail
    with pytest.raises(ValueError):
        lf.set_param_rule("length", lower=2, upper=0)
