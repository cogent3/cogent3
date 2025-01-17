import pytest

from cogent3 import load_aligned_seqs, load_tree
from cogent3.evolve import bootstrap, substitution_model

seqnames = ["Chimpanzee", "Rhesus", "Orangutan", "Human"]

REPLICATES = 2


def float_ge_zero(num, epsilon=1e-6):
    """compare whether a floating point value is >= zero with epsilon
    tolerance."""
    return num >= 0.0 or abs(num - 0.0) < epsilon


@pytest.fixture
def alignobj(DATA_DIR):
    aln = load_aligned_seqs(
        DATA_DIR / "brca1.fasta",
        moltype="dna",
    )
    aln = aln.take_seqs(seqnames)
    return aln[:1000]


@pytest.fixture
def treeobj(DATA_DIR):
    treeobj = load_tree(filename=DATA_DIR / "murphy.tree")
    return treeobj.get_sub_tree(seqnames)


@pytest.fixture
def f81mod():
    return substitution_model.TimeReversibleNucleotide(model_gaps=True)


def create_null_controller(align, tree, submod):
    """A null model controller creator.
    We constrain the human chimp branches to be equal."""
    lf = submod.make_likelihood_function(tree)
    lf.set_alignment(align)
    # we are setting a local molecular clock for human/chimp
    lf.set_local_clock("Human", "Chimpanzee")
    return lf


def create_alt_controller(align, tree, submod):
    """An alternative model controller. Chimp/Human
    branches are free to vary."""
    lf = submod.make_likelihood_function(tree)
    lf.set_alignment(align)
    return lf


def test_conf_int(alignobj, treeobj, f81mod):
    """testing estimation of confidence intervals."""
    lf = create_null_controller(alignobj, treeobj, f81mod)

    bstrap = bootstrap.EstimateConfidenceIntervals(
        lf,
        lambda x: x.get_param_value("length", "Human"),
        alignobj,
    )
    bstrap.set_num_replicates(REPLICATES)
    bstrap.set_seed(1984)
    bstrap.run(local=True, show_progress=False)
    samplelnL = bstrap.get_sample_lnL()
    assert all(lnL < 0.0 for lnL in samplelnL)

    observed_stat = bstrap.get_observed_stats()
    assert float_ge_zero(observed_stat)

    samplestats = bstrap.getSampleStats()

    assert all(float_ge_zero(stat) for stat in samplestats)

    assert len(samplelnL) == REPLICATES
    assert len(samplestats) == REPLICATES


def test_prob(alignobj, treeobj, f81mod):
    """testing estimation of probability."""

    prob_bstrap = bootstrap.EstimateProbability(
        create_null_controller(alignobj, treeobj, f81mod),
        create_alt_controller(alignobj, treeobj, f81mod),
        alignobj,
    )
    prob_bstrap.set_num_replicates(REPLICATES)
    prob_bstrap.set_seed(1984)
    prob_bstrap.run(local=True, show_progress=False)

    assert len(prob_bstrap.get_sample_LR_list()) == REPLICATES

    assert float_ge_zero(prob_bstrap.get_observed_LR())

    # check the returned sample LR's for being > 0.0
    assert all(
        float_ge_zero(sample_LR) for sample_LR in prob_bstrap.get_sample_LR_list()
    )

    # check the returned observed lnL fulfill this assertion too, really
    # testing their order
    null, alt = prob_bstrap.get_observed_lnL()
    assert float_ge_zero(2 * (alt - null))

    # now check the structure of the returned sample
    assert all(
        float_ge_zero(2 * (alt - null)) for null, alt in prob_bstrap.get_sample_lnL()
    )

    # be sure we get something back from getprob if proc rank is 0
    assert float_ge_zero(prob_bstrap.get_estimated_prob())
