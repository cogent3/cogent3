.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _rate-heterogeneity:

Analysis of rate heterogeneity
==============================

.. sectionauthor:: Gavin Huttley

A simple example for analyses involving rate heterogeneity among sites. In this case we will simulate an alignment with two rate categories and then try to recover the rates from the alignment.

.. jupyter-execute::

    from cogent3 import load_tree
    from cogent3.evolve.substitution_model import TimeReversibleNucleotide

Make an alignment with equal split between rates 0.6 and 0.2, and then concatenate them to create a new alignment.

.. jupyter-execute::

    model = TimeReversibleNucleotide(equal_motif_probs=True)
    tree = load_tree("data/test.tree")
    lf = model.make_likelihood_function(tree)
    lf.set_param_rule("length", value=0.6, is_constant=True)
    aln1 = lf.simulate_alignment(sequence_length=10000)
    lf.set_param_rule("length", value=0.2, is_constant=True)
    aln2 = lf.simulate_alignment(sequence_length=10000)
    aln3 = aln1 + aln2

Start from scratch, optimising only rates and the rate probability ratio.

.. jupyter-execute::

    model = TimeReversibleNucleotide(
        equal_motif_probs=True, ordered_param="rate", distribution="free"
    )
    lf = model.make_likelihood_function(tree, bins=2, digits=2, space=3)
    lf.set_alignment(aln3)
    lf.optimise(local=True, max_restarts=2, show_progress=False)

We want to know the bin probabilities and the posterior probabilities.

.. jupyter-execute::

    bprobs = [t for t in lf.get_statistics() if "bin" in t.title][0]
    bprobs

We'll now use a gamma distribution on the sample alignment, specifying the number of bins as 4. We specify that the bins have equal density using the ``lf.set_param_rule('bprobs', is_constant=True)`` command.

.. jupyter-execute::

    model = TimeReversibleNucleotide(
        equal_motif_probs=True, ordered_param="rate", distribution="gamma"
    )
    lf = model.make_likelihood_function(tree, bins=4)
    lf.set_param_rule("bprobs", is_constant=True)
    lf.set_alignment(aln3)
    lf.optimise(local=True, max_restarts=2, show_progress=False)
