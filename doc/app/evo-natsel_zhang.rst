.. jupyter-execute::
    :hide-code:

    import set_working_directory

``natsel_zhang`` – a branch-site test
=====================================

This is the hypothesis test presented in `Zhang et al <https://www.ncbi.nlm.nih.gov/pubmed/16107592>`__. This test evaluates the hypothesis that a set of sites have undergone positive natural selection on a pre-specified set of lineages.

For this model class, there are groups of branches for which all positions are evolving neutrally but some proportion of those neutrally evolving sites change to adaptively evolving on so-called foreground edges. For the current example, we’ll define the Chimpanzee and Human branches as foreground and everything else as background. The following table defines the parameter scopes.

========== ========== ================ ================
Site Class Proportion Background Edges Foreground Edges
========== ========== ================ ================
0          p0         0 < omega0 < 1   0 < omega0 < 1
1          p1         omega1 = 1       omega1 = 1
2a         p2         0 < omega0 < 1   0 < omega2 > 1
2b         p3         omega1 = 1       0 < omega0 < 1
========== ========== ================ ================

.. note:: Our implementation is not as parametrically succinct as that of Zhang et al, we have 1 additional bin probability.

.. jupyter-execute::

    from cogent3.app import io, evo

    loader = io.load_aligned(format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

    zhang_test = evo.natsel_zhang(
        "GNC",
        tree="data/primate_brca1.tree",
        optimise_motif_probs=False,
        tip1="Human",
        tip2="Chimpanzee",
    )

    result = zhang_test(aln)
    result

.. jupyter-execute::

    result.alt.lf

Getting the posterior probabilities of site-class membership
------------------------------------------------------------

.. jupyter-execute::

    bprobs = result.alt.lf.get_bin_probs()
    bprobs[:, :20]

Getting all the statistics in tabular form
------------------------------------------

.. jupyter-execute::

    tab = evo.tabulate_stats()
    stats = tab(result.alt)
    stats

.. jupyter-execute::

    stats["edge bin params"][:10]  # truncating the table