.. jupyter-execute::
    :hide-code:

    import set_working_directory

``natsel_neutral`` – a test for selective neutrality
----------------------------------------------------

We employ codon models to test hypotheses regarding the mode of natural selection that has operated on a gene.

Noting that ω (omega) is the ratio of nonsynonymous substitutions to synonymous substitutions, ω=1 is indicative a gene is evolving neutrally.

.. warning:: I’m setting ``optimise_motif_probs=False`` to speed up execution of the examples, not because it’s a good idea!

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

    omega_eq_1 = get_app("natsel_neutral",
        "GNC", tree="data/primate_brca1.tree", optimise_motif_probs=False
    )
    result = omega_eq_1(aln)
    type(result)

.. jupyter-execute::

    result

.. jupyter-execute::

    result.alt.lf
