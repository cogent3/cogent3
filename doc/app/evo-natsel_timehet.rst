.. jupyter-execute::
    :hide-code:

    import set_working_directory

``natsel_timehet`` – a test of branch heterogeneity
===================================================

We employ codon models to test whether the mode of natural selection affecting human and chimpanzee lineages is distinctive. This is done by specifying the edges of interest (`Yang 1998 <https://www.ncbi.nlm.nih.gov/pubmed/9580986>`__). (Note I’m setting ``optimise_motif_probs=False`` to speed up execution of the examples, not because it’s a good idea!)

.. jupyter-execute::

    from cogent3.app import evo, io

    loader = io.load_aligned(format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

    hc_differ = evo.natsel_timehet(
        "GNC",
        tree="data/primate_brca1.tree",
        optimise_motif_probs=False,
        tip1="Human",
        tip2="Chimpanzee",
    )
    result = hc_differ(aln)
    result

.. jupyter-execute::

    result.alt.lf
