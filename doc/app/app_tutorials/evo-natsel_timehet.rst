.. jupyter-execute::
    :hide-code:

    import set_working_directory

``natsel_timehet`` – a test of branch heterogeneity
---------------------------------------------------

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

We employ codon models to test whether the mode of natural selection affecting human and chimpanzee lineages is distinctive. This is done by specifying the edges of interest (`Yang 1998 <https://www.ncbi.nlm.nih.gov/pubmed/9580986>`__).

.. warning:: I’m setting ``optimise_motif_probs=False`` to speed up execution of the examples, not because it’s a good idea!

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

    hc_differ = get_app("natsel_timehet",
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
