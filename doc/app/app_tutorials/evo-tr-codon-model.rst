.. jupyter-execute::
    :hide-code:

    import set_working_directory

Applying a time-reversible codon model
--------------------------------------

We display the full set of codon models available.

.. jupyter-execute::

    from cogent3 import available_models

    available_models("codon")

Using the conditional nucleotide form codon model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CNFGTR model (`Yap et al <https://www.ncbi.nlm.nih.gov/pubmed/19815689>`__) is the most robust of the time-reversible codon models available (`Kaehler et al <https://www.ncbi.nlm.nih.gov/pubmed/28175284>`__). By default, this model does not optimise the codon frequencies but uses the average estimated from the alignment. We configure the model to optimise the root motif probabilities.

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")
    model = get_app("model",
        "CNFGTR",
        tree="data/primate_brca1.tree",
        optimise_motif_probs=True,
    )
    result = model(aln)
    result

.. jupyter-execute::

    result.lf
