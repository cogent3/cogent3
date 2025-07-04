.. jupyter-execute::
    :hide-code:

    import set_working_directory

Applying a time-reversible nucleotide model
-------------------------------------------

We display the available set of nucleotide substitution models.

.. jupyter-execute::

    from cogent3 import available_models

    available_models("nucleotide")

Using the GTR model
^^^^^^^^^^^^^^^^^^^

We specify the general time-reversible model (`Lanave et al <https://www.ncbi.nlm.nih.gov/pubmed/6429346>`__) by its abbreviation. By default, this model does not optimise the codon frequencies but uses the average estimated from the alignment. We configure the model to optimise the root motif probabilities.

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")
    model = get_app("model",
        "GTR", tree="data/primate_brca1.tree", optimise_motif_probs=True
    )
    result = model(aln)
    result

.. jupyter-execute::

    result.lf
