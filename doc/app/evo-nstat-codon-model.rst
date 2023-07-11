.. jupyter-execute::
    :hide-code:

    import set_working_directory

Applying GNC, a non-stationary codon model
------------------------------------------

See `Kaehler et al <https://www.ncbi.nlm.nih.gov/pubmed/28175284>`__ for the formal description of this model. Note that perform hypothesis testing using this model elsewhere.

We apply this to a sample alignment.

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

The model is specified using it’s abbreviation.

.. jupyter-execute::

    model = get_app("model", "GNC", tree="data/primate_brca1.tree")
    result = model(aln)
    result

.. jupyter-execute::

    result.lf

We can obtain the tree with branch lengths as ENS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If this tree is written to newick (using the ``write()`` method), the lengths will now be ENS.

.. jupyter-execute::

    tree = result.tree
    fig = tree.get_figure()
    fig.scale_bar = "top right"
    fig.show(width=500, height=500)
