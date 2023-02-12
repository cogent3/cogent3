.. jupyter-execute::
    :hide-code:

    import set_working_directory

Apply a non-stationary nucleotide model to an alignment with 3 sequences
========================================================================

We load some sample data first and select just 3 sequences.

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    select_seqs = get_app("take_named_seqs", "Human", "Rhesus", "Galago")
    process = loader + select_seqs
    aln = process("data/primate_brca1.fasta")
    aln.names

We analyses these using the general Markov nucleotide, GN, model.
Because we analyse just 3 sequences, there is only one possible unrooted
tree, hence it is not required to specify the tree in this instance.

.. jupyter-execute::

    gn = get_app("model", "GN")
    gn

We apply this to ``aln``.

.. jupyter-execute::

    fitted = gn(aln)
    type(fitted)

``model_result``
----------------

As the output above indicates, ``fitted`` is a ``model_result`` object.

This object provides an interface for accessing attributes of a fitted model. The representation display (below), a styled table in a jupyter notebook, presents a summary view with the log-likelihood (``lnL``), number of free parameters (``nfp``) and whether all matrices satisfied the identifiability conditions diagonal largest in column (DLC) and a unique mapping of Q to P. (For description of these quantities and why they matter see `Chang 1996 <https://www.ncbi.nlm.nih.gov/pubmed/?term=8854662>`__ and `Kaehler et al <https://www.ncbi.nlm.nih.gov/pubmed/25503772>`__.)

``model_result`` has dictionary behaviour, hence the ``key`` column. This will be demonstrated below.

.. jupyter-execute::

    fitted

More detail on the fitted model are available via attributes. For instance, display the maximum likelihood estimates via the likelihood function attribute

.. jupyter-execute::

    fitted.lf

.. jupyter-execute::

    fitted.lnL, fitted.nfp

.. jupyter-execute::

    fitted.source

The ``model_result.tree`` attribute is an "annotated tree". Maximum likelihood estimates from the model have been assigned to the tree. Of particular significance, the "length" attribute corresponds to the expected number of substitutions (or ENS). For a non-stationary model, like GN, this can be different to the conventional length (`Kaehler et al <https://www.ncbi.nlm.nih.gov/pubmed/25503772>`__).

.. jupyter-execute::

    fitted.tree, fitted.alignment

We can access the sum of all branch lengths. Either as "ENS" or "paralinear" using the ``total_length()`` method.

.. jupyter-execute::

    fitted.total_length(length_as="paralinear")

Fitting a separate nucleotide model to each codon position
----------------------------------------------------------

Controlled by setting ``split_codons=True``.

.. jupyter-execute::

    gn = get_app("model", "GN", split_codons=True)

    fitted = gn(aln)
    fitted

The model fit statistics, ``lnL`` and ``nfp`` are now sums of the corresponding values from the fits to the individual positions. The ``DLC`` and ``unique_Q`` are also a summary across all models. These only achieve the value ``True`` when all matrices, from all models, satisfy the condition.

We get access to the likelihood functions of the individual positions via the indicated dict keys.

.. jupyter-execute::

    fitted[3]
