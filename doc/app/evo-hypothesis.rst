.. jupyter-execute::
    :hide-code:

    import set_working_directory

Testing a hypothesis – non-stationary or time-reversible
========================================================

We test the hypothesis that the GTR model is sufficient for a data set, compared with the GN (non-stationary general nucleotide model).

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

    tree = "data/primate_brca1.tree"

    null = get_app("model", "GTR", tree=tree, optimise_motif_probs=True)
    alt = get_app("model", "GN", tree=tree, optimise_motif_probs=True)
    hyp = get_app("hypothesis", null, alt)
    result = hyp(aln)
    type(result)

``result`` is a ``hypothesis_result`` object. The ``repr()`` displays the likelihood ratio test statistic, degrees of freedom and associated p-value>

.. jupyter-execute::

    result

In this case, we accept the null given the p-value is > 0.05. We use this object to demonstrate the properties of a ``hypothesis_result``.

``hypothesis_result`` has attributes and keys
---------------------------------------------

Accessing the test statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    result.LR, result.df, result.pvalue

The null hypothesis
~~~~~~~~~~~~~~~~~~~

This model is accessed via the ``null`` attribute.

.. jupyter-execute::

    result.null

.. jupyter-execute::

    result.null.lf

The alternate hypothesis
~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    result.alt.lf

Saving hypothesis results
-------------------------

You are advised to save these results as serialised data since this provides maximum flexibility for downstream analyses.

The following would write the result into a ``sqlitedb``.

.. code-block:: python

    from cogent3 import get_app, open_data_store
    
    output = open_data_store("path/to/myresults.sqlitedb", mode="w")
    writer = get_app("write_db", data_store=output)
    writer(result)
