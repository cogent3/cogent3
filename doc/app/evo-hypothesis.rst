.. jupyter-execute::
    :hide-code:

    import set_working_directory

Testing a hypothesis – non-stationary or time-reversible
========================================================

We evaluate whether the GTR model is sufficient for a data set, compared with the GN (non-stationary general nucleotide model).

.. jupyter-execute::

    from cogent3.app import io, evo, sample

    loader = io.load_aligned(format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")

.. jupyter-execute::

    tree = "data/primate_brca1.tree"
    sm_args = dict(optimise_motif_probs=True)

    null = evo.model("GTR", tree=tree, sm_args=sm_args)
    alt = evo.model("GN", tree=tree, sm_args=sm_args)
    hyp = evo.hypothesis(null, alt)
    result = hyp(aln)
    type(result)

``result`` is a ``hypothesis_result`` object. The ``repr()`` displays the likelihood ratio test statistic, degrees of freedom and associated p-value>

.. jupyter-execute::

    result

In this case, we accept the null given the p-value is > 0.05. We still use this object to demonstrate the properties of a ``hypothesis_result``.

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

You are advised to save these results as json using the standard json writer, or the db writer.

This following would write the result into a ``tinydb``.

.. code-block:: python

    from cogent3.app.io import write_db

    writer = write_db("path/to/myresults.tinydb", create=True, if_exists="overwrite")
    writer(result)