.. jupyter-execute::
    :hide-code:

    import set_working_directory

Extracting maximum likelihood estimates from a ``model_result``
===============================================================

If you want to get the stats from a fitted model, use the ``tabulate_stats`` app.

We demonstrate this by first fitting a model.

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")
    model = get_app("model", "GN", tree="data/primate_brca1.tree")
    result = model(aln)

Create and apply ``tabulate_stats`` app
---------------------------------------

.. jupyter-execute::

    tabulator = get_app("tabulate_stats")
    tabulated = tabulator(result)
    tabulated

``tabulated`` is a ``tabular_result`` instance which, like other result types, has ``dict`` like behaviour. It also contains key/value pairs for each model parameter type.

Edge parameters
---------------

These are all parameters that differ between edges. Since the current model is time-homogeneous (a single rate matrix), the table only has entries for the branch scalar (denoted “length”).

.. jupyter-execute::

    tabulated["edge params"]

.. note:: Unless the model is time-reversible, the lengths in that table are not ENS (`Kaehler et al <https://www.ncbi.nlm.nih.gov/pubmed/28175284>`__). As we used a non-stationary nucleotide model in this example, the length values are a scalar used to adjust the matrices during optimisation.

Global parameters
-----------------

These are the elements of the rate matrix.

.. jupyter-execute::

    tabulated["global params"]

Motif parameters
----------------

These are estimates of the nucleotide probabilities in the unobserved ancestor.

.. jupyter-execute::

    tabulated["motif params"]
