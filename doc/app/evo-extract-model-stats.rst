.. jupyter-execute::
    :hide-code:

    import set_working_directory

Extracting maximum likelihood estimates from a ``model_result``
===============================================================

If you want to get the stats out-of a fitted model, use the ``evo.tabulate_stats()`` app.

We first fit a model.

.. jupyter-execute::

    from cogent3.app import evo, io

    loader = io.load_aligned(format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")
    model = evo.model("GN", tree="data/primate_brca1.tree")
    result = model(aln)

Create and apply ``tabulate_stats`` app
---------------------------------------

.. jupyter-execute::

    tabulator = evo.tabulate_stats()
    tabulated = tabulator(result)
    tabulated

``tabulated`` is a ``tabular_result`` instance which, like other result types, has ``dict`` like behaviour. It also contains key/value pairs for each model parameter type.

Edge parameters
---------------

These are all parameters that differ between edges. Since the current model is time-homogeneous (a single rate matrix), only the table only has entries for the branch scalar (denoted “length”).

.. jupyter-execute::

    tabulated["edge params"]

.. note:: Unless the model is time-reversible, the lengths in that table are not ENS (`Kaehler et al <https://www.ncbi.nlm.nih.gov/pubmed/28175284>`__). As we used a non-stationary nucleotide model in this example, the length values are a scalar used to adjust the matrices during optimisation.

Global parameters
-----------------

In this example, these are the elements of the rate matrix.

.. jupyter-execute::

    tabulated["global params"]

Motif parameters
----------------

In the current example, these are estimates of the nucleotide probabilities in the unobserved ancestor.

.. jupyter-execute::

    tabulated["motif params"]
