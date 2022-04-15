.. jupyter-execute::
    :hide-code:

    import set_working_directory

Applying a discrete-time, non-stationary nucleotide model
=========================================================

We fit a discrete-time Markov nucleotide model. This corresponds to a Barry and Hartigan 1987 model.

.. jupyter-execute::

    from cogent3.app import evo, io

    loader = io.load_aligned(format="fasta", moltype="dna")
    aln = loader("data/primate_brca1.fasta")
    model = evo.model("BH", tree="data/primate_brca1.tree")
    result = model(aln)
    result

**NOTE:** DLC stands for diagonal largest in column and the value is a
check on the identifiability of the model. ``unique_Q`` is not
applicable to a discrete-time model and so remains as ``None``.

Looking at the likelihood function, you will

.. jupyter-execute::

    result.lf

Get a tree with branch lengths as paralinear
--------------------------------------------

This is the only possible length metric for a discrete-time process.

.. jupyter-execute::

    tree = result.tree
    fig = tree.get_figure()
    fig.scale_bar = "top right"
    fig.show(width=500, height=500)

Getting parameter estimates
---------------------------

For a discrete-time model, aside from the root motif probabilities, everything is edge specific. But note that the ``tabular_result`` has different keys from the continuous-time case, as demonstrated below.

.. jupyter-execute::

    tabulator = evo.tabulate_stats()
    stats = tabulator(result)
    stats

.. jupyter-execute::

    stats["edge motif motif2 params"]
