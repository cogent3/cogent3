.. jupyter-execute::
    :hide-code:

    import set_working_directory

Reconstructing ancestral states
-------------------------------

This app takes a ``model_result`` and returns a ``tabular_result`` consisting of the posterior probabilities of ancestral states for each node of a tree. These probabilities are computed using the marginal reconstruction algorithm.

We first fit a model to the sample data.

.. jupyter-execute::

    from cogent3 import get_app

    loader = get_app("load_aligned", format="fasta")
    aln = loader("data/primate_brca1.fasta")
    gn = get_app("model", "GN", tree="data/primate_brca1.tree")
    result = gn(aln)

Define the ``ancestral_states`` app
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    reconstuctor = get_app("ancestral_states")
    states_result = reconstuctor(result)
    states_result

The ``tabular_result`` is keyed by the node name. Each value is a ``DictArray``, with header corresponding to the states and rows corresponding to alignment position.

.. jupyter-execute::

    states_result["edge.0"]

If not included in the newick tree file, the internal node names are automatically generated open loading. You can establish what those are by interrogating the tree bound to the likelihood function object. (If you move your mouse cursor over the nodes, their names will appear as hover text.)

.. jupyter-execute::

    result.tree.get_figure(contemporaneous=True).show(width=500, height=500)
