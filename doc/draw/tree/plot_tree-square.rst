.. jupyter-execute::
    :hide-code:

    import set_working_directory

Square Dendrogram Style
=======================

We use a tree saved in ``json`` format from a likelihood function analysis of a non-stationary model. The tree was derived such that the the branch lengths are now "ENS". 

.. note:: I change the scale bar placement. Valid values are ``"top left"``, ``"top right"``, ``"bottom left"`` (default), or ``"bottom right"``.

.. jupyter-execute::    

    from cogent3.app import io

    reader = io.load_json()

    ens_tree = reader("data/GN-tree.json")
    fig = ens_tree.get_figure(width=600, height=600)
    fig.scale_bar = "top right"
    fig.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_tree-square.png"

    fig.write(outpath)

Colouring a set of edges
------------------------

.. jupyter-execute::

    fig.style_edges(
        "AfricanEl",
        tip2="Manatee",
        legendgroup="Afrotheria",
        line=dict(color="magenta"),
    )
    fig.show()

With Contemporaneous Tips
-------------------------

.. jupyter-execute::

    fig.contemporaneous = True
    fig.show()
