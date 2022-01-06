.. jupyter-execute::
    :hide-code:

    import set_working_directory

Circular Dendrogram Style
=========================


.. jupyter-execute::


    from cogent3.app import io


    reader = io.load_json()

    ens_tree = reader("data/GN-tree.json")
    fig = ens_tree.get_figure("circular", width=600, height=600)
    fig.show()

Colouring a set of edges
------------------------

.. jupyter-execute::

    fig.style_edges("AfricanEl", tip2="Manatee", legendgroup="Afrotheria",
                line=dict(color="magenta", width=2))
    fig.show()


With Contemporaneous Tips
-------------------------

.. jupyter-execute::

    fig.contemporaneous = True
    fig.label_pad = 0.23
    fig.show(width=650, height=600)

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_tree-circular.png"

    fig.write(outpath)
