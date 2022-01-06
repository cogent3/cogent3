.. jupyter-execute::
    :hide-code:

    import set_working_directory

Radial Dendrogram Style
=======================

.. jupyter-execute::

    from cogent3.app import io

    reader = io.load_json()

    ens_tree = reader("data/GN-tree.json")
    fig = ens_tree.get_figure("radial", width=600, height=600)
    fig.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_tree-radial.png"

    fig.write(outpath)

With Contemporaneous Tips
-------------------------

.. jupyter-execute::

    fig.contemporaneous = True
    fig.label_pad = 0.23
    fig.show()
