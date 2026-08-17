.. jupyter-execute::
    :hide-code:

    import set_working_directory

Angular Dendrogram Style
========================

This is a left-right style. You'll note that there's overlap of edges at the bottom -- a known issue with this display style.

.. jupyter-execute::

    import cogent3

    reader = cogent3.get_app("load_json")

    ens_tree = reader("data/GN-tree.json")
    fig = ens_tree.get_figure(style="angular", width=600, height=600)
    fig.show()

With Contemporaneous Tips
-------------------------

.. jupyter-execute::

    fig.contemporaneous = True
    fig.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_tree-angular.png"

    fig.write(outpath)
