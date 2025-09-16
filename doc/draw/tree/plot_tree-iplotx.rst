.. jupyter-execute::
    :hide-code:

    import set_working_directory
    import matplotlib.pyplot as plt

Using iplotx with cogent3
=========================

.. jupyter-execute::

    import cogent3
    import iplotx as ipx

    reader = cogent3.get_app("load_json")

    ens_tree = reader("data/GN-tree.json")

    ipx.tree(
      ens_tree,
      layout="radial",
      leaf_deep=True,
      leaf_labels=True,
    )

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_tree-radial.png"

    plt.savefig(outpath)
