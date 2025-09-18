.. jupyter-execute::
    :hide-code:

    import set_working_directory

Using iplotx with cogent3
=========================
`iplotx <https://https://iplotx.readthedocs.io>`__ is a library to visualise trees and networks
using ``matplotlib`` (``cogent3`` uses ``plotly`` internally). It supports dozens of options
to style the appearance of trees and can produce static images, interactive plots, animations,
and so on. It can also be used to combine trees and other types of charts (e.g. bar charts,
scatter plots, networks, heatmaps, etc.) in the same figure.

Below is a simple example on how to combine ``cogent3`` with ``iplotx``:

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

    import matplotlib.pyplot as plt
    plt.savefig(outpath)
