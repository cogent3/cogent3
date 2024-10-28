.. _available_apps:

Displaying installed apps
-------------------------

Apps perform functions ranging from multiple sequence alignment (e.g. ``progressive_align``), to excluding alignment columns containing non-nucleotide characters (e.g. ``omit_degenerates``) to performing maximum-likelihood evolutionary analyses (e.g. ``model``).

.. jupyter-execute::

    from cogent3 import available_apps

    available_apps()

The ``name_filter`` argument can be used to display only the apps that match a string. For example, to display all the loader apps

.. jupyter-execute::

    available_apps(name_filter="load")


.. note:: The first column lists the package providing the app.