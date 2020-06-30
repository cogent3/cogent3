What apps are there?
--------------------

Apps perform functions ranging from multiple sequence alignment (e.g. ``progressive_align``), to excluding alignment columns containing non-nucleotide characters (e.g. ``omit_degenerates``) to performing maximum-likelihood evolutionary analyses (e.g. ``model``).

Apps that are identified as “composable” (with a value of ``True`` under that column) can be combined into a single function by addition.

.. jupyter-execute::

    from cogent3 import available_apps

    available_apps()
