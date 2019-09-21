*******
Warning
*******

Be warned that cogent3_ is a significantly changed library from the original `PyCogent <http://www.pycogent.org>`_. The renaming has been done to emphasise these differences and to make the project name and import statement consistent (``cogent`` was always the import name, originating in the ``pyevolve`` project from 2004).

Most of the changes from PyCogent involved elimination of modules, using `black <https://github.com/psf/black>`_ and `isort <https://github.com/timothycrosley/isort>`_ for coding style, rationalisation of interfaces and the addition of new features. For instance, we have an experimental ``cogent3.app`` module (documentation still being written) that is intended to present a functional programming style interface to ``cogent3`` capabilities. We now also use `Plotly <https://plot.ly/python/>`_ for all visualisation.

The rewrite has been a massive amount of work and unfortunately the changes to the API are only indirectly documented by virtue of having the documentation match the library state. Thus, the best way to get older scripts working is to check the Library documentation related to your code. More explicitly, you can also search in the `repository history <https://github.com/cogent3/cogent3>`_.

.. _cogent3: https://github.com/cogent3/cogent3