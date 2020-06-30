****************************************
How does ``cogent3`` relate to PyCogent?
****************************************

cogent3_ is a significantly changed library from the original PyCogent. The renaming has been done to emphasise these differences and to make the project name and import statement consistent (``cogent`` was always the import name, originating in the ``pyevolve`` project from 2004).

Most of the changes from PyCogent involved elimination of modules, using `black <https://github.com/psf/black>`_ and `isort <https://github.com/timothycrosley/isort>`_ for coding style, rationalisation of interfaces and the addition of new features. For instance, we have an experimental ``cogent3.app`` module (documentation still being written) that is intended to present a functional programming style interface to ``cogent3`` capabilities. We now also use `Plotly <https://plotly.com/python/>`_ for all visualisation.

The rewrite has been a massive amount of work and unfortunately the changes to the API are only indirectly documented by virtue of having the documentation match the library state. Thus, the best way to get older scripts working is to check the Library documentation related to your code. More explicitly, you can also search in the `repository history <https://github.com/cogent3/cogent3>`_.

``cogent3`` no longer includes module ``x``, what do I do?
==========================================================

.. glossary::
    ``cogent.app``
        The ``PyCogent`` module was concerned with wrapping external applications. There are multiple 3rd party alternatives to this, for example ``click``, ``burrito``, etc.. The ``cogent3.app`` module is very different being focussed on providing a functional style interface to ``cogent3`` capabilities.

    ``cogent.db.ensembl``
        This has been turned into a standalone project, `EnsemblDb3 <https://github.com/cogent3/ensembldb3>`_.

    ``cogent.db.eutils``
        For querying NCBI via their EUtils interface. A replacement is `biocommons/eutils <https://github.com/biocommons/eutils>`_.

    ``cogent.struct``
        For manipulating 3D structures. We know of no replacements.

    ``cogent.motif``
        A subset of k-mer analyses. No direct replacements.

    ``cogent.seqsim``
        See the sequence simulation capabilities in ``cogent3.evolve`` and ``cogent3.app.evo``.

    ``cogent.maths.unifrac``
        The microbiome related functions are now in `scikit-bio <https://scikit-bio.org>`_.

.. _cogent3: https://github.com/cogent3/cogent3
