.. _custom_apps:

Writing your own apps
=====================

The app infrastructure is provided by the ``scinexus`` package. See the |scinexus| documentation for full details on |define_app|, |app_types|, handling |not_completed|, and |citations|.

Below are cogent3-specific examples showing how to write apps that work with cogent3 data types.

Available cogent3 types
-----------------------

You can use existing type hints if your function takes or returns ``cogent3`` types.

.. jupyter-execute::

    from cogent3.app.typing import defined_types

    defined_types()

.. note:: You don't have to use cogent3 types. You can also use standard Python types.

Defining a cogent3 app from a function
--------------------------------------

We write a function that takes a cogent3 alignment and returns the first *n* positions.

.. jupyter-execute::

    from scinexus.composable import define_app
    from cogent3.app.typing import AlignedSeqsType

    @define_app
    def n_positions(val: AlignedSeqsType, n=2) -> AlignedSeqsType:
        return val[:n]

The critical elements are:

1. The ``define_app`` decorator is used.
2. Type hints are specified for the function's first argument and its return type.

.. dropdown:: Using the custom app

    We create an app instance for a specific value of ``n``.

    .. jupyter-execute::

        first4 = n_positions(n=4)
        first4

    The instance's ``repr()`` indicates the wrapped function and the argument values. You use ``first4()`` like all composable apps, e.g.

    .. jupyter-execute::

        from cogent3 import make_aligned_seqs

        aln = make_aligned_seqs(
            dict(a="GCAAGCGTTTAT", b="GCTTTTGTCAAT"), moltype="dna"
        )
        result = first4(aln)
        result

Defining a cogent3 app from a class
------------------------------------

.. jupyter-execute::

    from scinexus.composable import define_app
    from cogent3.app.typing import AlignedSeqsType

    @define_app
    class n_positions:
        def __init__(self, n=2):
            self.n = n

        def main(self, val: AlignedSeqsType) -> AlignedSeqsType:
            return val[:self.n]

The critical elements are:

1. The ``define_app`` decorator is used.
2. The class has a ``main()`` method.
3. Type hints are specified for the ``main()`` method's first argument and its return type.

.. dropdown:: Using the custom app

    This is identical to what we did above.

    .. jupyter-execute::

        first4 = n_positions(n=4)

        # we use the alignment defined above

        result = first4(aln)
        result

App naming conventions
----------------------

Use words in lower case separated by underscores (e.g. ``lower_case``) to name your apps. Apps are callable, just like functions, and the `PEP8 guidelines <https://peps.python.org/pep-0008/#function-and-variable-names>`_ specify this naming style.

If you will make your app available on the Python package index, we recommend prefixing each app with your package name. For example, the `piqtree2 <https://pypi.org/project/piqtree2>`_ library distributes apps with names such as ``piqtree_phylo``.

.. _app_citations:

See ``scinexus`` documentation for how to add a |citation| to your app.

How to get the citations for the apps you use
---------------------------------------------

Correctly attributing the authors of algorithms and software is a requirement of good scientific practice. The ``.citations`` property on an app instance returns its citations as a tuple.

.. jupyter-execute::

    app = strict_filter()
    app.citations

You can get the BibTeX string directly via the ``.bib`` property.

.. jupyter-execute::

    print(app.bib)

.. dropdown:: Citations in a composed pipeline

    When apps are composed into a pipeline, ``.citations`` collects unique citations from all apps in the chain.

    .. jupyter-execute::

        from cogent3 import get_app

        loader = get_app("load_aligned", moltype="dna", format_name="fasta")
        pipeline = loader + strict_filter()
        pipeline.citations

    The ``.bib`` property gives the combined BibTeX for the whole pipeline.

    .. jupyter-execute::

        print(pipeline.bib)

.. note::

    When a composed pipeline is run via ``apply_to()``, citations are
    automatically saved in the output data store. See |dstore_cites|
    for how to inspect and export them.
