.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _apps:

Overview of using apps
======================

What are apps?
--------------

Apps are ready-made "functions" that you can apply to your data without needing to know all the technical details. They are easy to use, even if you're not an expert programmer. Multiple apps can be naturally composed into "pipelines", which are fully equipped for robust and reproducible application to data. In fact, apps and app pipelines can be applied to a single, or thousands, of data file(s) without writing loops or conditionals.

The composable app infrastructure is provided by the ``scinexus`` package. See the |scinexus| documentation for details on composability rules, type checking, batch processing, parallel execution, and progress tracking.

.. _app_start:

How do I start to use ``cogent3`` apps?
---------------------------------------

Three top-level functions are very useful:

- ``available_apps()`` identifies what apps are installed (see :ref:`available_apps`)
- ``app_help()`` shows what a given app can do (see :ref:`app_help`)
- ``get_app()`` returns an app instance for you to use (see :ref:`get_app`)

Two other crucial concepts concern:

- |data_stores|
- |track_failures|

Types of apps
-------------

There are 3 |app_types|:

#. loaders (by convention, names starts with ``load_<data type>``)
#. writers (by convention, names starts with ``write_<data type>``)
#. generic (no naming convention)

As their names imply, loaders load, writers write and generic apps do other operations on data.

An example
----------

.. jupyter-execute::
    :hide-code:


    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

I illustrate the general approach for a simple example -- extracting third codon positions. As I'm defining a writer, I also need to define the destination (a directory in this case) where it will write to.

.. jupyter-execute::

    from cogent3 import get_app, open_data_store

    out_dstore = open_data_store(path_to_dir, suffix="fa", mode="w")

    loader = get_app("load_aligned", format_name="fasta", moltype="dna")
    cpos3 = get_app("take_codon_positions", 3)
    writer = get_app("write_seqs", out_dstore, format_name="fasta")

There are two ways in which I can apply the three above apps to data:

1. Using apps sequentially like functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    data = loader("data/primate_brca1.fasta")
    just3rd = cpos3(data)
    m = writer(just3rd)

The resulting alignment ``just3rd`` will be written into the ``out_dstore`` directory in fasta format with the same filename as the original data (``"primate_brca1.fasta"``).

.. note::

    ``m`` is a |data_member| of ``out_dstore``.

2. Composing several apps into a multi-step "process"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can make this simpler by creating a single composed function.

.. jupyter-execute::

    process = loader + cpos3 + writer
    m = process("data/primate_brca1.fasta")

Applying a process to multiple data records
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To apply a composed function to multiple files requires a |data_stores|. Using ``open_data_store()`` we identify all data files in a directory that we want to analyse, in the following case, all fasta file in the data directory. ``process`` can be then applied to all records in the data store without having to loop.

.. jupyter-execute::

    dstore = open_data_store("data", suffix="fasta", mode="r")
    result = process.apply_to(dstore)

.. note:: ``result`` is ``out_dstore``.

.. jupyter-execute::
    :hide-code:

    import shutil

    shutil.rmtree(path_to_dir, ignore_errors=True)
