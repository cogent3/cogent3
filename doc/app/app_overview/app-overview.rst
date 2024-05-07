.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _apps:

Overview of using apps
======================

What are apps?
--------------

Apps are ready-made "functions" that you can apply to your data without needing to know all the technical details. They are easy to use, even if you're not an expert programmer. Multiple apps can be naturally composed into "pipelines", which are fully equipped for robust and reproducible application to data. In fact, apps and app pipelines can be applied to a single, or thousands, of data file(s) without writing loops or conditionals.

Apps have several key features, they:

#. trap errors
#. check the validity of input data
#. can be used by themselves or combined into a "composed" app (aka pipeline)
#. automatically track the relationship between an input data record and its output record

Where to start?
---------------

Three top-level functions are very useful:

- ``available_apps()`` identifies what apps are installed (see :ref:`available_apps`)
- ``app_help()`` shows what a given app can do (see :ref:`app_help`)
- ``get_app()`` returns an app instance for you to use (see :ref:`get_app`)

Two other crucial concepts concern: 

- :ref:`data stores <data_stores>` and 
- :ref:`tracking failures <not_completed>`

.. _app_types:

Types of apps
-------------

There are 3 types of apps:

#. loaders (by convention, names starts with ``load_<data type>``)
#. writers (by convention, names starts with ``write_<data type>``)
#. generic (no naming convention)

As their names imply, loaders load, writers write and generic apps do other operations on data.

.. _app_composability:

Composability
-------------

Most ``cogent3`` apps are "composable", meaning that multiple apps can be combined into a single function by addition. For example, say we have an app (``fit_model``) that performs a molecular evolutionary analysis on an alignment, and another app (``extract_stats``) that gets the statistics from the result. We could perform these steps sequentially as follows

.. code-block:: python

    fitted = fit_model(alignment)
    stats = extract_stats(fitted)

Composability allows us to simplify this as follows

.. code-block:: python

    app = fit_model + extract_stats
    stats = app(fitted)

We can have many more apps in a composed function than just the two shown here.

.. _composability_rules:

Composability rules
^^^^^^^^^^^^^^^^^^^

There are rules around app composition, starting with app types. Loaders and writers are special cases. If included, a loader must always be first, e.g.

.. code-block:: python

    app = a_loader + a_generic

If included, a writer must always be last, e.g.

.. code-block:: python

    app = a_generic + a_writer

Changing the order for either of the above will result in a ``TypeError``.

The next constraint on app composition are the input and output types of the apps involved. Specifically, apps define the type of input they work on and the type of output they produce. For two apps to be composed, the output (or return) type of app on the left (e.g. ``a_loader``) must overlap with the input type of the app on the right (e.g. ``a_generic``). If they don't match, a ``TypeError`` is raised.

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

    loader = get_app("load_aligned", format="fasta", moltype="dna")
    cpos3 = get_app("take_codon_positions", 3)
    writer = get_app("write_seqs", out_dstore, format="fasta")

There are two ways in which I can apply the three above apps to data:

1. Using apps sequentially like functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    data = loader("data/primate_brca1.fasta")
    just3rd = cpos3(data)
    m = writer(just3rd)

The resulting alignment ``just3rd`` will be written into the ``out_dstore`` directory in fasta format with the same filename as the original data (``"primate_brca1.fasta"``).

.. note::

    ``m`` is a :ref:`DataMember <data_member>` of ``out_dstore``.

2. Composing several apps into a multi-step "process"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can make this simpler by creating a single composed function.

.. jupyter-execute::

    process = loader + cpos3 + writer
    m = process("data/primate_brca1.fasta")

Applying a process to multiple data records
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To apply a composed function to multiple files requires a :ref:`data store <data_stores>`. Using ``open_data_store()`` we identify all data files in a directory that we want to analyse, in the following case, all fasta file in the data directory. ``process`` can be then applied to all records in the data store without having to loop.

.. jupyter-execute::

    dstore = open_data_store("data", suffix="fasta", mode="r")
    result = process.apply_to(dstore)

.. note:: ``result`` is ``out_dstore``.

Other important features
------------------------

The settings and data analysed will be logged
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A log file will be written into the same data store as the output. The log includes information on the conditions under which the analysis was run and fingerprint all input and output files.

.. jupyter-execute::

    out_dstore.summary_logs

Failures are recorded
^^^^^^^^^^^^^^^^^^^^^

Any "failures" (see :ref:`not_completed`) are saved. The data store class provides methods for interrogating those. First, a general summary of the output data store indicates we have 6 records that did not complete.

.. jupyter-execute::

    out_dstore.describe

These occur for this example primarily because some of the files contain sequences that are not aligned

.. jupyter-execute::

    out_dstore.summary_not_completed

You can track progress
^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    result = process.apply_to(dstore, show_progress=True)

You can do parallel computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    result = process.apply_to(dstore, parallel=True)

By default, this will use all available processors on your machine. (See :ref:`parallel` for more details plus how to take advantage of multiple machines using MPI.)

All of the above
^^^^^^^^^^^^^^^^

.. code-block:: python

    process.apply_to(dstore, parallel=True, show_progress=True)

.. jupyter-execute::
    :hide-code:

    import shutil

    shutil.rmtree(path_to_dir, ignore_errors=True)
