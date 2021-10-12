.. _apps:

**********************
Overview of using apps
**********************

The general idea is you create the flavour of function you want from an existing app. What you have created is then callable. It can be applied to a single data object (like an alignment file), or a series of them  (like a directory of alignment files). Alternatively, an app can be combined with other apps to make a pipeline.

I illustrate the general approach for a simple example -- extracting third codon positions.

Define the apps
---------------

.. code-block:: python

    from cogent3.app import io, sample

    reader = io.load_aligned(format="fasta", moltype="dna")
    cpos3 = sample.take_codon_positions(3)
    writer = io.write_seqs("path/to/write/thirdpos", format="fasta")

Using apps like functions
-------------------------

.. code-block:: python

    data = reader("some/path/to/seqs.fasta")
    just3rd = cpos3(data)
    m = writer(just3rd, identifier="3rdpos_data.fasta")

In the above, ``m`` is a ``DataStoreMember``. The result will be written into the directory specified in constructing the ``writer``.

Composing a multi-step process from several apps
------------------------------------------------

The above can be simplified by creation of a single composed function. Executing this with the path argument will generate the same output.

.. code-block:: python

    process = reader + cpos3 + writer
    m = process("some/path/to/seqs.fasta")

Applying a process to multiple data records
-------------------------------------------

We use a data store to identify all data files in a directory that we want to analyse. ``process`` can be then applied to all records in the data store.

.. code-block:: python

    dstore = io.get_data_store("path/to/dir", suffix="fasta")
    r = process.apply_to(dstore)

Here ``r`` is a list of all the ``DataStoreMember`` instances.

Other important features
------------------------

You can track progress
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    process.apply_to(dstore, show_progress=True)

You can do parallel computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    process.apply_to(dstore, parallel=True)

By default, this will use all available processors on your machine. If you are running in an mpi environment, you can add the argument ``par_kw=dict(use_mpi=True)``. For more details, see :ref:`parallel`.

You can log the settings and data analysed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    process.apply_to(dstore, logger=True)

All of the above
^^^^^^^^^^^^^^^^

.. code-block:: python

    process.apply_to(dstore, parallel=True, logger=True, show_progress=True)

If you use the ``json`` based output formats (either explicitly, or via using the tinydb data store type), any "failures" (see :ref:`not_completed`) will be written to file also and a convenient interface is provided for interrogating those.
