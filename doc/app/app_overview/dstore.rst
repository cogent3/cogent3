.. _data_stores:

Data stores -- collections of data records
==========================================

Data stores are provided by the `scinexus <https://pypi.org/project/scinexus/>`_ package. A :index:`data store` is a collection of data members of the same type (e.g. all ``.fasta`` files in a directory). Data stores allow you to apply an app or composed pipeline to many data records without writing loops.

.. todo::

    Link to scinexus data store documentation once available.

..
    Placeholder: <scinexus data store docs URL>

Using data stores with cogent3
------------------------------

Use the ``open_data_store()`` function to open a data store and a :ref:`loader <app_types>` app to read member data.

.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. jupyter-execute::

    from cogent3 import get_app, open_data_store

    dstore = open_data_store("data/raw.zip", suffix="fa", mode="r")
    print(dstore)

.. jupyter-execute::

    loader = get_app("load_unaligned", moltype="dna")
    seqs = loader(dstore[0])
    seqs

.. _data_member:

Applying a pipeline to a data store
------------------------------------

.. jupyter-execute::
    :hide-code:

    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

.. jupyter-execute::

    out_dstore = open_data_store(path_to_dir, suffix="fa", mode="w")
    loader = get_app("load_aligned", moltype="dna", format_name="fasta")
    take3 = get_app("take_codon_positions", 3)
    writer = get_app("write_seqs", data_store=out_dstore, format_name="fasta")
    app = loader + take3 + writer
    result = app.apply_to(dstore)
    result.describe

.. jupyter-execute::
    :hide-code:

    import shutil
    shutil.rmtree(path_to_dir, ignore_errors=True)

.. _data_store_citations:

See the scinexus documentation for full details on data store types, structure, operations, locking, logging, and citations.
