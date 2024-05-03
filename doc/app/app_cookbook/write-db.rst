.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _write_db:

Writing a database to file
--------------------------

``write_db`` can be used to write serialised objects to a database instance. In the below example, we use it to write the output of a composed process to disk. 

.. jupyter-execute::
    :hide-code:

    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    dstore = open_data_store("data", suffix="fasta", limit=2)

    reader = get_app("load_aligned", format="fasta", moltype="dna")
    min_length = get_app("min_length", 300)
    out_dstore = open_data_store(f"{path_to_dir}.sqlitedb", mode="w")
    writer = get_app("write_db", out_dstore) 

    process = reader + min_length + writer
    result = process.apply_to(dstore)
    result.describe