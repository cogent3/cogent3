.. jupyter-execute::
    :hide-code:

    import set_working_directory

Writing tabular data
--------------------

With the ``write_tabular`` app, ``cogent3`` "TabularTypes" (``Table``, ``DictArray``, ``DistanceMatrix``) are supported for writing to disk. 

Let's generate a ``cogent3`` ``Table`` to use in the examples below. One way to do that is by applying the ``tabulate_stats`` app to a model result. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    # load alignment
    load_aligned_app = get_app("load_aligned", moltype="dna")
    aln = load_aligned_app("data/primate_brca1.fasta")

    # fit GN model
    gn_model_app = get_app("model", "GN", tree="data/primate_brca1.tree")
    model_result = gn_model_app(aln)

    # tabulate the model result
    tabulator = get_app("tabulate_stats")
    model_result_tab = tabulator(model_result)
    motif_params = model_result_tab["motif params"]
    print(type(motif_params))
    motif_params

Writing a CSV file
^^^^^^^^^^^^^^^^^^

To write in CSV format, we create the ``write_tabular`` app with ``format="csv"``. 

.. jupyter-execute::
    :hide-code:
    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name


.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    out_dstore = open_data_store(path_to_dir, mode="w", suffix="csv")

    write_tabular_app = get_app("write_tabular", data_store=out_dstore, format="csv")
    write_tabular_app(motif_params, identifier="gn_model_results.csv")


Writing a TSV file
^^^^^^^^^^^^^^^^^^

To write in TSV format, we create the ``write_tabular`` app with ``format="tsv"``. 

.. jupyter-execute::
    :hide-code:

    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    out_dstore = open_data_store(path_to_dir, mode="w", suffix="tsv")

    write_tabular_app = get_app("write_tabular", data_store=out_dstore, format="tsv")
    write_tabular_app(motif_params, identifier="gn_model_results.tsv")

Using ``write_tabular`` in a composed process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of applying the apps sequentially as above, we can add apps into a composed process, and apply the process to a data store. In this example, we define a process that calculates an unaligned distance measure between sequences, writing these estimated distances to a tsv file. 

.. jupyter-execute::
    :hide-code:

    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    loader = get_app("load_unaligned", moltype="dna")
    jdist = get_app("jaccard_dist")
    out_dstore = open_data_store(path_to_dir, mode="w", suffix="tsv")
    writer = get_app("write_tabular", data_store=out_dstore, format="tsv")

    process = loader + jdist + writer

    in_dstore = open_data_store("data", suffix="fasta", mode="r", limit=2)

    result = process.apply_to(in_dstore)
    result.describe

.. tip:: When running this code on your machine, remember to replace ``path_to_dir`` with an actual directory path.