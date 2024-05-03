.. jupyter-execute::
    :hide-code:

    import set_working_directory

Loading tabular data
--------------------

How to load data from a tsv file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can load a file containing tab seperated values with the ``load_tabular`` app, providing it with the appropriate seperator, ``sep="\t"``. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    load_tsv_app = get_app("load_tabular", sep="\t")
    tsv_data = load_tsv_app("data/stats.tsv")
    tsv_data

How to load data from a csv file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can load a file containing comma seperated values with the ``load_tabular`` app, providing it with the seperator, ``sep=","``. 

In the above example, ``data`` is a ``cogent3.util.Table`` object, if we write it to disk in csv format, we can re-use it for this example too. 

.. jupyter-execute::
    :hide-code:

    
    from tempfile import NamedTemporaryFile

    tmpdir = NamedTemporaryFile(dir=".")
    path_to_csv_file = tmpdir.name

.. tip:: If you are executing this code on your machine, replace ``path_to_csv_file`` with the path containing your csv file!

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    tsv_data.write(path_to_csv_file, sep=",")

    load_csv_app = get_app("load_tabular", sep=",")
    csv_data = load_csv_app(path_to_csv_file)
    csv_data