.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _write_json:

Writing JSON Serialised Objects
-------------------------------

Using JSON, we can serialise ``cogent3`` objects to a file for easy storage and retrieval.

Create an example object to serialise
"""""""""""""""""""""""""""""""""""""

Lets create a ``LikelihoodFunction`` object to use in this example. It is generated from fitting the General Nucleotide (GN) model to an alignment of BRCA1 in primates. 

.. jupyter-execute::
    :hide-code:

    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    # Load the alignment
    load_aligned_app = get_app("load_aligned", moltype="dna", format="fasta")
    aln = load_aligned_app("data/primate_brca1.fasta")

    # Fit the GN model
    gn_model_app = get_app("model", "GN", tree="data/primate_brca1.tree")
    model_result = gn_model_app(aln)
    model_result.lf

``write_json`` -  writing JSON-serialised object to file 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Using the ``write_json`` app, we can write out the likelihood function as a JSON-serialised object, making it easy to retrieve the model parameters for future reference if required. 

We need to provide the ``write_json`` app with a data store to which it will write to. Optionally when we apply the app we can specify a identifier for the data, which will name the file. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    # Initialise the write_json app with a data store
    out_dstore = open_data_store(path_to_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)

    write_json_app(model_result.lf, identifier="gn_params.json")

.. note:: Learn how to load a JSON serialised object in the :ref:`loading JSON serialised objects <load_json>` section!

.. tip:: When running this code on your machine, remember to replace ``path_to_dir`` with an actual directory path.
