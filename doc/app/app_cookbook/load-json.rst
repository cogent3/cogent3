.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _load_json:

Loading JSON serialised objects
-------------------------------

The ``load_json`` app loads JSON serialised ``cogent3`` objects from a file, returning whatever was stored.  

In the :ref:`writing JSON section <write_json>`, we wrote a likelihood function object to disk in JSON serialised format, now let's load it back! All we need is the path to which it was written. Since ``write_json`` writes to a data store, this will be ``<path_to_dstore>/<identifier>``. 

.. tip::  To follow along with this example, first head to the :ref:`writing JSON section <write_json>` to write out this data, or replace the path with a JSON file on your machine. 

.. jupyter-execute::
    :hide-code:
    
    from cogent3 import get_app, open_data_store
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

    # create lf object to write
    load_aligned_app = get_app("load_aligned", moltype="dna", format="fasta")
    aln = load_aligned_app("data/primate_brca1.fasta")
    gn_model_app = get_app("model", "GN", tree="data/primate_brca1.tree")
    lf = gn_model_app(aln).lf

    out_dstore = open_data_store(path_to_dir, mode="w", suffix="json")

    # write with the same identifier as the write_json example
    write_json_app = get_app("write_json", data_store=out_dstore)
    data_member = write_json_app(lf, identifier="gn_params.json")
    path_to_dstore = data_member.data_store.source

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    load_json_app = get_app("load_json")
    lf = load_json_app(f"{path_to_dstore}/gn_params.json")
    lf
