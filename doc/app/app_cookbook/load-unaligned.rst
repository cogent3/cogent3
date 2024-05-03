.. jupyter-execute::
    :hide-code:

    import set_working_directory

Loading unaligned sequence data
-------------------------------

We can load unaligned sequence data using the ``load_unaligned`` app, this will return a ``SequenceCollection``. 

Loading unaligned DNA sequences from a single fasta file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we load unaligned DNA sequences from a single fasta file using the load_unaligned app. We specify the molecular type ``(moltype="protein")`` and the file format ``(format="fasta")``.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    load_unaligned_app = get_app("load_unaligned", format="fasta", moltype="protein")
    seqs = load_unaligned_app("data/inseqs_protein.fasta")
    seqs

Loading unaligned DNA sequences from multiple fasta files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load unaligned DNA sequences from multiple fasta files, we need two things, a data store that identifies the files we are interested in and a process composed of our apps of interest. 

1. A data store that identifies the files we are interested in 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Here we open a read-only (``mode="r"``) data store that identifies all fasta files in the data directory, limiting the data store to two members as a minimum example.

.. jupyter-execute::
    :hide-code:

    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    fasta_seq_dstore = open_data_store("data", suffix="fasta", mode="r", limit=2)

2. A composed process that defines our workflow 
"""""""""""""""""""""""""""""""""""""""""""""""

In this example our process loads the unaligned sequences using ``load_unaligned``, then applies ``jaccard_dist`` to estimate a kmer based genetic distance, which we write out to a data store using ``write_tabular``. 

.. note:: Apps that are "writers" typically require a data store to write to, learn more about writers :ref:`here! <writers>`. 

.. jupyter-execute::
    :raises:

    out_dstore = open_data_store(path_to_dir, suffix="tsv", mode="w")
    
    load_unaligned_app = get_app("load_unaligned", format="fasta", moltype="dna")
    jdist = get_app("jaccard_dist")
    writer = get_app("write_tabular", out_dstore, format="tsv") 

    process = load_unaligned_app + jdist + writer

.. tip:: When running this code on your machine, remember to replace ``path_to_dir`` with an actual directory path.

Now we're good to go! We can apply ``process`` to our data store of fasta sequences. ``result`` is a data store, which you can index to see individual data members. We can inspect a given data member look using the ``.read()`` on data members. 

.. jupyter-execute::
    :raises:

    result = process.apply_to(fasta_seq_dstore)
    print(result[1].read())