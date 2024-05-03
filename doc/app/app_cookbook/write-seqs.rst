.. jupyter-execute::
    :hide-code:

    import set_working_directory

Writing sequences and sequence alignments
-----------------------------------------

Writing a sequence alignment to disk
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Writing an alignment to disk can be achieved easily with the ``write_seqs`` app. First let's load the alignment we want to write. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    load_aligned_app = get_app("load_aligned", moltype="dna", format="fasta")
    aln = load_aligned_app("data/primate_brca1.fasta")
    aln

When creating the ``write_seqs`` app, we need to provide a data store to which we want the data to be written, and optionally, we can specify the format we want the sequences to be written in. 

.. jupyter-execute::
    :hide-code:

    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, open_data_store

    seq_dstore = open_data_store(path_to_dir, suffix="phylip", mode="w")

    write_seqs_app = get_app("write_seqs", data_store=seq_dstore, format="phylip")
    result = write_seqs_app(aln)

    result.read()[:50]

Writing many sequence collections to a data store
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically, the final step of a data processing pipeline is writing out the filtered data. When ``write_seqs`` is composed into a process, the process will write out multiple sequence collections to a data store. 

.. jupyter-execute::
    :hide-code:

    
    from tempfile import TemporaryDirectory

    tmpdir = TemporaryDirectory(dir=".")
    path_to_dir = tmpdir.name

We can create our input data store containing all the files with the ".fasta" suffix in the data directory using ``open_data_store``. 

.. jupyter-execute::
    :raises:

    from cogent3 import open_data_store

    fasta_seq_dstore = open_data_store("data", suffix="fasta", mode="r")

Let's define a process. In this example, our process loads the sequences, filters the sequences to keep only those which are translateable, translates the sequences, and then writes the filtered sequences to a data store. 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app, open_data_store

    out_dstore = open_data_store(path_to_dir, suffix="fa", mode="w")

    loader = get_app("load_unaligned", format="fasta", moltype="dna")
    keep_translatable = get_app("select_translatable")
    translate = get_app("translate_seqs")
    writer = get_app("write_seqs", out_dstore, format="fasta")

    process = loader + keep_translatable + translate + writer

.. tip:: When running this code on your machine, remember to replace ``path_to_dir`` with an actual directory path.

We apply ``process`` to our input data store, and assign the resulting data store to ``result``. 

.. jupyter-execute::
    :raises:

    result = process.apply_to(fasta_seq_dstore)

Acessing an overview of our process
"""""""""""""""""""""""""""""""""""

We can interrogate ``result`` to see an overview of the process. 

.. jupyter-execute::
    :raises:

    result.describe

There were 10 data files to which the process was successfully applied. However, there were three to which the process did not complete. We can see a summary of the failures by acessing the ``summary_not_completed`` property. 

.. jupyter-execute::
    :raises:

    result.summary_not_completed

Looks like the first two failed because they are protein sequences and ``load_unaligned`` expected DNA sequences. 

Interestingly, another file failed in the ``keep_translatable`` step. By design, these failures did not stop the rest of the pipeline from being run. In fact, the data store collects the :ref:`NotCompleted objects <not_completed>`, which store traceback information, allowing you to interrogate any failings. 
