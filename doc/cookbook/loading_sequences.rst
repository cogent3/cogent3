.. jupyter-execute::
    :hide-code:

    import os

    import set_working_directory

.. _load_seq:

Loading a single sequence from a file
-------------------------------------

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

In this case, the filename suffix is used to infer the data format.

.. jupyter-execute::

    from cogent3 import load_seq

    seq = load_seq("data/mycoplasma-genitalium.fa", moltype="dna")
    seq

.. warning:: If a file has more than one sequence, only the first one is loaded.

.. jupyter-execute::

    seq = load_seq("data/brca1-bats.fasta", moltype="dna")
    seq

.. note:: It's also possible to load a sequence from a :ref:`url <load_url>`.

.. _load-seqs:

Loading an alignment from a file or url
---------------------------------------

.. author, Gavin Huttley, Tony Walters, Tom Elliott

Loading aligned sequences
^^^^^^^^^^^^^^^^^^^^^^^^^

Any file in which the sequences have exactly the same length can be loaded as an alignment.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    type(aln)

.. note:: The load functions record the origin of the data in the ``info`` attribute under a `"source"` key.

.. jupyter-execute::

    aln.info.source

.. todo:: add cross ref for description of Info class

Loading unaligned sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Files containing sequences that may differ in length can be loaded using ``load_unaligned_seqs()``, which returns a sequence collection.

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs

    seqs = load_unaligned_seqs("data/long_testseqs.fasta", moltype="dna")
    type(seqs)

.. _load_url:

Loading from a url
^^^^^^^^^^^^^^^^^^

The ``cogent3`` load functions support loading from a url. We load the above fasta file directly from GitHub.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs(
        "https://raw.githubusercontent.com/cogent3/cogent3/develop/doc/data/long_testseqs.fasta",
        moltype="dna",
    )

Specifying the file format
^^^^^^^^^^^^^^^^^^^^^^^^^^

The loading functions use the filename suffix to infer the file format. This can be overridden using the ``format`` argument.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna", format="fasta")
    aln

Specifying the sequence molecular type
--------------------------------------

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    protein_seqs = {"seq1": "DEKQL-RG", "seq2": "DDK--SRG"}
    proteins_loaded = make_aligned_seqs(protein_seqs, moltype="protein")
    proteins_loaded.moltype
    proteins_loaded

Making an alignment from standard python objects
------------------------------------------------

From a dict of strings
^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    seqs = {"seq1": "AATCG-A", "seq2": "AATCGGA"}
    seqs_loaded = make_aligned_seqs(seqs, moltype="dna")

From a dict of numpy arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_aligned_seqs
    from numpy import array, uint8

    seqs = {
        "seq1": array([2, 2, 0, 1, 3, 9, 2], dtype=uint8),
        "seq2": array([2, 2, 0, 1, 3, 3, 2], dtype=uint8),
    }
    seqs_loaded = make_aligned_seqs(seqs, moltype="dna")

From a series of strings
^^^^^^^^^^^^^^^^^^^^^^^^

The sequence names will be automatically created.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = ["AATCG-A", "AATCGGA"]
    coll = make_aligned_seqs(data, moltype="dna", new_type=True)
    coll

Changing sequence labels on loading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Load a list of aligned nucleotide sequences, while specifying the DNA molecule type and stripping the comments from the label. In this example, we rename sequences by passing a function that removes everything after the first whitespace to the ``label_to_name`` parameter.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = {
        "sample1 Mus musculus": "AACCTGC--C",
        "sample2 Gallus gallus": "AAC-TGCAAC",
    }
    loaded_seqs = make_aligned_seqs(
        data, moltype="dna", label_to_name=lambda x: x.split()[0]
    )
    loaded_seqs

Making a sequence collection from standard python objects
---------------------------------------------------------

This is done using ``make_unaligned_seqs()``, which returns a ``SequenceCollection`` instance. The function arguments match those of ``make_aligned_seqs()``. We demonstrate only for the case where the input data is a ``dict``.

.. jupyter-execute::

    from cogent3 import make_unaligned_seqs

    seqs = {"seq1": "AATCA", "seq2": "AATCGGA"}
    seqs = make_unaligned_seqs(data=seqs, moltype="dna")
    seqs

Loading sequences using format parsers
--------------------------------------

``load_aligned_seqs()`` and ``load_unaligned_seqs()`` are just convenience interfaces to format parsers. It can sometimes be more effective to use the parsers directly, say when you don't want to load everything into memory.

Loading FASTA sequences from an open file or list of lines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load FASTA formatted sequences directly, you can use ``iter_fasta_records``. This parser returns data as python strings.

.. note:: This returns the sequences as strings.

.. jupyter-execute::

    from cogent3.parse.fasta import iter_fasta_records

    seqs = list(iter_fasta_records("data/long_testseqs.fasta"))
    seqs

Handling overloaded FASTA sequence labels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The FASTA label field is frequently overloaded, with different information fields present in the field and separated by some delimiter. This can be flexibly addressed using the ``LabelParser``. By creating a custom label parser, we can decide which part we use as the sequence name. We show how to convert a field into something specific.

.. jupyter-execute::

    from cogent3.parse.fasta import LabelParser

    def latin_to_common(latin):
        return {"Homo sapiens": "human", "Pan troglodtyes": "chimp"}[latin]

    label_parser = LabelParser(
        "%(species)s", [[1, "species", latin_to_common]], split_with=":"
    )
    for label in ">abcd:Homo sapiens:misc", ">abcd:Pan troglodtyes:misc":
        label = label_parser(label)
        print(label, type(label))

``RichLabel`` objects have an ``Info`` object as an attribute, allowing specific reference to all the specified label fields.

.. jupyter-execute::

    from cogent3.parse.fasta import LabelParser, iter_fasta_records

    fasta_data = [
        ">gi|10047090|ref|NP_055147.1| small muscle protein, X-linked [Homo sapiens]",
        "MNMSKQPVSNVRAIQANINIPMGAFRPGAGQPPRRKECTPEVEEGVPPTSDEEKKPIPGAKKLPGPAVNL",
        "SEIQNIKSELKYVPKAEQ",
        ">gi|10047092|ref|NP_037391.1| neuronal protein [Homo sapiens]",
        "MANRGPSYGLSREVQEKIEQKYDADLENKLVDWIILQCAEDIEHPPPGRAHFQKWLMDGTVLCKLINSLY",
        "PPGQEPIPKISESKMAFKQMEQISQFLKAAETYGVRTTDIFQTVDLWEGKDMAAVQRTLMALGSVAVTKD",
    ]
    label_to_name = LabelParser(
        "%(ref)s",
        [[1, "gi", str], [3, "ref", str], [4, "description", str]],
        split_with="|",
    )
    for name, seq in iter_fasta_records(fasta_data, label_to_name=label_to_name):
        print(name)
        print(name.info.gi)
        print(name.info.description)

.. _storage-plugin:

Using a third-party plugin for sequence storage
-----------------------------------------------

Sequence collections and alignments have a ``.storage`` attribute which holds the underlying sequence data and provides basic functions for obtaining it. Users can install a third-party plugin which is customized for different types of sequence data. The following examples require you install the ``cogent3-h5seqs`` plugin. This project provides alternative storage for both unaligned sequences and for alignments.

.. code-block:: shell

    $ pip install cogent3-h5seqs

Selecting an alternate storage backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specify the storage using the ``storage_backend`` argument.

.. jupyter-execute::
    :raises:

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs(
        "data/long_testseqs.fasta", moltype="dna", storage_backend="h5seqs_aligned"
    )
    aln

That's it!

.. jupyter-execute::

    type(aln.storage)

For the ``cogent3-h5seqs`` package you specify a different storage backend for unaligned sequences.

.. jupyter-execute::
    :raises:

    from cogent3 import load_unaligned_seqs

    seqs = load_unaligned_seqs(
        "data/long_testseqs.fasta", moltype="dna", storage_backend="h5seqs_unaligned"
    )
    type(seqs.storage)

Set the default storage
^^^^^^^^^^^^^^^^^^^^^^^

You can set the default storage process-wide, so you don't need to use the ``storage_backend`` argument.

.. jupyter-execute::
    :raises:

    import cogent3

    cogent3.set_storage_defaults(
        unaligned_seqs="h5seqs_unaligned", aligned_seqs="h5seqs_aligned"
    )

    aln = cogent3.get_dataset("brca1")
    type(aln.storage)

When you apply operations, the new backend storage setting is applied.

.. jupyter-execute::
    :raises:

    coll = aln.degap()
    type(coll.storage)


.. note:: To revert to the ``cogent3`` defaults, use the ``reset`` argument.

    .. jupyter-execute::
        :raises:

        cogent3.set_storage_defaults(reset=True)
