.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _load_seq:

Loading a sequence from a file
------------------------------

It's also possible to load a sequence from a :ref:`url <load_url>`.

.. jupyter-execute::

    from cogent3 import load_seq
    
    seq = load_seq("data/mycoplasma-genitalium.fa", moltype="dna")
    seq

.. warning:: If a file has more than one sequence, only the first one is loaded.

.. jupyter-execute::

    seq = load_seq("data/brca1-bats.fasta", moltype="dna")
    seq

.. note:: The filename suffix is used to infer the data format.

.. _load-seqs:

Loading an alignment from a file or url
---------------------------------------

.. author, Gavin Huttley, Tony Walters, Tom Elliott

Loading aligned sequences
^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    type(aln)

The load functions record the origin of the data in the ``info`` attribute under a `"source"` key.

.. jupyter-execute::

    aln.info.source

.. note:: The function ``load_aligned_seqs()`` returns an ``ArrayAlignment`` by default. If you set the argument ``array_align=False``, you will get an ``Alignment``. (That class can be annotated.)

.. todo:: add cross ref for description of Info class

Loading unaligned sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``load_unaligned_seqs()`` function returns a sequence collection.

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


    aln = load_aligned_seqs("https://raw.githubusercontent.com/cogent3/cogent3/develop/doc/data/long_testseqs.fasta", moltype="dna")

Specifying the file format
^^^^^^^^^^^^^^^^^^^^^^^^^^

The loading functions use the filename suffix to infer the file format. This can be overridden using the ``format`` argument.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna", format="fasta")
    aln

Specifying the sequence molecular type
--------------------------------------

Simple case of loading a ``list`` of aligned amino acid sequences in FASTA format, with and without ``moltype`` specification. When ``moltype`` is not specified it defaults to ``BYTES`` for the ``ArrayAlignment`` class, ``ASCII`` for the ``Alignment`` class.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    protein_seqs = [">seq1", "DEKQL-RG", ">seq2", "DDK--SRG"]
    proteins_loaded = make_aligned_seqs(protein_seqs)
    proteins_loaded.moltype
    print(proteins_loaded)
    proteins_loaded = make_aligned_seqs(protein_seqs, moltype="protein")
    print(proteins_loaded)

.. note:: This applies to both the ``load_*`` or ``make_*`` functions.

Making an alignment from standard python objects
------------------------------------------------

From a dict of strings
^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    seqs = {"seq1": "AATCG-A", "seq2": "AATCGGA"}
    seqs_loaded = make_aligned_seqs(seqs)

From a series of strings
^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    seqs = [">seq1", "AATCG-A", ">seq2", "AATCGGA"]
    seqs_loaded = make_aligned_seqs(seqs)
    print(seqs_loaded)

Stripping label characters on loading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Load a list of aligned nucleotide sequences, while specifying the DNA molecule type and stripping the comments from the label. In this example, stripping is accomplished by passing a function that removes everything after the first whitespace to the ``label_to_name`` parameter.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    DNA_seqs = [
        ">sample1 Mus musculus",
        "AACCTGC--C",
        ">sample2 Gallus gallus",
        "AAC-TGCAAC",
    ]
    loaded_seqs = make_aligned_seqs(
        DNA_seqs, moltype="dna", label_to_name=lambda x: x.split()[0]
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

To load FASTA formatted sequences directly, you can use the ``MinimalFastaParser``.

.. note:: This returns the sequences as strings.

.. jupyter-execute::

    from cogent3.parse.fasta import MinimalFastaParser

    f = open("data/long_testseqs.fasta")
    seqs = [(name, seq) for name, seq in MinimalFastaParser(f)]
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

    from cogent3.parse.fasta import LabelParser, MinimalFastaParser

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
    for name, seq in MinimalFastaParser(fasta_data, label_to_name=label_to_name):
        print(name)
        print(name.info.gi)
        print(name.info.description)
