.. _load-seqs:

Loading an alignment from a file
--------------------------------

.. author, Gavin Huttley, Tony Walters, Tom Elliott

Loading aligned sequences
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/long_testseqs.fasta', moltype="dna")
    >>> type(aln)
    <class 'cogent3.core.alignment.ArrayAlignment'>

The load functions record the origin of the data in the ``info`` attribute under a `"source"` key.

.. doctest::
    
    >>> aln.info.source
    'data/long_testseqs.fasta'

.. note:: The function ``load_aligned_seqs()`` returns an ``ArrayAlignment`` by default. If you set the argument ``array_align=False``, you will get an ``Alignment``. (That class can be annotated.)

.. todo:: add cross ref for description of Info class

Loading unaligned sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``load_unaligned_seqs()`` function returns a sequence collection.

.. doctest::

    >>> from cogent3 import load_unaligned_seqs
    >>> seqs = load_unaligned_seqs('data/long_testseqs.fasta', moltype="dna")
    >>> type(seqs)
    <class 'cogent3.core.alignment.SequenceCollection'>

Specifying the file format
^^^^^^^^^^^^^^^^^^^^^^^^^^

The loading functions use the filename suffix to infer the file format. This can be overridden using the ``format`` argument.

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/long_testseqs.fasta', moltype="dna",
    ...                  format='fasta')
    ...
    >>> aln
    5 x 2532 dna alignment: Human[TGTGGCACAAA...

Specifying the sequence molecular type
--------------------------------------

Simple case of loading a ``list`` of aligned amino acid sequences in FASTA format, with and without ``moltype`` specification. When ``moltype`` is not specified it defaults to ``BYTES`` for the ``ArrayAlignment`` class, ``ASCII`` for the ``Alignment`` class.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> protein_seqs = ['>seq1','DEKQL-RG','>seq2','DDK--SRG']
    >>> proteins_loaded = make_aligned_seqs(protein_seqs)
    >>> proteins_loaded.moltype
    MolType(('\x00', '\x01', '\x02', '\x03'...
    >>> print(proteins_loaded)
    >seq1
    DEKQL-RG
    >seq2
    DDK--SRG
    <BLANKLINE>
    >>> proteins_loaded = make_aligned_seqs(protein_seqs, moltype="protein")
    >>> print(proteins_loaded)
    >seq1
    DEKQL-RG
    >seq2
    DDK--SRG
    <BLANKLINE>

.. note:: This applies to both the ``load_*`` or ``make_*`` functions.

Making an alignment from standard python objects
------------------------------------------------

From a series of strings
^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> seqs = ['>seq1','AATCG-A','>seq2','AATCGGA']
    >>> seqs_loaded = make_aligned_seqs(seqs)
    >>> print(seqs_loaded)
    >seq1
    AATCG-A
    >seq2
    AATCGGA
    <BLANKLINE>

From a dict of strings
^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> seqs = {'seq1': 'AATCG-A','seq2': 'AATCGGA'}
    >>> seqs_loaded = make_aligned_seqs(seqs)

Stripping label characters on loading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Load a list of aligned nucleotide sequences, while specifying the DNA molecule type and stripping the comments from the label. In this example, stripping is accomplished by passing a function that removes everything after the first whitespace to the ``label_to_name`` parameter.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> DNA_seqs = ['>sample1 Mus musculus','AACCTGC--C','>sample2 Gallus gallus','AAC-TGCAAC']
    >>> loaded_seqs = make_aligned_seqs(DNA_seqs, moltype="dna", label_to_name=lambda x: x.split()[0])
    >>> print(loaded_seqs)
    >sample1
    AACCTGC--C
    >sample2
    AAC-TGCAAC
    <BLANKLINE>

Loading sequences using format parsers
--------------------------------------

``load_aligned_seqs()`` and ``load_unaligned_seqs()`` are just convenience interfaces to format parsers. It can sometimes be more effective to use the parsers directly, say when you don't want to load everything into memory.

Loading FASTA sequences from an open file or list of lines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load FASTA formatted sequences directly, you can use the ``MinimalFastaParser``.

.. note:: This returns the sequences as strings.

.. doctest::

    >>> from cogent3.parse.fasta import MinimalFastaParser
    >>> f=open('data/long_testseqs.fasta')
    >>> seqs = [(name, seq) for name, seq in MinimalFastaParser(f)]
    >>> seqs
    [('Human', 'TGTGGCACAAATAC...

Handling overloaded FASTA sequence labels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The FASTA label field is frequently overloaded, with different information fields present in the field and separated by some delimiter. This can be flexibly addressed using the ``LabelParser``. By creating a custom label parser, we can decide which part we use as the sequence name. We show how to convert a field into something specific.

.. doctest::

    >>> from cogent3.parse.fasta import LabelParser
    >>> def latin_to_common(latin):
    ...     return {'Homo sapiens': 'human',
    ...             'Pan troglodtyes': 'chimp'}[latin]
    >>> label_parser = LabelParser("%(species)s",
    ...             [[1, "species", latin_to_common]], split_with=':')
    >>> for label in ">abcd:Homo sapiens:misc", ">abcd:Pan troglodtyes:misc":
    ...     label = label_parser(label)
    ...     print(label, type(label))
    human <class 'cogent3.parse.fasta.RichLabel'>
    chimp <class 'cogent3.parse.fasta.RichLabel'>

``RichLabel`` objects have an ``Info`` object as an attribute, allowing specific reference to all the specified label fields.

.. doctest::

    >>> from cogent3.parse.fasta import MinimalFastaParser, LabelParser
    >>> fasta_data = ['>gi|10047090|ref|NP_055147.1| small muscle protein, X-linked [Homo sapiens]',
    ...  'MNMSKQPVSNVRAIQANINIPMGAFRPGAGQPPRRKECTPEVEEGVPPTSDEEKKPIPGAKKLPGPAVNL',
    ... 'SEIQNIKSELKYVPKAEQ',
    ... '>gi|10047092|ref|NP_037391.1| neuronal protein [Homo sapiens]',
    ... 'MANRGPSYGLSREVQEKIEQKYDADLENKLVDWIILQCAEDIEHPPPGRAHFQKWLMDGTVLCKLINSLY',
    ... 'PPGQEPIPKISESKMAFKQMEQISQFLKAAETYGVRTTDIFQTVDLWEGKDMAAVQRTLMALGSVAVTKD']
    ...
    >>> label_to_name = LabelParser("%(ref)s",
    ...                              [[1,"gi", str],
    ...                               [3, "ref", str],
    ...                               [4, "description", str]],
    ...                               split_with="|")
    ...
    >>> for name, seq in MinimalFastaParser(fasta_data, label_to_name=label_to_name):
    ...     print(name)
    ...     print(name.info.gi)
    ...     print(name.info.description)
    NP_055147.1
    10047090
     small muscle protein, X-linked [Homo sapiens]
    NP_037391.1
    10047092
     neuronal protein [Homo sapiens]

