DNA and RNA sequences
---------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott

DNA sequence objects
^^^^^^^^^^^^^^^^^^^^

Creating a DNA sequence from a string
"""""""""""""""""""""""""""""""""""""

All sequence and alignment objects have a molecular type, or ``MolType`` which provides key properties for validating sequence characters. Here we use the ``DNA`` ``MolType`` to create a DNA sequence.

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> my_seq
    DnaSequence(AGTACAC... 11)
    >>> print my_seq
    AGTACACTGGT
    >>> str(my_seq)
    'AGTACACTGGT'

Converting to FASTA format
""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('AGTACACTGGT')
    >>> print my_seq.toFasta()
    >0
    AGTACACTGGT

Creating a named sequence
"""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('AGTACACTGGT','my_gene')
    >>> my_seq
    DnaSequence(AGTACAC... 11)
    >>> type(my_seq)
    <class 'cogent.core.sequence.DnaSequence'>

Setting or changing the name of a sequence
""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence('AGTACACTGGT')
    >>> my_seq.Name = 'my_gene'
    >>> print my_seq.toFasta()
    >my_gene
    AGTACACTGGT

Inverting a DNA sequence
""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> print my_seq.complement()
    TCATGTGACCA
    >>> print my_seq.reversecomplement()
    ACCAGTGTACT

The ``rc`` method name is easier to type

.. doctest::

    >>> print my_seq.rc()
    ACCAGTGTACT

.. _translation:

Converting a DnaSequence object to protein
""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT",'X')
    >>> pep = my_seq.getTranslation()
    >>> type(pep)
    <class 'cogent.core.sequence.ProteinSequence'>
    >>> print pep.toFasta()
    >X
    STL

Converting a DNA sequence to RNA
""""""""""""""""""""""""""""""""

.. doctest::

    >>> print my_seq.toRna()
    AGUACACUGGU

Testing complementarity
"""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> a = DNA.makeSequence("AGTACACTGGT")
    >>> a.canPair(a.complement())
    False
    >>> a.canPair(a.reversecomplement())
    True

Joining two DNA sequences
"""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> extra_seq = DNA.makeSequence("CTGAC")
    >>> long_seq = my_seq + extra_seq
    >>> long_seq
    DnaSequence(AGTACAC... 16)
    >>> str(long_seq)
    'AGTACACTGGTCTGAC'

Slicing DNA sequences
"""""""""""""""""""""

.. doctest::

    >>> my_seq[1:6]
    DnaSequence(GTACA)

Getting 3rd positions from codons
"""""""""""""""""""""""""""""""""

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. doctest::

    >>> from cogent import DNA
    >>> seq = DNA.makeSequence('ATGATGATGATG')

Creating the position indices, note that we start at the 2nd index (the 'first' codon's 3rd position) indicate each position as a *span* (``i -- i+1``).

.. doctest::

    >>> indices = [(i, i+1) for i in range(len(seq))[2::3]]

Create the sequence feature and use it to slice the sequence.

.. doctest::

    >>> pos3 = seq.addFeature('pos3', 'pos3', indices)
    >>> pos3 = pos3.getSlice()
    >>> assert str(pos3) == 'GGGG'

Getting 1st and 2nd positions from codons
"""""""""""""""""""""""""""""""""""""""""

The only difference here to above is that our spans cover 2 positions.

.. doctest::

    >>> from cogent import DNA
    >>> seq = DNA.makeSequence('ATGATGATGATG')
    >>> indices = [(i, i+2) for i in range(len(seq))[::3]]
    >>> pos12 = seq.addFeature('pos12', 'pos12', indices)
    >>> pos12 = pos12.getSlice()
    >>> assert str(pos12) == 'ATATATAT'

Creating a general sequence object
""""""""""""""""""""""""""""""""""

.. note:: the import statement

.. doctest::

    >>> from cogent.core.sequence import Sequence
    >>> seq = Sequence('ABCDEF','Name')
    >>> print seq.toFasta()
    >Name
    ABCDEF
    >>> print type(seq)
    <class 'cogent.core.sequence.Sequence'>

Loading sequences from a file
"""""""""""""""""""""""""""""

For loading collections of unaligned or aligned sequences see :ref:`load-seqs`.

Loading FASTA sequences from an open file or list of lines
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To load sequences from a fasta file directly, you can use the ``MinimalFastaParser``.

.. note:: This returns the sequences as strings.

.. doctest::

    >>> from cogent.parse.fasta import MinimalFastaParser
    >>> f=open('data/long_testseqs.fasta')
    >>> seqs = [(name, seq) for name, seq in MinimalFastaParser(f)]
    >>> print seqs
    [('Human', 'TGTGGCACAAATAC...

Handling overloaded FASTA sequence labels
+++++++++++++++++++++++++++++++++++++++++

The FASTA label field is frequently overloaded, with different information fields present in the field and separated by some delimiter. This can be flexibly addressed using the ``LabelParser``. By creating a custom label parser, we can decided which part we use as the sequence name. We show how convert a field into something specific.

.. doctest::
    
    >>> from cogent.parse.fasta import LabelParser
    >>> def latin_to_common(latin):
    ...     return {'Homo sapiens': 'human',
    ...             'Pan troglodtyes': 'chimp'}[latin]
    >>> label_parser = LabelParser("%(species)s",
    ...             [[1, "species", latin_to_common]], split_with=':')
    >>> for label in ">abcd:Homo sapiens:misc", ">abcd:Pan troglodtyes:misc":
    ...     label = label_parser(label)
    ...     print label, type(label)
    human <class 'cogent.parse.fasta.RichLabel'>
    chimp <class 'cogent.parse.fasta.RichLabel'>

The ``RichLabel`` objects have an ``Info`` object as an attribute, allowing specific reference to all the specified label fields.

.. doctest::
    
    >>> from cogent.parse.fasta import MinimalFastaParser, LabelParser
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
    ...     print name
    ...     print name.Info.gi
    ...     print name.Info.description
    NP_055147.1
    10047090
     small muscle protein, X-linked [Homo sapiens]
    NP_037391.1
    10047092
     neuronal protein [Homo sapiens]

Loading DNA sequences from a GenBank file
+++++++++++++++++++++++++++++++++++++++++

.. todo:: get sample data for this

*To be written.*
