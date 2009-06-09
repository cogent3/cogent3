****************************
Manipulating biological data
****************************

.. sectionauthor:: Gavin Huttley, Kristian Rother, Patrick Yannul

Genetic code
============

Translate DNA sequences
-----------------------

.. doctest::
    
    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> standard_code.translate('TTTGCAAAC')
    'FAN'

Translate a codon
-----------------

.. doctest::
    
    >>> from cogent.core.genetic_code import DEFAULT as standard_code
    >>> standard_code['TTT']
    'F'

or get the codons for a single amino acid

.. doctest::
    
    >>> standard_code['A']
    ['GCT', 'GCC', 'GCA', 'GCG']

Sequences
=========

DNA
---

Create a DNA sequence object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> my_seq
    DnaSequence(AGTACAC... 11)
    >>> print my_seq
    AGTACACTGGT
    >>> str(my_seq)
    'AGTACACTGGT'

Inverting a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> print my_seq.complement()
    TCATGTGACCA
    >>> print my_seq.reversecomplement()
    ACCAGTGTACT

or reverse complement using the shorter ``rc`` method.

.. doctest::

    >>> print my_seq.rc()
    ACCAGTGTACT

Transcribing a DNA sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> print my_seq.toRna()
    AGUACACUGGU

Comparing complementarity:
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import DNA
    >>> a = DNA.makeSequence("AGTACACTGGT")
    >>> a.canPair(a.complement())
    False
    >>> a.canPair(a.reversecomplement())
    True

Joining two DNA sequences
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> my_seq = DNA.makeSequence("AGTACACTGGT")
    >>> a = DNA.makeSequence("AGTACACTGGT")
    >>> print my_seq + a
    AGTACACTGGTAGTACACTGGT

Slicing DNA sequences
^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> my_seq[1:6]
    DnaSequence(GTACA)

Converting to Fasta format
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> print my_seq.toFasta()
    >0
    AGTACACTGGT

Changing the Name of a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> my_seq.Name = 'my_gene'
    >>> print my_seq.toFasta()
    >my_gene
    AGTACACTGGT

Convert to codon series
^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> my_seq.CodonAlphabet()
    ('CTT', 'ACC', 'ACA', 'ACG', 'ATC', 'ATA',...

Creating a sequence with a name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import PROTEIN
    >>> p = PROTEIN.makeSequence('THISISAPRQTEIN','myProtein')

Creating a general sequence object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.core.sequence import Sequence
    >>> seq = Sequence('ACDEF','Name')
    >>> print seq.toFasta()
    >Name
    ACDEF

Parsing files with many sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reading a FASTA file with DNA sequences
"""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs
    >>> seq = LoadSeqs('data/test.paml', aligned=True)
    >>> print seq
    >NineBande
    GCAAGGCGCCAACAGAGCAGATGGGCTGAAAGTAAGGAAACATGTAATGATAGGCAGACT
    >Mouse
    GCAGTGAGCCAGCAGAGCAGATGGGCTGCAAGTAAAGGAACATGTAACGACAGGCAGGTT
    >Human
    GCAAGGAGCCAACATAACAGATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACT
    >HowlerMon
    GCAAGGAGCCAACATAACAGATGGGCTGAAAGTGAGGAAACATGTAATGATAGGCAGACT
    >DogFaced
    GCAAGGAGCCAGCAGAACAGATGGGTTGAAACTAAGGAAACATGTAATGATAGGCAGACT
    <BLANKLINE>

Loading protein sequences in a Phylip file
""""""""""""""""""""""""""""""""""""""""""
.. doctest::

    >>> seq = LoadSeqs('data/abglobin_aa.phylip', moltype=PROTEIN,
    ...              aligned=True)

RNA with modifications from a dict
""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs, RNA
    >>> rna = {'seq1': '--ACGU--GU---',
    ...        'seq2': '--ACGUA-GU---',
    ...        'seq3': '--ACGUA-GU---'}
    >>> seqs = LoadSeqs(data=rna, moltype=RNA)

Loading FASTA sequences from an open file or list of lines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent.parse.fasta import FastaParser
    >>> f=open('data/long_testseqs.fasta')
    >>> seqs = [(name, seq) for name, seq in FastaParser(f)]
    >>> print seqs
    [('Human', ByteSequence(TGTGGCA... 2532)), ('HowlerMon',...

Loading DNA sequences from a GenBank file
"""""""""""""""""""""""""""""""""""""""""

.. todo:: get sample data for this

Converting a SequenceCollection to FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs
    >>> seq = LoadSeqs('data/test.paml', aligned=False)
    >>> fasta_data = seq.toFasta()

Removing some sequences from the alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/test.paml')
    >>> aln.Names
    ['NineBande', 'Mouse', 'Human', 'HowlerMon', 'DogFaced']
    >>> new = aln.takeSeqs(['Human', 'HowlerMon'])

Calculating gap fractions for each column in an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta')
    >>> for column in aln[113:150].iterPositions():
    ...     ungapped = filter(lambda x:x=='-', column)
    ...     gap_fraction = len(ungapped)*1.0/len(column)
    ...     print gap_fraction
    0.0
    0.666666666667
    0.0
    0.0...

Getting all variable positions from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> just_variable_aln = aln.filtered(lambda x: len(set(x)) > 1)
    >>> print just_variable_aln[:10]
    >Human
    AAGCAAAACT
    >HowlerMon
    AAGCAAGACT
    >Mouse
    GGGCCCAGCT
    >NineBande
    AAATAAAACT
    >DogFaced
    AAACAAAATA
    <BLANKLINE>

Getting all variable codons from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> variable_codons = aln.filtered(lambda x: len(set(x)) > 1,
    ...                                  motif_length=3)
    >>> print just_variable_aln[:9]
    >Human
    AAGCAAAAC
    >HowlerMon
    AAGCAAGAC
    >Mouse
    GGGCCCAGC
    >NineBande
    AAATAAAAC
    >DogFaced
    AAACAAAAT
    <BLANKLINE>

Remove all gaps from an alignment in FASTA format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs("data/primate_cdx2_promoter.fasta")
    >>> degapped = aln.degap()

Getting the third sequence from an Alignment as a Sequence object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/test.paml')
    >>> seq = aln.getSeq(aln.Names[2])


Getting 3rd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Getting 1st and 2nd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Translation
^^^^^^^^^^^

RNA
---

Protein
-------

Arbitrary
---------

Alignments
==========

Creating an Alignment object from a SequenceCollection
------------------------------------------------------

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml', aligned=False)
    >>> ali = Alignment(seq)
    >>> fasta_1 = seq.toFasta()
    >>> fasta_2 = ali.toFasta()
    >>> fasta_1 == fasta_2
    True

Converting an alignment to FASTA format
---------------------------------------

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/long_testseqs.fasta')
    >>> aln = Alignment(seq)
    >>> fasta_align = aln.toFasta()

Converting an alignment into Phylip format
------------------------------------------

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> phylip_file, name_dictionary = aln.toPhylip()

Convert an alignment to a list of strings
-----------------------------------------

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml')
    >>> ali = Alignment(seq)
    >>> string_list = ali.todict().values()

Slicing an alignment
--------------------

Getting a single column from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml')
    >>> ali = Alignment(seq)
    >>> column_four = ali[3]

Getting a region of columns
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> region = aln[50:70]

Filtering positions
-------------------

Filtering sequences
-------------------

Translating
-----------

Trees
=====

Selecting subtrees
------------------

Drawing trees
-------------

.. pdf, asciiArt

Tabular data
============

SQL like capabilities
---------------------

Reading large files
-------------------

Formatting
----------

.. columns for display, digits, spaces

Getting raw data
----------------

Filtering results
-----------------

Sorting
-------

Exporting
---------

Structure data
==============

2D
--

3D
--

Visualisation
-------------

Sequence metadata
=================

Annotations with coordinates
----------------------------

Automated introduction from reading genbank files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Manipulating annotated regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotation display
^^^^^^^^^^^^^^^^^^

Introducing your own
^^^^^^^^^^^^^^^^^^^^

Displaying them
^^^^^^^^^^^^^^^

Generic metadata
----------------

Info object
^^^^^^^^^^^

