****************************
Manipulating biological data
****************************

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

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

Computing motif probabilities from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For nucleotides.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta', moltype=DNA)
    >>> motif_probs = aln.getMotifProbs()
    >>> print motif_probs
    {'A': 0.24...

For dinucleotides, we need to pass in a dinucleotide alphabet.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> dinuc_alphabet = DNA.Alphabet.getWordAlphabet(2)
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta', moltype=DNA)
    >>> motif_probs = aln.getMotifProbs(alphabet=dinuc_alphabet)
    >>> print motif_probs
    {'AA': 0.078222...

The same holds for codons or other arbitrary alphabets, as long as they match the alignment ``MolType``.

Some calculations in cogent require there to be no-zero's in the motif probabilities, in which case we use a pseudo-count. We illustrate that here with a simple example where T is missing. Without the pseudo-count, the frequency of T is 0.0, with the pseudo-count defined as 1e-6 then the frequency of T will be slightly less than that.

.. doctest::
    
    >>> aln = LoadSeqs(data=[('a', 'AACAAC'),('b', 'AAGAAG')], moltype=DNA)
    >>> motif_probs = aln.getMotifProbs()
    >>> assert motif_probs['T'] == 0.0
    >>> motif_probs = aln.getMotifProbs(pseudocount=1e-6)
    >>> assert 0 < motif_probs['T'] <= 1e-6

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

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. doctest::

    >>> from cogent import DNA
    >>> seq = DNA.makeSequence('ATGATGATGATG')

Creating the position indices, note that we start at the 2nd index (the 'first' codon's 3rd position) indicate each position as a *span* (``i -- i+1``).

.. doctest::

    >>> indices = [(i,i+1) for i in range(len(seq))[2::3]]

Create the sequence feature and use it to slice the sequence.

.. doctest::

    >>> pos3 = seq.addFeature('pos3', 'pos3', indices)
    >>> pos3 = pos3.getSlice()
    >>> assert str(pos3) == 'GGGG'

Getting 1st and 2nd positions from codons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only difference here to above is our span's cover 2 positions.

.. doctest::

    >>> from cogent import DNA
    >>> seq = DNA.makeSequence('ATGATGATGATG')
    >>> indices = [(i,i+2) for i in range(len(seq))[::3]]
    >>> pos12 = seq.addFeature('pos12', 'pos12', indices)
    >>> pos12 = pos12.getSlice()
    >>> assert str(pos12) == 'ATATATAT'

Translation
^^^^^^^^^^^

*To be written.*

RNA
---

*To be written.*

Protein
-------

*To be written.*

Arbitrary
---------

*To be written.*

Alignments
==========

Creating an Alignment object from a SequenceCollection
------------------------------------------------------

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml', aligned=False)
    >>> aln = Alignment(seq)
    >>> fasta_1 = seq.toFasta()
    >>> fasta_2 = aln.toFasta()
    >>> assert fasta_1 == fasta_2

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
    >>> aln = Alignment(seq)
    >>> string_list = aln.todict().values()

Slicing an alignment
--------------------

Getting a single column from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> column_four = aln[3]

Getting a region of columns
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> region = aln[50:70]

Getting codon 3rd positions from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs(data={'seq1': 'ATGATGATG---', 'seq2': 'ATGATGATGATG'})
    >>> indices = [(i,i+1) for i in range(len(aln))[2::3]]
    >>> pos3 = aln.addFeature('pos3', 'pos3', indices)
    >>> pos3 = pos3.getSlice()
    >>> print pos3
    >seq2
    GGGG
    >seq1
    GGG-
    <BLANKLINE>

Filtering positions
-------------------

We sometimes want to eliminate ambiguous or gap data from our alignments. We show how to exclude alignment columns by the characters they contain. In the first instance we do this just for single nucleotide columns, then for trinucleotides (equivalent for handling codons).

.. doctest::
    
    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs(data= [('seq1', 'ATGAAGGTG---'),
    ...                       ('seq2', 'ATGAAGGTGATG'),
    ...                       ('seq3', 'ATGAAGGNGATG')], moltype=DNA)

We now just define a one-line function that returns ``True`` if the passed data contains only nucleotide characters, ``False`` otherwise. The function works by converting the aligned column into a ``set`` and checking it is equal to, or a subset of, all nucleotides. This function, which works for nucleotides or codons, has the effect of eliminating the (nucleotide/trinucleotide) columns with the 'N' and '-' characters.

.. doctest::
    
    >>> just_nucs = lambda x: set(''.join(x)) <= set('ACGT')

We apply to nucleotides,

.. doctest::
    
    >>> nucs = aln.filtered(just_nucs)
    >>> print nucs
    >seq1
    ATGAAGGG
    >seq2
    ATGAAGGG
    >seq3
    ATGAAGGG
    <BLANKLINE>

and trinucleotides (specified by setting ``motif_length=3``).

.. doctest::
    
    >>> trinucs = aln.filtered(just_nucs, motif_length=3)
    >>> print trinucs
    >seq1
    ATGAAG
    >seq2
    ATGAAG
    >seq3
    ATGAAG
    <BLANKLINE>

Filtering sequences
-------------------

*To be written.*

Translating
-----------

*To be written.*

Trees
=====

Selecting subtrees
------------------

*To be written.*

Drawing trees
-------------

*To be written.*

.. pdf, asciiArt

Tabular data
============

.. doctest::
    :hide:

    >>> # just saving some tabular data for subsequent data
    >>> from cogent import LoadTable
    >>> rows = (('NP_003077', 'Con', 2.5386013224378985),
    ... ('NP_004893', 'Con', 0.12135142635634111e+06),
    ... ('NP_005079', 'Con', 0.95165949788861326e+07),
    ... ('NP_005500', 'NonCon', 0.73827030202664901e-07),
    ... ('NP_055852', 'NonCon', 1.0933217708952725e+07))
    >>> table = LoadTable(header=['Locus', 'Region', 'Ratio'], rows=rows)
    >>> table.writeToFile('stats.txt', sep=',')

Loading delimited formats
-------------------------

We load a comma separated data file using the generic ``LoadTable`` function.

.. doctest::

    >>> from cogent import LoadTable
    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

Reading large files
-------------------

For really large files the automated conversion used by the standard read mechanism can be quite slow. If the data within a column is consistently of one type, set the ``LoadTable`` argument ``static_column_types=True``. This causes the ``Table`` object to create a custom reader.

.. doctest::

    >>> table = LoadTable('stats.txt', static_column_types=True)
    >>> print table
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

Formatting
----------

Changing displayed numerical precision
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We change the ``Ratio`` column to using scientific notation.

.. doctest::

    >>> table.setColumnFormat('Ratio', '%.1e')
    >>> print table
    ==============================
        Locus    Region      Ratio
    ------------------------------
    NP_003077       Con    2.5e+00
    NP_004893       Con    1.2e+05
    NP_005079       Con    9.5e+06
    NP_005500    NonCon    7.4e-08
    NP_055852    NonCon    1.1e+07
    ------------------------------

Change digits or column spacing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be done on table loading,

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',', digits=1, space=2)
    >>> print table
    =============================
        Locus  Region       Ratio
    -----------------------------
    NP_003077     Con         2.5
    NP_004893     Con    121351.4
    NP_005079     Con   9516595.0
    NP_005500  NonCon         0.0
    NP_055852  NonCon  10933217.7
    -----------------------------

or, for spacing at least, by modifying the attributes

.. doctest::

    >>> table.Space = '    '
    >>> print table
    =================================
        Locus    Region         Ratio
    ---------------------------------
    NP_003077       Con           2.5
    NP_004893       Con      121351.4
    NP_005079       Con     9516595.0
    NP_005500    NonCon           0.0
    NP_055852    NonCon    10933217.7
    ---------------------------------

Changing column headings
------------------------

The table ``Header`` is immutable. Changing column headings is done as follows.

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table.Header
    ['Locus', 'Region', 'Ratio']
    >>> table = table.withNewHeader('Ratio', 'Stat')
    >>> print table.Header
    ['Locus', 'Region', 'Stat']

Creating new columns from existing ones
---------------------------------------

This can be used to take a single, or multiple columns and generate a new column of values. Here we'll take 2 columns and return True/False based on a condition.

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> table = table.withNewColumn('LargeCon',
    ...                     lambda (r,v): r == 'Con' and v>10.0,
    ...                     columns=['Region', 'Ratio'])
    >>> print table
    ================================================
        Locus    Region            Ratio    LargeCon
    ------------------------------------------------
    NP_003077       Con           2.5386       False
    NP_004893       Con      121351.4264        True
    NP_005079       Con     9516594.9789        True
    NP_005500    NonCon           0.0000       False
    NP_055852    NonCon    10933217.7090       False
    ------------------------------------------------

Appending tables
----------------

Can be done without specifying a new column. Here we simply use the same table data.

.. doctest::

    >>> table1 = LoadTable('stats.txt', sep=',')
    >>> table2 = LoadTable('stats.txt', sep=',')
    >>> table = table1.appended(None, table2)
    >>> print table
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_005500    NonCon           0.0000
    NP_055852    NonCon    10933217.7090
    ------------------------------------

or with a new column

.. doctest::

    >>> table1.Title = 'Data1'
    >>> table2.Title = 'Data2'
    >>> table = table1.appended('Data#', table2, title='')
    >>> print table
    =============================================
    Data#        Locus    Region            Ratio
    ---------------------------------------------
    Data1    NP_003077       Con           2.5386
    Data1    NP_004893       Con      121351.4264
    Data1    NP_005079       Con     9516594.9789
    Data1    NP_005500    NonCon           0.0000
    Data1    NP_055852    NonCon    10933217.7090
    Data2    NP_003077       Con           2.5386
    Data2    NP_004893       Con      121351.4264
    Data2    NP_005079       Con     9516594.9789
    Data2    NP_005500    NonCon           0.0000
    Data2    NP_055852    NonCon    10933217.7090
    ---------------------------------------------

.. note:: We assigned an empty string to ``title``, otherwise the resulting table has the same ``Title`` attribute as that of ``table1``.

Summing a single column
-----------------------

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> table.summed('Ratio')
    20571166.652847398

Summing multiple columns or rows - strictly numerical data
----------------------------------------------------------

We define a strictly numerical table,

.. doctest::

    >>> all_numeric = LoadTable(header=['A', 'B', 'C'], rows=[range(3),
    ...                                 range(3,6), range(6,9), range(9,12)])
    >>> print all_numeric
    =============
    A     B     C
    -------------
    0     1     2
    3     4     5
    6     7     8
    9    10    11
    -------------

and sum all columns (default condition)

.. doctest::

    >>> all_numeric.summed()
    [18, 22, 26]

and all rows

.. doctest::

    >>> all_numeric.summed(col_sum=False)
    [3, 12, 21, 30]

Summing multiple columns or rows with mixed non-numeric/numeric data
--------------------------------------------------------------------

We define a table with mixed data, like a distance matrix.

.. doctest::

    >>> mixed = LoadTable(header=['A', 'B', 'C'], rows=[['*',1,2], [3,'*', 5],
    ...                                                 [6,7,'*']])
    >>> print mixed
    ===========
    A    B    C
    -----------
    *    1    2
    3    *    5
    6    7    *
    -----------

and sum all columns (default condition), ignoring non-numerical data

.. doctest::

    >>> mixed.summed(strict=False)
    [9, 8, 7]

and all rows

.. doctest::

    >>> mixed.summed(col_sum=False, strict=False)
    [3, 8, 13]


Filtering table rows
--------------------

We can do this by providing a reference to an external function

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> sub_table = table.filtered(lambda x: x < 10.0, columns='Ratio')
    >>> print sub_table
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_003077       Con    2.5386
    NP_005500    NonCon    0.0000
    -----------------------------

or using valid python syntax within a string, which is executed

.. doctest::

    >>> sub_table = table.filtered("Ratio < 10.0")
    >>> print sub_table
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_003077       Con    2.5386
    NP_005500    NonCon    0.0000
    -----------------------------

You can also filter for values in multiple columns

.. doctest::

    >>> sub_table = table.filtered("Ratio < 10.0 and Region == 'NonCon'")
    >>> print sub_table
    =============================
        Locus    Region     Ratio
    -----------------------------
    NP_005500    NonCon    0.0000
    -----------------------------

Filtering table columns
-----------------------

We select only columns that have a sum > 20 from the ``all_numeric`` table constructed above.

.. doctest::
    
    >>> big_numeric = all_numeric.filteredByColumn(lambda x: sum(x)>20)
    >>> print big_numeric
    ========
     B     C
    --------
     1     2
     4     5
     7     8
    10    11
    --------

Sorting
-------

Standard sorting
^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table.sorted(columns='Ratio')
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_005500    NonCon           0.0000
    NP_003077       Con           2.5386
    NP_004893       Con      121351.4264
    NP_005079       Con     9516594.9789
    NP_055852    NonCon    10933217.7090
    ------------------------------------

Reverse sorting
^^^^^^^^^^^^^^^

.. doctest::

    >>> print table.sorted(columns='Ratio', reverse='Ratio')
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_055852    NonCon    10933217.7090
    NP_005079       Con     9516594.9789
    NP_004893       Con      121351.4264
    NP_003077       Con           2.5386
    NP_005500    NonCon           0.0000
    ------------------------------------

Sorting involving multiple columns, one reversed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> print table.sorted(columns=['Region', 'Ratio'], reverse='Ratio')
    ====================================
        Locus    Region            Ratio
    ------------------------------------
    NP_005079       Con     9516594.9789
    NP_004893       Con      121351.4264
    NP_003077       Con           2.5386
    NP_055852    NonCon    10933217.7090
    NP_005500    NonCon           0.0000
    ------------------------------------

Getting raw data
----------------

For a single column
^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> raw = table.getRawData('Region')
    >>> print raw
    ['Con', 'Con', 'Con', 'NonCon', 'NonCon']

For multiple columns
^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> raw = table.getRawData(['Locus', 'Region'])
    >>> print raw
    [['NP_003077', 'Con'], ['NP_004893', 'Con'], ...

Iterating over table rows
-------------------------

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> for row in table:
    ...     print row['Locus']
    ...
    NP_003077
    NP_004893
    NP_005079
    NP_005500
    NP_055852

Table slicing
-------------

Using column names
^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table[:2, :'Region']
    =========
        Locus
    ---------
    NP_003077
    NP_004893
    ---------

Using column indices
^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table[:2,: 1]
    =========
        Locus
    ---------
    NP_003077
    NP_004893
    ---------

SQL like capabilities
---------------------

Distinct values
^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> assert table.getDistinctValues('Region') == set(['NonCon', 'Con'])

Counting
^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> assert table.count("Region == 'NonCon' and Ratio > 1") == 1

Joining tables
^^^^^^^^^^^^^^

SQL like join operations requires tables have different ``Title`` attributes which are not ``None``. We do a standard inner join here for a restricted subset. We must specify the columns that will be used for the join. Here we just use ``Locus`` but multiple columns can be used, and their names can be different between the tables. Note that the second table's title becomes a part of the column names.

.. doctest::

    >>> rows = [['NP_004893', True], ['NP_005079', True],
    ...         ['NP_005500', False], ['NP_055852', False]]
    >>> region_type = LoadTable(header=['Locus', 'LargeCon'], rows=rows,
    ...                 title='RegionClass')
    >>> stats_table = LoadTable('stats.txt', sep=',', title='Stats')
    >>> new = stats_table.joined(region_type, columns_self='Locus')
    >>> print new
    ============================================================
        Locus    Region            Ratio    RegionClass_LargeCon
    ------------------------------------------------------------
    NP_004893       Con      121351.4264                    True
    NP_005079       Con     9516594.9789                    True
    NP_005500    NonCon           0.0000                   False
    NP_055852    NonCon    10933217.7090                   False
    ------------------------------------------------------------

Exporting
---------

Writing delimited formats
^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> table.writeToFile('stats_tab.txt', sep='\t')

Writing latex format
^^^^^^^^^^^^^^^^^^^^

It is also possible to specify column alignment, table caption and other arguments.

.. doctest::

    >>> table = LoadTable('stats.txt', sep=',')
    >>> print table.tostring(format='latex')
    \begin{longtable}[htp!]{ r r r }
    \hline
    \bf{Locus} & \bf{Region} & \bf{Ratio} \\
    \hline
    \hline
    NP_003077 &    Con &        2.5386 \\
    NP_004893 &    Con &   121351.4264 \\
    NP_005079 &    Con &  9516594.9789 \\
    NP_005500 & NonCon &        0.0000 \\
    NP_055852 & NonCon & 10933217.7090 \\
    \hline
    \end{longtable}

.. we remove the table data

.. doctest::
    :hide:

    >>> import os
    >>> os.remove('stats.txt')
    >>> os.remove('stats_tab.txt')

Structure data
==============

2D
--

*To be written.*

3D
--

*To be written.*

Visualisation
-------------

*To be written.*

Sequence metadata
=================

Annotations with coordinates
----------------------------

*To be written.*

Automated introduction from reading genbank files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*To be written.*

Manipulating annotated regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*To be written.*

Annotation display
^^^^^^^^^^^^^^^^^^

*To be written.*

Introducing your own
^^^^^^^^^^^^^^^^^^^^

*To be written.*

Displaying them
^^^^^^^^^^^^^^^

*To be written.*

Generic metadata
----------------

*To be written.*

Info object
^^^^^^^^^^^

*To be written.*

