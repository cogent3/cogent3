Collections and Alignments
--------------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott

Loading sequences from a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As an alignment
"""""""""""""""

The function ``LoadSeqs()`` creates either a sequence collection or an alignment depending on the keyword argument ``aligned`` (the default is ``True``).

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/long_testseqs.fasta', moltype=DNA)
    >>> type(aln)
    <class 'cogent.core.alignment.Alignment'>

This example and those following used the file named `long_testseqs.fasta` available in the `data` directory. You can find it here :download:`long_testseqs.fasta <../data/long_testseqs.fasta>`.

As a sequence collection (unaligned)
""""""""""""""""""""""""""""""""""""

Setting the ``LoadSeqs()`` function keyword argument ``aligned=False`` returns a sequence collection.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> seqs = LoadSeqs('data/long_testseqs.fasta', moltype=DNA, aligned=False)
    >>> print type(seqs)
    <class 'cogent.core.alignment.SequenceCollection'>

.. note:: An alignment can be sliced, but a ``SequenceCollection`` can not.

Specifying the file format
""""""""""""""""""""""""""

``LoadSeqs`` uses the filename suffix to infer the file format. This can be overridden using the format argument.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/long_testseqs.fasta', moltype=DNA,
    ...                  format='fasta')
    ...
    >>> aln
    5 x 2532 dna alignment: Human[TGTGGCACAAA...

Basic Collection objects
^^^^^^^^^^^^^^^^^^^^^^^^

.. _load-seqs:

Constructing a SequenceCollection or Alignment object from strings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> dna  = {'seq1': 'ATGACC',
    ...         'seq2': 'ATCGCC'}
    >>> seqs = LoadSeqs(data=dna, moltype=DNA)
    >>> print type(seqs)
    <class 'cogent.core.alignment.Alignment'>
    >>> seqs = LoadSeqs(data=dna, moltype=DNA, aligned=False)
    >>> print type(seqs)
    <class 'cogent.core.alignment.SequenceCollection'>

To recover a single DNA sequence from the collection or alignment by name
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/long_testseqs.fasta', moltype=DNA, aligned=True)
    >>> aln.Names
    ['Human', 'HowlerMon', 'Mouse', 'NineBande', 'DogFaced']
    >>> seq = aln.getSeq('Human')
    >>> seq.Name
    'Human'
    >>> seq
    DnaSequence(TGTGGCA... 2532)
    >>> type(seq)
    <class 'cogent.core.sequence.DnaSequence'>

One can also slice the sequences from an alignment like a list
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> aln.Seqs[0]
    [0:2532]/2532 of DnaSequence(TGTGGCA... 2532)

An alignment can be sliced "vertically"
"""""""""""""""""""""""""""""""""""""""

Alignments are organised with sequences as 'rows' and aligned residues in 'columns'. Hence, vertical slicing returns columns.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/long_testseqs.fasta', moltype=DNA, aligned=True)
    >>> print aln[:24]
    >Human
    TGTGGCACAAATACTCATGCCAGC
    >HowlerMon
    TGTGGCACAAATACTCATGCCAGC
    >Mouse
    TGTGGCACAGATGCTCATGCCAGC
    >NineBande
    TGTGGCACAAATACTCATGCCAAC
    >DogFaced
    TGTGGCACAAATACTCATGCCAAC
    <BLANKLINE>

A SequenceCollection cannot be sliced (it's unaligned)
""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> seqs = LoadSeqs('data/long_testseqs.fasta', moltype=DNA, aligned=False)
    >>> try:
    ...     print seqs[:24]
    ... except TypeError, e:
    ...     print e
    ...
    'SequenceCollection' object is unsubscriptable


Converting a SequenceCollection to FASTA format
"""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs
    >>> seq = LoadSeqs('data/test.paml', aligned=False)
    >>> fasta_data = seq.toFasta()
    >>> print fasta_data
    >DogFaced
    GCAAGGAGCCAGCAGAACAGATGGGTTGAAACTAAGGAAACATGTAATGATAGGCAGACT
    >HowlerMon
    GCAAGGAGCCAACATAACAGATGGGCTGAAAGTGAGGAAACATGTAATGATAGGCAGACT
    >Human
    GCAAGGAGCCAACATAACAGATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACT
    >Mouse
    GCAGTGAGCCAGCAGAGCAGATGGGCTGCAAGTAAAGGAACATGTAACGACAGGCAGGTT
    >NineBande
    GCAAGGCGCCAACAGAGCAGATGGGCTGAAAGTAAGGAAACATGTAATGATAGGCAGACT

The elements of a collection or alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Accessing individual sequences by name
""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/long_testseqs.fasta', moltype=DNA, aligned=True)
    >>> aln.Names
    ['Human', 'HowlerMon', 'Mouse', 'NineBande', 'DogFaced']
    >>> seq = aln.getSeq('Human')
    >>> seq.Name
    'Human'
    >>> seq
    DnaSequence(TGTGGCA... 2532)
    >>> type(seq)
    <class 'cogent.core.sequence.DnaSequence'>

Accessing individual sequences by position
""""""""""""""""""""""""""""""""""""""""""

The usual approach is to access a ``SequenceCollection`` or ``Alignment`` object as a dictionary, obtaining the individual sequences using the titles as "keys" (above).  However, one can also iterate through the collection like a list.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> fn = 'data/long_testseqs.fasta'
    >>> seqs = LoadSeqs(fn, moltype=DNA, aligned=False)
    >>> my_seq = seqs.Seqs[0]
    >>> my_seq[:24]
    DnaSequence(TGTGGCA... 24)
    >>> str(my_seq[:24])
    'TGTGGCACAAATACTCATGCCAGC'
    >>> type(my_seq)
    <class 'cogent.core.sequence.DnaSequence'>
    >>> aln = LoadSeqs(fn, moltype=DNA, aligned=True)
    >>> aln.Seqs[0][:24]
    [0:24]/2532 of DnaSequence(TGTGGCA... 2532)
    >>> print aln.Seqs[0][:24]
    TGTGGCACAAATACTCATGCCAGC

Keeping a subset of sequences from the alignment
""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/test.paml', moltype=DNA)
    >>> aln.Names
    ['NineBande', 'Mouse', 'Human', 'HowlerMon', 'DogFaced']
    >>> new = aln.takeSeqs(['Human', 'HowlerMon'])
    >>> new.Names
    ['Human', 'HowlerMon']

Note the subset contain references to the original sequences, not copies.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/test.paml', moltype=DNA)
    >>> seq = aln.getSeq('Human')
    >>> new = aln.takeSeqs(['Human', 'HowlerMon'])
    >>> id(new.getSeq('Human')) == id(aln.getSeq('Human'))
    True

Alignments
^^^^^^^^^^

Creating an Alignment object from a SequenceCollection
""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml', aligned=False)
    >>> aln = Alignment(seq)
    >>> fasta_1 = seq.toFasta()
    >>> fasta_2 = aln.toFasta()
    >>> assert fasta_1 == fasta_2

Handling gaps
"""""""""""""

Remove all gaps from an alignment in FASTA format
+++++++++++++++++++++++++++++++++++++++++++++++++

This necessarily returns a ``SequenceCollection``.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs("data/primate_cdx2_promoter.fasta")
    >>> degapped = aln.degap()
    >>> print type(degapped)
    <class 'cogent.core.alignment.SequenceCollection'>

.. TODO the following should be preceded by a section describing the writeToFile method and format argument

Writing sequences to file
"""""""""""""""""""""""""

Both collection and alignment objects have a ``writeToFile`` method. The output format is inferred from the filename suffix,

.. doctest::
    
    >>> from cogent import LoadSeqs, DNA
    >>> dna  = {'seq1': 'ATGACC',
    ...         'seq2': 'ATCGCC'}
    >>> aln = LoadSeqs(data=dna, moltype=DNA)
    >>> aln.writeToFile('sample.fasta')

or by the ``format`` argument.

.. doctest::
    
    >>> aln.writeToFile('sample', format='fasta')

.. now clean the files up

.. doctest::
    :hide:
    
    >>> from cogent.util.misc import remove_files
    >>> remove_files(['sample', 'sample.fasta'], error_on_missing=False)

Converting an alignment to FASTA format
"""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/long_testseqs.fasta')
    >>> aln = Alignment(seq)
    >>> fasta_align = aln.toFasta()

Converting an alignment into Phylip format
""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> phylip_file, name_dictionary = aln.toPhylip()

Converting an alignment to a list of strings
""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> string_list = aln.todict().values()

Slicing an alignment
^^^^^^^^^^^^^^^^^^^^

By rows (sequences)
"""""""""""""""""""

An ``Alignment`` can be sliced

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> fn = 'data/long_testseqs.fasta'
    >>> aln = LoadSeqs(fn, moltype=DNA, aligned=True)
    >>> print aln[:24]
    >Human
    TGTGGCACAAATACTCATGCCAGC
    >HowlerMon
    TGTGGCACAAATACTCATGCCAGC
    >Mouse
    TGTGGCACAGATGCTCATGCCAGC
    >NineBande
    TGTGGCACAAATACTCATGCCAAC
    >DogFaced
    TGTGGCACAAATACTCATGCCAAC
    <BLANKLINE>

but a ``SequenceCollection`` cannot be sliced

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> fn = 'data/long_testseqs.fasta'
    >>> seqs = LoadSeqs(fn, moltype=DNA, aligned=False)
    >>> try:
    ...     print seqs[:24]
    ... except TypeError, e:
    ...     print e
    ...
    'SequenceCollection' object is unsubscriptable

Getting a single column from an Alignment
"""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> seq = LoadSeqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> column_four = aln[3]

Getting a region of contiguous columns
""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent.core.alignment import Alignment
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> region = aln[50:70]

Iterating over alignment positions
""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta')
    >>> col = aln[113:115].iterPositions()
    >>> type(col)
    <type 'generator'>
    >>> list(col)
    [['A', 'A', 'A'], ['T', '-', '-']]

Getting codon 3rd positions from an alignment
"""""""""""""""""""""""""""""""""""""""""""""

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs(data={'seq1': 'ATGATGATG---',
    ...                      'seq2': 'ATGATGATGATG'})
    >>> range(len(aln))[2::3]
    [2, 5, 8, 11]
    >>> indices = [(i, i+1) for i in range(len(aln))[2::3]]
    >>> indices
    [(2, 3), (5, 6), (8, 9), (11, 12)]
    >>> pos3 = aln.addFeature('pos3', 'pos3', indices)
    >>> pos3 = pos3.getSlice()
    >>> print pos3
    >seq2
    GGGG
    >seq1
    GGG-
    <BLANKLINE>

Filtering positions
"""""""""""""""""""

Eliminating columns with non-nucleotide characters
++++++++++++++++++++++++++++++++++++++++++++++++++

We sometimes want to eliminate ambiguous or gap data from our alignments. We show how to exclude alignment columns by the characters they contain. In the first instance we do this just for single nucleotide columns, then for trinucleotides (equivalent for handling codons).

.. doctest::

    >>> from cogent import LoadSeqs, DNA
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

We can also do this in a more longwinded but clearer fashion with a named multi-line function:

.. doctest::

    >>> def just_nucs(x, allowed = 'ACGT'):
    ...     for char in ''.join(x): # ensure char is a str with length 1
    ...         if not char in allowed:
    ...             return False
    ...     return True
    ...
    >>> nucs = aln.filtered(just_nucs)
    >>> nucs
    3 x 8 dna alignment: seq1[ATGAAGGG], seq2[ATGAAGGG], seq3[ATGAAGGG]
    >>> print nucs
    >seq1
    ATGAAGGG
    >seq2
    ATGAAGGG
    >seq3
    ATGAAGGG
    <BLANKLINE>

Applying the same filter to trinucleotides (specified by setting ``motif_length=3``).

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

Getting all variable positions from an alignment
++++++++++++++++++++++++++++++++++++++++++++++++

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

Getting all constant positions from an alignment
++++++++++++++++++++++++++++++++++++++++++++++++

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> just_constant_aln = aln.filtered(lambda x: len(set(x)) == 1)
    >>> print just_constant_aln[:10]
    >Human
    TGTGGCACAA
    >HowlerMon
    TGTGGCACAA
    >Mouse
    TGTGGCACAA
    >NineBande
    TGTGGCACAA
    >DogFaced
    TGTGGCACAA
    <BLANKLINE>

Getting all variable codons from an alignment
+++++++++++++++++++++++++++++++++++++++++++++

This is exactly the same as before, with a new keyword argument

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> variable_codons = aln.filtered(lambda x: len(set(x)) > 1,
    ...                                motif_length=3)
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

Filtering sequences
"""""""""""""""""""

Extracting sequences by sequence identifier into a new alignment object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can use ``takeSeqs()`` to extract some sequences by sequence identifier from an alignment to a new alignment object:

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/long_testseqs.fasta')
    >>> aln.takeSeqs(['Human','Mouse'])
    2 x 2532 text alignment: Human[TGTGGCACAAA...], Mouse[TGTGGCACAGA...]

Alternatively, you can extract only the sequences which are not specified by passing ``negate=True``:

.. doctest::

    >>> aln.takeSeqs(['Human','Mouse'],negate=True)
    3 x 2532 text alignment: NineBande[TGTGGCACAAA...], HowlerMon[TGTGGCACAAA...], DogFaced[TGTGGCACAAA...]

Extracting sequences using an arbitrary function into a new alignment object
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can use ``takeSeqsIf()`` to extract sequences into a new alignment object based on whether an arbitrary function applied to the sequence evaluates to True. For example, to extract sequences which don't contain any N bases you could do the following:

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs(data= [('seq1', 'ATGAAGGTG---'),
    ...                       ('seq2', 'ATGAAGGTGATG'),
    ...                       ('seq3', 'ATGAAGGNGATG')], moltype=DNA)
    >>> def no_N_chars(s):
    ...     return 'N' not in s
    >>> aln.takeSeqsIf(no_N_chars)
    2 x 12 text alignment: seq1[ATGAAGGTG--...], seq2[ATGAAGGTGAT...]

You can additionally get the sequences where the provided function evaluates to False:

.. doctest::

    >>> aln.takeSeqsIf(no_N_chars,negate=True)
    1 x 12 text alignment: seq3[ATGAAGGNGAT...]

Computing alignment statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Computing motif probabilities from an alignment
"""""""""""""""""""""""""""""""""""""""""""""""

The method ``getMotifProbs()`` of ``Alignment`` objects returns the probabilities for all motifs of a given length. For individual nucleotides:

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta', moltype=DNA)
    >>> motif_probs = aln.getMotifProbs()
    >>> print motif_probs
    {'A': 0.24...

For dinucleotides or longer, we need to pass in an ``Alphabet`` with the appropriate word length. Here is an example with trinucleotides:

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> trinuc_alphabet = DNA.Alphabet.getWordAlphabet(3)
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta', moltype=DNA)
    >>> motif_probs = aln.getMotifProbs(alphabet=trinuc_alphabet)
    >>> for m in sorted(motif_probs, key=lambda x: motif_probs[x],
    ...                 reverse=True):
    ...     print m, motif_probs[m]
    ...
    CAG 0.0374581939799
    CCT 0.0341137123746
    CGC 0.0301003344482...

The same holds for other arbitrary alphabets, as long as they match the alignment ``MolType``.

Some calculations in cogent require all non-zero values in the motif probabilities, in which case we use a pseudo-count. We illustrate that here with a simple example where T is missing. Without the pseudo-count, the frequency of T is 0.0, with the pseudo-count defined as 1e-6 then the frequency of T will be slightly less than 1e-6.

.. doctest::

    >>> aln = LoadSeqs(data=[('a', 'AACAAC'),('b', 'AAGAAG')], moltype=DNA)
    >>> motif_probs = aln.getMotifProbs()
    >>> assert motif_probs['T'] == 0.0
    >>> motif_probs = aln.getMotifProbs(pseudocount=1e-6)
    >>> assert 0 < motif_probs['T'] <= 1e-6

It is important to notice that motif probabilities are computed by treating sequences as non-overlapping tuples. Below is a very simple pair of identical sequences where there are clearly 2 'AA' dinucleotides per sequence but only the first one is 'in-frame' (frame width = 2).

We then create a dinucleotide ``Alphabet`` object and use this to get dinucleotide probabilities. These frequencies are determined by breaking each aligned sequence up into non-overlapping dinucleotides and then doing a count. The expected value for the 'AA' dinucleotide in this case will be 2/8 = 0.25.

.. doctest::

    >>> seqs = [('a', 'AACGTAAG'), ('b', 'AACGTAAG')]
    >>> aln = LoadSeqs(data=seqs, moltype=DNA)
    >>> dinuc_alphabet = DNA.Alphabet.getWordAlphabet(2)
    >>> motif_probs = aln.getMotifProbs(alphabet=dinuc_alphabet)
    >>> assert motif_probs['AA'] == 0.25

What about counting the total incidence of dinucleotides including those not in-frame?  A naive application of the Python string object's count method will not work as desired either because it "returns the number of non-overlapping occurrences".

.. doctest::

    >>> seqs = [('my_seq', 'AAAGTAAG')]
    >>> aln = LoadSeqs(data=seqs, moltype=DNA)
    >>> my_seq = aln.getSeq('my_seq')
    >>> my_seq.count('AA')
    2
    >>> 'AAA'.count('AA')
    1
    >>> 'AAAA'.count('AA')
    2

To count all occurrences of a given dinucleotide in a DNA sequence, one could use a standard Python approach such as list comprehension:

.. doctest::

    >>> from cogent import Sequence, DNA
    >>> seq = Sequence(moltype=DNA, seq='AAAGTAAG')
    >>> seq
    DnaSequence(AAAGTAAG)
    >>> di_nucs = [seq[i:i+2] for i in range(len(seq)-1)]
    >>> sum([nn == 'AA' for nn in di_nucs])
    3

Calculating gap fractions for each column in an alignment
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Filtering extracted columns for the gap character
+++++++++++++++++++++++++++++++++++++++++++++++++

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta')
    >>> col = aln[113:115].iterPositions()
    >>> c1, c2 = list(col)
    >>> c1, c2
    (['A', 'A', 'A'], ['T', '-', '-'])
    >>> filter(lambda x: x == '-', c1)
    []
    >>> filter(lambda x: x == '-', c2)
    ['-', '-']

Calculating the gap fraction
++++++++++++++++++++++++++++

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta')
    >>> for column in aln[113:150].iterPositions():
    ...     ungapped = filter(lambda x: x == '-', column)
    ...     gap_fraction = len(ungapped) * 1.0 / len(column)
    ...     print gap_fraction
    0.0
    0.666666666667
    0.0
    0.0...

