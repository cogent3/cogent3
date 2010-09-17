Collections and Alignments
--------------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott

For loading collections of unaligned or aligned sequences see :ref:`load-seqs`.

Basic Collection objects
^^^^^^^^^^^^^^^^^^^^^^^^

Constructing a ``SequenceCollection`` or ``Alignment`` object from strings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

Converting a ``SequenceCollection`` to FASTA format
"""""""""""""""""""""""""""""""""""""""""""""""""""

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

Accessing individual sequences from a collection or alignment by name
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Using the ``getSeq`` method allows for extracting an unaligned sequence from a collection or alignment by name.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs(data= [('seq1', 'ATGAA------'),
    ...                       ('seq2', 'ATG-AGTGATG'),
    ...                       ('seq3', 'AT--AG-GATG')], moltype=DNA)
    >>> seq = aln.getSeq('seq1')
    >>> seq.Name
    'seq1'
    >>> type(seq)
    <class 'cogent.core.sequence.DnaSequence'>
    >>> seq.isGapped()
    False

Alternatively, if you want to extract the aligned (i.e., gapped) sequence from an alignment, you can use ``getGappedSeq``.

.. doctest::

    >>> seq = aln.getGappedSeq('seq1')
    >>> seq.isGapped()
    True
    >>> print seq
    ATGAA------

To see the names of the sequences in a sequence collection, you can use either the ``Names`` attribute or ``getSeqNames`` method.

.. doctest::

    >>> aln.Names
    ['seq1', 'seq2', 'seq3']
    >>> aln.getSeqNames()
    ['seq1', 'seq2', 'seq3']

Slice the sequences from an alignment like a list
"""""""""""""""""""""""""""""""""""""""""""""""""

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

Getting a subset of sequences from the alignment
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

Creating an ``Alignment`` object from a ``SequenceCollection``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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
    >>> print seqs[:24]
    Traceback (most recent call last):
    TypeError: 'SequenceCollection' object is unsubscriptable

Getting a single column from an alignment
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

.. _filter-positions:

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

You can use ``takeSeqs`` to extract some sequences by sequence identifier from an alignment to a new alignment object:

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

You can use ``takeSeqsIf`` to extract sequences into a new alignment object based on whether an arbitrary function applied to the sequence evaluates to True. For example, to extract sequences which don't contain any N bases you could do the following:

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

The method ``getMotifProbs`` of ``Alignment`` objects returns the probabilities for all motifs of a given length. For individual nucleotides:

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

Working with alignment gaps
"""""""""""""""""""""""""""

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

Extracting maps of aligned to unaligned positions (i.e., gap maps)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

It's often important to know how an alignment position relates to a position in one or more of the sequences in the alignment. The ``gapMaps`` method of the individual sequences is useful for this. To get a map of sequence to alignment positions for a specific sequence in your alignment, do the following:

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs(data= [('seq1', 'ATGAAGG-TG--'),
    ...                       ('seq2', 'ATG-AGGTGATG'),
    ...                       ('seq3', 'ATGAAG--GATG')], moltype=DNA)
    >>> seq_to_aln_map = aln.getGappedSeq('seq1').gapMaps()[0]

It's now possible to look up positions in the ``seq1``, and find out what they map to in the alignment:

.. doctest::

    >>> seq_to_aln_map[3]
    3
    >>> seq_to_aln_map[8]
    9

This tells us that in position 3 in ``seq1`` corresponds to position 3 in ``aln``, and that position 8 in ``seq1`` corresponds to position 9 in ``aln``.

Notice that we grabbed the first result from the call to ``gapMaps``. This is the sequence position to alignment position map. The second value returned is the alignment position to sequence position map, so if you want to find out what sequence positions the alignment positions correspond to (opposed to what alignment positions the sequence positions correspond to) for a given sequence, you would take the following steps:

.. doctest::

    >>> aln_to_seq_map = aln.getGappedSeq('seq1').gapMaps()[1]
    >>> aln_to_seq_map[3]
    3
    >>> aln_to_seq_map[8]
    7

If an alignment position is a gap, and therefore has no corresponding sequence position, you'll get a ``KeyError``.

.. doctest::

   >>> seq_pos = aln_to_seq_map[7]
   Traceback (most recent call last):
   KeyError: 7

.. note:: The first position in alignments and sequences is always numbered position 0.

Filtering alignments based on gaps
++++++++++++++++++++++++++++++++++

.. note:: An alternate, computationally faster, approach to removing gaps is to use the ``filtered`` method as discussed in :ref:`filter-positions`.

The ``omitGapRuns`` method can be applied to remove long stretches of gaps in an alignment. In the following example, we remove sequences that have more than two adjacent gaps anywhere in the aligned sequence.

.. doctest::

    >>> aln = LoadSeqs(data= [('seq1', 'ATGAA---TG-'),
    ...                       ('seq2', 'ATG-AGTGATG'),
    ...                       ('seq3', 'AT--AG-GATG')], moltype=DNA)
    >>> print aln.omitGapRuns(2).toFasta()
    >seq2
    ATG-AGTGATG
    >seq3
    AT--AG-GATG

If instead, we just wanted to remove positions from the alignment which are gaps in more than a certain percentage of the sequences, we could use the ``omitGapPositions`` function. For example:

.. doctest::

    >>> aln = LoadSeqs(data= [('seq1', 'ATGAA---TG-'),
    ...                       ('seq2', 'ATG-AGTGATG'),
    ...                       ('seq3', 'AT--AG-GATG')], moltype=DNA)
    >>> print aln.omitGapPositions(0.40).toFasta()
    >seq1
    ATGA--TG-
    >seq2
    ATGAGGATG
    >seq3
    AT-AGGATG

You'll notice that the 4th and 7th columns of the alignment have been removed because they contained 66% gaps -- more than the allowed 40%. 

If you wanted to remove sequences which contain more than a certain percent gap characters, you could use the ``omitGapSeqs`` method. This is commonly applied to filter partial sequences from an alignment. 

    >>> aln = LoadSeqs(data= [('seq1', 'ATGAA------'),
    ...                       ('seq2', 'ATG-AGTGATG'),
    ...                       ('seq3', 'AT--AG-GATG')], moltype=DNA)
    >>> filtered_aln = aln.omitGapSeqs(0.50)
    >>> print filtered_aln.toFasta()
    >seq2
    ATG-AGTGATG
    >seq3
    AT--AG-GATG

Note that following this call to ``omitGapSeqs``, the 4th column of ``filtered_aln`` is 100% gaps. This is generally not desirable, so a call to ``omitGapSeqs`` is frequently followed with a call to ``omitGapPositions`` with no parameters -- this defaults to removing positions which are all gaps:

    >>> print filtered_aln.omitGapPositions().toFasta()
    >seq2
    ATGAGTGATG
    >seq3
    AT-AG-GATG

