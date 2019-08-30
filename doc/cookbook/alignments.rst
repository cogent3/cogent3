Collections and Alignments
--------------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott, Jan Kosinski

For loading collections of unaligned or aligned sequences see :ref:`load-seqs`.

Basic Collection objects
^^^^^^^^^^^^^^^^^^^^^^^^

Constructing a ``SequenceCollection`` or ``Alignment`` object from strings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import make_aligned_seqs, make_unaligned_seqs
    >>> dna  = {'seq1': 'ATGACC',
    ...         'seq2': 'ATCGCC'}
    >>> seqs = make_aligned_seqs(data=dna, moltype="dna")
    >>> print(type(seqs))
    <class 'cogent3.core.alignment.ArrayAlignment'>
    >>> seqs = make_unaligned_seqs(dna, moltype="dna")
    >>> print(type(seqs))
    <class 'cogent3.core.alignment.SequenceCollection'>

Constructing a ``ArrayAlignment`` using ``make_aligned_seqs``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> dna  = {'seq1': 'ATGACC',
    ...         'seq2': 'ATCGCC'}
    >>> seqs = make_aligned_seqs(data=dna, moltype="dna", array_align=True)
    >>> print(type(seqs))
    <class 'cogent3.core.alignment.ArrayAlignment'>
    >>> print(seqs)
    >seq1
    ATGACC
    >seq2
    ATCGCC
    <BLANKLINE>


Converting a ``SequenceCollection`` to FASTA format
"""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_unaligned_seqs
    >>> seqs = load_unaligned_seqs('data/test.paml')
    >>> print(seqs)  # doctest: +SKIP
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

Adding new sequences to an existing collection or alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New sequences can be either appended or inserted using the ``add_seqs`` method. More than one sequence can be added at the same time. Note that ``add_seqs`` does not modify the existing collection/alignment, it creates a new one.

Appending the sequences
"""""""""""""""""""""""

``add_seqs`` without additional parameters will append the sequences to the end of the collection/alignment.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs([('seq1', 'ATGAA------'),
    ...                       ('seq2', 'ATG-AGTGATG'),
    ...                       ('seq3', 'AT--AG-GATG')], moltype="dna")
    >>> print(aln)
    >seq1
    ATGAA------
    >seq2
    ATG-AGTGATG
    >seq3
    AT--AG-GATG
    <BLANKLINE>
    >>> new_seqs = make_aligned_seqs([('seq0', 'ATG-AGT-AGG'),
    ...                           ('seq4', 'ATGCC------')], moltype="dna")
    >>> new_aln = aln.add_seqs(new_seqs)
    >>> print(new_aln)
    >seq1
    ATGAA------
    >seq2
    ATG-AGTGATG
    >seq3
    AT--AG-GATG
    >seq0
    ATG-AGT-AGG
    >seq4
    ATGCC------
    <BLANKLINE>

.. note:: The order is not preserved if you use ``to_fasta`` method, which sorts sequences by name.

Inserting the sequences
"""""""""""""""""""""""

Sequences can be inserted into an alignment at the specified position using either the ``before_name`` or ``after_name`` arguments.

.. doctest::

   >>> new_aln = aln.add_seqs(new_seqs, before_name='seq2')
   >>> print(new_aln)
   >seq1
   ATGAA------
   >seq0
   ATG-AGT-AGG
   >seq4
   ATGCC------
   >seq2
   ATG-AGTGATG
   >seq3
   AT--AG-GATG
   <BLANKLINE>
   >>> new_aln = aln.add_seqs(new_seqs, after_name='seq2')
   >>> print(new_aln)
   >seq1
   ATGAA------
   >seq2
   ATG-AGTGATG
   >seq0
   ATG-AGT-AGG
   >seq4
   ATGCC------
   >seq3
   AT--AG-GATG
   <BLANKLINE>

Inserting sequence(s) based on their alignment to a reference sequence
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Already aligned sequences can be added to an existing ``Alignment`` object and aligned at the same time using the ``add_from_ref_aln`` method. The alignment is performed based on their alignment to a reference sequence (which must be present in both alignments). The method assumes the first sequence in ``ref_aln.names[0]`` is the reference.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs([('seq1', 'ATGAA------'),
    ...                      ('seq2', 'ATG-AGTGATG'),
    ...                      ('seq3', 'AT--AG-GATG')], moltype="dna")
    >>> ref_aln = make_aligned_seqs([('seq3', 'ATAGGATG'),
    ...                          ('seq0', 'ATG-AGCG'),
    ...                          ('seq4', 'ATGCTGGG')], moltype="dna")
    >>> new_aln = aln.add_from_ref_aln(ref_aln)
    >>> print(new_aln)
    >seq1
    ATGAA------
    >seq2
    ATG-AGTGATG
    >seq3
    AT--AG-GATG
    >seq0
    AT--G--AGCG
    >seq4
    AT--GC-TGGG
    <BLANKLINE>

``add_from_ref_aln`` has the same arguments as ``add_seqs`` so ``before_name`` and ``after_name`` can be used to insert the new sequences at the desired position.

.. note:: This method does not work with the ``ArrayAlignment`` class.

Removing all columns with gaps in a named sequence
++++++++++++++++++++++++++++++++++++++++++++++++++

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs([('seq1', 'ATGAA---TG-'),
    ...                      ('seq2', 'ATG-AGTGATG'),
    ...                      ('seq3', 'AT--AG-GATG')], moltype="dna")
    >>> new_aln = aln.get_degapped_relative_to('seq1')
    >>> print(new_aln)
    >seq1
    ATGAATG
    >seq2
    ATG-AAT
    >seq3
    AT--AAT
    <BLANKLINE>

The elements of a collection or alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Accessing individual sequences from a collection or alignment by name
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Using the ``get_seq`` method allows for extracting an unaligned sequence from a collection or alignment by name.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs([('seq1', 'ATGAA------'),
    ...                      ('seq2', 'ATG-AGTGATG'),
    ...                      ('seq3', 'AT--AG-GATG')],
    ...                 moltype="dna", array_align=False)
    >>> seq = aln.get_seq('seq1')
    >>> seq.name
    'seq1'
    >>> type(seq)
    <class 'cogent3.core.sequence.DnaSequence'>
    >>> seq.is_gapped()
    False

Alternatively, if you want to extract the aligned (i.e., gapped) sequence from an alignment, you can use ``get_gapped_seq``.

.. doctest::

    >>> seq = aln.get_gapped_seq('seq1')
    >>> seq.is_gapped()
    True
    >>> print(seq)
    ATGAA------

To see the names of the sequences in a sequence collection, you can use either the ``Names`` attribute or ``get_seq_names`` method.

.. doctest::

    >>> aln.names
    ['seq1', 'seq2', 'seq3']
    >>> aln.names
    ['seq1', 'seq2', 'seq3']

Slice the sequences from an alignment like a list
"""""""""""""""""""""""""""""""""""""""""""""""""

The usual approach is to access a ``SequenceCollection`` or ``Alignment`` object as a dictionary, obtaining the individual sequences using the titles as "keys" (above).  However, one can also iterate through the collection like a list.

.. doctest::

    >>> from cogent3 import load_unaligned_seqs, load_aligned_seqs
    >>> fn = 'data/long_testseqs.fasta'
    >>> seqs = load_unaligned_seqs(fn, moltype="dna")
    >>> my_seq = seqs.seqs[0]
    >>> my_seq[:24]
    DnaSequence(TGTGGCA... 24)
    >>> str(my_seq[:24])
    'TGTGGCACAAATACTCATGCCAGC'
    >>> type(my_seq)
    <class 'cogent3.core.sequence.DnaSequence'>
    >>> aln = load_aligned_seqs(fn, moltype="dna")
    >>> aln.seqs[0][:24]
    DnaSequence(TGTGGCA... 24)
    >>> print(aln.seqs[0][:24])
    TGTGGCACAAATACTCATGCCAGC

Getting a subset of sequences from the alignment
""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/test.paml', moltype="dna")
    >>> aln.names
    ['NineBande', 'Mouse', 'Human', 'HowlerMon', 'DogFaced']
    >>> new = aln.take_seqs(['Human', 'HowlerMon'])
    >>> new.names
    ['Human', 'HowlerMon']

Note, if you set ``array_align=False``, then the subset contain references to the original sequences, not copies.

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/test.paml', array_align=False, moltype="dna")
    >>> seq = aln.get_seq('Human')
    >>> new = aln.take_seqs(['Human', 'HowlerMon'])
    >>> id(new.get_seq('Human')) == id(aln.get_seq('Human'))
    True

Alignments
^^^^^^^^^^

Creating an ``Alignment`` object from a ``SequenceCollection``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_unaligned_seqs
    >>> from cogent3.core.alignment import Alignment
    >>> seq = load_unaligned_seqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> fasta_1 = seq
    >>> fasta_2 = aln
    >>> assert fasta_1 == fasta_2

Convert alignment to DNA, RNA or PROTEIN moltypes
"""""""""""""""""""""""""""""""""""""""""""""""""

This is useful if you've loaded a sequence alignment without specifying the moltype and later need to convert it.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> data = [('a', 'ACG---'), ('b', 'CCTGGG')]
    >>> aln = make_aligned_seqs(data=data)
    >>> dna = aln.to_dna()
    >>> dna
    2 x 6 dna alignment: a[ACG---], b[CCTGGG]

To RNA

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> data = [('a', 'ACG---'), ('b', 'CCUGGG')]
    >>> aln = make_aligned_seqs(data=data)
    >>> rna = aln.to_rna()
    >>> rna
    2 x 6 rna alignment: a[ACG---], b[CCUGGG]

To PROTEIN

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> data = [('x', 'TYV'), ('y', 'TE-')]
    >>> aln = make_aligned_seqs(data=data)
    >>> prot = aln.to_protein()
    >>> prot
    2 x 3 protein alignment: x[TYV], y[TE-]

Handling gaps
"""""""""""""

Remove all gaps from an alignment in FASTA format
+++++++++++++++++++++++++++++++++++++++++++++++++

This necessarily returns a ``SequenceCollection``.

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta")
    >>> degapped = aln.degap()
    >>> print(type(degapped))
    <class 'cogent3.core.alignment.SequenceCollection'>

.. TODO the following should be preceded by a section describing the write method and format argument

Writing sequences to file
"""""""""""""""""""""""""

Both collection and alignment objects have a ``write`` method. The output format is inferred from the filename suffix,

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> dna  = {'seq1': 'ATGACC',
    ...         'seq2': 'ATCGCC'}
    >>> aln = make_aligned_seqs(data=dna, moltype="dna")
    >>> aln.write('sample.fasta')

or by the ``format`` argument.

.. doctest::

    >>> aln.write('sample', format='fasta')

.. now clean the files up

.. doctest::
    :hide:

    >>> from cogent3.util.misc import remove_files
    >>> remove_files(['sample', 'sample.fasta'], error_on_missing=False)

Converting an alignment to FASTA format
"""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> from cogent3.core.alignment import Alignment
    >>> seq = load_aligned_seqs('data/long_testseqs.fasta')
    >>> aln = Alignment(seq)
    >>> fasta_align = aln

Converting an alignment into Phylip format
""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> from cogent3.core.alignment import Alignment
    >>> seq = load_aligned_seqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> got = aln.to_phylip()
    >>> print(got)
    5  60
    NineBande GCAAGGCGCCAACAGAGCAGATGGGCTGAAAGTAAGGAAACATGTAATGATAGGCAGACT
    Mouse     GCAGTGAGCCAGCAGAGCAGATGGGCTGCAAGTAAAGGAACATGTAACGACAGGCAGGTT
    Human     GCAAGGAGCCAACATAACAGATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACT
    HowlerMon GCAAGGAGCCAACATAACAGATGGGCTGAAAGTGAGGAAACATGTAATGATAGGCAGACT
    DogFaced  GCAAGGAGCCAGCAGAACAGATGGGTTGAAACTAAGGAAACATGTAATGATAGGCAGACT
    <BLANKLINE>

Converting an alignment to a list of strings
""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> from cogent3.core.alignment import Alignment
    >>> seq = load_aligned_seqs('data/test.paml')
    >>> aln = Alignment(seq)
    >>> string_list = aln.to_dict().values()

Slicing an alignment
^^^^^^^^^^^^^^^^^^^^

By rows (sequences)
"""""""""""""""""""

An ``Alignment`` can be sliced

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> fn = 'data/long_testseqs.fasta'
    >>> aln = load_aligned_seqs(fn, moltype="dna")
    >>> print(aln[:24])
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

    >>> from cogent3 import load_unaligned_seqs
    >>> fn = 'data/long_testseqs.fasta'
    >>> seqs = load_unaligned_seqs(fn)
    >>> print(seqs[:24])
    Traceback (most recent call last):
    TypeError: 'SequenceCollection' object...

Getting a single column from an alignment
"""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> seq = load_aligned_seqs('data/test.paml')
    >>> column_four = aln[3]

Getting a region of contiguous columns
""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/long_testseqs.fasta')
    >>> region = aln[50:70]

Iterating over alignment positions
""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/primate_cdx2_promoter.fasta')
    >>> col = aln[113:115].iter_positions()
    >>> type(col)
    <class 'generator'>
    >>> list(col)
    [[ByteSequence(A), ByteSequence(A), ByteSequence(A)], [ByteSequence(T)...

Getting codon 3rd positions from ``Alignment``
""""""""""""""""""""""""""""""""""""""""""""""

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data={'seq1': 'ATGATGATG---',
    ...                      'seq2': 'ATGATGATGATG'}, array_align=False)
    >>> list(range(len(aln))[2::3])
    [2, 5, 8, 11]
    >>> indices = [(i, i+1) for i in range(len(aln))[2::3]]
    >>> indices
    [(2, 3), (5, 6), (8, 9), (11, 12)]
    >>> pos3 = aln.add_feature('pos3', 'pos3', indices)
    >>> pos3 = pos3.get_slice()
    >>> print(pos3)  # doctest: +SKIP
    >seq2
    GGGG
    >seq1
    GGG-
    <BLANKLINE>

Getting codon 3rd positions from ``ArrayAlignment``
"""""""""""""""""""""""""""""""""""""""""""""""""""

We can use more conventional slice notation in this instance. Note, because Python counts from 0, the 3rd position starts at index 2.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data={'seq1': 'ATGATGATG---',
    ...                      'seq2': 'ATGATGATGATG'}, array_align=True)
    >>> pos3 = aln[2::3]
    >>> print(pos3)  # doctest: +SKIP
    >seq1
    GGG-
    >seq2
    GGGG
    <BLANKLINE>

.. _filter-positions:

Filtering positions
"""""""""""""""""""

Trim terminal stop codons
+++++++++++++++++++++++++

For evolutionary analyses that use codon models we need to exclude terminating stop codons. For the case where the sequences are all of length divisible by 3.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data={'seq1': 'ACGTAA---',
    ...                      'seq2': 'ACGACA---',
    ...                      'seq3': 'ACGCAATGA'}, moltype="dna")
    ...
    >>> new = aln.trim_stop_codons()
    >>> print(new)  # doctest: +SKIP
    >seq3
    ACGCAA
    >seq2
    ACGACA
    >seq1
    ACG---
    <BLANKLINE>

If the alignment contains sequences not divisible by 3, use the ``allow_partial`` argument.

.. doctest::

    >>> aln = make_aligned_seqs(data={'seq1': 'ACGTAA---',
    ...                      'seq2': 'ACGAC----', # terminal codon incomplete
    ...                      'seq3': 'ACGCAATGA'}, moltype="dna")
    ...
    >>> new = aln.trim_stop_codons(allow_partial=True)
    >>> print(new)  # doctest: +SKIP
    >seq3
    ACGCAA
    >seq2
    ACGAC-
    >seq1
    ACG---
    <BLANKLINE>


Eliminating columns with non-nucleotide characters
++++++++++++++++++++++++++++++++++++++++++++++++++

We sometimes want to eliminate ambiguous or gap data from our alignments. We show how to exclude alignment columns by the characters they contain. In the first instance we do this just for single nucleotide columns, then for trinucleotides (equivalent for handling codons). Both are done using the ``no_degenerates`` method.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data= [('seq1', 'ATGAAGGTG---'),
    ...                       ('seq2', 'ATGAAGGTGATG'),
    ...                       ('seq3', 'ATGAAGGNGATG')], moltype="dna")

We apply to nucleotides,

.. doctest::

    >>> nucs = aln.no_degenerates()
    >>> print(nucs)
    >seq1
    ATGAAGGG
    >seq2
    ATGAAGGG
    >seq3
    ATGAAGGG
    <BLANKLINE>

Applying the same filter to trinucleotides (specified by setting ``motif_length=3``).

.. doctest::

    >>> trinucs = aln.no_degenerates(motif_length=3)
    >>> print(trinucs)
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

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/long_testseqs.fasta')
    >>> pos = aln.variable_positions()
    >>> just_variable_aln = aln.take_positions(pos)
    >>> print(just_variable_aln[:10])
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

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/long_testseqs.fasta')
    >>> pos = aln.variable_positions()
    >>> just_constant_aln = aln.take_positions(pos, negate=True)
    >>> print(just_constant_aln[:10])
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

This is done using the ``filtered`` method using the ``motif_length`` argument. We demonstrate this first for the ``ArrayAlignment``.

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/long_testseqs.fasta')
    >>> variable_codons = aln.filtered(lambda x: len(set(map(tuple, x))) > 1,
    ...                                motif_length=3)
    >>> print(just_variable_aln[:9])
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

Then for the standard ``Alignment`` by first converting the ``ArrayAlignment``.

.. doctest::
    
    >>> aln = aln.to_type(array_align=False)
    >>> variable_codons = aln.filtered(lambda x: len(set(''.join(x))) > 1,
    ...                                motif_length=3)
    >>> print(just_variable_aln[:9])
    >Human
    AAGCAAAAC...

Filtering sequences
"""""""""""""""""""

Extracting sequences by sequence identifier into a new alignment object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can use ``take_seqs`` to extract some sequences by sequence identifier from an alignment to a new alignment object:

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/long_testseqs.fasta')
    >>> aln.take_seqs(['Human','Mouse'])
    2 x 2532 bytes alignment: Human[TGTGGCACAAA...], Mouse[TGTGGCACAGA...]

Alternatively, you can extract only the sequences which are not specified by passing ``negate=True``:

.. doctest::

    >>> aln.take_seqs(['Human','Mouse'], negate=True)  # doctest: +SKIP
    3 x 2532 bytes alignment: NineBande[TGTGGCACAAA...], HowlerMon[TGTGGCACAAA...], DogFaced[TGTGGCACAAA...]

Extracting sequences using an arbitrary function into a new alignment object
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can use ``take_seqs_if`` to extract sequences into a new alignment object based on whether an arbitrary function applied to the sequence evaluates to True. For example, to extract sequences which don't contain any N bases you could do the following:

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data= [('seq1', 'ATGAAGGTG---'),
    ...                       ('seq2', 'ATGAAGGTGATG'),
    ...                       ('seq3', 'ATGAAGGNGATG')], moltype="dna")
    >>> def no_N_chars(s):
    ...     return 'N' not in s
    >>> aln.take_seqs_if(no_N_chars)
    2 x 12 dna alignment: seq1[ATGAAGGTG--...], seq2[ATGAAGGTGAT...]

You can additionally get the sequences where the provided function evaluates to False:

.. doctest::

    >>> aln.take_seqs_if(no_N_chars,negate=True)
    1 x 12 dna alignment: seq3[ATGAAGGNGAT...]

Computing alignment statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Getting motif counts
""""""""""""""""""""

We state the motif length we want and whether to allow gap or ambiguous characters. The latter only has meaning for IPUAC character sets (the DNA, RNA or PROTEIN moltypes). We illustrate this for the DNA moltype with motif lengths of 1 and 3.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data= [('seq1', 'ATGAAGGTG---'),
    ...                       ('seq2', 'ATGAAGGTGATG'),
    ...                       ('seq3', 'ATGAAGGNGATG')], moltype="dna")
    >>> counts = aln.counts()
    >>> print(counts) # doctest: +SKIP
    Counter({'G': 14, 'A': 11, 'T': 7})
    >>> counts = aln.counts(motif_length=3)
    >>> print(counts) # doctest: +SKIP
    Counter({'ATG': 5, 'AAG': 3, 'GTG': 2})
    >>> counts = aln.counts(include_ambiguity=True)
    >>> print(counts) # doctest: +SKIP
    Counter({'G': 14, 'A': 11, 'T': 7, 'N': 1})

.. note::
    
    Only the observed motifs are returned, rather than all defined by the alphabet.

Computing motif probabilities from an alignment
"""""""""""""""""""""""""""""""""""""""""""""""

The method ``get_motif_probs`` of ``Alignment`` objects returns the probabilities for all motifs of a given length. For individual nucleotides:

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/primate_cdx2_promoter.fasta', moltype="dna")
    >>> motif_probs = aln.get_motif_probs()
    >>> print(motif_probs) # doctest: +SKIP
    {'A': 0.24...

For dinucleotides or longer, we need to pass in an ``Alphabet`` with the appropriate word length. Here is an example with trinucleotides:

.. doctest::

    >>> from cogent3 import load_aligned_seqs, DNA
    >>> trinuc_alphabet = DNA.alphabet.get_word_alphabet(3)
    >>> aln = load_aligned_seqs('data/primate_cdx2_promoter.fasta', moltype="dna")
    >>> motif_probs = aln.get_motif_probs(alphabet=trinuc_alphabet)
    >>> for m in sorted(motif_probs, key=lambda x: motif_probs[x],
    ...                 reverse=True):
    ...     print("%s  %.3f" % (m, motif_probs[m]))
    ...
    CAG  0.037
    CCT  0.034
    CGC  0.030...

The same holds for other arbitrary alphabets, as long as they match the alignment ``MolType``.

Some calculations in ``cogent3`` require all non-zero values in the motif probabilities, in which case we use a pseudo-count. We illustrate that here with a simple example where T is missing. Without the pseudo-count, the frequency of T is 0.0, with the pseudo-count defined as 1e-6 then the frequency of T will be slightly less than 1e-6.

.. doctest::

    >>> aln = make_aligned_seqs(data=[('a', 'AACAAC'),('b', 'AAGAAG')], moltype="dna")
    >>> motif_probs = aln.get_motif_probs()
    >>> assert motif_probs['T'] == 0.0
    >>> motif_probs = aln.get_motif_probs(pseudocount=1e-6)
    >>> assert 0 < motif_probs['T'] <= 1e-6

It is important to notice that motif probabilities are computed by treating sequences as non-overlapping tuples. Below is a very simple pair of identical sequences where there are clearly 2 'AA' dinucleotides per sequence but only the first one is 'in-frame' (frame width = 2).

We then create a dinucleotide ``Alphabet`` object and use this to get dinucleotide probabilities. These frequencies are determined by breaking each aligned sequence up into non-overlapping dinucleotides and then doing a count. The expected value for the 'AA' dinucleotide in this case will be 2/8 = 0.25.

.. doctest::

    >>> seqs = [('a', 'AACGTAAG'), ('b', 'AACGTAAG')]
    >>> aln = make_aligned_seqs(data=seqs, moltype="dna")
    >>> dinuc_alphabet = DNA.alphabet.get_word_alphabet(2)
    >>> motif_probs = aln.get_motif_probs(alphabet=dinuc_alphabet)
    >>> assert motif_probs['AA'] == 0.25

What about counting the total incidence of dinucleotides including those not in-frame?  A naive application of the Python string object's count method will not work as desired either because it "returns the number of non-overlapping occurrences".

.. doctest::

    >>> seqs = [('my_seq', 'AAAGTAAG')]
    >>> aln = make_aligned_seqs(data=seqs, moltype="dna")
    >>> my_seq = aln.get_seq('my_seq')
    >>> my_seq.count('AA')
    2
    >>> 'AAA'.count('AA')
    1
    >>> 'AAAA'.count('AA')
    2

To count all occurrences of a given dinucleotide in a DNA sequence, one could use a standard Python approach such as list comprehension:

.. doctest::

    >>> from cogent3 import make_seq
    >>> seq = make_seq(moltype="dna", seq='AAAGTAAG')
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

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/primate_cdx2_promoter.fasta')
    >>> col = aln[113:115].iter_positions()
    >>> c1, c2 = list(col)
    >>> c1, c2
    ([ByteSequence(A), ByteSequence(A), ByteSequence(A)], [ByteSequence(T),...
    >>> list(filter(lambda x: x == '-', c1))
    []
    >>> list(filter(lambda x: x == '-', c2))
    [ByteSequence(-), ByteSequence(-)]

Calculating the gap fraction
++++++++++++++++++++++++++++

.. doctest::

    >>> from cogent3 import load_aligned_seqs
    >>> aln = load_aligned_seqs('data/primate_cdx2_promoter.fasta')
    >>> for column in aln[113:150].iter_positions():
    ...     ungapped = list(filter(lambda x: x == '-', column))
    ...     gap_fraction = len(ungapped) * 1.0 / len(column)
    ...     print(gap_fraction)
    0.0
    0.66666...

Extracting maps of aligned to unaligned positions (i.e., gap maps)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

It's often important to know how an alignment position relates to a position in one or more of the sequences in the alignment. The ``gap_maps`` method of the individual sequences is useful for this. To get a map of sequence to alignment positions for a specific sequence in your alignment, do the following:

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data=[('seq1', 'ATGAAGG-TG--'),
    ...                      ('seq2', 'ATG-AGGTGATG'),
    ...                      ('seq3', 'ATGAAG--GATG')], moltype="dna")
    >>> seq_to_aln_map = aln.get_gapped_seq('seq1').gap_maps()[0]

It's now possible to look up positions in the ``seq1``, and find out what they map to in the alignment:

.. doctest::

    >>> seq_to_aln_map[3]
    3
    >>> seq_to_aln_map[8]
    9

This tells us that in position 3 in ``seq1`` corresponds to position 3 in ``aln``, and that position 8 in ``seq1`` corresponds to position 9 in ``aln``.

Notice that we grabbed the first result from the call to ``gap_maps``. This is the sequence position to alignment position map. The second value returned is the alignment position to sequence position map, so if you want to find out what sequence positions the alignment positions correspond to (opposed to what alignment positions the sequence positions correspond to) for a given sequence, you would take the following steps:

.. doctest::

    >>> aln_to_seq_map = aln.get_gapped_seq('seq1').gap_maps()[1]
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

The ``omit_gap_runs`` method can be applied to remove long stretches of gaps in an alignment. In the following example, we remove sequences that have more than two adjacent gaps anywhere in the aligned sequence.

.. doctest::

    >>> aln = make_aligned_seqs(data=[('seq1', 'ATGAA---TG-'),
    ...                      ('seq2', 'ATG-AGTGATG'),
    ...                      ('seq3', 'AT--AG-GATG')], moltype="dna")
    >>> print(aln.omit_gap_runs(2))  # doctest: +SKIP
    >seq2
    ATG-AGTGATG
    >seq3
    AT--AG-GATG

If instead, we just wanted to remove positions from the alignment which are gaps in more than a certain percentage of the sequences, we could use the ``omit_gap_pos`` function. For example:

.. doctest::

    >>> aln = make_aligned_seqs(data=[('seq1', 'ATGAA---TG-'),
    ...                      ('seq2', 'ATG-AGTGATG'),
    ...                      ('seq3', 'AT--AG-GATG')], moltype="dna")
    >>> print(aln.omit_gap_pos(0.40))  # doctest: +SKIP
    >seq1
    ATGA--TG-
    >seq2
    ATGAGGATG
    >seq3
    AT-AGGATG

You'll notice that the 4th and 7th columns of the alignment have been removed because they contained 66% gaps -- more than the allowed 40%.

If you wanted to remove sequences which contain more than a certain percent gap characters, you could use the ``omit_gap_seqs`` method. This is commonly applied to filter partial sequences from an alignment.

    >>> aln = make_aligned_seqs(data=[('seq1', 'ATGAA------'),
    ...                      ('seq2', 'ATG-AGTGATG'),
    ...                      ('seq3', 'AT--AG-GATG')], moltype="dna")
    >>> filtered_aln = aln.omit_gap_seqs(0.50)
    >>> print(filtered_aln)  # doctest: +SKIP
    >seq2
    ATG-AGTGATG
    >seq3
    AT--AG-GATG

Note that following this call to ``omit_gap_seqs``, the 4th column of ``filtered_aln`` is 100% gaps. This is generally not desirable, so a call to ``omit_gap_seqs`` is frequently followed with a call to ``omit_gap_pos`` with no parameters -- this defaults to removing positions which are all gaps:

    >>> print(filtered_aln.omit_gap_pos())  # doctest: +SKIP
    >seq2
    ATGAGTGATG
    >seq3
    AT-AG-GATG
