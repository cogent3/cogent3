.. jupyter-execute::
    :hide-code:

    import set_working_directory


Sequence Collections and Alignments
-----------------------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott, Jan Kosinski

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"


For loading collections of unaligned or aligned sequences see :ref:`load-seqs`.

Basic Collection objects
^^^^^^^^^^^^^^^^^^^^^^^^

Constructing a ``SequenceCollection`` or ``Alignment`` object from strings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = {"seq1": "ATGACC", "seq2": "ATCGCC"}
    # for an alignment, sequences must be the same length
    seqs = make_aligned_seqs(data=data, moltype="dna")
    type(seqs)

.. jupyter-execute::

    from cogent3 import make_unaligned_seqs

    data = {"seq1": "ATGCC", "seq2": "ATCG"}
    seqs = make_unaligned_seqs(data, moltype="dna")
    type(seqs)


Converting a ``SequenceCollection`` to FASTA format
"""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs

    seqs = load_unaligned_seqs("data/test.paml", moltype="dna")
    print(seqs.to_fasta())

Adding new sequences to an existing collection or alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

More than one sequence can be added to a collection simultaneously. Note that ``add_seqs()`` does not modify the existing collection/alignment, it creates a new one.

The elements of a collection or alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Accessing individual sequences from a collection or alignment by name
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

We can get a sequence by name by indexing the ``.seqs`` attribute.

.. jupyter-execute::

    from cogent3 import make_unaligned_seqs

    seqcoll = make_unaligned_seqs(
        [("seq1", "ATGAA"), ("seq2", "ATGAGTGATG"), ("seq3", "ATAGGATG")],
        moltype="dna",
    )
    seq = seqcoll.seqs["seq1"]
    seq

For an alignment, the result is an `Aligned` instance.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        [("seq1", "ATGAA------"), ("seq2", "ATG-AGTGATG"), ("seq3", "AT--AG-GATG")],
        moltype="dna",
    )
    aligned = aln.seqs["seq1"]
    aligned


For the alignment case, you can get the ungapped sequence by accessing the ``.seq`` attribute of the aligned instance.

.. jupyter-execute::

    aligned.seq


.. jupyter-execute::

    seq = aln.get_seq("seq1")
    seq.name
    type(seq)
    seq.is_gapped()

Alternatively, if you want to extract the aligned (i.e., gapped) sequence from an alignment, you can use ``get_gapped_seq``.

.. jupyter-execute::

    seq = aln.get_gapped_seq("seq1")
    seq.is_gapped()
    seq

To see the names of the sequences in a sequence collection, use the ``names`` attribute.

.. jupyter-execute::

    aln.names

Slice the sequences from a collection like a list
"""""""""""""""""""""""""""""""""""""""""""""""""

Use the ``.seqs`` attribute. We can index a single sequence

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs

    fn = "data/long_testseqs.fasta"
    seqs = load_unaligned_seqs(fn, moltype="dna")
    my_seq = seqs.seqs[0]
    my_seq

but you cannot index a slice (Use ``.take_seqs()`` for that).

.. jupyter-execute::
    :raises: TypeError

    seqs.seqs[:2]

Getting a subset of sequences from the alignment
""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    aln.names

.. jupyter-execute::

    new = aln.take_seqs(["Human", "HowlerMon"])
    new.names

Alternatively, you can extract only the sequences which are not specified by passing ``negate=True``:

.. jupyter-execute::

    new = aln.take_seqs(["Human", "HowlerMon"], negate=True)
    new.names

.. note:: The subset contains references to the original sequences, not copies.


Writing sequences to file
"""""""""""""""""""""""""

Both collection and alignment objects have a ``write()`` method. The output format is inferred from the filename suffix,

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    dna = {"seq1": "ATGACC", "seq2": "ATCGCC"}
    aln = make_aligned_seqs(data=dna, moltype="dna")
    aln.write("sample.fasta")

or by the ``format`` argument.

.. jupyter-execute::

    aln.write("sample", format="fasta")

.. now clean the files up

.. jupyter-execute::

    from cogent3.util.io import remove_files

    remove_files(["sample", "sample.fasta"], error_on_missing=False)

Alignments
^^^^^^^^^^

Creating an ``Alignment`` object from a ``SequenceCollection``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs
    from cogent3.core.alignment import Alignment

    seq = load_unaligned_seqs("data/test.paml", moltype="dna")
    seq

.. jupyter-execute::

    aln = Alignment(seq.seqs)
    aln

Convert alignment to DNA, RNA or PROTEIN moltypes
"""""""""""""""""""""""""""""""""""""""""""""""""

This is useful if you've loaded a sequence alignment without specifying the moltype and later need to convert it using the dedicated method

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = [("a", "ACG---"), ("b", "CCTGGG")]
    aln = make_aligned_seqs(data=data, moltype="text")
    dna = aln.to_dna()
    dna

Or using the generic ``to_moltype()`` method

.. jupyter-execute::

    aln.to_moltype("dna")

To RNA

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = [("a", "ACG---"), ("b", "CCUGGG")]
    aln = make_aligned_seqs(data=data, moltype="text")
    rna = aln.to_rna()
    rna

To PROTEIN

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = [("x", "TYV"), ("y", "TE-")]
    aln = make_aligned_seqs(data=data, moltype="text")
    prot = aln.to_moltype("protein")
    prot

Handling gaps
"""""""""""""

Remove all gaps from an alignment
+++++++++++++++++++++++++++++++++

This necessarily returns a ``SequenceCollection``.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    degapped = aln.degap()
    degapped

Removing all columns with gaps in a named sequence
++++++++++++++++++++++++++++++++++++++++++++++++++

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        [("seq1", "ATGAA---TG-"), ("seq2", "ATG-AGTGATG"), ("seq3", "AT--AG-GATG")],
        moltype="dna",
    )
    new_aln = aln.get_degapped_relative_to("seq1")
    new_aln


.. TODO the following should be preceded by a section describing the write method and format argument

Converting an alignment to FASTA format
"""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    fasta_align = aln.to_fasta()

Converting an alignment to Phylip format
""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    got = aln.to_phylip()
    print(got)

Converting an alignment to a list of strings
""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    data = [str(s) for s in aln.seqs]

Converting an alignment to a list of arrays
"""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    import numpy
    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    data = [numpy.array(s) for s in aln.seqs]
    data

Getting an alignment as an array
""""""""""""""""""""""""""""""""

The rows are sequences and their order corresponds to that of ``aln.names``.

.. jupyter-execute::

    import numpy
    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    aln.array_seqs

Slicing an alignment
""""""""""""""""""""

Alignments can be thought of as a matrix, with sequences along the rows and alignment positions as the columns. However, all slicing is only along positions.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    fn = "data/long_testseqs.fasta"
    aln = load_aligned_seqs(fn, moltype="dna")
    aln[:24]

.. warning:: A ``SequenceCollection`` cannot be sliced!

Getting a single column from an alignment
"""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    column_four = aln[3]
    column_four

Getting a region of contiguous columns
""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    region = aln[50:70]

Iterating over alignment positions
""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    col = list(aln[113:115].iter_positions())
    col

Getting codon 3rd positions
"""""""""""""""""""""""""""

We can use conventional slice notation. Note, because Python counts from 0, the 3rd position starts at index 2.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data={"seq1": "ATGATGATG---", "seq2": "ATGATGATGATG"},
        moltype="dna",
    )
    aln[2::3]

.. _filter-positions:

Filtering positions
"""""""""""""""""""

Trim terminal stop codons
+++++++++++++++++++++++++

For evolutionary analyses that use codon models we need to exclude terminating stop codons. For the case where the sequences are all of length divisible by 3.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data={"seq1": "ACGTAA---", "seq2": "ACGACA---", "seq3": "ACGCAATGA"},
        moltype="dna",
    )
    new = aln.trim_stop_codons()
    new

To detect if the alignment contains sequences not divisible by 3, use the ``strict`` argument. This argument covers both allowing partial terminating codons / not divisible by 3.

.. jupyter-execute::
    :raises:

    aln = make_aligned_seqs(
        data={
            "seq1": "ACGTAA---",
            "seq2": "ACGAC----",  # terminal codon incomplete
            "seq3": "ACGCAATGA",
        },
        moltype="dna",
    )
    new = aln.trim_stop_codons(strict=True)

Eliminating columns with non-nucleotide characters
++++++++++++++++++++++++++++++++++++++++++++++++++

We sometimes want to eliminate ambiguous or gap data from our alignments. We demonstrate how to exclude alignment columns based on the characters they contain. In the first instance, we do this just for single nucleotide columns, then for trinucleotides (equivalent for handling codons). Both are done using the ``no_degenerates()`` method.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[
            ("seq1", "ATGAAGGTG---"),
            ("seq2", "ATGAAGGTGATG"),
            ("seq3", "ATGAAGGNGATG"),
        ],
        moltype="dna",
    )

We apply to nucleotides,

.. jupyter-execute::

    nucs = aln.no_degenerates()
    nucs

Applying the same filter to trinucleotides (specified by setting ``motif_length=3``).

.. jupyter-execute::

    trinucs = aln.no_degenerates(motif_length=3)
    trinucs

Getting all variable positions from an alignment
++++++++++++++++++++++++++++++++++++++++++++++++

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    pos = aln.variable_positions()
    just_variable_aln = aln.take_positions(pos)
    just_variable_aln[:10]

Getting all constant positions from an alignment
++++++++++++++++++++++++++++++++++++++++++++++++

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    pos = aln.variable_positions()
    just_constant_aln = aln.take_positions(pos, negate=True)
    just_constant_aln[:10]

Getting all variable codons from an alignment
+++++++++++++++++++++++++++++++++++++++++++++

This is done using the ``filtered`` method using the ``motif_length`` argument.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    pos = aln.variable_positions(motif_length=3)
    variable_codons = aln.take_positions(pos)
    just_variable_aln[:9]

Filtering sequences
"""""""""""""""""""

Extracting sequences using an arbitrary function into a new alignment object
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can use ``take_seqs_if`` to extract sequences into a new alignment object based on whether an arbitrary function applied to the sequence evaluates to True. For example, to extract sequences which don't contain any N bases you could do the following:

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[
            ("seq1", "ATGAAGGTG---"),
            ("seq2", "ATGAAGGTGATG"),
            ("seq3", "ATGAAGGNGATG"),
        ],
        moltype="dna",
    )

    def no_N_chars(s):
        return "N" not in str(s)

    aln.take_seqs_if(no_N_chars)

You can additionally get the sequences where the provided function evaluates to False:

.. jupyter-execute::

    aln.take_seqs_if(no_N_chars, negate=True)

Computing alignment statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Getting motif counts
""""""""""""""""""""

We state the motif length we want and whether to allow gap or ambiguous characters. The latter only has meaning for IPUAC character sets (the DNA, RNA or PROTEIN moltypes). We illustrate this for the DNA moltype with motif lengths of 1 and 3.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[
            ("seq1", "ATGAAGGTG---"),
            ("seq2", "ATGAAGGTGATG"),
            ("seq3", "ATGAAGGNGATG"),
        ],
        moltype="dna",
    )
    counts = aln.counts()
    counts

.. jupyter-execute::

    counts = aln.counts(motif_length=3)
    counts

.. jupyter-execute::

    counts = aln.counts(include_ambiguity=True)
    counts

.. note::

    Only the observed motifs are returned, rather than all defined by the alphabet.

Getting motif counts per sequence
"""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[
            ("seq1", "ATGAAGGTG---"),
            ("seq2", "ATGAAGGTGATG"),
            ("seq3", "ATGAAGGNGATG"),
        ],
        moltype="dna",
    )
    counts = aln.counts_per_seq()
    counts

.. note:: There are also ``.probs_per_seq()`` and ``.entropy_per_seq()`` methods.

Getting motif counts per position
"""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[
            ("seq1", "ATGAAGGTG---"),
            ("seq2", "ATGAAGGTGATG"),
            ("seq3", "ATGAAGGNGATG"),
        ],
        moltype="dna",
    )
    counts = aln.counts_per_pos()
    counts

.. note:: There are also ``.probs_per_pos()`` and ``.entropy_per_pos()`` methods.

Computing motif probabilities from an alignment
"""""""""""""""""""""""""""""""""""""""""""""""

The method ``get_motif_probs`` of ``Alignment`` objects returns the probabilities for all motifs of a given length. For individual nucleotides:

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    motif_probs = aln.get_motif_probs()
    motif_probs

For dinucleotides or longer, we need to pass in a ``KmerAlphabet`` with the appropriate word length. Here is an example with trinucleotides:

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, get_moltype, make_table

    trinuc_alphabet = get_moltype("dna").alphabet.get_kmer_alphabet(3)
    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    motif_probs = aln.get_motif_probs(alphabet=trinuc_alphabet)
    table = make_table(header=["motif", "freq"], data=list(motif_probs.items()))
    table

Some calculations in ``cogent3`` require all non-zero values in the motif probabilities, in which case we use a pseudo-count. We illustrate that here with a simple example where T is missing. Without the pseudo-count, the frequency of T is 0.0, with the pseudo-count defined as 1e-6 then the frequency of T will be slightly less than 1e-6.

.. jupyter-execute::

    aln = make_aligned_seqs(data=[("a", "AACAAC"), ("b", "AAGAAG")], moltype="dna")
    motif_probs = aln.get_motif_probs()
    assert motif_probs["T"] == 0.0
    motif_probs = aln.get_motif_probs(pseudocount=1e-6)
    motif_probs

.. note:: For alignments, motif probabilities are computed by treating sequences as non-overlapping tuples. To get all possible k-mers, use the ``iter_kmers()`` method on the sequence classes.

Working with alignment gaps
"""""""""""""""""""""""""""

Filtering extracted columns for the gap character
+++++++++++++++++++++++++++++++++++++++++++++++++

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    col = aln[113:115].iter_positions()
    c1, c2 = list(col)
    c1, c2
    list(filter(lambda x: x == "-", c1))
    list(filter(lambda x: x == "-", c2))

Calculating the gaps per position
+++++++++++++++++++++++++++++++++

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    gap_counts = aln.count_gaps_per_pos()
    gap_counts # this is a DictArray

To turn that into grap fraction

.. jupyter-execute::

    gap_frac = gap_counts.array / aln.num_seqs


Filtering alignments based on gaps
++++++++++++++++++++++++++++++++++

If we want to remove positions from the alignment which are gaps in more than a certain percentage of the sequences, use the ``omit_gap_pos()``.

.. jupyter-execute::

    aln = make_aligned_seqs(
        data=[
            ("seq1", "ATGAA---TG-"),
            ("seq2", "ATG-AGTGATG"),
            ("seq3", "AT--AG-GATG"),
        ],
        moltype="dna",
    )
    filtered_aln = aln.omit_gap_pos(allowed_gap_frac=0.40)
    filtered_aln

.. note:: The default for ``filtered_aln.omit_gap_pos()`` is to remove columns with gaps in all the sequences. This can occur after sequences have been removed from the alignment.
