.. jupyter-execute::
    :hide-code:

    import set_working_directory

Sequence Collections and Alignments
-----------------------------------

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul, Tom Elliott, Jan Kosinski

.. note:: **Alpha Release of the New SequenceCollection API**

   We are pleased to announce an alpha release of our new ``SequenceCollection`` API. This version can be accessed by specifying the argument ``new_type=True`` in the ``load_unaligned_seqs()`` or ``make_unaligned_seqs()`` functions. The new API is designed for increased efficiency, offering access to the underlying data in multiple formats, including numpy arrays, strings, and bytes (via ``array(seqs)``, ``str(seqs)`` and ``bytes(seqs)`` respectively). 

   Please be aware that this alpha release has not been fully integrated with the library. Users are encouraged to explore its capabilities but should proceed with caution!

For loading collections of unaligned or aligned sequences see :ref:`load-seqs`.

What's the difference between ``Alignment`` and ``ArrayAlignment``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``Alignment`` class can be annotated, meaning you can add annotations to an Alignment or it's member sequences and you can slice the alignment via those objects. This capability is achieved, under the hood, by having the individual sequences represent gaps as a "span", rather than explicitly as a "-" character in the sequence itself. This representation is also efficient for very long sequences.

The ``ArrayAlignment`` class cannot be annotated. The class represents its sequences as a ``numpy.ndarray`` instance. Thus, the gaps are included as a state in the array. This class is better at handling a lot of sequences and should typically be faster. This is the default class returned by the ``load_aligned_seqs`` and ``make_aligned_seqs`` functions. (See :ref:`load-seqs` for details.)

You can change alignment types using the ``to_type()`` method.

Basic Collection objects
^^^^^^^^^^^^^^^^^^^^^^^^

Constructing a ``SequenceCollection`` or ``Alignment`` object from strings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_aligned_seqs, make_unaligned_seqs

    dna = {"seq1": "ATGACC", "seq2": "ATCGCC"}
    seqs = make_aligned_seqs(data=dna, moltype="dna")
    type(seqs)

.. jupyter-execute::

    seqs = make_unaligned_seqs(dna, moltype="dna")
    type(seqs)

Constructing a ``ArrayAlignment`` using ``make_aligned_seqs``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    dna = {"seq1": "ATGACC", "seq2": "ATCGCC"}
    seqs = make_aligned_seqs(data=dna, moltype="dna", array_align=True)
    print(type(seqs))
    seqs

Converting a ``SequenceCollection`` to FASTA format
"""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs

    seqs = load_unaligned_seqs("data/test.paml", moltype="dna")
    seqs

Adding new sequences to an existing collection or alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New sequences can be either appended or inserted using the ``add_seqs`` method. More than one sequence can be added at the same time. Note that ``add_seqs`` does not modify the existing collection/alignment, it creates a new one.

Appending the sequences
"""""""""""""""""""""""

``add_seqs`` without additional parameters will append the sequences to the end of the collection/alignment.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        [("seq1", "ATGAA------"), ("seq2", "ATG-AGTGATG"), ("seq3", "AT--AG-GATG")],
        moltype="dna",
    )
    aln

.. jupyter-execute::

    new_seqs = make_aligned_seqs(
        [("seq0", "ATG-AGT-AGG"), ("seq4", "ATGCC------")], moltype="dna"
    )
    new_aln = aln.add_seqs(new_seqs)
    new_aln

.. note:: The order is not preserved if you use the ``to_fasta()`` method, which sorts sequences by name.

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

The elements of a collection or alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Accessing individual sequences from a collection or alignment by name
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Using the ``get_seq`` method allows for extracting an unaligned sequence from a collection or alignment by name.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        [("seq1", "ATGAA------"), ("seq2", "ATG-AGTGATG"), ("seq3", "AT--AG-GATG")],
        moltype="dna",
        array_align=False,
    )
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

Slice the sequences from an alignment like a list
"""""""""""""""""""""""""""""""""""""""""""""""""

The usual approach is to access a ``SequenceCollection`` or ``Alignment`` object as a dictionary, obtaining the individual sequences using the titles as "keys" (above).  However, one can also iterate through the collection like a list.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, load_unaligned_seqs

    fn = "data/long_testseqs.fasta"
    seqs = load_unaligned_seqs(fn, moltype="dna")
    my_seq = seqs.seqs[0]
    my_seq[:24]

.. jupyter-execute::

    type(my_seq)

.. jupyter-execute::

    aln = load_aligned_seqs(fn, moltype="dna")
    aln.seqs[0][:24]


Getting a subset of sequences from the alignment
""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    aln.names

.. jupyter-execute::

    new = aln.take_seqs(["Human", "HowlerMon"])
    new.names

.. note:: The ``Alignment`` class (which you get if you set ``array_align=False``) is more memory efficient. The subset contain references to the original sequences, not copies.

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
    aln = make_aligned_seqs(data=data, moltype="dna")
    dna = aln.to_dna()
    dna

Or using the generic ``to_moltype()`` method

.. jupyter-execute::

    aln.to_moltype("dna")

To RNA

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    data = [("a", "ACG---"), ("b", "CCUGGG")]
    aln = make_aligned_seqs(data=data, moltype="dna")
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

Remove all gaps from an alignment in FASTA format
+++++++++++++++++++++++++++++++++++++++++++++++++

This necessarily returns a ``SequenceCollection``.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    degapped = aln.degap()
    print(type(degapped))

.. TODO the following should be preceded by a section describing the write method and format argument

Writing sequences to file
"""""""""""""""""""""""""

Both collection and alignment objects have a ``write`` method. The output format is inferred from the filename suffix,

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

Converting an alignment to FASTA format
"""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.core.alignment import Alignment

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    fasta_align = aln.to_fasta()

Converting an alignment into Phylip format
""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.core.alignment import Alignment

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    got = aln.to_phylip()
    print(got)

Converting an alignment to a list of strings
""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.core.alignment import Alignment

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    string_list = aln.to_dict().values()

Slicing an alignment
^^^^^^^^^^^^^^^^^^^^

By rows (sequences)
"""""""""""""""""""

An ``Alignment`` can be sliced

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    fn = "data/long_testseqs.fasta"
    aln = load_aligned_seqs(fn, moltype="dna")
    aln[:24]

but a ``SequenceCollection`` cannot be sliced

.. jupyter-execute::
    :raises: TypeError

    from cogent3 import load_unaligned_seqs

    fn = "data/long_testseqs.fasta"
    seqs = load_unaligned_seqs(fn, moltype="dna")
    seqs[:24]

Getting a single column from an alignment
"""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/test.paml", moltype="dna")
    column_four = aln[3]

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
    col = aln[113:115].iter_positions()
    type(col)
    list(col)

Getting codon 3rd positions from ``Alignment``
""""""""""""""""""""""""""""""""""""""""""""""

We'll do this by specifying the position indices of interest, creating a sequence ``Feature`` and using that to extract the positions.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data={"seq1": "ATGATGATG---", "seq2": "ATGATGATGATG"},
        array_align=False,
        moltype="dna",
    )
    list(range(len(aln))[2::3])
    indices = [(i, i + 1) for i in range(len(aln))[2::3]]
    indices

.. jupyter-execute::

    pos3 = aln.add_feature(biotype="pos3", name="pos3", spans=indices)
    pos3 = pos3.get_slice()
    pos3

Getting codon 3rd positions from ``ArrayAlignment``
"""""""""""""""""""""""""""""""""""""""""""""""""""

We can use more conventional slice notation in this instance. Note, because Python counts from 0, the 3rd position starts at index 2.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data={"seq1": "ATGATGATG---", "seq2": "ATGATGATGATG"},
        array_align=True,
        moltype="dna",
    )
    pos3 = aln[2::3]
    pos3

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

We sometimes want to eliminate ambiguous or gap data from our alignments. We show how to exclude alignment columns by the characters they contain. In the first instance we do this just for single nucleotide columns, then for trinucleotides (equivalent for handling codons). Both are done using the ``no_degenerates`` method.

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
    variable_codons = aln.filtered(lambda x: len({str(e) for e in x}) > 1, motif_length=3)
    just_variable_aln[:9]

Filtering sequences
"""""""""""""""""""

Extracting sequences by sequence identifier into a new alignment object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can use ``take_seqs`` to extract some sequences by sequence identifier from an alignment to a new alignment object:

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")
    aln.take_seqs(["Human", "Mouse"])

Alternatively, you can extract only the sequences which are not specified by passing ``negate=True``:

.. jupyter-execute::

    aln.take_seqs(["Human", "Mouse"], negate=True)

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
        return "N" not in s

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

Computing motif probabilities from an alignment
"""""""""""""""""""""""""""""""""""""""""""""""

The method ``get_motif_probs`` of ``Alignment`` objects returns the probabilities for all motifs of a given length. For individual nucleotides:

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    motif_probs = aln.get_motif_probs()
    motif_probs

For dinucleotides or longer, we need to pass in an ``Alphabet`` with the appropriate word length. Here is an example with trinucleotides:

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, get_moltype

    trinuc_alphabet = get_moltype("dna").alphabet.get_word_alphabet(3)
    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    motif_probs = aln.get_motif_probs(alphabet=trinuc_alphabet)
    for m in sorted(motif_probs, key=lambda x: motif_probs[x], reverse=True):
        print("%s  %.3f" % (m, motif_probs[m]))

The same holds for other arbitrary alphabets, as long as they match the alignment ``MolType``.

Some calculations in ``cogent3`` require all non-zero values in the motif probabilities, in which case we use a pseudo-count. We illustrate that here with a simple example where T is missing. Without the pseudo-count, the frequency of T is 0.0, with the pseudo-count defined as 1e-6 then the frequency of T will be slightly less than 1e-6.

.. jupyter-execute::

    aln = make_aligned_seqs(data=[("a", "AACAAC"), ("b", "AAGAAG")], moltype="dna")
    motif_probs = aln.get_motif_probs()
    assert motif_probs["T"] == 0.0
    motif_probs = aln.get_motif_probs(pseudocount=1e-6)
    assert 0 < motif_probs["T"] <= 1e-6

It is important to notice that motif probabilities are computed by treating sequences as non-overlapping tuples. Below is a very simple pair of identical sequences where there are clearly 2 'AA' dinucleotides per sequence but only the first one is 'in-frame' (frame width = 2).

We then create a dinucleotide ``Alphabet`` object and use this to get dinucleotide probabilities. These frequencies are determined by breaking each aligned sequence up into non-overlapping dinucleotides and then doing a count. The expected value for the 'AA' dinucleotide in this case will be 2/8 = 0.25.

.. jupyter-execute::

    from cogent3 import get_moltype

    dna = get_moltype("dna")
    seqs = [("a", "AACGTAAG"), ("b", "AACGTAAG")]
    aln = make_aligned_seqs(data=seqs, moltype="dna")
    dinuc_alphabet = dna.alphabet.get_word_alphabet(2)
    motif_probs = aln.get_motif_probs(alphabet=dinuc_alphabet)
    assert motif_probs["AA"] == 0.25

What about counting the total incidence of dinucleotides including those not in-frame?  A naive application of the Python string object's count method will not work as desired either because it "returns the number of non-overlapping occurrences".

.. jupyter-execute::

    seqs = [("my_seq", "AAAGTAAG")]
    aln = make_aligned_seqs(data=seqs, moltype="dna")
    my_seq = aln.get_seq("my_seq")
    my_seq.count("AA")
    "AAA".count("AA")
    "AAAA".count("AA")

To count all occurrences of a given dinucleotide in a DNA sequence, one could use a standard Python approach such as list comprehension:

.. jupyter-execute::

    from cogent3 import make_seq

    seq = make_seq(moltype="dna", seq="AAAGTAAG")
    seq
    di_nucs = [seq[i : i + 2] for i in range(len(seq) - 1)]
    sum([nn == "AA" for nn in di_nucs])

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

Calculating the gap fraction
++++++++++++++++++++++++++++

.. jupyter-execute::

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta", moltype="dna")
    for column in aln[113:150].iter_positions():
        ungapped = list(filter(lambda x: x == "-", column))
        gap_fraction = len(ungapped) * 1.0 / len(column)
        print(gap_fraction)

Filtering alignments based on gaps
++++++++++++++++++++++++++++++++++

If we want to remove positions from the alignment which are gaps in more than a certain percentage of the sequences, we could use the ``omit_gap_pos`` function. For example:

.. jupyter-execute::

    aln = make_aligned_seqs(
        data=[
            ("seq1", "ATGAA---TG-"),
            ("seq2", "ATG-AGTGATG"),
            ("seq3", "AT--AG-GATG"),
        ],
        moltype="dna",
    )
    aln.omit_gap_pos(0.40)

.. note:: The default for ``filtered_aln.omit_gap_pos()`` is to remove columns with gaps in all the sequences. This can occur after sequences have been removed from the alignment.
