*******************
Building alignments
*******************

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Using the cogent3 aligners
==========================

Running a pairwise Needleman-Wunsch-Alignment
---------------------------------------------

.. TODO look at the singapore workshop usage of cogent3.align.align

Running a progressive aligner
-----------------------------

We import useful functions and then load the sequences to be aligned.

.. doctest::

    >>> from cogent3 import load_unaligned_seqs, make_tree
    >>> seqs = load_unaligned_seqs('data/test2.fasta', moltype="dna")

For nucleotides
^^^^^^^^^^^^^^^

We load a canned nucleotide substitution model and the progressive aligner ``TreeAlign`` function.

.. doctest::

    >>> from cogent3.evolve.models import HKY85
    >>> from cogent3.align.progressive import TreeAlign

We first align without providing a guide tree. The ``TreeAlign`` algorithm builds pairwise alignments and estimates the substitution model parameters and pairwise distances. The distances are used to build a neighbour joining tree and the median value of substitution model parameters are provided to the substitution model for the progressive alignment step.

.. doctest::

    >>> aln, tree = TreeAlign(HKY85(), seqs)
    >>> aln # doctest: +SKIP
    5 x 60 bytes alignment: NineBande[-C-----GCCA...], Mouse[GCAGTGAGCCA...], DogFaced[GCAAGGAGCCA...], ...

We then align using a guide tree (pre-estimated) and specifying the ratio of transitions to transversions (kappa).

.. doctest::

    >>> tree = make_tree('(((NineBande:0.0128202449453,Mouse:0.184732725695):0.0289459522137,DogFaced:0.0456427810916):0.0271363715538,Human:0.0341320714654,HowlerMon:0.0188456837006)root;')
    >>> params={'kappa': 4.0}
    >>> aln, tree = TreeAlign(HKY85(), seqs, tree=tree, param_vals=params)
    >>> aln # doctest: +SKIP
    5 x 60 bytes alignment: NineBande[-C-----GCCA...], Mouse[GCAGTGAGCCA...], DogFaced[GCAAGGAGCCA...], ...

For codons
^^^^^^^^^^

We load a canned codon substitution model and use a pre-defined tree and parameter estimates.

.. doctest::

    >>> from cogent3.evolve.models import MG94HKY
    >>> tree = make_tree('((NineBande:0.0575781680031,Mouse:0.594704139406):0.078919659556,DogFaced:0.142151930069,(HowlerMon:0.0619991555435,Human:0.10343006422):0.0792423439112)')
    >>> params={'kappa': 4.0, 'omega': 1.3}
    >>> aln, tree = TreeAlign(MG94HKY(), seqs, tree=tree, param_vals=params)
    >>> aln # doctest: +SKIP
    5 x 60 bytes alignment: NineBande[------CGCCA...], Mouse[GCAGTGAGCCA...], DogFaced[GCAAGGAGCCA...], ...

Converting gaps from aa-seq alignment to nuc seq alignment
==========================================================

We load some unaligned DNA sequences and show their translation.

.. doctest::

    >>> from cogent3 import make_unaligned_seqs
    >>> seqs = [('hum', 'AAGCAGATCCAGGAAAGCAGCGAGAATGGCAGCCTGGCCGCGCGCCAGGAGAGGCAGGCCCAGGTCAACCTCACT'),
    ...         ('mus', 'AAGCAGATCCAGGAGAGCGGCGAGAGCGGCAGCCTGGCCGCGCGGCAGGAGAGGCAGGCCCAAGTCAACCTCACG'),
    ...         ('rat', 'CTGAACAAGCAGCCACTTTCAAACAAGAAA')]
    >>> unaligned_DNA = make_unaligned_seqs(seqs, moltype="dna")
    >>> print(unaligned_DNA)  # doctest: +SKIP
    >hum
    AAGCAGATCCAGGAAAGCAGCGAGAATGGCAGCCTGGCCGCGCGCCAGGAGAGGCAGGCCCAGGTCAACCTCACT
    >mus
    AAGCAGATCCAGGAGAGCGGCGAGAGCGGCAGCCTGGCCGCGCGGCAGGAGAGGCAGGCCCAAGTCAACCTCACG
    >rat
    CTGAACAAGCAGCCACTTTCAAACAAGAAA
    >>> print(unaligned_DNA.get_translation())  # doctest: +SKIP
    >hum
    KQIQESSENGSLAARQERQAQVNLT
    >mus
    KQIQESGESGSLAARQERQAQVNLT
    >rat
    LNKQPLSNKK
    <BLANKLINE>

We load an alignment of these protein sequences.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aligned_aa_seqs = [('hum', 'KQIQESSENGSLAARQERQAQVNLT'),
    ...                    ('mus', 'KQIQESGESGSLAARQERQAQVNLT'),
    ...                    ('rat', 'LNKQ------PLS---------NKK')]
    >>> aligned_aa = make_aligned_seqs(aligned_aa_seqs, moltype="protein")

We then obtain an alignment of the DNA sequences from the alignment of their translation.

.. doctest::

    >>> aligned_DNA = aligned_aa.replace_seqs(unaligned_DNA, aa_to_codon=True)
    >>> print(aligned_DNA)  # doctest: +SKIP
    >hum
    AAGCAGATCCAGGAAAGCAGCGAGAATGGCAGCCTGGCCGCGCGCCAGGAGAGGCAGGCCCAGGTCAACCTCACT
    >mus
    AAGCAGATCCAGGAGAGCGGCGAGAGCGGCAGCCTGGCCGCGCGGCAGGAGAGGCAGGCCCAAGTCAACCTCACG
    >rat
    CTGAACAAGCAG------------------CCACTTTCA---------------------------AACAAGAAA
    <BLANKLINE>

Setting the argument ``aa_to_codons=False`` is only useful when the sequences have exactly the length. One use case is to allow introducing the gaps onto another copy of the alignment where there are annotations.
