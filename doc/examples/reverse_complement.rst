Getting the reverse complement
==============================

This is a property of DNA, and hence alignments need to be created with the appropriate ``MolType``. In the following example, the alignment is truncated to just 50 bases for the sake of simplifying the presentation.

.. doctest::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs("data/long_testseqs.fasta", moltype=DNA)[:50]

The original alignment looks like this.

.. doctest::

    >>> print aln
    >FlyingFox
    TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAAC---AGTTT
    >DogFaced
    TGTGGCACAAATACTCATGCCAACTCATTACAGCATGAGAACAGCAGTTT
    <BLANKLINE>

We do reverse complement very simply.

.. doctest::

    >>> naln = aln.rc()

The reverse complemented alignment looks like this.

.. doctest::

    >>> print naln
    >FlyingFox
    AAACT---GTTCTCATGCTGTAAAGAGCTGGCATGAGCATTTGTGCCACA
    >DogFaced
    AAACTGCTGTTCTCATGCTGTAATGAGTTGGCATGAGTATTTGTGCCACA
    <BLANKLINE>
