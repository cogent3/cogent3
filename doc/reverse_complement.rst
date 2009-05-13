Getting the reverse complement
==============================

This is a property of DNA, and hence alignments need to be created with the appropriate ``MolType``. In the following example, the alignment is truncated to just 100 bases for the sake of simplifying the presentation.

.. pycode::

    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs("data/long_testseqs.fasta", moltype=DNA)[:100]

The original alignment looks like this.

.. pycode::

    >>> print aln
    >FlyingFox
    TGTGGCACAAATGCTCATGCCAGCTCTTTACAGCATGAGAAC---AGTTTATTATACACTAAAGACAGAATGAATGTAGAAAAGACTGACTTCTGTAATA
    >DogFaced
    TGTGGCACAAATACTCATGCCAACTCATTACAGCATGAGAACAGCAGTTTATTATACACTAAAGACAGAATGAATGTAGAAAAGACTGACTTCTGTAATA

We do reverse complement very simply.

.. pycode::

    >>> naln = aln.rc()

The reverse complemented alignment looks like this.

.. pycode::

    >>> print naln
    >FlyingFox
    TATTACAGAAGTCAGTCTTTTCTACATTCATTCTGTCTTTAGTGTATAATAAACT---GTTCTCATGCTGTAAAGAGCTGGCATGAGCATTTGTGCCACA
    >DogFaced
    TATTACAGAAGTCAGTCTTTTCTACATTCATTCTGTCTTTAGTGTATAATAAACTGCTGTTCTCATGCTGTAATGAGTTGGCATGAGTATTTGTGCCACA
    <BLANKLINE>
