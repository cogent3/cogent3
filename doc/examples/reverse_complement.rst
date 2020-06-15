Getting the reverse complement
==============================

.. sectionauthor:: Gavin Huttley

This is a property of DNA, and hence alignments need to be created with the appropriate ``MolType``. In the following example, the alignment is truncated to just 50 bases for the sake of simplifying the presentation.

.. jupyter-execute::
    :linenos:

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs("data/long_testseqs.fasta", moltype="dna")[:50]

The original alignment looks like this.

.. jupyter-execute::
    :linenos:

    print(aln)

We do reverse complement very simply.

.. jupyter-execute::
    :linenos:

    naln = aln.rc()

The reverse complemented alignment looks like this.

.. jupyter-execute::
    :linenos:

    print(naln)
