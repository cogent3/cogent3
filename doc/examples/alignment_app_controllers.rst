Using alignment application controllers to align unaligned sequences
====================================================================

.. sectionauthor:: Daniel McDonald

This document provides examples of how to align sequences using the alignment application controllers. Each alignment application controller module provides the support method ``align_unaligned_seqs``. This method takes as input a ``SequenceCollection`` object or a dict mapping sequence ids to sequences, the ``MolType`` of the sequences, and an option dict containing specific parameter settings. As output, the method returns an ``Alignment`` object.

First, lets import all of the ``align_unaligned_seqs`` methods:

.. doctest::
    
    >>> from cogent.app.clustalw import align_unaligned_seqs as clustalw_align_unaligned_seqs
    >>> from cogent.app.muscle import align_unaligned_seqs as muscle_align_unaligned_seqs
    >>> from cogent.app.mafft import align_unaligned_seqs as mafft_align_unaligned_seqs

Next, we'll load our test data. We will be using DNA sequences for this example:

.. doctest::

    >>> from cogent.core.moltype import DNA
    >>> from cogent import LoadSeqs
    >>> unaligned_seqs = LoadSeqs(filename='data/test2.fasta', aligned=False)

Lets align some sequences using default parameters!

.. note:: Output is truncated for document formatting

.. doctest::
    
    >>> clustalw_aln = clustalw_align_unaligned_seqs(unaligned_seqs, DNA) 
    >>> muscle_aln = muscle_align_unaligned_seqs(unaligned_seqs, DNA) 
    >>> mafft_aln = mafft_align_unaligned_seqs(unaligned_seqs, DNA) 
    >>> clustalw_aln
    5 x 60 dna alignment: NineBande[------CGCCA...], Mouse[GCAGTGAGCCA...], ...
    >>> muscle_aln
    5 x 60 dna alignment: NineBande[------CGCCA...], Mouse[GCAGTGAGCCA...], ...
    >>> mafft_aln
    5 x 60 dna alignment: NineBande[------CGCCA...], Mouse[GCAGTGAGCCA...], ...

To change specific parameters, simply specify the parameters in a dict and pass it in:

.. note:: Output is truncated for document formatting

.. doctest::
    
    >>> clustalw_params = {'-gapopen':-3, '-quicktree':True}
    >>> clustalw_aln = clustalw_align_unaligned_seqs(unaligned_seqs, DNA, params=clustalw_params)
    >>> clustalw_aln
    5 x 60 dna alignment: NineBande[------CGCCA...], Mouse[GCAGTGAGCCA...], ...
