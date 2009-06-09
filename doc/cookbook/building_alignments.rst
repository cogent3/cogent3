*******************
Building alignments
*******************

.. sectionauthor:: Kristian Rother, Patrick Yannul, Gavin Huttley

Using the cogent aligners
=========================

Running a pairwise Needleman-Wunsch-Alignment
---------------------------------------------

.. doctest::
    
    >>> from cogent.align.algorithm import nw_align
    >>> seq1 = 'AKSAMITNY'
    >>> seq2 = 'AKHSAMMIT'
    >>> print nw_align(seq1,seq2)
    ('AK-SAM-ITNY', 'AKHSAMMIT--')

3rd-party apps
==============

Converting gaps from aa-seq alignment to nuc seq alignment
==========================================================

