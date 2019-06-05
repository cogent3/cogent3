Drawing a dotplot
=================

.. sectionauthor:: Gavin Huttley

.. doctest::

    >>> from cogent3 import LoadSeqs, DNA
    >>> from cogent3.core import annotation
    >>> from cogent3.draw import dotplot

Load the alignment for illustrative purposes, I'll make one sequence a different length than the other and introduce a custom sequence annotation for a miscellaneous feature. Normally, those annotations would be on the unaligned sequences.

.. doctest::

    >>> aln = LoadSeqs("data/test.paml", moltype=DNA, array_align=False)
    >>> feature = aln.add_annotation(annotation.Feature, "misc_feature",
    ...                             "pprobs", [(38, 55)])
    >>> seq1 = aln.get_seq('NineBande')[10:-3]
    >>> seq2 = aln.get_seq('DogFaced')

Write out the dotplot as a pdf file in the current directory note that seq1 will be the x-axis, and seq2 the y-axis.

.. doctest::

    >>> dp = dotplot.Dotplot(seq1,seq2)
    >>> filename = 'dotplot_example.pdf'
    >>> dp.write_pdf(filename)

.. clean up

.. doctest::
    :hide:
    
    >>> import os
    >>> os.remove('dotplot_example.pdf')
