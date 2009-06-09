Drawing a dotplot
=================

.. sectionauthor:: Gavin Huttley

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.core import annotation
    >>> from cogent.draw import dotplot

Load the alignment for illustrative purposes, I'll make one sequence a different length than the other and introduce a custom sequence annotation for a miscellaneous feature. Normally, those annotations would be on the unaligned sequences.

.. doctest::

    >>> aln = LoadSeqs("data/test.paml")
    >>> feature = aln.addAnnotation(annotation.Feature, "misc_feature",
    ...                             "pprobs", [(38, 55)])
    >>> seq1 = aln.getSeq('NineBande')[10:-3]
    >>> seq2 = aln.getSeq('DogFaced')

Write out the dotplot as a pdf file in the current directory note that seq1 will be the x-axis, and seq2 the y-axis.

.. doctest::

    >>> dp = dotplot.Display2D(seq1,seq2)
    >>> filename = 'dotplot_example.pdf'
    >>> dp.drawToPDF(filename)

.. clean up

.. doctest::
    :hide:
    
    >>> import os
    >>> os.remove('dotplot_example.pdf')
