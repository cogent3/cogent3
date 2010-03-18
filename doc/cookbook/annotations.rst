Annotations
^^^^^^^^^^^

In its simplest usage, annotations provide a way of labeling a particular region or segment of a ``DnaSequence`` or ``Alignment`` object as containing a feature.  Different segments of the same type (cds or exon or any new feature you define), can then grouped together for analysis.  The Python slice operator ([ ], ``__getitem__``) can be used to access a particular feature.

All the features of a given type can be obtained in a list by using the method ``getAnnotationsMatching()``

.. doctest::

    >>> from cogent import DNA
    >>> s = DNA.makeSequence('ATGACCCTGTAAAAAATGTGTTAACCC',
    ...    Name='a')
    >>> cds1 = s.addFeature('cds','cds1', [(0,12)])
    >>> cds2 = s.addFeature('cds','cds2', [(15,24)])
    >>> s
    DnaSequence(ATGACCC... 27)
    >>> cds1
    cds "cds1" at [0:12]/27
    >>> type(cds1)
    <class 'cogent.core.annotation.AnnotatableFeature'>
    >>> s[cds2]
    DnaSequence(ATGTGTTAA)
    >>> all_cds = s.getAnnotationsMatching('cds')
    >>> all_cds
    [cds "cds1" at [0:12]/27, cds "cds2" at [15:24]/27]

In addition, ``getRegionCoveringAll()`` and ``getShadow()`` can be used to grab all the coding sequences or non-coding sequences in a ``DnaSequence`` object.

    >>> coding_seqs = s.getRegionCoveringAll(all_cds)
    >>> coding_seqs
    region "cds" at [0:12, 15:24]/27
    >>> type(coding_seqs)
    <class 'cogent.core.annotation._Feature'>
    >>> coding_seqs.getSlice()
    DnaSequence(ATGACCC... 21)
    >>> noncoding_seqs = coding_seqs.getShadow()
    >>> noncoding_seqs
    region "not cds" at [12:15, 24:27]/27
    >>> print noncoding_seqs.getSlice()
    AAACCC

For more extensive documentation about annotations see :ref:`seq-annotations`.

Annotations with coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*To be written.*

Automated introduction from reading genbank files
"""""""""""""""""""""""""""""""""""""""""""""""""

*To be written.*

Manipulating annotated regions
""""""""""""""""""""""""""""""

*To be written.*

Annotation display
""""""""""""""""""

*To be written.*

Introducing your own
""""""""""""""""""""

*To be written.*

Displaying them
"""""""""""""""

*To be written.*

Generic metadata
^^^^^^^^^^^^^^^^

*To be written.*

Info object
"""""""""""

*To be written.*

