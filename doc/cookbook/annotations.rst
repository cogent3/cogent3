Annotations
^^^^^^^^^^^

Annotations with coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For more extensive documentation about annotations see :ref:`seq-annotations`.

Automated introduction from reading genbank files
"""""""""""""""""""""""""""""""""""""""""""""""""

We load a sample genbank file with plenty of features and grab those corresponding to CDS.

.. doctest::

    >>> from cogent.parse.genbank import RichGenbankParser
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'))
    >>> for accession, seq in parser:
    ...     print accession
    ...
    AE006468
    >>> cds = seq.getAnnotationsMatching('CDS')
    >>> print cds
    [CDS "thrL" at [189:255]/10020, CDS "thrA" at ...

Masking annotated regions
"""""""""""""""""""""""""

We mask the CDS regions.

.. doctest::

    >>> from cogent.parse.genbank import RichGenbankParser
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'))
    >>> seq = [seq for accession, seq in parser][0]
    >>> no_cds = seq.withMaskedAnnotations('CDS')
    >>> print no_cds[150:400]
    CAAGACAGACAAATAAAAATGACAGAGTACACAACATCC?????????...

The above sequence could then have positions filtered so no position with the ambiguous character '?' was present.

.. note:: the same method exists on ``Alignment`` objects.

What features of a certain type are available?
""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> s = DNA.makeSequence('ATGACCCTGTAAAAAATGTGTTAACCC',
    ...    Name='a')
    >>> cds1 = s.addFeature('cds','cds1', [(0,12)])
    >>> cds2 = s.addFeature('cds','cds2', [(15,24)])
    >>> all_cds = s.getAnnotationsMatching('cds')
    >>> all_cds
    [cds "cds1" at [0:12]/27, cds "cds2" at [15:24]/27]

Getting all features of a type, or everything but that type
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The annotation methods ``getRegionCoveringAll()`` and ``getShadow()`` can be used to grab all the coding sequences or non-coding sequences in a ``DnaSequence`` object.

.. doctest::

    >>> from cogent.parse.genbank import RichGenbankParser
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'))
    >>> seq = [seq for accession, seq in parser][0]
    >>> all_cds = seq.getAnnotationsMatching('CDS')
    >>> coding_seqs = seq.getRegionCoveringAll(all_cds)
    >>> coding_seqs
    region "CDS" at [189:255, 336:2799, 2800:3730, 3733...
    >>> coding_seqs.getSlice()
    DnaSequence(ATGAACC... 9063)
    >>> noncoding_seqs = coding_seqs.getShadow()
    >>> noncoding_seqs
    region "not CDS" at [0:189, 255:336, 2799:2800, ...
    >>> noncoding_seqs.getSlice()
    DnaSequence(AGAGATT... 957)

Getting sequence features when you have an alignment object
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sequence features can be accessed via a containing ``Alignment``:

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln = LoadSeqs(data=[['x','-AAAAAAAAA'], ['y','TTTT--TTTT']])
    >>> print aln
    >x
    -AAAAAAAAA
    >y
    TTTT--TTTT
    <BLANKLINE>
    >>> exon = aln.getSeq('x').addFeature('exon', '1', [(3,8)])
    >>> aln_exons = aln.getAnnotationsFromSequence('x', 'exon')
    >>> aln_exons = aln.getAnnotationsFromAnySequence('exon')
    >>> aln_exons
    [exon "1" at [4:9]/10]

Introducing your own annotations
""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> s = DNA.makeSequence('ATGACCCTGTAAAAAATGTGTTAACCC',
    ...    Name='a')
    >>> cds = s.addFeature('cds','cds1', [(0,12)])
    >>> cds
    cds "cds1" at [0:12]/27

Annotation display on sequences
"""""""""""""""""""""""""""""""

We can display annotations on sequences, writing to file.

.. note:: This requires `matplotlib <http://matplotlib.sourceforge.net>`_ be installed.

We first make a sequence and add some annotations.

.. doctest::

    >>> from cogent import DNA
    >>> seq = DNA.makeSequence('aaaccggttt' * 10)
    >>> v = seq.addFeature('exon', 'exon', [(20,35)])
    >>> v = seq.addFeature('repeat_unit', 'repeat_unit', [(39,49)])
    >>> v = seq.addFeature('repeat_unit', 'rep2', [(49,60)])

We then make a ``Display`` instance and write to file. This will use standard feature policy for colouring and shape of feature types.

.. doctest::

    >>> from cogent.draw.linear import Display
    >>> seq_display = Display(seq, colour_sequences=True)
    >>> fig = seq_display.makeFigure()
    >>> fig.savefig('annotated_1.png')

Annotation display on alignments
""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA, LoadSeqs
    >>> from cogent.core.annotation import Variable
    >>> from cogent.draw.linear import Display
    >>> aln = LoadSeqs('data/primate_cdx2_promoter.fasta', moltype=DNA)[:150]
    >>> annot = aln.addAnnotation(Variable, 'redline', 'align', [((0,15),1),((15,30),2),((30,45),3)])
    >>> annot = aln.addAnnotation(Variable, 'blueline', 'align', [((0,15),1.5),((15,30),2.5),((30,45),3.5)])
    >>> align_display = Display(aln, colour_sequences=True)
    >>> fig = align_display.makeFigure(width=25, left=1, right=1)
    >>> fig.savefig('annotated_2.png')

Annotation display of a custom variable
"""""""""""""""""""""""""""""""""""""""

We just show a series of spans.

.. doctest::

    >>> from cogent import DNA
    >>> from cogent.draw.linear import Display
    >>> from cogent.core.annotation import Variable
    >>> seq = DNA.makeSequence('aaaccggttt' * 10)
    >>> annot = seq.addAnnotation(Variable, 'redline', 'align',
    ...     [((0,15),1),((15,30),2),((30,45),3)])
    ...
    >>> seq_display = Display(seq, colour_sequences=True)
    >>> fig = seq_display.makeFigure()
    >>> fig.savefig('annotated_3.png')

Generic metadata
^^^^^^^^^^^^^^^^

*To be written.*

Info object
"""""""""""

*To be written.*

.. following cleans up files

.. doctest::
    :hide:

    >>> from cogent.util.misc import remove_files
    >>> remove_files(['annotated_%d.png' % i for i in range(1,4)],
    ...               error_on_missing=False)

