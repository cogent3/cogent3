Annotations
^^^^^^^^^^^

.. Gavin Huttley, Tom Elliot

Annotations with coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For more extensive documentation about annotations see :ref:`seq-annotations`.

Automated introduction from reading genbank files
"""""""""""""""""""""""""""""""""""""""""""""""""

We load a sample genbank file with plenty of features and grab the CDS features.

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

Customising annotation construction from reading a genbank file
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can write your own code to construct annotation objects. One reason you might do this is some genbank files do not have a ``/gene`` tag on gene related features, instead only possessing a ``/locus_tag``. For illustrating the approach we only create annotations for ``CDS`` features. We write a custom callback function that uses the ``locus_tag`` as the ``Feature`` name.

.. doctest::
    
    >>> from cogent.core.annotation import Feature
    >>> def add_annotation(seq, feature, spans):
    ...     type_ = feature['type']
    ...     if type_ != 'CDS':
    ...         return
    ...     name = feature.get('locus_tag', None)
    ...     if name and not isinstance(name, basestring):
    ...         name = ' '.join(name)
    ...     seq.addAnnotation(Feature, type_, name, spans)
    ... 
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'),
    ...          add_annotation=add_annotation)
    >>> for accession, seq in parser: # just reading one accession,sequence
    ...     break
    ...  
    >>> genes = seq.getAnnotationsMatching('CDS')
    >>> print genes
    [CDS "STM0001" at [189:255]/10020, CDS "STM0002" at [336:2799]/10020...

Creating directly on a sequence
"""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> from cogent.core.annotation import Feature
    >>> s1 = DNA.makeSequence("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      Name="seq1")
    ...
    >>> print s1[10:15] # this will be exon 1
    CCCCC
    >>> print s1[30:40] # this will be exon 2
    TTTTTAAAAA
    >>> print s1[45:48] # this will be exon 3
    CCC
    >>> s2 = DNA.makeSequence("CGAAACGTTT", Name="seq2")
    >>> s3 = DNA.makeSequence("CGAAACGTTT", Name="seq3")

Via
"""

``addAnnotation``
+++++++++++++++++

.. doctest::

    >>> from cogent import DNA
    >>> from cogent.core.annotation import Feature
    >>> s1 = DNA.makeSequence("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      Name="seq1")
    ...
    >>> exon1 = s1.addAnnotation(Feature, 'exon', 'A', [(10,15)])
    >>> exon2 = s1.addAnnotation(Feature, 'exon', 'B', [(30,40)])

``addFeature``
++++++++++++++

.. doctest::

    >>> from cogent import DNA
    >>> s1 = DNA.makeSequence("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      Name="seq1")
    ...
    >>> exon3 = s1.addFeature('exon', 'C', [(45, 48)])

*There are other annotation types.*

Adding as a series or item-wise
"""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> s2 = DNA.makeSequence("CGAAACGTTT", Name="seq2")
    >>> cpgs_series = s2.addFeature('cpgsite', 'cpg', [(0,2), (5,7)])
    >>> s3 = DNA.makeSequence("CGAAACGTTT", Name="seq3")
    >>> cpg1 = s3.addFeature('cpgsite', 'cpg', [(0,2)])
    >>> cpg2 = s3.addFeature('cpgsite', 'cpg', [(5,7)])

Taking the union of annotations
"""""""""""""""""""""""""""""""

Construct a pseudo-feature (``cds``) that's a union of other features (``exon1``, ``exon2``, ``exon3``).

.. doctest::
    
    >>> from cogent import DNA
    >>> s1 = DNA.makeSequence("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      Name="seq1")
    ...
    >>> exon1 = s1.addFeature('exon', 'A', [(10,15)])
    >>> exon2 = s1.addFeature('exon', 'B', [(30,40)])
    >>> exon3 = s1.addFeature('exon', 'C', [(45, 48)])
    >>> cds = s1.getRegionCoveringAll([exon1, exon2, exon3])

Getting annotation coordinates
""""""""""""""""""""""""""""""

These are useful for doing custom things, e.g. you could construct intron features using the below.

.. doctest::
    
    >>> cds.getCoordinates()
    [(10, 15), (30, 40), (45, 48)]

Annotations have shadows
""""""""""""""""""""""""

A shadow is a span representing everything but the annotation.

.. doctest::

    >>> not_cds = cds.getShadow()
    >>> not_cds
    region "not exon" at [0:10, 15:30, 40:45, 48:49]/49

Compare to the coordinates of the original.

.. doctest::

    >>> cds
    region "exon" at [10:15, 30:40, 45:48]/49

Adding to a sequence member of an alignment
"""""""""""""""""""""""""""""""""""""""""""

The following annotation is directly applied onto the sequence and so is in ungapped sequence coordinates.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> aln1 = LoadSeqs(data=[['x','-AAACCCCCA'],
    ...                       ['y','TTTT--TTTT']])
    >>> seq_exon = aln1.getSeq('x').addFeature('exon', 'A', [(3,8)])

Adding to an alignment
""""""""""""""""""""""

We add an annotation directly onto an alignment. In this example we add a ``Variable`` that can be displayed as a red line on a drawing. The resulting annotation (``red_data`` here) is in **alignment coordinates**!

.. doctest::

    >>> from cogent.core.annotation import Variable
    >>> red_data = aln1.addAnnotation(Variable, 'redline', 'align',
    ...              [((0,15),1),((15,30),2),((30,45),3)])
    ...

Slicing sequences and alignments by annotations
"""""""""""""""""""""""""""""""""""""""""""""""

By a feature or coordinates returns same sequence span

.. doctest::

    >>> from cogent import DNA
    >>> s1 = DNA.makeSequence("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      Name="seq1")
    ...
    >>> exon1 = s1.addFeature('exon', 'A', [(10,15)])
    >>> exon2 = s1.addFeature('exon', 'B', [(30,40)])
    >>> s1[exon1]
    DnaSequence(CCCCC)
    >>> s1[10:15]
    DnaSequence(CCCCC)

Using the annotation object ``getSlice`` method returns the same thing.

.. doctest::

    >>> s1[exon2]
    DnaSequence(TTTTTAAAAA)
    >>> exon2.getSlice()
    DnaSequence(TTTTTAAAAA)

Slicing by pseudo-feature or feature series
"""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> s1 = DNA.makeSequence("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      Name="seq1")
    ...
    >>> exon1 = s1.addFeature('exon', 'A', [(10,15)])
    >>> exon2 = s1.addFeature('exon', 'B', [(30,40)])
    >>> exon3 = s1.addFeature('exon', 'C', [(45, 48)])
    >>> cds = s1.getRegionCoveringAll([exon1, exon2, exon3])
    >>> print s1[cds]
    CCCCCTTTTTAAAAACCC
    >>> print s1[exon1, exon2, exon3]
    CCCCCTTTTTAAAAACCC

.. warning:: Slices are applied in order!

.. doctest::

    >>> print s1
    AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT
    >>> print s1[exon1, exon2, exon3]
    CCCCCTTTTTAAAAACCC
    >>> print s1[exon2]
    TTTTTAAAAA
    >>> print s1[exon3]
    CCC
    >>> print s1[exon1, exon3, exon2]
    CCCCCCCCTTTTTAAAAA

Slice series must not be overlapping
""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> s1[1:10, 9:15]
    Traceback (most recent call last):
    ValueError: Uninvertable. Overlap: 9 < 10
    >>> s1[exon1, exon1]
    Traceback (most recent call last):
    ValueError: Uninvertable. Overlap: 10 < 15

But ``getRegionCoveringAll`` resolves this, ensuring no overlaps.

.. doctest::

    >>> print s1.getRegionCoveringAll([exon3, exon3]).getSlice()
    CCC

You can slice an annotation itself
""""""""""""""""""""""""""""""""""

.. doctest::

    >>> print s1[exon2]
    TTTTTAAAAA
    >>> ex2_start = exon2[0:3]
    >>> print s1[ex2_start]
    TTT
    >>> ex2_end = exon2[-3:]
    >>> print s1[ex2_end]
    AAA

Sequence vs Alignment slicing
"""""""""""""""""""""""""""""

You can't slice an alignment using an annotation from a sequence.

.. doctest::

    >>> aln1[seq_exon]
    Traceback (most recent call last):
    ValueError: Can't map exon "A" at [3:8]/9 onto 2 x 10 text alignment: x[-AAACCCCCA], y[TTTT--TTTT] via []

Copying annotations
"""""""""""""""""""

You can copy annotations onto sequences with the same name, even if the length differs

.. doctest::

    >>> aln2 = LoadSeqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']])
    >>> seq = DNA.makeSequence('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> match_exon = seq.addFeature('exon', 'A', [(3,8)])
    >>> aln2.getSeq('x').copyAnnotations(seq)
    >>> copied = list(aln2.getAnnotationsFromSequence('x', 'exon'))
    >>> copied
    [exon "A" at [4:9]/10]

but if the feature lies outside the sequence being copied to, you get a lost span

.. doctest::

    >>> aln2 = LoadSeqs(data=[['x', '-AAAA'], ['y', 'TTTTT']])
    >>> seq = DNA.makeSequence('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> match_exon = seq.addFeature('exon', 'A', [(5,8)])
    >>> aln2.getSeq('x').copyAnnotations(seq)
    >>> copied = list(aln2.getAnnotationsFromSequence('x', 'exon'))
    >>> copied
    [exon "A" at [5:5, -4-]/5]
    >>> copied[0].getSlice()
    2 x 4 text alignment: x[----], y[----]

You can copy to a sequence with a different name, in a different alignment if the feature lies within the length

.. doctest::

    >>> # new test
    >>> aln2 = LoadSeqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']])
    >>> seq = DNA.makeSequence('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> match_exon = seq.addFeature('exon', 'A', [(5,8)])
    >>> aln2.getSeq('y').copyAnnotations(seq)
    >>> copied = list(aln2.getAnnotationsFromSequence('y', 'exon'))
    >>> copied
    [exon "A" at [7:10]/10]

If the sequence is shorter, again you get a lost span.

.. doctest::

    >>> aln2 = LoadSeqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']])
    >>> diff_len_seq = DNA.makeSequence('CCCCCCCCCCCCCCCCCCCCCCCCCCCC', 'x')
    >>> nonmatch = diff_len_seq.addFeature('repeat', 'A', [(12,14)])
    >>> aln2.getSeq('y').copyAnnotations(diff_len_seq)
    >>> copied = list(aln2.getAnnotationsFromSequence('y', 'repeat'))
    >>> copied
    [repeat "A" at [10:10, -6-]/10]

Querying
""""""""

You need to get a corresponding annotation projected into alignment coordinates via a query.

.. doctest::

    >>> aln_exon = aln1.getAnnotationsFromAnySequence('exon')
    >>> print aln1[aln_exon]
    >x
    CCCCC
    >y
    --TTT
    <BLANKLINE>

Querying produces objects only valid for their source
"""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> cpgsite2 = s2.getAnnotationsMatching('cpgsite')
    >>> print s2[cpgsite2]
    CGCG
    >>> cpgsite3 = s3.getAnnotationsMatching('cpgsite')
    >>> s2[cpgsite3]
    Traceback (most recent call last):
    ValueError: Can't map cpgsite "cpg" at [0:2]/10 onto DnaSequence(CGAAACGTTT) via []

Querying for absent annotation
""""""""""""""""""""""""""""""

You get back an empty list, and slicing with this returns an empty sequence.

.. doctest::

    >>> # this test is new
    >>> dont_exist = s2.getAnnotationsMatching('dont_exist')
    >>> dont_exist
    []
    >>> s2[dont_exist]
    DnaSequence()

Querying features that span gaps in alignments
""""""""""""""""""""""""""""""""""""""""""""""

If you query for a feature from a sequence, it's alignment coordinates may be discontinuous.

.. doctest::

    >>> aln3 = LoadSeqs(data=[['x', 'C-CCCAAAAA'], ['y', '-T----TTTT']])
    >>> exon = aln3.getSeq('x').addFeature('exon', 'ex1', [(0,4)])
    >>> print exon.getSlice()
    CCCC
    >>> aln_exons = list(aln3.getAnnotationsFromSequence('x', 'exon'))
    >>> print aln_exons
    [exon "ex1" at [0:1, 2:5]/10]
    >>> print aln3[aln_exons]
    >x
    CCCC
    >y
    ----
    <BLANKLINE>

.. note:: The ``T`` opposite the gap is missing since this approach only returns positions directly corresponding to the feature.

``asOneSpan`` unifies features with discontinuous alignment coordinates
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To get positions spanned by a feature, including gaps, use ``asOneSpan``.

.. doctest::

    >>> unified = aln_exons[0].asOneSpan()
    >>> print aln3[unified]
    >x
    C-CCC
    >y
    -T---
    <BLANKLINE>

Behaviour of annotations on nucleic acid sequences
""""""""""""""""""""""""""""""""""""""""""""""""""

Reverse complementing a sequence **does not** reverse annotations, that is they retain the reference to the frame for which they were defined.

.. doctest::

    >>> plus = DNA.makeSequence("CCCCCAAAAAAAAAATTTTTTTTTTAAAGG")
    >>> plus_rpt = plus.addFeature('blah', 'a', [(5,15), (25, 28)])
    >>> print plus[plus_rpt]
    AAAAAAAAAAAAA
    >>> minus = plus.rc()
    >>> print minus
    CCTTTAAAAAAAAAATTTTTTTTTTGGGGG
    >>> minus_rpt = minus.getAnnotationsMatching('blah')
    >>> print minus[minus_rpt]
    AAAAAAAAAAAAA

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

Masking annotated regions on alignments
"""""""""""""""""""""""""""""""""""""""

We mask exon's on an alignment.

.. doctest::
    
    >>> from cogent import LoadSeqs, DNA
    >>> aln = LoadSeqs(data=[['x', 'C-CCCAAAAAGGGAA'],
    ...                       ['y', '-T----TTTTG-GTT']], moltype=DNA)
    >>> exon = aln.getSeq('x').addFeature('exon', 'norwegian', [(0,4)])
    >>> print aln.withMaskedAnnotations('exon', mask_char='?')
    >x
    ?-???AAAAAGGGAA
    >y
    -T----TTTTG-GTT
    <BLANKLINE>

These also persist through reverse complement operations.

.. doctest::
    
    >>> rc = aln.rc()
    >>> print rc
    >x
    TTCCCTTTTTGGG-G
    >y
    AAC-CAAAA----A-
    <BLANKLINE>
    >>> print rc.withMaskedAnnotations('exon', mask_char='?')
    >x
    TTCCCTTTTT???-?
    >y
    AAC-CAAAA----A-
    <BLANKLINE>

You can take mask of the shadow
"""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import DNA
    >>> s = DNA.makeSequence('CCCCAAAAAGGGAA', 'x')
    >>> exon = s.addFeature('exon', 'norwegian', [(0,4)])
    >>> rpt = s.addFeature('repeat', 'norwegian', [(9, 12)])
    >>> rc = s.rc()
    >>> print s.withMaskedAnnotations('exon', shadow=True)
    CCCC??????????
    >>> print rc.withMaskedAnnotations('exon', shadow=True)
    ??????????GGGG
    >>> print s.withMaskedAnnotations(['exon', 'repeat'], shadow=True)
    CCCC?????GGG??
    >>> print rc.withMaskedAnnotations(['exon', 'repeat'], shadow=True)
    ??CCC?????GGGG

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

The annotation methods ``getRegionCoveringAll`` and ``getShadow`` can be used to grab all the coding sequences or non-coding sequences in a ``DnaSequence`` object.

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

Sequence features can be accessed via a containing ``Alignment``.

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
    >>> annot = aln.addAnnotation(Variable, 'redline', 'align',
    ...                          [((0,15),1),((15,30),2),((30,45),3)])
    >>> annot = aln.addAnnotation(Variable, 'blueline', 'align',
    ...                          [((0,15),1.5),((15,30),2.5),((30,45),3.5)])
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

