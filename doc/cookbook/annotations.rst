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

    >>> from cogent3.parse.genbank import RichGenbankParser
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'))
    >>> for accession, seq in parser:
    ...     print(accession)
    ...
    AE006468
    >>> cds = seq.get_annotations_matching('CDS')
    >>> print(cds)
    [CDS "thrL" at [189:255]/10020, CDS "thrA" at ...

Customising annotation construction from reading a genbank file
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can write your own code to construct annotation objects. One reason you might do this is some genbank files do not have a ``/gene`` tag on gene related features, instead only possessing a ``/locus_tag``. For illustrating the approach we only create annotations for ``CDS`` features. We write a custom callback function that uses the ``locus_tag`` as the ``Feature`` name.

.. doctest::

    >>> from cogent3.core.annotation import Feature
    >>> def add_annotation(seq, feature, spans):
    ...     type_ = feature['type']
    ...     if type_ != 'CDS':
    ...         return
    ...     name = feature.get('locus_tag', None)
    ...     if name and not isinstance(name, str):
    ...         name = ' '.join(name)
    ...     seq.add_annotation(Feature, type_, name, spans)
    ...
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'),
    ...          add_annotation=add_annotation)
    >>> for accession, seq in parser: # just reading one accession,sequence
    ...     break
    ...
    >>> genes = seq.get_annotations_matching('CDS')
    >>> print(genes)
    [CDS "STM0001" at [189:255]/10020, CDS "STM0002" at [336:2799]/10020...

Creating directly on a sequence
"""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import DNA
    >>> from cogent3.core.annotation import Feature
    >>> s1 = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      name="seq1")
    ...
    >>> print(s1[10:15]) # this will be exon 1
    CCCCC
    >>> print(s1[30:40]) # this will be exon 2
    TTTTTAAAAA
    >>> print(s1[45:48]) # this will be exon 3
    CCC
    >>> s2 = DNA.make_seq("CGAAACGTTT", name="seq2")
    >>> s3 = DNA.make_seq("CGAAACGTTT", name="seq3")

Via
"""

``add_annotation``
++++++++++++++++++

.. doctest::

    >>> from cogent3 import DNA
    >>> from cogent3.core.annotation import Feature
    >>> s1 = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      name="seq1")
    ...
    >>> exon1 = s1.add_annotation(Feature, 'exon', 'A', [(10,15)])
    >>> exon2 = s1.add_annotation(Feature, 'exon', 'B', [(30,40)])

``add_feature``
+++++++++++++++

.. doctest::

    >>> from cogent3 import DNA
    >>> s1 = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      name="seq1")
    ...
    >>> exon3 = s1.add_feature('exon', 'C', [(45, 48)])

*There are other annotation types.*

Adding as a series or item-wise
"""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import DNA
    >>> s2 = DNA.make_seq("CGAAACGTTT", name="seq2")
    >>> cpgs_series = s2.add_feature('cpgsite', 'cpg', [(0,2), (5,7)])
    >>> s3 = DNA.make_seq("CGAAACGTTT", name="seq3")
    >>> cpg1 = s3.add_feature('cpgsite', 'cpg', [(0,2)])
    >>> cpg2 = s3.add_feature('cpgsite', 'cpg', [(5,7)])

Taking the union of annotations
"""""""""""""""""""""""""""""""

Construct a pseudo-feature (``cds``) that's a union of other features (``exon1``, ``exon2``, ``exon3``).

.. doctest::

    >>> from cogent3 import DNA
    >>> s1 = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      name="seq1")
    ...
    >>> exon1 = s1.add_feature('exon', 'A', [(10,15)])
    >>> exon2 = s1.add_feature('exon', 'B', [(30,40)])
    >>> exon3 = s1.add_feature('exon', 'C', [(45, 48)])
    >>> cds = s1.get_region_covering_all([exon1, exon2, exon3])

Getting annotation coordinates
""""""""""""""""""""""""""""""

These are useful for doing custom things, e.g. you could construct intron features using the below.

.. doctest::

    >>> cds.get_coordinates()
    [(10, 15), (30, 40), (45, 48)]

Annotations have shadows
""""""""""""""""""""""""

A shadow is a span representing everything but the annotation.

.. doctest::

    >>> not_cds = cds.get_shadow()
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

    >>> from cogent3 import make_aligned_seqs
    >>> aln1 = make_aligned_seqs(data=[['x','-AAACCCCCA'],
    ...                       ['y','TTTT--TTTT']], array_align=False)
    >>> seq_exon = aln1.get_seq('x').add_feature('exon', 'A', [(3,8)])

Adding to an alignment
""""""""""""""""""""""

We add an annotation directly onto an alignment. In this example we add a ``Variable`` that can be displayed as a red line on a drawing. The resulting annotation (``red_data`` here) is in **alignment coordinates**!

.. doctest::

    >>> from cogent3.core.annotation import Variable
    >>> red_data = aln1.add_annotation(Variable, 'redline', 'align',
    ...              [((0,15),1),((15,30),2),((30,45),3)])
    ...

Slicing sequences and alignments by annotations
"""""""""""""""""""""""""""""""""""""""""""""""

By a feature or coordinates returns same sequence span

.. doctest::

    >>> from cogent3 import DNA
    >>> s1 = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      name="seq1")
    ...
    >>> exon1 = s1.add_feature('exon', 'A', [(10,15)])
    >>> exon2 = s1.add_feature('exon', 'B', [(30,40)])
    >>> s1[exon1]
    DnaSequence(CCCCC)
    >>> s1[10:15]
    DnaSequence(CCCCC)

Using the annotation object ``get_slice`` method returns the same thing.

.. doctest::

    >>> s1[exon2]
    DnaSequence(TTTTTAAAAA)
    >>> exon2.get_slice()
    DnaSequence(TTTTTAAAAA)

Slicing by pseudo-feature or feature series
"""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import DNA
    >>> s1 = DNA.make_seq("AAGAAGAAGACCCCCAAAAAAAAAA"\
    ...                      "TTTTTTTTTTAAAAAGGGAACCCT",
    ...                      name="seq1")
    ...
    >>> exon1 = s1.add_feature('exon', 'A', [(10,15)])
    >>> exon2 = s1.add_feature('exon', 'B', [(30,40)])
    >>> exon3 = s1.add_feature('exon', 'C', [(45, 48)])
    >>> cds = s1.get_region_covering_all([exon1, exon2, exon3])
    >>> print(s1[cds])
    CCCCCTTTTTAAAAACCC
    >>> print(s1[exon1, exon2, exon3])
    CCCCCTTTTTAAAAACCC

.. warning:: Slices are applied in order!

.. doctest::

    >>> print(s1)
    AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT
    >>> print(s1[exon1, exon2, exon3])
    CCCCCTTTTTAAAAACCC
    >>> print(s1[exon2])
    TTTTTAAAAA
    >>> print(s1[exon3])
    CCC
    >>> print(s1[exon1, exon3, exon2])
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

But ``get_region_covering_all`` resolves this, ensuring no overlaps.

.. doctest::

    >>> print(s1.get_region_covering_all([exon3, exon3]).get_slice())
    CCC

You can slice an annotation itself
""""""""""""""""""""""""""""""""""

.. doctest::

    >>> print(s1[exon2])
    TTTTTAAAAA
    >>> ex2_start = exon2[0:3]
    >>> print(s1[ex2_start])
    TTT
    >>> ex2_end = exon2[-3:]
    >>> print(s1[ex2_end])
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

    >>> aln2 = make_aligned_seqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']],
    ...                 array_align=False)
    >>> seq = DNA.make_seq('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> match_exon = seq.add_feature('exon', 'A', [(3,8)])
    >>> aln2.get_seq('x').copy_annotations(seq)
    >>> copied = list(aln2.get_annotations_from_seq('x', 'exon'))
    >>> copied
    [exon "A" at [4:9]/10]

but if the feature lies outside the sequence being copied to, you get a lost span

.. doctest::

    >>> aln2 = make_aligned_seqs(data=[['x', '-AAAA'], ['y', 'TTTTT']], array_align=False)
    >>> seq = DNA.make_seq('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> match_exon = seq.add_feature('exon', 'A', [(5,8)])
    >>> aln2.get_seq('x').copy_annotations(seq)
    >>> copied = list(aln2.get_annotations_from_seq('x', 'exon'))
    >>> copied
    [exon "A" at [5:5, -4-]/5]
    >>> copied[0].get_slice()
    2 x 4 text alignment: x[----], y[----]

You can copy to a sequence with a different name, in a different alignment if the feature lies within the length

.. doctest::

    >>> # new test
    >>> aln2 = make_aligned_seqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']],
    ...                  array_align=False)
    >>> seq = DNA.make_seq('CCCCCCCCCCCCCCCCCCCC', 'x')
    >>> match_exon = seq.add_feature('exon', 'A', [(5,8)])
    >>> aln2.get_seq('y').copy_annotations(seq)
    >>> copied = list(aln2.get_annotations_from_seq('y', 'exon'))
    >>> copied
    [exon "A" at [7:10]/10]

If the sequence is shorter, again you get a lost span.

.. doctest::

    >>> aln2 = make_aligned_seqs(data=[['x', '-AAAAAAAAA'], ['y', 'TTTT--TTTT']],
    ...                 array_align=False)
    >>> diff_len_seq = DNA.make_seq('CCCCCCCCCCCCCCCCCCCCCCCCCCCC', 'x')
    >>> nonmatch = diff_len_seq.add_feature('repeat', 'A', [(12,14)])
    >>> aln2.get_seq('y').copy_annotations(diff_len_seq)
    >>> copied = list(aln2.get_annotations_from_seq('y', 'repeat'))
    >>> copied
    [repeat "A" at [10:10, -6-]/10]

Querying
""""""""

You need to get a corresponding annotation projected into alignment coordinates via a query.

.. doctest::

    >>> aln_exon = aln1.get_annotations_from_any_seq('exon')
    >>> print(aln1[aln_exon])
    >x
    CCCCC
    >y
    --TTT
    <BLANKLINE>

Querying produces objects only valid for their source
"""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> cpgsite2 = s2.get_annotations_matching('cpgsite')
    >>> print(s2[cpgsite2])
    CGCG
    >>> cpgsite3 = s3.get_annotations_matching('cpgsite')
    >>> s2[cpgsite3]
    Traceback (most recent call last):
    ValueError: Can't map cpgsite "cpg" at [0:2]/10 onto DnaSequence(CGAAACGTTT) via []

Querying for absent annotation
""""""""""""""""""""""""""""""

You get back an empty list, and slicing with this returns an empty sequence.

.. doctest::

    >>> # this test is new
    >>> dont_exist = s2.get_annotations_matching('dont_exist')
    >>> dont_exist
    []
    >>> s2[dont_exist]
    DnaSequence()

Querying features that span gaps in alignments
""""""""""""""""""""""""""""""""""""""""""""""

If you query for a feature from a sequence, it's alignment coordinates may be discontinuous.

.. doctest::

    >>> aln3 = make_aligned_seqs(data=[['x', 'C-CCCAAAAA'], ['y', '-T----TTTT']],
    ...                 array_align=False)
    >>> exon = aln3.get_seq('x').add_feature('exon', 'ex1', [(0,4)])
    >>> print(exon.get_slice())
    CCCC
    >>> aln_exons = list(aln3.get_annotations_from_seq('x', 'exon'))
    >>> print(aln_exons)
    [exon "ex1" at [0:1, 2:5]/10]
    >>> print(aln3[aln_exons])
    >x
    CCCC
    >y
    ----
    <BLANKLINE>

.. note:: The ``T`` opposite the gap is missing since this approach only returns positions directly corresponding to the feature.

``as_one_span`` unifies features with discontinuous alignment coordinates
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To get positions spanned by a feature, including gaps, use ``as_one_span``.

.. doctest::

    >>> unified = aln_exons[0].as_one_span()
    >>> print(aln3[unified])
    >x
    C-CCC
    >y
    -T---
    <BLANKLINE>

Behaviour of annotations on nucleic acid sequences
""""""""""""""""""""""""""""""""""""""""""""""""""

Reverse complementing a sequence **does not** reverse annotations, that is they retain the reference to the frame for which they were defined.

.. doctest::

    >>> plus = DNA.make_seq("CCCCCAAAAAAAAAATTTTTTTTTTAAAGG")
    >>> plus_rpt = plus.add_feature('blah', 'a', [(5,15), (25, 28)])
    >>> print(plus[plus_rpt])
    AAAAAAAAAAAAA
    >>> minus = plus.rc()
    >>> print(minus)
    CCTTTAAAAAAAAAATTTTTTTTTTGGGGG
    >>> minus_rpt = minus.get_annotations_matching('blah')
    >>> print(minus[minus_rpt])
    AAAAAAAAAAAAA

Masking annotated regions
"""""""""""""""""""""""""

We mask the CDS regions.

.. doctest::

    >>> from cogent3.parse.genbank import RichGenbankParser
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'))
    >>> seq = [seq for accession, seq in parser][0]
    >>> no_cds = seq.with_masked_annotations('CDS')
    >>> print(no_cds[150:400])
    CAAGACAGACAAATAAAAATGACAGAGTACACAACATCC?????????...

The above sequence could then have positions filtered so no position with the ambiguous character '?' was present.

Masking annotated regions on alignments
"""""""""""""""""""""""""""""""""""""""

We mask exon's on an alignment.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data=[['x', 'C-CCCAAAAAGGGAA'],
    ...                      ['y', '-T----TTTTG-GTT']],
    ...                moltype="dna", array_align=False)
    >>> exon = aln.get_seq('x').add_feature('exon', 'norwegian', [(0,4)])
    >>> print(aln.with_masked_annotations('exon', mask_char='?'))
    >x
    ?-???AAAAAGGGAA
    >y
    -T----TTTTG-GTT
    <BLANKLINE>

These also persist through reverse complement operations.

.. doctest::

    >>> rc = aln.rc()
    >>> print(rc)
    >x
    TTCCCTTTTTGGG-G
    >y
    AAC-CAAAA----A-
    <BLANKLINE>
    >>> print(rc.with_masked_annotations('exon', mask_char='?'))
    >x
    TTCCCTTTTT???-?
    >y
    AAC-CAAAA----A-
    <BLANKLINE>

You can take mask of the shadow
"""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import DNA
    >>> s = DNA.make_seq('CCCCAAAAAGGGAA', 'x')
    >>> exon = s.add_feature('exon', 'norwegian', [(0,4)])
    >>> rpt = s.add_feature('repeat', 'norwegian', [(9, 12)])
    >>> rc = s.rc()
    >>> print(s.with_masked_annotations('exon', shadow=True))
    CCCC??????????
    >>> print(rc.with_masked_annotations('exon', shadow=True))
    ??????????GGGG
    >>> print(s.with_masked_annotations(['exon', 'repeat'], shadow=True))
    CCCC?????GGG??
    >>> print(rc.with_masked_annotations(['exon', 'repeat'], shadow=True))
    ??CCC?????GGGG

What features of a certain type are available?
""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import DNA
    >>> s = DNA.make_seq('ATGACCCTGTAAAAAATGTGTTAACCC',
    ...    name='a')
    >>> cds1 = s.add_feature('cds','cds1', [(0,12)])
    >>> cds2 = s.add_feature('cds','cds2', [(15,24)])
    >>> all_cds = s.get_annotations_matching('cds')
    >>> all_cds
    [cds "cds1" at [0:12]/27, cds "cds2" at [15:24]/27]

Getting all features of a type, or everything but that type
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The annotation methods ``get_region_covering_all`` and ``get_shadow`` can be used to grab all the coding sequences or non-coding sequences in a ``DnaSequence`` object.

.. doctest::

    >>> from cogent3.parse.genbank import RichGenbankParser
    >>> parser = RichGenbankParser(open('data/ST_genome_part.gb'))
    >>> seq = [seq for accession, seq in parser][0]
    >>> all_cds = seq.get_annotations_matching('CDS')
    >>> coding_seqs = seq.get_region_covering_all(all_cds)
    >>> coding_seqs
    region "CDS" at [189:255, 336:2799, 2800:3730, 3733...
    >>> coding_seqs.get_slice()
    DnaSequence(ATGAACC... 9063)
    >>> noncoding_seqs = coding_seqs.get_shadow()
    >>> noncoding_seqs
    region "not CDS" at [0:189, 255:336, 2799:2800, ...
    >>> noncoding_seqs.get_slice()
    DnaSequence(AGAGATT... 957)

Getting sequence features when you have an alignment object
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sequence features can be accessed via a containing ``Alignment``.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data=[['x','-AAAAAAAAA'], ['y','TTTT--TTTT']],
    ...                array_align=False)
    >>> print(aln)
    >x
    -AAAAAAAAA
    >y
    TTTT--TTTT
    <BLANKLINE>
    >>> exon = aln.get_seq('x').add_feature('exon', '1', [(3,8)])
    >>> aln_exons = aln.get_annotations_from_seq('x', 'exon')
    >>> aln_exons = aln.get_annotations_from_any_seq('exon')
    >>> aln_exons
    [exon "1" at [4:9]/10]

Annotation display on sequences
"""""""""""""""""""""""""""""""

We can display annotations on sequences, writing to file.

.. note:: This requires `matplotlib <http://matplotlib.org/>`_ be installed.

We first make a sequence and add some annotations.

.. doctest::

    >>> from cogent3 import DNA
    >>> seq = DNA.make_seq('aaaccggttt' * 10)
    >>> v = seq.add_feature('exon', 'exon', [(20,35)])
    >>> v = seq.add_feature('repeat_unit', 'repeat_unit', [(39,49)])
    >>> v = seq.add_feature('repeat_unit', 'rep2', [(49,60)])

.. todo document info attribute
.. todo document how to visualise annotations

.. following cleans up files

.. doctest::
    :hide:

    >>> from cogent3.util.misc import remove_files
    >>> remove_files(['annotated_%d.png' % i for i in range(1,4)],
    ...               error_on_missing=False)
