.. jupyter-execute::
    :hide-code:

    import set_working_directory

Annotations
^^^^^^^^^^^

.. Gavin Huttley, Tom Elliot

Annotations with coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For more extensive documentation about annotations see :ref:`seq-annotations`.

Automated introduction from reading genbank files
"""""""""""""""""""""""""""""""""""""""""""""""""

We load a sample genbank file with plenty of features

.. jupyter-execute::

    from cogent3 import load_seq

    seq = load_seq("data/ST_genome_part.gb", moltype="dna")
    seq

.. jupyter-execute::

    seq.annotation_db

and grab the CDS features.

.. jupyter-execute::

    cds = list(seq.get_features(biotype="CDS"))
    cds[:3]

Creating directly on a sequence
"""""""""""""""""""""""""""""""

Create some sequences to use first.

.. jupyter-execute::

    from cogent3 import make_seq

    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT", name="seq1", moltype="dna"
    )
    s2 = make_seq("CGAAACGTTT", name="seq2", moltype="dna")
    s3 = make_seq("CGAAACGTTT", name="seq3", moltype="dna")

These correspond to the different exons

.. jupyter-execute::

    s1[10:15]  # this will be exon 1

.. jupyter-execute::

    s1[30:40]  # this will be exon 2

.. jupyter-execute::

    s1[45:48]  # this will be exon 3

Add features via
""""""""""""""""

An annotation db
++++++++++++++++

Directly into an annotation db and assigning it to the sequence attribute.

.. jupyter-execute::

    from cogent3 import make_seq
    from cogent3.core.annotation_db import BasicAnnotationDb

    db = BasicAnnotationDb()

    db.add_feature(seqid="seq1", biotype="exon", name="C", spans=[(45, 48)])
    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT", name="seq1", moltype="dna"
    )
    s1.annotation_db = db

``add_feature``
+++++++++++++++

Using the sequence method.

.. jupyter-execute::

    from cogent3 import make_seq

    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT", name="seq1", moltype="dna"
    )
    exon1 = s1.add_feature(biotype="exon", name="A", spans=[(10, 15)])
    exon2 = s1.add_feature(biotype="exon", name="B", spans=[(30, 40)])
    exon3 = s1.add_feature(biotype="exon", name="C", spans=[(45, 48)])

Adding as a series
""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    s2 = make_seq("CGAAACGTTT", name="seq2", moltype="dna")
    cpgs_series = s2.add_feature(biotype="cpgsite", name="cpg", spans=[(0, 2), (5, 7)])
    cpgs_series

Adding item-wise
""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    s3 = make_seq("CGAAACGTTT", name="seq3", moltype="dna")
    cpg1 = s3.add_feature(biotype="cpgsite", name="cpg", spans=[(0, 2)])
    cpg2 = s3.add_feature(biotype="cpgsite", name="cpg", spans=[(5, 7)])
    cpg1, cpg2

Taking the union of annotations
"""""""""""""""""""""""""""""""

Construct a pseudo-feature (``cds``) that's a union of other features (``exon1``, ``exon2``, ``exon3``).

.. jupyter-execute::

    from cogent3 import make_seq

    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT",
        name="seq1",
        moltype="dna",
    )
    exon1 = s1.add_feature(biotype="exon", name="A", spans=[(10, 15)])
    exon2 = s1.add_feature(biotype="exon", name="B", spans=[(30, 40)])
    exon3 = s1.add_feature(biotype="exon", name="C", spans=[(45, 48)])
    cds = exon1.union([exon2, exon3])
    cds

Getting annotation coordinates
""""""""""""""""""""""""""""""

These are useful for doing custom things, e.g. you could construct intron features using the below.

.. jupyter-execute::

    cds.get_coordinates()

Annotations have shadows
""""""""""""""""""""""""

A shadow is a span representing everything but the annotation.

.. jupyter-execute::

    not_cds = cds.shadow()
    not_cds

Compare to the coordinates of the original.

.. jupyter-execute::

    cds

Adding to a sequence member of an alignment
"""""""""""""""""""""""""""""""""""""""""""

The following annotation is for the sequence.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln1 = make_aligned_seqs(
        data=[["x", "-AAACCCCCA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    aln1.add_feature(
        seqid="x", biotype="exon", name="A", spans=[(3, 8)], on_alignment=False
    )

Adding to an alignment
""""""""""""""""""""""

We add an annotation directly onto an alignment. The resulting annotation (``shared`` here) is in **alignment coordinates**!

.. note:: issues to document, that on the db we must define ``seqid is None``. If on the alignment, the coordinates must lie within the alignment.

.. jupyter-execute::

    aln1.annotation_db.add_feature(
        seqid=None,
        biotype="shared",
        name="demo",
        spans=[(0, 15), (15, 30), (30, 45)],
        on_alignment=True,
    )

Slicing sequences and alignments by annotations
"""""""""""""""""""""""""""""""""""""""""""""""

By a feature or coordinates returns same sequence span

.. jupyter-execute::

    from cogent3 import make_seq

    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT",
        name="seq1",
        moltype="dna",
    )
    exon1 = s1.add_feature(biotype="exon", name="A", spans=[(10, 15)])
    exon2 = s1.add_feature(biotype="exon", name="B", spans=[(30, 40)])
    s1[exon1]

.. jupyter-execute::

    s1[10:15]

Using the annotation object ``get_slice`` method returns the same thing.

.. jupyter-execute::

    exon1.get_slice()

Slicing by pseudo-feature or feature series
"""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT",
        name="seq1",
        moltype="dna",
    )
    exon1 = s1.add_feature(biotype="exon", name="A", spans=[(10, 15)])
    exon2 = s1.add_feature(biotype="exon", name="B", spans=[(30, 40)])
    exon3 = s1.add_feature(biotype="exon", name="C", spans=[(45, 48)])
    cds = exon1.union([exon2, exon3])
    s1[cds]

Copying annotations
"""""""""""""""""""

You can copy annotations onto sequences with the same name. Note that the db instance bound to alignment and its member sequences is the same.

.. jupyter-execute::

    aln2 = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]],
        array_align=False,
        moltype="dna",
    )
    x, y = aln2.get_seq("x"), aln2.get_seq("y")
    x.annotation_db is y.annotation_db is aln2.annotation_db

.. warning:: Despite this, it is possible for the attributes to get out-of-sync. So, any copy annotations should be done using ``alignment.copy_annotations()``, **not** ``alignment.get_seq("x").copy_annotations()``.

.. jupyter-execute::

    seq = make_seq("CCCCCCCCCCCCCCCCCCCC", name="x", moltype="dna")
    match_exon = seq.add_feature(biotype="exon", name="A", spans=[(3, 8)])
    aln2.copy_annotations(seq.annotation_db)
    aln2.annotation_db

but if the feature lies outside the sequence being copied to, you get a lost span

.. jupyter-execute::

    copied = list(aln2.get_features(seqid="x", biotype="exon"))
    copied

Querying on alignment gives an alignment feature
""""""""""""""""""""""""""""""""""""""""""""""""

Sequence coordinates are projected into alignment coordinates.

.. jupyter-execute::

    aln_exon = list(aln1.get_features(biotype="exon"))[0]
    aln_exon.get_slice()

Querying produces objects only valid for their source
"""""""""""""""""""""""""""""""""""""""""""""""""""""

To get a sequence annotation via the alignment, we must get the sequence itself.

.. jupyter-execute::

    x = aln1.get_seq("x")
    x_exon = list(x.get_features(biotype="exon"))[0]
    x_exon

As the representation indicates, the ``Feature`` is bound to the sequence and not the alignment. So trying to use it on the alignment raises an exception.

.. jupyter-execute::
    :raises: ValueError

    aln1[x_exon]

Querying for absent annotation
""""""""""""""""""""""""""""""

You get back an empty list, and slicing with this returns an empty sequence.

.. jupyter-execute::

    # this test is new
    dont_exist = list(s2.get_features(biotype="dont_exist"))
    dont_exist

Querying features that span gaps in alignments
""""""""""""""""""""""""""""""""""""""""""""""

If you query for a feature from a sequence (i.e. the feature is in sequence coordinates), it's alignment coordinates may be discontinuous. leading to omission of data from other sequences

.. jupyter-execute::

    aln3 = make_aligned_seqs(
        data=[["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]],
        array_align=False,
        moltype="dna",
    )
    exon = aln3.add_feature(seqid="x", biotype="exon", name="ex1", spans=[(0, 4)], on_alignment=False)
    exon.get_slice()

.. jupyter-execute::

    aln_exons = list(aln3.get_features(seqid="x", biotype="exon"))[0]
    aln_exons

.. note:: The ``T`` opposite the gap is missing since this approach only returns positions directly corresponding to the feature.

To include the gaps, use the ``allow_gaps`` argument

.. jupyter-execute::

    exon.get_slice(allow_gaps=True)


``as_one_span`` unifies features with discontinuous alignment coordinates
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To get positions spanned by a feature, including gaps, use ``as_one_span``.

.. jupyter-execute::

    unified = aln_exons.as_one_span()
    aln3[unified]

Behaviour of annotations on nucleic acid sequences
""""""""""""""""""""""""""""""""""""""""""""""""""

Reverse complementing a sequence **does not** reverse annotations, that is they retain the reference to the frame for which they were defined.

.. jupyter-execute::

    plus = make_seq("CCCCCAAAAAAAAAATTTTTTTTTTAAAGG", moltype="dna")
    plus_rpt = plus.add_feature(biotype="blah", name="a", spans=[(5, 15), (25, 28)])
    plus[plus_rpt]

.. jupyter-execute::

    minus = plus.rc()
    minus

.. jupyter-execute::

    minus_rpt = list(minus.get_features(biotype="blah"))[0]
    minus[minus_rpt]

Masking annotated regions
"""""""""""""""""""""""""

We mask the CDS regions.

.. jupyter-execute::

    from cogent3 import load_seq

    seq = load_seq("data/ST_genome_part.gb", moltype="dna")
    no_cds = seq.with_masked_annotations("CDS")
    no_cds[150:400]

The above sequence could then have positions filtered so no position with the ambiguous character '?' was present.

Masking annotated regions on alignments
"""""""""""""""""""""""""""""""""""""""

We mask exon's on an alignment.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[["x", "C-CCCAAAAAGGGAA"], ["y", "-T----TTTTG-GTT"]],
        moltype="dna",
        array_align=False,
    )
    exon = aln.add_feature(
        seqid="x", biotype="exon", name="norwegian", spans=[(0, 4)], on_alignment=False
    )
    aln.with_masked_annotations("exon", mask_char="?")

After a reverse complement operation

.. jupyter-execute::

    rc = aln.rc()
    rc

these persist.

.. jupyter-execute::

    rc.with_masked_annotations("exon", mask_char="?")

You can take mask of the shadow
"""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import make_seq

    s = make_seq("CCCCAAAAAGGGAA", name="x", moltype="dna")
    exon = s.add_feature(biotype="exon", name="norwegian", spans=[(0, 4)])
    rpt = s.add_feature(biotype="repeat", name="norwegian", spans=[(9, 12)])
    s.with_masked_annotations("exon", shadow=True)

.. jupyter-execute::

    rc = s.rc()
    rc.with_masked_annotations("exon", shadow=True)

.. jupyter-execute::

    s.with_masked_annotations(["exon", "repeat"], shadow=True)

.. jupyter-execute::

    rc.with_masked_annotations(["exon", "repeat"], shadow=True)

What features of a certain type are available?
""""""""""""""""""""""""""""""""""""""""""""""

Query the database directly

.. jupyter-execute::

    from cogent3 import load_seq

    s = load_seq("data/ST_genome_part.gb", moltype="dna")
    s.annotation_db.describe

There are no user created features, but 31 features from the genbank annotations.

Getting all features of a type
""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_seq

    seq = load_seq("data/ST_genome_part.gb", moltype="dna")
    all_cds = list(seq.get_features(biotype="CDS"))
    coding_seqs = all_cds[0].union(all_cds[1:])
    coding_seqs

.. jupyter-execute::

    coding_seqs.get_slice()

.. jupyter-execute::

    noncoding_seqs = coding_seqs.shadow()
    noncoding_seqs

.. jupyter-execute::

    noncoding_seqs.get_slice()

Annotation display on sequences
"""""""""""""""""""""""""""""""

We can display annotations on sequences, writing to file. Using ``seq`` defined above

.. jupyter-execute::

    from cogent3 import load_seq

    seq = load_seq("data/ST_genome_part.gb", moltype="dna")
    fig = seq.get_drawable()
    fig.show()

.. following cleans up files

.. jupyter-execute::
    :hide-code:

    from cogent3.util.io import remove_files

    remove_files(["annotated_%d.png" % i for i in range(1, 4)], error_on_missing=False)
