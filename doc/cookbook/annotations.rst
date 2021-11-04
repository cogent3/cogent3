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

We load a sample genbank file with plenty of features and grab the CDS features.

.. jupyter-execute::

    from cogent3.parse.genbank import RichGenbankParser

    parser = RichGenbankParser(open("data/ST_genome_part.gb"))
    for accession, seq in parser:
        print(accession)

.. jupyter-execute::

    cds = seq.get_annotations_matching("CDS")
    print(cds)

Customising annotation construction from reading a genbank file
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can write your own code to construct annotation objects. One reason you might do this is some genbank files do not have a ``/gene`` tag on gene related features, instead only possessing a ``/locus_tag``. For illustrating the approach we only create annotations for ``CDS`` features. We write a custom callback function that uses the ``locus_tag`` as the ``Feature`` name.

.. jupyter-execute::

    from cogent3.core.annotation import Feature

    def add_annotation(seq, feature, spans):
        type_ = feature["type"]
        if type_ != "CDS":
            return
        name = feature.get("locus_tag", None)
        if name and not isinstance(name, str):
            name = " ".join(name)
        seq.add_annotation(Feature, type_, name, spans)

    parser = RichGenbankParser(
        open("data/ST_genome_part.gb"), add_annotation=add_annotation
    )
    for accession, seq in parser:  # just reading one accession,sequence
        break
    genes = seq.get_annotations_matching("CDS")
    print(genes)

Creating directly on a sequence
"""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import DNA
    from cogent3.core.annotation import Feature

    s1 = DNA.make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAA" "TTTTTTTTTTAAAAAGGGAACCCT", name="seq1"
    )
    print(s1[10:15])  # this will be exon 1
    print(s1[30:40])  # this will be exon 2
    print(s1[45:48])  # this will be exon 3
    s2 = DNA.make_seq("CGAAACGTTT", name="seq2")
    s3 = DNA.make_seq("CGAAACGTTT", name="seq3")

Via
"""

``add_annotation``
++++++++++++++++++

.. jupyter-execute::

    from cogent3 import DNA
    from cogent3.core.annotation import Feature

    s1 = DNA.make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAA" "TTTTTTTTTTAAAAAGGGAACCCT", name="seq1"
    )
    exon1 = s1.add_annotation(Feature, "exon", "A", [(10, 15)])
    exon2 = s1.add_annotation(Feature, "exon", "B", [(30, 40)])

``add_feature``
+++++++++++++++

.. jupyter-execute::

    from cogent3 import DNA

    s1 = DNA.make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAA" "TTTTTTTTTTAAAAAGGGAACCCT", name="seq1"
    )
    exon3 = s1.add_feature("exon", "C", [(45, 48)])

*There are other annotation types.*

Adding as a series or item-wise
"""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import DNA

    s2 = DNA.make_seq("CGAAACGTTT", name="seq2")
    cpgs_series = s2.add_feature("cpgsite", "cpg", [(0, 2), (5, 7)])
    s3 = DNA.make_seq("CGAAACGTTT", name="seq3")
    cpg1 = s3.add_feature("cpgsite", "cpg", [(0, 2)])
    cpg2 = s3.add_feature("cpgsite", "cpg", [(5, 7)])

Taking the union of annotations
"""""""""""""""""""""""""""""""

Construct a pseudo-feature (``cds``) that's a union of other features (``exon1``, ``exon2``, ``exon3``).

.. jupyter-execute::

    from cogent3 import DNA

    s1 = DNA.make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAA" "TTTTTTTTTTAAAAAGGGAACCCT", name="seq1"
    )
    exon1 = s1.add_feature("exon", "A", [(10, 15)])
    exon2 = s1.add_feature("exon", "B", [(30, 40)])
    exon3 = s1.add_feature("exon", "C", [(45, 48)])
    cds = s1.get_region_covering_all([exon1, exon2, exon3])

Getting annotation coordinates
""""""""""""""""""""""""""""""

These are useful for doing custom things, e.g. you could construct intron features using the below.

.. jupyter-execute::

    cds.get_coordinates()

Annotations have shadows
""""""""""""""""""""""""

A shadow is a span representing everything but the annotation.

.. jupyter-execute::

    not_cds = cds.get_shadow()
    not_cds

Compare to the coordinates of the original.

.. jupyter-execute::

    cds

Adding to a sequence member of an alignment
"""""""""""""""""""""""""""""""""""""""""""

The following annotation is directly applied onto the sequence and so is in ungapped sequence coordinates.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln1 = make_aligned_seqs(
        data=[["x", "-AAACCCCCA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    seq_exon = aln1.get_seq("x").add_feature("exon", "A", [(3, 8)])

Adding to an alignment
""""""""""""""""""""""

We add an annotation directly onto an alignment. In this example we add a ``Variable`` that can be displayed as a red line on a drawing. The resulting annotation (``red_data`` here) is in **alignment coordinates**!

.. jupyter-execute::

    from cogent3.core.annotation import Variable

    red_data = aln1.add_annotation(
        Variable, "redline", "align", [((0, 15), 1), ((15, 30), 2), ((30, 45), 3)]
    )

Slicing sequences and alignments by annotations
"""""""""""""""""""""""""""""""""""""""""""""""

By a feature or coordinates returns same sequence span

.. jupyter-execute::

    from cogent3 import DNA

    s1 = DNA.make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAA" "TTTTTTTTTTAAAAAGGGAACCCT", name="seq1"
    )
    exon1 = s1.add_feature("exon", "A", [(10, 15)])
    exon2 = s1.add_feature("exon", "B", [(30, 40)])
    s1[exon1]
    s1[10:15]

Using the annotation object ``get_slice`` method returns the same thing.

.. jupyter-execute::

    s1[exon2]
    exon2.get_slice()

Slicing by pseudo-feature or feature series
"""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import DNA

    s1 = DNA.make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAA" "TTTTTTTTTTAAAAAGGGAACCCT", name="seq1"
    )
    exon1 = s1.add_feature("exon", "A", [(10, 15)])
    exon2 = s1.add_feature("exon", "B", [(30, 40)])
    exon3 = s1.add_feature("exon", "C", [(45, 48)])
    cds = s1.get_region_covering_all([exon1, exon2, exon3])
    print(s1[cds])
    print(s1[exon1, exon2, exon3])

.. warning:: Slices are applied in order!

.. jupyter-execute::

    print(s1)
    print(s1[exon1, exon2, exon3])
    print(s1[exon2])
    print(s1[exon3])
    print(s1[exon1, exon3, exon2])

Slice series must not be overlapping
""""""""""""""""""""""""""""""""""""

.. jupyter-execute::
    :raises: ValueError

    s1[1:10, 9:15]
    s1[exon1, exon1]

But ``get_region_covering_all`` resolves this, ensuring no overlaps.

.. jupyter-execute::

    print(s1.get_region_covering_all([exon3, exon3]).get_slice())

You can slice an annotation itself
""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    print(s1[exon2])
    ex2_start = exon2[0:3]
    print(s1[ex2_start])
    ex2_end = exon2[-3:]
    print(s1[ex2_end])

Sequence vs Alignment slicing
"""""""""""""""""""""""""""""

You can't slice an alignment using an annotation from a sequence.

.. jupyter-execute::
    :raises: ValueError

    aln1[seq_exon]

Copying annotations
"""""""""""""""""""

You can copy annotations onto sequences with the same name, even if the length differs

.. jupyter-execute::

    aln2 = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
    match_exon = seq.add_feature("exon", "A", [(3, 8)])
    aln2.get_seq("x").copy_annotations(seq)
    copied = list(aln2.get_annotations_from_seq("x", "exon"))
    copied

but if the feature lies outside the sequence being copied to, you get a lost span

.. jupyter-execute::

    aln2 = make_aligned_seqs(data=[["x", "-AAAA"], ["y", "TTTTT"]], array_align=False)
    seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
    match_exon = seq.add_feature("exon", "A", [(5, 8)])
    aln2.get_seq("x").copy_annotations(seq)
    copied = list(aln2.get_annotations_from_seq("x", "exon"))
    copied
    copied[0].get_slice()

You can copy to a sequence with a different name, in a different alignment if the feature lies within the length

.. jupyter-execute::

    # new test
    aln2 = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCC", "x")
    match_exon = seq.add_feature("exon", "A", [(5, 8)])
    aln2.get_seq("y").copy_annotations(seq)
    copied = list(aln2.get_annotations_from_seq("y", "exon"))
    copied

If the sequence is shorter, again you get a lost span.

.. jupyter-execute::

    aln2 = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    diff_len_seq = DNA.make_seq("CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "x")
    nonmatch = diff_len_seq.add_feature("repeat", "A", [(12, 14)])
    aln2.get_seq("y").copy_annotations(diff_len_seq)
    copied = list(aln2.get_annotations_from_seq("y", "repeat"))
    copied

Querying
""""""""

You need to get a corresponding annotation projected into alignment coordinates via a query.

.. jupyter-execute::

    aln_exon = aln1.get_annotations_from_any_seq("exon")
    print(aln1[aln_exon])

Querying produces objects only valid for their source
"""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::
    :raises: ValueError

    cpgsite2 = s2.get_annotations_matching("cpgsite")
    print(s2[cpgsite2])
    cpgsite3 = s3.get_annotations_matching("cpgsite")
    s2[cpgsite3]

Querying for absent annotation
""""""""""""""""""""""""""""""

You get back an empty list, and slicing with this returns an empty sequence.

.. jupyter-execute::

    # this test is new
    dont_exist = s2.get_annotations_matching("dont_exist")
    dont_exist
    s2[dont_exist]

Querying features that span gaps in alignments
""""""""""""""""""""""""""""""""""""""""""""""

If you query for a feature from a sequence, it's alignment coordinates may be discontinuous.

.. jupyter-execute::

    aln3 = make_aligned_seqs(
        data=[["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]], array_align=False
    )
    exon = aln3.get_seq("x").add_feature("exon", "ex1", [(0, 4)])
    print(exon.get_slice())
    aln_exons = list(aln3.get_annotations_from_seq("x", "exon"))
    print(aln_exons)
    print(aln3[aln_exons])

.. note:: The ``T`` opposite the gap is missing since this approach only returns positions directly corresponding to the feature.

``as_one_span`` unifies features with discontinuous alignment coordinates
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To get positions spanned by a feature, including gaps, use ``as_one_span``.

.. jupyter-execute::

    unified = aln_exons[0].as_one_span()
    print(aln3[unified])

Behaviour of annotations on nucleic acid sequences
""""""""""""""""""""""""""""""""""""""""""""""""""

Reverse complementing a sequence **does not** reverse annotations, that is they retain the reference to the frame for which they were defined.

.. jupyter-execute::

    plus = DNA.make_seq("CCCCCAAAAAAAAAATTTTTTTTTTAAAGG")
    plus_rpt = plus.add_feature("blah", "a", [(5, 15), (25, 28)])
    print(plus[plus_rpt])
    minus = plus.rc()
    print(minus)
    minus_rpt = minus.get_annotations_matching("blah")
    print(minus[minus_rpt])

Masking annotated regions
"""""""""""""""""""""""""

We mask the CDS regions.

.. jupyter-execute::

    from cogent3.parse.genbank import RichGenbankParser

    parser = RichGenbankParser(open("data/ST_genome_part.gb"))
    seq = [seq for accession, seq in parser][0]
    no_cds = seq.with_masked_annotations("CDS")
    print(no_cds[150:400])

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
    exon = aln.get_seq("x").add_feature("exon", "norwegian", [(0, 4)])
    print(aln.with_masked_annotations("exon", mask_char="?"))

These also persist through reverse complement operations.

.. jupyter-execute::

    rc = aln.rc()
    print(rc)
    print(rc.with_masked_annotations("exon", mask_char="?"))

You can take mask of the shadow
"""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import DNA

    s = DNA.make_seq("CCCCAAAAAGGGAA", "x")
    exon = s.add_feature("exon", "norwegian", [(0, 4)])
    rpt = s.add_feature("repeat", "norwegian", [(9, 12)])
    rc = s.rc()
    print(s.with_masked_annotations("exon", shadow=True))
    print(rc.with_masked_annotations("exon", shadow=True))
    print(s.with_masked_annotations(["exon", "repeat"], shadow=True))
    print(rc.with_masked_annotations(["exon", "repeat"], shadow=True))

What features of a certain type are available?
""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import DNA

    s = DNA.make_seq("ATGACCCTGTAAAAAATGTGTTAACCC", name="a")
    cds1 = s.add_feature("cds", "cds1", [(0, 12)])
    cds2 = s.add_feature("cds", "cds2", [(15, 24)])
    all_cds = s.get_annotations_matching("cds")
    all_cds

Getting all features of a type, or everything but that type
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The annotation methods ``get_region_covering_all`` and ``get_shadow`` can be used to grab all the coding sequences or non-coding sequences in a ``DnaSequence`` object.

.. jupyter-execute::

    from cogent3.parse.genbank import RichGenbankParser

    parser = RichGenbankParser(open("data/ST_genome_part.gb"))
    seq = [seq for accession, seq in parser][0]
    all_cds = seq.get_annotations_matching("CDS")
    coding_seqs = seq.get_region_covering_all(all_cds)
    coding_seqs
    coding_seqs.get_slice()
    noncoding_seqs = coding_seqs.get_shadow()
    noncoding_seqs
    noncoding_seqs.get_slice()

Getting sequence features when you have an alignment object
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sequence features can be accessed via a containing ``Alignment``.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    print(aln)
    exon = aln.get_seq("x").add_feature("exon", "1", [(3, 8)])
    aln_exons = aln.get_annotations_from_seq("x", "exon")
    aln_exons = aln.get_annotations_from_any_seq("exon")
    aln_exons

Annotation display on sequences
"""""""""""""""""""""""""""""""

We can display annotations on sequences, writing to file.

We first make a sequence and add some annotations.

.. jupyter-execute::

    from cogent3 import DNA

    seq = DNA.make_seq("aaaccggttt" * 10)
    v = seq.add_feature("exon", "exon", [(20, 35)])
    v = seq.add_feature("repeat_unit", "repeat_unit", [(39, 49)])
    v = seq.add_feature("repeat_unit", "rep2", [(49, 60)])

.. todo:: document info attribute
.. todo:: document how to visualise annotations

.. following cleans up files

.. jupyter-execute::
    :hide-code:

    from cogent3.util.io import remove_files

    remove_files(["annotated_%d.png" % i for i in range(1, 4)], error_on_missing=False)