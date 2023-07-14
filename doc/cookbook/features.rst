.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _howto-features:

Features
--------

This guide provides instructions on creating, querying, and utilising features to manipulate biological sequence data.

How to create a custom ``Feature``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Via ``add_feature``
""""""""""""""""""""

We can dynamically create a ``Feature`` via the ``add_feature()`` method, providing the key information about the feature.

The new feature will be added to the ``annotation_db`` attribute of the ``Sequence`` and/or ``Alignment``. A ``Feature`` instance will be returned.

On a ``Sequence``
+++++++++++++++++

We use ``add_feature`` to add a ``Feature`` to a ``Sequence``, providing the biotype, name/id, and a list of start and stop indices.

.. jupyter-execute::
    :raises:

    from cogent3 import make_seq

    # create an example sequence
    seq = make_seq(
        "GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA",
        name="seq1",
        moltype="dna",
    )

    # Add a feature
    seq.add_feature(biotype="gene", name="a_gene", spans=[(12, 48)], seqid="seq1")

A feature can also be added as a series of non-overlapping genomic segments:

.. jupyter-execute::
    :raises:

    # Add a feature as a series
    seq.add_feature(
        biotype="cpgsite",
        name="a_cpgsite",
        spans=[(12, 18), (21, 29), (35, 48)],
        seqid="seq1",
    )

On a ``Sequence`` within an ``Alignment``
+++++++++++++++++++++++++++++++++++++++++

We use ``add_feature`` and provide a value for ``seqid`` to add a feature to a sequence within an ``Alignment`` or ``SequenceCollection``. 

Note the difference between the value provided to ``spans``, and the value of ``map`` in the returned ``Feature``... the resulting feature is in **alignment coordinates**!

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    # make demo alignment
    aln1 = make_aligned_seqs(
        data=[["seq1", "-AAACCCCCA"], ["seq2", "TTTT--TTTT"]], array_align=False
    )
    # add feature to seq1
    aln1.add_feature(
        seqid="seq1", biotype="exon", name="A", spans=[(3, 8)], on_alignment=False
    )

The ``Feature`` specifies that ``seqid="from 'seq1'"``, indicating that it "belongs" to seq1.

On an ``Alignment``
+++++++++++++++++++

We use ``add_feature`` to add a ``Feature`` to an ``Alignment``.

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    # make demo alignment
    aln1 = make_aligned_seqs(
        data=[["seq1", "-AAACCCCCA"], ["seq2", "TTTT--TTTT"]], array_align=False
    )

    aln1.add_feature(
        biotype="exon",
        name="aligned_exon",
        spans=[(0, 8)],
        on_alignment=True,
    )

The ``Feature`` specifies that ``seqid=None``, indicating that it belongs to the alignment

Via an ``AnnotationDb``
+++++++++++++++++++++++

We can use ``add_feature`` to add a feature directly into an ``AnnotationDb``, and assign it to the ``annotation_db`` attribute of a ``Sequence`` or ``Alignment``. For extensive documentation on handling features directly via an ``AnnotationDb`` see :ref:`anno_db`.

.. jupyter-execute::

    from cogent3 import make_seq
    from cogent3.core.annotation_db import BasicAnnotationDb

    # init empty db and add feature
    db = BasicAnnotationDb()
    db.add_feature(seqid="seq1", biotype="exon", name="C", spans=[(45, 48)])

    # make demo seq
    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT", name="seq1", moltype="dna"
    )

    # assign db to sequence
    s1.annotation_db = db
    s1.annotation_db

How to load bulk Features from a File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically, we want to load bulk features from a genomic annotation file, such as a GFF or Genbank file. For the following examples, we will use *Caenorhabditis elegans* chromosome I.

.. note:: See the list of :ref:`data_links` to download the data used in the following examples.

To load features from a genomic annotation file along with the corresponding sequence, we can use the ``load_seq`` function. The features are stored in a ``AnnotationDb`` and assigned to the ``annotation_db`` attribute of the sequence.

From a Genbank file
"""""""""""""""""""

How to load features and sequence data
++++++++++++++++++++++++++++++++++++++

To load the sequence and all 40,578 features from *C. elegans* Chromosome 1, we use the ``load_seq`` function. The loading process takes approximately 1.4 seconds ⚡️

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq
    
    %timeit load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")

.. jupyter-execute::
    :raises:
    :hide-code:

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")

The features are stored in the ``annotation_db`` attribute.

.. jupyter-execute::
    :raises:

    seq.annotation_db

Now that the ``Sequence`` is annotated, we can query it for specific features. For more details on querying, skip to :ref:`Querying for Features <query_for_features>`.

From a GFF file
"""""""""""""""

How to load features and sequence data
++++++++++++++++++++++++++++++++++++++

Given a FASTA file with sequence data and a GFF file with annotations, we can use ``load_seq`` to load both the sequence and its corresponding features.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq(
        filename="data/C-elegans-chromosome-I.fa",
        annotation_path="data/C-elegans-chromosome-I.gff",
        moltype="dna",
    )
    seq.annotation_db

.. warning:: ``total_records=0``? 🤔 This is because ``load_seq`` assumes the sequence names match exactly between files! If the names are different, you need to provide function to the ``label_to_name`` argument.

Because the names above are different, for FASTA its ``"I dna:chromosome chromosome:WBcel235:I:1:15072434:1 REF"`` and for GFF its ``"I"``, we need a ``label_to_name`` argument. We provide a lambda function.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq(
        "data/C-elegans-chromosome-I.fa",
        annotation_path="data/C-elegans-chromosome-I.gff",
        label_to_name=lambda x: x.split()[0],
        moltype="dna",
    )
    seq.annotation_db

How to load features and associate them with an existing sequence
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We can use the ``annotate_from_gff()`` method to associate the features from a GFF file with the existing ``Sequence``.

If we know that the features and the sequence share the same coordinate space, then we only need to provide the path to the annotation file.

.. jupyter-execute::
    :hide-code:

    from cogent3 import load_seq

    loaded_seq = load_seq(
        "data/C-elegans-chromosome-I.fa",
        label_to_name=lambda x: x.split()[0],
        moltype="dna",
    )

.. jupyter-execute::
    :raises:

    # loaded_seq = < loaded / created the seq>
    loaded_seq.annotate_from_gff("data/C-elegans-chromosome-I.gff")
    loaded_seq.annotation_db

If the feature coordinates precede the sequence, for example, a sequence corresponds to a gene that starts 600 base pairs from the beginning of chromosome, but the annotation file is for the entire chromosome, we need to provide an offset to the ``annotate_from_gff()`` method.

.. jupyter-execute::
    :hide-code:

    from cogent3 import make_seq

    sub_seq = make_seq(
        "GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAAGACTAAGCCTAAGA",
        name="I",
        moltype="dna",
    )

.. jupyter-execute::
    :raises:

    # sub_seq = <genomic region starting at the 600th nt>
    sub_seq.annotate_from_gff("data/C-elegans-chromosome-I.gff", offset=600)
    sub_seq.annotation_db

How to load features and associate them with sequences in an existing alignment
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To annotate one or more ``Sequence`` in an ``Alignment``, call ``annotate_from_gff()`` on the ``Alignment`` instance, passing in the path to the GFF annotation file and a list of sequence names to annotate to the ``seq_ids`` argument.

For example, first we load an alignment of the brca1 gene in primates.

.. jupyter-execute::
    :raises:

    from cogent3 import load_aligned_seqs

    brca1_aln = load_aligned_seqs(
        "data/primate_brca1.fasta", array_align=False, moltype="dna"
    )
    brca1_aln

Next, we annotate with a GFF file that contains features specific to the human gene.
    
.. jupyter-execute::
    :raises:

    brca1_aln.annotate_from_gff("data/brca1_hsa_shortened.gff", seq_ids=["Human"])
    brca1_aln.annotation_db

Note that the ``AnnotationDb`` is accessible via the ``Alignment`` (above) and ``Sequence`` (below) attribute.

.. jupyter-execute::
    :raises:

    brca1_aln.get_seq("Human").annotation_db

.. note:: ``Alignment.annotate_from_gff()`` does not support setting an offset. If you need to set the offset for a sequence within an alignment, you can do so directly using the ``Sequence.annotation_offset`` attribute.

.. _query_for_features:

How to query a Sequence or Alignment for Features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The method ``get_features`` yields all features that match the given arguments. You can provide conditions for the name, biotype, and start/stop location of a feature.

Querying a ``Sequence`` for Features
""""""""""""""""""""""""""""""""""""

Querying via Feature Name
+++++++++++++++++++++++++

We can search for a gene given its name (AKA its unique ID). For example we can search for a gene with ``name="WBGene00021661"``.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")

    # note we wrap `get_features` in `list` as generator is returned
    gene = list(seq.get_features(name="WBGene00021661", biotype="gene"))
    gene

Querying via Feature Biotype
++++++++++++++++++++++++++++

We can search for features with a certain biotype, for example, all coding sequences (CDS):

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    cds = list(seq.get_features(biotype="CDS"))
    cds[:3]

We can also provide combinations of argument to search, for example, all CDS with a given name:

.. jupyter-execute::
    :raises:

    cds = list(seq.get_features(biotype="CDS", name="WBGene00021661"))
    cds

Querying via region of interest
+++++++++++++++++++++++++++++++

We can provide ``start`` and ``end`` arguments to ``get_features()`` and all features within the coordinates will be returned.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    region_features = list(seq.get_features(start=10148, stop=26732))
    region_features[:3]

We can again provide a combination of conditions, for example, querying for all features with ``biotype="mRNA"`` within a certain range, and returning the first match.

.. jupyter-execute::
    :raises:

    mRNA = list(seq.get_features(start=10148, stop=29322, biotype="mRNA"))[0]
    mRNA

Querying a Sequence (via an Alignment) for Features
"""""""""""""""""""""""""""""""""""""""""""""""""""

To query for a particular ``Sequence`` within an ``Alignment`` or ``SequenceCollection``, we can use ``get_features`` as shown above for a ``Sequence``, but providing the seqid for the sequence of interest.

For example, given an alignment of primates, we can search for features that are just on the human sequence as follows:

.. jupyter-execute::
    :raises:

    from cogent3 import load_aligned_seqs

    # first load alignment and annotate the human seq
    aln = load_aligned_seqs(
        "data/primate_brca1.fasta", array_align=False, moltype="dna"
    )
    aln.annotate_from_gff("data/brca1_hsa_shortened.gff", seq_ids=["Human"])

    # query alignment providing seqid of interest
    human_exons = list(aln.get_features(biotype="exon", seqid="Human"))
    human_exons

Note that ``seqid="from'Human'"`` indicated this feature belongs to this particular sequence.

Querying an Alignment for Features
""""""""""""""""""""""""""""""""""

Querying for features on any ``Sequence`` in an ``Alignment``
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

todo: ``on_alignment=False`` and dont provide seqid

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    # add a feature to the alignment we created above on difference sequence
    aln.add_feature(biotype="gene", name="gene:101", spans=[(40, 387)], seqid="Rhesus")

    any_feature = list(aln.get_features(on_alignment=False))
    any_feature

Note there are features from both Rhesus, which we just added, and Human, which we annotated above

Querying for features on an ``Alignment``
+++++++++++++++++++++++++++++++++++++++++

todo: ``on_alignment=True`` and dont provide seqid

Using ``add_feature`` we add a feature to the brca1 alignment we have been using above, by specifying ``on_alignment=True`` this feature will be on the ``Alignment``. 

To query for features on the alignment, we use ``get_features``, again specifying ``on_alignment=True``. 

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    # first we add the feature to the alignment
    aln.add_feature(
        biotype="pseudogene", name="pseudogene1", spans=[(420, 666)], on_alignment=True
    )

    # query for features on the alignment
    aln_features = list(aln.get_features(on_alignment=True))
    aln_features

Note how even though we annotated the Human and Rhesus sequences in the above examples, only the pseudogene we added to ``Alignment`` is returned by this query.  

Querying features that span gaps in alignments
++++++++++++++++++++++++++++++++++++++++++++++

If you query for a ``Feature`` from a ``Sequence`` (i.e. the feature is in sequence coordinates), its alignment coordinates may be discontinuous. This will lead to an omission of data from other sequences!

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln3 = make_aligned_seqs(
        data=[["x", "C-CCCAAAAA"], ["y", "-T----TTTT"]],
        array_align=False,
        moltype="dna",
    )
    exon = aln3.add_feature(
        seqid="x", biotype="exon", name="ex1", spans=[(0, 4)], on_alignment=False
    )
    exon.get_slice()

.. jupyter-execute::
    :raises:

    aln_exons = list(aln3.get_features(seqid="x", biotype="exon"))[0]
    aln_exons

.. note:: In the above, the ``T`` in sequence Y opposite the gap is missing since this approach only returns positions directly corresponding to the feature.

To include the gaps, use the ``allow_gaps`` argument

.. jupyter-execute::
    :raises:

    exon.get_slice(allow_gaps=True)

Examples using the methods available on Features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``Feature`` has many methods to manipulate the sequence or alignment that they are bound to.

How to slice a ``Sequence`` or ``Alignment`` by its features
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Given a ``Feature``, we can directly slice its parent sequence to return its sequence information

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq(
        "data/C-elegans-chromosome-I.fa",
        annotation_path="data/C-elegans-chromosome-I.gff",
        label_to_name=lambda x: x.split()[0],
        moltype="dna",
    )
    pseudogene = list(seq.get_features(start=10148, stop=26732, biotype="pseudogene"))[0]
    seq[pseudogene]

.. note:: This only works for the ``Sequence`` that the ``Feature`` "belongs" to.

We can also achieve this via ``get_slice()``

.. jupyter-execute::
    :raises:

    pseudogene.get_slice()

How to display the features of a Sequence
"""""""""""""""""""""""""""""""""""""""""

We can display all the features on a sequence using ``.get_drawable()``. We show it for only the first 50,000 base pairs. The plotly figure returned, as displayed below, is interactive! 🤩

.. jupyter-execute::

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    subseq = seq[25000:35000]
    fig = subseq.get_drawable()
    fig.show()

How to find the coordinates of a feature
""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::
    :raises:

    pseudogene.get_coordinates()

These are useful for doing custom things, e.g. if the introns are not annotated for a gene, we can generate the introns from the coordinates of the exons as follows:

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    cds = list(seq.get_features(biotype="CDS"))[0]
    exon_coords = cds.get_coordinates()

    exon_coords

We generate the intron coordinates from the second element of the first tuple, and the first element of the second tuple and so on:

.. jupyter-execute::
    :raises:

    intron_coords = []

    for i in range(len(exon_coords) - 1):
        intron_coords.append((exon_coords[i][1], exon_coords[i + 1][0]))

    intron_coords

We can then add the introns as a ``Feature`` to the sequence!

.. jupyter-execute::
    :raises:

    intron = seq.add_feature(
        biotype="intron", name="intron:Y74C9A.3.1", seqid="I", spans=intron_coords
    )
    intron

How to take the union of features
"""""""""""""""""""""""""""""""""

We can create a feature that is the union of all coding sequence.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    cds = list(seq.get_features(biotype="CDS"))
    union_cds = cds[0].union(cds[1:])

How to get the shadow of a Feature
""""""""""""""""""""""""""""""""""

The "shadow" of a feature is a new feature containing all of the sequence **except the feature**!

How to use the shadow of a Feature to return the intergenic sequence
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We first need to query our sequence for all genes. Using the ``union()`` method we combine all genes into a single feature.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    genes = list(seq.get_features(biotype="gene"))
    genes = genes[0].union(genes[1:])
    genes

Taking the "shadow" of all genes will return the intergenic region as a valid ``Feature``

.. jupyter-execute::
    :raises:

    intergenic = genes.shadow()
    intergenic

We can slice the sequence by this new Feature to return the complete intergenic sequence!

.. jupyter-execute::
    :raises:

    intergenic.get_slice()

How to mask annotated regions
"""""""""""""""""""""""""""""

Masking annotated regions on a sequence
+++++++++++++++++++++++++++++++++++++++

We can mask a certain annotation using ``with_masked_annotations()``

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    no_cds = seq.with_masked_annotations("CDS")
    no_cds[2575800:2575900]

The above sequence could then have positions filtered so no position with the ambiguous character '?' was present.

Masking annotated regions on an Alignment
+++++++++++++++++++++++++++++++++++++++++

We can mask exons on an alignment.

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=[["x", "C-CCCAAAAAGGGAA"], ["y", "-T----TTTTG-GTT"]],
        moltype="dna",
        array_align=False,
    )
    exon = aln.add_feature(
        seqid="x",
        biotype="exon",
        name="exon-be-gone",
        spans=[(0, 4)],
        on_alignment=False,
    )
    aln.with_masked_annotations("exon", mask_char="?")

After a reverse complement operation

.. jupyter-execute::
    :raises:

    rc = aln.rc()
    rc

these persist.

.. jupyter-execute::
    :raises:

    rc.with_masked_annotations("exon", mask_char="?")

How to find the "children" of a Feature
"""""""""""""""""""""""""""""""""""""""

To find the "children" of a feature, we can use the ``get_children()`` method. A "child" refers to a feature that is nested within or contained by another "parent" feature. For example, a child feature could be an exon contained within a gene or a CDS contained within a transcript.

This method returns a generator that yields all the child features of the specified feature.

For example, let's find the children of the gene "WBGene00021661":

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    seq = load_seq(
        "data/C-elegans-chromosome-I.fa",
        annotation_path="data/C-elegans-chromosome-I.gff",
        label_to_name=lambda x: x.split()[0],
        moltype="dna",
    )
    gene = list(seq.get_features(name="gene:WBGene00022276", biotype="gene"))[0]
    children = list(gene.get_children())
    children

How to find the "parent" of a Feature
"""""""""""""""""""""""""""""""""""""

To find the "parent" of a feature, we can use the ``get_parent()`` method, which achieves the inverse of the above method.

For example, we can use the first "child" we returned above, ``"transcript:Y74C9A.2a.3"``, to find the original parent gene!

.. jupyter-execute::
    :raises:

    child = list(seq.get_features(name="transcript:Y74C9A.2a.3", biotype="mRNA"))[0]
    parent = list(child.get_parent())
    parent

How to copy features
""""""""""""""""""""

We can copy features onto sequences with the same name. Note that the ``AnnotationDb`` instance bound to the alignment and its member sequences is the **same**.

.. jupyter-execute::
    :raises:

    aln2 = make_aligned_seqs(
        data=[["x", "-AAAAAAAAA"], ["y", "TTTT--TTTT"]],
        array_align=False,
        moltype="dna",
    )
    x, y = aln2.get_seq("x"), aln2.get_seq("y")
    x.annotation_db is y.annotation_db is aln2.annotation_db

.. warning:: Despite this, it is possible for the attributes to get out-of-sync. So, any copy annotations should be done using ``alignment.copy_annotations()``, **not** ``alignment.get_seq("x").copy_annotations()``.

.. jupyter-execute::
    :raises:

    seq = make_seq("CCCCCCCCCCCCCCCCCCCC", name="x", moltype="dna")
    match_exon = seq.add_feature(biotype="exon", name="A", spans=[(3, 8)])
    aln2.copy_annotations(seq.annotation_db)
    aln2.annotation_db

However, if the feature lies outside the sequence being copied to, you get a lost span

.. jupyter-execute::
    :raises:

    copied = list(aln2.get_features(seqid="x", biotype="exon"))
    copied

How to get the positions of a feature as one span
"""""""""""""""""""""""""""""""""""""""""""""""""

``as_one_span`` unifies features with discontinuous alignment coordinates and returns positions spanned by a feature, including gaps.

.. jupyter-execute::
    :raises:

    unified = aln_exons.as_one_span()
    aln3[unified]

Behaviour of annotations on nucleic acid sequences
""""""""""""""""""""""""""""""""""""""""""""""""""

Reverse complementing a sequence **does not** reverse features. Features are considered to have strand specific meaning (.e.g CDS, exons) and so they retain the reference to the frame for which they were defined.

.. jupyter-execute::
    :raises:

    plus = make_seq("CCCCCAAAAAAAAAATTTTTTTTTTAAAGG", moltype="dna")
    plus_rpt = plus.add_feature(biotype="blah", name="a", spans=[(5, 15), (25, 28)])
    plus[plus_rpt]

.. jupyter-execute::
    :raises:

    minus = plus.rc()
    minus

.. jupyter-execute::
    :raises:

    minus_rpt = list(minus.get_features(biotype="blah"))[0]
    minus[minus_rpt]
