.. jupyter-execute::
    :hide-code:

    import set_working_directory

.. _intro_annotations:

Annotations
-----------

This guide provides instructions on creating, querying, and utilising Features to manipulate biological sequence data. For more extensive documentation about annotations, see :ref:`seq-annotations`.

Creating custom Features
^^^^^^^^^^^^^^^^^^^^^^^^

Via ``add_features``
""""""""""""""""""""

One way to dynamically create a ``Feature`` is via the ``add_feature()`` method, providing the biotype, name/id, and a list of start and stop indices. 

The new feature will be added to the ``annotation_db`` attribute of the ``Sequence`` and or ``Alignment``. A ``Feature`` instance will be returned.

On a ``Sequence``
+++++++++++++++++

We use ``add_feature`` to add a feature to a sequence.

.. jupyter-execute::
    :raises:

    from cogent3 import make_seq

    seq = make_seq(
        "GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA",
        name="seq1",
        moltype="dna",
    )

    # Add a feature item-wise
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

We use ``add_feature`` to add a feature to a sequence within an alignment. The resulting feature is in **alignment coordinates**!

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln1 = make_aligned_seqs(
        data=[["x", "-AAACCCCCA"], ["y", "TTTT--TTTT"]], array_align=False
    )
    aln1.add_feature(
        seqid="x", biotype="exon", name="A", spans=[(3, 8)], on_alignment=False
    )

On an ``Alignment``
+++++++++++++++++++

We use ``add_feature`` to add a feature to an alignment. The resulting feature is in **alignment coordinates**!

.. jupyter-execute::
    :raises:

    aln1.add_feature(
        seqid=None,
        biotype="shared",
        name="demo",
        spans=[(0, 8)],
        on_alignment=True,
    )

Via an ``AnnotationDb``
+++++++++++++++++++++++

We use ``add_feature`` to add a feature directly into into an annotation database, and assign it to the ``annotation_db`` attribute of the sequence. We could also assign it to the ``annotation_db`` attribute of an alignment.

.. jupyter-execute::

    from cogent3 import make_seq
    from cogent3.core.annotation_db import BasicAnnotationDb

    db = BasicAnnotationDb()

    db.add_feature(seqid="seq1", biotype="exon", name="C", spans=[(45, 48)])
    s1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT", name="seq1", moltype="dna"
    )
    s1.annotation_db = db

Loading Features from a File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically, we want to load features from a genomic annotation file, such as a GFF or Genbank file. For the following examples, we will use data from *Caenorhabditis elegans* chromosome I.

.. note:: See the list of :ref:`data_links` to download the data used in the following examples.

To load features from a genomic annotation file along with the corresponding sequence, we can use the ``load_seq`` function. The features are stored in a ``AnnotationDb`` and assigned to the ``annotation_db`` attribute of the sequence. 

From a Genbank file
"""""""""""""""""""

To load the sequence and all 40,578 features from *C. elegans* Chromosome 1, use the ``load_seq`` function. The loading process takes approximately 1.5 seconds.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq
    
    %timeit load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")

.. jupyter-execute::
    :hide-code:

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")

The features are stored in the ``annotation_db`` attribute.

.. jupyter-execute::
    :raises:

    seq.annotation_db

We can query the sequence for specific features (more details to follow).

From a GFF file
"""""""""""""""

How to load features and sequence data
++++++++++++++++++++++++++++++++++++++

If you have the FASTA file for the sequence, you can use ``load_seq`` and provide the GFF file to the ``annotation_path`` argument.

.. jupyter-execute::
    :raises:

    seq = load_seq(
        "data/C-elegans-chromosome-I.fa",
        annotation_path="data/C-elegans-chromosome-I.gff",
        moltype="dna",
    )
    seq.annotation_db

.. warning:: ``total_records=0``? 🧐 Note that this method assumes an exact match between the names in the files! If the names are different, you need to provide a ``label_to_name`` argument. 

As the names are different in our example (``"I dna:chromosome chromosome:WBcel235:I:1:15072434:1 REF"`` in the FASTA file and ``"I"`` in the gff file) we need to provide a ``label_to_name`` argument as follows:

.. jupyter-execute::
    :raises:

    seq = load_seq(
        "data/C-elegans-chromosome-I.fa",
        annotation_path="data/C-elegans-chromosome-I.gff",
        label_to_name=lambda x: x.split()[0],
        moltype="dna",
    )
    seq.annotation_db

How to load features and associate them with an existing sequence
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If we know that the features lie within the sequence coordinates, we can use the ``annotate_from_gff()`` method to associate the features with the existing sequence.

.. jupyter-execute::
    :hide-code:

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

If the features precede the sequence, we can still use the ``annotate_from_gff()`` method, but we need to provide the offset value. For example, given a sequence that starts 600 base pairs from the beginning of chromosome 1, we can adjust the features as follows:

.. jupyter-execute::
    :hide-code:

    from cogent3 import make_seq

    # sub_seq starts at position 600, so we need to provide an offset.
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

How to load features and associate them with an existing alignment
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To annotate one or more sequences in an alignment, call ``annotate_from_gff()`` on the ``Alignment`` instance, passing in the path to the GFF annotation file and a list of sequence names to annotate:

This will create an annotation database that is accessible via the ``annotation_db`` attribute on both the ``Alignment`` and named ``Sequence`` instances. 

For example, we can grab the sequence alignment of the brca1 gene in primates, and an annotation file for human chromosome 17. 

todo: check for bug in the below.

.. jupyter-execute::
    :raises:

    from cogent3 import load_aligned_seqs
    
    brca1_aln = load_aligned_seqs("data/primate_brca1.fasta", array_align=False, moltype="dna")
    brca1_aln.annotate_from_gff("data/human_annotations.gff", seq_ids=["Human"])
    brca1_aln.annotation_db

.. note:: ``Alignment.annotate_from_gff()`` does not support setting an offset. If you need to set the offset for a sequence within an alignment, you can do so directly using the ``Sequence.annotation_offset`` attribute.


Our annotation file includes features for the whole of chromosome 17, so we need to set the offset. The start location of the brca1 gene in humans is position 43,044,295 so we set this as the annotation offset.

.. jupyter-execute::
    :raises:
    
    brca1_aln.get_seq("Human").annotation_offset = 43044295
    brca1_aln
    
Querying for Features
^^^^^^^^^^^^^^^^^^^^^

The method ``get_features`` yields all features that match the given arguments. You can provide conditions for the name, biotype, and start/stop location of a feature. 

Querying a Sequence for Features
""""""""""""""""""""""""""""""""

Querying via Feature Name
+++++++++++++++++++++++++

We can search for the a gene given its unique ID ``"WBGene00021661"``. We wrap the method call in a ``list()`` as it is a generator object.

.. jupyter-execute::
    :raises:

    mbtr_1 = list(seq.get_features(name="WBGene00021661", biotype="gene"))
    mbtr_1
    
Querying via Feature Biotype
++++++++++++++++++++++++++++

We can search for all CDS 

.. jupyter-execute::
    :raises:

    cds = list(seq.get_features(biotype="CDS"))
    cds[:3]

We can search for all CDS for a given gene

.. jupyter-execute::
    :raises:

    cds = list(seq.get_features(biotype="CDS", name="CDS:Y74C9A.3.1"))
    cds


Querying via region of interest
+++++++++++++++++++++++++++++++

We can provide ``start`` and ``end`` arguments to ``get_features()`` and all features within the coordinates will be returned.

.. jupyter-execute::
    :raises:

    region_features = list(seq.get_features(start=10148, stop=26732))
    region_features[:3]


We can also provide a combination of conditions, for example, querying for all pseudogenes within a certain range, and returning the first match.

.. jupyter-execute::
    :raises:

    pseudogene = list(seq.get_features(start=10148, stop=26732, biotype="pseudogene"))[0]
    pseudogene
    
Querying a Sequence (via an Alignment) for Features
"""""""""""""""""""""""""""""""""""""""""""""""""""

#need data C elegans aln of gene

To query for a particular sequence within and Alignment, we can use the ``get_features`` as above, but providing the seqid for the sequence of interest.

We can query for CDS features on sequence "x" in our example alignment:

.. jupyter-execute::
    :raises:

    seq_cds = list(aln1.get_features(seqid="x", biotype="CDS"))
    seq_cds


Querying an Alignment for Features
""""""""""""""""""""""""""""""""""

todo

.. jupyter-execute::
    :raises:

    from cogent3 import load_aligned_seqs
    
    brca1 = load_aligned_seqs("data/primate_brca1.fasta", array_align=False, moltype="dna")
    brca1.get_seq("Human").annotate_from_gff("/Users/u6096186/repos/cogent3/doc/data/human_annotations.gff")
    
    # Recall, human_annotations.gff includes annotation for chromosome 17, so we need to set the offset. 
    # the start location of brca1 on the chromosome is 43,044,295 so we set this as the annotation offset.
    brca1.get_seq("Human").annotation_offset = 43044295
    
    brca_exons = list(brca1.get_seq("Human").get_features(biotype="exon"))
    brca_exons


Example functionalities of Features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``Feature`` has many methods that can manipulate the sequence or alignment that they are bound to. 

How to slice a ``Sequence`` or ``Alignment`` by its features
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Given the pseudogene feature that we queried for above, we can directly slice its parent sequence to return its sequence information 

.. jupyter-execute::
    :raises:
    
    seq[pseudogene]

.. note:: This only works for the ``Sequence`` that the ``Feature`` "belongs" to.

We can also achieve this via ``get_slice()``

.. jupyter-execute::
    :raises:

    pseudogene.get_slice()

How to find the coordinates of a feature
""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::
    :raises:

    pseudogene.get_coordinates()

These are useful for doing custom things, e.g. If the introns are not annotated for a gene, we can generate the introns from the coordinates of the exons

.. jupyter-execute::
    :raises:

    cds = list(seq.get_features(biotype="CDS"))[0]
    exon_coords = cds.get_coordinates()
    
    exon_coords
    
We generate the intron coordinates from the second element of the first tuple, and first element of the second tuple.   

.. jupyter-execute::
    :raises:

    intron_coords = []
    
    for i in range(len(exon_coords)-1):
        intron_coords.append((exon_coords[i][1], exon_coords[i+1][0]))
        
    intron_coords

    
We can then generate a custom feature for the introns

.. jupyter-execute::
    :raises:

    seq.add_feature(biotype="intron", name="intron:Y74C9A.3.1", seqid="I", spans=intron_coords)

How to get the shadow of a Feature
""""""""""""""""""""""""""""""""""

We first need to query the Sequence for all genes. Using the ``union()`` method we combine all genes into a single feature. Finally, taking the "shadow" of all genes will return the intergenic region as a valid ``Feature``

.. jupyter-execute::
    :raises:

    genes = list(seq.get_features(biotype="gene"))
    genes = genes[0].union(genes[1:])
    intergenic = genes.shadow()
    intergenic

We can slice the sequence by this new Feautre

.. jupyter-execute::
    :raises:

    intergenic.get_slice()

Interrogating Features on a Sequence or Alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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



Methods on sequence to manipulate Features

``copy_annotations()``


``with_masked_annotations()``


``is_annotated()``


Methods on Features to manipulate Sequences


``without_lost_spans()``



``get_children()``

``get_parent()``


``get_drawable()``

``to_dict()``
