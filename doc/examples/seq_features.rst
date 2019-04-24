.. _seq-annotations:

Advanced sequence handling
==========================

.. sectionauthor:: Gavin Huttley

Individual sequences and alignments can be manipulated by annotations. Most value in the genome sequences arises from sequence annotations regarding specific sequence feature types, e.g. genes with introns / exons, repeat sequences. These can be applied to an alignment either using data formats available from genome portals (e.g. GFF, or GenBank annotation formats) or by custom assignments.

Annotations can be added in two ways: using either the ``add_annotation`` or the ``add_feature`` method. The distinction between these two is that ``add_features`` is more specialised. Features can be thought of as a type of annotation representing standard sequence properties eg introns/exons. Annotations are the more general case, such as a computed property which has, say a numerical value and a span.

For illustrative purposes we define a sequence with 2 exons and grab the 1\ :sup:`st` \ exon:

.. doctest::

    >>> from cogent3 import DNA
    >>> s = DNA.make_seq("aagaagaagacccccaaaaaaaaaattttttttttaaaaaaaaaaaaa",
    ... name="Orig")
    >>> exon1 = s.add_feature('exon', 'exon1', [(10,15)])
    >>> exon2 = s.add_feature('exon', 'exon2', [(30,40)])

Here, '``exon``' is the feature type, and '``exon#``' the feature name. The feature type is used for the display formatting, which won't be illustrated here, and also for selecting all features of the same type, shown below.

We could also have created an annotation using the ``add_annotation`` method:

.. doctest::

    >>> from cogent3.core.annotation import Feature
    >>> s2=DNA.make_seq("aagaagaagacccccaaaaaaaaaattttttttttaaaaaaaaaaaaa",
    ... name="Orig2")
    >>> exon3 = s2.add_annotation(Feature, 'exon', 'exon1', [(35,40)])

We can use the features (eg ``exon1``) to get the corresponding sequence region.

.. doctest::

    >>> s[exon1]
    DnaSequence(CCCCC)

You can query annotations by type and optionally by label, receiving a list of features:

.. doctest::

    >>> exons = s.get_annotations_matching('exon')
    >>> print(exons)
    [exon "exon1" at [10:15]/48, exon "exon2" at [30:40]/48]

We can use this list to construct a pseudo-feature covering (or excluding) multiple features using ``get_region_covering_all``. For instance, getting all exons,

.. doctest::

    >>> print(s.get_region_covering_all(exons))
    region "exon" at [10:15, 30:40]/48
    >>> s.get_region_covering_all(exons).get_slice()
    DnaSequence(CCCCCTT... 15)

or not exons (the exon *shadow*):

.. doctest::

    >>> print(s.get_region_covering_all(exons).get_shadow().get_slice())
    AAGAAGAAGAAAAAAAAAAATTTTTAAAAAAAA

The first of these essentially returns the CDS of the gene.

Features are themselves sliceable:

.. doctest::

    >>> exon1[0:3].get_slice()
    DnaSequence(CCC)

This approach to sequence / alignment handling allows the user to manipulate them according to things they know about such as genes or repeat elements. Most of this annotation data can be obtained from genome portals.

The toolkit can perform standard sequence / alignment manipulations such as getting a subset of sequences or aligned columns, translating sequences, reading and writing standard formats.
