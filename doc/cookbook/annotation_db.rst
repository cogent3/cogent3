.. jupyter-execute::
    :hide-code:

    import set_working_directory


Annotation Databases
--------------------

This guide shows you how to use ``cogent3``'s annotation databases, which are in-memory SQLite databases, to store, query and manipulate the features (also known as annotations) of one or more biological sequences.

For more extensive documentation about features see :ref:`intro_annotations` and :ref:`seq-annotations`.

What are the different types of ``AnnotationDb``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three types of databases that are available: ``BasicAnnotationDb``, ``GffAnnotationDb``, and ``GenbankAnnotationDb``. A ``BasicAnnotationDb`` provides a "user" table where custom features can be added. The latter two, as hinted in their names, are distinguished by the type of data file used to generate them. Both ``GffAnnotationDb`` and ``GenbankAnnotationDb`` have a "user" table for custom features, which allows them to be merged with a ``BasicAnnotationDb``. However, they cannot be merged with each other!

How to create a standalone ``BasicAnnotationDb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Achieved by creating an ``BasicAnnotationDb`` instance. This is an empty database to which we can add features.

.. jupyter-execute::
    :raises:

    from cogent3.core.annotation_db import BasicAnnotationDb

    anno_db = BasicAnnotationDb()
    anno_db
    
How to load an standalone ``AnnotationDb`` from a data file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically, we want to load a collection of features from a genomic annotation file, such as a GFF or Genbank file. For the following examples, we will use data from the bacterium *Mycoplasma genitalium*.

.. note:: See the list of :ref:`data_links` to download the data used in the following examples.

From a GFF file
"""""""""""""""

To load features from a GFF file, you can use the ``load_annotations`` function and provide the path to the GFF file. The function automatically determines the file type based on the file extension, and it returns a ``GffAnnotationDb`` object.

.. jupyter-execute::
    :raises:

    from cogent3 import load_annotations

    gff_db = load_annotations(path="data/mycoplasma-genitalium.gff")
    gff_db

From a Genbank file
"""""""""""""""""""

To load features from a Genbank file, you can once again use the ``load_annotations`` function and provide the path to the Genbank file. The function detects the file type based on the file extension, and it returns a ``GenbankAnnotationDb`` object.

.. jupyter-execute::
    :raises:

    from cogent3 import load_annotations

    gb_db = load_annotations(path="data/mycoplasma-genitalium.gb")
    gb_db

How to generate a summary of an ``AnnotationDb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To generate a summary of an ``AnnotationDb``, you can access the ``describe`` attribute of the database. This attribute returns a ``cogent3.util.table.Table`` instance that shows the number of records for each seqid, the count for each biotype, and the number of rows in each table (in this example there is a "gff" table with 1,169 rows and an empty "user" table).

.. jupyter-execute::
    :raises:

    summary = gff_db.describe
    summary

How to add custom features to an ``AnnotationDb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is achieved via the ``add_features`` method and for all three types of ``AnnotationDb`` it will be added to the "user" table. The method requires information about the feature, such as its biotype, name, genomic location (spans), and the seqid. The seqid is necessary when linking an ``AnnotationDb`` to a ``Sequence`` object, see :ref:`How to assign an AnnotationDb to a sequence <assign_db_to_seq>` for more information.

We can add a feature to the empty ``BasicAnnotationDb`` we created above. Now the database has one record!

.. jupyter-execute::
    :raises:

     anno_db.add_feature(
               seqid="NC_000908",
               biotype="gene",
               name="interesting_gene",
               spans=[(1, 4)],
               strand="+",
                )
    anno_db.describe

We can also add a feature to our ``GffAnnotationDb`` or ``GenbankAnnotationDb``. Below, the previously empty "user" table now has a row count of one, indicating that our feature has been successfully added to the database.

.. jupyter-execute::
    :raises:

    gff_db.add_feature(
        seqid="seq1",
        biotype="gene",
        name="interesting_gene",
        spans=[(1, 4)],
        strand="+",
    )
    gff_db.describe[-2:, :] # showing just last two rows

How to write an ``AnnotationDb`` to disk for efficient re-loading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the above examples, all databases indicate that ``source=":memory:"``, i.e. they are in-memory databases. We can write any database to disk using the ``write()`` method and providing an outpath.

.. code-block:: python

    # write to disk
    gb_db.write("data/m-genitalium-database.gbdb")

    # do something

    # re-load from disk
    quick_load_gb_db = GenbankAnnotationDb(source="data/m-genitalium-database.gbdb")

.. note:: The suffix of the outpath (".gbdb" in the above example) can be arbitrarily chosen, however, this behaviour may change in the future to only accept registered suffixes! 👀

How to query an ``AnnotationDb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note, there are two methods with the same interface available to query an ``AnnotationDb``:

1. ``get_features_matching()``. A generator that yields all features that matched the query. The **minimal information** required to create a ``cogent3`` ``Feature`` object is provided in the returned dictionary. For more information on Features see :ref:`seq-intro_annotations` and :ref:`seq-annotations`.

2. ``get_records_matching()``. A generator that yields all features that matched the query. The **complete record** for each matching feature is provided in the returned dictionary.

Put simply, a "feature" is a subset of a "record".

Querying via Feature Name
"""""""""""""""""""""""""

To query a database for a feature by its name, provide the name of the feature as an argument to either ``get_features_matching()`` or ``get_records_matching()``. Since an ``AnnotationDb`` can contain records for more than one sequence, it is best practice to also include the seqid of the sequence of interest.

For example, querying the ``GenbankAnnotationDb`` for the 16s rRNA gene:

.. jupyter-execute::
    :raises:

    mg_16s = list(
        gb_db.get_features_matching(
            name="MG_RS00775", biotype="gene", seqid="NC_000908"
        )
    )
    mg_16s

Querying via Feature Biotype
""""""""""""""""""""""""""""

Similarly, ``get_features_matching()`` and ``get_records_matching()`` can be used to query the database for all features that match a given biotype.

For example, querying the ``GffAnnotationDb`` for all pseudogenes:

.. jupyter-execute::
    :raises:

    pseudogenes = list(gff_db.get_features_matching(biotype="pseudogene"))
    pseudogenes[:2] # showing just the first two

Querying via region of interest
"""""""""""""""""""""""""""""""

We can provide ``start`` and ``end`` arguments to ``get_features_matching()`` and ``get_records_matching()`` and all features within the coordinates will be returned.

For example, the adhesin protein of *M. genitalium* is organised in an operon between positions 220600 to 229079, so we can query for genes in that region to return all operon genes:

.. jupyter-execute::
    :raises:

    operon_cds = list(
        gff_db.get_features_matching(start=220600, end=229067, biotype="CDS")
    )
    operon_cds

Querying via the extended attributes field
""""""""""""""""""""""""""""""""""""""""""

A particularly useful functionality of a ``GffAnnotationDb`` is the ability to search the extended attributes field. This allows you to query for records that have matches to a specific string provided to the ``attributes`` argument within their extended attributes field.

For example, you can query for all CDS related to replication:

.. jupyter-execute::
    :raises:

    replication_records = list(
        gff_db.get_records_matching(attributes="replication", biotype="CDS")
    )
    replication_records[0] # showing just the first match

.. note:: Extended attribute querying only works for GFF databases!

How to interrogate an ``AnnotationDb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An ``AnnotationDb`` can be interrogated to explore the properties of a sequence without needing the sequence information.

How many unique genes are in a given genome?
""""""""""""""""""""""""""""""""""""""""""""

*Mycoplasma genitalium* has the smallest bacterial genome, so the number of genes in the loaded database represents the approximate minimal set of genes required for bacterial life! We can see the total number of genes by using the ``num_matches()`` method and specifying the condition we want to be matched is that the biotype is "gene".

.. jupyter-execute::
    :raises:

    gb_db.num_matches(biotype="gene")

The count is 563, however, this may include genes with more than one copy. To determine the number of distinct genes we can use the ``count_distinct()`` method and specify ``biotype="gene"`` and ``name=True`` to indicate we are interested in genes with distinct names.

.. jupyter-execute::
    :raises:

    total_genes = gb_db.count_distinct(biotype="gene", name=True)
    single_copy = total_genes[total_genes.columns["count"] == 1, :]
    len(single_copy)

The count of unique genes is 561. This means that almost every gene is present only once in the genome, very little redundancy here!

Just for fun, let's try this with the GFF database... (downloaded from the exact same source)

.. jupyter-execute::
    :raises:

    total_genes = gff_db.num_matches(biotype="gene")
    print("total genes: ", total_genes)
    genes = gff_db.count_distinct(biotype="gene", name=True)
    single_copy = genes[genes.columns["count"] == 1, :]
    print("single copy genes: ", len(single_copy))

What? 🤯

How to find the "children" of a Feature
"""""""""""""""""""""""""""""""""""""""

Achieved via ``get_feature_children``

.. jupyter-execute::
    :raises:

    children = list(gff_db.get_feature_children(name="gene-MG_RS00035"))
    children

How to find the "parent" of a Feature
"""""""""""""""""""""""""""""""""""""

Achieved via ``get_feature_parent``

.. jupyter-execute::
    :raises:

    parents = list(gff_db.get_feature_parent(name="cds-WP_009885556.1"))
    parents

How to combine two ``AnnotationDb`` instances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Checking the compatibility of two ``AnnotationDb`` instances
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Combining data requires compatibility of the databases, this can be checked via the ``compatible()`` method. Below we check whether a ``GffAnnotationDb`` is compatible with a ``BasicAnnotationDb``.

.. jupyter-execute::
    :raises:

    gff_db.compatible(anno_db)

The method evaluates to ``True``, indicating that the data of the two databases can be merged.

What about merging a ``GffAnnotationDb`` and ``GenbankAnnotationDb``?

.. jupyter-execute::
    :raises:

    gff_db.compatible(gb_db)

The method evaluates to ``False``. Merging a ``GffAnnotationDb`` and ``GenbankAnnotationDb`` is not possible.

Taking the union of two ``AnnotationDb`` instances
""""""""""""""""""""""""""""""""""""""""""""""""""

The ``union()`` method will return a **new instance** with merged records.

.. jupyter-execute::
    :raises:

    union_db = gb_db.union(anno_db)
    union_db.describe[-2:, :]

In the new merged database, there is now content in both the "user" and "gff" table.

Updating an ``AnnotationDb`` with the record from another database
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The ``update()`` method will update records of a given database with another and return the **same instance** of the database.

.. jupyter-execute::
    :raises:

    gff_db.update(anno_db)
    gff_db.describe[-2:, :]

Initialise a ``AnnotationDb`` with another database
"""""""""""""""""""""""""""""""""""""""""""""""""""

You can assign a compatible database to the ``db`` argument in the ``AnnotationDb`` constructor. If its the same class, it's db will be bound to self and directly modified.

.. jupyter-execute::
    :raises:

    from cogent3.core.annotation_db import GenbankAnnotationDb
    
    new_gb_db = GenbankAnnotationDb(source="m-genitalium-database.gbdb", db=anno_db)
    new_gb_db

.. _assign_db_to_seq:

How to assign an ``AnnotationDb`` to a sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For more extensive documentation about annotations see :ref:`intro_annotations` and :ref:`seq-annotations`.

Directly assign an ``AnnotationDb`` to a Sequence
"""""""""""""""""""""""""""""""""""""""""""""""""

Assign the AnnotationDb to the ``annotation_db`` attribute of a Sequence

.. jupyter-execute::
    :raises:

    from cogent3 import make_seq

    seq1 = make_seq(
        "AAGAAGAAGACCCCCAAAAAAAAAATTTTTTTTTTAAAAAGGGAACCCT",
        name="NC_000908",
        moltype="dna",
    )

    seq1.annotation_db = anno_db
    seq1.annotation_db

Loading an ``AnnotationDb`` and ``Sequence`` using the ``load_seq()`` function
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For a single sequence Genbank file
++++++++++++++++++++++++++++++++++

Loading a sequence from a Genbank file will automatically create a database instance containing all features present in the file. This database instance will be bound to the ``Sequence`` instance via the ``.annotation_db`` attribute, accessing this attribute displays a representation of the bound annotations.

.. jupyter-execute::
    :raises:

    from cogent3 import load_seq

    gb_seq = load_seq("data/mycoplasma-genitalium.gb")
    gb_seq.annotation_db

.. note:: Only single sequence Genbank files are supported. To load multiple sequences with annotations, first load the sequences (using ``load_aligned_seqs`` or ``load_unaligned_seqs``), then annotate from a gff file.

For a single sequence FASTA file and an associated GFF annotation file
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Data can be loaded by providing the path to the gff file to the ``annotation_path`` argument of ``load_seq()``.

.. jupyter-execute::
    :raises:

    gff_seq = load_seq(
        "data/mycoplasma-genitalium.fa",
        annotation_path="data/mycoplasma-genitalium.gff",
    )
    gff_seq.annotation_db

.. note:: This assumes an exact match of the sequence name between files!

In the above example, the sequence name in the fasta file does not match any records in the gff3 file (it is ``"NC_000908.2 Mycoplasmoides genitalium G37, complete sequence"`` in the former, and ``"NC_000908.2"`` in the latter). However, if you are confident that they are related, then you can use the ``label_to_name`` argument of ``load_seq()`` to change the sequence name as follows:

.. jupyter-execute::
    :raises:

    seq = load_seq(
        "data/mycoplasma-genitalium.fa",
        annotation_path="data/mycoplasma-genitalium.gff",
        label_to_name=lambda x: x.split()[0],
    )
    seq.annotation_db
