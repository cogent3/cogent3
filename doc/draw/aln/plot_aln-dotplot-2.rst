.. jupyter-execute::
    :hide-code:

    import set_working_directory

Dotplot with annotated sequences
================================

.. note:: These docs now use the ``new_type`` core objects via the following setting.

    .. jupyter-execute::

        import os

        # using new types without requiring an explicit argument
        os.environ["COGENT3_NEW_TYPE"] = "1"

If sequences in a dotplot have been annotated, the ``dotplot()`` method returns an ``AnnotatedDrawable``.

Reloading sequences and annotations
-----------------------------------

The data file, ``tp53-annotations.json``, was created from a query of Ensembl for one-to-one orthologs of human TP53 between Human, Macaque, Orangutan and Marmoset. The resulting sequences were annotated with the location of the CDS for the canonical transcript, then the ``annotation_db`` was saved as json using ``annotation_db.to_json()``.

.. jupyter-execute::

    import cogent3

    ann_db = cogent3.load_annotations(
        path="data/tp53-annotations.json"
    )
    seqs = cogent3.load_unaligned_seqs(
        "data/tp53.fa", moltype="dna"
    )
    seqs.annotation_db = ann_db
    dp = seqs.dotplot(name1="Macaque", name2="Marmoset", width=600)
    dp.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-dotplot-2.png"

    dp.write(outpath)

Removing annotation tracks
--------------------------

.. jupyter-execute::

    help(dp.remove_track)

Thus we could remove the left annotation track, for instance with

.. jupyter-execute::

    dp.remove_track(left_track=True)
    dp.show()
