.. jupyter-execute::
    :hide-code:

    import set_working_directory

Dotplot with annotated sequences
================================

If sequences in a dotplot have been annotated, the `dotplot()` method returns an `AnnotatedDrawable`.

Reloading from json
-------------------

The data file, ``tp53.json``, was created from a query of ensembl for one-to-one orthologs of human TP53 between Human, Macaque, Orangutan and Marmoset. The resulting sequences were annotated with the location of the CDS for the canonical transcript, then the ``SequenceCollection`` was saved as json using ``cogent3.app.write_json``.

.. jupyter-execute::

    from cogent3.app.io import get_data_store, load_json

    loader = load_json()
    seqs = loader("data/tp53.json")
    dp = seqs.dotplot(name1="Macaque", name2="Marmoset", width=600)
    dp.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-dotplot-2.png"

    fig.write(outpath)

Removing annotation tracks
--------------------------

.. jupyter-execute::

    help(dp.remove_track)

Thus we could remove the left annotation track, for instance with

.. jupyter-execute::

    dp.remove_track(left_track=True)
