.. jupyter-execute::
    :hide-code:

    import set_working_directory

Annotation DB Features
======================

Genomic annotations can be visualised directly from an annotation database, without needing a sequence object. This is useful for exploring annotations loaded from GFF or GenBank files.

Drawing features from an annotation database
---------------------------------------------

We load annotations for *Caenorhabditis elegans* chromosome I from a GFF file.

.. jupyter-execute::

    from cogent3 import load_annotations, draw_annotations

    db = load_annotations(path="data/C-elegans-chromosome-I.gff")
    db

Drawing all features for a region
----------------------------------

We can draw features for a specific region by providing coordinate bounds.

.. jupyter-execute::

    fig = draw_annotations(db, seqid="I", biotype=("gene", "CDS"), start=25000, stop=35000)
    fig.show(height=400, width=700)

.. note:: This produces the same visual output as ``seq.get_drawable(biotype=...)``, but without loading the sequence data.

Using the range slider
-----------------------

By default, figures include an interactive range slider below the main plot. Drag the handles or select a region on the slider to zoom into a portion of the sequence. The x-axis label indicates the slider.

.. note:: Hover over filled shapes to see feature details (name, biotype, coordinates, strand).

Filtering by biotype
---------------------

We can select specific feature types to display.

.. jupyter-execute::

    fig = draw_annotations(db, seqid="I", biotype="gene", start=25000, stop=35000)
    fig.show(height=300, width=650)

Disabling the range slider for static images
----------------------------------------------

The range slider is included when saving figures as images. To produce a clean static image, disable the controls with ``show_controls=False``.

.. jupyter-execute::

    fig = draw_annotations(db, seqid="I", biotype="gene", start=25000, stop=35000, show_controls=False)
    fig.show(height=300, width=650)

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_annotation-db.png"

    fig.write(outpath)
