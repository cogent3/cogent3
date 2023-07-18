.. jupyter-execute::
    :hide-code:

    import set_working_directory

Sequence Features
=================

A sequence "feature" is an annotated segment, with the annotations being generated either computationally, e.g. repeat classification, or experimentally, e.g. single nucleotide polymorphism. In this example, we just load both sequence and features from a GenBank record.

Drawing all features on a sequence segment
------------------------------------------

We load chromosome I of *Caenorhabditis elegans*.

.. jupyter-execute::

    from cogent3 import load_seq

    seq = load_seq("data/C-elegans-chromosome-I.gb", moltype="dna")
    seq

As you can see it's quite large. It doesn't make sense to try and display all the features, so we will slice it down to a 10kbp segment.

.. jupyter-execute::

    seq = seq[25000:35000]

Drawing features is then limited to features within that segment.

.. jupyter-execute::

    fig = seq.get_drawable()
    fig.show(height=400, width=700)

.. note:: If a feature extends outside the displayed segment, it's hover text indicates it as "(incomplete)".

Drawing selected feature biotypes
---------------------------------

We specify what biotypes we want to display.

.. jupyter-execute::

    fig = seq.get_drawable(biotype=("gene", "CDS", "mRNA"))
    fig.show(height=300, width=650)
