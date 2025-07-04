.. jupyter-execute::
    :hide-code:

    import set_working_directory

Dotplot basics
==============

A technique (`Gibbs and McIntyre <https://www.ncbi.nlm.nih.gov/pubmed/5456129>`_) for comparing sequences. All ``cogent3`` sequence collections classes (``SequenceCollection``, ``Alignment``) have a dotplot method.

.. todo:: Change dotplot ref to a citation

The method returns a drawable, as demonstrated below between unaligned sequences.

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs

    seqs = load_unaligned_seqs("data/SCA1-cds.fasta", moltype="dna")
    draw = seqs.dotplot()
    draw.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-dotplot-1.png"

    draw.write(outpath)

If sequence names are not provided, two randomly chosen sequences are selected (see below). The plot title reflects the parameter values for defining a match. ``window`` is the size of the sequence segments being compared. ``threshold`` is the number of exact matches within ``window`` required for the two sequence segments to be considered a match. ``gap`` is the size of a gap between adjacent matches before merging.

Modifying the matching parameters
---------------------------------

If we set window and threshold to be equal, this is equivalent to an exact match approach.

.. jupyter-execute::

    draw = seqs.dotplot(name1="Human", name2="Mouse", window=8, threshold=8)
    draw.show()

Displaying dotplot for the reverse complement
---------------------------------------------

.. jupyter-execute::

    draw = seqs.dotplot(name1="Human", name2="Mouse", rc=True)
    draw.show()

.. note:: clicking on an entry in the legend turns it off

Setting plot attributes
-----------------------

I'll modify the title and figure width.

.. jupyter-execute::

    draw = seqs.dotplot(name1="Human", name2="Mouse", rc=True, title="SCA1", width=400)
    draw.show()
