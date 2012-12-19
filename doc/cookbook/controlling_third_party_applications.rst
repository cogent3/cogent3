************************************
Controlling third party applications
************************************

Existing supported apps
=======================

Alignment apps
--------------

clustalw and muscle
^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> from cogent import LoadSeqs, DNA
    >>> from cogent.app.clustalw import align_unaligned_seqs as clustal_aln
    >>> from cogent.app.muscle import align_unaligned_seqs as muscle_aln
    >>> seqs = LoadSeqs(filename='data/test2.fasta', aligned=False)
    >>> aln1 = clustal_aln(seqs, DNA)
    >>> aln2 = muscle_aln(seqs, DNA)
    >>> aln1 == aln2
    True
    >>> from cogent.app.fasttree import build_tree_from_alignment
    >>> tr = build_tree_from_alignment(aln1,moltype=DNA)
    >>> print tr.asciiArt()
              /-Mouse
             |
    ---------|--NineBande
             |
             |          /-DogFaced
              \0.508---|
                       |          /-HowlerMon
                        \0.752---|
                                  \-Human

And if you have matplotlib installed you can draw the tree (see :ref:`draw-trees`).

.. note:: Tree output based on v2.0.1
.. TODO add in cross-ref to drawing usage example

BLAST
-----

See :ref:`blast-usage`.

Phylo apps
----------

*To be written.*

RNA structure apps
------------------

*To be written.*

Visualisation
-------------

*To be written.*

PyMol, MAGE
^^^^^^^^^^^

Adding new apps
===============

*To be written.*

Pipelining 3rd party apps
=========================

*To be written.*

.. integrating with cogent features

.. grab seqs from genbank, align, build tree, cogent evolutionary analysis

