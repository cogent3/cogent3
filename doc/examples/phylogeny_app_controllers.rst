Using phylogeny application controllers to construct phylogenetic trees from alignments
=======================================================================================

.. sectionauthor:: Daniel McDonald

This document provides a few use case examples of how to use the phylogeny application controllers available in PyCogent. Each phylogeny application controller provides the support method ``build_tree_from_alignment``. This method takes as input an ``Alignment`` object, a ``SequenceColleciton`` object or a dict mapping sequence IDs to sequences. The ``MolType`` must also be specified. Optionally, you can indicate if you would like the "best_tree$", as well as any additional application parameters. These methods return a ``PhyloNode`` object.

To start, lets import all of our ``build_tree_from_alignment`` methods and our ``MolType``:

.. doctest::

    >>> from cogent.core.moltype import DNA
    >>> from cogent.app.clearcut import build_tree_from_alignment as clearcut_build_tree
    >>> from cogent.app.clustalw import build_tree_from_alignment as clustalw_build_tree
    >>> from cogent.app.fasttree import build_tree_from_alignment as fasttree_build_tree
    >>> from cogent.app.muscle import build_tree_from_alignment as muscle_build_tree
    >>> from cogent.app.raxml import build_tree_from_alignment as raxml_build_tree

Next, we'll load up a test set of sequences and construct an ``Alignment``:

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.app.muscle import align_unaligned_seqs
    >>> unaligned = LoadSeqs(filename='data/test2.fasta', aligned=False)
    >>> aln = align_unaligned_seqs(unaligned, DNA)

Now, let's construct some trees with default parameters!

.. note:: We are explicitly seeding Clearcut and RAxML to ensure reproducible results, and FastTree's output depends slightly on which version of FastTree is installed

.. doctest::

    >>> clearcut_tree = clearcut_build_tree(aln, DNA, params={'-s':42})
    >>> clustalw_tree = clustalw_build_tree(aln, DNA)
    >>> fasttree_tree = fasttree_build_tree(aln, DNA)
    >>> muscle_tree = muscle_build_tree(aln, DNA)
    >>> raxml_tree = raxml_build_tree(aln, DNA, params={'-p':42})
    >>> clearcut_tree
    Tree("(Mouse,(((HowlerMon,Human),DogFaced),NineBande));")
    >>> clustalw_tree
    Tree("((DogFaced,(HowlerMon,Human)),Mouse,NineBande);")
    >>> muscle_tree
    Tree("(Mouse,(DogFaced,(Human,(HowlerMon,NineBande))));")
    >>> raxml_tree
    Tree("((HowlerMon,Human),(DogFaced,Mouse),NineBande);")

.. code-block:: python

    >>> fasttree_tree
    Tree("(Mouse,NineBande,(DogFaced,(HowlerMon,Human)0.752)0.508);")

These methods allow the programmer to specify any of the applications parameters. Let's look at an example where we tell Clearcut to use traditional neighbor-joining, shuffle the distance matrix, use Kimura distance correction and explicitly seed the random number generator:

.. doctest::

    >>> clearcut_params = {'-N':True,'-k':True,'-S':True,'-s':42}
    >>> clearcut_tree = clearcut_build_tree(aln, DNA, params=clearcut_params)
    >>> clearcut_tree
    Tree("(((HowlerMon,Human),(NineBande,Mouse)),DogFaced);")

