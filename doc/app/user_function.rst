Turn your functions into composable apps
========================================

This is super easy -- just use the ``define_app`` decorator! This generates a wrapper class that has a reference to your functionand can then become part of a composed function.

You need two things, your function and type hints on the first argument and the function return type.

A simple example
----------------

Let's make an app that returns the elements of an alignment up to a specified index, with the index being a keyword argument. We now define a decorated function ``up_to`` and import the type hints and the decorator.

.. jupyter-execute::

    from cogent3.app.composable import define_app
    from cogent3.app.typing import AlignedSeqsType

    @define_app
    def up_to(val: AlignedSeqsType, index=2) -> AlignedSeqsType:
        return val[:index]

We create an app instance for a specific value of ``index``

.. jupyter-execute::

    first4 = up_to(index=4)
    first4

The ``repr()`` of the instance indicates the wrapped function and the module it's in.

You use ``first4()`` like all composable apps, e.g.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        data=dict(a="GCAAGCGTTTAT", b="GCTTTTGTCAAT"), array_align=False, moltype="dna"
    )
    result = first4(aln)
    result

Renaming sequences
------------------

This time we wrap a method call on a ``SequenceCollection`` (and the alignment sub-classes) for renaming sequences. We also illustrate here that to support both aligned and unaligned data types as input/output, we have to include these in the construction of the custom function.

.. note:: The ``SerialisableType`` indicates the data has the ability to be converted to ``json``.

.. jupyter-execute::
    
    from typing import Union

    from cogent3.app.composable import define_app
    from cogent3.app.typing import SeqsCollectionType, SerialisableType
    
    T = Union[SeqsCollectionType, SerialisableType]
    
    @define_app
    def rename_seqs(seqs: SeqsCollectionType) -> T:
        """upper case names"""
        return seqs.rename_seqs(lambda x: x.upper())

    renamer = rename_seqs()
    result = renamer(aln)
    result

A user app with a different output type
---------------------------------------

In this example, we make an function that returns a ``DistanceMatrix`` from an alignment.

.. jupyter-execute::

    from cogent3.app.composable import define_app
    from cogent3.app.typing import AlignedSeqsType, PairwiseDistanceType

    @define_app
    def get_dists(aln: AlignedSeqsType, calc="hamming") -> PairwiseDistanceType:
        return aln.distance_matrix(calc=calc, show_progress=False)

    percent_dist = get_dists(calc="percent")
    result = percent_dist(aln)
    result
