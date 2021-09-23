Turn your functions into composable apps
========================================

This is super easy -- just use the ``appify`` decorator! This generates a ``user_function`` wrapper class that takes a reference to your function and the input, output and data types. The resulting app can then become part of a composed function.

You need four things.

``func``
    A function to decorate ... duh!

``input_types``
    A type, or collection of type that your function can handle. This setting dictates what other apps have an output that is a compatable input for your function.

``output_types``
    A type, or collection of type that your function produces. This setting dictates what other apps can have yours as input.

``data_types``
    The data class names, as strings, that your function can handle. Not required, but useful.

A simple example
----------------

Let's make an app that returns the elements of an alignment up to a specified index, with the index being a keyword argument. We now define a decorated function ``up_to()``

.. jupyter-execute::

    from cogent3.app.composable import ALIGNED_TYPE, appify

    @appify(ALIGNED_TYPE, ALIGNED_TYPE, data_types="Alignment")
    def up_to(val, index=4):
        return val[:index]

Now we define a ``user_function`` instance that takes and ret
The ``repr()`` of your ``user_function`` instance indicates the wrapped function and the module it's in.

.. jupyter-execute::

    up_to

We create an app instance for a specific value of ``index``

.. jupyter-execute::

    first4 = up_to(index=4)
    first4

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

.. note:: The ``SERIALISABLE_TYPE`` indicates the data has the ability to be converted to ``json``.

.. jupyter-execute::

    from cogent3.app.composable import (
        ALIGNED_TYPE,
        SEQUENCE_TYPE,
        SERIALISABLE_TYPE,
        appify,
    )

    @appify((ALIGNED_TYPE, SEQUENCE_TYPE), SERIALISABLE_TYPE)
    def rename_seqs(aln):
        """upper case names"""
        return aln.rename_seqs(lambda x: x.upper())

    renamer = rename_seqs()
    result = renamer(aln)
    result

A user app with a different output type
---------------------------------------

In this example, we make an function that returns ``DistanceMatrix`` of an alignment.

.. jupyter-execute::

    from cogent3.app.composable import (
        ALIGNED_TYPE,
        PAIRWISE_DISTANCE_TYPE,
        appify,
    )

    @appify(ALIGNED_TYPE, PAIRWISE_DISTANCE_TYPE)
    def get_dists(aln, calc="hamming"):
        return aln.distance_matrix(calc=calc, show_progress=False)

    percent_dist = get_dists(calc="percent")
    result = percent_dist(aln)
    result

.. note:: We omitted the ``data_types`` argument just for demonstration purposes.
