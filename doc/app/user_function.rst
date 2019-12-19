Custom composable apps
======================

You can make a simple customised app using the ``user_function`` app. This is a wrapper class that takes a reference to your function and the input, output and data types. The resulting app can then become part of a composed function.

Defining a ``user_function`` requires you consider four things.

``func``
    A function you have written. This is required.

``input_types``
    A type, or collection of type that your function can handle. This setting dictates what other apps have an output that is a compatable input for your function.

``output_types``
    A type, or collection of type that your function produces. This setting dictates what other apps can have yours as input.

``data_types``
    The data class names, as strings, that your function can handle. Not required, but useful.

A simple example
----------------

We make a very simple function ``first4``, that returns the first 4 elements of an alignment.

.. doctest::

    >>> def first4(val):
    ...     return val[:4]

Now we define a ``user_function`` instance that takes and returns an ``ALIGNED_TYPE``.

.. doctest::

    >>> from cogent3.app.composable import user_function, ALIGNED_TYPE
    >>> just4 = user_function(first4, input_types=ALIGNED_TYPE, output_types=ALIGNED_TYPE, data_types="Alignment")

The ``repr()`` of your ``user_function`` instance indicates the wrapped function and the module it's in.

.. doctest::

    >>> just4
    user_function(name='first4', module='__main__')

You use it like all composable apps which we demonstrate using a small sample alignment.

.. doctest::

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data=dict(a="GCAAGCGTTTAT", b="GCTTTTGTCAAT"), array_align=False)
    >>> result = just4(aln)
    >>> result
    2 x 4 text alignment: a[GCAA], b[GCTT]

Renaming sequences
------------------

This time we wrap a method call on a ``SequenceCollection`` (and the alignment sub-classes) for renaming sequences. We also illustrate here that to support both aligned and unaligned data types as input/output, we have to include these in the construction of the custom function.

.. note:: The ``SERIALISABLE_TYPE`` indicates the data has the ability to be converted to ``json``.

.. doctest::
    
    >>> from cogent3.app.composable import user_function, ALIGNED_TYPE, SEQUENCE_TYPE, SERIALISABLE_TYPE
    >>> def renamer(aln):
    ...     """upper case names"""
    ...     return aln.rename_seqs(lambda x: x.upper())
    >>> rename_seqs = user_function(renamer, 
    ...                             input_types=(ALIGNED_TYPE, SEQUENCE_TYPE),
    ...                             output_types=SERIALISABLE_TYPE,
    ...                             data_types=("SequenceCollection", "Alignment", "ArrayAlignment"))
    >>> result = rename_seqs(aln)
    >>> result.names
    ['A', 'B']

A user function for with a different output type
------------------------------------------------

In this example, we make an function that returns ``DistanceMatrix`` of an alignment.

.. doctest::

    >>> from cogent3.app.composable import user_function, ALIGNED_TYPE, PAIRWISE_DISTANCE_TYPE
    >>> def _get_dist(aln):
    ...     return aln.distance_matrix(calc="hamming", show_progress=False)
    >>> get_dist = user_function(_get_dist, input_types=ALIGNED_TYPE,
    ...                          output_types=PAIRWISE_DISTANCE_TYPE,
    ...                          data_types=("Alignment", "ArrayAlignment"))
    >>> result = get_dist(aln)
    >>> result
    =====================
              a         b
    ---------------------
    a    0.0000    6.0000
    b    6.0000    0.0000
    ---------------------
