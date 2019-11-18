######################
Custom composable apps
######################

Overview
--------

``user_function`` is a wrapper class for user specified function

Configuration
-------------
``user_function`` is supposed to wrap a function with associated ``input_types``, ``output_types`` and ``data_types``

:``input_types``: Allowed input types

                  is a String or a Collection of String
:``output_types``: The types of output

                  is a String or a Collection of String
:``data_types``: Allowed data types

                  is a String or a Collection of String

                  the default value is ``None``

Example:

.. code-block:: python

    from cogent3.app.composable import user_function
    u_function = user_function(func, input_types="allowed_input_types", output_types="corresponding_output_types", data_types="allowed_data_types")

Usages
----------------------

Example 1
~~~~~~~~~

We first make a very simple function ``foo``.

.. doctest::

    >>> def foo(self, val, *args, **kwargs):
    ...     return val[:4]

Now we define ``user_function``:

.. code-block:: python

    >>> from cogent3.app.composable import user_function
    >>> u_function = user_function(foo, input_types="aligned", output_types="aligned")

Here we specify that the user specified function ``foo`` will have ``aligned`` input type and ``aligned`` output type

Next, we make an ``ArrayAlignment`` data type by invoking function ``make_aligned_seqs`` and feed it into defined ``u_function``

.. code-block:: python

    >>> from cogent3 import make_aligned_seqs
    >>> aln = make_aligned_seqs(data=[("a", "GCAAGCGTTTAT"), ("b", "GCTTTTGTCAAT")])
    >>> result = u_function(aln).to_dict()
    >>> result

    {"a": "GCAA", "b": "GCTT"}

Example 2
~~~~~~~~~

In this example, we make another function ``bar`` that returns DistanceMatrix of an alignment

.. doctest::

    >>> def bar(self, val, *args, **kwargs):
    ...     return return val.distance_matrix(calc="hamming", show_progress=False)

Noticeably, function ``bar`` is supposed to have ``"aligned"`` as its input type and ``"pairwise_distances"`` as its output type, and we will be feeding it an ``Alignment`` data type

.. code-block:: python

    >>> from cogent3.app.composable import user_function
    >>> from cogent3.core.alignment import Alignment
    >>> data = dict([("s1", "ACGTACGTA"), ("s2", "GTGTACGTA")])
    >>> aln = Alignment(data=data, moltype="dna")
    >>> u_function = user_function(bar, input_types="aligned", output_types="pairwise_distances", data_types="Alignment")
    >>> result = u_function(aln)
    >>> result

    {("s1", "s2"): 2.0, ("s2", "s1"): 2.0}

Example 3
~~~~~~~~~

This time we wrap a function ``rename_seqs`` which takes the responsibility of capitalising all the names of sequences in a ``SequenceCollection``

Define function rename_seqs:

.. code-block:: python

        >>> _renamer = lambda x: x.upper()
        >>> def rename_seqs(seqs_collection):
        ...     new = {}
        ...     for name, seq in seqs_collection.named_seqs.items():
        ...         new_name = _renamer(name)
        ...         new_seq = seqs_collection.moltype.make_seq(seq, new_name)
        ...         new[new_name] = new_seq
        ...     result = seqs_collection.__class__(data=new, info=seqs_collection.info, moltype=seqs.moltype)
        ...     return result


Next, we make a ``SequenceCollection``

.. code-block:: python

    >>> from cogent3.core.alignment import SequenceCollection
    >>> from cogent3.core.moltype import DNA
    >>> data = {"seq1": "ACGTACGTA", "seq2": "ACCGAA---", "seq3": "ACGTACGTT"}
    >>> seqs = SequenceCollection(data, moltype=DNA)

We then feed this ``SequenceCollection`` into a defined ``user_function`` wrapping the user specified function ``rename_seqs`` and return the result

.. code-block:: python

    >>> from cogent3.app.composable import user_function
    >>> u_function = user_function(rename_seqs, input_types="aligned", output_types="aligned", data_types="SequenceCollection")
    >>> result = u_function(seqs)

    >SEQ1\nACGTACGTA\n>SEQ2\nACCGAA---\n>SEQ3\nACGTACGTT\n