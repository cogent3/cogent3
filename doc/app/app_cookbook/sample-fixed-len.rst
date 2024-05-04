Sample an alignment to a fixed length
-------------------------------------

Let's load in an alignment of rodents to use in the examples. 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app

    loader = get_app("load_aligned", moltype="protein", format="phylip")
    aln = loader("data/abglobin_aa.phylip")
    aln

How to sample the first ``n`` positions of an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can use the ``fixed_length`` app to sample an alignment to a fixed length. By default, it will sample from the beginning of an alignment, the argument ``length=20`` specifies how many positions to sample. 


.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    first_20 = get_app("fixed_length", length=20)
    first_20(aln)


How to sample ``n`` positions from within an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating the ``fixed_length`` app with the argument ``start=x`` specifies that the sampled sequence should begin ``x`` positions into the alignment. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    skip_10_take_20 = get_app("fixed_length", length=20, start=10)
    skip_10_take_20(aln)


How to sample ``n`` positions randomly from within an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The start position can be selected at random with ``random=True``. An optional ``seed`` can be provided to ensure the same start position is used when the app is called.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    random_20 = get_app("fixed_length", length=20, random=True)
    random_20(aln)
