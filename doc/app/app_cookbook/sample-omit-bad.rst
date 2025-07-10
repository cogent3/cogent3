Remove problem sequences from an alignment
------------------------------------------

Using ``omit_bad_seqs`` we can eliminate sequences from an ``Alignment`` based on their gap fraction and/or the number of gaps they uniquely introduce. 

Let's create a sample alignment with some gaps. 

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        {
            "s1": "---ACC---TT-",
            "s2": "---ACC---TT-",
            "s3": "---ACC---TT-",
            "s4": "--AACCG-GTT-",
            "s5": "--AACCGGGTTT",
            "s6": "AGAACCGGGTT-",
            "s7": "------------",
        },
        moltype="dna",
    )

Removing sequences with more than X% gaps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating the ``omit_bad_seqs`` app with the argument ``gap_fraction=0.5`` will omit sequences that contain 50% or more gaps.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_frac_05 = get_app("omit_bad_seqs", gap_fraction=0.5)
    omit_frac_05(aln)

Removing sequences that contribute many gaps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``quantile=0.8`` argument omits sequences that are ranked above the specified quantile with respect to the number of gaps uniquely introduced into the alignment. In the following example, sequence ``s6`` is omitted, as it uniquely introduces gaps in the first two positions of the alignment.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_quant_08 = get_app("omit_bad_seqs", quantile=0.8)
    omit_quant_08(aln)
