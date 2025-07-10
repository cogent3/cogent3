Remove duplicated sequences from an alignment
---------------------------------------------

The ``omit_duplicated`` app removes redundant sequences from a sequence collection (aligned or unaligned).

Let's create sample data with duplicated sequences.

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    data = {
        "a": "ACGT",
        "b": "ACG-",  # identical to 'a' except has a gap
        "c": "ACGG",  # duplicate
        "d": "ACGG",  # duplicate
        "e": "AGTC",  # unique
    }
    aln = make_aligned_seqs(data, moltype="dna")
    aln

Creating the ``omit_duplicated`` app with the argument ``choose="longest"`` selects the duplicated sequence with the least number of gaps and ambiguous characters. In the above example, only one of ``c`` and ``d`` will be retained.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_shorter_duplicate = get_app("omit_duplicated", moltype="dna", choose="longest")
    omit_shorter_duplicate(aln)

Creating the ``omit_duplicated`` app with the argument  ``choose=None`` means only unique sequences are retained.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_all_duplicates = get_app("omit_duplicated", moltype="dna", choose=None)
    omit_all_duplicates(aln)

The ``mask_degen`` argument specifies how to treat matches between sequences with degenerate characters. 

Let's create sample data that has a DNA ambiguity code.

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs(
        {
            "s1": "ATCG",
            "s2": "ATYG",  # matches s1 with ambiguity
            "s3": "GGTA",
        },
        moltype="dna",
    )

Since "Y" represents pyrimidines where the site can be either "C" or "T", s1 indeed matches s2 and one of them will be removed. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    app_dna = get_app("omit_duplicated", mask_degen=True, choose="longest", moltype="dna")
    app_dna(aln)