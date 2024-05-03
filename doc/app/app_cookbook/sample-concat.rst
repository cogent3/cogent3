Concatenating alignments
------------------------

The ``concat`` app provides a mechanism to concatenate alignments. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    concat_alns_app = get_app("concat", moltype="dna")

Let's create sample alignments with matching sequence names to use in the below examples. 

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln1 = make_aligned_seqs({"s1": "AAA", "s2": "CAA", "s3": "AAA"}, moltype="dna")
    aln2 = make_aligned_seqs({"s1": "GCG", "s2": "GGG", "s3": "GGT"}, moltype="dna")
    aln1

.. jupyter-execute::
    :raises:

    aln2


How to concatenate alignments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, sequences without matching names in the corresponding alignment are omitted (``intersect=True``).

.. jupyter-execute::
    :raises:

    result = concat_alns_app([aln1, aln2])
    result

How to concatenate alignments with missing sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By providing the argument ``intersect=False``, the ``concat`` app will includes missing sequences across alignments. Missing sequences are replaced by a sequence of ``"?"``.

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs, get_app

    concat_missing = get_app("concat", moltype="dna", intersect=False)
    aln3 = make_aligned_seqs({"s4": "GCG", "s5": "GGG"}, moltype="dna")
    result = concat_missing([aln1, aln3])
    result

How to concatenated alignments with a deliminator ``"N"``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can insert an ``"N"`` character in between the concatenated sequences. 

.. jupyter-execute::
    :raises:
    
    from cogent3 import get_app

    concat_delim = get_app("concat", join_seq="N", moltype="dna")
    result = concat_delim([aln1, aln2])
    result