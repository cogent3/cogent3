Removing degenerate characters
------------------------------

Degenerate IUPAC base symbols represent a site position that can have multiple possible characters. For a DNA example, "Y" represents pyrimidines where the site can be either "C" or "T".

.. note:: In many molecular evolutionary and phylogenetic analyses, the gap character "-" is treated "N", meaning any base.

Let's create sample data with degenerate characters

.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs({"s1": "ACGA-GACG", "s2": "GATGATGYT"}, moltype="dna")
    aln

Omit aligned columns containing a degenerate character
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_degens = get_app("omit_degenerates", moltype="dna")
    result = omit_degens(aln)
    result


Omit all degenerate characters except gaps from an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we create the app with the argument ``gap_is_degen=False``, we can omit degenerate characters but retain gaps. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_degens_keep_gaps = get_app("omit_degenerates", moltype="dna", gap_is_degen=False)
    result = omit_degens_keep_gaps(aln)
    result

Omit k-mers which contain degenerate characters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we create ``omit_degenerates`` with the argument ``motif_length``, it will split sequences into non-overlapping tuples of the specified length and exclude any tuple that contains a degenerate character. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_degenerates_app = get_app("omit_degenerates", moltype="dna", motif_length=2)
    result = omit_degenerates_app(aln)

    result
