
Removing highly gapped positions
--------------------------------

Using the ``omit_gap_pos`` app, we can remove position from an alignment which exceed a specified proportions of gaps. 

Let's create a sample alignment with gaps. 
        
.. jupyter-execute::
    :raises:

    from cogent3 import make_aligned_seqs

    aln = make_aligned_seqs({"s1": "ACGA-GA-CG", "s2": "GATGATG-AT"}, moltype="dna")
    aln

Removing highly gapped nucleotide positions
"""""""""""""""""""""""""""""""""""""""""""

Sites with over 99% gaps are excluded by default.

.. jupyter-execute::
    :raises:

    from cogent3 import get_app

    omit_gap_pos_app = get_app("omit_gap_pos", moltype="dna")
    result = omit_gap_pos_app(aln)
    result

We can alter the threshold for the allowed fraction of gaps with the ``allowed_frac`` argument. Let's create an app that excludes all aligned sites with over 49% gaps.

.. jupyter-execute::
    :raises:

    omit_gap_pos_app = get_app("omit_gap_pos", allowed_frac=0.49, moltype="dna")
    result = omit_gap_pos_app(aln)
    result

Removing highly gapped codon positions
""""""""""""""""""""""""""""""""""""""

To eliminate any codon columns (where a column is a triple of nucleotides) that contain a gap character, we use the ``motif_length`` argument.

.. jupyter-execute::
    :raises:

    omit_gap_pos_app = get_app(
        "omit_gap_pos", allowed_frac=0, motif_length=3, moltype="dna"
    )
    result = omit_gap_pos_app(aln)
    result