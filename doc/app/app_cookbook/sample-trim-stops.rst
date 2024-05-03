Remove trailing stop codons from sequences in an alignment
----------------------------------------------------------

For evolutionary analyses that use codon models we need to exclude terminating stop codons. To account for the difference in stop codons in different genetic codes, we can provide an argument for the genetic code. For a list of all genetic codes, see :ref:`here <genetic-codes>`. 

.. jupyter-execute::
    :raises:

    from cogent3 import get_app, make_aligned_seqs

    aln = make_aligned_seqs(
        data={"seq1": "ACGTAA---", "seq2": "ACGACA---", "seq3": "ACGCAATGA"},
        moltype="dna",
    )
    no_stops = get_app("trim_stop_codons", gc=1)
    no_stops(aln)