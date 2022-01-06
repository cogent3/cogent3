.. jupyter-execute::
    :hide-code:

    import set_working_directory

Sequence logos
==============

Sequence logo's display sequence information. They're extensively applied to transcription factor binding site (TFBS) display. They can also be applied to sequence alignments more generally.

Drawing logo for a TFBS
-----------------------

We use the TFBS for the TATA box binding protein.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs
    from cogent3.parse import jaspar

    _, pwm = jaspar.read("data/tbp.jaspar")
    freqarr = pwm.to_freq_array()
    freqarr[:5]  # illustrating the contents of the MotifFreqsArray

.. jupyter-execute::

    logo = freqarr.logo()
    logo.show(height=250, width=500)

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-seqlogo.png"

    logo.write(outpath)

Drawing a sequence logo from a multiple sequence alignment
----------------------------------------------------------

This can be done for an entire alignment, but bear in mind it can take some time to render. Note that we include gap characters in the display.

.. jupyter-execute::

    aln = load_aligned_seqs("data/brca1-bats.fasta", moltype="dna")
    l = aln[:311].seqlogo(height=300, width=500, wrap=60, vspace=0.05)
    l.show()

Sequence logo of protein alignment
----------------------------------

No difference here except it uses the built-in colour scheme from the protein ``MolType``.

.. jupyter-execute::

    aa = aln.get_translation(incomplete_ok=True)[:120]
    logo = aa.seqlogo(width=500, height=300, wrap=50, vspace=0.1)
    logo.show()
