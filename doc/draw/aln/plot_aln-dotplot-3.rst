.. jupyter-execute::
    :hide-code:

    import set_working_directory

Assess alignment quality via dotplots
=====================================

The dotplot algorithm within ``cogent3`` offers a powerful tool to visualise the quality of a sequence alignment, as it shows both the matching segments between two sequences and the path chosen by the alignment algorithm. 

Visualising an alignment with poor parameter choices
----------------------------------------------------

Let's begin by importing the functions and loading sequences for this example. Note that although we load aligned sequences, we use the ``degap`` method to remove all gap characters from the sequences. This produces a set of unaligned sequences, for which we can align with different parameters to demonstrate the use of dotplots to visualise alignment quality. 

.. jupyter-execute::
    :raises:

    from cogent3 import load_aligned_seqs, get_app

    seqs = load_aligned_seqs("data/brca1.fasta", moltype="dna").degap()
    seqs = seqs.take_seqs(
        ["LeafNose", "FalseVamp", "RoundEare", "Sloth", "Anteater", "HairyArma"]
    )

The choice of parameters, such as gap insertion and gap extension parameters, can dramatically impact the quality of an alignment. An example of a **bad choice** of parameters is a high probability of opening a gap, but a low probability of extending a gap. 

Let's align using such parameters, and take a look at a dotplot between two of the sequences. 

.. jupyter-execute::
    :raises:

    aligner = get_app("progressive_align", "nucleotide", indel_rate=1e-2, indel_length=1e-9)
    aln = aligner(seqs)
    dp = aln[2200:2500].dotplot("HairyArma", "RoundEare")
    dp.show()

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-dotplot-3.png"

    dp.write(outpath)

The dotplot clearly shows the misalignment between the two sequences. Towards the second half of the alignment, there is a section of the alignment path (the dashed line), which does not line up with any matching segments between the sequences (the blue line). 

Taking a closer look at the poorly aligned section, we can see multiple small gaps in many of the sequences. The likelihood of this reflecting the true history of the sequences is low. It requires many indel events in many species and instead is likely to be an artefact of the alignment algorithm.

.. jupyter-execute::
    :raises:

    aln[2300:2500] 

Visualising an alignment with good parameter choices
----------------------------------------------------

Let's align the same sequences, but with more biologically realistic parameters. The default parameters for the progressive aligner are a better choice as they have a low probability of opening a gap, but a high probability of extending a gap.

.. note:: The default parameters (of any program) are not always the best choice, but they are a good starting point for many alignments.

.. jupyter-execute::
    :raises:

    aligner = get_app("progressive_align", "nucleotide")
    aln = aligner(seqs)
    aln[2200:2500].dotplot("HairyArma", "RoundEare").show()

The dotplot shows a much better alignment between the two sequences. The alignment path (the dashed line) now follows the matching segments between the sequences (the blue line). 

Taking a closer look at the same section, we can see a single large gap in two of the sequences. The likelihood of this reflecting the true history of the sequences is much higher than in the previous example, as it requires far fewer independent indel events to explain the alignment.

.. jupyter-execute::
    :raises:

    aln[2300:2500]
