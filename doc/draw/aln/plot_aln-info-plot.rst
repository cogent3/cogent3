.. jupyter-execute::
    :hide-code:

    import set_working_directory

Information analysis of an alignment
====================================

Information here is in the formal sense -- maximum entropy minus the entropy at a position. This is fast to compute and is an indicator of the variability at a position.

Illustrated with a simple example
---------------------------------

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, make_aligned_seqs, make_seq

    s1 = make_seq("TGATGTAAGGTAGTT", name="s1", moltype="dna")
    s2 = make_seq("--CTGGAAGGGT---", name="s2", moltype="dna")

    seqs = make_aligned_seqs([s1, s2], array_align=False, moltype="dna")
    draw = seqs.information_plot(window=2, include_gap=True)
    draw.show(width=500, height=400)

On a sample data set
--------------------

Clicking on any of the legend items causes that to disappear from the plot.

.. jupyter-execute::

    aln = load_aligned_seqs("data/brca1.fasta", moltype="protein")

    fig = aln.information_plot(stat="median")
    fig.show(width=500, height=400)

.. jupyter-execute::
    :hide-code:

    outpath = set_working_directory.get_thumbnail_dir() / "plot_aln-info-plot.png"

    fig.write(outpath)
