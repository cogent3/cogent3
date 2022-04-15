.. jupyter-execute::
    :hide-code:

    import set_working_directory

The simplest script
===================

.. sectionauthor:: Gavin Huttley

We use a canned nucleotide substitution model (the ``HKY85`` model) on just three primate species. As there is only one unrooted tree possible, the sequence names are all that's required to make the tree.

.. jupyter-execute::

    from cogent3 import load_aligned_seqs, make_tree
    from cogent3.evolve.models import get_model

    model = get_model("HKY85")
    aln = load_aligned_seqs("data/primate_cdx2_promoter.fasta")
    tree = make_tree(tip_names=aln.names)
    lf = model.make_likelihood_function(tree)
    lf.set_alignment(aln)
    lf.optimise(show_progress=False)
    lf
