.. jupyter-execute::
    :hide-code:

    import set_working_directory

*******************
Building alignments
*******************

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Using a ``cogent3`` progressive aligner for nucleotides
=======================================================

We load a canned nucleotide substitution model and the progressive aligner ``tree_align`` function.

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs, make_tree
    from cogent3.align.progressive import tree_align

We first align without providing a guide tree. The ``tree_align`` algorithm builds pairwise alignments and estimates the substitution model parameters and pairwise distances. The distances are used to build a neighbour joining tree and the median value of substitution model parameters are provided to the substitution model for the progressive alignment step.

.. jupyter-execute::

    seqs = load_unaligned_seqs("data/test2.fasta", moltype="dna")
    aln, tree = tree_align("HKY85", seqs, show_progress=False)
    aln

We then align using a guide tree (pre-estimated) and specifying the ratio of transitions to transversions (kappa).

.. jupyter-execute::

    tree = make_tree(
        "(((NineBande:0.013,Mouse:0.185):0.023,DogFaced:0.046):0.027,Human:0.034,HowlerMon:0.019)"
    )
    params = {"kappa": 4.0}
    aln, tree = tree_align(
        "HKY85", seqs, tree=tree, param_vals=params, show_progress=False
    )
    aln

Using a ``cogent3`` progressive aligner for codons
==================================================

We load a canned codon substitution model and use a pre-defined tree and parameter estimates.

.. jupyter-execute::

    from cogent3 import load_unaligned_seqs, make_tree
    from cogent3.align.progressive import tree_align

    seqs = load_unaligned_seqs("data/test2.fasta", moltype="dna")
    tree = make_tree(
        "((NineBande:0.058,Mouse:0.595):0.079,DogFaced:0.142,(HowlerMon:0.062,Human:0.103):0.079)"
    )
    params = {"kappa": 4.0, "omega": 1.3}
    aln, tree = tree_align(
        "MG94HKY", seqs, tree=tree, param_vals=params, show_progress=False
    )
    aln

Converting gaps from aa-seq alignment to nuc seq alignment
==========================================================

We load some unaligned DNA sequences and show their translation.

.. jupyter-execute::

    from cogent3 import make_unaligned_seqs
    from cogent3.align.progressive import tree_align
    from cogent3.evolve.models import get_model

    seqs = [
        (
            "hum",
            "AAGCAGATCCAGGAAAGCAGCGAGAATGGCAGCCTGGCCGCGCGCCAGGAGAGGCAGGCCCAGGTCAACCTCACT",
        ),
        (
            "mus",
            "AAGCAGATCCAGGAGAGCGGCGAGAGCGGCAGCCTGGCCGCGCGGCAGGAGAGGCAGGCCCAAGTCAACCTCACG",
        ),
        ("rat", "CTGAACAAGCAGCCACTTTCAAACAAGAAA"),
    ]
    unaligned_DNA = make_unaligned_seqs(seqs, moltype="dna")
    unaligned_DNA

.. jupyter-execute::

    unaligned_DNA.get_translation()

We load an alignment of these protein sequences.

.. jupyter-execute::

    from cogent3 import make_aligned_seqs

    aligned_aa_seqs = [
        ("hum", "KQIQESSENGSLAARQERQAQVNLT"),
        ("mus", "KQIQESGESGSLAARQERQAQVNLT"),
        ("rat", "LNKQ------PLS---------NKK"),
    ]
    aligned_aa = make_aligned_seqs(aligned_aa_seqs, moltype="protein")

We then obtain an alignment of the DNA sequences from the alignment of their translation.

.. jupyter-execute::

    aligned_DNA = aligned_aa.replace_seqs(unaligned_DNA, aa_to_codon=True)
    aligned_DNA
