*******************
Building alignments
*******************

.. authors, Gavin Huttley, Kristian Rother, Patrick Yannul

Using the cogent3 aligners
==========================

Running a progressive aligner
-----------------------------

We import useful functions and then load the sequences to be aligned.

.. jupyter-execute::
    :linenos:

    from cogent3 import load_unaligned_seqs, make_tree

    seqs = load_unaligned_seqs("data/test2.fasta", moltype="dna")

For nucleotides
^^^^^^^^^^^^^^^

We load a canned nucleotide substitution model and the progressive aligner ``TreeAlign`` function.

.. jupyter-execute::
    :linenos:

    from cogent3.evolve.models import HKY85
    from cogent3.align.progressive import TreeAlign

We first align without providing a guide tree. The ``TreeAlign`` algorithm builds pairwise alignments and estimates the substitution model parameters and pairwise distances. The distances are used to build a neighbour joining tree and the median value of substitution model parameters are provided to the substitution model for the progressive alignment step.

.. jupyter-execute::
    :linenos:

    aln, tree = TreeAlign(HKY85(), seqs)
    aln

We then align using a guide tree (pre-estimated) and specifying the ratio of transitions to transversions (kappa).

.. jupyter-execute::
    :linenos:

    tree = make_tree(
        "(((NineBande:0.0128202449453,Mouse:0.184732725695):0.0289459522137,DogFaced:0.0456427810916):0.0271363715538,Human:0.0341320714654,HowlerMon:0.0188456837006)root;"
    )
    params = {"kappa": 4.0}
    aln, tree = TreeAlign(HKY85(), seqs, tree=tree, param_vals=params)
    aln

For codons
^^^^^^^^^^

We load a canned codon substitution model and use a pre-defined tree and parameter estimates.

.. jupyter-execute::
    :linenos:

    from cogent3.evolve.models import MG94HKY

    tree = make_tree(
        "((NineBande:0.0575781680031,Mouse:0.594704139406):0.078919659556,DogFaced:0.142151930069,(HowlerMon:0.0619991555435,Human:0.10343006422):0.0792423439112)"
    )
    params = {"kappa": 4.0, "omega": 1.3}
    aln, tree = TreeAlign(MG94HKY(), seqs, tree=tree, param_vals=params)
    aln

Converting gaps from aa-seq alignment to nuc seq alignment
==========================================================

We load some unaligned DNA sequences and show their translation.

.. jupyter-execute::
    :linenos:

    from cogent3 import make_unaligned_seqs

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
    print(unaligned_DNA)
    print(unaligned_DNA.get_translation())

We load an alignment of these protein sequences.

.. jupyter-execute::
    :linenos:

    from cogent3 import make_aligned_seqs

    aligned_aa_seqs = [
        ("hum", "KQIQESSENGSLAARQERQAQVNLT"),
        ("mus", "KQIQESGESGSLAARQERQAQVNLT"),
        ("rat", "LNKQ------PLS---------NKK"),
    ]
    aligned_aa = make_aligned_seqs(aligned_aa_seqs, moltype="protein")

We then obtain an alignment of the DNA sequences from the alignment of their translation.

.. jupyter-execute::
    :linenos:

    aligned_DNA = aligned_aa.replace_seqs(unaligned_DNA, aa_to_codon=True)
    print(aligned_DNA)
