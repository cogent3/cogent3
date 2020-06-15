Simulate an alignment
=====================

.. sectionauthor:: Gavin Huttley

How to  simulate an alignment. For this example we just create a simple model using a four taxon tree with very different branch lengths, a Felsenstein model with very different nucleotide frequencies and a long alignment.

See the other examples for how to define other substitution models.

.. jupyter-execute::
    :linenos:

    import sys
    from cogent3 import make_tree
    from cogent3.evolve import substitution_model

Specify the 4 taxon tree,

.. jupyter-execute::
    :linenos:

    t = make_tree("(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1);")

Define our Felsenstein 1981 substitution model.

.. jupyter-execute::
    :linenos:

    sm = substitution_model.TimeReversibleNucleotide(
        motif_probs={"A": 0.5, "C": 0.2, "G": 0.2, "T": 0.1}, model_gaps=False
    )
    lf = sm.make_likelihood_function(t)
    lf.set_constant_lengths()
    lf.set_name("F81 model")
    print(lf)

We'll now create a simulated alignment of length 1000 nucleotides.

.. jupyter-execute::
    :linenos:

    simulated = lf.simulate_alignment(sequence_length=1000)
