Simulate an alignment
=====================

.. sectionauthor:: Gavin Huttley

For this example we just create a simple model using a four taxon tree with different branch lengths and a Felsenstein model.

.. jupyter-execute::

    import sys

    from cogent3 import make_tree
    from cogent3.evolve.models import get_model

Specify the 4 taxon tree,

.. jupyter-execute::

    t = make_tree("(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1);")

Define our Felsenstein 1981 substitution model.

.. jupyter-execute::

    sm = get_model("F81")
    lf = sm.make_likelihood_function(t)
    lf.set_constant_lengths()
    lf.set_motif_probs(dict(A=0.1, C=0.2, G=0.3, T=0.4))
    lf

We'll now create a simulated alignment of length 1000 nucleotides.

.. jupyter-execute::

    simulated = lf.simulate_alignment(sequence_length=1000)
    simulated
