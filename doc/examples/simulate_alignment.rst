Simulate an alignment
=====================

.. sectionauthor:: Gavin Huttley

How to  simulate an alignment. For this example we just create a simple model using a four taxon tree with very different branch lengths, a Felsenstein model with very different nucleotide frequencies and a long alignment.

See the other examples for how to define other substitution models.

.. doctest::

    >>> import sys
    >>> from cogent3 import make_tree
    >>> from cogent3.evolve import substitution_model

Specify the 4 taxon tree,

.. doctest::

    >>> t = make_tree('(a:0.4,b:0.3,(c:0.15,d:0.2)edge.0:0.1);')

Define our Felsenstein 1981 substitution model.

.. doctest::

    >>> sm = substitution_model.TimeReversibleNucleotide(motif_probs = {'A': 0.5, 'C': 0.2,
    ... 'G': 0.2, 'T': 0.1}, model_gaps=False)
    >>> lf = sm.make_likelihood_function(t)
    >>> lf.set_constant_lengths()
    >>> lf.set_name('F81 model')
    >>> print(lf)
    F81 model
    number of free parameters = 0
    ==========================
      edge    parent    length
    --------------------------
         a      root    0.4000
         b      root    0.3000
         c    edge.0    0.1500
         d    edge.0    0.2000
    edge.0      root    0.1000
    --------------------------
    ====================================
         A         C         G         T
    ------------------------------------
    0.5000    0.2000    0.2000    0.1000
    ------------------------------------

We'll now create a simulated alignment of length 1000 nucleotides.

.. doctest::

    >>> simulated = lf.simulate_alignment(sequence_length=1000)

The result is a normal ``Cogent`` alignment object, which can be used in the same way as any other alignment object.
