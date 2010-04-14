Seqsim Simple Tree Simulation
=============================

.. sectionauthor:: Julia Goodrich

This is an example of how to use the birth-death model in cogent to simulate
a tree.

.. doctest::

    >>> from cogent.seqsim.birth_death import BirthDeathModel, ExtinctionError,\
    ... TooManyTaxaError

Create a model with specific death probabilities per timestep using 
``BirthProb`` and ``DeathProb``, and ``TimePerStep``. The desired maximum 
number of taxa on the tree can be set using ``MaxTaxa``.


.. doctest::

    >>> b = BirthDeathModel(BirthProb=0.2, DeathProb=0.1, TimePerStep=0.03, 
    ... MaxTaxa=20)

To simulate a tree with an exact number of taxa, use the following loop. 
The ``exact`` flag raises a ``TooManyTaxaError`` exception if the call produces
the wrong number of taxa (e.g. because too many things died in the same 
timestep). An ``ExtinctionError`` is raised if the ancestral node dies off, or 
all the nodes die off.


.. doctest::

    >>> while True:
    ...     try:
    ...         t = b(exact=True)
    ...     except (ExtinctionError, TooManyTaxaError):
    ...         pass
    ...     else:
    ...         break

t now contains a ``RangeNode`` tree object that was returned, and can be used 
the same ways any other tree object in cogent is used. For instance, ``seqsim`` 
can be used to evolve alignments.
