Estimate parameter values using a sampling from a dataset
=========================================================

This script uses the sample method of the alignment class to provide an
estimate for a two stage optimisation.
This allows rapid optimisation of long alignments and complex models with
a good chance of arriving at the global maximum for the model and data.
Local optimisation of the full alignment may end up in local maximum and
for this reason results from this strategy my be inaccurate.

From cogent import all the components we need.

.. doctest::

    >>> from cogent import LoadSeqs, LoadTree
    >>> from cogent.evolve import  substitution_model

Load your alignment, note that if your file ends with a suffix that is the same as it's format (assuming it's a supported format) then you can just give the filename. Otherwise you can specify the format using the format argument.

.. doctest::

    >>> al = LoadSeqs(filename = "data/test.paml")

Get your tree

.. doctest::

    >>> t = LoadTree(filename = "data/test.tree")

Get the raw substitution model

.. doctest::

    >>> sm = substitution_model.Nucleotide()

Make a likelihood function from a sample of the alignment the .sample method selects the chosen number of bases at random because we set motif probabilities earlier the motif probabilities of the whole alignment rather than the sample will be used for the calculator

.. doctest::

    >>> lf = sm.makeLikelihoodFunction(t)
    >>> lf.setMotifProbsFromData(al)
    >>> lf.setAlignment(al.sample(20))

Optimise with the slower but more accurate simulated annealing optimiser

.. doctest::

    >>> lf.optimise()
    Outer loop = 0...

Next use the whole alignment

.. doctest::

    >>> lf.setAlignment(al)

and the faster Powell optimiser that will only find the best result near the provided starting point

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> lf.optimise(local=True)
        Number of function evaluations = 1; current F = ...


Print the result using ``print lf``.
