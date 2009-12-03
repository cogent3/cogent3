Perform Principal Coordinates Analysis
======================================

.. sectionauthor:: Cathy Lozupone

An example of how to calculate the pairwise distances for a set of sequences.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.phylo import distance
    >>> from cogent.cluster.metric_scaling import PCoA

Import a substitution model (or create your own)

.. doctest::

    >>> from cogent.evolve.models import HKY85

Load the alignment.

.. doctest::

    >>> al = LoadSeqs("data/test.paml")

Create a pairwise distances object calculator for the alignment, providing a substitution model instance.

.. doctest::

    >>> d = distance.EstimateDistances(al, submodel= HKY85())
    >>> d.run(show_progress=False)

Now use this matrix to perform principal coordinates analysis.

.. doctest::

    >>> PCoA_result = PCoA(d.getPairwiseDistances())
    >>> print PCoA_result
    ======================================================================================
            Type              Label  vec_num-0  vec_num-1  vec_num-2  vec_num-3  vec_num-4
    --------------------------------------------------------------------------------------
    Eigenvectors          NineBande      -0.02       0.01      -0.04       0.01       0.00
    Eigenvectors           DogFaced      -0.04      -0.06       0.01       0.00       0.00
    Eigenvectors          HowlerMon      -0.07       0.01      -0.01      -0.02       0.00
    Eigenvectors              Mouse       0.20       0.01       0.01      -0.00       0.00
    Eigenvectors              Human      -0.07       0.04       0.03       0.01       0.00
     Eigenvalues        eigenvalues       0.05       0.01       0.00       0.00      -0.00
     Eigenvalues  var explained (%)      85.71       9.60       3.73       0.95      -0.00
    --------------------------------------------------------------------------------------

We can save these results to a file in a delimited format (we'll use tab here) that can be opened up in any data analysis program, like R or Excel. Here the principal coordinates can be plotted against each other for visualization.

.. doctest::

    >>> PCoA_result.writeToFile('PCoA_results.txt',sep='\t')

Cleanup
-------

We don't actually want to keep that file now, so I'm importing the ``os`` module to delete it.

.. doctest::

    >>> import os
    >>> os.remove('PCoA_results.txt')
