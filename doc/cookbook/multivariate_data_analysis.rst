**************************
Multivariate data analysis
**************************

Principal Coordinates Analysis
==============================

Principal Coordinates Analysis works on a matrix of pairwise distances. In this example we start by calculating the pairwise distances for a set of aligned sequences.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.phylo import distance
    >>> from cogent.cluster.metric_scaling import PCoA

Import a substitution model (or create your own).

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
    Eigenvectors          NineBande      -0.02       0.01       0.04       0.01       0.00
    Eigenvectors           DogFaced      -0.04      -0.06      -0.01       0.00       0.00
    Eigenvectors          HowlerMon      -0.07       0.01       0.01      -0.02       0.00
    Eigenvectors              Mouse       0.20       0.01      -0.01      -0.00       0.00
    Eigenvectors              Human      -0.07       0.04      -0.03       0.01       0.00
     Eigenvalues        eigenvalues       0.05       0.01       0.00       0.00      -0.00
     Eigenvalues  var explained (%)      85.71       9.60       3.73       0.95      -0.00
    --------------------------------------------------------------------------------------

We can save these results to a file in a delimited format (we'll use tab here) that can be opened up in any data analysis program, like R or Excel. Here the principal coordinates can be plotted against each other for visualization.

.. doctest::

    >>> PCoA_result.writeToFile('PCoA_results.txt',sep='\t')

NMDS
====
NMDS (Non-metric MultiDimensional Scaling) works on a matrix of pairwise distances. In this example, we generate a matrix based on the euclidean distances of an abundance matrix.

.. doctest::

    >>> from cogent.cluster.nmds import NMDS
    >>> from cogent.maths.distance_transform import dist_euclidean
    >>> from numpy import array

We start with an abundance matrix, samples (rows) by sequences/species (cols)

.. doctest::

    >>> abundance = array(
    ...        [[7,1,0,0,0,0,0,0,0],
    ...        [4,2,0,0,0,1,0,0,0],
    ...        [2,4,0,0,0,1,0,0,0],
    ...        [1,7,0,0,0,0,0,0,0],
    ...        [0,8,0,0,0,0,0,0,0],
    ...        [0,7,1,0,0,0,0,0,0],#idx 5
    ...        [0,4,2,0,0,0,2,0,0],
    ...        [0,2,4,0,0,0,1,0,0],
    ...        [0,1,7,0,0,0,0,0,0],
    ...        [0,0,8,0,0,0,0,0,0],
    ...        [0,0,7,1,0,0,0,0,0],#idx 10
    ...        [0,0,4,2,0,0,0,3,0],
    ...        [0,0,2,4,0,0,0,1,0],
    ...        [0,0,1,7,0,0,0,0,0],
    ...        [0,0,0,8,0,0,0,0,0],
    ...        [0,0,0,7,1,0,0,0,0],#idx 15
    ...        [0,0,0,4,2,0,0,0,4],
    ...        [0,0,0,2,4,0,0,0,1],
    ...        [0,0,0,1,7,0,0,0,0]], 'float')

Then compute a distance matrix of your choosing, and perform nmds on that matrix

.. doctest::

    >>> distmtx = dist_euclidean(abundance)
    >>> nm = NMDS(distmtx, verbosity=0)

The NMDS object provides a list of points, which can be plotted if desired

.. doctest::

    >>> pts = nm.getPoints()
    >>> stress = nm.getStress()

With matplotlib installed, we could then do ``plt.plot(pts[:,0], pts[:,1])``

Hierarchical clustering (UPGMA, NJ)
===================================

Hierarchical clustering techniques work on a matrix of pairwise distances. In this case, we generate the distance matrix between a set of sequences.

**NOTE:** UPGMA should not be used for phylogenetic reconstruction.

.. doctest::

    >>> from cogent import LoadSeqs
    >>> from cogent.phylo import distance
    >>> from cogent.cluster.UPGMA import upgma

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

Now use this matrix to build a UPGMA cluster.

.. doctest::

    >>> mycluster = upgma(d.getPairwiseDistances())
    >>> print mycluster.asciiArt()
                                  /-NineBande
                        /edge.1--|
                       |         |          /-HowlerMon
              /edge.0--|          \edge.2--|
             |         |                    \-Human
    -root----|         |
             |          \-DogFaced
             |
              \-Mouse

We demonstrate saving this UPGMA cluster to a file.

.. doctest::

    >>> mycluster.writeToFile('test_upgma.tree')

..
    We don't actually want to keep that file now, so I'm importing the ``os`` module to delete it.

.. doctest::
    :hide:
    
    >>> import os
    >>> os.remove('test_upgma.tree')

