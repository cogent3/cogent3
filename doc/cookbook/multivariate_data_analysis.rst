.. _multivariate-analysis:

**************************
Multivariate data analysis
**************************

.. sectionauthor Justin Kuczynski, Catherine Lozupone, Andreas Wilm

Principal Coordinates Analysis
==============================

Principal Coordinates Analysis works on a matrix of pairwise distances. In this example we start by calculating the pairwise distances for a set of aligned sequences, though any distance matrix can be used with PCoA, relating any objects, not only sequences.

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
    >>> print PCoA_result # doctest: +SKIP
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


Fast-MDS
========

The eigendecomposition step in Principal Coordinates Analysis (PCoA)
doesn't scale very well. And for thousands of objects the computation
of all pairwise distance alone can get very slow, because it scales
quadratically. For a huge number of objects this might even pose a
memory problem. Fast-MDS methods approximate an MDS/PCoA solution and
do not suffer from these problems.

First, let's simulate a big data sample by creating 4000 objects living
in 10 dimension. Then compute their pairwise distances and perform a
principal coordinates analysis on it. Note that the last two steps might take
already a couple of minutes.

.. doctest::

    >>> from cogent.maths.distance_transform import dist_euclidean
    >>> from cogent.cluster.metric_scaling import principal_coordinates_analysis
    >>> from numpy import random
    >>> objs = random.random((4000, 10))
    >>> distmtx = dist_euclidean(objs);
    >>> full_pcoa = principal_coordinates_analysis(distmtx)


PyCogent implements two fast MDS approximations called
Split-and-Combine MDS (SCMDS, still in development) and Nystrom (also known as
Landmark-MDS). Both can easily handle many thousands objects. One
reason is that they don't require all distances to be computed.
Instead you pass down the distance function and only required
distances are calculated.

Nystrom works by using a so called seed-matrix, which contains (only) k by
n distances, where n is the total number of objects and k<<n. The
bigger k, the more exact the approximation will be and the longer the
computation will take. One further difference to normal Principal
Coordinates Analysis is, that no eigenvalues, but only approximate
eigenvectors of length dim will be returned.

.. doctest::

   >>> from cogent.cluster.approximate_mds import nystrom
   >>> n_seeds = 100
   >>> dim = 4
   >>> nystrom_frontend(distmtx[:n_seeds], dim)

A good rule of thumb for picking n_seeds is log(n), log(n)**2 or
sqrt(n).


SCMDS works by dividing the pairwise distance matrix into chunks of
certain size and overlap. MDS is performed on each chunk individually
and the resulting solutions are progressively joined. As in the case
of Nystrom not all distances will be computed, but only those of the
overlapping tiles. The size and overlap of the tiles determine the
quality of the approximation as well as the run-time.

.. doctest::

   >>> from cogent.cluster.approximate_mds import scmds
   >>> n_full = 4000
   >>> tile_size = 500
   >>> tile_overlap = 50
   >>> dim = 4
   >>> dist_func = lambda x, y: distmtx[x, y]
   >>> scmds(n_full, tile_size, tile_overlap, dim, dist_func)


If you want to know how good the returned approximations are, you will
have to perform principal_coordinates_analysis() on a smallish
submatrix and perform a goodness_of_fit analysis.







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

Then compute a distance matrix using euclidean distance, and perform nmds on that matrix

.. doctest::

    >>> euc_distmtx = dist_euclidean(abundance)
    >>> nm = NMDS(euc_distmtx, verbosity=0)

The NMDS object provides a list of points, which can be plotted if desired

.. doctest::

    >>> pts = nm.getPoints()
    >>> stress = nm.getStress()

With matplotlib installed, we could then do ``plt.plot(pts[:,0], pts[:,1])``

Hierarchical clustering (UPGMA, NJ)
===================================

Hierarchical clustering techniques work on a matrix of pairwise distances. In this case, we use the distance matrix from the NMDS example, relating samples of species to one another using UPGMA (NJ below).

.. note:: UPGMA should not be used for phylogenetic reconstruction.

.. doctest::

    >>> from cogent.cluster.UPGMA import upgma

we start with the distance matrix and list of sample names:

.. doctest::

    >>> sample_names = ['sample'+str(i) for i in range(len(euc_distmtx))]

make 2d dict:

.. doctest::

    >>> euc_distdict = {}
    >>> for i in range(len(sample_names)):
    ...    for j in range(len(sample_names)):
    ...        euc_distdict[(sample_names[i],sample_names[j])]=euc_distmtx[i,j]

e.g.: ``euc_distdict[('sample6', 'sample5')] == 3.7416573867739413``

Now use this matrix to build a UPGMA cluster.

.. doctest::

    >>> mycluster = upgma(euc_distdict)
    >>> print mycluster.asciiArt()
                                                      /-sample10
                                            /edge.3--|
                                  /edge.2--|          \-sample8
                                 |         |
                                 |          \-sample9
                        /edge.1--|
                       |         |                    /-sample12
                       |         |          /edge.5--|
                       |         |         |          \-sample11
                       |          \edge.4--|
                       |                   |          /-sample6
                       |                    \edge.6--|
              /edge.0--|                              \-sample7
             |         |
             |         |                                        /-sample15
             |         |                              /edge.10-|
             |         |                    /edge.9--|          \-sample14
             |         |                   |         |
             |         |          /edge.8--|          \-sample13
             |         |         |         |
             |          \edge.7--|          \-sample16
    -root----|                   |
             |                   |          /-sample17
             |                    \edge.11-|
             |                              \-sample18
             |
             |                              /-sample5
             |                    /edge.14-|
             |          /edge.13-|          \-sample4
             |         |         |
             |         |          \-sample3
              \edge.12-|
                       |                    /-sample2
                       |          /edge.16-|
                        \edge.15-|          \-sample1
                                 |
                                  \-sample0

We demonstrate saving this UPGMA cluster to a file.

.. doctest::

    >>> mycluster.writeToFile('test_upgma.tree')

..
    We don't actually want to keep that file now, so I'm importing the ``os`` module to delete it.

.. doctest::
    :hide:

    >>> import os
    >>> os.remove('test_upgma.tree')

We can use neighbor joining (NJ) instead of UPGMA:

.. doctest::

    >>> from cogent.phylo.nj import nj
    >>> njtree = nj(euc_distdict)
    >>> print njtree.asciiArt()
              /-sample16
             |
             |                    /-sample12
             |          /edge.2--|
             |         |         |          /-sample13
             |         |          \edge.1--|
             |         |                   |          /-sample14
             |         |                    \edge.0--|
             |         |                              \-sample15
             |         |
             |         |                              /-sample7
             |-edge.14-|                    /edge.5--|
             |         |                   |         |          /-sample8
             |         |                   |          \edge.4--|
             |         |          /edge.6--|                   |          /-sample10
             |         |         |         |                    \edge.3--|
             |         |         |         |                              \-sample9
    -root----|         |         |         |
             |         |         |          \-sample11
             |         |         |
             |          \edge.13-|                    /-sample6
             |                   |                   |
             |                   |                   |                              /-sample4
             |                   |          /edge.10-|                    /edge.7--|
             |                   |         |         |          /edge.8--|          \-sample3
             |                   |         |         |         |         |
             |                   |         |          \edge.9--|          \-sample5
             |                    \edge.12-|                   |
             |                             |                    \-sample2
             |                             |
             |                             |          /-sample0
             |                              \edge.11-|
             |                                        \-sample1
             |
             |          /-sample18
              \edge.15-|
                        \-sample17
