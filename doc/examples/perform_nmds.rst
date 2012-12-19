Perform Nonmetric Multidimensional Scaling
==========================================

.. sectionauthor:: Justin Kuczynski

An example of how to use nmds.

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
