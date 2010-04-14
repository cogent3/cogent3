******************
Community analysis
******************

alpha diversity
===============

Rarefaction
-----------

*To be written.*

Parametric methods
------------------

*To be written.*

Nonparametric methods
---------------------

*To be written.*

beta diversity
==============

Unifrac
-------

*To be written.*

Taxon-based
-----------

Computing a distance matrix between samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pycogent provides many different ways to compute pairwise distances between objects.  cogent/maths/distance_transform.py provides a set of functions to calculate dissimilarities/distances between samples, given an abundance matrix.  Here is one example:
.. doctest::

>>> from cogent.maths.distance_transform import dist_euclidean
    
    .. note:: see distance_transform.py for other metrics than euclidean

>>> from numpy import array

>>> abundance_data = array([[1, 3],
...                        [5, 2],
...                        [0.1, 22]],'float')
    
We now have 3 samples, and the abundance of each column (e.g.: species) in that sample.  The first sample has 1 individual of species 1, 3 individuals of species 2.  We now compute the relatedness between these samples, using euclidean distance between the rows:

>>> dists = dist_euclidean(abundance_data)

.. doctest::

    >>> print str(dists.round(2)) # doctest: +SKIP
    [[  0.        ,   4.12,  19.02]
    [  4.12,   0.        ,  20.59 ]
    [ 19.02,  20.59 ,   0.        ]]
    
    
this distance matrix can be visualized via multivariate reduction techniques such as `PCoA or NMDS <./multivariate_data_analysis.html>`_.

Taxonomy
========

*To be written.*

.. need to decide on methods here

