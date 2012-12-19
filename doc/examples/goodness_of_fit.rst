The ``goodness_of_fit`` module
==============================

.. sectionauthor:: Andreas Wilm

This is a short example of how to use the goodness_of_fit module.

``goodness_of_fit`` will measure the degree of correspondence (or goodness of fit) between a multidimensional scaling (i.e. a lower dimensional representation of distances between objects) and the input distance matrix.

Let's setup some example data:

.. doctest::
    
    >>> import numpy
    >>> import cogent.cluster.goodness_of_fit as goodness_of_fit
    >>> distmat = numpy.array([
    ... [ 0.      ,  0.039806,  0.056853,  0.21595 ,  0.056853,  0.0138  ,
    ...   0.203862,  0.219002,  0.056853,  0.064283],
    ... [ 0.039806,  0.      ,  0.025505,  0.203862,  0.0208  ,  0.039806,
    ...   0.194917,  0.21291 ,  0.0208  ,  0.027869],
    ... [ 0.056853,  0.025505,  0.      ,  0.197887,  0.018459,  0.056853,
    ...   0.191958,  0.203862,  0.018459,  0.025505],
    ... [ 0.21595 ,  0.203862,  0.197887,  0.      ,  0.206866,  0.206866,
    ...   0.07956 ,  0.066935,  0.203862,  0.206866],
    ... [ 0.056853,  0.0208  ,  0.018459,  0.206866,  0.      ,  0.056853,
    ...   0.203862,  0.21595 ,  0.0138  ,  0.0208  ],
    ... [ 0.0138  ,  0.039806,  0.056853,  0.206866,  0.056853,  0.      ,
    ...   0.197887,  0.209882,  0.056853,  0.064283],
    ... [ 0.203862,  0.194917,  0.191958,  0.07956 ,  0.203862,  0.197887,
    ...   0.      ,  0.030311,  0.200869,  0.206866],
    ... [ 0.219002,  0.21291 ,  0.203862,  0.066935,  0.21595 ,  0.209882,
    ...   0.030311,  0.      ,  0.21291 ,  0.219002],
    ... [ 0.056853,  0.0208  ,  0.018459,  0.203862,  0.0138  ,  0.056853,
    ...   0.200869,  0.21291 ,  0.      ,  0.011481],
    ... [ 0.064283,  0.027869,  0.025505,  0.206866,  0.0208  ,  0.064283,
    ...   0.206866,  0.219002,  0.011481,  0.      ]])
    >>> mds_coords = numpy.array([
    ... [ 0.065233,  0.035019,  0.015413],
    ... [ 0.059604,  0.00168 , -0.003254],
    ... [ 0.052371, -0.010959, -0.014047],
    ... [-0.13804 , -0.036031,  0.031628],
    ... [ 0.063703, -0.015483, -0.00751 ],
    ... [ 0.056803,  0.031762,  0.021767],
    ... [-0.135082,  0.023552, -0.021006],
    ... [-0.150323,  0.011935, -0.010013],
    ... [ 0.06072 , -0.01622 , -0.007721],
    ... [ 0.065009, -0.025254, -0.005257]])

You are now good to compute stress values as shown in the following.

``Stress.calcKruskalStress()``
------------------------------

This computes Kruskal's Stress AKA Stress Formula 1 (Kruskal, 1964).
Kruskal gives the following numbers as guideline: 


==========  ===============
Stress [%]  Goodness of fit
==========  ===============
20          Poor
10          Fair
5           Good
2.5         Excellent
0           Perfect
==========  ===============

Our example data gives a very good stress value:

.. doctest::
    
    >>> stress = goodness_of_fit.Stress(distmat, mds_coords)
    >>> print stress.calcKruskalStress()
    0.02255...

``Stress.calcSstress()``
------------------------

This computes S-Stress (Takane, 1977). S-Stress emphasises larger dissimilarities more than smaller ones. S-Stress values behave nicer then Stress-1 values:

- The value of S-Stress is always between 0 and 1
- Values less then 0.1 mean good representation. 

Using the example data again:

.. doctest::
    
    >>> stress = goodness_of_fit.Stress(distmat, mds_coords)
    >>> print stress.calcSstress()
    0.00883...
