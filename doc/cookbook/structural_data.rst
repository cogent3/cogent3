***************
Structural data
***************

.. sectionauthor:: Kristian Rother, Patrick Yannul

Protein structures
==================

Calculating euclidean distances between atoms
---------------------------------------------

.. doctest::
    
    >>> from numpy import array
    >>> from cogent.maths.geometry import distance
    >>> a1 = array([1.0, 2.0, 3.0])
    >>> a2 = array([1.0, 4.0, 9.0])
    >>> distance(a1,a2)
    6.324...