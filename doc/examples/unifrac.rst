Run a Fast Unifrac community analysis
=====================================

.. sectionauthor:: Justin Kuczynski

Below is a simple example of using the fast unifrac function.

first we import some tools

.. doctest::

    >>> from cogent.parse.tree import DndParser
    >>> from cogent.maths.unifrac.fast_unifrac import fast_unifrac
    >>> from cogent.maths.unifrac.fast_tree import UniFracTreeNode

then we make a small example tree with tips B, C, D representing the relationship
between species B, C, and D

.. doctest::

    >>> tree_str = "(B:0.2,(C:0.3,D:0.4)E:0.6)F;"
    >>> tr = DndParser(tree_str, UniFracTreeNode)
    >>> print tr.asciiArt() # doctest: +SKIP
              /-B
    -F-------|
             |          /-C
              \E-------|
                        \-D

here's what the sample (rows) by sequence (cols) abundance matrix looks like::

    ...    [10,11,0]
    ...    [2,0,9]
    ...    [2,2,2]

and here it is in dict format for unifrac

.. doctest::

    >>> envs = {'B':{'sample1':10, 'sample2':2, 'sample3':2},
    ...        'C':{'sample1':11,'sample2':0, 'sample3':2},
    ...        'D':{'sample1':0, 'sample2':9, 'sample3':2}
    ...        }
    

now we run unifrac::


    
    >>> res = fast_unifrac(tr, envs)
    >>> print res['distance_matrix'] # doctest: +SKIP
    
    (array([[ 0.        ,  0.46666667,  0.26666667],
           [ 0.46666667,  0.        ,  0.2       ],
           [ 0.26666667,  0.2       ,  0.        ]]),
           ['sample1', 'sample2', 'sample3'])
    

the pcoa results are misleading for such a small dataset, but the distance
matrix is accurate
