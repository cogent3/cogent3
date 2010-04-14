******************
Community analysis
******************

alpha diversity
===============

Phylogenetic Diversity (PD)
---------------------------

For each environment (i.e. sample), calculates the amount of branch length in a phylogenetic tree that lead to its sequences.  

First we will load in a Newick formatted tree.

.. doctest::

    >>> from cogent.parse.tree import DndParser
    >>> from cogent.maths.unifrac.fast_tree import UniFracTreeNode
    >>> tree_in = open("data/Crump_example_tree_newick.txt")
    >>> tree = DndParser(tree_in, UniFracTreeNode)
    
Next we will load information on which sequences in the tree come from which environment.

.. doctest::

    >>> from cogent.maths.unifrac.fast_unifrac import count_envs
    >>> envs_in = open("data/Crump_et_al_example_env_file.txt")
    >>> envs = count_envs(envs_in)

Finally, we can calculate the PD values for each environment in the tree

.. doctest:: 

    >>> from cogent.maths.unifrac.fast_unifrac import PD_whole_tree
    >>> envs, PD_values = PD_whole_tree(tree, envs)
    >>> print envs
    ['E_FL', 'E_PA', 'O_FL', 'O_UN', 'R_FL', 'R_PA']
    >>> print PD_vals#doctest: +SKIP
    [ 5.85389  7.60352  2.19215  2.81821  3.93728  3.7534 ]

PD_vals is a numpy array with the values representing each environment in envs.

Rarefaction
-------------

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

The Fast UniFrac implementation in PyCogent is the source code for the Fast UniFrac web tool available at http://bmf2.colorado.edu/fastunifrac and the QIIME pipeline for Microbial community analysis http://qiime.sourceforge.net

Calculate a UniFrac Distance Matrix and apply PCoA and UPGMA
------------------------------------------------------------

The UniFrac analysis is run on open tree and environment file objects. The resulting dictionary has a distance matrix of pairwise UniFrac values ('distance_matrix'), a Newick string representing the results of performing UPGMA clustering on this distance matrix ('cluster_envs') and the results of running Principal Coordinates Analysis on the distance matrix ('pcoa'). One can specify weighted UniFrac with weighted=True. Here we run an unweighted analysis.

.. doctest::

    >>> from cogent.maths.unifrac.fast_unifrac import fast_unifrac_file
    >>> tree_in = open("data/Crump_example_tree_newick.txt")
    >>> envs_in = open("data/Crump_et_al_example_env_file.txt")
    >>> result = fast_unifrac_file(tree_in, envs_in, weighted=False)
    >>> print result['cluster_envs']
    ((((('E_FL':0.339607103063,'R_FL':0.339607103063):0.0279991540511,'R_PA':0.367606257114):0.0103026524101,'E_PA':0.377908909524):0.0223322024492,'O_UN':0.400241111973):0.00976759866402,'O_FL':0.410008710637);
    >>> print result['pcoa']
    =================================================================================================
            Type              Label  vec_num-0  vec_num-1  vec_num-2  vec_num-3  vec_num-4  vec_num-5
    -------------------------------------------------------------------------------------------------
    Eigenvectors               E_FL       0.05       0.22      -0.09      -0.26      -0.29      -0.00
    Eigenvectors               E_PA      -0.36       0.24       0.21      -0.08       0.18      -0.00
    Eigenvectors               O_FL       0.32      -0.26       0.30      -0.13       0.05      -0.00
    Eigenvectors               O_UN      -0.28      -0.40      -0.24      -0.04       0.01      -0.00
    Eigenvectors               R_FL       0.29       0.18      -0.28       0.09       0.22      -0.00
    Eigenvectors               R_PA      -0.02       0.02       0.11       0.42      -0.17      -0.00
     Eigenvalues        eigenvalues       0.40       0.36       0.29       0.27       0.19      -0.00
     Eigenvalues  var explained (%)      26.34      23.84      19.06      18.02      12.74      -0.00
    -------------------------------------------------------------------------------------------------

Perform pairwise tests of whether samples are significantly different with UniFrac
----------------------------------------------------------------------------------

The analysis is run on open tree and environment file objects. In this example, we use unweighted unifrac (weighted=False), we permute the environment assignments on the tree 50 times (num_iters=50) and we perform UniFrac on all pairs of environments (test_on="Pairwise"). A list is returned with a tuple for each pairwise comparison with items: 0 - the first environment, 1 - the second environment, 2- the uncorrected p-value and 3 - the p-value after correcting for multiple comparisons with the Bonferroni correction.

.. doctest::

    >>> from cogent.maths.unifrac.fast_unifrac import fast_unifrac_permutations_file
    >>> tree_in = open("data/Crump_example_tree_newick.txt")
    >>> envs_in = open("data/Crump_et_al_example_env_file.txt")
    >>> result = fast_unifrac_permutations_file(tree_in, envs_in, weighted=False, num_iters=50, test_on="Pairwise")
    >>> print result[0]#doctest: +SKIP
    ('E_FL', 'E_PA', 0.17999999999999999, 1.0)

Perform a single UniFrac significance test on the whole tree
------------------------------------------------------------

The analysis is run on open tree and environment file objects. In this example, we use weighted unifrac (weighted=True), we permute the environment assignments on the tree 50 times (num_iters=50) and we perform a unifrac significance test on the whole tree (test_on="Tree"). The resulting list has only one item since a single test was performed. It is a 3 item tuple where the second and third values are the p-value.

.. doctest::

    >>> from cogent.maths.unifrac.fast_unifrac import fast_unifrac_permutations_file
    >>> tree_in = open("data/Crump_example_tree_newick.txt")
    >>> envs_in = open("data/Crump_et_al_example_env_file.txt")
    >>> result = fast_unifrac_permutations_file(tree_in, envs_in, weighted=True, num_iters=50, test_on="Tree")
    >>> print result#doctest: +SKIP
    [('whole tree', 0.56000000000000005, 0.56000000000000005)]

P-test
-------

Perform pairwise tests of whether samples are significantly different with the P-test (Martin, 2002)
----------------------------------------------------------------------------------------------------

The analysis is run on open tree and environment file objects. In this example, we permute the environment assignments on the tree 50 times (num_iters=50) and perform the p test for all pairs of environments (test_on="Pairwise"). A list is returned with a tuple for each pairwise comparison with items: 0 - the first environment, 1 - the second environment, 2- the uncorrected p-value and 3 - the p-value after correcting for multiple comparisons with the Bonferroni correction.

.. doctest::

    >>> from cogent.maths.unifrac.fast_unifrac import fast_p_test_file
    >>> tree_in = open("data/Crump_example_tree_newick.txt")
    >>> envs_in = open("data/Crump_et_al_example_env_file.txt")
    >>> result = fast_p_test_file(tree_in, envs_in, num_iters=50, test_on="Pairwise")
    >>> print result[0]#doctest: +SKIP
    ('E_FL', 'E_PA', 0.040000000000000001, 0.59999999999999998)


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

