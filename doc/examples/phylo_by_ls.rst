Phylogenetic reconstruction by least squares
============================================

.. sectionauthor:: Gavin Huttley

We will load some pre-computed pairwise distance data. To see how that data was computed see the :ref:`calculating-pairwise-distances` example. That data is saved in a format called ``pickle`` which is native to python. As per usual, we import the basic components we need.

.. recompute the data matrix and then delete file at end

.. doctest::
    :hide:
    
    >>> from cogent import LoadSeqs
    >>> from cogent.phylo import distance
    >>> from cogent.evolve.models import HKY85
    >>> al = LoadSeqs("data/long_testseqs.fasta")
    >>> d = distance.EstimateDistances(al, submodel= HKY85())
    >>> d.run(show_progress=False)
    >>> import cPickle
    >>> f = open('dists_for_phylo.pickle', "w")
    >>> cPickle.dump(d.getPairwiseDistances(), f)
    >>> f.close()

.. doctest::

    >>> import cPickle
    >>> from cogent.phylo import distance, least_squares

Now load the distance data.

.. doctest::

    >>> f = file('dists_for_phylo.pickle', 'r')
    >>> dists = cPickle.load(f)
    >>> f.close()

If there are extremely small distances, they can cause an error in the least squares calculation. Since such estimates are between extremely closely related sequences we could simply drop all distances for one of the sequences. We won't do that here, we'll leave that as exercise.

We make the ls calculator.

.. doctest::

    >>> ls = least_squares.WLS(dists)

We will search tree space for the collection of best trees using the advanced stepwise addition algorithm (hereafter *asaa*).

Look for the single best tree
-----------------------------

In this use case we are after just 1 tree. We specify up to what taxa size all possible trees for the sample will be computed. Here we are specifying ``a=5``. This means 5 sequences will be picked randomly and all possible trees relating them will be evaluated. ``k=1`` means only the best tree will be kept at the end of each such round of evaluation. For every remaining sequence it is grafted onto every possible branch of this tree. The best ``k`` results are then taken to the next round, when another sequence is randomly selected for addition. This proceeds until all sequences have been added. The result with following arguments is a single wls score and a single ``Tree`` which can be saved etc ..

.. doctest::
    
    >>> score, tree = ls.trex(a = 5, k = 1, show_progress = True)
    3 trees of size 4 at start
    15 trees of size 5 ... done
    >>> assert score < 1e-4

We won't display this tree, because we are doing more below.

A more rigourous tree space search
----------------------------------

We change the asaa settings, so we keep more trees and then look at the distribution of the statistics for the last collection of trees. We could also change ``a`` to be larger, but in the current case we just adjust ``k``. We also set the argument ``return_all = True``, the effect of which is to return the complete set of 145 trees. These, and their support statistic, can then be inspected.

.. doctest::

    >>> trees = ls.trex(a = 5, k = 5, show_progress = True, return_all = True)
    3 trees of size 4 at start
    15 trees of size 5 ... done

Remember the sum-of-squares statistic will be smaller for 'good' trees. The order of the trees returned is from good to bad. As indicated above, ``trees`` has length 145.

.. doctest::

    >>> print len(trees)
    15

Lets inspect the resulting statistics. First, the object ``trees`` is a list of (wls, Tree) tuples. We will therefore loop over the list to generate a separate list of just the wls statistics. The following syntax is called a list comprehension - basically just a very succinct ``for`` loop.

.. doctest::

    >>> wls_stats = [tree[0] for tree in trees]

The ``wls_stats`` is a list which, if printed, looks like

.. code-block:: python
    
    [1.3308768548934439e-05, 0.0015588630350439783, ...

From this you'll see that the first 5 results are very similar to each other and would probably reasonably be considered equivalently supported topologies. I'll just print the first two of the these trees after balancing them (in order to make their representations as equal as possible).

.. doctest::

    >>> t1 = trees[0][1].balanced()
    >>> t2 = trees[1][1].balanced()
    >>> print t1.asciiArt()
                        /-Human
              /edge.0--|
             |          \-HowlerMon
             |
    -root----|--Mouse
             |
             |          /-NineBande
              \edge.1--|
                        \-DogFaced
    >>> print t2.asciiArt()
              /-DogFaced
             |
             |          /-Human
    -root----|-edge.0--|
             |          \-HowlerMon
             |
             |          /-NineBande
              \edge.1--|
                        \-Mouse

You can see the difference involves the Jackrabbit, TreeShrew, Gorilla, Rat clade.

Assessing the fit for a pre-specified tree topology
---------------------------------------------------

In some instances we may have a tree from the literature or elsewhere whose fit to the data we seek to evaluate. In this case I'm going load a tree as follows.

.. doctest::

    >>> from cogent import LoadTree
    >>> query_tree = LoadTree(
    ... treestring="((Human:.2,DogFaced:.2):.3,(NineBande:.1, Mouse:.5):.2,HowlerMon:.1)")

We now just use the ``ls`` object created above. The following evaluates the query using it's associated branch lengths, returning only the wls statistic.

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> ls.evaluateTree(query_tree)
    2.8...

We can also evaluate just the tree's topology, returning both the wls statistic and the tree with best fit branch lengths.

.. doctest::

    >>> wls, t = ls.evaluateTopology(query_tree)
    >>> assert "%.4f" % wls == '0.0084'

Using maximum likelihood for measuring tree fit
-----------------------------------------------

This is a much slower algorithm and the interface largely mirrors that for the above. The difference is you import ``maximum_likelihood`` instead of ``least_squares``, and use the ``ML`` instead of ``WLS`` classes. The ``ML`` class requires a substitution model (like a HKY85 for DNA or JTT92 for protein), and an alignment. It also optionally takes a distance matrix, such as that used here, computed for the same sequences. These distances are then used to obtain estimates of branch lengths by the WLS method for each evaluated tree topology which are then used as starting values for the likelihood optimisation.

.. clean up

.. doctest::
    :hide:
    
    >>> import os
    >>> os.remove('dists_for_phylo.pickle')
