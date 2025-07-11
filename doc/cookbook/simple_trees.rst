.. jupyter-execute::
    :hide-code:

    import set_working_directory

Trees
-----

.. authors, Gavin Huttley, Tom Elliott

Loading a tree from a file
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr)

The path to the original file is stored in the ``.source`` attribute.

.. jupyter-execute::

    tr.source

Visualising a tree with ``ascii_art()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.ascii_art())

.. note:: See the :ref:`tree-gallery` for interactive graphical display of dendrograms.

Writing a tree to a file
^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    tr.write("data/temp.tree")

Getting a ``dict`` nodes keyed by their name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    names_nodes = tr.get_nodes_dict()
    names_nodes["Human"]

Getting the name of a node
^^^^^^^^^^^^^^^^^^^^^^^^^^

The root node name defaults to ``"root"``.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    tr.name

.. jupyter-execute::

    hu = tr.get_node_matching_name("Human")
    hu.name

You can ensure internal nodes get named

.. jupyter-execute::

    tr.name_unnamed_nodes()

The object type of a tree and its nodes is the same
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    type(tr)

.. jupyter-execute::

    nodes = tr.get_nodes_dict()
    hu = tr.get_node_matching_name("Human")
    type(hu)

Working with the nodes of a tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get all the nodes, tips and edges as a ``dict``.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    nodes = tr.get_nodes_dict()
    for n in nodes.items():
        print(n)

As a list.

.. jupyter-execute::

    nodes = tr.get_edge_vector()

Only the tip (terminal) nodes as a list.

.. jupyter-execute::

    tips = tr.tips()

Iterate the tip nodes.

.. jupyter-execute::

    for n in tr.iter_tips():
        print(n.name)

Get just the internal nodes as a list

.. jupyter-execute::

    non_tips = tr.nontips()

or iteratively.

.. jupyter-execute::

    for n in tr.iter_nontips():
        print(n.name)

Getting the path between two tips or edges (connecting nodes)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    edges = tr.get_connecting_edges("edge.1", "Human")
    for edge in edges:
        print(edge.name)

Get tip-to-root distances
^^^^^^^^^^^^^^^^^^^^^^^^^

The sum of all lengths on nodes connecting tips to the root node.

.. jupyter-execute::

    from cogent3 import make_tree

    tr = make_tree("(B:3,(C:2,D:4):5);")
    tr.tip_to_root_distances()

Can also be done for a subset of tips.

.. jupyter-execute::

    tr.tip_to_root_distances(names=["B", "D"])

Get tip-to-tip distances
^^^^^^^^^^^^^^^^^^^^^^^^

Get a distance matrix between all pairs of tips
and a list of the tip nodes.

.. jupyter-execute::

    from cogent3 import make_tree

    tr = make_tree("(B:3,(C:2,D:4)F:5)G;")
    dmat = tr.tip_to_tip_distances()
    dmat

.. note:: ``tip_to_tip_distances()`` is an alias for ``get_distances()``.

Getting the distance between two nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Via pairwise distances, which returns a ``DistanceMatrix`` instance.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    dists = tr.get_distances(names=["Human", "Mouse"])
    dists

Or directly between the node objects.

.. jupyter-execute::

    tr = load_tree("data/test.tree")
    nodes = tr.get_nodes_dict()
    hu = nodes["Human"]
    mu = nodes["Mouse"]
    hu.distance(mu)

Get sum of all branch lengths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_tree

    tr = make_tree("(B:3,(C:2,D:4)F:5)G;")
    tr.total_length()

Getting the last common ancestor (LCA) for two nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    nodes = tr.get_nodes_dict()
    hu = nodes["Human"]
    mu = nodes["Mouse"]
    lca = hu.last_common_ancestor(mu)
    lca.name, lca

Getting all the ancestors for a node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A list of all nodes to the tree root.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    hu = tr.get_node_matching_name("Human")
    for a in hu.ancestors():
        print(a.name)

Getting all the children for a node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    node = tr.get_node_matching_name("edge.1")
    children = list(node.iter_tips()) + list(node.iter_nontips())
    for child in children:
        print(child.name)

Getting all the distances for a tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On a ``TreeNode``, each branh has a weight of 1 so the distances represent the number of connected nodes. On a ``PhyloNode`` the measure is the sum of branch lengths.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    dists = tr.get_distances()
    dists

Getting the two nodes that are farthest apart
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    tr.max_tip_tip_distance()

Get the nodes within a given distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    hu = tr.get_node_matching_name("Human")
    tips = hu.tips_within_distance(0.2)
    for t in tips:
        print(t)

Rerooting trees
^^^^^^^^^^^^^^^

Reorienting a tree at a named node
""""""""""""""""""""""""""""""""""

The method name is a bit misleading. If ``tr`` is an unrooted tree (loosely, this is a tree whose root node has > 2 children) then the result is more a re-orientation of the tree rather than true root.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.rooted_at("edge.0").ascii_art())

At the midpoint
"""""""""""""""

This does produce a rooted tree.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.root_at_midpoint().ascii_art())

Root at a named edge
""""""""""""""""""""

The edge can be either a tip or an internal node.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.ascii_art())

.. jupyter-execute::

    print(tr.rooted("Mouse").ascii_art())

Tree representations
^^^^^^^^^^^^^^^^^^^^

Newick format
"""""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    tr.get_newick()

.. jupyter-execute::

    tr.get_newick(with_distances=True)

.. jupyter-execute::

    tr.get_newick(with_distances=True, with_node_names=True)

Tree traversal
^^^^^^^^^^^^^^

Here is the example tree for reference:

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.ascii_art())

Preorder
""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    for t in tr.preorder():
        print(t.get_newick())

Postorder
"""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    for t in tr.postorder():
        print(t.get_newick())

Selecting subtrees
^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import make_tree

    tr = make_tree("((a,b),((c,d),(e,f),(g,h)));")
    print(tr.ascii_art(show_internal=False))

Provide the names of nodes you want the subtree for. The  default behaviour is to force the subtree to have the same number of children at the root as the original tree, in this case 2.

.. jupyter-execute::

    subtree = tr.get_sub_tree(["c", "e", "g"])
    print(subtree.ascii_art(show_internal=False))

Use the ``as_rooted`` argument to ensure the selected subtree topology is as it existed on the original tree.

.. jupyter-execute::

    subtree = tr.get_sub_tree(["c", "e", "g"], as_rooted=True)
    print(subtree.ascii_art(show_internal=False))

..
    We do some file clean up

.. jupyter-execute::
    :hide-code:

    from cogent3.util.io import remove_files

    remove_files(["data/temp.tree", "data/temp.pdf"], error_on_missing=False)

Tree manipulation methods
^^^^^^^^^^^^^^^^^^^^^^^^^

Pruning the tree
""""""""""""""""

Remove internal nodes with only one child. Create new connections
and branch lengths (if tree is a PhyloNode) to reflect the change.

.. jupyter-execute::

    from cogent3 import make_tree

    simple_tree = make_tree("(B:0.2,(D:0.4)E:0.5);")
    print(simple_tree.ascii_art())

The ``prune()`` modifies the tree in place.

.. jupyter-execute::

    simple_tree.prune()
    print(simple_tree.ascii_art())

Create a full unrooted copy of the tree
"""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr1 = load_tree("data/test.tree")
    print(tr1.get_newick())

.. jupyter-execute::

    tr2 = tr1.unrooted_deepcopy()
    print(tr2.get_newick())

Transform tree into a bifurcating tree
""""""""""""""""""""""""""""""""""""""

Add internal nodes so that every node has 2 or fewer children.

.. jupyter-execute::

    from cogent3 import make_tree

    tree_string = "(B:0.2,H:0.2,(C:0.3,D:0.4,E:0.1)F:0.5)G;"
    tr = make_tree(tree_string)
    print(tr.ascii_art())

.. jupyter-execute::

    print(tr.bifurcating().ascii_art())

Transform tree into a balanced tree
"""""""""""""""""""""""""""""""""""

Using a balanced tree can substantially improve performance of
likelihood calculations for time-reversible models. Note that
the resulting tree has a different orientation with the effect
that specifying clades or stems for model parameterisation
should be done using the "outgroup_name" argument.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.ascii_art())

.. jupyter-execute::

    print(tr.balanced().ascii_art())

Test two trees for same topology
""""""""""""""""""""""""""""""""

Branch lengths don't matter.

.. jupyter-execute::

    from cogent3 import make_tree

    tr1 = make_tree("(B:0.2,(C:0.2,D:0.2)F:0.2)G;")
    tr2 = make_tree("((C:0.1,D:0.1)F:0.1,B:0.1)G;")
    tr1.same_topology(tr2)

Measure topological distances between two trees
"""""""""""""""""""""""""""""""""""""""""""""""

A number of topological tree distance metrics are available. They include:

* The Robinson-Foulds Distance for rooted trees.
* The Matching Cluster Distance for rooted trees.
* The Robinson-Foulds Distance for unrooted trees.
* The Lin-Rajan-Moret Distance for unrooted trees.

There are several variations of the Robinson-Foulds metric in the literature. The definition used by ``cogent3`` is the
cardinality of the symmetric difference of the sets of clades/splits in the two rooted/unrooted trees. Other definitions sometimes
divide this by two, or normalise it to the unit interval. 

The Robinson-Foulds distance is quick to compute, but is known to saturate quickly. Moving a single leaf in a tree can maximise this metric.

The Matching Cluster and Lin-Rajan-Moret are two matching-based distances that are more statistically robust. 
Unlike the Robinson-Foulds distance which counts how many of the splits/clades are not exactly same, the matching-based distances
measures the degree by which the splits/clades are different. The matching-based distances solve a min-weight matching problem,
which for large trees may take longer to compute.

.. jupyter-execute::

    # Distance metrics for rooted trees
    from cogent3 import make_tree

    tr1 = make_tree(treestring="(a,(b,(c,(d,e))));")
    tr2 = make_tree(treestring="(e,(d,(c,(b,a))));")
    
    mc_distance = tr1.tree_distance(tr2, method="matching_cluster") # or method="mc" or method="matching"
    rooted_rf_distance = tr1.tree_distance(tr2, method="rooted_robinson_foulds") # or method="rrf" or method="rf"

    print("Matching Cluster Distance:", mc_distance)
    print("Rooted Robinson Foulds Distance:", rooted_rf_distance)

.. jupyter-execute::

    # Distance metrics for unrooted trees
    from cogent3 import make_tree
    
    tr1 = make_tree(treestring="(a,b,(c,(d,e)));")
    tr2 = make_tree(treestring="((a,c),(b,d),e);")
    
    lrm_distance = tr1.tree_distance(tr2, method="lin_rajan_moret") # or method="lrm" or method="matching"
    unrooted_rf_distance = tr1.tree_distance(tr2, method="unrooted_robinson_foulds") # or method="urf" or method="rf"
    
    print("Lin-Rajan-Moret Distance:", lrm_distance)
    print("Unrooted Robinson Foulds Distance:", unrooted_rf_distance)
