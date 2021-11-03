.. jupyter-execute::
    :hide-code:

    import set_working_directory

Trees
-----

.. authors, Gavin Huttley, Tom Elliott

Loading a tree from a file and visualizing it with ``ascii_art()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.ascii_art())

Writing a tree to a file
^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    tr.write("data/temp.tree")

Getting the individual nodes of a tree by name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    names = tr.get_node_names()
    names[:4]

.. jupyter-execute::

    names[4:]
    names_nodes = tr.get_nodes_dict()
    names_nodes["Human"]

.. jupyter-execute::

    tr.get_node_matching_name("Mouse")

Getting the name of a node (or a tree)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    hu = tr.get_node_matching_name("Human")
    tr.name

.. jupyter-execute::

    hu.name

The object type of a tree and its nodes is the same
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    nodes = tr.get_nodes_dict()
    hu = nodes["Human"]
    type(hu)

.. jupyter-execute::

    type(tr)

Working with the nodes of a tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get all the nodes, tips and edges

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    nodes = tr.get_nodes_dict()
    for n in nodes.items():
        print(n)

only the terminal nodes (tips)

.. jupyter-execute::

    for n in tr.iter_tips():
        print(n)

for internal nodes (edges) we can use Newick format to simplify the output

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    for n in tr.iter_nontips():
        print(n.get_newick())

Getting the path between two tips or edges (connecting edges)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    edges = tr.get_connecting_edges("edge.1", "Human")
    for edge in edges:
        print(edge.name)

Getting the distance between two nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    nodes = tr.get_nodes_dict()
    hu = nodes["Human"]
    mu = nodes["Mouse"]
    hu.distance(mu)
    hu.is_tip()

Getting the last common ancestor (LCA) for two nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    nodes = tr.get_nodes_dict()
    hu = nodes["Human"]
    mu = nodes["Mouse"]
    lca = hu.last_common_ancestor(mu)
    lca

.. jupyter-execute::

    type(lca)

Getting all the ancestors for a node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    dists = tr.get_distances()

We also show how to select a subset of distances involving just one species.

.. jupyter-execute::

    human_dists = [names for names in dists if "Human" in names]
    for dist in human_dists:
        print(dist, dists[dist])

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

At a named node
"""""""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.rooted_at("edge.0").ascii_art())

At the midpoint
"""""""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.root_at_midpoint().ascii_art())

.. jupyter-execute::

    print(tr.ascii_art())

Near a given tip
""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    print(tr.ascii_art())

.. jupyter-execute::

    print(tr.rooted_with_tip("Mouse").ascii_art())

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

XML format
""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    xml = tr.get_xml()
    for line in xml.splitlines():
        print(line)

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

One way to do it
""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = load_tree("data/test.tree")
    for tip in tr.iter_nontips():
        tip_names = tip.get_tip_names()
        print(tip_names)
        sub_tree = tr.get_sub_tree(tip_names)
        print(sub_tree.ascii_art())

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

    simple_tree_string = "(B:0.2,(D:0.4)E:0.5)F;"
    simple_tree = make_tree(simple_tree_string)
    print(simple_tree.ascii_art())

.. jupyter-execute::

    simple_tree.prune()
    print(simple_tree.ascii_art())

.. jupyter-execute::

    print(simple_tree)

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

    from cogent3 import load_tree

    tree_string = "(B:0.2,H:0.2,(C:0.3,D:0.4,E:0.1)F:0.5)G;"
    tr = make_tree(tree_string)
    print(tr.ascii_art())

.. jupyter-execute::

    print(tr.bifurcating().ascii_art())

Transform tree into a balanced tree
"""""""""""""""""""""""""""""""""""

Using a balanced tree can substantially improve performance of
likelihood calculations. Note that the resulting tree has a
different orientation with the effect that specifying clades or
stems for model parameterisation should be done using the
"outgroup_name" argument.

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

    from cogent3 import load_tree

    tr1 = make_tree("(B:0.2,(C:0.2,D:0.2)F:0.2)G;")
    tr2 = make_tree("((C:0.1,D:0.1)F:0.1,B:0.1)G;")
    tr1.same_topology(tr2)

Calculate each node's maximum distance to a tip
"""""""""""""""""""""""""""""""""""""""""""""""

Sets each node's "TipDistance" attribute to be
the distance from that node to its most distant tip.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = make_tree("(B:0.2,(C:0.3,D:0.4)F:0.5)G;")
    print(tr.ascii_art())

.. jupyter-execute::

    tr.set_tip_distances()
    for t in tr.preorder():
        print(t.name, t.TipDistance)

Scale branch lengths in place to integers for ascii output
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. jupyter-execute::

    from cogent3 import load_tree

    tr = make_tree("(B:0.2,(C:0.3,D:0.4)F:0.5)G;")
    print(tr)

.. jupyter-execute::

    tr.scale_branch_lengths()
    print(tr)

Get tip-to-tip distances
""""""""""""""""""""""""
Get a distance matrix between all pairs of tips
and a list of the tip nodes.

.. jupyter-execute::

    from cogent3 import load_tree

    tr = make_tree("(B:3,(C:2,D:4)F:5)G;")
    d, tips = tr.tip_to_tip_distances()
    for i, t in enumerate(tips):
        print(t.name, d[i])

Compare two trees using tip-to-tip distance matrices
""""""""""""""""""""""""""""""""""""""""""""""""""""

Score ranges from 0 (minimum distance) to 1 (maximum
distance). The default is to use Pearson's correlation,
in which case a score of 0 means that the Pearson's
correlation was perfectly good (1), and a score of 1
means that the Pearson's correlation was perfectly bad (-1).

Note: automatically strips out the names that don't match.

.. jupyter-execute::

    from cogent3 import load_tree

    tr1 = make_tree("(B:2,(C:3,D:4)F:5)G;")
    tr2 = make_tree("(C:2,(B:3,D:4)F:5)G;")
    tr1.compare_by_tip_distances(tr2)