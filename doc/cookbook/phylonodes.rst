Phylonodes Methods
------------------

.. authors, Dan Knights

Basics
^^^^^^

Loading a tree from a file and visualizing it with ``ascii_art()``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> print(tr.ascii_art())
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced

Writing a tree to a file
""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> tr.write('data/temp.tree')

Getting the individual nodes of a tree by name
""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> names = tr.get_node_names()
    >>> names[:4]
    ['root', 'edge.1', 'edge.0', 'Human']
    >>> names[4:]
    ['HowlerMon', 'Mouse', 'NineBande', 'DogFaced']
    >>> names_nodes = tr.get_nodes_dict()
    >>> names_nodes['Human']
    Tree("Human;")
    >>> tr.get_node_matching_name('Mouse')
    Tree("Mouse;")

Getting the name of a node (or a tree)
""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> hu = tr.get_node_matching_name('Human')
    >>> tr.name
    'root'
    >>> hu.name
    'Human'

The object type of a tree and its nodes is the same
"""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> nodes = tr.get_nodes_dict()
    >>> hu = nodes['Human']
    >>> type(hu)
    <class 'cogent3.core.tree.PhyloNode'>
    >>> type(tr)
    <class 'cogent3.core.tree.PhyloNode'>

Working with the nodes of a tree
""""""""""""""""""""""""""""""""

Get all the nodes, tips and edges

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> nodes = tr.get_nodes_dict()
    >>> for n in nodes.items():  # doctest: +SKIP
    ...     print(n)
    ...
    ('NineBande', Tree("NineBande;"))
    ('edge.1', Tree("((Human,HowlerMon),Mouse);"))
    ('root', Tree("(((Human,HowlerMon),Mouse),NineBande,DogFaced);"))
    ('DogFaced', Tree("DogFaced;"))
    ('Human', Tree("Human;"))
    ('edge.0', Tree("(Human,HowlerMon);"))
    ('Mouse', Tree("Mouse;"))
    ('HowlerMon', Tree("HowlerMon;"))

only the terminal nodes (tips)

.. doctest::

    >>> for n in tr.iter_tips():
    ...     print(n)
    ...
    Human:0.0311054096183;
    HowlerMon:0.0415847131449;
    Mouse:0.277353608988;
    NineBande:0.0939768158209;
    DogFaced:0.113211053859;

for internal nodes (edges) we can use Newick format to simplify the output

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> for n in tr.iter_nontips():
    ...     print(n.get_newick())
    ...
    ((Human,HowlerMon),Mouse);
    (Human,HowlerMon);

Getting the path between two tips or edges (connecting edges)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> edges = tr.get_connecting_edges('edge.1','Human')
    >>> for edge in edges:
    ...    print(edge.name)
    ...
    edge.1
    edge.0
    Human

Getting the distance between two nodes
""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> nodes = tr.get_nodes_dict()
    >>> hu = nodes['Human']
    >>> mu = nodes['Mouse']
    >>> hu.distance(mu)
    0.3467553...
    >>> hu.is_tip()
    True

Getting the last common ancestor (LCA) for two nodes
""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> nodes = tr.get_nodes_dict()
    >>> hu = nodes['Human']
    >>> mu = nodes['Mouse']
    >>> lca = hu.last_common_ancestor(mu)
    >>> lca
    Tree("((Human,HowlerMon),Mouse);")
    >>> type(lca)
    <class 'cogent3.core.tree.PhyloNode'>

Getting all the ancestors for a node
""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> hu = tr.get_node_matching_name('Human')
    >>> for a in hu.ancestors():
    ...     print(a.name)
    ...
    edge.0
    edge.1
    root

Getting all the children for a node
"""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> node = tr.get_node_matching_name('edge.1')
    >>> children = list(node.iter_tips()) + list(node.iter_nontips())
    >>> for child in children:
    ...     print(child.name)
    ...
    Human
    HowlerMon
    Mouse
    edge.0

Getting all the distances for a tree
""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> dists = tr.get_distances()

We also show how to select a subset of distances involving just one species.

.. doctest::

    >>> human_dists = [names for names in dists if 'Human' in names]
    >>> for dist in human_dists:
    ...     print(dist, dists[dist])  # doctest: +SKIP
    ...
    ('Human', 'NineBande') 0.183106418165
    ('DogFaced', 'Human') 0.202340656203
    ('NineBande', 'Human') 0.183106418165
    ('Human', 'DogFaced') 0.202340656203
    ('Mouse', 'Human') 0.346755361094
    ('HowlerMon', 'Human') 0.0726901227632
    ('Human', 'Mouse') 0.346755361094
    ('Human', 'HowlerMon') 0.0726901227632


Getting the two nodes that are farthest apart
"""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> tr.max_tip_tip_distance()
    (0.4102925130849, ('Mouse', 'DogFaced'))


Get the nodes within a given distance
"""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> hu = tr.get_node_matching_name('Human')
    >>> tips = hu.tips_within_distance(0.2)
    >>> for t in tips:
    ...     print(t)
    ...
    HowlerMon:0.0415847131449;
    NineBande:0.0939768158209;

Rerooting trees
^^^^^^^^^^^^^^^

At a named node
"""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> print(tr.rooted_at('edge.0').ascii_art())
              /-Human
             |
    -root----|--HowlerMon
             |
             |          /-Mouse
              \edge.0--|
                       |          /-NineBande
                        \edge.1--|
                                  \-DogFaced


At the midpoint
"""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> print(tr.root_at_midpoint().ascii_art())
              /-Mouse
             |
    -root----|                    /-Human
             |          /edge.0--|
             |         |          \-HowlerMon
              \edge.0.2|
                       |          /-NineBande
                        \edge.1--|
                                  \-DogFaced
    >>> print(tr.ascii_art())
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-------- /-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced

Tree representations
^^^^^^^^^^^^^^^^^^^^

Newick format
"""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> tr.get_newick()
    '(((Human,HowlerMon),Mouse),NineBande,DogFaced);'
    >>> tr.get_newick(with_distances=True)
    '(((Human:0.0311054096183,HowlerMon:0.0415847131449)...

XML format
""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> xml = tr.get_xml()
    >>> for line in xml.splitlines():
    ...    print(line)
    ...
    <?xml version="1.0"?>
    <clade>
      <clade>
         <param><name>length</name><value>0.0197278502379</value></param>
        <clade>
           <param><name>length</name><value>0.0382963424874</value></param>
          <clade>
             <name>Human</name>...

.. todo write docs for display of dendrograms

Tree traversal
^^^^^^^^^^^^^^

Here is the example tree for reference:

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> print(tr.ascii_art())
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced

Preorder
""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> for t in tr.preorder():
    ...     print(t.get_newick())
    ...
    (((Human,HowlerMon),Mouse),NineBande,DogFaced);
    ((Human,HowlerMon),Mouse);
    (Human,HowlerMon);
    Human;
    HowlerMon;
    Mouse;
    NineBande;
    DogFaced;

Postorder
"""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> for t in tr.postorder():
    ...     print(t.get_newick())
    ...
    Human;
    HowlerMon;
    (Human,HowlerMon);
    Mouse;
    ((Human,HowlerMon),Mouse);
    NineBande;
    DogFaced;
    (((Human,HowlerMon),Mouse),NineBande,DogFaced);

Selecting subtrees
^^^^^^^^^^^^^^^^^^

One way to do it
""""""""""""""""

.. doctest::

    >>> from cogent3 import load_tree
    >>> tr = load_tree('data/test.tree')
    >>> for tip in tr.iter_nontips():
    ...     tip_names = tip.get_tip_names()
    ...     print(tip_names)
    ...     sub_tree = tr.get_sub_tree(tip_names)
    ...     print(sub_tree.ascii_art())
    ...     print()
    ...
    ['Human', 'HowlerMon', 'Mouse']
              /-Human
             |
    -root----|--HowlerMon
             |
              \-Mouse
    <BLANKLINE>
    ['Human', 'HowlerMon']
              /-Human
    -root----|
              \-HowlerMon
    <BLANKLINE>

..
    We do some file clean up

.. doctest::
    :hide:

    >>> from cogent3.util.misc import remove_files
    >>> remove_files(['data/temp.tree', 'data/temp.pdf'],
    ...                 error_on_missing=False)
