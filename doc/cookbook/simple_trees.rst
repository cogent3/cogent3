Trees
-----

.. authors, Gavin Huttley, Tom Elliott

Basics
^^^^^^

Loading a tree from a file and visualizing it with ``asciiArt()``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> print tr.asciiArt()
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

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> tr.writeToFile('data/temp.tree')

Getting the individual nodes of a tree by name
""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> names = tr.getNodeNames()
    >>> names[:4]
    ['root', 'edge.1', 'edge.0', 'Human']
    >>> names[4:]
    ['HowlerMon', 'Mouse', 'NineBande', 'DogFaced']
    >>> names_nodes = tr.getNodesDict()
    >>> names_nodes['Human']
    Tree("Human;")
    >>> tr.getNodeMatchingName('Mouse')
    Tree("Mouse;")

Getting the name of a node (or a tree)
""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> hu = tr.getNodeMatchingName('Human')
    >>> tr.Name
    'root'
    >>> hu.Name
    'Human'

The object type of a tree and its nodes is the same
"""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> nodes = tr.getNodesDict()
    >>> hu = nodes['Human']
    >>> type(hu)
    <class 'cogent.core.tree.PhyloNode'>
    >>> type(tr)
    <class 'cogent.core.tree.PhyloNode'>

Working with the nodes of a tree
""""""""""""""""""""""""""""""""

Get all the nodes, tips and edges

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> nodes = tr.getNodesDict()
    >>> for n in nodes.items():
    ...     print n
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

    >>> for n in tr.iterTips():
    ...     print n
    ...
    Human:0.0311054096183;
    HowlerMon:0.0415847131449;
    Mouse:0.277353608988;
    NineBande:0.0939768158209;
    DogFaced:0.113211053859;

for internal nodes (edges) we can use Newick format to simplify the output

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> for n in tr.iterNontips():
    ...     print n.getNewick()
    ...
    ((Human,HowlerMon),Mouse);
    (Human,HowlerMon);

Getting the path between two tips or edges (connecting edges)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> edges = tr.getConnectingEdges('edge.1','Human')
    >>> for edge in edges:
    ...    print edge.Name
    ...
    edge.1
    edge.0
    Human

Getting the distance between two nodes
""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> nodes = tr.getNodesDict()
    >>> hu = nodes['Human']
    >>> mu = nodes['Mouse']
    >>> hu.distance(mu)
    0.34675536109369998
    >>> hu.isTip()
    True

Getting the last common ancestor (LCA) for two nodes
""""""""""""""""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> nodes = tr.getNodesDict()
    >>> hu = nodes['Human']
    >>> mu = nodes['Mouse']
    >>> lca = hu.lastCommonAncestor(mu)
    >>> lca
    Tree("((Human,HowlerMon),Mouse);")
    >>> type(lca)
    <class 'cogent.core.tree.PhyloNode'>

Getting all the ancestors for a node
""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> hu = tr.getNodeMatchingName('Human')
    >>> for a in hu.ancestors():
    ...     print a.Name
    ...
    edge.0
    edge.1
    root

Getting all the children for a node
"""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> node = tr.getNodeMatchingName('edge.1')
    >>> children = list(node.iterTips()) + list(node.iterNontips())
    >>> for child in children:
    ...     print child.Name
    ...
    Human
    HowlerMon
    Mouse
    edge.0

Getting all the distances for a tree
""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> dists = tr.getDistances()

We also show how to select a subset of distances involving just one species.

.. doctest::

    >>> human_dists = [names for names in dists if 'Human' in names]
    >>> for dist in human_dists:
    ...     print dist, dists[dist]
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

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> tr.maxTipTipDistance()
    (0.4102925130849, ('Mouse', 'DogFaced'))


Get the nodes within a given distance
"""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> hu = tr.getNodeMatchingName('Human')
    >>> tips = hu.tipsWithinDistance(0.2)
    >>> for t in tips:
    ...     print t
    ...
    HowlerMon:0.0415847131449;
    NineBande:0.0939768158209;

Rerooting trees
^^^^^^^^^^^^^^^

At a named node
"""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> print tr.rootedAt('edge.0').asciiArt()
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

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> print tr.rootAtMidpoint().asciiArt()
              /-Mouse
             |
    -root----|                    /-Human
             |          /edge.0--|
             |         |          \-HowlerMon
              \edge.0.2|
                       |          /-NineBande
                        \edge.1--|
                                  \-DogFaced
    >>> print tr.asciiArt()
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-------- /-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced

Near a given tip
""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> print tr.asciiArt()
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-------- /-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced
    >>> print tr.rootedWithTip("Mouse").asciiArt()
              /-Mouse
             |
    -root----|                    /-Human
             |          /edge.0--|
             |         |          \-HowlerMon
              \edge.0.2|
                       |          /-NineBande
                        \edge.1--|
                                  \-DogFaced
        
                       |          /-NineBande
                        \edge.1--|
                                  \-DogFaced


Tree representations
^^^^^^^^^^^^^^^^^^^^

Newick format
"""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> tr.getNewick()
    '(((Human,HowlerMon),Mouse),NineBande,DogFaced);'
    >>> tr.getNewick(with_distances=True)
    '(((Human:0.0311054096183,HowlerMon:0.0415847131449)...

XML format
""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> xml = tr.getXML()
    >>> for line in xml.splitlines():
    ...    print line
    ...
    <?xml version="1.0"?>
    <clade>
      <clade>
         <param><name>length</name><value>0.0197278502379</value></param>
        <clade>
           <param><name>length</name><value>0.0382963424874</value></param>
          <clade>
             <name>Human</name>...

Write to PDF
""""""""""""

.. note:: This requires ``matplotlib``. It will bring up a ``matplotlib`` window if run from the command line. But in any case, it will write the pdf file to the data directory.

.. doctest::

    >>> from cogent import LoadTree
    >>> from cogent.draw import dendrogram
    >>> tr = LoadTree('data/test.tree')
    >>> h, w = 500, 500
    >>> np = dendrogram.ContemporaneousDendrogram(tr)
    >>> np.drawToPDF('temp.pdf', w, h, font_size=14)

.. doctest::
    :hide:
    
    >>> from cogent.util.misc import remove_files
    >>> remove_files('temp.pdf', error_on_missing=False)


Tree traversal
^^^^^^^^^^^^^^

Here is the example tree for reference:

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> print tr.asciiArt()
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

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> for t in tr.preorder():
    ...     print t.getNewick()
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

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> for t in tr.postorder():
    ...     print t.getNewick()
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

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> for tip in tr.iterNontips():
    ...     tip_names = tip.getTipNames()
    ...     print tip_names
    ...     sub_tree = tr.getSubTree(tip_names)
    ...     print sub_tree.asciiArt()
    ...     print
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

    >>> from cogent.util.misc import remove_files
    >>> remove_files(['data/temp.tree', 'data/temp.pdf'],
    ...                 error_on_missing=False)


Tree manipulation methods
^^^^^^^^^^^^^

Pruning the tree
""""""""""""""""

Remove internal nodes with only one child. Create new connections
and branch lengths (if tree is a PhyloNode) to reflect the change. 

.. doctest::

    >>> from cogent import LoadTree
    >>> simple_tree_string="(B:0.2,(D:0.4)E:0.5)F;"
    >>> simple_tree=LoadTree(treestring=simple_tree_string)
    >>> print simple_tree.asciiArt()
              /-B
    -F-------|
              \E------- /-D
    >>> simple_tree.prune()
    >>> print simple_tree.asciiArt()
              /-B
    -F-------|
              \-D
    >>> print simple_tree
    (B:0.2,D:0.9)F;


Create a full unrooted copy of the tree
"""""""""""""""""""""""""""""""""""""""

.. doctest::

    >>> from cogent import LoadTree
    >>> tr1 = LoadTree('data/test.tree')
    >>> print tr1.getNewick()
    (((Human,HowlerMon),Mouse),NineBande,DogFaced);
    >>> tr2 = tr1.unrootedDeepcopy()
    >>> print tr2.getNewick()
    (((Human,HowlerMon),Mouse),NineBande,DogFaced);


Transform tree into a bifurcating tree
""""""""""""""""""""""""""""""""""""""

Add internal nodes so that every node has 2 or fewer children.

.. doctest::

    >>> from cogent import LoadTree
    >>> tree_string="(B:0.2,H:0.2,(C:0.3,D:0.4,E:0.1)F:0.5)G;"
    >>> tr = LoadTree(treestring=tree_string)
    >>> print tr.asciiArt()
              /-B
             |
             |--H
    -G-------|
             |          /-C
             |         |
              \F-------|--D
                       |
                        \-E
    >>> print tr.bifurcating().asciiArt()
              /-B
    -G-------|
             |          /-H
              \root.2--|
                       |          /-C
                        \F-------|
                                 |          /-D
                                  \root----|
                                            \-E

    
Transform tree into a balanced tree
"""""""""""""""""""""""""""""""""""

Using a balanced tree can substantially improve performance of 
likelihood calculations. Note that the resulting tree has a 
different orientation with the effect that specifying clades or 
stems for model parameterization should be done using the 
"outgroup_name" argument.

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree('data/test.tree')
    >>> print tr.asciiArt()
                                  /-Human
                        /edge.0--|
              /edge.1--|          \-HowlerMon
             |         |
             |          \-Mouse
    -root----|
             |--NineBande
             |
              \-DogFaced
    >>> print tr.balanced().asciiArt()
                        /-Human
              /edge.0--|
             |          \-HowlerMon
             |
    -root----|--Mouse
             |
             |          /-NineBande
              \edge.1--|
                        \-DogFaced
    

Test two trees for same topology
""""""""""""""""""""""""""""""""

Branch lengths don't matter.

.. doctest::

    >>> from cogent import LoadTree
    >>> tr1 = LoadTree(treestring="(B:0.2,(C:0.2,D:0.2)F:0.2)G;")
    >>> tr2 = LoadTree(treestring="((C:0.1,D:0.1)F:0.1,B:0.1)G;")
    >>> tr1.sameTopology(tr2)
    True
    


Calculate each node's maximum distance to a tip
"""""""""""""""""""""""""""""""""""""""""""""""

Sets each node's "TipDistance" attribute to be
the distance from that node to its most distant tip.

.. doctest::

    >>> from cogent import LoadTree
    >>> tr = LoadTree(treestring="(B:0.2,(C:0.3,D:0.4)F:0.5)G;")
    >>> print tr.asciiArt()
              /-B
    -G-------|
             |          /-C
              \F-------|
                        \-D
    >>> tr.setTipDistances()
    >>> for t in tr.preorder():
    ...     print t.Name, t.TipDistance
    ... 
    G 0.9
    B 0
    F 0.4
    C 0
    D 0
    

scaleBranchLengths()
""""""""""""""""""""
Scales BranchLengths in place to integers for ascii output.

        Warning: tree might not be exactly the length you specify.

        Set ultrametric=True if you want all the root-tip distances to end
        up precisely the same.


tipToTipDistances()
"""""""""""""""""""
           """Returns distance matrix between all pairs of tips, and a tip order.
            
        Warning: .__start and .__stop added to self and its descendants.

        tip_order contains the actual node objects, not their names (may be
        confusing in some cases).
        """


compareByTipDistances()
"""""""""""""""""""""""
        """Compares self to other using tip-to-tip distance matrices.

        Value returned is dist_f(m1, m2) for the two matrices. Default is
        to use the Pearson correlation coefficient, with +1 giving a distance
        of 0 and -1 giving a distance of +1 (the madimum possible value).
        Depending on the application, you might instead want to use
        distance_from_r_squared, which counts correlations of both +1 and -1
        as identical (0 distance).
        
        Note: automatically strips out the names that don't match (this is
        necessary for this method because the distance between non-matching 
        names and matching names is undefined in the tree where they don't 
        match, and because we need to reorder the names in the two trees to 
        match up the distance matrices).
        """
