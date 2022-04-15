Manipulation of Tree Node Objects
=================================

.. sectionauthor:: Tony Walters

Examples of how to initialize and manipulate various tree node objects.

.. jupyter-execute::

    from cogent3 import load_tree, make_tree
    from cogent3.core.tree import PhyloNode
    from cogent3.parse.tree import DndParser

The general method to initialize a tree is ``load_tree``, however, for exceptionally large trees or if one needs to specify the node objects (``TreeNode``, ``PhyloNode``, or ``RangeNode``), ``DndParser`` should be used.  ``load_tree`` uses ``PhyloNode`` objects by default.

The basic properties of the tree node objects are:

- ``TreeNode`` objects are general purpose in nature, and lack phylogenetic distance values.
- ``PhyloNode`` objects inherit the methods of the ``TreeNode`` class and in addition contain phylogenetic distances.
- ``RangeNode`` objects contain evolution simulation methods in addition to the standard features of a ``PhyloNode``.

The following demonstrates the two methods for initializing a phylogenetic tree object.

.. jupyter-execute::

    simple_tree_string = "(B:0.2,(C:0.3,D:0.4)E:0.5)F;"
    complex_tree_string = "(((363564 AB294167.1 Alkalibacterium putridalgicola:0.0028006,55874 AB083411.1 Marinilactibacillus psychrotolerans:0.0022089):0.40998,(15050 Y10772.1 Facklamia hominis:0.32304,(132509 AY707780.1 Aerococcus viridans:0.58815,((143063 AY879307.1 Abiotrophia defectiva:0.5807,83619 AB042060.1 Bacillus schlegelii:0.23569):0.03586,169722 AB275483.1 Fibrobacter succinogenes:0.38272):0.06516):0.03492):0.14265):0.63594,(3589 M62687.1 Fibrobacter intestinalis:0.65866,314063 CP001146.1 Dictyoglomus thermophilum:0.38791):0.32147,276579 EU652053.1 Thermus scotoductus:0.57336);"
    simple_tree = make_tree(simple_tree_string)
    complex_tree = DndParser(complex_tree_string, PhyloNode)

Now to displaying, creating, deleting, and inserting a node in simple_tree.  Note that simple_tree has three tips, one internal node 'E', and the root 'F.'  For this example, we will create a node named 'A', with a distance of 0.1, delete the node 'C' through its parent, the internal node 'E', and finally we will insert 'A' where 'C' once was.

Display the original tree.

.. jupyter-execute::

    print(simple_tree.ascii_art())

Create a new node object.

.. jupyter-execute::

    A_node = PhyloNode(name="A", Length=0.1)

Display the children of the root node, one of which is the parent of the tip we wish to alter.  To add or remove a node, we need to use the parent of the target node, which in this case is the internal node 'E.'

.. jupyter-execute::

    print(simple_tree.children)

Remove the 'C' tip.  **Note:** ``remove()`` and ``remove_node()`` return 'True' if a node is removed, 'False' if they cannot remove a node.

.. jupyter-execute::

    simple_tree.children[1].remove("C")

Insert the new 'A' tip where 'C' was previously.

.. jupyter-execute::

    simple_tree.children[1].insert(0, A_node)

Finally, display the modified tree.

.. jupyter-execute::

    print(simple_tree.ascii_art())

When deleting tree nodes, it is often desirable to clean up any unbranched internal nodes that may have resulted from removal of tips.  For example, if we wanted to delete the node 'A' that was previously added, the resulting tree would have an unbranched internal node 'E.'

.. jupyter-execute::

    simple_tree.children[1].remove("A")
    print(simple_tree.ascii_art())

With the ``prune()`` method, internal nodes with only a single branch are removed.

.. jupyter-execute::

    simple_tree.prune()
    print(simple_tree.ascii_art())

An Example of Conditional Tree Node Modifications
=================================================

Now to look at the more complex and realistic tree.  In complex_tree, there are no internal nodes or a defined root.  In order to display this tree in a more succinct manner, we can rename these tips to only contain the genus and species names.  To step through the tips only, we can use the ``iter_tips()`` iterator, and rename each node.  The ``ascii_art()`` function, by default, will attempt to display internal nodes; this can be suppressed by the parameter ``show_internal=False``.

First, let's split the ungainly name string for each tip and only preserve the genus and species component, separated by a space.

.. jupyter-execute::

    for n in complex_tree.iter_tips():
        n.name = n.name.split()[2] + " " + n.name.split()[3]

Now we display the tree with ``ascii_art()``.

.. jupyter-execute::

    print(complex_tree.ascii_art(show_internal=False))

For another example of manipulating a phylogenetic tree, let us suppose that we want to remove any species in the tree that are not closely related to *Aerococcus viridans*.  To do this, we will delete any nodes that have a greater phylogenetic distance than 1.8 from *Aerococcus viridans*.  The best method to remove a large number of nodes from a tree is to first create a list of nodes to delete, followed by the actual removal process.  It is important that the ``prune()`` function be called after deletion of each node to ensure that internal nodes whose tips are deleted are removed instead of becoming tips.  Alternatively, one could test for internal nodes whose children are deleted in the procedure and flag these nodes to be deleted as well.

First, generate a list of tip nodes.

.. jupyter-execute::

    tips = complex_tree.tips()

Next, iterate through this list, compare the distances to *Aerococcus*, and append to the deletion list if greater than 1.8.

.. jupyter-execute::

    tips_to_delete = []
    AEROCOCCUS_INDEX = 3
    for n in tips:
        if tips[AEROCOCCUS_INDEX].distance(n) > 1.8:
            tips_to_delete.append(n)

Now for the actual deletion process.  We can simply use the parent of each node in the deletion list to remove itself.  Pruning is necessary to prevent internal nodes from being left as tips.  **Note:** ``remove()`` and ``remove_node()`` return 'True' if a node is successfully removed, 'False' otherwise.

.. jupyter-execute::

    for n in tips_to_delete:
        n.parent.remove(n)
        complex_tree.prune()

Finally, print the modified ``complex_tree``.

.. jupyter-execute::

    print(complex_tree.ascii_art(show_internal=False))
