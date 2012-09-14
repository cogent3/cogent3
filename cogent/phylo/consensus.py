#! /usr/bin/env python
"""This module implements methods for generating consensus trees from a list of trees"""

from cogent.core.tree import TreeBuilder
from cogent import LoadTree

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Matthew Wakefield", "Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"


def majorityRule(trees, strict=False):
    """Determines the consensus tree from a list of rooted trees using the 
     majority rules method of Margush and McMorris 1981
    Arguments:
        - trees: A list of cogent.evolve.tree objects
        - strict: A boolean flag for strict majority rule tree
          construction when true only nodes occurring >50% will be used
          when false the highest scoring < 50% node will be used if
          there is more than one node with the same score this will be
          arbitrarily chosen on sort order
    
    Returns:
        a list of cogent.evolve.tree objects
    """
    trees = [(1, tree) for tree in trees]
    return weightedMajorityRule(trees, strict, "count")

def weightedMajorityRule(weighted_trees, strict=False, attr="support"):
    cladecounts = {}
    edgelengths = {}
    total = 0
    for (weight, tree) in weighted_trees:
        total += weight
        edges = tree.getEdgeVector()
        for edge in edges:
            tips = edge.getTipNames(includeself=True)
            tips = frozenset(tips)
            if tips not in cladecounts:
                cladecounts[tips] = 0
            cladecounts[tips] += weight
            length = edge.Length and edge.Length * weight
            if edgelengths.get(tips, None):
                edgelengths[tips] += length
            else:
                edgelengths[tips] = length
    cladecounts = [(count, clade) for (clade, count) in cladecounts.items()]
    cladecounts.sort()
    cladecounts.reverse()
    
    if strict:
        # Remove any with support < 50%
        for index, (count, clade) in enumerate(cladecounts):
            if count <= 0.5 * total:
                cladecounts = cladecounts[:index]
                break
    
    # Remove conflicts
    accepted_clades = set()
    counts = {}
    for (count, clade) in cladecounts:
        for accepted_clade in accepted_clades:
            if clade.intersection(accepted_clade) and not (
                    clade.issubset(accepted_clade) or
                    clade.issuperset(accepted_clade)):
                        break
        else:
            accepted_clades.add(clade)
            counts[clade] = count
            weighted_length = edgelengths[clade]
            edgelengths[clade] = weighted_length and weighted_length / total
    
    nodes = {}
    queue = []
    tree_build = TreeBuilder().createEdge    
    for clade in accepted_clades:
        if len(clade) == 1:
            tip_name = iter(clade).next()
            params = {'length':edgelengths[clade], attr:counts[clade]}
            nodes[tip_name] = tree_build([], tip_name, params)
        else:
            queue.append(((len(clade), clade)))
            
    while queue:
        queue.sort()
        (size, clade) = queue.pop(0)
        new_queue = []
        for (size2, ancestor) in queue:
            if clade.issubset(ancestor):
                new_ancestor = (ancestor - clade) | frozenset([clade])
                counts[new_ancestor] = counts.pop(ancestor)
                edgelengths[new_ancestor] = edgelengths.pop(ancestor)
                ancestor = new_ancestor
            new_queue.append((len(ancestor), ancestor))
        children = [nodes.pop(c) for c in clade]
        assert len([children])
        nodes[clade] = tree_build(children, None, 
            {attr:counts[clade], 'length':edgelengths[clade]})
        queue = new_queue
    
    for root in nodes.values():
        root.Name = 'root' # Yuk
    
    return [root for root in nodes.values()]


if __name__ == "__main__":
    import sys
    trees = []
    for filename in sys.argv[1:]:
        for tree in open(filename):
            trees.append(LoadTree(treestring=tree))
    print "Consensus of %s trees from %s" % (len(trees),sys.argv[1:])
    outtrees = majorityRule(trees, strict=True)
    for tree in outtrees:
            print tree.asciiArt(compact=True, show_internal=False)
