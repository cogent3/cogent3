#!/usr/bin/env python
"""Neighbour joining phylogenetic tree estimation.

Note that negative branch lengths are reset to 0.0 during the calculations.

This is the algorithm of Studier and Keppler, as described in the book
Biological sequence analysis by Durbin et al.
"""

import numpy
from cogent.core.tree import TreeBuilder
from util import distanceDictTo2D

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def nj(dists, no_negatives=True):
    """Arguments:
        - dists: dict of (name1, name2): distance
        - no_negatives: negative branch lengths will be set to 0
    """
    
    constructor = TreeBuilder(mutable=True).createEdge
    (names, d) = distanceDictTo2D(dists)
    
    nodes = [constructor([], name, {}) for name in names]
    
    while len(nodes) > 2:
        # Eliminate one node per iteration until 2 left
        L = len(nodes)
        
        # Find neighbours i and j
        r = numpy.sum(d, 0) * 1.0/(len(nodes)-2)
        D = d - numpy.add.outer(r, r)
        D += numpy.identity(L) * max(numpy.absolute(D.ravel())) * 2
        (i, j) = divmod(numpy.argmin(D.ravel()), L)
        assert i != j, (i, j, D[i, j])
        
        # Branch lengths from i and j to new node
        nodes[i].Length = 0.5 * (d[i,j] + r[i] - r[j])
        nodes[j].Length = 0.5 * (d[i,j] + r[j] - r[i])
        
        # no negative branch lengths
        if no_negatives:
            nodes[i].Length = max(0.0, nodes[i].Length)
            nodes[j].Length = max(0.0, nodes[j].Length)
        
        # Join i and k to make new node
        new_node = constructor([nodes[i], nodes[j]], None, {})
        
        # Store new node at i
        new_dists = 0.5 * (d[i] + d[j] - d[i,j])
        d[:, i] = new_dists
        d[i, :] = new_dists
        d[i, i] = 0.0
        nodes[i] = new_node
        
        # Eliminate j
        d[j, :] = d[L-1, :]
        d[:, j] = d[:, L-1]
        assert d[j, j] == 0.0, d
        d = d[0:L-1, 0:L-1]
        nodes[j] = nodes[L-1]
        nodes.pop()
    
    # no negative branch lengths
    if len(nodes[0].Children) < len(nodes[1].Children):
        nodes.reverse()
    
    # 2 remaining nodes will be [root, extra_child]
    nodes[1].Length = d[0,1]
    if no_negatives:
        nodes[1].Length = max(0.0, nodes[1].Length)
    
    #Need to replace nodes[0] with new root
    nodes[1].Parent = nodes[0]
    return constructor(nodes[0].Children, 'root', {}).deepcopy()
