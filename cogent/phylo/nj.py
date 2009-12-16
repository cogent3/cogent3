#!/usr/bin/env python
"""Generalised Neighbour Joining phylogenetic tree estimation.

By default negative branch lengths are reset to 0.0 during the calculations.

This is based on the algorithm of Studier and Keppler, as described in the book
Biological sequence analysis by Durbin et al

Generalised as described by Pearson, Robins & Zhang, 1999.
"""

import numpy
from cogent.core.tree import TreeBuilder
from util import distanceDictTo2D

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class LightweightTreeTip(str):
    def convert(self, constructor, length):
        node = constructor([], str(self), {})
        node.Length = max(0.0, length)
        return node
        
class LightweightTreeNode(frozenset):
    """Set of (length, child node) tuples"""
    def convert(self, constructor=None, length=None):
        if constructor is None:
            constructor = TreeBuilder().createEdge
        children = [child.convert(constructor, clength) 
                for (clength, child) in self]
        node = constructor(children, None, {})
        if length is not None:
            node.Length = max(0.0, length)
        return node
        
    def __or__(self, other):
        return type(self)(frozenset.__or__(self, other))
        
class PartialTree(object):
    # A candidate tree stored as
    #  (distance matrix, list of subtrees, list of tip sets, set of partitions, score)
    #  at each iteration (ie: call of the join method) the number of subtrees is 
    #  reduced as 2 of them are joined, while
    #  the number of partitions is increased as a new edge is introduced.
    def __init__(self, d, nodes, tips, topology, score, encode_partition):
        self.d = d
        self.nodes = nodes
        self.tips = tips
        self.topology = topology
        self.score = score
        self._encode_partition = encode_partition
        
    def getDistSavedJoinScoreMatrix(self):
        d = self.d
        L = len(d)
        r = numpy.sum(d, 0)
        Q = d - numpy.add.outer(r, r)/(L-2.0)
        return Q/2.0 + sum(r)/(L-2.0)/2 + self.score

    def join(self, i, j, topologies=None):
        tips = self.tips[:]
        new_tip_set = tips[i] | tips[j]
        # check is new topology
        topology = self.topology | frozenset(
                [self._encode_partition(new_tip_set)])
        if topologies is not None and topology in topologies:
            return None
        
        nodes = self.nodes[:]
        d = self.d.copy()

        # Branch lengths from i and j to new node
        L = len(nodes)
        r = numpy.sum(d, axis=0)
        ij_dist_diff = (r[i]-r[j]) / (L-2.0)
        left_length = 0.5 * (d[i,j] + ij_dist_diff)
        right_length = 0.5 * (d[i,j] - ij_dist_diff)

        score = self.score + d[i,j]
        
        left_length = max(0.0, left_length)
        right_length = max(0.0, right_length)
        
        # Join i and k to make new node
        new_node = LightweightTreeNode(
                [(left_length, nodes[i]), (right_length, nodes[j])])
        
        # Store new node at i
        new_dists = 0.5 * (d[i] + d[j] - d[i,j])
        d[:, i] = new_dists
        d[i, :] = new_dists
        d[i, i] = 0.0
        nodes[i] = new_node
        tips[i] = new_tip_set
        
        # Eliminate j
        d[j, :] = d[L-1, :]
        d[:, j] = d[:, L-1]
        assert d[j, j] == 0.0, d
        d = d[0:L-1, 0:L-1]
        nodes[j] = nodes[L-1]
        nodes.pop()
        tips[j] = tips[L-1]
        tips.pop()
        
        return type(self)(d, nodes, tips, topology, score, self._encode_partition)
    
    def asScoreTreeTuple(self):
        assert len(self.nodes) == 2
        nodes = self.nodes
        if not isinstance(nodes[0], frozenset):
            nodes = nodes[::-1]
        length = max(self.d[0,1], 0.0)
        root = nodes[0] | LightweightTreeNode([(length, nodes[1])])
        tree = root.convert()
        tree.Name = "root"
        return (self.score + length, tree)

def gnj(dists, keep=None, show_progress=False):
    """Arguments:
        - dists: dict of (name1, name2): distance
        - keep: maximum number of partial trees to keep at each iteration, and 
          therefore to return.  Same as Q parameter in original GNJ paper.
    Result:
        - a sorted list of (tree length, tree) tuples
    """
        
    (names, d) = distanceDictTo2D(dists)

    if keep is None:
        keep = len(names) * 5
        
    # For recognising duplicate topologies, encode partitions (ie: edges) as 
    # frozensets of tip names, which should be quickly comparable.
    arbitrary_anchor = names[0]
    all_tips = frozenset(names)
    def encode_partition(tips):
        included = frozenset(tips)
        if arbitrary_anchor not in included:
            included = all_tips - included
        # could also convert to long int, would be faster?
        return included
        
    tips = [frozenset([n]) for n in names]
    star_topology = frozenset(encode_partition(tip) for tip in tips)
    nodes = [LightweightTreeTip(name) for name in names]
    star_tree = PartialTree(d, nodes, tips, star_topology, 0.0, encode_partition)
    trees = [star_tree]
    
    for L in range(len(names), 2, -1):
        scores = numpy.zeros([len(trees), L, L])
        # scores[k,i,j] is score of the kth tree with its ith and jth nodes joined.
        for (k, tree) in enumerate(trees):
            scores[k] = tree.getDistSavedJoinScoreMatrix()
        
        next_trees = []
        topologies = set()
        order = numpy.argsort(scores.flat)
        for index in order:
            if len(next_trees) == keep:
                break
            (k, ij) = divmod(index, L*L)
            (i, j) = divmod(ij, L)
            if i == j:
                continue
            
            tree = trees[k].join(i,j, topologies)
            if tree is None: # ie: a duplicate topology
                continue
            
            topologies.add(tree.topology)
            next_trees.append(tree)
            #if len(next_trees) > 1:
            #    print len(next_trees[0].topology ^ tree.topology)
        
        trees = next_trees
        
        if show_progress:
            print L
        
    result = []
    for tree in trees:
        result.append(tree.asScoreTreeTuple())
    result.sort()
    
    return result

def nj(dists, no_negatives=True):
    """Arguments:
        - dists: dict of (name1, name2): distance
        - no_negatives: negative branch lengths will be set to 0
    """
    assert no_negatives, "no_negatives=False is deprecated"
    (result,) = gnj(dists, keep=1)
    (score, tree) = result
    return tree

