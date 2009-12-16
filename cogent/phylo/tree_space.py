#!/usr/bin/env python
import heapq
import numpy
from cogent.core.tree import TreeBuilder
from cogent.util import parallel, checkpointing

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# Trees are represented as "ancestry" matricies in which A[i,j] iff j is an
# ancestor of i.  For LS calculations the ancestry matrix is converted
# to a "paths" matrix or "split metric" in which S[p,j] iff the path between
# the pth pair of tips passes through edge j.  For ML calculations the
# ancestry matrix is converted back into an ordinary cogent tree object.

def tree2ancestry(tree, order=None):
    nodes = tree.unrooted().getEdgeVector()[:-1]
    if order is not None:
        lookup = dict([(k,i) for (i,k) in enumerate(order)])
        def _ordered_tips_first(n):
            if n.Children:
                return len(order)
            else:
                return lookup[n.Name]
        nodes.sort(key=_ordered_tips_first)

    n = len(nodes)
    A = numpy.zeros([n, n], int)
    seen = {}
    for (i, node) in enumerate(nodes):
        A[i, i] = 1
        seen[id(node)] = i
        for c in node.Children:
            A[:,i] |= A[:,seen[id(c)]]
    names = [n.Name for n in nodes if not n.Children]
    lengths = [n.Length for n in nodes]
    return (A, names, lengths)

def ancestry2tree(A, lengths, tip_names):
    """Convert edge x edge ancestry matrix to a cogent Tree object"""
    tips = {}
    tip = 0
    for i in range(len(A)):
        if numpy.sum(A[:,i]) == 1:
            tips[i] = tip_names[tip]
            tip += 1
    assert tip == len(tip_names)
    
    constructor = TreeBuilder().createEdge
    free = {}
    for i in numpy.argsort(numpy.sum(A, axis=0)):
        children = [j for j in range(len(A)) if A[j, i] and j != i]
        child_nodes = [free.pop(j) for j in children if j in free]
        if child_nodes:
            name = None
        else:
            name = tips[i]
        if lengths is None:
            params = {}
        else:
            params = {'length':lengths[i]}
        node = constructor(child_nodes, name, params)
        free[i] = node
    return constructor(free.values(), 'root', {})

def grown(B, split_edge):
    """Ancestry matrix 'B' with one extra leaf added at 'split_edge'.
    Row/column order within the matrix is independent of the topology it 
    represents. The added leaf will be the last one in the matrix, which keeps 
    the leaf node order the same as the order in which they are added, which is 
    what is assumed by ancestry2tree and ancestry2paths"""
    n = len(B)
    A = numpy.zeros([n+2, n+2], int)
    A[:n, :n] = B
    (sibling, parent) = (n, n + 1)
    A[sibling] = A[parent] = A[split_edge]
    A[:,parent] = A[:,split_edge]
    A[sibling,split_edge] = 0
    A[parent, split_edge] = 0
    A[sibling,sibling] = 1
    A[parent,parent] = 1
    A[sibling,parent] = 1
    A[split_edge,parent] = 1
    return A

class TreeEvaluator(object):
    """Subclass must provide makeTreeScorer and result2output"""
    
    def _result_iter(self, trees, names):
        heapq.heapify(trees)
        while trees:
            (err, lengths, ancestry) = heapq.heappop(trees)
            yield self.result2output(err, ancestry, lengths, names)
    
    def results2output(self, results):
        return list(results)
        
    def evaluateTopology(self, tree):
        """Optimal (score, tree) for the one topology 'tree'"""
        (ancestry, names, lengths) = tree2ancestry(tree)
        evaluate = self.makeTreeScorer(names)
        (err, lengths) = evaluate(ancestry)
        return self.result2output(err, ancestry, lengths, names)
    
    def evaluateTree(self, tree):
        """score for 'tree' with lengths as-is"""
        (ancestry, names, lengths) = tree2ancestry(tree)
        evaluate = self.makeTreeScorer(names)
        (err, result) = evaluate(ancestry, lengths=lengths)
        return err
    
    def _consistentNameOrder(self, fixed_names, ordered_names=None):
        """fixed_names followed by ordered_names without duplicates"""
        all_names = set(self.names)

        fixed_names_set = set(fixed_names)
        assert fixed_names_set.issubset(all_names)

        if ordered_names:
            assert set(ordered_names).issubset(all_names)
        else:
            ordered_names = self.names
        names = list(fixed_names) + [n for n in ordered_names 
                if n not in fixed_names_set]
        return names
    
    def trex(self, a=8, k=1000, start=None, order=None, show_progress=False,
            return_all=False, checkpoint_filename=None):
        """TrexML policy for tree sampling - all trees up to size 'a' and
        then keep no more than 'k' best trees at each tree size.
        'order' is an optional list of tip names.  
        'start' is an optional list of initial trees.  Each of the trees must
        contain the same tips.
        
        Advanced step-wise addition algorithm
        M. J. Wolf, S. Easteal, M. Kahn, B. D. McKay, and L. S. Jermiin.
        Trexml: a maximum-likelihood approach for extensive tree-space
        exploration.
        Bioinformatics, 16(4):383 94, 2000."""
        
        printing_cpu = parallel.getCommunicator().Get_rank() == 0
        
        checkpointer = checkpointing.Checkpointer(checkpoint_filename)
        if checkpointer.available():
            (init_tree_size, fixed_names, trees) = checkpointer.load()
            names = self._consistentNameOrder(fixed_names, order)
        elif start is not None:
            if not isinstance(start, list):
                start = [start]
            fixed_names = start[0].getTipNames()
            names = self._consistentNameOrder(fixed_names, order)
            trees = []
            for tree in start:
                # check the start tree represents a subset of tips
                assert set(tree.getTipNames()) < set(self.names), \
                    "Starting tree names not a subset of the sequence names"
                
                (ancestry, fixed_names2, lengths) = tree2ancestry(
                        tree, order=fixed_names)
                assert fixed_names2 == fixed_names
                trees.append((None, None, ancestry))
            init_tree_size = len(fixed_names)
        else:
            trees = [(None, None, numpy.identity(3, int))]
            names = self._consistentNameOrder([], order)
            init_tree_size = 3
            
        tree_size = len(names)
        assert tree_size > 3
        if a > tree_size:
            a = tree_size
        if a < 4:
            a = 4

        # All trees of size a-1, no need to compare them
        for n in range(init_tree_size+1, a):
            trees2 = []
            for (err2, lengths2, ancestry) in trees:
                for split_edge in range(len(ancestry)):
                    ancestry2 = grown(ancestry, split_edge)
                    trees2.append((None, None, ancestry2))
            trees = trees2
            init_tree_size = n
        
        if show_progress and printing_cpu:
            print len(trees), 'trees of size', init_tree_size, 'at start'
        
        for n in range(init_tree_size+1, tree_size+1):
            # parallel - combine tree heaps up a tree of cpus
            evals = len(trees) * len(trees[0][2])
            if show_progress and printing_cpu:
                print evals, 'trees of size', n, '...',
            
            (comm, leftover) = parallel.getSplitCommunicators(evals)
            parallel.push(leftover)
            (cpu_count, cpu_rank) = (comm.Get_size(), comm.Get_rank())
            try:
                # tree of cpus
                parent = (cpu_rank+1) // 2 - 1
                left = (cpu_rank+1) * 2 - 1
                children = [child for child in [left, left+1]
                        if child < cpu_count]
                
                trees2 = []
                evaluate = self.makeTreeScorer(names[:n])
                evaluation_count = 0
                for (err2, lengths2, ancestry) in trees:
                    for split_edge in range(len(ancestry)):
                        evaluation_count += 1
                        if evaluation_count % cpu_count != cpu_rank:
                            continue
                        ancestry2 = grown(ancestry, split_edge)
                        (err, lengths) = evaluate(ancestry2)
                        trees2.append((err, lengths, ancestry2))
                
                for child in children:
                    trees2.extend(comm.recv(source=child))
                
                if len(trees2) > k and (n < tree_size or not return_all):
                    trees2 = heapq.nsmallest(k, trees2)
                
                if show_progress and printing_cpu:
                    if n < tree_size:
                        print 'kept', len(trees2), 'trees'
                    else:
                        print 'done'
                
                if cpu_rank > 0:
                    comm.send(trees2, parent)
                    trees2 = comm.recv(source=parent)
                
                for child in children:
                    comm.send(trees2, child)
            finally:
                parallel.pop(leftover)
            trees = trees2
            if printing_cpu:
                checkpointer.record((n, names[:n], trees))
        
        results = self._result_iter(trees, names)
        if return_all:
            result = self.results2output(results)
        else:
            result = results.next()
        return result
    
