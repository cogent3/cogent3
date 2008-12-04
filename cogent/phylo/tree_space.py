#!/usr/bin/env python
import heapq
import numpy
from cogent.core.tree import TreeBuilder
from cogent.util import parallel, checkpointing

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# Trees are represented as "ancestry" matricies in which A[i,j] iff j is an
# ancestor of i.  For LS calculations the ancestry matrix is converted
# to a "paths" matrix or "split metric" in which S[p,j] iff the path between
# the pth pair of tips passes through edge j.  For ML calculations the
# ancestry matrix is converted into an ordinary cogent tree object.

def tree2ancestry(tree):
    nodes = tree.getEdgeVector()[:-1]
    n = len(nodes)
    A = numpy.zeros([n, n], int)
    seen = []
    for (i, node) in enumerate(nodes):
        seen.append((i, node))
        A[i, i] = 1
        for c in node.Children:
            for (j,o) in seen:
                if o is c:
                    break
            else:
                raise RuntimeError
            A[:,i] |= A[:,j]
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
    """Given 'tree' B return tree with one extra edge"""
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
    
    def trex(self, a=8, k=1000, show_progress=False,
            return_all=False, checkpoint_filename=None):
        """TrexML policy for tree sampling - all trees up to size 'a' and
        then keep no more than 'k' best trees at each tree size
        
        Advanced step-wise addition algorithm
        M. J. Wolf, S. Easteal, M. Kahn, B. D. McKay, and L. S. Jermiin.
        Trexml: a maximum-likelihood approach for extensive tree-space
        exploration.
        Bioinformatics, 16(4):383 94, 2000."""
        
        printing_cpu = parallel.getCommunicator().rank == 0
        
        names = self.names
        tree_size = len(names)
        assert tree_size > 3
        
        checkpointer = checkpointing.Checkpointer(checkpoint_filename)
        if checkpointer.available():
            (a, old_names, trees) = checkpointer.load()
            assert sorted(old_names) == sorted(names[:a])
            names = old_names + names[a:]
        else:
            # All trees of size a-1, no need to compare them
            if a > tree_size:
                a = tree_size
            if a < 4:
                a = 4
            trees = [(None, None, numpy.array([[1,0,0], [0,1,0], [0,0,1]]))]
            for n in range(4, a):
                trees2 = []
                for (err2, lengths2, ancestry) in trees:
                    for split_edge in range(len(ancestry)):
                        ancestry2 = grown(ancestry, split_edge)
                        trees2.append((None, None, ancestry2))
                trees = trees2
        
        if show_progress and printing_cpu:
            print len(trees), 'trees of size', a-1, 'at start'
        
        for n in range(a, tree_size+1):
            # parallel - combine tree heaps up a tree of cpus
            evals = len(trees) * len(trees[0][2])
            if show_progress and printing_cpu:
                print evals, 'trees of size', n, '...',
            
            (comm, leftover) = parallel.getSplitCommunicators(evals)
            parallel.push(leftover)
            try:
                # tree of cpus
                parent = (comm.rank+1) // 2 - 1
                left = (comm.rank+1) * 2 - 1
                children = [child for child in [left, left+1]
                        if child < comm.size]
                
                trees2 = []
                evaluate = self.makeTreeScorer(names[:n])
                evaluation_count = 0
                for (err2, lengths2, ancestry) in trees:
                    for split_edge in range(len(ancestry)):
                        evaluation_count += 1
                        if evaluation_count % comm.size != comm.rank:
                            continue
                        ancestry2 = grown(ancestry, split_edge)
                        (err, lengths) = evaluate(ancestry2)
                        trees2.append((err, lengths, ancestry2))
                
                for child in children:
                    trees2.extend(comm.receive_obj(child))
                
                if len(trees2) > k and (n < tree_size or not return_all):
                    trees2 = heapq.nsmallest(k, trees2)
                
                if show_progress and printing_cpu:
                    if n < tree_size:
                        print 'kept', len(trees2), 'trees'
                    else:
                        print 'done'
                
                if comm.rank > 0:
                    comm.send_obj(trees2, parent)
                    trees2 = comm.receive_obj(parent)
                
                for child in children:
                    comm.send_obj(trees2, child)
            finally:
                parallel.pop(leftover)
            trees = trees2
            if printing_cpu:
                checkpointer.record((n, names[:n], trees))
        
        results = self._result_iter(trees, names)
        if return_all:
            result = list(results)
        else:
            result = results.next()
        return result
    
