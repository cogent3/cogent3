#!/usr/bin/env python
"""Leaf and Edge classes that can calculate their likelihoods.
Each leaf holds a sequence.  Used by a likelihood function."""

from __future__ import division
from cogent.util.modules import importVersionedModule, ExpectedImportError
from cogent.util.parallel import MPI
from cogent import LoadTable

import numpy

numpy.seterr(all='ignore')

numerictypes = numpy.core.numerictypes.sctype2char

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

try:
    pyrex = importVersionedModule('_likelihood_tree', globals(), 
            (2, 1), "pure Python/NumPy likelihoodihood tree")
except ExpectedImportError:
    pyrex = None
        
class _LikelihoodTreeEdge(object):
    def __init__(self, children, edge_name, alignment=None):
        self.edge_name = edge_name
        self.alphabet = children[0].alphabet
        self.comm = None  # for MPI 
        
        M = children[0].shape[-1]
        for child in children:
            assert child.shape[-1] == M
        
        # Unique positions are unique combos of input positions
        if alignment is None:
            # The children are pre-aligned gapped sequences
            assignments = [c.index for c in children]
        else:
            self.alignment = alignment  #XXX preserve through MPI split?
            # The children are ungapped sequences, 'alignment'
            # indicates where gaps need to go.
            assignments = []
            for (i, c) in enumerate(children):
                a = []
                for align_index in alignment:
                    col = align_index[i]
                    if col is None:
                        u = len(c.uniq)-1 # gap
                    else:
                        u = c.index[col]
                        assert 0 <= u < len(c.uniq)-1, (
                                u, len(c.uniq), c.uniq[-1], align_index)
                    a.append(u)
                assignments.append(a)
        (uniq, counts, self.index) = _indexed(zip(*assignments))
        
        # extra column for gap
        uniq.append(tuple([len(c.uniq)-1 for c in children]))
        counts.append(0)
        
        self.uniq = numpy.asarray(uniq, self.integer_type)
        
        # For faster math, a contiguous index array for each child
        self.indexes = [
                numpy.array(list(ch), self.integer_type)
                for ch in numpy.transpose(self.uniq)]

        # If this is the root it will need to weight the total
        # log likelihoods by these counts:
        self.counts = numpy.array(counts, self.float_type)
        
        # For product of child likelihoods
        self._indexed_children = zip(self.indexes, children)
        self.shape = [len(self.uniq), M]

        # Derive per-column degree of ambiguity from children's
        ambigs = [child.ambig[index] for (index, child) in self._indexed_children]
        self.ambig = numpy.product(ambigs, axis=0)
                
    def getSitePatterns(self, cols):
        # Recursive lookup of Site Patterns aka Alignment Columns
        child_motifs = [child.getSitePatterns(index[cols])
                for (index, child) in self._indexed_children]
        return [''.join(child[u] for child in child_motifs)
                for u in range(len(cols))]
    
    def restrictMotif(self, input_likelihoods, fixed_motif):
        # for reconstructAncestralSeqs
        mask = numpy.zeros([input_likelihoods.shape[-1]], self.float_type)
        mask[fixed_motif] = 1.0
        input_likelihoods *= mask
            
    def parallelShare(self, comm):
        """A local version of self for a single CPU in an MPI group"""
        if comm is None or comm.Get_size() == 1:
            return self
        assert self.comm is None
        U = len(self.uniq) - 1 # Gap column
        (size, rank) = (comm.Get_size(), comm.Get_rank())
        (share, remainder) = divmod(U, size)
        if share == 0:
            return self # not enough to share
        share_sizes = [share+1]*remainder + [share]*(size-remainder)
        assert sum(share_sizes) == U
        (lo,hi) = [sum(share_sizes[:i]) for i in (rank, rank+1)]
        local_cols = [i for (i,u) in enumerate(self.index) 
                if lo <= u < hi]
        local = self.selectColumns(local_cols)
        
        # Attributes for reconstructing/refinding the global arrays.
        # should maybe make a wrapping class instead.
        local.share_sizes = share_sizes
        local.comm = comm
        local.full_length_version = self
        return local
    
    def selectColumns(self, cols):
        children = []
        for (index, child) in self._indexed_children:
            child = child.selectColumns(cols)
            children.append(child)
        return self.__class__(children, self.edge_name)
    
    def parallelReconstructColumns(self, likelihoods):
        """Recombine full uniq array (eg: likelihoods) from MPI CPUs"""
        if self.comm is None:
            return (self, likelihoods)
        result = numpy.empty([sum(self.share_sizes)+1], likelihoods.dtype)
        mpi_type = None  # infer it from the likelihoods/result array dtype
        send = (likelihoods[:-1], mpi_type)  # drop gap column
        recv = (result, self.share_sizes, None, mpi_type)
        self.comm.Allgatherv(send, recv)
        result[-1] = likelihoods[-1]  # restore gap column
        return (self.full_length_version, result)
    
    def getFullLengthLikelihoods(self, likelihoods):
        (self, likelihoods) = self.parallelReconstructColumns(likelihoods)
        return likelihoods[self.index]
    
    def calcGStatistic(self, likelihoods, return_table=False):
        # A Goodness-of-fit statistic
        (self, likelihoods) = self.parallelReconstructColumns(likelihoods)
        
        unambig = (self.ambig == 1.0).nonzero()[0]
        observed = self.counts[unambig].astype(int)
        expected = likelihoods[unambig] * observed.sum()
        #chisq = ((observed-expected)**2 / expected).sum()
        G = 2 * observed.dot(numpy.log(observed/expected))
        
        if return_table:
            motifs = self.getSitePatterns(unambig)
            rows = zip(motifs, observed, expected)
            rows.sort(key=lambda row:(-row[1], row[0]))
            table = LoadTable(header=['Pattern', 'Observed', 'Expected'], rows=rows, row_ids=True)
            return (G, table)
        else:
            return G
    
    def getEdge(self, name):
        if self.edge_name == name:
            return self
        else:
            for (i,c) in self._indexed_children:
                r = c.getEdge(name)
                if r is not None:
                    return r
        return None
    
    def makePartialLikelihoodsArray(self):
        return numpy.ones(self.shape, self.float_type)
    
    def sumInputLikelihoods(self, *likelihoods):
        result = numpy.ones(self.shape, self.float_type)
        self.sumInputLikelihoodsR(result, *likelihoods)
        return result

    def asLeaf(self, likelihoods):
        (self, likelihoods) = self.parallelReconstructColumns(likelihoods)
        assert len(likelihoods) == len(self.counts)
        return LikelihoodTreeLeaf(likelihoods, likelihoods, 
                self.counts, self.index, self.edge_name, self.alphabet, None)

class _PyLikelihoodTreeEdge(_LikelihoodTreeEdge):
    # Should be a subclass of regular tree edge?

    float_type = numerictypes(float)
    integer_type = numerictypes(int)

    # For scaling very very small numbers
    BASE = 2.0 ** 100
    LOG_BASE = numpy.log(BASE)

    def sumInputLikelihoodsR(self, result, *likelihoods):
        result[:] = 1.0
        for (i, index) in enumerate(self.indexes):
            result *= numpy.take(likelihoods[i], index, 0)
        return result
    
    # For root
    
    def logDotReduce(self, patch_probs, switch_probs, plhs):
        (self, plhs) = self.parallelReconstructColumns(plhs)
        exponent = 0
        state_probs = patch_probs.copy()
        for site in self.index:
            state_probs = numpy.dot(switch_probs, state_probs) * plhs[site]
            while max(state_probs) < 1.0:
                state_probs *= self.BASE
                exponent -= 1
        return numpy.log(sum(state_probs)) + exponent * self.LOG_BASE

    def getTotalLogLikelihood(self, input_likelihoods, mprobs):
        lhs = numpy.inner(input_likelihoods, mprobs)
        return self.getLogSumAcrossSites(lhs)

    def getLogSumAcrossSites(self, lhs):
        return numpy.inner(numpy.log(lhs), self.counts)

class _PyxLikelihoodTreeEdge(_LikelihoodTreeEdge):
    integer_type = numerictypes(int)   # match checkArrayInt1D
    float_type = numerictypes(float)   # match checkArrayDouble1D/2D
    
    def sumInputLikelihoodsR(self, result, *likelihoods):
        pyrex.sumInputLikelihoods(self.indexes, result, likelihoods)
        return result
    
    # For root
    
    def logDotReduce(self, patch_probs, switch_probs, plhs):
        (self, plhs) = self.parallelReconstructColumns(plhs)
        return pyrex.logDotReduce(self.index, patch_probs, switch_probs, plhs)
        
    def getTotalLogLikelihood(self, input_likelihoods, mprobs):
        return pyrex.getTotalLogLikelihood(self.counts, input_likelihoods, 
                mprobs)
    
    def getLogSumAcrossSites(self, lhs):
        return pyrex.getLogSumAcrossSites(self.counts, lhs)

if pyrex is None:
    LikelihoodTreeEdge = _PyLikelihoodTreeEdge
else:
    LikelihoodTreeEdge = _PyxLikelihoodTreeEdge
    
FLOAT_TYPE = LikelihoodTreeEdge.float_type
INTEGER_TYPE = LikelihoodTreeEdge.integer_type

def _indexed(values):
    # >>> _indexed(['a', 'b', 'c', 'a', 'a'])
    # (['a', 'b', 'c'], [3, 1, 1], [0, 1, 2, 0, 0])
    index = numpy.zeros([len(values)], INTEGER_TYPE)
    unique = []
    counts = []
    seen = {}
    for (c, key) in enumerate(values):
        if key in seen:
            i = seen[key]
            counts[i] += 1
        else:
            i = len(unique)
            unique.append(key)
            counts.append(1)
            seen[key] = i
        index[c] = i
    return unique, counts, index

def makeLikelihoodTreeLeaf(sequence, alphabet=None, seq_name=None):    
    if alphabet is None:
        alphabet = sequence.MolType.Alphabet
    if seq_name is None:
        seq_name = sequence.getName()
        
    motif_len = alphabet.getMotifLen()
    sequence2 = sequence.getInMotifSize(motif_len)
    
    # Convert sequence to indexed list of unique motifs
    (uniq_motifs, counts, index) = _indexed(sequence2)
    
    # extra column for gap
    uniq_motifs.append('?' * motif_len)
    counts.append(0)
    
    counts = numpy.array(counts, FLOAT_TYPE)
    
    # Convert list of unique motifs to array of unique profiles
    try:
        likelihoods = alphabet.fromAmbigToLikelihoods(
            uniq_motifs, FLOAT_TYPE)
    except alphabet.AlphabetError, detail:
        motif = str(detail)
        posn = list(sequence2).index(motif) * motif_len
        raise ValueError, '%s at %s:%s not in alphabet' % (
                repr(motif), seq_name, posn)
    
    return LikelihoodTreeLeaf(uniq_motifs, likelihoods, 
                counts, index, seq_name, alphabet, sequence)

class LikelihoodTreeLeaf(object):
    def __init__(self, uniq, likelihoods, counts, index, edge_name, 
            alphabet, sequence):
        if sequence is not None:
            self.sequence = sequence 
        self.alphabet = alphabet
        self.name = self.edge_name = edge_name
        self.uniq = uniq
        self.motifs = numpy.asarray(uniq)
        self.input_likelihoods = likelihoods
        self.counts = counts
        self.index = index
        self.shape = likelihoods.shape
        self.ambig = numpy.sum(self.input_likelihoods, axis=-1)
    
    def backward(self):
        index = numpy.array(self.index[::-1,...])
        result = self.__class__(self.uniq, self.input_likelihoods, self.counts, 
                index, self.edge_name, self.alphabet, None)
        return result

    def __len__(self):
        return len(self.index)
        
    def __getitem__(self, index):
        cols = range(*index.indices(len(self.index)))
        return self.selectColumns(cols)
        
    def getMotifCounts(self, include_ambiguity=False):
        weights = self.counts / self.ambig
        profile = self.input_likelihoods * weights[...,numpy.newaxis]
        if not include_ambiguity:
            unambig = self.ambig == 1.0
            profile = numpy.compress(unambig, profile, axis=0)
        return numpy.sum(profile, axis=0)
    
    def getAmbiguousPositions(self):
        ambig = {}
        for (i,u) in enumerate(self.index):
            if self.ambig[u] != 1.0:
                ambig[i] = self.uniq[u]
        return ambig
    
    def selectColumns(self, cols):
        sub_index = [self.index[i] for i in cols]
        (keep, counts, index) = _indexed(sub_index)
        keep.append(len(self.uniq)-1)  # extra column for gap
        counts.append(0)
        counts = numpy.array(counts, FLOAT_TYPE)
        uniq = [self.uniq[u] for u in keep]
        likelihoods = self.input_likelihoods[keep]
        return self.__class__(
                uniq, likelihoods, counts, index, self.edge_name, 
                self.alphabet, None)
        
    def getEdge(self, name):
        if self.edge_name == name:
            return self
        else:
            return None
            
    def getSitePatterns(self, cols):
        return numpy.asarray(self.uniq)[cols]
        
