#!/usr/bin/env python

"""Align two Alignables, each of which can be a sequence or a subalignment
produced by the same code."""

# How many cells before using linear space alignment algorithm.
# Should probably set to about half of physical memory / PointerEncoder.bytes
HIRSCHBERG_LIMIT = 10**8

import numpy

# setting global state on module load is bad practice and can be ineffective,
# and this setting matches the defaults anyway, but left here as reference 
# in case it needs to be put back in a more runtime way.
#numpy.seterr(all='ignore')

import warnings

from cogent.align.traceback import alignment_traceback
from cogent.evolve.likelihood_tree import LikelihoodTreeEdge
from indel_positions import leaf2pog
from cogent import LoadSeqs
from cogent.core.alignment import Aligned
from cogent.align.traceback import map_traceback
from cogent.util import parallel

from cogent.util.modules import importVersionedModule, ExpectedImportError
try:
    pyrex_align_module = importVersionedModule('_pairwise_pogs', globals(),
            (3, 1), "slow Python alignment implementation")
except ExpectedImportError:
    pyrex_align_module = None
try:
    pyrex_seq_align_module = importVersionedModule('_pairwise_seqs', globals(),
            (3, 1), "slow Python alignment implementation")
except ExpectedImportError:
    pyrex_seq_align_module = None

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class PointerEncoding(object):
    """Pack very small ints into a byte.  The last field, state, is assigned
    whatever bits are left over after the x and y pointers have claimed what 
    they need, which is expected to be only 2 bits each at most"""
    
    dtype = numpy.int8
    bytes = 1

    def __init__(self, x, y):
        assert x > 0 and y > 0, (x,y)
        (x, y) = (numpy.ceil(numpy.log2([x+1,y+1]))).astype(int)
        s = 8 * self.bytes - sum([x, y])
        assert s**2 >= 4+1, (x,y,s) # min states required
        self.widths = numpy.array([x,y,s]).astype(int)
        self.limits = 2 ** self.widths
        self.max_states = self.limits[-1]
        if DEBUG:
            print self.max_states, "states allowed in viterbi traceback"
        self.positions = numpy.array([0, x, x+y], int)
        #a.flags.writeable = False
    def encode(self, x, y, s):
        parts = numpy.asarray([x, y, s], int)
        assert all(parts < self.limits), (parts, self.limits)
        return (parts << self.positions).sum()
    def decode(self, coded):
        return (coded >> self.positions) % self.limits
    def getEmptyArray(self, shape):
        return numpy.zeros(shape, self.dtype)
    
DEBUG = False

def py_calc_rows(plan, x_index, y_index, i_low, i_high, j_low, j_high,
        preds, state_directions, T,
        xgap_scores, ygap_scores, match_scores, rows, track, track_enc, 
        viterbi, local=False, use_scaling=False, use_logs=False):
    """Pure python version of the dynamic programming algorithms
    Forward and Viterbi.  Works on sequences and POGs.  Unli"""
    if use_scaling:
        warnings.warn("Pure python version of DP code can suffer underflows")
        # because it ignores 'exponents', the Pyrex version doesn't.
    source_states = range(len(T))
    BEGIN = 0
    ERROR = len(T)
    (rows, exponents) = rows
    if use_logs:
        neutral_score = 0.0
        impossible = -numpy.inf
    else:
        neutral_score = 1.0
        impossible = 0.0
    best_score = impossible
    for i in range(i_low, i_high):
        x = x_index[i]
        i_sources = preds[0][i]
        current_row = rows[plan[i]]
        #current_row[:] = 0.0
        current_row[:,0] = impossible
        if i == 0 and not local:
            current_row[0,0] = neutral_score
        for j in range(j_low, j_high):
            y = y_index[j]
            j_sources = preds[1][j]
            for (state, bin, dx, dy) in state_directions:
                if (local and dx and dy):
                    cumulative_score = T[BEGIN, state]
                    pointer = (dx, dy, BEGIN)
                else:
                    cumulative_score = impossible
                    pointer = (0, 0, ERROR)
                for (a, prev_i) in enumerate([[i], i_sources][dx]):
                    source_row = rows[plan[prev_i]]
                    for (b, prev_j) in enumerate([[j], j_sources][dy]):
                        source_posn = source_row[prev_j]
                        for prev_state in source_states:
                            prev_value = source_posn[prev_state]
                            transition = T[prev_state, state]
                            if viterbi:
                                if use_logs:
                                    candidate = prev_value + transition
                                else:
                                    candidate = prev_value * transition
                                #if DEBUG:
                                #    print prev_state, prev_value, state
                                if candidate > cumulative_score:
                                    cumulative_score = candidate
                                    pointer = (a+dx, b+dy, prev_state)
                            else:
                                cumulative_score += prev_value * transition
                if dx and dy:
                    d_score = match_scores[bin, x, y]
                elif dx:
                    d_score = xgap_scores[bin, x]
                elif dy:
                    d_score = ygap_scores[bin, y]
                else:
                    d_score = neutral_score
                #if DEBUG:
                #    print (dx, dy), d_score, cumulative_score
                if use_logs:
                    current_row[j, state] = cumulative_score + d_score
                else:
                    current_row[j, state] = cumulative_score * d_score
                if track is not None:
                    track[i,j,state] = (numpy.array(pointer) << track_enc).sum()
                if (i==i_high-1 and j==j_high-1 and not local) or (
                        local and dx and dy and current_row[j, state] > best_score):
                    (best, best_score) = (((i, j), state), current_row[j, state])
    #if DEBUG:
    #    print i_low, i_high, j_low, j_high
    #    print 'best_score %5.1f  at  %s' % (numpy.log(best_score), best)
    if not use_logs:
        best_score = numpy.log(best_score)
    return best + (best_score,)

class TrackBack(object):
    def __init__(self, tlist):
        self.tlist = tlist
    
    def __str__(self):
        return ''.join('(%s,%s)%s' % 
            (x,y, '.xym'[dx+2*dy]) for (state, (x,y), (dx,dy))
            in self.tlist)
            
    def offset(self, X, Y):
        tlist = [(state, (x+X, y+Y), dxy) for (state, (x,y), dxy) in self.tlist]
        return TrackBack(tlist)
                
    def __add__(self, other):
        return TrackBack(self.tlist + other.tlist)
        
    #def asStatePosnTransTuples(self):
    #    return iter(self.tlist)
        
    def asStatePosnTuples(self):
        return [(s,p) for (s,p,d) in self.tlist]

    def asBinPosTuples(self, state_directions):
        bin_map = dict((state, bin) for (state, bin, dx, dy) in
            state_directions)
        result = []
        for (state, posn, (dx, dy)) in self.tlist:
            pos = [[None, i-1][d] for (i,d) in zip(
                    posn, [dx, dy])]
            result.append((bin_map.get(int(state), None), pos))
        return result
    

class Pair(object):
    def __init__(self, alignable1, alignable2, backward=False):
        alignables = [alignable1, alignable2]
        assert alignable1.alphabet == alignable2.alphabet
        self.alphabet = alignable1.alphabet
        
        for alignable in alignables:
            assert isinstance(alignable, _Alignable), type(alignable)
            if not isinstance(alignable, AlignableSeq):
                some_pogs = True
                break
        else:
            some_pogs = False
        
        if some_pogs and pyrex_align_module is not None:
            aligner = pyrex_align_module.calc_rows
        elif (not some_pogs) and pyrex_seq_align_module is not None:
            aligner =  pyrex_seq_align_module.calc_rows
        else:
            aligner = py_calc_rows
        
        self.both_seqs = not some_pogs
        self.aligner = aligner
        
        if backward:
            alignables = [a.backward() for a in alignables]
        
        self.children = [alignable1, alignable2] = alignables
        self.max_preds = [alignable.max_preds for alignable in alignables]
        self.pointer_encoding = PointerEncoding(*self.max_preds)

        self.size = [len(alignable1), len(alignable2)]
        self.uniq_size = [len(alignable1.plh), len(alignable2.plh)]
        self.plan = numpy.array(alignable1.getRowAssignmentPlan())
        self.x_index = alignable1.index
        self.y_index = alignable2.index
    
    def getSeqNamePairs(self):
        return [(a.leaf.edge_name, a.leaf.sequence) for a in self.children]
    
    def makeSimpleEmissionProbs(self, mprobs, psubs1):
        psubs2 = [numpy.identity(len(psub)) for psub in psubs1]
        bins = [PairBinData(mprobs, *ppsubs) for ppsubs in zip(
            psubs1, psubs2) ]
        return PairEmissionProbs(self, bins)
    
    def makeEmissionProbs(self, bins):
        bins = [PairBinData(*args) for args in bins]
        return PairEmissionProbs(self, bins)
    
    def makeReversibleEmissionProbs(self, bins, length):
        bins = [BinData(*bin) for bin in bins]
        return ReversiblePairEmissionProbs(self, bins, length)
    
    def backward(self):
        return Pair(*self.children, **dict(backward=True))
    
    def __getitem__(self, index):
        assert len(index) == 2, index
        children = [child[dim_index] for (child, dim_index) in zip(
                self.children, index)]
        return Pair(*children)
        
    def _decode_state(self, track, encoding, posn, pstate):
        coded = int(track[posn[0], posn[1], pstate])
        (a, b, state) = encoding.decode(coded)
        if state >= track.shape[-1]:
            raise ArithmeticError('Error state in traceback')
        (x, y) = posn
        if state == -1:
            next = (x, y)
        else:
            if a: x = self.children[0][x][a-1]
            if b: y = self.children[1][y][b-1]
            next = numpy.array([x,y], int)
        return (next, (a, b), state)
    
    def traceback(self, track, encoding, posn, state, skip_last=False):
        result = []
        started = False
        while 1:
            (nposn, (a, b), nstate) = self._decode_state(track, encoding, 
                    posn, state)
            if state:
                result.append((state, posn, (a>0, b>0)))
            if started and state == 0:
                break
            (posn, state) = (nposn, nstate)
            started = True
        result.reverse()
        if skip_last:
            result.pop()
        return TrackBack(result)
    
    def edge2plh(self, edge, plhs):
        bins = plhs[0].shape[0]
        plh = [edge.sumInputLikelihoods(*[p[bin][1:-1] for p in plhs])
                for bin in range(bins)]
        return plh
    
    def getPOG(self, aligned_positions):
        (pog1, pog2) = [child.getPOG() for child in self.children]
        return pog1.traceback(pog2, aligned_positions)
        
    def getPointerEncoding(self, n_states):
        assert n_states <= self.pointer_encoding.max_states, (
            n_states, self.pointer_encoding.max_states)
        return self.pointer_encoding

    def getScoreArraysShape(self):
        needed = max(self.plan) + 1
        N = self.size[1]
        return (needed, N)
    
    def getEmptyScoreArrays(self, n_states, dp_options):
        shape = self.getScoreArraysShape() + (n_states,)
        mantissas = numpy.zeros(shape, float)
        if dp_options.use_logs:
            mantissas[:] = numpy.log(0.0)
        if dp_options.use_scaling:
            exponents = numpy.ones(shape, int) * -10000
        else:
            exponents = None
        return (mantissas, exponents)
    
    def calcRows(self, i_low, i_high, j_low, j_high, state_directions,
            T, scores, rows, track, track_encoding, viterbi, **kw):
        (match_scores, (xscores, yscores)) = scores
        track_enc = track_encoding and track_encoding.positions
        #print T
        return self.aligner(self.plan, self.x_index, self.y_index,
                i_low, i_high, j_low, j_high, self.children, state_directions,
                T, xscores, yscores, match_scores, rows, track, 
                track_enc, viterbi, **kw)
    

class _Alignable(object):
    def __init__(self, leaf):
        self.leaf = leaf
        self.alphabet = leaf.alphabet
        (uniq, alphabet_size) = leaf.input_likelihoods.shape
        full = len(leaf.index)
        self.plh = numpy.zeros([uniq+2, alphabet_size], float)
        self.plh[1:-1] = leaf.input_likelihoods
        self.index = numpy.zeros([full+2], int)
        self.index[1:-1] = numpy.asarray(leaf.index) + 1
        self.index[0] = 0
        self.index[full+1] = uniq+1
    
    def _asCombinedArray(self):
        # POG in a format suitable for Pyrex code, two arrays
        # preds here means predecessor
        pred = []
        offsets = []
        for pre in self:
            offsets.append(len(pred))
            pred.extend(pre)
        offsets.append(len(pred))
        # provides the paths leading to a point (predecessors), and offsets
        # records index positions fdor each point (graph node)
        return (numpy.array(pred), numpy.array(offsets))
    
    def asCombinedArray(self):
        if not hasattr(self, '_combined'):
            self._combined = self._asCombinedArray()
        return self._combined
    
    def getRowAssignmentPlan(self):
        d = self.getOuterLoopDiscardPoints()
        free = set()
        top = 0
        assignments = []
        for i in range(len(d)):
            if free:
                assignments.append(free.pop())
            else:
                assignments.append(top)
                top = top + 1
            for j in d[i]:
                free.add(assignments[j])
        
        return assignments
    

class AlignablePOG(_Alignable):
    """Alignable wrapper of a Partial Object Graph, ie: subalignment"""
    
    def __init__(self, leaf, pog, children=None):
        assert len(leaf) == len(pog), (len(leaf), len(pog))
        _Alignable.__init__(self, leaf)
        self.pred = pog.asListOfPredLists()
        self.max_preds = max(len(pre) for pre in self.pred)
        self.pog = pog
        if children is not None:
            self.aligneds = self._calcAligneds(children)
        self.leaf = leaf
    
    def __repr__(self):
        return 'AlPOG(%s,%s)' % (self.pog.all_jumps, repr(self.leaf))
    
    def getAlignment(self):
        return LoadSeqs(data=self.aligneds)
    
    def _calcAligneds(self, children):
        word_length = self.alphabet.getMotifLen()
        (starts, ends, maps) = map_traceback(self.pog.getFullAlignedPositions())
        aligneds = []
        for (dim, child) in enumerate(children):
            for (seq_name, aligned) in child.aligneds:
                #aligned = aligned[(starts[dim]-1)*word_length:(ends[dim]-1)*word_length]
                aligned = aligned.remappedTo((maps[dim]*word_length).inverse())
                aligneds.append((seq_name, aligned))
        return aligneds
    
    def backward(self):
        return self.__class__(self.leaf.backward(), self.pog.backward())
    
    def getPOG(self):
        return self.pog
    
    def __len__(self):
        return len(self.pred)
    
    def __iter__(self):
        return iter(self.pred)
    
    def __getitem__(self, index):
        # XXX the int case should be a different method?
        if isinstance(index, int):
            return self.pred[index]
        else:
            pog = self.pog[index]
            leaf = self.leaf[index]
            return AlignablePOG(leaf, pog)
    
    def midlinks(self):
        return self.pog.midlinks()
    
    def getOuterLoopDiscardPoints(self):
        # for score row caching
        last_successor = {}
        discard_list = {}
        for (successor, ps) in enumerate(self):
            for i in ps:
                last_successor[i] = successor
            discard_list[successor] = []
        for (i, successor) in last_successor.items():
            discard_list[successor].append(i)
        return discard_list
    

class AlignableSeq(_Alignable):
    """Wrapper for a Sequence which gives it the same interface as an
    AlignablePOG"""
    
    def __init__(self, leaf):
        _Alignable.__init__(self, leaf)
        if hasattr(leaf, 'sequence'):
            self.seq = leaf.sequence
            aligned = Aligned([(0, len(self.seq))], self.seq, len(self.seq))
            self.aligneds = [(self.leaf.edge_name, aligned)]
        self.max_preds = 1
        self._pog = None
    
    def __repr__(self):
        return 'AlSeq(%s)' % (getattr(self, 'seq', '?'))
    
    def getPOG(self):
        if self._pog is None:
            self._pog = leaf2pog(self.leaf)
        return self._pog
    
    def __len__(self):
        return len(self.index)
    
    def backward(self):
        return self.__class__(self.leaf.backward())
    
    def __iter__(self):
        # empty list 1st since 0th position has no predecessor
        yield []
        for i in range(1, len(self.index)):
            yield [i-1]
    
    def __getitem__(self, index):
        # XXX the int case should be a different method?
        if isinstance(index, int):
            if index == 0:
                return []
            elif 0 < index < len(self.index):
                return [index-1]
            else:
                raise IndexError(index)
        #elif index == slice(None, None, None):
        #    return self
        else:
            return AlignableSeq(self.leaf[index])

    def midlinks(self):
        half = len(self.leaf) // 2
        return  [(half, half)]

    def getOuterLoopDiscardPoints(self):
        return [[]] + [[i] for i in range(len(self)-1)]
    

def adaptPairTM(pairTM, finite=False):
    # constructs state_directions
    if finite:
        # BEGIN and END already specified
        assert list(pairTM.Tags[0]) == list(pairTM.Tags[-1]) == []
        T = pairTM.Matrix
        assert not T[-1, ...] and not T[..., 0]
        for tag in pairTM.Tags[1:-1]:
            assert tag, 'silent state'
        state_directions_list = list(enumerate(pairTM.Tags[1:-1]))
    else:
        pairTM = pairTM.withoutSilentStates()
        stationary_probs = numpy.array(pairTM.StationaryProbs)
        T = pairTM.Matrix
        full_matrix = numpy.zeros([len(T)+2, len(T)+2], float)
        full_matrix[1:-1,1:-1] = T
        full_matrix[0,1:-1] = stationary_probs # from BEGIN
        full_matrix[:,-1] = 1.0  #  to END
        T = full_matrix
        state_directions_list = list(enumerate(pairTM.Tags))
    
    this_row_last = lambda (state, (dx,dy)):(not (dx or dy), not dx)
    state_directions_list.sort(key=this_row_last)
    # sorting into desirable order (sort may not be necessary)
    
    state_directions = numpy.zeros([len(state_directions_list), 4], int)
    for (i, (state, emit)) in enumerate(state_directions_list):
        (dx, dy) = emit
        assert dx==0 or dy==0 or dx==dy
        bin = max(dx, dy)-1
        state_directions[i] = (state+1, bin, dx>0, dy>0)
    return (state_directions, T)


class PairEmissionProbs(object):
    """A pair of sequences and the psubs that relate them, but no gap TM"""
    def __init__(self, pair, bins):
        self.pair = pair
        self.bins = bins
        self.scores = {}
    
    def makePartialLikelihoods(self, use_cost_function):
        # use_cost_function specifies whether eqn 2 of Loytynoja & Goldman 
        # is applied.  Without it insertions may be favored over deletions
        # because the emission probs of the insert aren't counted.
        plhs = [[], []]
        gap_plhs = [[], []]
        for bin in self.bins:
            for (dim, pred) in enumerate(self.pair.children):
                # first and last should be special START and END nodes
                plh = numpy.inner(pred.plh, bin.ppsubs[dim])
                gap_plh = numpy.inner(pred.plh, bin.mprobs)
                if use_cost_function:
                    plh /= gap_plh[..., numpy.newaxis]
                    gap_plh[:] = 1.0
                else:
                    gap_plh[0] = gap_plh[-1] = 1.0
                gap_plhs[dim].append(gap_plh)
                plhs[dim].append(plh)
        for dim in [0,1]:
            plhs[dim] = numpy.array(plhs[dim])
            gap_plhs[dim] = numpy.array(gap_plhs[dim])
        return (plhs, gap_plhs)
    
    def _makeEmissionProbs(self, use_cost_function):
        (plhs, gap_scores) = self.makePartialLikelihoods(use_cost_function)
        match_scores = numpy.zeros([len(self.bins)] + self.pair.uniq_size,
                                    float)
        for (b, (x, y, bin)) in enumerate(zip(plhs[0], plhs[1], self.bins)):
            match_scores[b] = numpy.inner(x*bin.mprobs, y)
        match_scores[:, 0, 0] = match_scores[:, -1, -1] = 1.0
        return (match_scores, gap_scores)
    
    def _getEmissionProbs(self, use_logs, use_cost_function):
        key = (use_logs, use_cost_function)
        if key not in self.scores:
            if use_logs:
                (M, (X, Y)) = self._getEmissionProbs(False, use_cost_function)
                (M, X, Y) = [numpy.log(a) for a in [M, X, Y]]
                self.scores[key] = (M, (X, Y))
            else:
                self.scores[key] = self._makeEmissionProbs(use_cost_function)
        return self.scores[key]
    
    def _calc_global_probs(self, pair, scores, kw, state_directions,
            T, rows, cells, backward=False):
        if kw['use_logs']:
            (impossible, inevitable) = (-numpy.inf, 0.0)
        else:
            (impossible, inevitable) = (0.0, 1.0)
        (M, N) = pair.size
        (mantissas, exponents) = rows
        mantissas[0,0,0] = inevitable
        if exponents is not None:
            exponents[0,0,0] = 0
        probs = []
        last_i = -1
        to_end = numpy.array([(len(T)-1, 0, 0, 0)])
        for (state, (i,j)) in cells:
            if i > last_i:
                rr = pair.calcRows(last_i+1, i+1, 0, N-1,
                    state_directions, T, scores, rows, None, None, **kw)
            else:
                assert i == last_i, (i, last_i)
            last_i = i
            T2 = T.copy()
            if backward:
                T2[:, -1] = T[:, state]
            else:
                T2[:, -1] = impossible
                T2[state, -1] = inevitable
            global DEBUG
            _d = DEBUG
            DEBUG = False
            (maxpos, state, score) = pair.calcRows(
                    i, i+1, j, j+1, to_end, T2, scores, rows, None, None, **kw)
            DEBUG = _d
            probs.append(score)
        return numpy.array(probs)
        
    def __getitem__(self, index):
        assert len(index) == 2, index
        return PairEmissionProbs(self.pair[index], self.bins)
    
    def hirschberg(self, TM, dp_options):
        """linear-space alignment algorithm
        A linear space algorithm for computing maximal common subsequences.
        Comm. ACM 18,6 (1975) 341-343.
        Dan Hirschberg
        """
        (states, T) = TM
        
        # This implementation is slightly complicated by the need to handle 
        # alignments of alignments, because a subalignment may have an indel
        # spanning the midpoint where we want to divide the problem in half.
        # That must be the sense in which the fatter and slower method used
        # in "Prank" (Loytynoja A, Goldman N. 2005) is "computationally more 
        # attractive": for them there is only one link in the list:
        links = self.pair.children[0].midlinks()

        def _half_row_scores(backward):
            T2 = T.copy()
            if backward:
                T2[0, 1:-1] = 1.0  # don't count the begin state transition twice
            else:
                T2[1:-1:,-1] = 1.0  # don't count the end state transition twice
            return self.scores_at_rows(
                (states, T2), dp_options, 
                last_row=[link[backward] for link in links],
                backward = not not backward)
        
        (last_row1, last_row2) = parallel.map(_half_row_scores, [0,1])
        middle_row = (last_row1 + last_row2)
        (link, anchor, anchor_state) = numpy.unravel_index(
            numpy.argmax(middle_row.flat), middle_row.shape)
        score = middle_row[link, anchor, anchor_state]
        (join1, join2) = links[link]
        
        def _half_solution(part):
            T2 = T.copy()
            if part == 0:
                T2[-1] = 0.0
                T2[anchor_state, -1] = 1.0  # Must end with the anchor's state
                part = self[:join1, :anchor]
            else:
                T2[0, :] = T[anchor_state, :]  # Starting from the anchor's state
                part = self[join2:, anchor:]
            return part.dp((states, T2), dp_options)
        
        [(s1, tb_a), (s2, tb_b)] = parallel.map(_half_solution, [0,1])
        tb = tb_a + tb_b.offset(join2, anchor)
        # Same return as for self.dp(..., tb=...)
        return score, tb
    
    def scores_at_rows(self, TM, dp_options, last_row, backward=False):
        """A score array shaped [rows, columns, states] but only for those
        row numbers requested.  Used by Hirschberg algorithm"""
        (M, N) = self.pair.size
        (state_directions, T) = TM
        reverse = bool(dp_options.backward) ^ bool(backward)
        cells = []
        p_rows = sorted(set(last_row))
        if reverse:
            p_rows.reverse()
        for i in p_rows:
            for j in range(0, N-1):
                for state in range(len(T)-1):
                    if reverse:
                        cells.append((state, (M-2-i, N-2-j)))
                    else:
                        cells.append((state, (i, j)))
        probs = self.dp(TM, dp_options, cells=cells, backward=backward)
        probs = numpy.array(probs)
        probs.shape = (len(p_rows), N-1, len(T)-1)
        result = numpy.array([
            probs[p_rows.index(i)] for i in last_row])
        return result
        
    def dp(self, TM, dp_options, cells=None, backward=False):
        """Score etc. from a Dynamic Programming function applied to this pair.
        
        TM - (state_directions, array) describing the Transition Matrix.
        dp_options - instance of DPFlags indicating algorithm etc.
        cells - List of (state, posn) for which posterior probs are requested.
        backward - run algorithm in reverse order.
        """
        (state_directions, T) = TM
        if dp_options.viterbi and cells is None:
            encoder = self.pair.getPointerEncoding(len(T))
            problem_dimensions = self.pair.size + [len(T)]
            problem_size = numpy.product(problem_dimensions)
            memory = problem_size * encoder.bytes / 10**6
            if dp_options.local:
                msg = 'Local alignment'
            elif cells is not None:
                msg = 'Posterior probs'
            elif self.pair.size[0]-2 >= 3 and not backward and (
                    problem_size > HIRSCHBERG_LIMIT or 
                    parallel.getCommunicator().Get_size() > 1):
                 return self.hirschberg(TM, dp_options)
            else:
                msg = 'dp'
            if memory > 500:
                warnings.warn('%s will use > %sMb.' % (msg, memory))
            track = encoder.getEmptyArray(problem_dimensions)
        else:
            track = encoder = None
        
        kw = dict(
                use_scaling=dp_options.use_scaling,
                use_logs=dp_options.use_logs,
                viterbi=dp_options.viterbi,
                local=dp_options.local)
        
        if dp_options.backward:
            backward = not backward
        
        if backward:
            pair = self.pair.backward()
            origT = T
            T = numpy.zeros(T.shape, float)
            T[1:-1,1:-1] = numpy.transpose(origT[1:-1,1:-1])
            T[0,:] = origT[:, -1]
            T[:,-1] = origT[0,:]
        else:
            pair = self.pair
        
        if dp_options.use_logs:
            T = numpy.log(T)
        
        scores = self._getEmissionProbs(
                dp_options.use_logs, dp_options.use_cost_function)
        
        rows = pair.getEmptyScoreArrays(len(T), dp_options)
        
        if cells is not None:
            assert not dp_options.local
            result = self._calc_global_probs(
                    pair, scores, kw, state_directions, T, rows, cells,
                    backward)
        else:
            (M, N) = pair.size
            if dp_options.local:
                (maxpos, state, score) =  pair.calcRows(1, M-1, 1, N-1,
                    state_directions, T, scores, rows, track, encoder, **kw)
            else:
                pair.calcRows(0, M-1, 0, N-1,
                    state_directions, T, scores, rows, track, encoder, **kw)
                end_state_only = numpy.array([(len(T)-1, 0, 1, 1)])
                (maxpos, state, score) = pair.calcRows(M-1, M, N-1, N,
                    end_state_only, T, scores, rows, track, encoder, **kw)
                    
            if track is None:
                result = score
            else:
                tb = self.pair.traceback(track, encoder, maxpos, state,
                    skip_last = not dp_options.local)
                result = (score, tb)
        return result
        
    def getAlignable(self, aligned_positions, ratio=None):
        assert ratio is None, "already 2-branched"
        children = self.pair.children # alignables
        leaves = [c.leaf for c in children]
        aligned_positions = [posn for (bin, posn) in aligned_positions]
        pog = self.pair.getPOG(aligned_positions)
        edge = LikelihoodTreeEdge(leaves, 'parent', pog.getAlignedPositions())
        (plhs, gapscores) = self.makePartialLikelihoods(use_cost_function=False)
        plh = self.pair.edge2plh(edge, plhs)
        assert len(plh) == 1, ('bins!', len(plh))
        leaf = edge.asLeaf(plh[0]) # like profile
        return AlignablePOG(leaf, pog, children)
    
    def makePairHMM(self, transition_matrix, finite=False):
        # whether TM includes Begin and End states
        return PairHMM(self, transition_matrix, finite=finite)
    

class BinData(object):
    def __init__(self, mprobs, Qd, rate=1.0):
        self.Qd = Qd
        self.mprobs = mprobs
        self.rate = rate
    
    def forLengths(self, length1, length2):
        psub1 = self.Qd(length1 * self.rate)
        psub2 = self.Qd(length2 * self.rate)
        return PairBinData(self.mprobs, psub1, psub2)
    

class PairBinData(object):
    def __init__(self, mprobs, psub1, psub2):
        self.mprobs = mprobs
        self.ppsubs = [psub1, psub2]
    

class ReversiblePairEmissionProbs(object):
    """A pair of sequences and the psubs that relate them, but no gap TM
    
    'Reversible' in the sense that how `length` is divided between the 2 edges
    shouldn't change the forward and viterbi results"""
    
    def __init__(self, pair, bins, length):
        self.pair = pair
        self.bins = bins
        self.length = length
        self.midpoint = self._makePairEmissionProbs(0.5)
    
    def dp(self, *args, **kw):
        return self.midpoint.dp(*args, **kw)
    
    def _makePairEmissionProbs(self, ratio):
        assert 0.0 <= ratio <= 1.0
        lengths = [self.length * ratio, self.length * (1.0-ratio)]
        pbins = [bin.forLengths(*lengths) for bin in self.bins]
        return PairEmissionProbs(self.pair, pbins)
    
    def getAlignable(self, a_p, ratio=None):
        # a_p alignment positions
        if ratio in [None, 0.5]:
            ep = self.midpoint
        else:
            ep = self._makePairEmissionProbs(ratio=ratio)
        return ep.getAlignable(a_p)
    
    def makePairHMM(self, transition_matrix):
        return PairHMM(self, transition_matrix)
    

class DPFlags(object):
    def __init__(self, viterbi, local=False, use_logs=None,
            use_cost_function=True, use_scaling=None, backward=False):
        if use_logs is None:
            use_logs = viterbi and not use_scaling
        if use_scaling is None:
            use_scaling = not use_logs
        if use_logs:
            assert viterbi and not use_scaling
        self.use_cost_function = use_cost_function
        self.local = local
        self.use_logs = use_logs
        self.use_scaling = use_scaling
        self.viterbi = viterbi
        self.backward = backward
        self.as_tuple = (local, use_logs, use_cost_function, use_scaling,
                viterbi, backward)
    
    def __hash__(self):
        return hash(self.as_tuple)
    
    def __eq__(self, other):
        return self.as_tuple == other.as_tuple
    

class PairHMM(object):
    def __init__(self, emission_probs, transition_matrix, finite=False):
        self.emission_probs = emission_probs
        self.transition_matrix = transition_matrix
        self._transition_matrix = adaptPairTM(transition_matrix,
                finite=finite)
        self.results = {}
    
    def _getDPResult(self, **kw):
        dp_options = DPFlags(**kw)
        if dp_options not in self.results:
            self.results[dp_options] = \
                self.emission_probs.dp(self._transition_matrix, dp_options)
        return self.results[dp_options]
    
    def getForwardScore(self, **kw):
        return self._getDPResult(viterbi=False, **kw)
    
    def _getPosteriorProbs(self, tb, **kw):
        cells = tb.asStatePosTuples()
        score = self.getForwardScore(**kw)
        dp_options = DPFlags(viterbi=False, **kw)
        fwd = self.emission_probs.dp(self._transition_matrix, dp_options, cells)
        (N, M) = self.emission_probs.pair.size
        cells = [(state, (N-x-2, M-y-2)) for (state, (x,y)) in cells]
        tb.reverse()
        bck = self.emission_probs.dp(self._transition_matrix, dp_options, cells,
                backward=True)[::-1]
        return fwd + bck - score
    
    def getViterbiPath(self, **kw):
        result = self._getDPResult(viterbi=True,**kw)
        return ViterbiPath(self, result, **kw)
    
    def getViterbiScoreAndAlignment(self, ratio=None, **kw):
        # deprecate
        vpath = self.getViterbiPath(**kw)
        return (vpath.getScore(), vpath.getAlignment(ratio=ratio))
    
    def getLocalViterbiScoreAndAlignment(self, posterior_probs=False, **kw):
        # Only for pairwise.  Merge with getViterbiScoreAndAlignable above.
        # Local and POGs doesn't mix well.
        (vscore, tb) = self._getDPResult(viterbi=True, local=True, **kw)
        (state_directions, T) = self._transition_matrix
        aligned_positions = tb.asBinPosTuples(state_directions)
        seqs = self.emission_probs.pair.getSeqNamePairs()
        aligned_positions = [posn for (bin, posn) in aligned_positions]
        word_length = self.emission_probs.pair.alphabet.getMotifLen()
        align = alignment_traceback(seqs, aligned_positions, word_length)
        if posterior_probs:
            pp = self._getPosteriorProbs(tb, use_cost_function=False)
            return (vscore, align, numpy.exp(pp))
        else:
            return (vscore, align)
    

class ViterbiPath(object):
    def __init__(self, pair_hmm, result, **kw):
        (self.vscore, self.tb) = result
        (state_directions, T) = pair_hmm._transition_matrix
        self.aligned_positions = self.tb.asBinPosTuples(state_directions)
        self.pair_hmm = pair_hmm
        self.kw = kw
    
    def getScore(self):
        return self.vscore
    
    def getAlignable(self, ratio=None):
        # Because the alignment depends on the total length (so long as the
        # model is reversable!) the same cached viterbi result can be re-used
        # to calculate the partial likelihoods even if the root of the 2-seq
        # tree is moved around.
        alignable = self.pair_hmm.emission_probs.getAlignable(
                self.aligned_positions, ratio=ratio)
        return alignable
    
    def getAlignment(self, **kw):
        alignable = self.getAlignable(**kw)
        return alignable.getAlignment()
    
    def getPosteriorProbs(self):
        pp = self.pair_hmm._getPosteriorProbs(self.tb, use_cost_function=True)
        return numpy.exp(pp)
    
