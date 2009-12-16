#!/usr/bin/env python

"""A data structure currently superseded by indel_positions.py.

Multiple Sequence Alignment Using Partial Order Graphs

Bioinformatics 2002 18:452-464
Christopher Lee, Catherine Grasso, & Mark Sharlow
"""

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def pog_traceback(pogs, aligned_positions):
    full_pogs = [getattr(pog, 'full', pog) for pog in pogs]
    upto = [0, 0]
    align_builder = POGBuilder(pogs)
    full_builder = POGBuilder(full_pogs)
    full_coords = []
    full_count = 0
    full_aligned_positions = []
    for posn in aligned_positions:
        # Find equivalent position in the full POG
        # which includes all of the original sites
        full_posn = []
        for (pog,pos) in zip(pogs, posn):
            if pos is not None:
                pos = pog.full_coords[pos]
            full_posn.append(pos)
        
        # Add all skipped over positions to the full POG
        for (dim, pos) in enumerate(full_posn):
            if pos is not None:
                for p in range(upto[dim]+1, pos+1):
                    full_builder.addAligned([(dim, p)])
                    full_count += 1
                    fp = [None, None]
                    fp[dim] = p - 1
                    full_aligned_positions.append(fp)
                upto[dim] = pos+1
        
        # Add this position to the full POG
        full_merged_origins = []
        fap = []
        for (dim, pos) in enumerate(full_posn):
            if pos is not None:
                full_merged_origins.append((dim,pos+1))
                fap.append(pos+1)
            else:
                fap.append(None)
        full_builder.addAligned(full_merged_origins)
        full_count += 1
        full_aligned_positions.append(fap)
        
        # Add this position to the minimal POG which does not
        # include sequence which has been called as an insert.
        merged_origins = []
        for (dim, pos) in enumerate(posn):
            if pos is not None:
                merged_origins.append((dim,pos+1))
        align_builder.addAligned(merged_origins)
        full_coords.append(full_count-1)
    
    result_pog = align_builder.getPOG()
    result_pog.aligned_positions = aligned_positions
    
    # Add leftover positions to the full POG
    for (dim, pos) in enumerate([len(p) for p in full_pogs]):
        if pos is not None:
            for p in range(upto[dim]+1, pos+1):
                full_builder.addAligned([(dim, p)])
                full_count += 1
                fp = [None, None]
                fp[dim] = p - 1
                full_aligned_positions.append(fp)
            upto[dim] = pos+1
    
    result_pog.full = full_builder.getPOG()
    result_pog.full_coords = full_coords
    
    result_pog.full.aligned_positions = full_aligned_positions
    
    return result_pog

class POGBuilder(object):
    def __init__(self, children):
        self.children = children
        self.remap = [{} for child in children]
        self.started = [False, False]
        self.last = [None, None]
        self.result = [[]]
        self.origins = [[]]
    
    def addAligned(self, nodes):
        pre_merged = set()
        origins_merged = []
        for (dim, pos) in nodes:
            if pos is None:
                continue
            if not self.started[dim]:
                pre_merged.add(0)
                self.started[dim] = True
            pog = self.children[dim]
            pre = pog.pred_sets[pos]
            origins_merged.append((dim,pos))
            pre = [self.remap[dim].get(n, None) for n in pre]
            pre_merged.update([p for p in pre if p is not None])
            self.remap[dim][pos] = len(self.result)
            self.last[dim] = pos
        self.result.append(pre_merged)
        self.origins.append(origins_merged)
    
    def getPOG(self):
        end = set(self.remap[dim][self.last[dim]] for dim in [0,1])
        self.result.append(end)
        self.origins.append([])
        return POG(self.result, self.origins)
    

class POG(object):
    def __init__(self, pog, origins):
        self.origins = origins
        if origins is not None:
            assert len(origins) == len(pog), (len(origins), len(pog))
        for (i, pred_set) in enumerate(pog):
            for pre in pred_set:
                assert 0 <= pre < i, pog
        self.pred_sets = pog
    
    def traceback(self, other, aligned_positions):
        return pog_traceback([self, other], aligned_positions)
    
    def asListOfPredLists(self):
        return [list(pre) for pre in self.pred_sets]
    
    def getAlignedPositions(self):
        return self.aligned_positions
    
    def getFullAlignedPositions(self):
        return self.full.aligned_positions
    
    def __len__(self):
        return len(self.pred_sets) - 2
    
    def __getitem__(self, index):
        if index.start is None:
            start = 1
        else:
            start = index.start + 1
        if index.stop is None:
            end = len(self.pred_sets)-2 + 1
        else:
            end = index.stop + 1
        pred = [list(pre) for pre in self.pred_sets[start:end]]
        origins = ['s'] + list(self.origins[start:end]) + ['e']
        before = self.pred_sets[1:start]
        after = self.pred_sets[end:]
        ends = []
        for (i,p) in enumerate(after):
            for pre in p:
                if pre < end:
                    ends.append(max(pre-start+1, 0))
        pog = [[]] + [[max(0, pre-start+1) for pre in p] for p in pred]
        pog.append(ends)
        return POG(pog, origins)
    
    def backward(self):
        # Switches predecessors / successors
        pog = [[]]
        for (i,p) in enumerate(self.pred_sets[1:]):
            pog.append([])
            for pre in p:
                pog[pre].append(len(self.pred_sets)-i-2)
        pog.reverse()
        return POG(pog, None)
    
    def writeToDot(self, dot):
        print >>dot, 'digraph POG {'
        for (i, preds) in enumerate(self.pred_sets):
            #print i, preds
            for pred in preds:
                print >>dot, '  ', ('node%s -> node%s' % (pred, i))
            if i == 0:
                label = 'START'
            elif i == len(self.pred_sets) - 1:
                label = 'END'
            else:
                label = ','.join(c for c in self.origins[i])
            print >>dot, '  ', ('node%s' % i), '[label="%s"]' % label
        print >>dot, '}'
        print >>dot, ''


class LeafPOG(POG):
    def __init__(self, length):
        pog = [set([i]) for i in range(length)]
        pog = [[]] + pog + [[len(pog)]]
        self.pred_sets = pog
        self.full = self
        self.full_coords = range(len(pog))
        self.length = length
        cap = [(0, None)]
        self.origins = [cap] + [[(0, i)]
                for i in range(length)] + [cap]
    
    def __len__(self):
        return len(self.pred_sets) - 2
    
    def backward(self):
        return LeafPOG(self.length)
    

def leaf2pog(leaf):
    return LeafPOG(len(leaf))
    