#!/usr/bin/env python'
from tree_space import TreeEvaluator, ancestry2tree
from least_squares import WLS
from math import exp
import consensus

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class ML(TreeEvaluator):
    """(err, best_tree) = ML(model, alignment).trex()
    
    If 'dists' is provided uses WLS to get initial values for lengths"""
    
    def __init__(self, model, alignment, dists=None, opt_args={}):
        self.opt_args = opt_args
        self.names = alignment.getSeqNames()
        self.alignment = alignment
        self.model = model
        if dists:
            self.wlsMakeTreeScorer = WLS(dists).makeTreeScorer
        else:
            fake_wls = lambda a:(None, None)
            self.wlsMakeTreeScorer = lambda n: fake_wls
    
    def evaluateTree(self, tree):
        names = tree.getTipNames()
        subalign = self.alignment.takeSeqs(names)
        lf = self.model.makeLikelihoodFunction(tree)
        lf.setAlignment(subalign)
        return lf.getLogLikelihood()
            
    def makeTreeScorer(self, names):
        subalign = self.alignment.takeSeqs(names)
        wls_eval = self.wlsMakeTreeScorer(names)
        def evaluate(ancestry, lengths=None):
            if lengths is None:
                (wls_err, init_lengths) = wls_eval(ancestry)
            else:
                init_lengths = lengths
            tree = ancestry2tree(ancestry, init_lengths, names)
            lf = self.model.makeLikelihoodFunction(tree)
            lf.setAlignment(subalign)
            if lengths is not None:
                lf.setParamRule('length', is_const=True)
            lf.optimise(show_progress=False, **self.opt_args)
            err = -1.0 * lf.getLogLikelihood()
            tree = lf.getAnnotatedTree()
            return (err, tree)
        return evaluate
    
    def result2output(self, err, ancestry, annotated_tree, names):
        return (-1.0*err, annotated_tree)

    def results2output(self, results):
        return LogLikelihoodScoredTreeCollection(results)
    

class _ScoredTreeCollection(list):
    def __getitem__(self, index):
        # Helpful to keep type after truncation like [self[:10]],
        # but not self[0] or self[x,y,-1]
        result = list.__getitem__(self, index)
        if isinstance(index, slice) and index.step is None:
            result = type(self)(result)
        return result

    def writeToFile(self, filename):
        f = open(filename, 'w')
        for (score, tree) in self:
            f.writelines(
                [str(score), '\t', tree.getNewick(with_distances=True),'\n'])
        f.close()
        
class WeightedTreeCollection(_ScoredTreeCollection):
    def getConsensusTree(self, strict=False):
        ctrees = consensus.weightedMajorityRule(self, strict)
        assert len(ctrees) == 1, len(ctrees)
        return ctrees[0]
        
class LogLikelihoodScoredTreeCollection(_ScoredTreeCollection):
    """An ordered list of scored trees all covering the same set of tip names"""
        
    def __init__(self, trees):
        list.__init__(self, trees)
        # Quick and very dirty check of order
        assert self[0][0] >= self[-1][0]
        
    def getConsensusTree(self, cutoff=None, strict=False):
        return self.getWeightedTrees(cutoff).getConsensusTree(strict)
        
    def getWeightedTrees(self, cutoff=None):
        if cutoff is None:
            cutoff = 0.99
        assert 0 <= cutoff <= 1.0
        max_lnL = self[0][0]
        weights = [exp(lnL-max_lnL) for (lnL, t) in self]
        # add from smallest end to avoid rounding errors
        weights.reverse()
        tail = (1.0-cutoff) * sum(weights)
        dropped = 0.0
        for (index, weight) in enumerate(weights):
            dropped += weight
            if dropped > tail:
                weights = weights[index:]
                break
        denominator = sum(weights)
        weights.reverse()
        return WeightedTreeCollection((weight/denominator, tree) 
                for (weight, (lnL, tree)) in zip(weights, self))
        

def LoadTrees(filename):
    """Parse a file of (score, tree) lines. Scores can be positive probabilities
    or negative log likelihoods."""
    from cogent import LoadTree
    infile = open(filename, 'r')
    trees = []
    klass = list
    for line in infile:
        line = line.split(None, 1)
        lnL = float(line[0])
        if lnL > 1:
            raise ValueError('likelihoods expected, not %s' % lnL)
        elif lnL > 0:
            assert klass in [list, WeightedTreeCollection]
            klass = WeightedTreeCollection
        else:
            assert klass in [list,  LogLikelihoodScoredTreeCollection]
            klass = LogLikelihoodScoredTreeCollection
        tree = LoadTree(treestring=line[1])
        trees.append((lnL, tree))
    trees.sort(reverse=True)
    return klass(trees)
    