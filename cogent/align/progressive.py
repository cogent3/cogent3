#!/usr/bin/env python

from cogent import LoadTree
from cogent.phylo import nj as NJ
from cogent.phylo.distance import EstimateDistances
from cogent.core.info import Info

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def TreeAlign(model, seqs, tree=None, indel_rate=0.01, indel_length=0.01,
    show_progress = True, ests_from_pairwise=True, param_vals=None):
    """Returns a multiple alignment and tree.
    
    Uses the provided substitution model and a tree for determining the
    progressive order. If a tree is not provided a Neighbour Joining tree is
    constructed from pairwise distances estimated from pairwise aligning the
    sequences. If running in parallel, only the distance estimation is
    parallelised and only the master CPU returns the alignment and tree, other
    CPU's return None, None.
    
    Arguments:
        - model: a substitution model
        - seqs: a sequence collection
        - indel_rate, indel_length: parameters for the progressive pair-HMM
        - ests_from_pairwise: if no tree provided and True, the median value
          of the substitution model parameters are used
        - param_vals: named key, value pairs for model parameters. These
          override ests_from_pairwise.
    """
    _exclude_params = ['mprobs', 'rate', 'bin_switch']
    if param_vals:
        param_vals = dict(param_vals)
    else:
        param_vals = {}
    if isinstance(seqs, dict):
        seq_names = seqs.keys()
    else:
        seq_names = seqs.getSeqNames()
    
    two_seqs = len(seq_names) == 2
    
    if tree:
        tip_names = tree.getTipNames()
        tip_names.sort()
        seq_names.sort()
        assert tip_names == seq_names, \
            "names don't match between seqs and tree: tree=%s; seqs=%s" % \
            (tip_names, seq_names)
        ests_from_pairwise = False
    elif two_seqs:
        tree = LoadTree(tip_names=seqs.getSeqNames())
        ests_from_pairwise = False
    else:
        if show_progress:
            print "\nEstimating pairwise distances to build NJ tree"
        if ests_from_pairwise:
            est_params = [param for param in model.getParamList() \
                                    if param not in _exclude_params]
        else:
            est_params = None
        dcalc = EstimateDistances(seqs, model, do_pair_align=True,
                                    est_params=est_params)
        dcalc.run(show_progress=show_progress)
        dists = dcalc.getPairwiseDistances()
        tree = NJ.nj(dists)
    
    LF = model.makeLikelihoodFunction(tree.bifurcating(), aligned=False)
    if ests_from_pairwise and not param_vals:
        # we use the Median to avoid the influence of outlier pairs
        param_vals = {}
        for param in est_params:
            numbers = dcalc.getParamValues(param)
            if show_progress:
                print "\nParam Estimate Summary Stats: %s" % param
                print numbers.summarize()
            param_vals[param] = numbers.Median
    
    for param, val in param_vals.items():
        LF.setParamRule(param, value=val, is_const=True)
    
    if show_progress:
        msg = "\nDoing %s alignment" % ["progressive", "pairwise"][two_seqs]
        print msg
    
    LF.setParamRule('indel_rate', value=indel_rate, is_const=True)
    LF.setParamRule('indel_length', value=indel_length, is_const=True)
    
    LF.setSequences(seqs)
    edge = LF.getLogLikelihood().edge
    (vtLnL, align) = edge.getViterbiScoreAndAlignment(0.5)
    info = Info()
    info["AlignParams"] = param_vals
    info["AlignParams"].update(dict(indel_length=indel_length, indel_rate=indel_rate))
    align.Info = info
    return align, tree
