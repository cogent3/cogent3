#!/usr/bin/env python


from cogent3 import LoadTree
from cogent3.phylo import nj as NJ
from cogent3.phylo.distance import EstimateDistances
from cogent3.core.info import Info
from cogent3.util import progress_display as UI

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


@UI.display_wrap
def TreeAlign(model, seqs, tree=None, indel_rate=0.01, indel_length=0.01,
              ui=None, ests_from_pairwise=True, param_vals=None):
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
        seq_names = list(seqs.keys())
    else:
        seq_names = seqs.get_seq_names()

    two_seqs = len(seq_names) == 2

    if tree:
        tip_names = tree.get_tip_names()
        tip_names.sort()
        seq_names.sort()
        assert tip_names == seq_names, \
            "names don't match between seqs and tree: tree=%s; seqs=%s" % \
            (tip_names, seq_names)
        ests_from_pairwise = False
    elif two_seqs:
        tree = LoadTree(tip_names=seqs.get_seq_names())
        ests_from_pairwise = False
    else:
        if ests_from_pairwise:
            est_params = [param for param in model.get_param_list()
                          if param not in _exclude_params]
        else:
            est_params = None

        dcalc = EstimateDistances(seqs, model, do_pair_align=True,
                                  est_params=est_params)
        dcalc.run()
        dists = dcalc.get_pairwise_distances()
        tree = NJ.nj(dists)

    LF = model.make_likelihood_function(
        tree.bifurcating(name_unnamed=True), aligned=False)
    if ests_from_pairwise and not param_vals:
        # we use the Median to avoid the influence of outlier pairs
        param_vals = {}
        for param in est_params:
            numbers = dcalc.get_param_values(param)
            print("Param Estimate Summary Stats: %s" % param)
            print(numbers.summarize())
            param_vals[param] = numbers.Median

    ui.display("Doing %s alignment" % ["progressive", "pairwise"][two_seqs])
    with LF.updates_postponed():
        for param, val in list(param_vals.items()):
            LF.set_param_rule(param, value=val, is_constant=True)
        LF.set_param_rule('indel_rate', value=indel_rate, is_constant=True)
        LF.set_param_rule('indel_length', value=indel_length, is_constant=True)
        LF.set_sequences(seqs)
    edge = LF.get_log_likelihood().edge
    align = edge.getViterbiPath().getAlignment()
    info = Info()
    info["AlignParams"] = param_vals
    info["AlignParams"].update(
        dict(indel_length=indel_length, indel_rate=indel_rate))
    align.info = info
    return align, tree
