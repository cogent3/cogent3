#!/usr/bin/env python

import numpy

from cogent3.align import indel_model, pairwise, pycompare
from cogent3.evolve.likelihood_tree import make_likelihood_tree_leaf


Float = numpy.core.numerictypes.sctype2char(float)


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def dotplot(seq1, seq2, window, threshold, min_gap_length=0, band=None, **kw):
    # warnings.warn("cogent3.align.align.dotplot moved to cogent3.align.compare.dotplot",
    #    DeprecationWarning)
    return pycompare.dotplot(seq1, seq2, window, threshold, min_gap_length, band, **kw)


def make_dna_scoring_dict(match, transition, transversion):
    DNA = {}
    for a in "ATCG":
        ar = a in "AG"
        for b in "ATCG":
            br = b in "AG"
            if a == b:
                score = match
            elif ar == br:
                score = transition
            else:
                score = transversion
            DNA[a, b] = score
    return DNA


def make_generic_scoring_dict(match, mtype):
    """returns scoring dict for alignment

    Parameters
    ----------
    match : int
        value for a match, mismatches default to -1
    mtype
        MolType instance or string that can be used to get_moltype
    """
    from cogent3 import get_moltype

    mtype = get_moltype(mtype)
    S = {}
    for a in mtype:
        for b in mtype:
            if a == b:
                score = match
            else:
                score = -1
            S[a, b] = score
    return S


def _align_pairwise(
    s1, s2, mprobs, psub, TM, local, return_alignment=True, return_score=False, **kw
):
    """Generic alignment with any substitution model and indel model"""
    [p1, p2] = [make_likelihood_tree_leaf(seq) for seq in [s1, s2]]
    [p1, p2] = [pairwise.AlignableSeq(leaf) for leaf in [p1, p2]]
    pair = pairwise.Pair(p1, p2)
    EP = pair.make_simple_emission_probs(mprobs, [psub])
    hmm = EP.make_pair_HMM(TM)
    vpath = hmm.get_viterbi_path(local=local, **kw)
    score = vpath.get_score()
    if return_alignment:
        alignment = vpath.get_alignment()
        if return_score:
            return (alignment, score)
        else:
            return alignment
    else:
        return score


def classic_align_pairwise(s1, s2, Sd, d, e, local, return_score=False, **kw):
    """Alignment specified by gap costs and a score matrix"""
    TM = indel_model.classic_gap_scores(d, e)
    a1 = s1.moltype.alphabet
    a2 = s2.moltype.alphabet
    S = numpy.zeros([len(a1), len(a2)], Float)
    for (i, m1) in enumerate(a1):
        for (j, m2) in enumerate(a2):
            S[i, j] = Sd[m1, m2]
    psub = numpy.exp(S)
    mprobs = numpy.ones(len(psub), Float) / len(psub)
    return _align_pairwise(
        s1, s2, mprobs, psub, TM, local, return_score=return_score, **kw
    )


# these can't do codon sequences
# they could be replaced with something more sophisticated, like the HMM


def local_pairwise(s1, s2, S, d, e, return_score=False):
    return classic_align_pairwise(s1, s2, S, d, e, True, return_score=return_score)


def global_pairwise(s1, s2, S, d, e, return_score=False):
    return classic_align_pairwise(s1, s2, S, d, e, False, return_score=return_score)
