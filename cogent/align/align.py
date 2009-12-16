#!/usr/bin/env python

import numpy
Float = numpy.core.numerictypes.sctype2char(float)
import logging
LOG = logging.getLogger('cogent.align')

from cogent.align import pairwise, indel_model
from cogent.evolve.likelihood_tree import makeLikelihoodTreeLeaf

try:
    from _compare import dotplot
except ImportError:
    from pycompare import dotplot
    LOG.warning('Using slow Python dotplot implementation')

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def make_dna_scoring_dict(match, transition, transversion):
    DNA = {}
    for a in 'ATCG':
        ar = a in 'AG'
        for b in 'ATCG':
            br = b in 'AG'
            if a == b:
                score = match
            elif ar == br:
                score = transition
            else:
                score = transversion
            DNA[a,b] = score
    return DNA

def _align_pairwise(s1, s2, mprobs, psub, TM, local, **kw):
    """Generic alignment with any substitution model and indel model"""
    [p1, p2] = [makeLikelihoodTreeLeaf(seq) for seq in [s1, s2]]
    [p1, p2] = [pairwise.AlignableSeq(leaf) for leaf in [p1, p2]]
    pair = pairwise.Pair(p1, p2)
    EP = pair.makeSimpleEmissionProbs(mprobs, [psub])
    hmm = EP.makePairHMM(TM)
    if local:
        (score, alignment) = hmm.getLocalViterbiScoreAndAlignment(**kw)
    else:
        (score, alignment) = hmm.getViterbiScoreAndAlignment(**kw)
    return alignment

def classic_align_pairwise(s1, s2, Sd, d, e, local, **kw):
    """Alignment specified by gap costs and a score matrix"""
    TM = indel_model.ClassicGapScores(d, e)
    a1 = s1.MolType.Alphabet
    a2 = s2.MolType.Alphabet
    S = numpy.zeros([len(a1), len(a2)], Float)
    for (i,m1) in enumerate(a1):
        for (j,m2) in enumerate(a2):
            S[i, j] = Sd[m1, m2]
    psub = numpy.exp(S)
    mprobs = numpy.ones(len(psub), Float) / len(psub)
    return _align_pairwise(s1, s2, mprobs, psub, TM, local, **kw)

# these can't do codon sequences
# they could be replaced with something more sophisticated, like the HMM
# may not give same answer's as algorithm
def local_pairwise(s1, s2, S, d, e):
    return classic_align_pairwise(s1, s2, S, d, e, True)

def global_pairwise(s1, s2, S, d, e):
    return classic_align_pairwise(s1, s2, S, d, e, False)
