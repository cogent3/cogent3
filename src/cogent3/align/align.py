import numpy

from cogent3.align import indel_model, pairwise
from cogent3.evolve.likelihood_tree import make_likelihood_tree_leaf


def make_dna_scoring_dict(
    match: float, transition: float, transversion: float
) -> dict[tuple[str, str], float]:
    score_mat = {}
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
            score_mat[a, b] = score
    return score_mat


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
            score = match if a == b else -1
            S[a, b] = score
    return S


def _align_pairwise(
    s1,
    s2,
    mprobs,
    psub,
    TM,
    local,
    return_alignment=True,
    return_score=False,
    **kw,
):
    """Generic alignment with any substitution model and indel model"""
    [p1, p2] = [
        make_likelihood_tree_leaf(seq, seq.moltype.alphabet, seq.name)
        for seq in [s1, s2]
    ]
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
        return alignment
    return score


def classic_align_pairwise(s1, s2, Sd, d, e, local, return_score=False, **kw):
    """Alignment specified by gap costs and a score matrix"""
    TM = indel_model.classic_gap_scores(d, e)
    a1 = s1.moltype.alphabet
    a2 = s2.moltype.alphabet
    S = numpy.zeros([len(a1), len(a2)], float)
    for i, m1 in enumerate(a1):
        for j, m2 in enumerate(a2):
            S[i, j] = Sd[m1, m2]
    psub = numpy.exp(S)
    mprobs = numpy.ones(len(psub), float) / len(psub)
    return _align_pairwise(
        s1,
        s2,
        mprobs,
        psub,
        TM,
        local,
        return_score=return_score,
        **kw,
    )


# these can't do codon sequences
# they could be replaced with something more sophisticated, like the HMM


def local_pairwise(s1, s2, S, d, e, return_score=False):
    return classic_align_pairwise(s1, s2, S, d, e, True, return_score=return_score)


def global_pairwise(s1, s2, S, d, e, return_score=False):
    return classic_align_pairwise(s1, s2, S, d, e, False, return_score=return_score)
