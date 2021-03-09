import math

from cogent3 import DNA, make_aligned_seqs


def logLikelihood(pi, rateMatrix, alignment):
    """Produced a log-likelihood of an unrooted tree of three DNA sequences

    Given a constant rate matrix and the ancestors ratios calculate the
    log-likelihood of three sequences unrooted tree.

    Paramters
    ---------
    pi: Sequence
        An array of length 4, with ratios of ancestor in order T,C,G,A

    rateMatrix: Sequence of Sequence
          A 4X4 sequence representing each possible transition and their probabilities
          ordered T,C,G,A where (2,4) entry represents P(C -> A)

    Returns
    -------
    log(probability): The log of the likelihood of the unrooted tree of the
                      three sequences

    Examples
    --------
    >>> likelihoodCalculator(pi, rateMatrix, alignment2) # See above
    -3.278666124
    """
    probability = 1
    length = len(alignment)
    for i in range(0, length):
        probability *= probabilityPosition(pi, rateMatrix, alignment[i])
    return math.log10(probability)


def probabilityPosition(pi, rm, alnCol):
    probability = 0
    for i in range(0, 4):
        probability += (
            pi[i]
            * rm[i][baseToPos(alnCol.seqs[0])]
            * rm[i][baseToPos(alnCol.seqs[1])]
            * rm[i][baseToPos(alnCol.seqs[2])]
        )
    return probability


def baseToPos(base):
    if base == "T":
        return 0
    if base == "C":
        return 1
    if base == "G":
        return 2
    if base == "A":
        return 3
