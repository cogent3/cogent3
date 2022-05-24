#!/usr/bin/env python


import numpy

from cogent3.maths.markov import TransitionMatrix


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def pair_transition_matrix(order, a):
    """A matrix for Pair HMMs with gap states X and Y, match state M,
    and optionally a silent wait state W"""
    size = len(order)
    assert a.shape == (size, size)
    directions = {"W": (0, 0), "X": (1, 0), "Y": (0, 1), "M": (1, 1)}
    emits = [directions[state.upper()] for state in order]
    return TransitionMatrix(a, emits).withoutSilentStates()


def classic_gap_scores(d, e):
    """gap open / gap extend costs.  No X to Y transitions"""
    _ = numpy.inf
    C = numpy.array([[e, _, 0], [_, e, 0], [d, d, 0]])
    T = numpy.exp(-1.0 * C)
    T = T / numpy.sum(T, axis=1)[..., numpy.newaxis]
    return pair_transition_matrix("XYM", T)


class _SimpleIndelParams(object):
    def __init__(self, indel_rate, indel_length):
        assert 0.0 < indel_length < 1.0, indel_length
        assert 0.0 < indel_rate < 1.0, indel_rate
        self.indel_rate = indel_rate
        self.indel_length = indel_length


class SimpleIndelModel(_SimpleIndelParams):
    """P(gap open), P(gap extend) with const P(extend match)"""

    def calc_transition_matrix(self, distance):
        d = 1.0 - numpy.exp(-self.indel_rate * distance)
        e = self.indel_length
        g = 0.0
        T = numpy.array(
            [[0, d, d, 1 - 2 * d], [1 - e, e, 0, 0], [1 - e, 0, e, 0], [1 - g, 0, 0, g]]
        )
        return pair_transition_matrix("WXYM", T).withoutSilentStates()


class KnudsenMiyamotoIndelModel(_SimpleIndelParams):
    """Sequence Alignments and Pair HMMs Using Evolutionary History
    Bjarne Knudsen and Michael M. Miyamoto
    Journal of Molecular Biology
    Volume 333, Issue 2 , 17 October 2003, Pages 453-460"""

    def calc_transition_matrix(self, distance):
        distance = distance * self.indel_rate
        extend = self.indel_length
        close = 1.0 - extend

        # First indel event
        indel = 1.0 - numpy.exp(-distance)
        insert = deletion = indel / 2.0

        # Second indel event
        if distance < 0.0001:
            secondary_indel = distance / 2
        else:
            x = numpy.exp(-distance - numpy.log(distance))
            secondary_indel = 1 - 1 / distance + x
        assert 0.0 <= secondary_indel <= 1.0
        secondary_deletion = secondary_insert = secondary_indel / 2.0
        same_len = (1 - extend) / (1 + extend)
        longer = shorted = (1.0 - same_len) / 2

        eMX = insert * (1 - secondary_deletion * (1 - longer))
        eMY = deletion + insert * secondary_deletion * longer
        eMZ = eMX + eMY

        eXM = extend * secondary_deletion * same_len + close * (
            deletion / 4 + (1.0 - indel)
        )
        eXX = (
            extend * (secondary_insert * extend / close + secondary_deletion * shorted)
            + close * insert
            + extend
        )
        # eXX = a + a**2/(1-a**2)*secondary_indel + (1-a)*insert
        eXY = extend * secondary_deletion * longer + close * deletion * 3 / 4

        # e = 1 + ( extend * secondary_indel/2 / close)
        # print e, (eXM + eXX + eXY)
        e = eXM + eXX + eXY

        # assert eMM + eMX + eMY == 1.0, (eMM + eMX + eMY)

        tMX = tMY = eMZ / 2
        tMM = 1.0 - eMZ
        tXM = tYM = eXM / e
        tXX = tYY = eXX / e
        tXY = tYX = eXY / e

        tm = numpy.array([[tMM, tMX, tMY], [tXM, tXX, tXY], [tYM, tYX, tYY]])

        if min(tm.flat) < 0:
            raise ValueError

        return pair_transition_matrix("MXY", tm)
