#!/usr/bin/env python

import bisect

import numpy


Float = numpy.core.numerictypes.sctype2char(float)

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


class TransitionMatrix(object):
    """The transition matrix for a Markov process.  Just a square
    numpy array plus a list of state 'tags', eg:

    >>> a = numpy.array([
        [.9, .0, .1],
        [.0, .9, .1],
        [.5, .5, .0],], Float)
    >>> T = TransitionMatrix(a, ['x', 'y', 'm'])
    """

    def __init__(self, matrix, tags, stationary_probs=None):
        self.Matrix = numpy.array(matrix, Float)
        self.Tags = list(tags)
        self.size = len(matrix)
        assert matrix.shape == (self.size, self.size)
        assert len(tags) == self.size
        if stationary_probs is None:
            self._stationary_probs = None
        else:
            # Could recalculate, but if it is provided then faster to
            # just trust it
            assert len(stationary_probs) == self.size
            self._stationary_probs = numpy.array(stationary_probs, Float)

    def _getStationaryProbs(self):
        if self._stationary_probs is None:
            matrix = self.Matrix
            for i in range(10):
                matrix = numpy.core.multiarray.dot(matrix, matrix)
            self._stationary_probs = matrix[0]
        return self._stationary_probs

    StationaryProbs = property(_getStationaryProbs)

    def emit(self, random_series):
        """Generates an infinite sequence of states"""
        partitions = numpy.add.accumulate(self.Matrix, axis=1)
        for (state, row) in enumerate(partitions[:-1]):
            assert abs(row[-1] - 1.0) < 1e-6, (state, self.Matrix[state])
        x = random_series.uniform(0.0, 1.0)
        state = bisect.bisect_left(numpy.add.accumulate(self.StationaryProbs), x)
        while 1:
            yield self.Tags[state]
            x = random_series.uniform(0.0, 1.0)
            state = bisect.bisect_left(partitions[state], x)

    def __repr__(self):
        from cogent3.util.table import Table

        labels = []
        for (i, label) in enumerate(self.Tags):
            if hasattr(label, "__len__") and not isinstance(label, str):
                label = ",".join(str(z) for z in label)
            # Table needs unique labels
            label = f"{label} ({i})"
            labels.append(label)
        heading = [""] + labels
        a = [[name] + list(row) for (name, row) in zip(labels, self.Matrix)]
        return str(Table(header=heading, data=a))

    def withoutSilentStates(self):
        """An equivalent matrix without any of the states that have a
        false tag value"""
        N = self.size
        silent = numpy.array([(not max(tag)) for tag in self.Tags], Float)
        matrix = numpy.zeros([N, N], Float)
        for i in range(N):
            row = numpy.zeros([N], Float)
            emt = numpy.zeros([N], Float)
            row[i] = 1.0
            while max((row + emt) - emt):
                row = numpy.dot(row, self.Matrix)
                nul = silent * row
                emt += row - nul
                row = nul
            matrix[i] = emt
        keep = [i for i in range(self.size) if max(matrix[:, i])]
        matrix = numpy.take(matrix, keep, axis=0)
        matrix = numpy.take(matrix, keep, axis=1)
        tags = numpy.take(self.Tags, keep, axis=0)
        return type(self)(matrix, tags)

    def getLikelihoodOfSequence(self, obs, backward=False):
        """Just for testing really"""
        profile = numpy.zeros([len(obs), self.size], Float)
        for (i, a) in enumerate(obs):
            # This is suspiciously alphabet-like!
            profile[i, self.Tags.index(obs[i])] = 1.0
        return self.getLikelihoodOfProfile(profile, backward=backward)

    def getLikelihoodOfProfile(self, obs, backward=False):
        """Just for testing really"""
        if not backward:
            state_probs = self.StationaryProbs.copy()
            for i in range(len(obs)):
                state_probs = numpy.dot(state_probs, self.Matrix) * obs[i]
            return sum(state_probs)
        else:
            state_probs = numpy.ones([self.size], Float)
            for i in range(len(obs) - 1, -1, -1):
                state_probs = numpy.dot(self.Matrix, state_probs * obs[i])
            return sum(state_probs * self.StationaryProbs)

    def get_posterior_probs(self, obs):
        """'obs' is a sequence of state probability vectors"""
        result = numpy.zeros([len(obs), self.size], Float)

        # Forward
        state_probs = self.StationaryProbs.copy()
        for i in range(len(obs)):
            state_probs = numpy.dot(state_probs, self.Matrix) * obs[i]
            state_probs /= sum(state_probs)
            result[i] = state_probs

        # and Backward
        state_probs = numpy.ones([self.size], Float)
        for i in range(len(obs) - 1, -1, -1):
            state_probs = numpy.dot(self.Matrix, state_probs)
            state_probs /= sum(state_probs)
            result[i] *= state_probs
            result[i] /= sum(result[i])
            state_probs *= obs[i]

        return result

    def nestTransitionMatricies(self, Ts, blended=None):
        """Useful for combining several X/Y/M Pair HMMs into
        one large HMM. The transition matricies 'Ts' end up
        along the diagonal blocks of the result, and the off
        diagonal values depend on the region switching probablities
        defined by 'self'"""

        if blended is None:
            blended = lambda a, b: (a + b) / 2.0
            # blended = lambda a,b: numpy.sqrt(a*b)
            # blended = lambda a,b: b

        R = self.Matrix
        n = len(R)
        assert len(Ts) == n
        result = None
        for (x, a) in enumerate(self.Tags):
            a = Ts[a - 1]
            if result is None:
                tags = a.Tags
                c = len(tags)
                result = numpy.zeros([c * n, c * n], Float)
            else:
                assert a.Tags == tags, (a.Tags, tags)
            a = a.Matrix
            for (y, b) in enumerate(self.Tags):
                b = Ts[b - 1].Matrix
                block = self.Matrix[x, y] * blended(a, b)
                result[x * c : (x + 1) * c, y * c : (y + 1) * c] = block
        all_tags = []
        for i in self.Tags:
            all_tags.extend([[i * e for e in tag] for tag in tags])
        return TransitionMatrix(result, all_tags)


def SiteClassTransitionMatrix(switch, probs):
    """TM defined by stationary probabilities and a 'switch' parameter,
    switch=0.0 gives an identity Matrix, switch=1.0 gives independence,
    ie: zero-order markov process"""
    probs = numpy.asarray(probs)
    assert numpy.allclose(sum(probs), 1.0), probs
    I = numpy.identity(len(probs), Float)
    switch_probs = (1.0 - I) * (probs * switch) + I * (1.0 - (1.0 - probs) * switch)
    tags = [i + 1 for i in range(len(switch_probs))]
    return TransitionMatrix(switch_probs, tags, stationary_probs=probs.copy())
