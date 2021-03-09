import math

from unittest import TestCase, main

from cogent3 import DNA, make_aligned_seqs
from cogent3.evolve.alignment_likelihood import logLikelihood


alignment = make_aligned_seqs(data=[("s1", "A"), ("s2", "G"), ("s3", "A")], moltype=DNA)
alignment2 = make_aligned_seqs(
    data=[("s1", "AT"), ("s2", "AT"), ("s3", "AC")], moltype=DNA
)
alignment3 = make_aligned_seqs(
    data=[("s1", "AGC"), ("s2", "ATA"), ("s3", "ACC")], moltype=DNA
)
rateMatrix = [
    [0.95, 0.01, 0.01, 0.03],
    [0.005, 0.98, 0.01, 0.005],
    [0.005, 0.02, 0.95, 0.015],
    [0.001, 0.003, 0.006, 0.99],
]
pi = [0.3, 0.2, 0.3, 0.2]

import pytest

from numpy.testing import assert_allclose, assert_almost_equal, assert_equal


def test_alignment():
    assert_allclose((logLikelihood(pi, rateMatrix, alignment2)), -3.278666124)


def test_alignment2():
    assert_almost_equal(logLikelihood(pi, rateMatrix, alignment3), -7.9027699)


def test_alignment3():
    assert_almost_equal(logLikelihood(pi, rateMatrix, alignment), -2.905530618)
