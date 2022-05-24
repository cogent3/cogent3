#!/usr/bin/env python

import time

import numpy

from cogent3 import DNA
from cogent3.align.align import classic_align_pairwise, make_dna_scoring_dict


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def _s2i(s):
    return numpy.array(["ATCG".index(c) for c in s])


def test(r=1, **kw):
    S = make_dna_scoring_dict(10, -1, -8)

    seq2 = DNA.make_seq("AAAATGCTTA" * r)
    seq1 = DNA.make_seq("AATTTTGCTG" * r)

    t0 = time.clock()
    try:
        # return_alignment is False in order to emphasise the quadratic part of
        # the work.
        aln = classic_align_pairwise(
            seq1, seq2, S, 10, 2, local=False, return_alignment=False, **kw
        )
    except ArithmeticError:
        return "*"
    else:
        t = time.clock() - t0
        return int((len(seq1) * len(seq2)) / t / 1000)


if __name__ == "__main__":
    d = 2
    e = 1
    options = [(False, False), (True, False), (False, True)]
    template = "%10s " * 4
    print("               1000s positions^2 per second")
    print(template % ("length", "simple", "logs", "scaled"))
    for r in [100, 200, 300, 400, 500]:
        times = [test(r, use_logs=l, use_scaling=s) for (l, s) in options]
        print(template % tuple([r * 10] + times))
