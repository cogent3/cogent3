#!/usr/bin/env python

import numpy
import time
from cogent import DNA
from cogent.align.align import classic_align_pairwise, make_dna_scoring_dict

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def _s2i(s):
    return numpy.array(['ATCG'.index(c) for c in s])

def test(r=1, **kw):   
    S = make_dna_scoring_dict(10, -1, -8)
    
    seq2 = DNA.makeSequence('AAAATGCTTA' * r)
    seq1 = DNA.makeSequence('AATTTTGCTG' * r)
    
    t0 = time.time()
    aln = classic_align_pairwise(seq1, seq2, S, 10, 2, local=False, **kw)
    t = time.time() - t0
    return (len(seq1)*len(seq2))/t
    
    print t 
        
if __name__ == '__main__':
    d = 2
    e = 1
    options = [(False, False), (True, False), (False, True)]
    template = "%10s " * 4
    print "                1000s positions per second"
    print template % ("size", "simple", "logs", "scaled")
    for r in [50, 100, 200, 500]:
        times = [test(r, use_logs=l, use_scaling=s) for (l,s) in options]
        print template % tuple([r*10] + [int(t/1000) for t in times])
