#!/usr/bin/env python

import numpy
import time
from cogent.align import *
from cogent.align.traceback import seq_traceback

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.2"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def _s2i(s):
    return numpy.array(['ATCG'.index(c) for c in s])

def test(f, decode, d, e, r=1):
            
    S = numpy.identity(4) * 7 - 4
    
    s2 = 'AAAATGCTT' * r
    s1 = 'AATTTTGCT' * r
    
    t0 = time.time()
    (score, states, maxpos, last_state) = f(_s2i(s1), _s2i(s2), S, d, e)
    print time.time() - t0
    
    alignment = seq_traceback(s1, s2, decode, states, maxpos, last_state, '.')
    print ' ', r*10, d, e, score
    print ' ', ''.join(alignment[0][:40])
    print ' ', ''.join(alignment[1][:40])
    
    #for a in gap_traceback(decode, states, maxpos, last_state):
    #    print a
        
        
if __name__ == '__main__':
    d = 2
    e = 1
    for r in [1]: #, 10, 20, 30, 50]:
        test(needle.needleman_wunsch, needle.decode_state, d, e, r)
