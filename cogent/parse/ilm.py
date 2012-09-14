#!/usr/bin/env python

from cogent.struct.rna2d      import Pairs
from cogent.struct.knots      import opt_single_random
from cogent.struct.pairs_util import adjust_base

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def ilm_parser(lines=None,pseudo=True):
    """Ilm format parser

    Takes lines as input and returns a list with Pairs object.
    Pseudo - if True returns pairs with possible pseudoknot
             if False removes pseudoknots       
    """
    pairs = []
    for line in lines:
        if line.startswith('Final') or len(line)==1:#skip these lines
            continue
        line = line.strip('\n')
        line = map(int,line.split(None,2))
        if line[1] == 0:
            continue #Skip this line, not a pair
        else:
            pairs.append(line) 

    pairs = adjust_base(pairs,-1)
    tmp = Pairs(pairs).directed()
    tmp.sort()
    if not pseudo:
        tmp = opt_single_random(tmp)
        tmp.sort()
    result = []
    result.append(tmp)

    return result
