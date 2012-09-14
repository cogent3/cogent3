#!/usr/bin/env python
#file: column_parser.py
"""Parser for column format

Works for the following column format:
; COL 1             label
; COL 2             residue
; COL 3             seqpos
; COL 4             alignpos
; COL 5             align_bp
; COL 6             certainty/seqpos_bp

Structure part separated by '; ------' and ends with '; ******' 
"""

from string import split
from cogent.struct.rna2d import Pairs
from cogent.struct.pairs_util import adjust_base

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def column_parser(lines):
    """Parser column format"""

    record = False
    result = []
    struct = []
    seq = ''
    for line in lines:
        if line.startswith('; ------'): #structure part beginns
            record = True
            continue
        if line.startswith('; ******'): #structure part ends
            record = False
            struct =  adjust_base(struct,-1)
            struct = Pairs(struct).directed()#remove duplicates
            struct.sort()

            result.append([seq,struct])
            struct = []
            seq = ''
            continue
        if record:
            sline = line.split()
            if sline[4] == '.': #skip not paired
                seq = ''.join([seq,sline[1]])
                continue
            seq = ''.join([seq,sline[1]])
            pair = (int(sline[3]),int(sline[4])) #(alignpos,align_bp)
            struct.append(pair)
        
    return result
