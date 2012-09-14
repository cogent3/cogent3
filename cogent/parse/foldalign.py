#!/usr/bin/env python

from string import split
from cogent.struct.rna2d  import Pairs,ViennaStructure
from cogent.struct.pairs_util import adjust_base
from cogent.parse.column import column_parser

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def foldalign_parser(lines,col=True):
    """Parser foldalign output"""

    data = lines
    
    if col:
        return column_parser(data)
    else:
        return find_struct(data)


def find_struct(lines):
    """Finds structures in output data"""

    struct = ''
    name1 = ''
    name2 = ''
    seq1 = ''
    seq2 = ''
    result = []
    for line in lines:
        if line.startswith('; ========'):
            break
        if line.startswith('; ALIGNING'):
            line = line.split()
            name1 = line[2]
            name2 = line[4]
            continue
        if line.startswith('; ALIGN               %s' % name1):
            line = line.split()[3:]
            line = ''.join(line)
            seq1 = ''.join([seq1,line])
            continue
        if line.startswith('; ALIGN               %s' % name2):
            line = line.split()[3:]
            line = ''.join(line)
            seq2 = ''.join([seq2,line])
            continue
        if line.startswith('; ALIGN               Structure'):
            line = line.split()[3:]
            line = ''.join(line)
            struct = ''.join([struct,line])
            continue
    struct = ViennaStructure(struct).toPairs()
    struct.sort()
    result.append([struct,seq1,seq2])
    return result


        
            

