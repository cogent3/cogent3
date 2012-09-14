#!/usr/bin/env python

from cogent.struct.rna2d   import ViennaStructure,wuss_to_vienna
from cogent.util.transform import make_trans

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

to_vienna_table = make_trans('><','()')

def consan_parser(lines):
    """
    Takes a series of lines as input.

    Returns a list containing alignment and structure
    ex: [{alignment},[structure]]
    """
    seqs = []
    struct = ''
    pairs = []
    alignment = {}
    for line in lines:
        if sequence(line):
            line = line.split()
            name = line[0].strip()
            seq = line[1].strip()
            #add sequence to alignment
            if name in alignment:
                alignment[name] += seq
            else:
                alignment[name] = seq
        elif line.startswith('#=GC SS_cons'):
            line = line.split()
            struct += line[2].strip()
    pairs = convert_to_pairs(struct)
    return [alignment, pairs]
    
def convert_to_pairs(data):
    """
    Converts format >< to () format, viennaformat. 
    """
    try:
        vienna = ViennaStructure(data.translate(to_vienna_table))
        return toPairs(vienna)
    except IndexError:
        return ''
def toPairs(vienna):
    """
    Converts a vienna structure to a pairs obejct
    """

    pairs = vienna.toPairs()
    return pairs
    
def sequence(line):
    """Determines if line is a sequence line
    """
    answer = False
    if not len(line) == 1 and not line.startswith('#') and not line.startswith('/') and not line.startswith('Using'):
        answer = True
    return answer
