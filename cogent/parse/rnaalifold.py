#!/usr/bin/env python

from string                import strip,split,atof
from cogent.struct.rna2d   import ViennaStructure,Pairs

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman","Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def MinimalRnaalifoldParser(lines):
    """MinimalRnaalifoldParser.
        returns lists of sequence, structure_string, energy
    """
    res = []
    if lines:
        for i in range(0,len(lines),2):
            seq = lines[i].strip()
            struct,energy = lines[i+1].split(" (")
            energy = float(energy.split('=')[0].strip(' \n)'))
            res.append([seq,struct,energy])
    return res
        

def rnaalifold_parser(lines=None):
    """Parser for rnaalifold stdout output

    Returns a list containing: sequence,structure(pairs object) and energy
    Ex: [seq,[struct],energy]
    """
    result = line_parser(lines)
    return result
    
def line_parser(lines=None):
    """Parses RNAalifold output line for line """
    s = False
    seq = ''
    energy = ''
    pairs = ''
    result = []
    for line in lines:
        if len(line)>1 and s==False:
            seq = line.strip()
            s = True
        elif s == True:
            s=False
            struct = line.split(None,2)[0].strip('\n')
            energy = atof(line.split(' (',1)[1].split(None,1)[0].strip())
            pairs = to_pairs(struct)
            pairs.sort()
            
            result.append([seq,pairs,energy])
    return result

def to_pairs(struct=None):
    """
    Converts a vienna structure into a pairs object
    Returns pairs object

    pairs functions tested in test for rna2d.py
    """
    struct = ViennaStructure(struct)
    pairs = struct.toPairs()

    return pairs
