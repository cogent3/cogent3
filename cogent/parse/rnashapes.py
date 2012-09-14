#!/usr/bin/env python
#file: rnashapes_parser.py

"""
Author: Shandy Wikman (ens01svn@cs.umu.se)

Status: Development. According to future requirements behavior will be changed

Parser to parse RNAshapes output and returns list of lists [Seq,Pairs,Ene]

Revision History:
2006 Shandy Wikman created file
"""

from string                import split,strip,atof
from cogent.util.transform import make_trans
from cogent.struct.rna2d   import Pairs,ViennaStructure

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def RNAshapes_parser(lines=None,order=True):
    """
    Returns a list containing tuples of (sequence,pairs object,energy) for
    every sequence
    [[Seq,Pairs,Ene],[Seq,Pairs,Ene],...]
    Structures will be ordered by the structure energy by default, of ordered 
    isnt desired set order to False
    """
    result = lineParser(lines)
    if order:
        result = order_structs(result)
    return result

def lineParser(list=None):
    """
    Parses Lines from output and returns a list of tuples
    """
    result = []
    s = False
    seq = ''
    energy = 0.0
    pairs = ''

    for line in list:
        if len(seq)>1 and energy != 0.0 and len(pairs)>1:
            result.append([seq,pairs,energy])
            seq = energy = pairs = ''
        if line.startswith('>'):
            name = line.strip('>\n')
            s = True #signals that sequence is next line
        elif s:
            seq = line.strip()
            s = False
        elif line.startswith('-'):
            struct = line.split(None,2)[1].strip('\n')
            energy = atof(line.split(None,2)[0].strip('\n'))
            pairs = to_Pairs(struct)

    return result      
        

def to_Pairs(struct=None):
    """
    Converts a vienna structure into a pairs object
    Returns pairs object
    """
    struct = ViennaStructure(struct)
    pairs = struct.toPairs()

    return pairs

def order_structs(result):
    """
    order structure so that the structure whit highetst MFE(most negative)
    will be placed first and so on to the lowest MFE structure.

    """
    for i in result:
        i.reverse()
    result.sort()
    #result.reverse()  #to test with the lowest negetiv value as the best struct
    for i in result:
        i.reverse()

    return result
