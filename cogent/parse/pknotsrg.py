#!/usr/bin/env python

"""Parser for  NUPACK output format 

If pseudoknotted first steam will be denoted by [] brackets and second steam
with {} brackets
"""

from string                import strip,split,atof
from cogent.util.transform import make_trans
from cogent.struct.rna2d   import Pairs,ViennaStructure
from cogent.struct.knots   import opt_single_random

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def pknotsrg_parser(lines=None,pseudo=True):
    """Parser for pknotsrg output format

    Returns a list containing: sequence, structure and energy
    ex: [[seq,[structure],energy]]

    pseudo - If True pairs will be returned with pseudoknots
             If False pairs will be returned without pseudoknots
    """
    result = []
    struct = str(lines[1]).strip('\n')
    seq = lines[0].strip()
    tmp_pairs,energy = to_pairs(struct)
    tmp_pairs.sort()
    if not pseudo:
        tmp_pairs = opt_single_random(tmp_pairs)
        tmp_pairs.sort()
    result.append([seq,tmp_pairs,energy])
    return result
          
primary_table = make_trans('{[]}','....')
first_table   = make_trans('({[]})','..()..')
second_table  = make_trans('([{}])','..()..')
    

def to_pairs(struct=None):
    """
    Converts structure string in to a pairs object.
    Starts by checking for pseudoknots if pseudoknotted it translates each 
    steam in to vienna notation and from there makes a pairs object. 
    Each pairs object is then joined to form the final pairs object of 
    the entire structure

    Returns a tuple of the pairs object and the energy
    """
    primary = first = second = struct.split(None,2)[0]
    energy = atof(struct.split(None,2)[1].strip('()'))

    if struct.__contains__('['): #Checks for first pseudoknot steam
        primary = ViennaStructure(primary.translate(primary_table))
        first = ViennaStructure(first.translate(first_table))
        pairs = primary.toPairs()
        pairs.extend(first.toPairs()) #Adds the first steam to pairs object
        if struct.__contains__('{'): #Checks for second pseudo steam
            second = ViennaStructure(second.translate(second_table))
            pairs.extend(second.toPairs())
    else: 
          primary = ViennaStructure(primary.translate(primary_table))
          pairs = primary.toPairs()

    return pairs,energy
    
