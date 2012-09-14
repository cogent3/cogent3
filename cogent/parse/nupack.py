#!/usr/bin/env python

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

def nupack_parser(lines=None,pseudo=True):
    """Parser for NUPACK output format

    pseudo - If True pseudoknot will be keept if False it will be removed
    """
    result = line_parser(lines,pseudo)
    return result
    
curly_to_dots_table = make_trans('{}','..')
bracket_to_dots_table = make_trans('()','..')
curly_to_bracket_table = make_trans('{}','()')

def line_parser(lines=None,pseudo=True):
    """Parser for nupack output format

    Returns list containing: sequence, paris and energy
    ex: [[seq,[struct],energy]]
    """
    record = False
    result = []
    SSEList = []    #Sequence,Structure,Energy
    for line in lines:
        if line.startswith('Error'):#Error no structure found
            result = [Pairs([])] #return empty pairs list
            return result
        if line.startswith('Sequence and a Minimum Energy Structure'):
            record = True
        elif record:
            line = line.strip('\n')
            SSEList.append(line)

    SSEList[1] = to_pairs(SSEList,pseudo) #pairs
    SSEList = SSEList[:3]
    SSEList[-1] = atof(SSEList[-1].split()[2]) #energy
    result.append(SSEList)
    return result

def to_pairs(list=None,pseudo=True):
    """
    Converts nupack structure string into pairs object
    pseudoknotted and not pseudoknotted.
    """
    tmp = list[1]
    pairs = []
 
    if list.__contains__('pseudoknotted!'):
        #since pseudoknotted is denoted by {} it divides into {} and () string
        #they are then turned into pairs lists seperatly and then the lists are 
        #joined to form the complete set of pairs
        first = second = tmp
        
        first = ViennaStructure(first.translate(curly_to_dots_table))
        second = second.translate(bracket_to_dots_table)
        second = ViennaStructure(second.translate(curly_to_bracket_table))
        
        pairs = first.toPairs()
        pairs.extend(second.toPairs())
        pairs.sort()

        if not pseudo:
            pairs = opt_single_random(pairs)
            pairs.sort()
    else:
        structure = ViennaStructure(tmp.translate(curly_to_bracket_table))
        pairs = structure.toPairs()
        pairs.sort()

    return pairs

