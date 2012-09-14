#!/usr/bin/env python

"""Parser for ct rna secondary structure format 

Works on ct files containing one or more structures..
supports: Carnac
          dynalign
          mfold
          sfold
          unafold
          knetfold

Should work on all ct formats conforming to format: 
header, structure, header, structure ...

Header is line beginning every structure, containing length,energy,input file: 
72 ENERGY = -23.4 trna_phe.fasta

currently only works on multiple structures files if header lines contain 
the word 'Structure', 'ENERGY' or 'dG'. Further support added as needed

Convention of Connect format(ct) is to include 'ENERGY = value' as header above 
(value left blank if not applicable)
"""

from string              import split,atof
from cogent.struct.rna2d import Pairs
from cogent.struct.pairs_util  import adjust_base

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def ct_parser(lines=None):
    """Ct format parser

    Takes lines from a ct file as input
    
    Returns a list containing sequence,structure and if available the energy.
    [[seq1,[struct1],energy1],[seq2,[struct2],energy2],...]
    """

    count = 0
    length = ''
    energy = None
    seq = ''
    struct = []
    result = []

    for line in lines:
        count+=1
        sline = line.split(None,6) #sline = split line
        if count==1 or new_struct(line):#first line or new struct line.
            if count > 1:
                struct = adjust_base(struct,-1)
                struct = Pairs(struct).directed()
                struct.sort()
                if energy is not None:
                    result.append([seq,struct,energy])
                    energy = None
                else:
                    result.append([seq,pairs])
                struct = []
                seq = ''
            #checks if energy for predicted struct is given
            if sline.__contains__('dG') or sline.__contains__('ENERGY'):
                energy = atof(sline[3])
            if sline.__contains__('Structure'):
                energy = atof(sline[2])
        else:
            seq = ''.join([seq,sline[1]])
            if not int(sline[4]) == 0:#unpaired base
                pair = ( int(sline[0]),int(sline[4]) )
                struct.append(pair) 
    #structs are one(1) based, adjust to zero based
    struct = adjust_base(struct,-1)
    struct = Pairs(struct).directed()
    struct.sort()

    if energy is not None:
        result.append([seq,struct,energy])
    else:
        result.append([seq,struct])
    return result 

def new_struct(line):
    """
    Determines if a new structure begins on line in question.

    Currently only works for multiple structure files containing these key
    words in their header. 

    Convention of Connect format (ct format) is to include 'ENERGY = value' 
    (value left blank if not applicable)

    Support for additional formats will be added as needed
    """
    answer=False
    if 'Structure' in line or 'dG' in line or 'ENERGY' in line:
        answer = True

    return answer
