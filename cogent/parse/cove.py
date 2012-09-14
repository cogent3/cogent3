#!/usr/bin/env python

from string                import split,strip
from cogent.util.transform import make_trans
from cogent.struct.rna2d   import ViennaStructure,wuss_to_vienna

__author__ = "Shandy Wikman"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Shandy Wikman"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Shandy Wikman"
__email__ = "ens01svn@cs.umu.se"
__status__ = "Development"

def coves_parser(data=None):
    """ Parser for coves output using option -m

    Takes lines as input.
    Returns a list of lists(containing structure and pairs data)
    ex: [[seq1,pairs1],[seq2,pairs2],...]
    """
    c = 1
    count = 0
    tmp_seq = ''
    tmp_struct = ''
    seq = ''
    struct = ''
    result = []
    seqNames = NameList(data)
    for name in seqNames: #parse each sequence in turn
        name = '%s ' % (name) #add blank to differ ex seq1 from seq10 
        count+=1
        if count != 1:
            struct, seq = remove_gaps(struct,seq)
            pairs = convert_to_vienna(struct)
            result.append([seq,pairs])
            seq = ''
            struct = ''
        for line in data:
            line = str(line)
            line = line.strip()
            sline = line.split(None,1)
            if c==1: #sequence line (every other line)
                if line.startswith(name):
                    c=0
                    tmp_seq = sline[-1]
                    seq = ''.join([seq,tmp_seq])
            elif c==0: #struct line (every other line)
                if line.startswith(name):
                    c=1
                    tmp_struct = sline[-1]
                    struct = ''.join([struct,tmp_struct])
    struct,seq = remove_gaps(struct,seq)
    pairs = convert_to_vienna(struct)
    result.append([seq,pairs])

    return result

cove_to_vienna_table = make_trans('><','()')

def remove_gaps(struct,seq):
    """Remove gaps function

    Some results comes with gaps that need to be removed
    """
    seq = seq.replace('-','')
    tmp_struct = struct.split() #removes gaps
    tmp = ''
    for i in range(len(tmp_struct)):
        tmp = ''.join([tmp,tmp_struct[i]]) #put struct parts together
    struct = tmp
    if len(struct) != len(seq): #check so that struct and seq match in length
        raise ValueError, 'Sequence length don\'t match structure length'
    return struct,seq

def NameList(data=None):
    """
    Takes coves results and retrieves the sequence names for further parsing 
    """
    nameList = []
    if not isinstance(data,list):
        data = open(data).readlines()
    for line in data:
        if line.__contains__('bits'): #every unique sequense begins with 'bits'
            line = line.split()
            nameList.append(line[-1])
    return nameList

def convert_to_vienna(data):
    """
    Converts into vienna dot bracket format, >< to () 
    """
    try:
        return toPairs(ViennaStructure(data.translate(cove_to_vienna_table)))
    except IndexError:
        return ''
def toPairs(vienna):
    """
    Converts a vienna structure to a pairs obejct
    """

    pairs = vienna.toPairs()
    return pairs
    
