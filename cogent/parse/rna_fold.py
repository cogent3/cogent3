#!/usr/bin/env python
"""Parses RNAfold dot plot output file.
"""
from string import strip
from cogent.parse.record_finder import LabeledRecordFinder

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def RnaFoldParser(lines):
    """Returns a tuple containing sequence and dot plot indices.
    
        (sequence, (index1, index2, pair probability))
    """
    sequence = ''
    indices = []
    #Make sure lines is not empty
    if lines:
        #Get the block of lines that starts with /sequence
        sequence_block = LabeledRecordFinder(\
            lambda x: x.startswith('/sequence'))
        #only care about the second element in the result
        seq_block = list(sequence_block(lines))[1]
        #Get the sequence from the block
        sequence = getSequence(seq_block)
        #Get the indices and pair probabilites from the block
        indices = getIndices(seq_block)
    return (sequence, indices)


def getSequence(lines):
    """Returns sequence that RNAfold dot plot represents, given lines.
    """
    sequence_pieces = []
    #For each line after the first, containing /sequence
    for line in map(strip, lines[1:]):
        #If the line which denotes the end of a sequence is not reached
        if line.startswith(')'):
            break
        else:
            #Strip of whitespace and add to list
            sequence_pieces.append(line.replace('\\',''))
    return ''.join(sequence_pieces)
    
def getIndices(lines):
    """Returns list of tuples: (index1,index2,pair probability).
    """
    index_list = []
    #For each line that ends with 'ubox' (which denotes lines with indices)
    for line in filter(lambda x: x.endswith('ubox'), lines):
            #split on whitespace
            up_index,down_index,probability,ubox = line.split()
            #Build tuple with indices and pair probability
            index_list.append((int(up_index), int(down_index), \
                float(probability)))
    return index_list
