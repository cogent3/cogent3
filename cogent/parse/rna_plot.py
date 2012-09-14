#!/usr/bin/env python
"""
Parser for RNAPlot postscript output.
"""

from cogent.parse.record_finder import LabeledRecordFinder

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def get_sequence(lines):
    """Returns sequence string.
    """
    if lines:
        seq_line = ''.join(lines)
        seq = seq_line.split('(')[1].split(') def')[0].strip('\n\\')
    else:
        seq = ''
    return seq
    
def get_coordinates(lines):
    """Returns list of coordinates.
    """
    coords = []
    for l in lines:
        if l.startswith('] def'):
            break
        elif l.startswith('['):
            coords.append(map(float,l.strip('\n[]').split()))
    return coords

def get_pairs(lines):
    """Returns list of pairing indices.
    
        - Changes indices to zero based rather than 1 based.
    """
    pairs = []
    for l in lines:
        if l.startswith('] def'):
            break
        elif l.startswith('['):
            first,second = l.strip('\n[]').split()
            pairs.append([int(first)-1,int(second)-1])
    return pairs

def RnaPlotParser(lines):
    """Returns sequence, coordinates, and pairing indices.
    """
    sequence = ''
    coordinates = []
    pairs = []
    if lines:
        #Split on sequence block
        sequence_finder = LabeledRecordFinder(is_label_line=\
            lambda x: x.startswith('/sequence'))
        prefix, seq_block = list(sequence_finder(lines))
        
        #split on coordinate block
        coordinate_finder = LabeledRecordFinder(is_label_line=\
            lambda x: x.startswith('/coor'))
        #sequence block is first item in list
        sequence_block, coord_block = list(coordinate_finder(seq_block))
        
        #split on pairs block
        pairs_finder = LabeledRecordFinder(is_label_line=\
            lambda x: x.startswith('/pairs'))
        #coordinate block is first item in list
        coordinate_block, pairs_block = list(pairs_finder(coord_block))
        
        sequence = get_sequence(sequence_block)
        coordinates = get_coordinates(coordinate_block)
        pairs = get_pairs(pairs_block)
    
    return sequence, coordinates, pairs
