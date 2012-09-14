#!/usr/bin/env python
"""
Writer for Clustal format.
"""
from cogent.core.alignment import SequenceCollection
from copy import copy

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def clustal_from_alignment(aln, interleave_len=None):
    """Returns a string in Clustal format.
    
        - aln: can be an Alignment object or a dict.
        - interleave_len: sequence line width.  Only available if sequences are
            aligned.
    """
    if not aln:
        return ''
    
     # get seq output order
    try:
        order = aln.RowOrder
    except:
        order = aln.keys()
        order.sort()
    
    seqs = SequenceCollection(aln)
    clustal_list = ["CLUSTAL\n"]
    
    if seqs.isRagged():
        raise ValueError,\
             "Sequences in alignment are not all the same length." +\
             "Cannot generate Clustal format."
    
    aln_len = seqs.SeqLen
    #Get all labels
    labels = copy(seqs.Names)
    
    #Find all label lengths in order to get padding.
    label_lengths = [len(l) for l in labels]
    label_max = max(label_lengths)
    max_spaces = label_max+4
    
    #Get ordered seqs
    ordered_seqs = [seqs.NamedSeqs[label] for label in order]
    
    if interleave_len is not None:
        curr_ix = 0
        while curr_ix < aln_len:
            clustal_list.extend(["%s%s%s"%(x,' '*(max_spaces-len(x)),\
                y[curr_ix:curr_ix+ \
                interleave_len]) for x,y in zip(order,ordered_seqs)])
            clustal_list.append("")
            curr_ix += interleave_len
    else:
        clustal_list.extend(["%s%s%s"%(x,' '*(max_spaces-len(x)),y) \
            for x,y in zip(order,ordered_seqs)])
        clustal_list.append("")
    
    return '\n'.join(clustal_list)    
        
        
