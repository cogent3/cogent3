#!/usr/bin/env python
"""
Writer for Stockholm format.
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

def stockholm_from_alignment(aln, interleave_len=None, GC_annotation=None):
    """Returns a string in Stockholm format.
    
        - aln: can be an Alignment object or a dict.
        - interleave_len: sequence line width.  Only available if sequences are
            aligned.
        - GC_annotation: dict containing Per-column annotation {<tag>:<s>},
            added to Stockholm file in the following format: #=GC <tag> <s>
            - <s> is an aligned text line of annotation type <tag>.
            - #=GC lines are associated with a sequence alignment block;
            - <s> is aligned to the residues in the alignment block, and has the 
            same length as the rest of the block. #=GC lines are
            placed at the end of each block. 
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
    stockholm_list = ["# STOCKHOLM 1.0\n"]
    
    if seqs.isRagged():
        raise ValueError,\
             "Sequences in alignment are not all the same length." +\
             "Cannot generate Stockholm format."
    
    aln_len = seqs.SeqLen
    #Get all labels
    labels = copy(seqs.Names)
    
    #Get ordered seqs
    ordered_seqs = [seqs.NamedSeqs[label] for label in order]

    if GC_annotation is not None:
        GC_annotation_list = \
            [(k,GC_annotation[k]) for k in sorted(GC_annotation.keys())]
        #Add GC_annotation to list of labels.
        labels.extend(['#=GC '+ k for k in GC_annotation.keys()])
        for k,v in GC_annotation.items():
            if len(v) != aln_len:
                raise ValueError, """GC annotation %s is not same length as alignment. Cannot generate Stockholm format."""%(k)
    
    #Find all label lengths in order to get padding.
    label_lengths = [len(l) for l in labels]
    label_max = max(label_lengths)
    max_spaces = label_max+4
        
    if interleave_len is not None:
        curr_ix = 0
        while curr_ix < aln_len:
            stockholm_list.extend(["%s%s%s"%(x,' '*(max_spaces-len(x)),\
                y[curr_ix:curr_ix+ \
                interleave_len]) for x,y in zip(order, ordered_seqs)])
            if GC_annotation is not None:
                stockholm_list.extend(["#=GC %s%s%s"%(x,\
                    ' '*(max_spaces-len(x)-5),\
                    y[curr_ix:curr_ix + interleave_len]) for x,y in\
                    GC_annotation_list])
            stockholm_list.append("")
            curr_ix += interleave_len
    else:
        stockholm_list.extend(["%s%s%s"%(x,' '*(max_spaces-len(x)),y) \
            for x,y in zip(order, ordered_seqs)])
        if GC_annotation is not None:
            stockholm_list.extend(["#=GC %s%s%s"%(x,' '*(max_spaces-len(x)-5),\
                y) for x,y in GC_annotation_list])
        stockholm_list.append("")
    
    return '\n'.join(stockholm_list)+'//'
        
        
