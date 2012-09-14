#!/usr/bin/env python

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def phylip_from_alignment(aln, generic_label=True, make_seqlabel=None):
    """returns a phylip formatted string and an ID map of new to original
    sequence names. Sequences are sequential.
    
    Fails if all sequences are not the same length.
    
    Arguments:
        - generic_label: if true then numbered seq labels are generated
        - make_seqlabel: callback function that takes the seq object and
          returns a label str. The user must ensure these are correct for Phylip
          format.
    """
    assert generic_label or make_seqlabel is not None
    if aln.isRagged():
        raise ValueError, "Sequences in alignment are not all the same " +\
                          "length. Cannot generate PHYLIP format."
    num_seqs = len(aln.Seqs)
    if not aln or not num_seqs:
        return ""
    
    phylip_out = ["%d %d" % (num_seqs, aln.SeqLen)]
    id_map = {}
    cur_seq_id = 1
    
    for seq_name, seq in zip(aln.Names, aln.Seqs):
        if make_seqlabel is not None:
            label = make_seqlabel(seq)
        elif generic_label:
            label = "seq%07d" % cur_seq_id
        id_map[label] = seq_name
        phylip_out.append("%s %s" % (label, seq))
        cur_seq_id += 1
    
    return '\n'.join(phylip_out), id_map
