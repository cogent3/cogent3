#!/usr/bin/env python

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0.alpha"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def nexus_from_alignment(aln, seq_type, interleave_len=50):
    """returns a nexus formatted string

    Arguments:
        - seq_type: dna, rna, or protein
        - interleave_len: the line width"""
    if aln.is_ragged():
        raise ValueError("Sequences in alignment are not all the same " +
                         "length. Cannot generate NEXUS format.")
    num_seq = len(aln.Seqs)
    if not aln or not num_seq:
        return ""
    aln_len = aln.seq_len
    nexus_out = ["#NEXUS\n\nbegin data;"]
    nexus_out.append("    dimensions ntax=%d nchar=%d;" % (num_seq,
                                                           aln_len))
    nexus_out.append("    format datatype=%s interleave=yes missing=? " %
                     seq_type + "gap=-;")
    nexus_out.append("    matrix")
    cur_ix = 0
    names_seqs = sorted(aln.named_seqs.items())
    while cur_ix < aln_len:
        nexus_out.extend(["    %s    %s" % (x, y[cur_ix:cur_ix +
                                                 interleave_len]) for x, y in names_seqs])
        nexus_out.append("")
        cur_ix += interleave_len
    nexus_out.append("    ;\nend;")

    return '\n'.join(nexus_out)
