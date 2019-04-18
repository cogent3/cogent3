#!/usr/bin/env python

from cogent3.format.formatter import _AlignmentFormatter

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La"]
__license__ = "GPL"
__version__ = "3.0a2"
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
    if aln.is_ragged():
        raise ValueError("Sequences in alignment are not all the same " +
                         "length. Cannot generate PHYLIP format.")
    num_seqs = len(aln.seqs)
    if not aln or not num_seqs:
        return ""

    phylip_out = ["%d %d" % (num_seqs, aln.seq_len)]
    id_map = {}
    cur_seq_id = 1

    for seq_name, seq in zip(aln.names, aln.seqs):
        if make_seqlabel is not None:
            label = make_seqlabel(seq)
        elif generic_label:
            label = "seq%07d" % cur_seq_id
        id_map[label] = seq_name
        phylip_out.append("%s %s" % (label, seq))
        cur_seq_id += 1

    return '\n'.join(phylip_out), id_map


def alignment_to_phylip(alignmentdict, block_size=60, order=[]):
    """Returns a Phylip string given an alignment.
    """
    return PhylipFormatter().format(alignmentdict, block_size, order)


class PhylipFormatter(_AlignmentFormatter):

    def format(self, alignmentdict, block_size, order):
        """Format the alignment to Phylip.

        Arguments:
            - alignmentdict: dict of seqname -> seqstring.
            - blocksize: the sequence length to write to each line,
              default is 60
            - order: optional list of sequence names, which order to
              print in.
        (Assumes complete and correct list of names)
        """
        # setup
        if not order:
            order = list(alignmentdict.keys())
        self.setaligninfo(alignmentdict, order)
        self.setblocksize(block_size)

        # header
        header = '%d  %d\n' % (self.number_sequences, self.align_length)

        seqs = []

        # sequences (pretty much as writ by Gavin)

        for seqname in self.align_order:
            seq = alignmentdict[seqname]
            for block in range(0, self.align_length, self.block_size):
                if not block:
                    # write the otu name
                    if len(seqname) > 9:
                        prefix = '%-10s' % seqname[:9]
                    else:
                        prefix = '%-10s' % seqname
                else:
                    prefix = ' ' * 10

                if block + self.block_size > self.align_length:
                    to = self.align_length
                else:
                    to = block + self.block_size

                seqs.append('%s%s\n' % (prefix, seq[block:to]))

        return header + ''.join(seqs)
