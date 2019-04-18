#!/usr/bin/env python
"""Writer for PAML sequence format
"""

from cogent3.format.formatter import _AlignmentFormatter

__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Gavin Huttley", "Thomas La"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Thomas La"


def alignment_to_paml(alignmentdict, block_size=60, order=[]):
    """Returns a Paml string given an alignment.
    """
    return PamlFormatter().format(alignmentdict, block_size, order)


class PamlFormatter(_AlignmentFormatter):

    def format(self, alignmentdict, block_size, order):
        """Format the alignment to Paml.

        Arguments:
            - alignmentdict: dict of seqname -> seqstring.
            - blocksize: the sequence length to write to each line,
              default is 60
            - order: optional list of sequence names, which order to
              print in.
        (Assumes order is a complete and correct list of names)
        """

        # setup
        if not order:
            order = list(alignmentdict.keys())
        self.setaligninfo(alignmentdict, order)
        self.setblocksize(block_size)

        header = '%d  %d\n' % (self.number_sequences, self.align_length)
        return header + ''.join(['%s\n%s' % (seq, self.wrapstringtoblocksize(
                alignmentdict[seq], altblocksize=block_size)) for seq in self.align_order])
