#!/usr/bin/env python
"""Writer for GDE sequence format
"""

from cogent3.format.formatter import _AlignmentFormatter

__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Gavin Huttley", "Thomas La"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Thomas La"


def alignment_to_gde(alignmentdict, block_size=60, order=[]):
    """Returns a GDE string given an alignment.
    """
    return GDEFormatter().format(alignmentdict, block_size, order)


class GDEFormatter(_AlignmentFormatter):

    def format(self, alignmentdict, block_size, order):
        """Format the alignment to GDE.

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

        return ''.join(['%s%s\n%s' % ("%", seq, self.wrapstringtoblocksize(alignmentdict[seq],
                        altblocksize=block_size)) for seq in self.align_order])


