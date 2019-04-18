#!/usr/bin/env python
"""Alignment formatter template class
"""

__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Gavin Huttley", "Thomas La"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Thomas La"


class _AlignmentFormatter:
    """A virtual class for formatting alignments."""

    # other utility function
    def setblocksize(self, size):
        """Set the length of the sequence to be printed to each line.

        Arguments:
            - size: the sequence length to put on each line."""

        self.block_size = size

    def setaligninfo(self, alignmentdict, order=[]):
        """Set alignment attributes for writing.

        Arguments:
            - alignmentdict: dictionary of seqname -> seqstring
            - order: a list of seqname's in the order for writing
        """

        self.number_sequences = len(alignmentdict)
        # supersede the use of alignment length
        self.align_length = len(alignmentdict[list(alignmentdict.keys())[0]])

        if order != [] and len(order) == len(alignmentdict):
            # not testing contents - possibly should.
            self.align_order = order
        else:
            self.align_order = list(alignmentdict.keys()).sort()

    def slicestringinblocks(self, seqstring, altblocksize=0):
        """Return a list of string slices of specified length. No line returns.

        Arguments:
            - seqstring: the raw sequence string
            - altblocksize: the length of sequence for writing to each
              line, default (0) means default value specified by blocksize
              will be used.
        """

        if altblocksize:
            block_size = altblocksize
        else:
            block_size = self.block_size

        blocklist = []
        seqlength = len(seqstring)
        for block in range(0, seqlength, block_size):
            if block + block_size < seqlength:
                blocklist.append(seqstring[block: block + block_size])
            else:
                blocklist.append(seqstring[block:])

        return blocklist

    def wrapstringtoblocksize(self, seqstring, altblocksize=0):
        """Return sequence slices with line returns inserted at the end
        of each slice.

        Arguments:
            - seqstring: the raw sequence string
            - altblocksize: the length of sequence for writing to each
              line, default (0) means default value specified by blocksize
              will be used.
        """

        if altblocksize:
            self.block_size = altblocksize

        strlist = self.slicestringinblocks(seqstring, self.block_size)
        return '\n'.join(strlist) + "\n"
