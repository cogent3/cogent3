#!/usr/bin/env python
"""Alignment formatter template class
"""

__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Gavin Huttley", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


class _AlignmentFormatter:
    """A virtual class for formatting alignments."""

    block_size = None
    number_sequences = None
    align_length = None
    align_order = None

    # other utility function
    def set_block_size(self, size):
        """Set the length of the sequence to be printed to each line.

        Parameters
        ----------
        size
            the sequence length to put on each line.

        """

        self.block_size = size

    def set_align_info(self, alignment_dict, order=None):
        """Set alignment attributes for writing.

        Parameters
        ----------
        alignment_dict
            dictionary of seq_name
        order
            a list of seq_name's in the order for writing

        """

        order = order or []

        self.number_sequences = len(alignment_dict)
        # supersede the use of alignment length
        if len(alignment_dict) > 0:
            self.align_length = len(alignment_dict[list(alignment_dict.keys())[0]])

        if order != [] and len(order) == len(alignment_dict):
            # not testing contents - possibly should.
            self.align_order = order
        elif order is None or order == []:
            self.align_order = []
        else:
            self.align_order = list(alignment_dict.keys()).sort()

    def slice_string_in_blocks(self, seq_string, alt_block_size=0):
        """Return a list of string slices of specified length. No line returns.

        Parameters
        ----------
        seqstring
            the raw sequence string
        altblocksize
            the length of sequence for writing to each
            line, default (0) means default value specified by blocksize
            will be used.

        """

        if alt_block_size:
            block_size = alt_block_size
        else:
            block_size = self.block_size

        block_list = []
        seq_length = len(seq_string)
        for block in range(0, seq_length, block_size):
            if block + block_size < seq_length:
                block_list.append(seq_string[block : block + block_size])
            else:
                block_list.append(seq_string[block:])

        return block_list

    def wrap_string_to_block_size(self, seq_string, alt_block_size=0):
        """Return sequence slices with line returns inserted at the end
        of each slice.

        Parameters
        ----------
        seq_string
            the raw sequence string
        alt_block_size
            the length of sequence for writing to each
            line, default (0) means default value specified by block size
            will be used.

        """

        if alt_block_size:
            self.block_size = alt_block_size

        str_list = self.slice_string_in_blocks(seq_string, self.block_size)
        return "\n".join(str_list) + "\n"
