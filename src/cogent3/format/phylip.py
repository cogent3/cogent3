#!/usr/bin/env python
"""Writer for PHYLIP sequence format
"""

from cogent3.format.util import _AlignmentFormatter


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def alignment_to_phylip(alignment_dict, block_size=60, order=None):
    """Returns a Phylip string given an alignment."""
    return PhylipFormatter().format(
        alignment_dict, block_size, [] if order is None else order
    )


class PhylipFormatter(_AlignmentFormatter):
    def format(self, alignment_dict, block_size, order):
        """Format the alignment to Phylip.

        Parameters
        ----------
        alignment_dict
            dict of seq_name
        block_size
            the sequence length to write to each line,
            default is 60
        order
            optional list of sequence names, which order to
            print in.
            (Assumes complete and correct list of names)

        """
        # setup
        if not order:
            order = list(alignment_dict.keys())
        self.set_align_info(alignment_dict, order)
        self.set_block_size(block_size)

        # header
        header = "%d  %d\n" % (self.number_sequences, self.align_length)

        seqs = []

        # sequences (pretty much as writ by Gavin)

        for seq_name in self.align_order:
            seq = alignment_dict[seq_name]
            for block in range(0, self.align_length, self.block_size):
                if not block:
                    # write the otu name
                    if len(seq_name) > 9:
                        prefix = "%-10s" % seq_name[:9]
                    else:
                        prefix = "%-10s" % seq_name
                else:
                    prefix = " " * 10

                if block + self.block_size > self.align_length:
                    to = self.align_length
                else:
                    to = block + self.block_size

                seqs.append(f"{prefix}{seq[block:to]}\n")

        return header + "".join(seqs)
