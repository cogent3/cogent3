#!/usr/bin/env python
"""Writer for FASTA sequence format
"""

from cogent3.format.util import _AlignmentFormatter


__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight", "Gavin Huttley", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"


def alignment_to_fasta(alignment_dict, block_size=60, order=None):
    """Returns a Fasta string given an alignment."""
    order = order or []
    return FastaFormatter().format(alignment_dict, block_size, order)


class FastaFormatter(_AlignmentFormatter):
    def format(self, alignment_dict, block_size, order):
        """Format the alignment to Fasta.

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

        if len(alignment_dict) == 0:
            return ""

        return "".join(
            [
                ">%s\n%s"
                % (
                    seq,
                    self.wrap_string_to_block_size(
                        alignment_dict[seq], alt_block_size=block_size
                    ),
                )
                for seq in self.align_order
            ]
        )
