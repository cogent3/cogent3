#!/usr/bin/env python
"""Writer for GDE sequence format
"""

from cogent3.format.util import _AlignmentFormatter


def alignment_to_gde(alignment_dict, block_size=60, order=None):
    """Returns a Gde string given an alignment."""
    return GDEFormatter().format(
        alignment_dict, block_size, [] if order is None else order
    )


class GDEFormatter(_AlignmentFormatter):
    def format(self, alignment_dict, block_size, order):
        """Format the alignment to GDE.

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

        return "".join(
            [
                "%s%s\n%s"
                % (
                    "%",
                    seq,
                    self.wrap_string_to_block_size(
                        alignment_dict[seq], alt_block_size=block_size
                    ),
                )
                for seq in self.align_order
            ]
        )
