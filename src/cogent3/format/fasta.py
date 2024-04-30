#!/usr/bin/env python
"""Writer for FASTA sequence format
"""

from typing import Optional

from cogent3.format.util import _AlignmentFormatter


def alignment_to_fasta(
    alignment_dict: dict[str, str],
    block_size: int = 60,
    order: Optional[list[str]] = None,
) -> str:
    """Returns a Fasta string given an alignment.

    Parameters
        ----------
        alignment_dict
            dict of seq_name
        block_size
            the sequence length to write to each line,
            by default 60
        order
            optional list of sequence names, which order to print in.
            Assumes complete and correct list of names,
            by default None

        Returns
        -------
        The alignment in the Fasta format.
    """
    order = order or []
    return FastaFormatter().format(alignment_dict, block_size, order)


class FastaFormatter(_AlignmentFormatter):
    def format(
        self,
        alignment_dict: dict[str, str],
        block_size: int = 60,
        order: Optional[list[str]] = None,
    ) -> str:
        """Format the alignment to Fasta.

        Parameters
        ----------
        alignment_dict
            dict of seq_name
        block_size
            the sequence length to write to each line,
            by default 60
        order
            optional list of sequence names, which order to print in.
            Assumes complete and correct list of names,
            by default None

        Returns
        -------
        The alignment in the Fasta format.
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
