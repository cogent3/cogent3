#!/usr/bin/env python
"""
Writer for Clustal format.
"""
from copy import copy

from cogent3.core.alignment import SequenceCollection


__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"


def clustal_from_alignment(aln, wrap=None):
    """
    Parameters
    ----------
    aln
        can be an Alignment object or a dict
    wrap
        sequence line width.  Only available if sequences are
        aligned.

    Returns
    -------
    Returns a string in Clustal format
    """
    if not aln:
        return ""

    # get seq output order
    try:
        order = aln.RowOrder
    except:
        order = list(aln.keys())
        order.sort()

    seqs = SequenceCollection(aln)
    clustal_list = ["CLUSTAL\n"]

    if seqs.is_ragged():
        raise ValueError(
            "Sequences in alignment are not all the same length."
            + "Cannot generate Clustal format."
        )

    aln_len = seqs.seq_len
    # Get all labels
    labels = copy(seqs.names)

    # Find all label lengths in order to get padding.
    label_lengths = [len(l) for l in labels]
    label_max = max(label_lengths)
    max_spaces = label_max + 4

    # Get ordered seqs
    ordered_seqs = [seqs.named_seqs[label] for label in order]

    if wrap is not None:
        curr_ix = 0
        while curr_ix < aln_len:
            clustal_list.extend(
                [
                    f"{x}{' ' * (max_spaces - len(x))}{y[curr_ix:curr_ix + wrap]}"
                    for x, y in zip(order, ordered_seqs)
                ]
            )
            clustal_list.append("")
            curr_ix += wrap
    else:
        clustal_list.extend(
            [
                f"{x}{' ' * (max_spaces - len(x))}{y}"
                for x, y in zip(order, ordered_seqs)
            ]
        )
        clustal_list.append("")

    return "\n".join(clustal_list)
