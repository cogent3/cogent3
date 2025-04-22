#!/usr/bin/env python
"""
Writer for Clustal format.
"""

import os

_NEW_TYPE = "COGENT3_NEW_TYPE" in os.environ


def clustal_from_alignment(
    data: dict[str, str],
    wrap: int | None = None,
    order: list[str] | None = None,
) -> str:
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
    if not data:
        return ""

    # get seq output order
    order = order or list(data)

    clustal_list = ["CLUSTAL\n"]
    aln_len = len(data[order[0]])
    for l, s in data.items():
        if len(s) != aln_len:
            msg = f"length of {l!r}={len(s)} != {aln_len}"
            raise ValueError(msg)

    # Find all label lengths in order to get padding.
    label_max = len(max(order, key=len))
    max_spaces = label_max + 4

    # Get ordered seqs
    ordered_seqs = [data[label] for label in order]

    if wrap is not None:
        curr_ix = 0
        while curr_ix < aln_len:
            clustal_list.extend(
                [
                    f"{x}{' ' * (max_spaces - len(x))}{y[curr_ix : curr_ix + wrap]}"
                    for x, y in zip(order, ordered_seqs, strict=False)
                ],
            )
            clustal_list.append("")
            curr_ix += wrap
    else:
        clustal_list.extend(
            [
                f"{x}{' ' * (max_spaces - len(x))}{y}"
                for x, y in zip(order, ordered_seqs, strict=False)
            ],
        )
        clustal_list.append("")

    return "\n".join(clustal_list)
