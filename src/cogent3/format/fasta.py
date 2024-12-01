"""Writer for FASTA sequence format"""

from collections.abc import Iterable
from typing import Optional


def _iter_in_block_size(series: str, block_size: int) -> Iterable[str]:
    """Yield chunks of a series of items"""
    for i in range(0, len(series), block_size):
        yield series[i : i + block_size]


def seqs_to_fasta(
    seqs: dict[str, str],
    block_size: int = 60,
    order: Optional[list[str]] = None,
) -> str:
    """Returns a Fasta string given an alignment.

    Parameters
    ----------
    seqs
        seq_name to sequence mapping
    block_size
        the sequence length to write to each line,
        by default 60
    order
        optional list of sequence names, which order to print in.
        Assumes complete and correct list of names. Defaults to
        iteration order of seqs.

    Returns
    -------
    The sequences in the Fasta format.
    """
    order = order or seqs
    result = []
    for name in order:
        result.append(f">{name}")
        seq = str(seqs[name])
        if len(seq) <= block_size:
            result.append(seq)
        else:
            result.extend(_iter_in_block_size(seq, block_size))
    if result:
        result.append("")
    return "\n".join(result)
