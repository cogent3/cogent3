"""Writer for FASTA sequence format"""

from typing import Iterable, Optional

from cogent3.format.util import _AlignmentFormatter
from cogent3.util import warning as c3warn


def _iter_in_block_size(series: str, block_size: int) -> Iterable[str]:
    """Yield chunks of a series of items"""
    for i in range(0, len(series), block_size):
        yield series[i : i + block_size]


@c3warn.deprecated_args(
    "2024.9", reason="better name", old_new=[("alignment_dict", "seqs")]
)
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


@c3warn.deprecated_callable("2024.9", reason="better name", new="seqs_to_fasta")
def alignment_to_fasta(
    alignment_dict: dict[str, str],
    block_size: int = 60,
    order: Optional[list[str]] = None,
) -> str:  # pragma: no cover
    """use seqs_to_fasta() instead"""
    return seqs_to_fasta(seqs=alignment_dict, block_size=block_size, order=order)


@c3warn.deprecated_callable(
    "2024.9", reason="inefficient", is_discontinued=True, new="seqs_to_fasta"
)
class FastaFormatter(_AlignmentFormatter):  # pragma: no cover
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
