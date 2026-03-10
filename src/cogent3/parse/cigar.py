import re
import typing

import numpy

import cogent3
from cogent3.core.location import IndelMap

if typing.TYPE_CHECKING:
    from cogent3.core.alignment import Alignment

_pattern = re.compile("([0-9]*)([M=XIDNSHP])")


def map_to_cigar(imap: IndelMap) -> str:
    """convert a IndelMap into a cigar string"""
    cigar = []
    if imap.num_gaps == 0:
        return f"{imap.parent_length}M"

    last_cumsum = 0
    last_seq_pos = 0
    for seq_pos, cum_gaps in zip(imap.gap_pos, imap.cum_gap_lengths, strict=False):
        gap_len = cum_gaps - last_cumsum
        last_cumsum = cum_gaps
        if seq_span := seq_pos - last_seq_pos:
            num_chars = "" if seq_span == 1 else seq_span
            cigar.append(f"{num_chars}M")

        num_gap = "" if gap_len == 1 else gap_len
        cigar.append(f"{num_gap}D")
        last_seq_pos = seq_pos

    if last_seq_pos < imap.parent_length:
        num_chars = (
            ""
            if imap.parent_length - last_seq_pos == 1
            else imap.parent_length - last_seq_pos
        )
        cigar.append(f"{num_chars}M")
    return "".join(cigar)


def cigar_to_map(cigar_text: str) -> IndelMap:
    """convert cigar string into IndelMap"""
    # TODO handle soft/hard clipping

    gpos = []
    glengths = []
    seq_pos = 0

    for n, c in _pattern.findall(cigar_text):
        n = int(n) if n else 1
        if c in {"M", "I", "=", "X"}:
            seq_pos += n
        else:
            gpos.append(seq_pos)
            glengths.append(n)

    gpos = numpy.array(gpos, dtype=numpy.int32)
    glengths = numpy.array(glengths, dtype=numpy.int32)

    return IndelMap(
        gap_pos=gpos,
        gap_lengths=glengths,
        parent_length=seq_pos,
    )


def aligned_from_cigar(cigar_text: str, seq: str) -> str:
    """returns an Aligned sequence from a cigar string, sequence and moltype"""
    seq = str(seq)
    imap = cigar_to_map(cigar_text)
    pos_length = imap.get_gap_coordinates()
    if not pos_length:
        return seq

    new_seq: list[str] = []
    last_pos = 0
    for pos, length in pos_length:
        new_seq.extend((seq[last_pos:pos], "-" * length))
        last_pos = pos
    new_seq.append(seq[last_pos : imap.parent_length])
    return "".join(new_seq)


def slice_cigar(
    cigar_text: str, start: int, end: int, by_align: bool = True
) -> tuple[IndelMap, list[int]]:
    """slices a cigar string as an alignment"""
    imap = cigar_to_map(cigar_text)
    if by_align:
        new_map = imap[start:end]
        location = [
            imap.get_seq_index(start),
            imap.get_seq_index(end),
        ]
    else:
        location = [
            imap.get_align_index(start),
            imap.get_align_index(end, slice_stop=True),
        ]
        new_map = imap[location[0] : location[1]]
    return new_map, location


def cigar_to_alignment(
    seqs: dict[str, str],
    cigars: dict[str, str],
    moltype: str = "dna",
) -> "Alignment":
    """return an alignment from raw sequences and cigar strings

    Parameters
    ----------
    seqs
        raw sequences as {seqname: seq}
    cigars
        corresponding cigar text as {seqname: cigar_text}
        cigars and seqs should have the same seqnames
    moltype
        default to 'dna'
    """
    moltype = cogent3.get_moltype(moltype)
    data = {}
    for seqname in list(seqs.keys()):
        aligned_seq = aligned_from_cigar(
            cigars[seqname],
            seqs[seqname],
        )
        data[seqname] = aligned_seq
    return cogent3.make_aligned_seqs(data, moltype=moltype)
