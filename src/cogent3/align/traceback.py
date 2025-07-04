"""Conversion of dynamic program results ("arrays of arrows") into gap vectors,
gapped sequences or Cogent Alignment objects"""

import typing

import numpy

from cogent3.core import alignment as c3_alignment
from cogent3.core.location import IndelMap

IntOrNone = int | None
IntListType = list[int]
CoordsListType = list[list[typing.Sequence[int]]]


def seq_traceback(s1, s2, aligned_positions, gap_value):
    """gapped sequences from state matrix and ending point
    gaps are signified by 'gap_value' inserted in the sequences.
    """
    seqs = [s1, s2]
    alignments = [[], []]

    for posn in aligned_positions:
        for dimension, pos in enumerate(posn):
            c = seqs[dimension][pos] if pos is not None else gap_value
            alignments[dimension].append(c)

    for dimension in [0, 1]:
        alignments[dimension].reverse()

    if isinstance(s1, str):
        alignments = ["".join(a) for a in alignments]

    return alignments


def gap_traceback(
    aligned_positions: list[list[IntOrNone, IntOrNone]],
) -> tuple[IntListType, IntListType, CoordsListType, int]:
    """gap Vectors from state matrix and ending point"""
    consuming = [False, False]
    starts = [None, None]
    ends = [None, None]
    gap_vectors = [[], []]
    for a, posn in enumerate(aligned_positions):
        for dimension in [0, 1]:
            delta = posn[dimension] is not None
            if delta:
                if starts[dimension] is None:
                    starts[dimension] = posn[dimension]
                ends[dimension] = posn[dimension] + 1
            if consuming[dimension] != delta:
                gap_vectors[dimension].append(a)
                consuming[dimension] = delta
    a += 1
    for dimension in [0, 1]:
        gv = gap_vectors[dimension]
        if consuming[dimension]:
            gv.append(a)
        gap_vectors[dimension] = [(gv[i], gv[i + 1]) for i in range(0, len(gv), 2)]
    return starts, ends, gap_vectors, a


def map_traceback(
    aligned_positions: list[list[IntOrNone, IntOrNone]],
) -> tuple[IntListType, IntListType, list[IndelMap]]:
    # using IndelMap's to keep track of gaps for indel alignment
    starts, ends, gap_vectors, alignment_len = gap_traceback(aligned_positions)
    maps = [
        IndelMap.from_aligned_segments(locations=gv, aligned_length=alignment_len)
        for gv in gap_vectors
    ]
    return starts, ends, maps


def alignment_traceback(seqs, aligned_positions, word_length) -> c3_alignment.Alignment:
    """Alignment object from state matrix and ending point."""
    (starts, ends, maps) = map_traceback(aligned_positions)
    ungapped_seqs = {}
    gaps = {}
    moltype = None
    for start, end, amap, (name, seq) in zip(starts, ends, maps, seqs, strict=False):
        if moltype is None:
            moltype = seq.moltype
        ungapped_seqs[name] = numpy.array(seq[start * word_length : end * word_length])
        gaps[name] = (amap * word_length).array

    asd = c3_alignment.AlignedSeqsData.from_seqs_and_gaps(
        seqs=ungapped_seqs, gaps=gaps, alphabet=moltype.most_degen_alphabet()
    )
    return c3_alignment.make_aligned_seqs(asd, moltype=moltype)
