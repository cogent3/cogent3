#!/usr/bin/env python
"""Conversion of dynamic program results ("arrays of arrows") into gap vectors,
gapped sequences or Cogent Alignment objects"""

from cogent3.core.alignment import Aligned, Alignment
from cogent3.core.annotation import Map


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Rob Knight", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


def seq_traceback(s1, s2, aligned_positions, gap_value):
    """gapped sequences from state matrix and ending point
    gaps are signified by 'gap_value' inserted in the sequences.
    """
    seqs = [s1, s2]
    alignments = [[], []]

    for posn in aligned_positions:
        for (dimension, pos) in enumerate(posn):
            if pos is not None:
                c = seqs[dimension][pos]
            else:
                c = gap_value
            alignments[dimension].append(c)

    for dimension in [0, 1]:
        alignments[dimension].reverse()

    if isinstance(s1, str):
        alignments = ["".join(a) for a in alignments]

    return alignments


def gap_traceback(aligned_positions):
    """gap Vectors from state matrix and ending point"""
    consuming = [False, False]
    starts = [None, None]
    ends = [None, None]
    gap_vectors = [[], []]
    for (a, posn) in enumerate(aligned_positions):
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
    return (starts, ends, gap_vectors, a)


def map_traceback(aligned_positions):
    # using Map's to keep track of gaps for indel alignment
    (starts, ends, gap_vectors, alignment_len) = gap_traceback(aligned_positions)
    # print 'gv', gap_vectors
    maps = [Map(gv, parent_length=alignment_len).inverse() for gv in gap_vectors]
    return (starts, ends, maps)


def alignment_traceback(seqs, aligned_positions, word_length):
    """Alignment object from state matrix and ending point."""
    (starts, ends, maps) = map_traceback(aligned_positions)
    aligneds = []
    for (start, end, amap, (name, seq)) in zip(starts, ends, maps, seqs):
        gs = Aligned(amap * word_length, seq[start * word_length : end * word_length])
        aligneds.append((name, gs))
    return Alignment(moltype=None, data=aligneds)
