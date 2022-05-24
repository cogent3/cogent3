"""Parsers for the cigar format

Cigar stands for Compact Idiosyncratic gapped Alignment Report and defines the sequence
of matches/mismatches and deletions (or gaps). Cigar line is used in Ensembl database
for both multiple and pairwise genomic alignment.

for example, this cigar line 2MD3M2D2M will mean that the alignment contains 2 matches/
mismatches, 1 deletion (number 1 is omitted in order to save some spaces), 3 matches/
mismatches, 2 deletion and 2 matches/mismatches.

if the original sequence is: AACGCTT
   the cigar line is: 2MD3M2D2M
the aligned sequence will be:
    M M D M M M D D M M
    A A - C G C - - T T
"""

import re

from cogent3 import DNA, make_aligned_seqs
from cogent3.core.location import LostSpan, Map, Span


__author__ = "Hua Ying"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Hua Ying"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Hua Ying"
__email__ = "hua.ying@anu.edu.au"
__status__ = "Production"

_pattern = re.compile("([0-9]*)([DM])")


def map_to_cigar(map):
    """convert a Map into a cigar string"""
    cigar = ""
    for span in map.spans:
        if isinstance(span, Span):
            num_chars = span.end - span.start
            char = "M"
        else:
            num_chars = span.length
            char = "D"
        cigar += char if num_chars == 1 else str(num_chars) + char
    return cigar


def cigar_to_map(cigar_text):
    """convert cigar string into Map"""
    assert "I" not in cigar_text
    spans, posn = [], 0
    for n, c in _pattern.findall(cigar_text):
        n = int(n) if n else 1
        if c == "M":
            spans.append(Span(posn, posn + n))
            posn += n
        else:
            spans.append(LostSpan(n))
    return Map(spans=spans, parent_length=posn)


def aligned_from_cigar(cigar_text, seq, moltype=DNA):
    """returns an Aligned sequence from a cigar string, sequence and moltype"""
    if isinstance(seq, str):
        seq = moltype.make_seq(seq)
    map = cigar_to_map(cigar_text)
    return seq.gapped_by_map(map)


def _slice_by_aln(map, left, right):
    slicemap = map[left:right]
    location = [slicemap.start, slicemap.end] if hasattr(slicemap, "start") else []
    return slicemap, location


def _slice_by_seq(map, start, end):
    re_map = map.inverse()
    slicemap = re_map[start:end]
    aln_start, aln_end = slicemap.start, slicemap.end
    new_map = map[aln_start:aln_end]
    return new_map, [aln_start, aln_end]


def _remap(map):
    start = map.start
    if start == 0:
        new_map = map
        new_map.parent_length = map.end
    else:
        spans = []
        length = None
        for span in map.spans:
            if not span.lost:
                span.start = span.start - start
                span.end = span.end - start
                length = span.end
            spans.append(span)
        new_map = Map(spans=spans, parent_length=length)
    return new_map


def slice_cigar(cigar_text, start, end, by_align=True):
    """slices a cigar string as an alignment"""
    map = cigar_to_map(cigar_text)
    if by_align:
        new_map, location = _slice_by_aln(map, start, end)
    else:
        new_map, location = _slice_by_seq(map, start, end)
    if hasattr(new_map, "start"):
        new_map = _remap(new_map)
    return new_map, location


def CigarParser(
    seqs, cigars, sliced=False, ref_seqname=None, start=None, end=None, moltype=DNA
):
    """return an alignment from raw sequences and cigar strings
    if sliced, will return an alignment correspondent to ref sequence start to end

    Parameters
    ----------
        seqs - raw sequences as {seqname: seq}
        cigars - corresponding cigar text as {seqname: cigar_text}
        cigars and seqs should have the same seqnames
        moltype - optional default to DNA

    """
    data = {}
    if not sliced:
        for seqname in list(seqs.keys()):
            aligned_seq = aligned_from_cigar(
                cigars[seqname], seqs[seqname], moltype=moltype
            )
            data[seqname] = aligned_seq
    else:
        ref_aln_seq = aligned_from_cigar(
            cigars[ref_seqname], seqs[ref_seqname], moltype=moltype
        )
        m, aln_loc = slice_cigar(cigars[ref_seqname], start, end, by_align=False)
        data[ref_seqname] = ref_aln_seq[aln_loc[0] : aln_loc[1]]
        for seqname in [
            seqname for seqname in list(seqs.keys()) if seqname != ref_seqname
        ]:
            m, seq_loc = slice_cigar(cigars[seqname], aln_loc[0], aln_loc[1])
            if seq_loc:
                seq = seqs[seqname]
                if isinstance(seq, str):
                    seq = moltype.make_seq(seq)
                data[seqname] = seq[seq_loc[0] : seq_loc[1]].gapped_by_map(m)
            else:
                data[seqname] = DNA.make_seq("-" * (aln_loc[1] - aln_loc[0]))
    return make_aligned_seqs(data)
