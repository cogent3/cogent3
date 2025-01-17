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

import cogent3
from cogent3.core.location import IndelMap, LostSpan, Span

_pattern = re.compile("([0-9]*)([DM])")


def map_to_cigar(map):
    """convert a IndelMap into a cigar string"""
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
    """convert cigar string into IndelMap"""
    assert "I" not in cigar_text
    spans, posn = [], 0
    for n, c in _pattern.findall(cigar_text):
        n = int(n) if n else 1
        if c == "M":
            spans.append(Span(posn, posn + n))
            posn += n
        else:
            spans.append(LostSpan(n))
    return IndelMap.from_spans(spans=spans, parent_length=posn)


def aligned_from_cigar(cigar_text, seq, moltype="dna"):
    """returns an Aligned sequence from a cigar string, sequence and moltype"""
    moltype = cogent3.get_moltype(moltype)
    if isinstance(seq, str):
        seq = moltype.make_seq(seq=seq)
    map = cigar_to_map(cigar_text)
    return seq.gapped_by_map(map)


def _slice_by_aln(map, left, right):
    slicemap = map[left:right]
    location = [
        map.get_seq_index(left),
        map.get_seq_index(right),
    ]
    return slicemap, location


def _slice_by_seq(map, start, end):
    # start, end are in sequence coords
    aln_start = map.get_align_index(start)
    aln_end = map.get_align_index(end, slice_stop=True)
    return map[aln_start:aln_end], [aln_start, aln_end]


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
        new_map = IndelMap(spans=spans, parent_length=length)
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
    seqs,
    cigars,
    sliced=False,
    ref_seqname=None,
    start=None,
    end=None,
    moltype="dna",
):
    """return an alignment from raw sequences and cigar strings
    if sliced, will return an alignment correspondent to ref sequence start to end

    Parameters
    ----------
        seqs - raw sequences as {seqname: seq}
        cigars - corresponding cigar text as {seqname: cigar_text}
        cigars and seqs should have the same seqnames
        moltype - optional default to 'dna'

    """
    moltype = cogent3.get_moltype(moltype)
    data = {}
    if not sliced:
        for seqname in list(seqs.keys()):
            aligned_seq = aligned_from_cigar(
                cigars[seqname],
                seqs[seqname],
                moltype=moltype,
            )
            data[seqname] = aligned_seq
    else:
        ref_aln_seq = aligned_from_cigar(
            cigars[ref_seqname],
            seqs[ref_seqname],
            moltype=moltype,
        )
        m, aln_loc = slice_cigar(cigars[ref_seqname], start, end, by_align=False)
        data[ref_seqname] = ref_aln_seq[aln_loc[0] : aln_loc[1]]
        for seqname in [
            seqname for seqname in list(seqs.keys()) if seqname != ref_seqname
        ]:
            m, seq_loc = slice_cigar(cigars[seqname], aln_loc[0], aln_loc[1])
            if seq_loc:
                start, stop = seq_loc
                seq = seqs[seqname][start:stop]
                if isinstance(seq, str):
                    seq = moltype.make_seq(seq=seq)
                data[seqname] = seq.gapped_by_map(m)
            else:
                data[seqname] = moltype.make_seq("-" * (aln_loc[1] - aln_loc[0]))
    return cogent3.make_aligned_seqs(data, moltype=moltype)
