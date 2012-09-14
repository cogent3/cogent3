#!/usr/bin/env python
"""Parsers for the cigar format

Cigar stands for Compact Idiosyncratic Gapped Alignment Report and defines the sequence 
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
from cogent.core.location import LostSpan, Span, Map, _LostSpan
from cogent import DNA, LoadSeqs

__author__ = "Hua Ying"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Hua Ying"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Hua Ying"
__email__ = "hua.ying@anu.edu.au"
__status__ = "Production"

pattern = re.compile('([0-9]*)([DM])')

def map_to_cigar(map):
    """convert a Map into a cigar string"""
    cigar = ''
    for span in map.spans:
        if isinstance(span, Span):
            num_chars = span.End-span.Start
            char = 'M'
        else:
            num_chars = span.length
            char = 'D'
        if num_chars == 1:
            cigar += char
        else:
            cigar += str(num_chars)+char
    return cigar

def cigar_to_map(cigar_text):
    """convert cigar string into Map"""
    assert 'I' not in cigar_text
    spans, posn = [], 0
    for n, c in pattern.findall(cigar_text):
        if n:
            n = int(n)
        else:
            n = 1
            
        if c == 'M':
            spans.append(Span(posn, posn+n))
            posn += n
        else:
            spans.append(LostSpan(n))
    map = Map(spans = spans, parent_length = posn)
    return map

def aligned_from_cigar(cigar_text, seq, moltype=DNA):
    """returns an Aligned sequence from a cigar string, sequence and moltype"""
    if isinstance(seq, str):
        seq = moltype.makeSequence(seq)
    map = cigar_to_map(cigar_text)
    aligned_seq = seq.gappedByMap(map)
    return aligned_seq

def _slice_by_aln(map, left, right):
    slicemap = map[left:right]
    if hasattr(slicemap, 'Start'):
        location = [slicemap.Start, slicemap.End]
    else:
        location = []
    return slicemap, location

def _slice_by_seq(map, start, end):
    re_map = map.inverse()
    slicemap = re_map[start:end]
    aln_start, aln_end = slicemap.Start, slicemap.End
    new_map = map[aln_start:aln_end]
    return new_map, [aln_start, aln_end]

def _remap(map):
    start = map.Start
    if start == 0:
        new_map = map
        new_map.parent_length = map.End
    else:
        spans = []
        for span in map.spans:
            if span.lost:
                spans.append(span)
            else:
                span.Start = span.Start - start
                span.End = span.End - start
                length = span.End
                spans.append(span)
        new_map = Map(spans = spans, parent_length = length)
    return new_map

def slice_cigar(cigar_text, start, end, by_align=True):
    """slices a cigar string as an alignment"""
    map = cigar_to_map(cigar_text)
    if by_align:
        new_map, location = _slice_by_aln(map, start, end)
    else:
        new_map, location = _slice_by_seq(map, start, end)
    if hasattr(new_map, 'Start'):
        new_map = _remap(new_map)
    return new_map, location

def CigarParser(seqs, cigars, sliced = False, ref_seqname = None, start = None, end = None, moltype=DNA):
    """return an alignment from raw sequences and cigar strings
    if sliced, will return an alignment correspondent to ref sequence start to end
    
    Arguments:
        seqs - raw sequences as {seqname: seq}
        cigars - corresponding cigar text as {seqname: cigar_text}
        cigars and seqs should have the same seqnames
        MolType - optional default to DNA
    """
    data = {}
    if not sliced:
        for seqname in seqs.keys():
            aligned_seq = aligned_from_cigar(cigars[seqname], 
                                            seqs[seqname], moltype=moltype)
            data[seqname] = aligned_seq
    else:
        ref_aln_seq = aligned_from_cigar(cigars[ref_seqname], 
                                        seqs[ref_seqname], moltype=moltype)
        m, aln_loc = slice_cigar(cigars[ref_seqname], start, end, by_align = False)
        data[ref_seqname] = ref_aln_seq[aln_loc[0]:aln_loc[1]]
        for seqname in [seqname for seqname in seqs.keys() if seqname != ref_seqname]:
            m, seq_loc = slice_cigar(cigars[seqname], aln_loc[0], aln_loc[1])
            if seq_loc:
                seq = seqs[seqname]
                if isinstance(seq, str):
                    seq = moltype.makeSequence(seq)
                data[seqname] = seq[seq_loc[0]:seq_loc[1]].gappedByMap(m)
            else:
                data[seqname] = DNA.makeSequence('-'*(aln_loc[1] - aln_loc[0]))
    aln = LoadSeqs(data = data, aligned = True)
    return aln

