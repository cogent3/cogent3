import pytest

import cogent3
from cogent3.parse.cigar import (
    aligned_from_cigar,
    cigar_to_alignment,
    cigar_to_map,
    map_to_cigar,
    slice_cigar,
)

DNA = cogent3.get_moltype("dna")


@pytest.fixture
def cigar_text():
    return "3D2M3D6MDM2D3MD"


@pytest.fixture
def aln_seq():
    return DNA.make_seq(seq="---AA---GCTTAG-A--CCT-")


@pytest.fixture
def aln_seq1():
    return DNA.make_seq(seq="CCAAAAAA---TAGT-GGC--G")


@pytest.fixture
def map_and_seq(aln_seq):
    return aln_seq.parse_out_gaps()


@pytest.fixture
def map1_and_seq1(aln_seq1):
    return aln_seq1.parse_out_gaps()


@pytest.fixture(params=[(1, 4), (0, 8), (7, 12), (0, 1), (3, 5)])
def start_end(request):
    return request.param


@pytest.fixture
def aln(aln_seq, aln_seq1):
    return cogent3.make_aligned_seqs(
        {"FAKE01": aln_seq, "FAKE02": aln_seq1},
        moltype="dna",
    )


@pytest.fixture
def cigars(cigar_text, map1_and_seq1):
    map1, _ = map1_and_seq1
    return {"FAKE01": cigar_text, "FAKE02": map_to_cigar(map1)}


@pytest.fixture
def seqs(map_and_seq, map1_and_seq1):
    _, seq = map_and_seq
    _, seq1 = map1_and_seq1
    return {"FAKE01": str(seq), "FAKE02": str(seq1)}


def test_map_to_cigar(map_and_seq, cigar_text):
    """convert a Map to cigar string"""
    map_, _ = map_and_seq
    assert map_to_cigar(map_) == cigar_text


def test_cigar_to_map(cigar_text, map_and_seq):
    """test generating a Map from cigar"""
    map_, _ = map_and_seq
    imap = cigar_to_map(cigar_text)
    assert str(imap) == str(map_)


def test_cigar_to_map_with_equals():
    """= symbol should advance sequence position like M"""
    imap = cigar_to_map("5=")
    assert imap.parent_length == 5
    assert imap.num_gaps == 0


def test_cigar_to_map_with_mismatch():
    """X symbol should advance sequence position like M"""
    imap = cigar_to_map("3X2M")
    assert imap.parent_length == 5
    assert imap.num_gaps == 0


def test_cigar_to_map_with_insertion():
    """I symbol should advance sequence position"""
    imap = cigar_to_map("3M2I3M")
    assert imap.parent_length == 8
    assert imap.num_gaps == 0


def test_cigar_to_map_with_skipped():
    """N symbol should create gaps like D"""
    imap = cigar_to_map("3M2N3M")
    assert imap.parent_length == 6
    assert imap.num_gaps == 1


def test_cigar_to_map_with_soft_clip():
    """S symbol should create gaps"""
    imap = cigar_to_map("2S4M2S")
    assert imap.parent_length == 4
    assert imap.num_gaps == 2


def test_cigar_to_map_with_hard_clip():
    """H symbol should create gaps"""
    imap = cigar_to_map("2H4M2H")
    assert imap.parent_length == 4
    assert imap.num_gaps == 2


def test_cigar_to_map_with_padding():
    """P symbol should create gaps"""
    imap = cigar_to_map("3M1P3M")
    assert imap.parent_length == 6
    assert imap.num_gaps == 1


def test_cigar_to_map_mixed_symbols():
    """Test CIGAR with multiple symbol types"""
    imap = cigar_to_map("2=3X2M1D2I")
    assert imap.parent_length == 9  # 2+3+2+2
    assert imap.num_gaps == 1  # 1D


def test_cigar_to_map_implicit_count():
    """CIGAR ops without numbers default to count of 1"""
    imap = cigar_to_map("MDMDM")
    assert imap.parent_length == 3  # M + M + M
    assert imap.num_gaps == 2  # D + D


def test_aligned_from_cigar(cigar_text, map_and_seq, aln_seq):
    """test generating aligned seq from cigar"""
    _, seq = map_and_seq
    aligned_seq = aligned_from_cigar(cigar_text, seq)
    assert aligned_seq == aln_seq


def test_slice_cigar(start_end, cigar_text, aln_seq, map_and_seq):
    """test slicing cigars"""
    _, seq = map_and_seq
    start, end = start_end
    # test by_align = True
    map1, loc1 = slice_cigar(cigar_text, start, end)
    ori1 = aln_seq[start:end]
    if loc1:
        slicealn1 = seq[loc1[0] : loc1[1]].gapped_by_map(map1)
        assert ori1 == slicealn1
    else:
        assert len(map1) == len(ori1)

    # test by_align = False
    map2, loc2 = slice_cigar(cigar_text, start, end, by_align=False)
    slicealn2 = seq[start:end].gapped_by_map(map2)
    ori2 = aln_seq[loc2[0] : loc2[1]]
    assert slicealn2 == ori2


def test_cigar_parser(seqs, cigars, aln):
    """test without slice"""
    result = cigar_to_alignment(seqs, cigars)
    assert result == aln
