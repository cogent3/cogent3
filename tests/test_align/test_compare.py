from collections import deque

import pytest

from cogent3 import make_seq
from cogent3.align.pycompare import (
    Kmer,
    MatchedSeqPaths,
    SeqKmers,
    _calc_seed_size,
    _extend_from_position,
    _extend_left,
    find_matched_paths,
    segment,
)


def _brute_force(
    seq1,
    seq2,
    window,
    threshold,
):
    """an exhaustive comparison of all windows between the two sequences"""
    mp = MatchedSeqPaths()
    for s1 in range(len(seq1) - window + 1):
        subseq1 = seq1[s1 : s1 + window]
        for s2 in range(len(seq2) - window + 1):
            subseq2 = seq2[s2 : s2 + window]
            total = sum(b1 == b2 for b1, b2 in zip(subseq1, subseq2, strict=False))
            if total >= threshold:
                mp[s2 - s1].append((segment(s1, s1 + window), segment(s2, s2 + window)))

    for y_intercept in mp.paths:
        if len(mp.paths[y_intercept]) == 1:
            continue
        merged = [mp.paths[y_intercept][0]]
        for x2, y2 in mp.paths[y_intercept][1:]:
            x1, y1 = merged[-1]
            if x1.overlap(x2):
                merged[-1] = (x1 | x2, y1 | y2)
            else:
                merged.append((x2, y2))
        mp.paths[y_intercept] = merged
    return mp


@pytest.fixture
def smallseq():
    return "ACCGGTT"


def test_find_matched_k_eq_1():
    s1 = make_seq(seq="TGATGTAAGGTAGTT", name="1")
    s2 = make_seq(seq="CTGGAAGGGT", name="2")
    expect = _brute_force(s1, s2, window=5, threshold=3)
    sk = SeqKmers(s1, k=1, canonical=set("ACGT"))
    got = find_matched_paths(seq_kmers=sk, seq1=s1, seq2=s2, window=5, threshold=3)
    assert got.paths == expect.paths


def test_calc_seed_size():
    """default seed size should not be larger than than threshold"""
    for i in range(2, 30):
        x = _calc_seed_size(30, i)
        assert x <= i


def test_calc_seed_size_singleton_window():
    x = _calc_seed_size(1, 1)
    assert x == 1


def test_segment():
    c = segment(2, 4)
    s, e = c
    assert (s, e) == (2, 4)
    c2 = segment(4, 6)
    assert not c.overlap(c2) and not c2.overlap(c)
    assert c.overlap(c)
    c3 = segment(1, 3)
    c4 = segment(3, 4)
    c5 = segment(3, 5)
    for other in (c3, c4, c5):
        assert c.overlap(other)

    c = segment(0, 2) - segment(4, 6)
    assert c == segment(2, 4)
    c = segment(4, 6) - segment(0, 2)
    assert c == segment(2, 4)
    assert c3 - c4 == segment(0, 0)
    false = segment(0, 0)
    assert not false


def test_segment_or():
    c2 = segment(3, 7)
    c3 = segment(4, 5)
    got = c2 | c3
    assert got == segment(3, 7)


def test_segment_merge():
    c1 = segment(1, 3)
    c2 = segment(3, 7)
    with pytest.raises(AssertionError):
        c1.merge(c2, strict=True)

    got = c1.merge(c2, strict=False)
    assert got == segment(1, 7)


def test_segment_nonzero():
    c = segment(1, 3)
    assert c
    c = segment(0, 0)
    assert not c


def test_segment_adjacent():
    c1 = segment(1, 3)
    c2 = segment(3, 7)
    c3 = segment(4, 5)

    assert c1.adjacent(c2) and c2.adjacent(c1)
    assert not c1.adjacent(c3) and not c3.adjacent(c1)
    assert not c2.adjacent(c3) and not c3.adjacent(c2)


def test_segment_rc():
    """should correctly transform for reverse complemented sequence"""
    seq = make_seq(seq="AACCCTTTTT", moltype="dna")
    c = segment(2, 5)
    assert seq[c.start : c.end] == "CCC"
    r = c.for_rc(len(seq))
    assert seq.rc()[r.start : r.end] == "GGG"


def test_get_segments():
    seq1 = "ACCGCTT"
    seq2 = "TTCCGCTTA"
    r = _extend_from_position(seq1, 1, seq2, 2, 2, 2, "ACGT")
    assert seq1[r[0].start : r[0].end] == seq2[r[1].start : r[1].end]


def test_kmer_one(smallseq):
    seq = make_seq(seq=smallseq, name="seq1")
    kmers = {Kmer(e, seq.name, i) for i, e in enumerate(seq.iter_kmers(k=2))}
    assert kmers == {smallseq[i : i + 2] for i in range(len(smallseq) - 1)}

    # with k==1, there's duplicates
    kmers = {}
    for i, kmer in enumerate(seq.iter_kmers(k=1)):
        if kmer in kmers:
            kmers[kmer].add_location(seq.name, i)
        else:
            kmer = Kmer(kmer, seq.name, i)
            kmers[kmer] = kmer

    assert len(kmers) == 4
    assert kmers["A"].indices[seq.name] == [0]
    assert kmers["C"].indices[seq.name] == [1, 2]
    assert kmers["G"].indices[seq.name] == [3, 4]
    assert kmers["T"].indices[seq.name] == [5, 6]

    # add a kmer from a different seq
    kmers["A"].add_location("seq2", 23)
    assert kmers["A"].indices["seq2"] == [23]
    with pytest.raises(NotImplementedError):
        kmers["A"].add_location("seq3", 23)


def test_seqkmers_1seq(smallseq):
    seq1 = make_seq(seq=smallseq, name="seq1", moltype="dna")
    sk = SeqKmers(seq1, 1, canonical=set(seq1.moltype))
    assert sk.num_seqs == 1
    # iter k-mers skips single occurrence in seq1
    for kmer in sk.iter_matching_kmers():
        assert kmer.kmer != "A"


@pytest.mark.parametrize("k,expect", [(1, 3), (2, 5), (7, 0)])
def test_seqkmers_1seq_degenerate(k, expect):
    # k-mers with degenerate characters are not stored
    seq1 = make_seq(seq="NCCGGTT", name="seq1", moltype="dna")
    sk = SeqKmers(seq1, k, canonical=set(seq1.moltype))
    assert len(sk.kmers) == expect
    assert "N" not in sk.kmers
    for kmer in sk.kmers:
        assert "N" not in kmer.kmer


def test_seqkmers_2seqs(smallseq):
    seq1 = make_seq(seq=smallseq, name="seq1", moltype="dna")
    seq2 = make_seq(seq=smallseq[1:], name="seq2", moltype="dna")
    sk = SeqKmers(seq1, 2, canonical=set(seq1.moltype))
    sk.add_seq(seq2)
    # the k-mer indices for "AC" should not have an entry for seq2, but should
    # have an entry for all others that is -1 the entry for seq1
    for kmer in sk.kmers:
        if kmer == "AC":
            assert seq2.name not in kmer.indices
        else:
            assert kmer.indices[seq1.name] == [e + 1 for e in kmer.indices[seq2.name]]


def test_seqkmers_1seq_iter_matched_indices():
    from itertools import product

    seq1 = make_seq(seq="", name="seq1", moltype="dna")
    sk = SeqKmers(seq1, 2, canonical=set("ACGT"))
    kmer = Kmer("AA", "seq1", 0)
    kmer.add_location("seq1", 1)
    sk.kmers[kmer] = kmer
    kmer = Kmer("CC", "seq1", 2)
    kmer.add_location("seq1", 4)
    sk.kmers[kmer] = kmer
    got = list(sk.iter_matched_indices())

    expect = list(product([0, 1], [0, 1])) + list(product([2, 4], [2, 4]))
    assert got == expect


def test_seqkmers_2seq_iter_matched_indices():
    from itertools import product

    seq1 = make_seq(seq="", name="seq1", moltype="dna")
    sk = SeqKmers(seq1, 2, canonical=set("ACGT"))
    sk.num_seqs += 1
    sk.other_name = "seq2"
    kmer = Kmer("AA", "seq1", 0)
    kmer.add_location("seq1", 1)
    kmer.add_location("seq2", 2)
    kmer.add_location("seq2", 4)
    sk.kmers[kmer] = kmer

    got = list(sk.iter_matched_indices())
    expect = list(product([0, 1], [2, 4]))
    assert got == expect


def test_seqkmers_2seqs_dropseq(smallseq):
    seq1 = make_seq(seq=smallseq, name="seq1", moltype="dna")
    seq2 = make_seq(seq=smallseq[1:], name="seq2", moltype="dna")
    sk = SeqKmers(seq1, 2, canonical=set(seq1.moltype))
    sk.add_seq(seq2)
    assert sk.ref_name == seq1.name
    assert sk.other_name == seq2.name
    assert sk.num_seqs == 2

    sk.drop_seq(seq2.name)
    assert sk.ref_name == seq1.name
    assert sk.other_name is None
    assert sk.num_seqs == 1


def test_dropseq_default(smallseq):
    seq1 = make_seq(seq=smallseq, name="seq1", moltype="dna")
    seq2 = make_seq(seq=smallseq[1:], name="seq2", moltype="dna")
    sk = SeqKmers(seq1, 2, canonical=set(seq1.moltype))
    sk.add_seq(seq2)
    sk.drop_seq()
    assert sk.other_name is None


def test_one_seq_dropseq(smallseq):
    seq1 = make_seq(seq=smallseq, name="seq1", moltype="dna")
    sk = SeqKmers(seq1, 2, canonical=set(seq1.moltype))
    assert sk.other_name is None
    assert sk.num_seqs == 1
    sk.drop_seq()
    assert sk.other_name is None
    assert sk.num_seqs == 1


def test_seqkmers_iter_matching_kmers(smallseq):
    seq1 = make_seq(seq=smallseq, name="seq1", moltype="dna")
    seq2 = make_seq(seq=smallseq[1:], name="seq2", moltype="dna")
    sk = SeqKmers(seq1, 2, canonical=set(seq1.moltype))
    sk.add_seq(seq2)
    for kmer in sk.iter_matching_kmers():
        assert kmer != "AC"
        assert len(kmer.indices) == 2


def test_matchedpaths():
    mp = MatchedSeqPaths()
    mp.append(segment(1, 3), segment(2, 5))
    assert mp.last_on_path(1, 1) == (segment(0, 0), segment(0, 0))
    assert mp.last_on_path(1, 2) == mp[1][0]
    mp.append(segment(10, 30), segment(11, 31))
    assert mp.last_on_path(1, 2) == mp[1][-1]
    assert mp.last_on_path(1, 2) != mp[1][0]


@pytest.fixture
def aseq1():
    return make_seq(seq="ATACGT", name="seq1", moltype="dna")


@pytest.fixture
def aseq2():
    return make_seq(seq="AGTTTACGTTTGACGTAA", name="seq2", moltype="dna")


def test_bruteforce(aseq1, aseq2):
    got = _brute_force(aseq1, aseq2, 3, 3)
    expect = MatchedSeqPaths()
    expect[3].append((segment(1, 6), segment(4, 9)))
    expect[10].append((segment(2, 6), segment(12, 16)))
    assert got.paths == expect.paths


def test_find_matched_paths_2seq(aseq1, aseq2):
    expect = _brute_force(aseq1, aseq2, 3, 3)
    sk = SeqKmers(aseq1, k=3, canonical=set("ACGT"))
    got = find_matched_paths(
        seq_kmers=sk,
        seq1=aseq1,
        seq2=aseq2,
        window=3,
        threshold=3,
    )
    assert got.paths == expect.paths


def test_matched_paths_rc():
    """adjusting y-coordinate for reverse complement"""
    path = MatchedSeqPaths()
    path[4].append((segment(start=6, end=14), segment(start=10, end=18)))
    xrc, yrc = path.get_coords(rc=True, length=24)
    # slope must be negative
    assert xrc[0] < xrc[1] and yrc[0] > yrc[1]


def test_matched_paths_min_gap():
    """correctly handle merging segments via gap position"""
    path = MatchedSeqPaths()
    path[4].extend(
        [
            (segment(start=6, end=8), segment(start=10, end=12)),
            (segment(start=10, end=12), segment(start=14, end=16)),
        ],
    )
    got = path._get_segments(min_gap=1)
    assert got == path.paths
    got2 = path._get_segments(min_gap=2)
    assert got2 != path.paths

    got3 = path._get_segments(min_gap=3)
    assert got3 == got2

    assert got2[4] == [(segment(6, 12), segment(10, 16))]

    # everything is hooked up
    trace = path.plotly_trace(min_gap=3)
    assert trace["x"] == [6, 12]
    assert trace["y"] == [10, 16]


@pytest.mark.parametrize("moltype", ["text", "rna", "bytes", "protein"])
def test_find_matched_paths_moltype(aseq1, aseq2, moltype):
    s1 = aseq1.to_moltype(moltype)
    s2 = aseq2.to_moltype(moltype)
    expect = _brute_force(s1, s2, 3, 3)
    sk = SeqKmers(aseq1, k=3, canonical="ACGT")
    got = find_matched_paths(
        seq_kmers=sk,
        seq1=s1,
        seq2=s2,
        window=3,
        threshold=3,
    )
    assert got.paths == expect.paths


def test_find_matched_with_rc():
    s = make_seq(seq="CACACCACTGCAGTCGGATAGACC", moltype="dna", name="s1")
    r = s.rc()
    r.name = "rc1"
    k = _calc_seed_size(4, 4)
    expect = _brute_force(s, r, 4, 4)
    sk = SeqKmers(s, k=k, canonical="ACGT")
    # note there's a bug here, it creates a y_intercept at -11
    got = find_matched_paths(seq_kmers=sk, seq1=s, seq2=r, window=4, threshold=4)
    assert got.paths == expect.paths
    x, y = got.paths[4][0]
    yrc = y.for_rc(len(s))
    # the rev complemented y-coord == x
    assert x == yrc


def test_find_matched_with_small_seed():
    s1 = make_seq(seq="CACACCACTGCAGTCGGATAGACC", moltype="dna", name="s1")
    s2 = make_seq(seq="GGTCTATCCGACTGCAGTGGTGTG", moltype="dna", name="s2")
    k = 2
    expect = _brute_force(s1, s2, 4, 4)
    sk = SeqKmers(s1, k=k, canonical="ACGT")
    got = find_matched_paths(seq_kmers=sk, seq1=s1, seq2=s2, window=4, threshold=4)
    assert got.paths == expect.paths


@pytest.mark.parametrize("w,t", [(4, 4), (3, 3)])
def test_find_matched_1seq(w, t):
    s = make_seq(seq="CACACCACTGCAGTCGGATAGACC", moltype="dna", name="s1")
    expect = _brute_force(s, s, w, t)
    sk = SeqKmers(s, k=w, canonical=set("ACGT"))
    got = find_matched_paths(seq_kmers=sk, seq1=s, window=w, threshold=t)
    assert got.paths == expect.paths


def test_plotly_trace(aseq1, aseq2):
    sk = SeqKmers(aseq1, k=3, canonical=set("ACGT"))
    got = find_matched_paths(
        seq_kmers=sk,
        seq1=aseq1,
        seq2=aseq2,
        window=3,
        threshold=3,
    )
    trace = got.plotly_trace()
    assert isinstance(trace, dict)
    assert trace["type"] == "scatter"
    assert len(trace["x"]) == len(trace["y"]) and len(trace["x"]) > 0


def _construct_matches(s1, s2, window):
    return [a == b for a, b in zip(s1, s2, strict=False)][:window]


@pytest.mark.parametrize("a,b", [(3, 5), (3, 0), (0, 4), (0, 0)])
def test_extend_left_no_shifts(a, b):
    window, threshold = 4, 3
    # preceeding base is mismatch, so starts should be unchanged
    #        | * mismatch within window
    s1 = "TTACGCAGCA"
    s2 = "TTTCCCGTAGCA"
    #          |
    matches = deque(_construct_matches(s1[a:], s2[b:], window), maxlen=4)
    total = sum(matches)
    expect = deque(tuple(matches), maxlen=4)

    # if either seq index is 0, just returns original values
    start1, start2, m = _extend_left(matches, s1, a, s2, b, 3, threshold, total)
    assert start1 == a and start2 == b and m == expect


_input_data = [("TTACGTTGCA", 1), ("ATCCGTCGCA", 2), ("TCCCGTCGCA", 3)]


@pytest.mark.parametrize("s1,diff", _input_data)
def test_extend_left_shifted(s1, diff):
    window, threshold = 4, 3
    #          ****
    s2 = "TTTCCCGTCGCA"
    a, b = 3, 5
    matches = deque(_construct_matches(s1[a:], s2[b:], window), maxlen=4)
    start1, start2, matches = _extend_left(
        matches,
        s1,
        a,
        s2,
        b,
        3,
        threshold,
        sum(matches),
    )
    assert start1 == a - diff
    assert start2 == b - diff
    assert (
        list(matches)
        == [a == b for a, b in zip(s1[a - diff :], s2[b - diff :], strict=False)][
            :window
        ]
    )


@pytest.mark.parametrize("left_limit", [3, 4, 5, 6])
def test_extend_left_truncated(left_limit):
    # cannot be beyond the start of a sequence
    window, threshold = 4, 3
    # preceeding base is mismatch, start shifted left by 1 due to
    #        |  * mismatch end of window
    s1 = "TTACGTTGCA"
    s2 = "TTTCCCGTCGCA"
    #          |
    a, b = 3, 5
    matches = deque(_construct_matches(s1[a:], s2[b:], window), maxlen=4)
    start1, start2, matches = _extend_left(
        matches,
        s1,
        a,
        s2,
        b,
        left_limit,
        threshold,
        sum(matches),
    )
    assert start1 == a - 1
    assert start2 == b - 1
    assert (
        list(matches)
        == [a == b for a, b in zip(s1[a - 1 :], s2[b - 1 :], strict=False)][:window]
    )
