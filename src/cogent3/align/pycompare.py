from __future__ import annotations

from collections import deque
from dataclasses import InitVar, dataclass, field
from itertools import product
from typing import Generator, Optional

from cogent3.core.sequence import Sequence


@dataclass
class segment:
    """coordinates of a segment of a linear series"""

    start: int
    end: int

    def __contains__(self, index: int):
        return self.start <= index < self.end

    def __len__(self):
        return abs(self.end - self.start)

    def __sub__(self, other):
        a, b = (self, other) if self.start < other.start else (other, self)
        if a.overlap(b) or a.end == b.start:
            return segment(0, 0)

        return segment(a.end, b.start)

    def __nonzero__(self):
        return self.start != self.end

    def __getitem__(self, index: int):
        return (self.start, self.end)[index]

    def overlap(self, other):
        """whether coordinates overlap"""
        return (self.start < other.end <= self.end) or (
            self.start <= other.start < self.end
        )

    def adjacent(self, o):
        """whether coordinates are adjacent"""
        if not all([self, o]):
            return False
        return (o.start == self.end < o.end) or (self.start == o.end < self.end)

    def __or__(self, other):
        return self.merge(other, strict=True)

    def merge(self, other, strict: bool = True):
        """merge segments

        Parameters
        ----------
        other : segment
            instance to be merged
        strict
            if True, start
        Returns
        -------

        """
        if strict and not self.overlap(other):
            raise AssertionError(
                f"failed strict merge as {self} and {other} do not overlap"
            )
        return segment(min(self.start, other.start), max(self.end, other.end))

    def for_rc(self, length: int):
        """returns segment for reverse complement

        Parameters
        ----------
        length : int
            length of the sequence
        """
        return segment(length - self.end, length - self.start)


def _extend_left(
    matches: deque,
    seq1: str,
    idx1: int,
    seq2: str,
    idx2: int,
    left_limit: int,
    threshold: int,
    total: int,
) -> tuple[int, int, deque]:
    """extend a matching seed to the left

    Parameters
    ----------
    matches : deque
        series of bool values satisfying threshold
    seq1 : str
        first sequence
    idx1 : int
        initial position within seq1
    seq2 : str
        second sequence
    idx2 : int
        initial position within seq2
    left_limit : int
        positive integer indicating maximum number of left-wards
        positions to evaluate
    threshold : int
        minimum number of matches
    total : int
        total number of True values in matches

    Returns
    -------
    tuple[int, int, deque]
        The left-most position start indices for seq1/seq2 that satisfy
        threshold. The deque is the values corresponding to that position.
    """
    if idx1 == 0 or idx2 == 0:
        return idx1, idx2, matches

    left_limit = min(idx1, idx2, left_limit)

    matches.reverse()  # because we seek start backwards from seed
    offset = 0
    for offset in range(-1, -left_limit - 1, -1):
        prev = matches[0]
        current = seq1[idx1 + offset] == seq2[idx2 + offset]
        total += current - prev
        if total < threshold:
            # revert to previous position
            idx1 += offset + 1
            idx2 += offset + 1
            break

        matches.append(current)
    else:
        idx1 += offset
        idx2 += offset

    matches.reverse()
    return idx1, idx2, matches


def _extend_from_position(
    seq1, idx1, seq2, idx2, window, threshold, canonical, left_limit=0
) -> tuple[segment, segment]:
    """extend a matching seed

    Parameters
    ----------
    seq1 : str
        first sequence
    idx1 : int
        initial position within seq1
    seq2 : str
        second sequence
    idx2 : int
        initial position within seq2
    left_limit : int
        positive integer indicating maximum number of left-wards
        positions to evaluate
    threshold : int
        minimum number of matches

    Returns
    -------
    tuple[int, int, deque]
        The left-most position start indices for seq1/seq2 that satisfy
        threshold. The deque is the values corresponding to that position.
    """
    canonical = set(canonical)
    if len(seq1) - idx1 < window or len(seq2) - idx2 < window:
        return segment(0, 0), segment(0, 0)
    matched = [seq1[idx1 + i] == seq2[idx2 + i] for i in range(window)]
    total = sum(matched)
    if total < threshold:
        return segment(0, 0), segment(0, 0)

    matches = deque(
        matched,
        maxlen=window,
    )
    idx1, idx2, matches = _extend_left(
        matches, seq1, idx1, seq2, idx2, left_limit, threshold, total
    )

    total = sum(matches)
    offset = -1
    for offset, (s1, s2) in enumerate(
        zip(seq1[idx1 + window :], seq2[idx2 + window :])
    ):
        prev = matches[0]
        matches.append(s1 == s2 if {s1, s2} < canonical else 0)
        total += matches[window - 1] - prev
        if total < threshold:
            break
    else:
        offset += 1

    return segment(idx1, idx1 + window + offset), segment(idx2, idx2 + window + offset)


@dataclass
class Kmer:
    """individual k-mer"""

    kmer: str
    ref_name: str
    index: InitVar[int]
    indices: dict = field(init=False)

    def __post_init__(self, index):
        self.kmer = str(self.kmer)
        self.indices = {self.ref_name: [index]}

    def __hash__(self) -> int:
        return hash(self.kmer)

    def __eq__(self, __o: object) -> bool:
        return getattr(__o, "kmer", __o) == self.kmer

    def add_location(self, seq_name: str, index: int):
        """add a location of self for seq_name

        Parameters
        ----------
        seq_name : str
            name of sequence
        index : int
            position of this k-mer

        Raises
        -------
        NotImplementedError if number of sequence names exceeds 2
        """
        indices = self.indices.get(seq_name, [])
        indices.append(index)
        self.indices[seq_name] = indices
        if len(self.indices) > 2:
            raise NotImplementedError(
                f"Maximum of 2 sequences allowed, have {list(self.indices.keys())}"
            )


@dataclass
class SeqKmers:
    """representing a sequence as Kmer instances"""

    seq: InitVar[Sequence]
    k: int
    canonical: set

    kmers: dict = field(init=False)
    num_seqs: int = field(init=False)
    ref_name: str = field(init=False)
    other_name: Optional[str] = field(init=False)

    def __post_init__(self, seq):
        self.canonical = set(self.canonical)
        canonical = self.canonical
        self.ref_name = name = seq.name

        kmers = {}
        for i, kmer in enumerate(seq.iter_kmers(k=self.k)):
            if not set(kmer) <= canonical:
                # exclude k-mers with  non-canonical characters
                continue

            if kmer in kmers:
                kmers[kmer].add_location(name, i)
            else:
                kmer = Kmer(kmer, name, i)
                kmers[kmer] = kmer

        self.kmers = kmers
        self.num_seqs = 1
        self.other_name = None

    def add_seq(self, seq: Sequence) -> None:
        """transforms seq into k-mers, adds indices of matches

        Parameters
        ----------
        seq : Sequence
            has a different name to current sequence

        Notes
        -----
        k-mers unique to this sequence are ignored.
        """
        if seq.name == self.ref_name:
            raise ValueError(f"seq name {seq.name} matches ref seq name")

        self.other_name = name = seq.name

        for i, kmer in enumerate(seq.iter_kmers(self.k)):
            if kmer in self.kmers:
                self.kmers[kmer].add_location(name, i)

        self.num_seqs += 1

    def drop_seq(self, seq_name: Optional[str] = None) -> None:
        """removes other seq from all k-mers"""
        seq_name = seq_name if seq_name else self.other_name
        if seq_name is None:
            return

        for kmer in self.kmers:
            kmer.indices.pop(seq_name, None)

        self.other_name = None
        self.num_seqs = 1

    def _iter_indices_1_seq(self) -> Generator:
        for kmer in self.iter_matching_kmers():
            indices = kmer.indices[self.ref_name]
            yield from product(indices, indices)

    def _iter_indices_2_seq(self) -> Generator:
        for kmer in self.iter_matching_kmers():
            indices = list(kmer.indices.values())
            yield from product(*indices)

    def iter_matched_indices(self) -> Generator:
        """generate (seq1 index, seq2 index) for matched k-mer seeds

        Notes
        -----

        The indices correspond to the start of occurrences of exact k-mer
        matches. If only one sequence, the indices are for matches off
        the main diagonal.
        """
        return (
            self._iter_indices_1_seq()
            if self.num_seqs == 1
            else self._iter_indices_2_seq()
        )

    def _iter_kmers_1_seq(self) -> Generator:
        for kmer in self.kmers:
            if len(kmer.indices[self.ref_name]) == 1:
                continue
            yield kmer

    def _iter_kmers_2_seq(self) -> Generator:
        for kmer in self.kmers:
            if len(kmer.indices) == 1:
                continue
            yield kmer

    def iter_matching_kmers(self) -> Generator:
        """generator of Kmer instances

        Notes
        -----
        If two sequences recorded, yields k-mers with entries in both
        sequences.
        If a single sequence, yields k-mers with multiple occurrences in the
        primary sequence.
        """
        return (
            self._iter_kmers_1_seq() if self.num_seqs == 1 else self._iter_kmers_2_seq()
        )


@dataclass
class MatchedSeqPaths:
    """holder of paths between two sequences"""

    name: str = ""
    paths: dict = field(init=False)

    def __post_init__(self):
        self.paths = {}

    def __getitem__(self, y_intercept: int):
        self.paths[y_intercept] = self.paths.get(y_intercept, [])
        return self.paths[y_intercept]

    def last_on_path(self, refseq_pos: int, otherseq_pos: int):
        """returns the last coordinates on path"""
        y_intercept = otherseq_pos - refseq_pos
        if y_intercept in self.paths:
            return self[y_intercept][-1]

        # if no entry, corresponds to origin for refseq, y_intercept for other
        return segment(0, 0), segment(y_intercept, y_intercept)

    def append(self, seq1_segment: segment, seq2_segment: segment) -> None:
        """
        add segment to path

        Parameters
        ----------
        seq1_segment,  seq2_segment
            spans for refseq and otherseq

        Notes
        -----
        If a segment has 0 distance to the last segment, they are merged
        """
        y_intercept = seq2_segment.start - seq1_segment.start
        last_x, last_y = self.last_on_path(seq1_segment.start, seq2_segment.start)
        if seq1_segment.adjacent(last_x):
            self[y_intercept][-1] = (
                last_x.merge(seq1_segment, strict=False),
                last_y.merge(seq2_segment, strict=False),
            )
        else:
            self[y_intercept].append((seq1_segment, seq2_segment))

    def position_covered(self, refseq_pos: int, otherseq_pos: int) -> bool:
        """evaluates whether the position already evaluated for y_intercept"""
        evaluated = self.paths.get(otherseq_pos - refseq_pos, [])
        # I'm reversing the order here in the expectation that the most
        # recent evaluation is most likely to cover a span
        return any(
            refseq_pos in refseq_segment for refseq_segment, _ in evaluated[::-1]
        )

    def _get_segments(self, min_gap):
        if not min_gap:
            return self.paths

        paths = {}
        for y_intercept in self.paths:
            result = self.paths[y_intercept][:]
            if len(result) == 1:
                paths[y_intercept] = result
                continue
            result = result[:1]
            for x, y in self.paths[y_intercept][1:]:
                last_x, last_y = result[-1]
                if len(x - last_x) <= min_gap:
                    result[-1] = (
                        last_x.merge(x, strict=False),
                        last_y.merge(y, strict=False),
                    )
                else:
                    result.append((x, y))
            paths[y_intercept] = result
        return paths

    def get_coords(
        self, rc: bool = False, length: Optional[int] = None, min_gap: int = 0
    ):
        """returns x, y coordinates for plotting

        Parameters
        ----------
        rc
            whether coordinates are for a reverse complemented sequence
        length
            length of the sequence
        min_gap
            connect disjoint segments separated by <= min_gap
        """
        paths = self._get_segments(min_gap)
        x = []
        y = []
        for y_intercept in paths.values():
            for a, b in y_intercept:
                if rc:
                    # adjust coordinate and reverse
                    s, e = b.for_rc(length)
                    b = segment(e, s)
                x.extend(a)
                x.append(None)
                y.extend(b)
                y.append(None)
        return x[:-1], y[:-1]  # discard last None

    def plotly_trace(
        self, rc: bool = False, length: Optional[int] = None, min_gap: int = 0
    ):
        x, y = self.get_coords(rc=rc, length=length, min_gap=min_gap)
        trace = {
            "type": "scatter",
            "x": x,
            "y": y,
            "mode": "lines",
        }
        if self.name:
            trace["name"] = self.name
            trace["showlegend"] = True
        return trace


def _calc_seed_size(w: int, t: int, min_val: int = 5) -> int:
    """computes k-mer size

    Parameters
    ----------
    w : int
        window size
    t : int
        threshold for minimum number of matches
    min_val : int
        minimum k-mer size

    Returns
    -------
    int
        k-mer size

    Notes
    -----
    This seeks to return the maximum value of k that would guarantee the
    detection of a w sized segment with t matches. Evenly spaced mismatches
    represent the hardest case. We set a min_val for k to ensure performance
    is maintained at the cost of some loss of sensitivity.
    """
    assert 0 < t <= w, f"threshold={t} > window size={w}"
    d = w - t
    if w == t:
        return w
    elif d > t:
        return t
    k = t // d
    return max(k, min_val)


def find_matched_paths(
    *,
    seq_kmers: SeqKmers,
    seq1: Sequence,
    seq2: Optional[Sequence] = None,
    window: int = 20,
    threshold: int = 17,
) -> MatchedSeqPaths:
    """determine all matches between seq1 and seq2

    Parameters
    ----------
    seq_kmers : SeqKmers
        instance with seq1 converted to k-mers
    seq1, seq2 : Sequence
        cogent3 Sequence instances. If seq2 not provided, compares seq1 to
        itself.
    window : int
        size of sequence segment to be compared
    threshold : int
        Minimum number of positions that must be equal within the window.
        Less than or equal to window.

    Returns
    -------
    MatchedSeqPaths
    """
    k = seq_kmers.k
    delta = (window // k) * (k - 1) + (window - threshold)
    canonical = set(seq1.moltype)
    paths = MatchedSeqPaths()

    if seq2 is None:
        seq2 = seq1
        paths[0].append((segment(0, len(seq1)), segment(0, len(seq1))))
    else:
        seq_kmers.add_seq(seq2)

    # for every matched k-mer, we determine the y-intercept.
    for seq1_idx, seq2_idx in seq_kmers.iter_matched_indices():
        if paths.position_covered(seq1_idx, seq2_idx):
            # seq1_idx already evaluated for y_intercept, skip
            continue

        # it hasn't been evaluated, compute the left most limit to start from
        s1_coord, _ = paths.last_on_path(seq1_idx, seq2_idx)
        left_limit = min(abs(s1_coord.end - seq1_idx), delta)
        ref_coord, other_coord = _extend_from_position(
            str(seq1),
            seq1_idx,
            str(seq2),
            seq2_idx,
            window=window,
            threshold=threshold,
            canonical=canonical,
            left_limit=left_limit,
        )

        if not ref_coord:
            continue

        paths.append(ref_coord, other_coord)

    for y_intercept in paths.paths:
        if len(paths.paths[y_intercept]) == 1:
            continue
        merged = [paths.paths[y_intercept][0]]
        for x2, y2 in paths.paths[y_intercept][1:]:
            x1, y1 = merged[-1]
            if x1.overlap(x2):
                merged[-1] = (x1 | x2, y1 | y2)
            else:
                merged.append((x2, y2))
        paths.paths[y_intercept] = merged

    return paths
