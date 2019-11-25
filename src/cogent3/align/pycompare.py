#!/usr/bin/env python
# Very slow.  See compare.pyx


import cogent3.util.progress_display as UI

from cogent3.util.modules import ExpectedImportError, importVersionedModule


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def py_segments_from_diagonal(seq1, seq2, window, threshold, min_gap_length, diagonal):
    d_segments = []
    was_high = False
    scores = [0] * window
    score = 0
    (i_lo, i_hi) = max(0, -diagonal), min(len(seq1), len(seq2) - diagonal)
    for i in range(i_lo, i_hi):
        j = i + diagonal
        k = i % window
        score -= scores[k]
        scores[k] = seq1[i] == seq2[j]
        score += scores[k]
        if score >= threshold:
            if not was_high:
                start = max(i_lo, i - window)
                if d_segments and start - d_segments[-1][1] < min_gap_length:
                    (start, jumped_end) = d_segments.pop()
                was_high = True
        else:
            if was_high:
                d_segments.append((start, i))
                was_high = False
    if was_high:
        d_segments.append((start, i_hi))

    return d_segments


try:
    from . import _compare

    # _compare = importVersionedModule('_compare', globals(),
    # (1, 3), "slow Python dotplot")
    segments_from_diagonal = _compare.segments_from_diagonal
except ImportError:
    _compare = None
    segments_from_diagonal = py_segments_from_diagonal


@UI.display_wrap
def dotplot(seq1, seq2, window, threshold, min_gap_length=0, band=None, ui=None):
    """A list of line segments covering the window-mers with identical matches > threshold

    gaps of size less than min_gap will be hidden, which saves on line segments.
    if 'band' is not None then it limits the searched area
    """

    def one_diagonal(dia):
        segs = segments_from_diagonal(
            seq1, seq2, window, threshold, min_gap_length, dia
        )
        return [((start, start + dia), (end, end + dia)) for (start, end) in segs]

    if band is None:
        band = max(len(seq1), len(seq2))

    if isinstance(seq1, str):
        seq1 = seq1.encode("utf8")

    if isinstance(seq2, str):
        seq2 = seq2.encode("utf8")

    diagonals = list(range(-min(len(seq1), band), min(len(seq2), band) + 1))
    result = []
    for diag_segments in ui.imap(one_diagonal, diagonals, noun="offset"):
        result.extend(diag_segments)
    return result
