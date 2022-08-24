import numpy as np


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Stephen Ma"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def segments_from_diagonal(
    seq1, seq2, window, threshold, min_gap_length, diagonal
):  # pragma: no cover
    from cogent3.util.warning import discontinued

    discontinued(
        "function",
        "segments_from_diagonal",
        "2023.5",
        "replaced by much faster code in cogent3.align.compare",
    )

    assert window < 100
    scores = np.zeros(window)

    len1 = len(seq1)
    len2 = len(seq2)
    result = []
    was_high = 0
    score = 0
    i_lo = max(0, 0 - diagonal)
    i_hi = min(len1, len2 - diagonal)
    prior_end = 0
    for i in range(i_lo, i_hi):
        j = i + diagonal
        k = i % window
        score -= scores[k]
        scores[k] = seq1[i] == seq2[j]
        score += scores[k]
        if score >= threshold:
            if not was_high:
                start = max(i_lo, i - window)
                if min_gap_length and prior_end:
                    if start < prior_end + min_gap_length:
                        (start, jumped_end) = result.pop()
                was_high = 1
        else:
            if was_high:
                result.append((start, i))
                prior_end = i
                was_high = 0
    if was_high:
        result.append((start, i_hi))
    return result
