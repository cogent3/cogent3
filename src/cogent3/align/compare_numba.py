import numba
import numpy as np

from numba import int64, njit
from numba.types import List
from numba.types.containers import Tuple


@njit(
    List(dtype=Tuple(types=(int64, int64)))(
        numba.typeof(b"seq"), numba.typeof(b"seq"), int64, int64, int64, int64
    ),
    cache=True,
)
def segments_from_diagonal(seq1, seq2, window, threshold, min_gap_length, diagonal):
    scores = np.zeros(window)
    assert window < 100

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
