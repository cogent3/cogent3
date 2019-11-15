"""50x speedup for dotplots, but sequences must be strings and scoring is based on identity only
"""

version_info = (1, 3)
__version__ = "('2019', '11', '15', 'a')"

cdef int cmax(int a, int b):
    if a > b:
        return a
    else:
        return b
            
cdef int cmin(int a, int b):
    if a < b:
        return a
    else:
        return b
            
def segments_from_diagonal(
        char seq1[],
        char seq2[],
        int window,
        int threshold,
        int min_gap_length,
        int diagonal):
    """List of ((x1,y1), (x2,y2)) for diagonal line segments scoring >= threshold/window"""
        
    cdef int was_high, score, i, i_lo, i_hi, j, k, start, prior_end
    cdef int len1, len2
    cdef int scores[100]   # yuk - use malloc?
    assert window < 100
    len1 = len(seq1)
    len2 = len(seq2)
    result = []
    was_high = 0
    for i from 0 <= i < window:
        scores[i] = 0
    score = 0
    i_lo = cmax(0, 0-diagonal)
    i_hi = cmin(len1, len2-diagonal)
    prior_end = 0
    for i from i_lo <=  i < i_hi:
        j = i + diagonal
        k = i % window
        score -= scores[k]
        scores[k] = (seq1[i] == seq2[j])
        score += scores[k]
        if score >= threshold:
            if not was_high:
                start = cmax(i_lo, i - window)
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
