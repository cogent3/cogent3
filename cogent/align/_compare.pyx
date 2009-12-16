"""50x speedup for dotplots, but sequences must be strings and scoring is based on identity only
"""

__version__ = "('1', '4')"

# Because this can be long-running check for interupts in the outer loop
cdef extern from "Python.h":
    PyErr_Occurred()
    int PyErr_CheckSignals()

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
            
def dotplot(
        char seq1[],
        char seq2[],
        int window,
        int threshold,
        min_gap_length,
        band,
        int progress):
    """List of ((x1,y1), (x2,y2)) for diagonal line segments scoring >= threshold/window"""
        
    cdef int was_high, score, i, i_lo, i_hi, diagonal, j, k, start, len1, len2, end
    cdef int scores[100]   # yuk - use malloc?
    assert window < 100
    len1 = len(seq1)
    len2 = len(seq2)
    result = []
    if band is None:
        band = max(len1, len2)
    for diagonal in range(-1*min(len1, band), min(len2, band)+1):
        if PyErr_CheckSignals():
            raise PyErr_Occurred()
        was_high = 0
        end = -1 * (min_gap_length+1)
        for i from 0 <= i < window:
            scores[i] = 0
        score = 0
        i_lo = cmax(0, 0-diagonal)
        i_hi = cmin(len1, len2-diagonal)
        if progress:
            if diagonal %100 == 0:
                print diagonal, (i_lo, i_lo+diagonal), (i_hi, i_hi+diagonal), len(result)
        for i from i_lo <=  i < i_hi:
            j = i + diagonal
            k = i % window
            score -= scores[k]
            scores[k] = (seq1[i] == seq2[j])
            score += scores[k]
            if score >= threshold:
                if not was_high:
                    start = cmax(i_lo, i - window)
                    was_high = 1
            else:
                if was_high:
                    if start - end  < min_gap_length:
                        ((start, forget1), (jumped_end, forget2)) = result.pop()
                    end = i - 1
                    result.append(((start, start+diagonal), (end, end+diagonal)))
                    was_high = 0
        if was_high:
            end = i - 1
            result.append(((start, start+diagonal), (end, end+diagonal)))
    return result
