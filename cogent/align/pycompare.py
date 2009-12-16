#!/usr/bin/env python
# Very slow.  See compare.pyx

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def dotplot(seq1, seq2, window, threshold, min_gap_length=0, band=None, progress=False):
    """A list of line segments covering the window-mers with identical matches > threshold
    
    Gaps of size less than min_gap will be hidden, which saves on line segments.
    if 'band' is not None then it limits the searched area
    """
    
    result = []
    if band is None:
        band = max(len(seq1), len(seq2))
    for diagonal in range(-min(len(seq1), band), min(len(seq2), band)+1):
        d_segments = []
        was_high = False
        scores = [0] * window
        score = 0
        (i_lo, i_hi) = max(0, -diagonal), min(len(seq1), len(seq2)-diagonal)
        if progress:
            if diagonal %100 == 0:
                print diagonal, (i_lo, i_lo+diagonal), (i_hi, i_hi+diagonal), len(result)
        for i in range(i_lo, i_hi):
            j = i + diagonal
            k = i % window
            score -= scores[k]
            scores[k] = seq1[i] == seq2[j]
            score += scores[k]
            if score >= threshold:
                if not was_high:
                    start = max(i_lo, i - window)
                    was_high = True
            else:
                if was_high:
                    end = i - 1
                    if d_segments and start-d_segments[-1][1] < min_gap_length:
                        (start, jumped_end) = d_segments.pop()
                    d_segments.append((start, end))
                    was_high = False
        if was_high:
            end = i - 1
            d_segments.append((start, end))
        
        for (start, end) in d_segments:
            result.append(((start, start+diagonal), (end, end+diagonal)))
    return result
