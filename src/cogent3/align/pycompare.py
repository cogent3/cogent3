#!/usr/bin/env python
# Very slow.  See compare.pyx


import cogent3.util.progress_display as UI

from . import compare_numba


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


segments_from_diagonal = compare_numba.segments_from_diagonal


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
