#!/usr/bin/env python
from __future__ import division

from matplotlib.path import Path
from matplotlib.patches import PathPatch

from cogent.util.warning import discontinued
from cogent.draw.linear import Display
from cogent.draw.rlg2mpl import Drawable, figureLayout
from cogent.align.align import dotplot

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def suitable_threshold(window, desired_probability):
    """Use cumulative binomial distribution to find the number of identical
    bases which we expect a nucleotide window-mer to have with the desired
    probability"""
    cumulative_p = 0.0
    for matches in range(window, 0, -1):
        mismatches = window - matches
        p = 0.75 ** mismatches
        for i in range(matches, 0, -1):   # n
            p *= (i + mismatches)
            p /= i
            p *= 0.25
        cumulative_p += p
        if cumulative_p > desired_probability:
            break
    return matches

def _reinchify(figsize, posn, *args):
    (fw, fh) = figsize
    (x,y,w,h) = posn
    return [fw*x, fh*y, fw*w, fh*h]
    
def comparison_display(seq1, seq2, left=.5, bottom=.5, **kw):
    """'Fat' annotated X and Y axes for a dotplot
    
    Returns a matplotlib axes object placed and scaled ready for plotting 
    a sequence vs sequence comparison between the sequences (or alignments) 
    seq1 and seq2, which are also displayed. The aspect ratio will depend on 
    the sequence lengths as the sequences are drawn to the same scale"""
    
    import matplotlib.pyplot as plt
    
    (x1, y1, w1, h1) = _reinchify(*seq1.figureLayout(
        labeled=True, bottom=bottom, margin=0))
    (x2, y2, w2, h2) = _reinchify(*seq2.figureLayout(
        labeled=False, bottom=left, margin=0))
    
    # equalize points-per-base scales to get aspect ratio 1.0
    ipb = min(w1/len(seq1), w2/len(seq2))
    (w1, w2) = ipb*len(seq1), ipb*len(seq2)
    
    # Figure with correct aspect
    # Indent enough for labels and/or vertical display
    (w,h), posn = figureLayout(width=w1, height=w2,
        left=max(x1,y2+h2), bottom=y1+h1, **kw)
    fig = plt.figure(figsize=(w,h), facecolor='white')
    
    fw = fig.get_figwidth()
    fh = fig.get_figheight()
    # 2 sequence display axes
    x = seq1.asAxes(fig, [posn[0], posn[1]-h1/fh, posn[2], h1/fh])
    y = seq2.asAxes(fig, [posn[0]-h2/fw, posn[1], h2/fw, posn[3]], 
        vertical=True, labeled=False)
    
    # and 1 dotplot axes
    d = fig.add_axes(posn, sharex=x, sharey=y)
    d.xaxis.set_visible(False)
    d.yaxis.set_visible(False)
    return d

class Display2D(Drawable):
    def __init__(self, seq1, seq2, **kw):
        if not isinstance(seq1, Display):
            seq1 = Display(seq1, **kw)
        if not isinstance(seq2, Display):
            seq2 = Display(seq2, **kw)
        self.seq1 = seq1.base
        self.seq1d = seq1
        self.seq2 = seq2.base
        self.seq2d = seq2
        self._cache = {}
        # Check inputs are sufficiently sequence-like
        assert len(self.seq1) == len(str(self.seq1))
        assert len(self.seq2) == len(str(self.seq2))
    
    def _calc_lines(self, window, threshold, min_gap):
        # Cache dotplot line segment coordinates as they can sometimes
        # be re-used at different resolutions, colours etc.
        (len1, len2) = (len(self.seq1), len(self.seq2))
        if threshold is None:
            universe = (len1-window) * (len2-window)
            acceptable_noise = min(len1, len2) / window
            threshold = suitable_threshold(window, acceptable_noise/universe)
            # print 'require %s / %s bases' % (threshold, window)
            # print 'expect %s / %s matching' % (acceptable_noise, universe)
        
        key = (min_gap, window, threshold)
        if not self._cache.has_key(key):
            fwd = dotplot(str(self.seq1), str(self.seq2),
                    window, threshold, min_gap, None)
            if hasattr(self.seq1, "reversecomplement"):
                rev = dotplot(str(self.seq1.reversecomplement()), 
                        str(self.seq2), window, threshold, min_gap, None)
                rev = [((len1-x1,y1),(len1-x2,y2)) for ((x1,y1),(x2,y2)) in rev]
            else:
                rev = []
            self._cache[key] = (fwd, rev)
        
        return self._cache[key]
                
    def makeFigure(self, window=20, join_gaps=None, min_gap=0, **kw):
        """Drawing of a line segment based dotplot with annotated axes"""
        # hard to pick min_gap without knowing pixels per base, and
        # matplotlib is reasonably fast anyway, so:
        if join_gaps is not None:
            discontinued('argument', 'join_gaps', '1.6')
        ax = comparison_display(self.seq1d, self.seq2d, **kw)
        (fwd, rev) = self._calc_lines(window, None, min_gap)
        for (lines, colour) in [(fwd, 'blue'), (rev, 'red')]:
            vertices = []
            for segment in lines:
                vertices.extend(segment)
            if vertices:
                ops = [Path.MOVETO, Path.LINETO] * (len(vertices)//2)
                path = Path(vertices, ops)
                patch = PathPatch(path, edgecolor=colour, fill=False)
                ax.add_patch(patch)
        return ax.get_figure()
    
    def simplerMakeFigure(self):
        """Drawing of a matrix style dotplot with annotated axes"""
        import numpy
        ax = comparison_display(self.seq1, self.seq2)
        alphabet = self.seq1.MolType.Alphabet
        seq1 = alphabet.toIndices(self.seq1)
        seq2 = alphabet.toIndices(self.seq2)
        ax.pcolorfast(numpy.equal.outer(seq2, seq1))
        return ax.get_figure()
        
        
        
        
