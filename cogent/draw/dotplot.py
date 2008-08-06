#!/usr/bin/env python
from reportlab.graphics import shapes, renderPDF
from reportlab.lib import colors

from cogent.draw import Scale, Display
from cogent.align.align import dotplot

import logging

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

LOG = logging.getLogger("cogent.draw")

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

def comparison_display(t1, t2, total_width, margin=10):
    """'Fat' annotated X and Y axes for a dotplot"""
    left = max(t1.label_width, t2.height)
    bottom = max(t2.label_width, t1.height)
    width = total_width-2*margin-left
    g1 = t1.asShape(width+t1.label_width)
    scale = 1.0*width/len(t1.base)
    height = scale * len(t2.base)
    g2 = t2.asShape(height, withTrackLabelColumn=False)
    
    g1.translate(left+margin-t1.label_width, bottom+margin-t1.height)
    
    g2.transform = (0.0, 1.0, -1.0, 0.0, left+margin, margin+bottom)
    
    D = shapes.Drawing(total_width, height+bottom+2*margin)
    D.add(g1)
    D.add(g2)
    
    D.add(shapes.Rect(margin+left, margin+bottom, width, height,
        strokeWidth=0.5, fillColor=None))
    return (D, left, bottom, scale)

class Display2D(object):
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.d1 = Display(self.seq1)
        self.d2 = Display(self.seq2)
        self.label_width = self.d1.label_width
        self._cache = {}
    
    def _calc_lines(self, scale, window, threshold, join_gaps):
        # Cache dotplot line segment coordinates as they can sometimes
        # be re-used at different resolutions, colours etc.
        (len1, len2) = (len(self.seq1), len(self.seq2))
        if threshold is None:
            universe = (len1-window) * (len2-window)
            acceptable_noise = 1.0*min(len1, len2)/window
            threshold = suitable_threshold(window, acceptable_noise/universe)
            LOG.info('require %s / %s bases' % (threshold, window))
            LOG.info('expect %s / %s matching' % (acceptable_noise, universe))
        
        if join_gaps:
            min_gap = int(1.0/scale)
        else:
            min_gap = 0
        
        key = (min_gap, window, threshold)
        if not self._cache.has_key(key):
            fwd = dotplot(self.seq1._seq, self.seq2._seq,
                    window, threshold, min_gap, None, False)
            if hasattr(self.seq1, "reversecomplement"):
                rev = dotplot(self.seq1.reversecomplement()._seq, self.seq2._seq,
                        window, threshold, min_gap, None, False)
                rev = [((len1-x1,y1), (len1-x2,y2)) for ((x1,y1), (x2,y2)) in rev]
            else:
                rev = []
            self._cache[key] = (fwd, rev)
        
        return self._cache[key]
    
    def asDrawing(self, total_width=600, margin=20, window=20, join_gaps=True):
        """Drawing of a dotplot with annotated axes"""
        (D, x, y, scale) = comparison_display(self.d1, self.d2, total_width=total_width, margin=margin)
        LOG.info('scale %s' %  scale)
        (fwd, rev) = self._calc_lines(scale, window, None, join_gaps)
        LOG.info('lines %s fwd, %s rev' % (len(fwd), len(rev)))
        g = shapes.Group()
        for lines, colour in [(fwd, colors.blue), (rev, colors.red)]:
            for ((x1,y1), (x2,y2)) in lines:
                g.add(shapes.Line(x1, y1, x2, y2,
                    strokeWidth=1.0/scale, strokeColor=colour))
        g.translate(margin+x, margin+y)
        g.scale(scale, scale)
        D.add(g)
        return D
    
    def drawToPDF(self, filename, *args, **kw):
        D = self.asDrawing(*args, **kw)
        renderPDF.drawToFile(D, filename)
    
