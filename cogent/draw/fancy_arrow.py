#!/usr/bin/env python
"""Extension of matplotlib's arrows to give additional control.

Contributed this back to the matplotlib developers, so may be obsolete.
"""
from matplotlib.patches import Polygon
from numpy import dot, array, sqrt, concatenate

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class FancyArrow(Polygon):
    """Like Arrow, but lets you set head width and head height independently."""

    def __init__(self, x, y, dx, dy, width=0.001, length_includes_head=False, \
        head_width=None, head_length=None, shape='full', overhang=0, \
        head_starts_at_zero=False,**kwargs):
        """Returns a new Arrow.

        length_includes_head: True if head is counted in calculating the length.
        
        shape: ['full', 'left', 'right']
        
        overhang: distance that the arrow is swept back (0 overhang means
        triangular shape).

        head_starts_at_zero: if True, the head starts being drawn at coordinate
        0 instead of ending at coordinate 0.
        """
        if head_width is None:
            head_width = 3 * width
        if head_length is None:
            head_length = 1.5 * head_width

        distance = sqrt(dx**2 + dy**2)
        if length_includes_head:
            length=distance
        else:
            length=distance+head_length
        if not length:
            verts = [] #display nothing if empty
        else:
            #start by drawing horizontal arrow, point at (0,0)
            hw, hl, hs, lw = head_width, head_length, overhang, width
            left_half_arrow = array([
                [0.0,0.0],                  #tip
                [-hl, -hw/2.0],             #leftmost
                [-hl*(1-hs), -lw/2.0], #meets stem
                [-length, -lw/2.0],          #bottom left
                [-length, 0],
            ])
            #if we're not including the head, shift up by head length
            if not length_includes_head:
                left_half_arrow += [head_length, 0]
            #if the head starts at 0, shift up by another head length
            if head_starts_at_zero:
                left_half_arrow += [head_length/2.0, 0]
            #figure out the shape, and complete accordingly
            if shape == 'left':
                coords = left_half_arrow
            else:
                right_half_arrow = left_half_arrow*[1,-1]
                if shape == 'right':
                    coords = right_half_arrow
                elif shape == 'full':
                    coords=concatenate([left_half_arrow,right_half_arrow[::-1]])
                else:
                    raise ValueError, "Got unknown shape: %s" % shape
            cx = float(dx)/distance
            sx = float(dy)/distance
            M = array([[cx, sx],[-sx,cx]])
            verts = dot(coords, M) + (x+dx, y+dy)
        
        Polygon.__init__(self, map(tuple, verts), **kwargs)

def arrow(axis, x, y, dx, dy, **kwargs):
    """Draws arrow on specified axis from (x,y) to (x+dx,y+dy)."""
    a = FancyArrow(x, y, dx, dy, **kwargs)
    axis.add_artist(a)
    return a
