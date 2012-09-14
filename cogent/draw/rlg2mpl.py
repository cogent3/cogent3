"""ReportLab Graphics -> Matplotlib helpers"""

#linear, dotplot and dendrogram were originally written to target ReportLab Graphics rather
#than Matplotlib.  This module can be slowly boiled away as they become more matplotlib native.

from __future__ import division
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.text import Text
import matplotlib.patches as mpatches
import matplotlib.artist
import matplotlib.transforms
import matplotlib.colors
import numpy

from cogent.util.warning import discontinued

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def line_options(strokeColor='black', fillColor='black', strokeWidth=1):
    """Map from RLD line option names"""
    return dict(edgecolor=strokeColor, facecolor=fillColor, linewidth=strokeWidth)

def Line(x1, y1, x2, y2, **kw):
    """Acts like the RLG shape class of the same name"""
    path = Path([(x1, y1), (x2, y2)], [Path.MOVETO, Path.LINETO])
    return mpatches.PathPatch(path, **line_options(**kw))

def Rect(x, y, width, height, **kw):
    """Acts like the RLG shape class of the same name"""
    return mpatches.Rectangle((x,y),width,height,**line_options(**kw))

def Polygon(vertices, **kw):
    """Acts like the RLG shape class of the same name"""
    return mpatches.Polygon(vertices, **line_options(**kw))
    
def String(x, y, text, textAnchor='start', fontName=None, fontSize=10, 
        fillColor='black', rotation=None):
    """Acts like the RLG shape class of the same name"""
    fontname = fontName
    fontsize = fontSize
    color = fillColor
    ha = {'start':'left', 'middle':'center', 'end':'right'}[textAnchor]
    va = 'baseline'
    mpl_kw = dict((n,v) for (n,v) in locals().items() if n.islower())
    return Text(**mpl_kw)

class Group(matplotlib.artist.Artist):
    """Acts like the RLG shape class of the same name
    
    Groups elements together.  May apply a transform to its contents."""

    def __init__(self, *elements):
        """Initial lists of elements may be provided to allow
        compact definitions in literal Python code.  May or
        may not be useful."""
        matplotlib.artist.Artist.__init__(self)
        self.contents = []
        self.group_transform = matplotlib.transforms.Affine2D()
        self.outer_transform = matplotlib.transforms.TransformWrapper(
            self.get_transform())
        self.combined_transform = self.group_transform + self.outer_transform
        for elt in elements:
            self.add(elt)

    def set_figure(self, fig):
        matplotlib.artist.Artist.set_figure(self, fig)
        for c in self.contents:
            c.set_figure(fig)
    
    def set_clip_path(self, patch):
        matplotlib.artist.Artist.set_clip_path(self, patch)
        for c in self.contents:
            c.set_clip_path(patch)
          
    def set_transform(self, transform):
        matplotlib.artist.Artist.set_transform(self, transform)
        self.outer_transform.set(self.get_transform())
        
    def add(self, node, name=None):
        """Appends non-None child node to the 'contents' attribute. In addition,
        if a name is provided, it is subsequently accessible by name
        """
        # propagates properties down
        node.set_transform(node.get_transform() + self.combined_transform)
        self.contents.append(node)
    
    def rotate(self, theta):
        """Convenience to help you set transforms"""
        self.group_transform.rotate_deg(theta)

    def translate(self, dx, dy):
        """Convenience to help you set transforms"""
        self.group_transform.translate(dx, dy)

    def scale(self, sx, sy):
        """Convenience to help you set transforms"""
        self.group_transform.scale(sx, sy)

    def draw(self, renderer, *args, **kw):
        for c in self.contents:
            c.draw(renderer, *args, **kw)         


def figureLayout(width=None, height=None, margin=0.25, aspect=None, 
        default_aspect=0.75, useful_width=None, leftovers=False, **margins):
    """Width and height of a figure, plus a bounding box that nearly fills it, derived
    from defaults or provided margins.  All input figures are in inches."""
    left = margins.pop('left', 0) + margin
    right = margins.pop('right', 0) + margin
    top = margins.pop('top', 0) + margin
    bottom = margins.pop('bottom', 0) + margin
    default_width = 6  # use rcParams here?
    if useful_width:
        default_width = min(useful_width, default_width)
    width = width or default_width
    if aspect is not None:
        assert not height
        height = aspect * width
    else:
        height = height or default_aspect * width
    total_height = height + top + bottom
    total_width = width + left + right 
    posn = [left/total_width, bottom/total_height, width/total_width, height/total_height]
    if leftovers:
        return (total_width, total_height), posn, margins
    else:   
        assert not margins, margins.keys()
        return (total_width, total_height), posn

class Drawable(object):
    # Superclass for objects which can generate a matplotlib figure, in order 
    # to supply consistent and convenient showFigure() and drawToFile() 
    # methods.
    # Subclasses must provide .makeFigure() which will make use of 
    # _makeFigure() matplotlib.pyplot import done at runtime to give the 
    # user every chance to change the matplotlib backend first
    def _makeFigure(self, width, height, **kw):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(width,height), facecolor='white')
        return fig
                
    def drawFigure(self, title=None, **kw):
        """Draw the figure.  
        Extra arguments are forwarded to self.makeFigure()"""
        import matplotlib.pyplot as plt
        fig = self.makeFigure(**kw)
        if title is not None:
            fig.suptitle(title)
        plt.draw_if_interactive()
        
    def showFigure(self, title=None, **kw):
        """Make the figure and immediately pyplot.show() it.  
        Extra arguments are forwarded to self.makeFigure()"""
        self.drawFigure(title, **kw)
        import matplotlib.pyplot as plt
        plt.show()
        
    def drawToFile(self, fname, **kw):
        """Save in a file named 'fname'
        Extra arguments are forwarded to self.makeFigure() unless 
        they are valid for savefig()"""
        makefig_kw = {}
        savefig_kw = {}
        for (k,v) in kw.items():
            if k in ['dpi', 'facecolor', 'edgecolor', 'orientation',
                    'papertype', 'format', 'transparent']:
                savefig_kw[k] = v
            else:
                makefig_kw[k] = v    
        fig = self.makeFigure(**makefig_kw)
        fig.savefig(fname, **savefig_kw)
        
    def drawToPDF(self, filename, total_width=None, height=None, **kw):
        # Matches, as far as possible, old ReportLab version
        
        if total_width is not None:
            kw['width'] = total_width / 72
        kw2 = {}
        for (k,v) in kw.items():
            if k in ['wraps', 'border', 'withTrackLabelColumn']:
                discontinued('argument', "%s" % k, '1.6')
            else:
                kw2[k] = v
        kw2['format'] = 'pdf'
        if height:
            kw2['height'] = height / 72
        return self.drawToFile(filename, **kw2)


# For sequence feature styles:
# Matplotlib has fancy_box and fancy_arrow.  The code below is 
# similar, except that the two ends of the box have independent 
# styles: open, square, rounded, pointy, or blunt
        
class PathBuilder(object):
    """Path, made up of straight lines and bezier curves."""
    # Only used by the _End classes below

    def __init__(self):
        self.points = []
        self.operators = []

    def asPath(self):
        return Path(self.points, self.operators)

    def moveTo(self, x, y):
        self.points.append((x, y))
        self.operators.append(Path.MOVETO)

    def lineTo(self, x, y):
        self.points.append((x, y))
        self.operators.append(Path.LINETO)

    def curveTo(self, x1, y1, x2, y2, x3, y3):
        self.points.extend([(x1, y1), (x2, y2), (x3, y3)])
        self.operators.extend([Path.CURVE4]*3)

    def closePath(self):
        self.points.append((0.0, 0.0)) # ignored
        self.operators.append(Path.CLOSEPOLY)

class _End(object):
    
    def __init__(self, x_near, x_far, y_first, y_second, **kw):
        self.x_near = x_near
        self.x_far = x_far
        self.y_first = y_first
        self.y_second = y_second
        for (n, v) in kw.items():
            setattr(self, n, v)
    
    def moveToStart(self, path):
        path.moveTo(*self.startPoint())
    
    def drawToStart(self, path):
        path.lineTo(*self.startPoint())
    
    def finish(self, path):
        path.closePath()
    
    def startPoint(self):
        return (self.x_near, self.y_first)
        
    def __add__(self, oppo):
        p = PathBuilder()
        self.moveToStart(p)
        self.drawEnd(p)
        oppo.drawToStart(p)
        oppo.drawEnd(p)
        self.finish(p)
        return p.asPath()        
    
class Open(_End):
    def finish(self, path):
        self.drawToStart(path)
    
    def drawEnd(self, path):
        path.moveTo(self.x_near, self.y_second)
    
class Square(_End):
    def drawEnd(self, path):
        path.lineTo(self.x_near, self.y_second)
    
class Rounded(_End):
    def startPoint(self):
        return (self.x_near + self.dx, self.y_first)
    
    def drawEnd(self, path):
        path.curveTo(self.x_near, self.y_first, self.x_near, self.y_first,
                self.x_near, self.y_first + self.dy)
        path.lineTo(self.x_near, self.y_second - self.dy)
        path.curveTo(self.x_near, self.y_second, self.x_near, self.y_second,
                self.x_near + self.dx, self.y_second)    

class Pointy(_End):
    def _effective_dx(self):
        return max(abs(self.dx), abs(self.dy))*self.dx/abs(self.dx)
    
    def startPoint(self):
        return (self.x_near + self._effective_dx(), self.y_first)
    
    def drawEnd(self, path):
        head_start = self.x_near + self._effective_dx()
        middle = (self.y_first + self.y_second) / 2
        if self.blunt:
            for (x, y) in [
                    (head_start, self.y_first + self.dy),
                    (self.x_near, self.y_first),
                    (self.x_near, self.y_second),
                    (head_start, self.y_second - self.dy),
                    (head_start, self.y_second)]:
                path.lineTo(x, y)
        else:
            for (x, y) in [
                    (head_start, self.y_first + self.dy),
                    (self.x_near, middle),
                    (head_start, self.y_second - self.dy),
                    (head_start, self.y_second)]:
                path.lineTo(x, y)

def _sign(x):
    return x and x/abs(x)

def End(x1, x2, y1, y2, closed=True, rounded=False, pointy=False, blunt=False, 
        min_width=0.5, proportion_of_track=0.6):
    inwards = _sign(x2 - x1)
    span = max(abs(x2 - x1), min_width)
    thickness = abs(y1-y2)
    if pointy:
        head_size = min(thickness/20, span)
        head_size = max(head_size, min_width/2)
        height = thickness / proportion_of_track
        spare = (height - thickness) / 2
        end = Pointy(x1, x2, y1, y2,
                dx=inwards*head_size/2, dy=_sign(y1-y2)*spare,
                blunt=blunt)
    elif rounded:
        ry = thickness/4
        rx = min(span, ry)
        end = Rounded(x1, x2, y1, y2, dx=rx*inwards, dy=ry*_sign(y2-y1))
    elif not closed:
        end = Open(x1, x2, y1, y2)
    else:
        end = Square(x1, x2, y1, y2)
    return end

    
