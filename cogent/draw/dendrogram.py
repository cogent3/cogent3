#!/usr/bin/env python
"""Drawing trees.

Draws horizontal trees where the vertical spacing between taxa is
constant.  Since dy is fixed dendrograms can be either:
 - square:     dx = distance
 - not square: dx = max(0, sqrt(distance**2 - dy**2))

Also draws basic unrooted trees.

For drawing trees use either:
 - SquareDendrogram
 - StraightDendrogram
 - ContemporaneousDendrogram
 - ContemporaneousStraightDendrogram
 - ShelvedDendrogram
 - UnrootedDendrogram
"""

# Future:
#  - font styles
#  - orientation switch
# Layout gets more complicated for rooted tree styles if dy is allowed to vary,
# and constant-y is suitable for placing alongside a sequence alignment anyway.
from cogent.core.tree import TreeNode
import rlg2mpl
import matplotlib.colors
from matplotlib.patches import PathPatch, Polygon
from matplotlib.path import Path
from matplotlib.text import Text
import numpy

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Zongzhi Liu", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

to_rgb = matplotlib.colors.colorConverter.to_rgb

def _sign(x):
    """Returns True if x is positive, False otherwise."""
    return x and x/abs(x)

def _first_non_none(values):
    for item in values:
        if item is not None:
            return item
            
def SimpleColormap(color0, color1, name=None):
    """Linear interpolation between any two colours"""
    c0 = to_rgb(color0)
    c1 = to_rgb(color1)
    cn = ['red', 'green', 'blue']
    d = dict((n,[(0,s,s),(1,e,e)]) for (n,s,e) in zip(cn, c0, c1))
    return matplotlib.colors.LinearSegmentedColormap(name, d)

class ScalarColormapShading(object):
    """Returns color interpolated between two colors based on scalar value.
    
    Parameters:
    shade_param: parameter to look at for shading purposes
    min_val: minimum value of the parameter over the whole tree
    max_val: maximum value of the parameter over the whole tree
    color0: color to use when the parameter is at its minimum
    color1: color to use when the parameter is at its maximum
    
    Note: this is just a convenience wrapper for coloring. You can color
    the tree using any arbitrary function of the individual nodes by passing
    edge_color_callback to makeColorCallback (or to other objects
    that delegate to it).
    """
    def __init__(self, shade_param, min_val, max_val, cmap):
        """Returns a new callback for SimpleScalarShading: f(obj) -> color"""
        assert max_val, 'Need to know the maximum before shading can be done'
        self.min_val = min_val
        self.max_val = max_val
        self.shade_param = shade_param
        self.cmap = cmap
    
    def __call__(self, edge):
        """Returns color given node."""
        value = edge.params.get(self.shade_param, None)
        if value is None:
            return "grey"
        else:
            value_to_show = max(min(value, self.max_val), self.min_val)
            normed = (value_to_show-self.min_val)/(self.max_val-self.min_val)
            color = self.cmap(normed)
            return color    

def makeColorCallback(shade_param=None, min_value=0.0, max_value=None,
        edge_color="black", highlight_color="red", edge_color_callback=None,
        callback_returns_name=None, cmap=None):
    """Makes a callback f(node)->color using several strategies.
    
    The possibilities are:
    
    1.  You want all the nodes to have the same color. Pass in edge_color as
        something other than "black".
    
    2.  You want to shade the nodes from one color to another based on the value
        of a parameter. Pass the name of the parameter as a string to shade_param
        (e.g. the parameter might be "GC" for GC content). Pass in the max and
        min values (e.g. calculated from the range actually on the tree), and
        the colors for the "normal" (low) and "special" (high) values. The
        renderer will automatically calculate what the colors are.

    3.  You want some nodes to be one color, and other nodes to be highlighted
        in a different color. Pass in these colors as edge_color (e.g. "blue")
        and highlight_color (e.g. "green"). Set a parameter to 0 (for normal) or
        1 (for highlight), and pass in the name of this parameter as
        shade_param (e.g. the parameter might be "is_mammal", which you would
        set to 0 (False) or 1 (True) to highlight the mammals.
    
    4.  You have f(node) -> color. Pass in f as edge_color_callback.
    
    Alternatively, set the Color attribute of the dendrogram edges.
    """
    if callback_returns_name is not None:
        pass # give deprecation warning?, no longer needed
    
    edge_color = to_rgb(edge_color)
    highlight_color = to_rgb(highlight_color)
    
    if edge_color_callback is not None:
        return lambda edge:edge_color_callback(edge)
    elif shade_param:
        if cmap is None:
            cmap = SimpleColormap(edge_color, highlight_color)
        return ScalarColormapShading(
                shade_param, min_value, max_value, cmap)
    else:
        return lambda edge:edge_color
        
    
class MatplotlibRenderer(object):
    """Returns a matplitlib render including font size, stroke width, etc.
    
    Note: see documentation for makeColorCallback above to figure
    out how to make it color things the way you want. Dynamically varying the
    stroke width is not yet implemented by should be.
    """
    def __init__(self, font_size=None, stroke_width=3, label_pad=None, **kw):
        self.calculated_edge_color = makeColorCallback(**kw)
        self.text_opts = {}
        if font_size is not None:
            self.text_opts['fontsize'] = font_size
        self.line_opts = {}
        if stroke_width is not None:
            self.line_opts['linewidth'] = stroke_width
        if label_pad is None:
            label_pad = 8
        self.labelPadDistance = label_pad
    
    def edge_color(self, edge):
        if edge.Color is None:
            return self.calculated_edge_color(edge.original)
        else:
            return edge.Color
        
    def line(self, x1, y1, x2, y2, edge=None):
        opts = self.line_opts.copy()
        if edge is not None:
            opts['edgecolor'] = self.edge_color(edge)
        path = Path([(x1, y1), (x2, y2)], [Path.MOVETO, Path.LINETO])
        return PathPatch(path,  **opts)
    
    def polygon(self, vertices, color):
        opts = self.line_opts.copy()
        opts['color'] = color
        return Polygon(vertices, **opts)
        
    def string(self, x, y, string, ha=None, va=None, rotation=None, color=None):
        opts = self.text_opts.copy()
        if ha is not None:
            opts['ha'] = ha
        if va is not None:
            opts['va'] = va
        if rotation is not None:
            opts['rotation'] = rotation
        if color is not None:
            opts['color'] = color
        return Text(x, y, string, **opts)
        
class DendrogramLabelStyle(object):
    """Label options"""
    
    def __init__(self, show_params=None, show_internal_labels=False,
            label_template=None, edge_label_callback=None):
        if edge_label_callback is None:
            if label_template is None:
                if hasattr(show_params, "__contains__"):
                    if len(show_params) == 1:
                        label_template = "%%(%s)s" % show_params[0]
                    else:
                        label_template = "\n".join(
                            ["%s: %%(%s)s" % (p,p) for p in show_params])
                elif show_params:
                    label_template = "%s"
                else:
                    label_template = ""
            
            def edge_label_callback(edge):
                try:
                    if hasattr(label_template, 'substitute'):
                        # A new style (as of Py 2.4?) string template
                        return label_template.substitute(edge.params)
                    else:
                        return label_template % edge.params
                except KeyError:
                    return ""  # param missing - probably just the root edge
                        
        self.edgeLabelCallback = edge_label_callback
        self.showInternalLabels = show_internal_labels
    
    def getEdgeLabel(self, edge):
        return self.edgeLabelCallback(edge)
    
    def getNodeLabel(self, edge):
        if edge.Name is not None:
            return edge.Name
        elif self.showInternalLabels or not edge.Children:
            return edge.original.Name
        else:
            return ""

def ValidColorProperty(real_name, doc='A color name or other spec'):
    """Can only be set to Null or a valid color"""
    def getter(obj):
        return getattr(obj, real_name, None)
    def setter(obj, value):
        if value is not None: to_rgb(value)
        setattr(obj, real_name, value)
    def deleter(obj):
        setattr(obj, real_name, None)
    return property(getter, setter, deleter, doc)

class _Dendrogram(rlg2mpl.Drawable, TreeNode):
    # One of these for each tree edge.  Extra attributes:
    #    depth - distance from root to bottom of edge
    #    height - max distance from a decendant leaf to top of edge
    #    width - number of decendant leaves
    # note these are named tree-wise, not geometricaly, so think
    # of a vertical tree (for this part anyway)
    #
    #   x1, y1, x2, y2 - coordinates
    # these are horizontal / vertical as you would expect
    #
    # The algorithm is split into 4 passes over the tree for easier
    # code reuse - vertical drawing, new tree styles, new graphics
    # libraries etc.
    
    aspect_distorts_lengths = True

    def __init__(self, edge, use_lengths=True):
        children = [type(self)(child) for child in edge.Children]
        TreeNode.__init__(self, Params=edge.params.copy(), Children=children, 
            Name=("" if children else edge.Name))
        self.Length = edge.Length
        self.original = edge  # for edge_color_callback
        self.Collapsed = False
        self.use_lengths_default = use_lengths
    
    # Colors are properties so that invalid color names are caught immediately
    Color = ValidColorProperty('_Color', 'Color of line segment')
    NameColor = ValidColorProperty('_NameColor', 'Color of node name')
    CladeColor = ValidColorProperty('_CladeColor', 'Color of collapsed descendants')
    
    def __repr__(self):
        return '%s %s %s %s' % (
                self.depth, self.length, self.height, self.Children)
    
    def updateGeometry(self, use_lengths, depth=None, track_coordinates=None):
        """Calculate tree node attributes such as height and depth.
        Despite the name this first pass is ignorant of issues like
        scale and orientation"""
        
        if self.Length is None or not use_lengths:
            if depth is None:
                self.length = 0
            else:
                self.length = 1
        else:
            self.length = self.Length
        
        self.depth = (depth or 0) + self.length
        
        children = self.Children
        if children:
            for c in children:
                c.updateGeometry(use_lengths, self.depth, track_coordinates)
            self.height = max([c.height for c in children]) + self.length
            self.leafcount  = sum([c.leafcount for c in children])
            self.edgecount  = sum([c.edgecount for c in children]) + 1
            self.longest_label = max([c.longest_label for c in children],
                    key=len)
        else:
            self.height = self.length
            self.leafcount = self.edgecount = 1
            self.longest_label = self.Name or ''
        
        if track_coordinates is not None and self.Name != "root":
            self.track_y = track_coordinates[self.Name]
        else:
            self.track_y = 0

    def coords(self, height, width):
        """Return list of [node_name, node_id, x, y, child_ids]"""
        self.asArtist(height, width)
        result = []
        for node in self.postorder(include_self=True):
            result.append([node.Name, id(node), node.x2, node.y2] + [map(id, node.Children)])
        return result
    
    def makeFigure(self, width=None, height=None, margin=.25, use_lengths=None, **kw):
        (width, height),posn,kw = rlg2mpl.figureLayout(width, height, margin=0,
            default_aspect=0.5, leftovers=True, **kw)
        fig = self._makeFigure(width, height)
        ax = fig.add_axes(posn, frameon=False)
        width = 72 * posn[2] * fig.get_figwidth()
        height = 72 * posn[3] * fig.get_figheight()
        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_xticks([])
        ax.set_yticks([])
        if use_lengths is None:
            use_lengths = self.use_lengths_default
        else:
            pass # deprecate setting use_lengths here?
        if use_lengths and self.aspect_distorts_lengths:
            ax.set_aspect('equal')
        g = self.asArtist(width, height, use_lengths=use_lengths, 
            margin=margin*72, **kw)
        ax.add_artist(g)
        return fig
    
    def asArtist(self, width, height, margin=20, use_lengths=None,
            scale_bar="left", show_params=None, show_internal_labels=False,
            label_template=None, edge_label_callback=None, shade_param=None, 
            max_value=None, font_size=None, **kw):
        
        if use_lengths is None:
            use_lengths = self.use_lengths_default
        self.updateGeometry(use_lengths=use_lengths)
        
        if width <= 2 * margin:
            raise ValueError('%spt not wide enough for %spt margins' %
                    (width, margin))
        if height <= 2 * margin:
            raise ValueError('%spt not high enough for %spt margins' %
                    (height, margin))
        width -= 2 * margin
        height -= 2 * margin

        label_length = len(self.longest_label)
        label_width = label_length * 0.8 * (font_size or 10) # not very accurate
        (left_labels, right_labels) = self.labelMargins(label_width)
        total_label_width = left_labels + right_labels
        if width < total_label_width:
            raise ValueError('%spt not wide enough for ""%s"' %
                    (width, self.longest_label))
    
        scale = self.updateCoordinates(width-total_label_width, height)

        if shade_param is not None and max_value is None:
            for edge in self.postorder(include_self=True):
                sp = edge.params.get(shade_param, None)
                if max_value is None or sp > max_value:
                    max_value = sp
        renderer = MatplotlibRenderer(shade_param=shade_param, 
                max_value=max_value, font_size=font_size, **kw)

        labelopts = {}
        for labelopt in ['show_params', 'show_internal_labels', 
                'label_template', 'edge_label_callback']:
            labelopts[labelopt] = locals()[labelopt]
        label_style = DendrogramLabelStyle(**labelopts)

        ss = self._draw(renderer, label_style)
        if use_lengths:
            # Placing the scale properly might take some work,
            # for now just always put it in a bottom corner.
            unit = 10**min(0.0, numpy.floor(numpy.log10(width/scale/2.0)))
            if scale_bar == "right":
                x1, x2 = (width-scale*unit, width)
            elif scale_bar == "left":
                x1, x2 = (-left_labels, scale*unit-left_labels)
            else:
                assert not scale_bar, scale_bar
            if scale_bar:
                ss.append(renderer.line(x1, 0.0, x2, 0.0))
                ss.append(renderer.string((x1+x2)/2, 5, str(unit), va='bottom', ha='center'))
        
        g = rlg2mpl.Group(*ss)
        g.translate(margin+left_labels, margin)
        return g
    
    def _draw(self, renderer, label_style):
        g = []
        g += self._draw_edge(renderer, label_style)
        if self.Collapsed:
            g += self._draw_collapsed_clade(renderer, label_style)
        else:
            g += self._draw_node(renderer, label_style)
            for child in self.Children:
                g += child._draw(renderer, label_style)
            g += self._draw_node_label(renderer, label_style)
        return g
    
    def _draw_node(self, renderer, label_style):
        g = []
        # Joining line for square form
        if self.Children:
            cys = [c.y1 for c in self.Children] + [self.y2]
            if max(cys) > min(cys):
                g.append(renderer.line(self.x2, min(cys), self.x2, max(cys), self))
        return g
    
    def _draw_edge(self, renderer, label_style):
        g = []
        if ((self.x1, self.y1) == (self.x2, self.y2)):
            # avoid labeling zero length line, eg: root
            return g

        # Main line
        g.append(renderer.line(self.x1, self.y1, self.x2, self.y2, self))
        
        # Edge Label
        text = label_style.getEdgeLabel(self)
        if text:
            midx, midy = (self.x1+self.x2)/2, (self.y1+self.y2)/2
            if self.x1 == self.x2:
                rot = 0
            else:
                rot = numpy.arctan((self.y2-self.y1)/(self.x2-self.x1))
            midx += numpy.cos(rot+numpy.pi/2)*3
            midy += numpy.sin(rot+numpy.pi/2)*3
            g.append(renderer.string(midx, midy, text, ha='center', va='bottom', 
                    rotation=180/numpy.pi*rot))
        return g
    
    def _draw_node_label(self, renderer, label_style):
        text = label_style.getNodeLabel(self)
        color = self.NameColor
        (x, ha, y, va) = self.getLabelCoordinates(text, renderer)
        return [renderer.string(x, y, text, ha=ha, va=va, color=color)]
        
    def _draw_collapsed_clade(self, renderer, label_style):
        text = label_style.getNodeLabel(self)
        color = _first_non_none([self.CladeColor, self.Color, 'black'])
        icolor = 'white' if sum(to_rgb(color))/3 < 0.5 else 'black'
        g = []
        if not self.Children:
            return g
        (l,r,t,b), vertices = self.wedgeVertices()
        g.append(renderer.polygon(vertices, color))
        if not b <= self.y2 <= t:
            # ShelvedDendrogram needs this extra line segment
            g.append(renderer.line(self.x2, self.y2, self.x2, b, self))
        (x, ha, y, va) = self.getLabelCoordinates(text, renderer)
        g.append(renderer.string(
                (self.x2+r)/2, (t+b)/2, str(self.leafcount), ha=ha, va=va,
                color=icolor))
        g.append(renderer.string(
                x-self.x2+r, y, text, ha=ha, va=va, color=self.NameColor))
        return g
    
    def setCollapsed(self, collapsed=True, label=None, color=None):
        if color is not None:
            self.CladeColor = color
        if label is not None:
            self.Name = label
        self.Collapsed = collapsed


class Dimensions(object):
    def __init__(self, xscale, yscale, total_tree_height):
        self.x = xscale
        self.y = yscale
        self.height = total_tree_height
    

class _RootedDendrogram(_Dendrogram):
    """_RootedDendrogram subclasses provide yCoords and xCoords, which examine
    attributes of a node (its length, coodinates of its children) and return
    a tuple for start/end of the line representing the edge."""
    def labelMargins(self, label_width):
        return (0, label_width)
            
    def widthRequired(self):
        return self.leafcount
    
    def xCoords(self, scale, x1):
        raise NotImplementedError
    
    def yCoords(self, scale, y1):
        raise NotImplementedError
        
    def updateCoordinates(self, width, height):
        xscale = width / self.height
        yscale = height / self.widthRequired()
        scale = Dimensions(xscale, yscale, self.height)
        
        # y coords done postorder, x preorder, y first.
        # so it has to be done in 2 passes.
        self.update_y_coordinates(scale)
        self.update_x_coordinates(scale)
        return xscale
    
    def update_y_coordinates(self, scale, y1=None):
        """The second pass through the tree.  Y coordinates only
        depend on the shape of the tree and yscale"""
        if y1 is None:
            y1 = self.widthRequired() * scale.y
        child_y = y1
        for child in self.Children:
            child.update_y_coordinates(scale, child_y)
            child_y -= child.widthRequired() * scale.y
        (self.y1, self.y2) = self.yCoords(scale, y1)
    
    def update_x_coordinates(self, scale, x1=0):
        """For non 'square' styles the x coordinates will depend
        (a bit) on the y coodinates, so they should be done first"""
        (self.x1, self.x2) = self.xCoords(scale, x1)
        for child in self.Children:
            child.update_x_coordinates(scale, self.x2)
    
    def getLabelCoordinates(self, text, renderer):
        return (self.x2+renderer.labelPadDistance, 'left', self.y2, 'center')

class SquareDendrogram(_RootedDendrogram):
    aspect_distorts_lengths = False
    
    def yCoords(self, scale, y1):
        cys = [c.y1 for c in self.Children]
        if cys:
            y2 = (cys[0]+cys[-1]) / 2.0
        else:
            y2 = y1 - 0.5 * scale.y
        return (y2, y2)
    
    def xCoords(self, scale, x1):
        dx = scale.x * self.length
        x2 = x1 + dx
        return (x1, x2)

    def wedgeVertices(self):
        tip_ys = [(c.y2 + self.y2)/2 for c in self.iterTips()]
        t,b = max(tip_ys), min(tip_ys)
        cxs = [c.x2 for c in self.iterTips()]
        l,r = min(cxs), max(cxs)
        return (l,r,t,b), [(self.x2, b), (self.x2, t), (l, t), (r, b)]


class StraightDendrogram(_RootedDendrogram):
    def yCoords(self, scale, y1):
        # has a side effect of adjusting the child y1's to meet nodes' y2's
        cys = [c.y1 for c in self.Children]
        if cys:
            y2 = (cys[0]+cys[-1]) / 2.0
            distances = [child.length for child in self.Children]
            closest_child = self.Children[distances.index(min(distances))]
            dy = closest_child.y1 - y2
            max_dy = 0.8*max(5, closest_child.length*scale.x)
            if abs(dy) > max_dy:
                # 'moved', node.Name, y2, 'to within', max_dy,
                # 'of', closest_child.Name, closest_child.y1
                y2 = closest_child.y1 - _sign(dy) * max_dy
        else:
            y2 = y1 - scale.y / 2.0
        y1 = y2
        for child in self.Children:
            child.y1 = y2
        return (y1, y2)
    
    def xCoords(self, scale, x1):
        dx = self.length * scale.x
        dy = self.y2 - self.y1
        dx = numpy.sqrt(max(dx**2 - dy**2, 1))
        return (x1, x1 + dx)

    def wedgeVertices(self):
        tip_ys = [(c.y2 + self.y2)/2 for c in self.iterTips()]
        t,b = max(tip_ys), min(tip_ys)
        cxs = [c.x2 for c in self.iterTips()]
        l,r = min(cxs), max(cxs)
        vertices = [(self.x2, self.y2), (l, t), (r, b)]
        return (l,r,t,b), vertices

class _ContemporaneousMixin(object):
    """A dendrogram with all of the tips lined up.  
    Tidy but not suitable for displaying evolutionary distances accurately"""

    # Overrides init to change default for use_lengths
    def __init__(self, edge, use_lengths=False):
        super(_ContemporaneousMixin, self).__init__(edge, use_lengths)
        
    def xCoords(self, scale, x1):
        return (x1, (scale.height-(self.height-self.length))*scale.x)

class ContemporaneousDendrogram(_ContemporaneousMixin, SquareDendrogram):
    pass
    
class ContemporaneousStraightDendrogram(_ContemporaneousMixin, StraightDendrogram):
    pass


class ShelvedDendrogram(ContemporaneousDendrogram):
    """A dendrogram in which internal nodes also get a row to themselves"""
    def widthRequired(self):
        return self.edgecount  # as opposed to tipcount
    
    def yCoords(self, scale, y1):
        cys = [c.y1 for c in self.Children]
        if cys:
            y2 = cys[-1] - 1.0 * scale.y
        else:
            y2 = y1 - 0.5 * scale.y
        return (y2, y2)

class AlignedShelvedDendrogram(ShelvedDendrogram):
    
    def update_y_coordinates(self, scale, y1=None):
        """The second pass through the tree.  Y coordinates only
        depend on the shape of the tree and yscale"""
        for child in self.Children:
            child.update_y_coordinates(scale, None)
        (self.y1, self.y2) = self.yCoords(scale, None)
    
    def yCoords(self, scale, y1):
        if hasattr(self, 'track_y'):
            return (self.track_y, self.track_y)
        else:
            raise RuntimeError, self.Name
    

class UnrootedDendrogram(_Dendrogram):
    aspect_distorts_lengths = True

    def labelMargins(self, label_width):
        return (label_width, label_width)
    
    def wedgeVertices(self):
        tip_dists = [(c.depth-self.depth)*self.scale for c in self.iterTips()]
        (near, far) = (min(tip_dists), max(tip_dists))
        a = self.angle - 0.25 * self.wedge
        (x1, y1) = (self.x2+near*numpy.sin(a), self.y2+near*numpy.cos(a))
        a = self.angle + 0.25 * self.wedge
        (x2, y2) = (self.x2+far*numpy.sin(a), self.y2+far*numpy.cos(a))
        vertices = [(self.x2, self.y2), (x1, y1), (x2, y2)]
        return (self.x2, (x1+x2)/2, self.y2, (y1+y2)/2), vertices

    def updateCoordinates(self, width, height):
        angle = 2*numpy.pi / self.leafcount
        # this loop is a horrible brute force hack
        # there are better (but complex) ways to find
        # the best rotation of the tree to fit the display.
        best_scale = 0
        for i in range(60):
            direction = i/60.0*numpy.pi
            points = self._update_coordinates(1.0, 0, 0, direction, angle)
            xs = [x for (x,y) in points]
            ys = [y for (x,y) in points]
            scale = min(float(width)/(max(xs)-min(xs)), float(height)/(max(ys)-min(ys)))
            scale *= 0.95 # extra margin for labels
            if scale > best_scale:
                best_scale = scale
                mid_x = width/2-((max(xs)+min(xs))/2)*scale
                mid_y = height/2-((max(ys)+min(ys))/2)*scale
                best_args = (scale, mid_x, mid_y, direction, angle)
        self._update_coordinates(*best_args)
        return best_scale
    
    def _update_coordinates(self, s, x1, y1, a, da):
        # Constant angle algorithm.  Should add maximim daylight step.
        (x2, y2) = (x1+self.length*s*numpy.sin(a), y1+self.length*s*numpy.cos(a))
        (self.x1, self.y1, self.x2, self.y2, self.angle) = (x1, y1, x2, y2, a)
        if self.Collapsed:
            self.wedge = self.leafcount * da
            self.scale = s
            (l,r,t,b), vertices = self.wedgeVertices()
            return vertices
            
        a -= self.leafcount * da / 2
        if not self.Children:
            points = [(x2, y2)]
        else:
            points = []
            for (i,child) in enumerate(self.Children):
                ca = child.leafcount * da
                points += child._update_coordinates(s, x2, y2, a+ca/2, da)
                a += ca
        return points
    
    def getLabelCoordinates(self, text, renderer):
        (dx, dy) = (numpy.sin(self.angle), numpy.cos(self.angle))
        pad = renderer.labelPadDistance
        (x, y) = (self.x2+pad*dx, self.y2+pad*dy)
        if dx > abs(dy):
            return (x, 'left', y, 'center')
        elif -dx > abs(dy):
            return (x, 'right', y, 'center')
        elif dy > 0:
            return (x, 'center', y, 'bottom')
        else:
            return (x, 'center', y, 'top')
    
