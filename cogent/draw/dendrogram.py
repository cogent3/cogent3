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
import rlg2mpl
import matplotlib.colors
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.text import Text
import numpy

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Zongzhi Liu", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def _sign(x):
    """Returns True if x is positive, False otherwise."""
    return x and x/abs(x)

def _treemaxof(param):
    """Returns maximum value of a parameter in a tree."""
    def _max(node, child_results):
        params = node.edge.params
        if param in params:
            return max([params[param]] + child_results)
        elif child_results:
            return max(child_results)
        else:
            return None
    return _max

def SimpleColormap(color0, color1, name=None):
    """Linear interpolation between any two colours"""
    to_rgb = matplotlib.colors.colorConverter.to_rgb
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
    """
    if callback_returns_name is not None:
        pass # give deprecation warning?, no longer needed
        
    if edge_color_callback is not None:
        return edge_color_callback
    elif shade_param:
        if cmap is None:
            cmap = SimpleColormap(edge_color, highlight_color)
        return ScalarColormapShading(shade_param, min_value, max_value, cmap)
    else:
        return lambda edge:edge_color
    
    
class MatplotlibRenderer(object):
    """Returns a matplitlib render including font size, stroke width, etc.
    
    Note: see documentation for makeColorCallback above to figure
    out how to make it color things the way you want. Dynamically varying the
    stroke width is not yet implemented by should be.
    """
    def __init__(self, font_size=None, stroke_width=3, label_pad=None, **kw):
        self.edge_color = makeColorCallback(**kw)
        self.text_opts = {}
        if font_size is not None:
            self.text_opts['fontsize'] = font_size
        self.line_opts = {}
        if stroke_width is not None:
            self.line_opts['linewidth'] = stroke_width
        if label_pad is None:
            label_pad = self.stringWidth(' ')
        self.labelPadDistance = label_pad
    
    def line(self, x1, y1, x2, y2, edge=None):
        path = Path([(x1, y1), (x2, y2)], [Path.MOVETO, Path.LINETO])
        opts = self.line_opts.copy()
        if edge is not None:
            opts['edgecolor'] = self.edge_color(edge)
        return PathPatch(path,  **opts)
    
    def string(self, x, y, string, ha=None, va=None, rotation=None):
        opts = self.text_opts.copy()
        if ha is not None:
            opts['ha'] = ha
        if va is not None:
            opts['va'] = va
        if rotation is not None:
            opts['rotation'] = rotation
        return Text(x, y, string, **opts)
    
    def stringWidth(self, text):
        #rlg2mpl.String(0, 0, text, **self.text_opts).getEast()
        return self.text_opts.get('fontSize',10) * .8 * len(text) # not very accurate!
    
    def stringHeight(self, text):
        return self.text_opts.get('fontSize', 10) * len(text.split('\n'))
    
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
        if self.showInternalLabels or not edge.Children:
            return edge.Name
        else:
            return ""
        
class _Dendrogram(rlg2mpl.Drawable):
    # One of these for each tree edge.  Attributes:
    #    edge - the real tree node
    #    children - Dendrograms of self.edge.Children
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
    
    def __init__(self, edge):
        if hasattr(edge, 'TreeRoot'):
            edge = edge.TreeRoot
        self.edge = edge
        self.children = [self.__class__(c) for c in edge.Children]
        self._up_to_date = False
    
    def __repr__(self):
        return '%s %s %s %s' % (
                self.depth, self.length, self.height, self.children)
    
    def postorder(self, f):
        rs = [child.postorder(f) for child in self.children]
        return f(self, rs)
    
    def updateGeometry(self, renderer, use_lengths=True,
            depth=None, track_coordinates=None):
        """Calculate tree node attributes such as height and depth.
        Despite the name this first pass is ignorant of issues like
        scale and orientation"""
        
        if self.edge.Length is None or not use_lengths:
            if depth is None:
                self.length = 0
            else:
                self.length = 1
        else:
            self.length = self.edge.Length
        
        self.depth = (depth or 0) + self.length
        
        if self.children:
            for c in self.children:
                c.updateGeometry(renderer, use_lengths, self.depth, track_coordinates)
            self.height = max([c.height for c in self.children]) + self.length
            self.leafcount  = sum([c.leafcount for c in self.children])
            self.edgecount  = sum([c.edgecount for c in self.children]) + 1
            self.labelwidth = max([c.labelwidth for c in self.children])
        else:
            self.height = self.length
            self.leafcount = self.edgecount = 1
            node_label = self.edge.Name or ''
            self.labelwidth = renderer.stringWidth(node_label)
        
        if track_coordinates is not None and self.edge.Name != "root":
            self.track_y = track_coordinates[self.edge.Name]
        else:
            self.track_y = 0

    def coords(self, height, width):
        """Return list of [node_name, node_id, x, y, child_ids]"""
        self.asArtist(height, width)
        result = []
        def _f(node, child_results):
            result.append([node.edge.Name, id(node), node.x2, node.y2] + [map(id, node.children)])
        self.postorder(_f)
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
        if use_lengths and self.aspect_distorts_lengths:
            ax.set_aspect('equal')
        g = self.asArtist(width, height, use_lengths=use_lengths, 
            margin=margin*72, **kw)
        ax.add_artist(g)
        return fig
    
    def asArtist(self, width, height, margin=20, use_lengths=None,
            scale_bar="left", show_params=None, show_internal_labels=False,
            label_template=None, edge_label_callback=None, **kw):
        """A reportlab drawing"""
        
        label_style = DendrogramLabelStyle(
                show_params = show_params,
                show_internal_labels = show_internal_labels,
                label_template = label_template,
                edge_label_callback = edge_label_callback,
                )
        if kw.get('shade_param', None) is not None and \
                kw.get('max_value', None) is None:
            kw['max_value'] = self.postorder(_treemaxof(kw['shade_param']))
        renderer = MatplotlibRenderer(**kw)
        self.updateGeometry(renderer, use_lengths=use_lengths)
        if width <= 2 * margin:
            raise ValueError('%spt not wide enough for %spt margins' %
                    (width, margin))
        if height <= 2 * margin:
            raise ValueError('%spt not high enough for %spt margins' %
                    (height, margin))
        width -= 2 * margin
        height -= 2 * margin
        scale = self.updateCoordinates(width, height)
        ss = self._draw(renderer, label_style)
        if use_lengths:
            # Placing the scale properly might take some work,
            # for now just always put it in a bottom corner.
            unit = 10**min(0.0, numpy.floor(numpy.log10(width/scale/2.0)))
            if scale_bar == "right":
                x1, x2 = (width-scale*unit, width)
            elif scale_bar == "left":
                x1, x2 = (0, scale*unit)
            else:
                assert not scale_bar, scale_bar
            if scale_bar:
                ss.append(renderer.line(x1, 0.0, x2, 0.0))
                ss.append(renderer.string((x1+x2)/2, 5, str(unit), va='bottom', ha='center'))
        
        g = rlg2mpl.Group(*ss)
        g.translate(margin, margin)
        return g
    
    def _draw(self, renderer, label_style):
        g = []
        g += self._draw_edge(renderer, label_style)
        for child in self.children:
            g += child._draw(renderer, label_style)
        g += self._draw_node_label(renderer, label_style)
        return g
    
    def _draw_edge(self, renderer, label_style):
        g = []
        # Joining line for square form
        if self.children:
            cys = [c.y1 for c in self.children] + [self.y2]
            if max(cys) > min(cys):
                g.append(renderer.line(self.x2, min(cys), self.x2, max(cys), self.edge))
        
        if ((self.x1, self.y1) == (self.x2, self.y2)):
            # avoid labeling zero length line, eg: root
            return g
            
        # Main line
        g.append(renderer.line(self.x1, self.y1, self.x2, self.y2, self.edge))
        
        # Edge Label
        text = label_style.getEdgeLabel(self.edge)
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
        text = label_style.getNodeLabel(self.edge)
        (x, ha, y, va) = self.getLabelCoordinates(text, renderer)
        return [renderer.string(x, y, text, ha=ha, va=va)]
    

class _RootedDendrogramForm(object):
    """A rooted dendrogram form defines how lengths get mapped to X and Y coodinates.
    
    _RootedDendrogramStyle subclasses provide yCoords and xCoords, which examine
    attributes of a node (its length, coodinates of its children) and return
    a tuple for start/end of the line representing the edge."""
    aspect_distorts_lengths = True
    
    def __init__(self, tree, width, height):
        self.yscale = 1.0
        self.yscale = height / self.widthRequiredFor(tree)
        self.xscale = width / tree.height
        self.total_tree_height = tree.height
    
    def widthRequiredFor(self, node):
        return node.leafcount * self.yscale
    
    def xCoords(self, node, x1):
        raise NotImplementedError
    
    def yCoords(self, node, x1):
        raise NotImplementedError
    

class SquareDendrogramForm(_RootedDendrogramForm):
    aspect_distorts_lengths = False
    
    def yCoords(self, node, y1):
        cys = [c.y1 for c in node.children]
        if cys:
            y2 = (cys[0]+cys[-1]) / 2.0
        else:
            y2 = y1 - self.yscale / 2.0
        return (y2, y2)
    
    def xCoords(self, node, x1):
        dx = self.xscale * node.length
        x2 = x1 + dx
        return (x1, x2)
    

class ContemporaneousDendrogramForm(SquareDendrogramForm):
    def xCoords(self, node, x1):
        return (x1, (self.total_tree_height-(node.height-node.length))*self.xscale)
    

class ShelvedDendrogramForm(ContemporaneousDendrogramForm):
    def widthRequiredFor(self, node):
        return node.edgecount * self.yscale
    
    def yCoords(self, node, y1):
        cys = [c.y1 for c in node.children]
        if cys:
            y2 = cys[-1] - 1.0 * self.yscale
        else:
            y2 = y1    - 0.5 * self.yscale
        return (y2, y2)
    

class AlignedShelvedDendrogramForm(ShelvedDendrogramForm):
    def __init__(self, tree, width, height):
        self.yscale = 1.0
        self.xscale = width / tree.height
        self.total_tree_height = tree.height
    
    def yCoords(self, node, y1):
        if hasattr(node, 'track_y'):
            return (node.track_y, node.track_y)
        else:
            raise RuntimeError, node.edge.Name
            return ShelvedDendrogramForm.yCoords(self, node, y1)
    
    def widthRequiredFor(self, node):
        raise RuntimeError, node.edge.Name
    

class StraightDendrogramForm(_RootedDendrogramForm):
    def yCoords(self, node, y1):
        # has a side effect of adjusting the child y1's to meet nodes' y2's
        cys = [c.y1 for c in node.children]
        if cys:
            y2 = (cys[0]+cys[-1]) / 2.0
            distances = [child.length for child in node.children]
            closest_child = node.children[distances.index(min(distances))]
            dy = closest_child.y1 - y2
            max_dy = 0.8*max(5, closest_child.length*self.xscale)
            if abs(dy) > max_dy:
                # 'moved', node.edge.Name, y2, 'to within', max_dy,
                # 'of', closest_child.edge.Name, closest_child.y1
                y2 = closest_child.y1 - _sign(dy) * max_dy
        else:
            y2 = y1 - self.yscale / 2.0
        y1 = y2
        for child in node.children:
            child.y1 = y2
        return (y1, y2)
    
    def xCoords(self, node, x1):
        dx = node.length * self.xscale
        dy = node.y2 - node.y1
        dx = numpy.sqrt(max(dx**2 - dy**2, 1))
        return (x1, x1 + dx)
    

class ContemporaneousStraightDendrogramForm(StraightDendrogramForm):
    def xCoords(self, node, x1):
        return (x1, (self.total_tree_height-(node.height-node.length))*self.xscale)
    

class _RootedDendrogram(_Dendrogram):
    def updateCoordinates(self, width, height):
        if width < self.labelwidth:
            raise ValueError('%spt not wide enough for %spt wide labels' %
                    (width, self.labelwidth))
        width -= self.labelwidth
        
        form = self.FormClass(self, width, height)
        
        # y coords done postorder, x preorder, y first.
        # so it has to be done in 2 passes.
        self.update_y_coordinates(form)
        self.update_x_coordinates(form)
        return form.xscale
    
    def update_y_coordinates(self, style, y1=None):
        """The second pass through the tree.  Y coordinates only
        depend on the shape of the tree and yscale"""
        if y1 is None:
            y1 = style.widthRequiredFor(self)
        child_y = y1
        for child in self.children:
            child.update_y_coordinates(style, child_y)
            child_y -= style.widthRequiredFor(child)
        (self.y1, self.y2) = style.yCoords(self, y1)
    
    def update_x_coordinates(self, style, x1=0):
        """For non 'square' styles the x coordinates will depend
        (a bit) on the y coodinates, so they should be done first"""
        (self.x1, self.x2) = style.xCoords(self, x1)
        for child in self.children:
            child.update_x_coordinates(style, self.x2)
    
    def getLabelCoordinates(self, text, renderer):
        return (self.x2+renderer.labelPadDistance, 'left', self.y2, 'center')

class SquareDendrogram(_RootedDendrogram):
    FormClass = SquareDendrogramForm
    aspect_distorts_lengths = FormClass.aspect_distorts_lengths
    use_lengths_default = True

class StraightDendrogram(_RootedDendrogram):
    FormClass = StraightDendrogramForm
    aspect_distorts_lengths = FormClass.aspect_distorts_lengths
    use_lengths_default = True

class ContemporaneousDendrogram(_RootedDendrogram):
    FormClass = ContemporaneousDendrogramForm
    aspect_distorts_lengths = FormClass.aspect_distorts_lengths
    use_lengths_default = False

class ShelvedDendrogram(_RootedDendrogram):
    FormClass = ShelvedDendrogramForm
    aspect_distorts_lengths = FormClass.aspect_distorts_lengths
    use_lengths_default = False

class AlignedShelvedDendrogram(_RootedDendrogram):
    FormClass = AlignedShelvedDendrogramForm
    aspect_distorts_lengths = FormClass.aspect_distorts_lengths
    use_lengths_default = False
    
    def update_y_coordinates(self, style, y1=None):
        """The second pass through the tree.  Y coordinates only
        depend on the shape of the tree and yscale"""
        for child in self.children:
            child.update_y_coordinates(style, None)
        (self.y1, self.y2) = style.yCoords(self, None)

class ContemporaneousStraightDendrogram(_RootedDendrogram):
    FormClass = ContemporaneousStraightDendrogramForm
    aspect_distorts_lengths = FormClass.aspect_distorts_lengths
    use_lengths_default = False

class UnrootedDendrogram(_Dendrogram):
    use_lengths_default = True
    aspect_distorts_lengths = True
    
    def updateCoordinates(self, width, height):
        if width < 2*self.labelwidth:
            raise ValueError('%spt not wide enough for %spt wide labels' %
                    (width, self.labelwidth))
        width -= 2*self.labelwidth
        
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
                best_args = (scale, mid_x+self.labelwidth, mid_y, direction, angle)
        self._update_coordinates(*best_args)
        return best_scale
    
    def _update_coordinates(self, s, x1, y1, a, da):
        # Constant angle algorithm.  Should add maximim daylight step.
        (x2, y2) = (x1+self.length*s*numpy.sin(a), y1+self.length*s*numpy.cos(a))
        (self.x1, self.y1, self.x2, self.y2, self.angle) = (x1, y1, x2, y2, a)
        a -= self.leafcount * da / 2
        if not self.children:
            points = [(x2, y2)]
        else:
            points = []
            for (i,child) in enumerate(self.children):
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
    
