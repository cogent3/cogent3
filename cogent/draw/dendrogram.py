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
from reportlab.graphics import shapes, renderPDF
from reportlab.lib import colors
import numpy

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Zongzhi Liu", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.1"
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
        else:
            return max(child_results)
    return _max


class SimpleScalarShading(object):
    """Returns color interpolated between two colors based on scalar value.
    
    Parameters:
    shade_param: parameter to look at for shading purposes
    min_val: minimum value of the parameter over the whole tree
    max_val: maximum value of the parameter over the whole tree
    color0: color to use when the parameter is at its minimum
    color1: color to use when the parameter is at its maximum
    
    Note: this is just a convenience wrapper for coloring. You can color
    the tree using any arbitrary function of the individual nodes by passing
    edge_color_callback to makeReportlabColorCallback (or to other objects
    that delegate to it).
    """
    def __init__(self, shade_param, min_val, max_val, color0, color1):
        """Returns a new callback for SimpleScalarShading: f(obj) -> color"""
        assert max_val, 'Need to know the maximum before shading can be done'
        self.color0 = color0
        self.color1 = color1
        self.min_val = min_val
        self.max_val = max_val
        self.shade_param = shade_param
    
    def __call__(self, edge):
        """Returns color given node."""
        value = edge.params.get(self.shade_param, None)
        if value is None:
            return colors.grey
        else:
            value_to_show = max(min(value, self.max_val), self.min_val)
            color = colors.linearlyInterpolatedColor(self.color0, self.color1,
                    self.min_val, self.max_val, value_to_show)
            return color
    

class ColorNameCallbackWrapper(object):
    """Returns simple wrapper for an arbitrary f(node) -> color_name.
    
    Takes on initialization an arbitrary function of a node that returns a
    reportlab color (e.g. 'red'). Returns the actual reportlab color
    corresponding to this color name.
    """
    
    def __init__(self, color_callback):
        self.color_callback = color_callback
    
    def __call__(self, edge):
        color_name = self.color_callback(edge)
        return getattr(colors, color_name)
    

def makeReportlabColorCallback(shade_param=None, min_value=0.0, max_value=None,
        edge_color="black", highlight_color="red", edge_color_callback=None,
        callback_returns_name=True):
    """Makes a callback f(node)->color using several strategies.
    
    The possibilities are:
    1.  You have f(node) -> color. Pass in f as edge_color_callback, and
        callback_returns_name as False.
    
    2.  You have f(node) -> color_name. Pass in f as edge_color_callback, and
        leave callback_returns_name as True (the default). The result will be
        a callabale object that returns the color given the node.
    
    3.  You want all the nodes to have the same color. Pass in edge_color as
        something other than "black".
    
    4.  You want some nodes to be one color, and other nodes to be highlighted
        in a different color. Pass in these colors as edge_color (e.g. "blue")
        and highlight_color (e.g. "green"). Set a parameter to 0 (for normal) or
        1 (for highlight), and pass in the name of this parameter as
        shade_param (e.g. the parameter might be "is_mammal", which you would
        set to 0 (False) or 1 (True) to highlight the mammals.
    
    5.  You want to shade the nodes from one color to another based on the value
        of a parameter. Pass the name of the parameter as a string to shade_param
        (e.g. the parameter might be "GC" for GC content). Pass in the max and
        min values (e.g. calculated from the range actually on the tree), and
        the colors for the "normal" (low) and "special" (high) values. The
        renderer will automatically calculate what the colors are.
    """
    
    if edge_color_callback is not None:
        if callback_returns_name:
            return ColorNameCallbackWrapper(edge_color_callback)
        else:
            return edge_color_callback
    if isinstance(edge_color, basestring):
        edge_color = getattr(colors, edge_color)
    if isinstance(highlight_color, basestring):
        highlight_color = getattr(colors, highlight_color)
    if shade_param:
        return SimpleScalarShading(shade_param,
                min_value, max_value, edge_color, highlight_color)
    else:
        return lambda edge:edge_color
    

class ReportlabRenderer(object):
    """Returns a reportlab render including font size, stroke width, etc.
    
    Note: see documentation for makeReportlabColorCallback above to figure
    out how to make it color things the way you want. Dynamically varying the
    stroke width is not yet implemented by should be.
    """
    def __init__(self, font_size=None, stroke_width=3, label_pad=None, **kw):
        self.edge_color = makeReportlabColorCallback(**kw)
        self.text_opts = {}
        if font_size is not None:
            self.text_opts['fontSize'] = font_size
        self.line_opts = {'strokeLineCap':1}
        if stroke_width is not None:
            self.line_opts['strokeWidth'] = stroke_width
        if label_pad is None:
            label_pad = self.stringWidth(' ')
        self.labelPadDistance = label_pad
    
    def line(self, x1, y1, x2, y2, edge=None):
        opts = self.line_opts.copy()
        if edge is not None:
            opts['strokeColor'] = self.edge_color(edge)
        return shapes.Line(x1, y1, x2, y2, **opts)
    
    def string(self, x, y, string):
        return shapes.String(x, y, string, **self.text_opts)
    
    def stringWidth(self, text):
        return shapes.String(0, 0, text, **self.text_opts).getEast()
    
    def stringHeight(self, text):
        return self.text_opts.get('fontSize', 10) * len(text.split('\n'))
    
    def makeDrawing(self, shape_list, width, height, margin):
        g = shapes.Group(*shape_list)
        g.translate(margin, margin)
        D = shapes.Drawing(width+2*margin, height+2*margin)
        D.add(g)
        return D
    

class DendrogramLabelStyle(object):
    """Label options"""
    
    def __init__(self, show_params=None, show_internal_labels=False,
            label_template=None):
        if label_template is None:
            if hasattr(show_params, "__contains__"):
                label_template = "\n".join(
                        ["%s: %%(%s)s" % (p,p) for p in show_params])
            elif show_params:
                label_template = "%s"
            else:
                label_template = ""
        
        self.labelTemplate = label_template
        self.showInternalLabels = show_internal_labels
    
    def getEdgeLabel(self, edge):
        try:
            text = self.labelTemplate % edge.params
        except KeyError:
            text = ''  # param missing - probably just the root edge
        return text
    
    def getNodeLabel(self, edge):
        if self.showInternalLabels or not edge.Children:
            return edge.Name or ""
        else:
            return ""
    

class _Dendrogram(object):
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
            self.leafcount    = sum([c.leafcount     for c in self.children])
            self.edgecount    = sum([c.edgecount     for c in self.children]) + 1
            self.labelwidth    = max([c.labelwidth for c in self.children])
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
        self.draw(ReportlabRenderer, height, width)
        result = []
        f = lambda node, ignored: result.append([node.edge.Name, id(node), node.x2, node.y2] + [map(id, node.children)])
        self.postorder(f)
        return result

    
    # ReportLab specific from here - maybe should split class
    
    def drawToPDF(self, filename, *args, **kw):
        D = self.asDrawing(*args, **kw)
        renderPDF.drawToFile(D, filename)
    
    def asDrawing(self, width, height, **kw):
        return self.draw(ReportlabRenderer, width, height, **kw)
    
    def draw(self, rendererClass, width, height, margin=10, use_lengths=True,
            scale_bar="left", show_params=None, show_internal_labels=False,
            label_template=None, **kw):
        """A reportlab drawing"""
        
        label_style = DendrogramLabelStyle(
                show_params=show_params,
                show_internal_labels=show_internal_labels,
                label_template=label_template)
        if kw.get('shade_param', None) is not None and \
                kw.get('max_value', None) is None:
            kw['max_value'] = self.postorder(_treemaxof(kw['shade_param']))
        renderer = rendererClass(**kw)
        if use_lengths is None:
            use_lengths = self.use_lengths_default
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
            unit = 10**numpy.floor(numpy.log10(width/scale/2.0))
            if scale_bar == "right":
                x1, x2 = (width-scale*unit, width)
            elif scale_bar == "left":
                x1, x2 = (0, scale*unit)
            else:
                assert not scale_bar, scale_bar
            if scale_bar:
                ss.append(renderer.line(x1, 0.0, x2, 0.0))
                ss.append(renderer.string((x1+x2)/2, 5, str(unit)))
        
        return renderer.makeDrawing(ss, width, height, margin)
    
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
        
        # Main line
        g.append(renderer.line(self.x1, self.y1, self.x2, self.y2, self.edge))
        
        # Edge Label
        text = label_style.getEdgeLabel(self.edge)
        if text:
            midx, midy = (self.x1+self.x2)/2, (self.y1+self.y2)/2
            label_width = renderer.stringWidth(text)
            midx -= label_width / 2
            midy -=  renderer.stringHeight(text) * 1.25
            g.append(renderer.string(midx, midy, text))
        return g
    
    def _draw_node_label(self, renderer, label_style):
        text = label_style.getNodeLabel(self.edge)
        (x, y) = self.getLabelCoordinates(text, renderer)
        label_height = renderer.stringHeight(text)
        return [renderer.string(x, y+label_height/3, text)]
    

class _RootedDendrogramForm(object):
    """A rooted dendrogram form defines how lengths get mapped to X and Y coodinates.
    
    _RootedDendrogramStyle subclasses provide yCoords and xCoords, which examine
    attributes of a node (its length, coodinates of its children) and return
    a tuple for start/end of the line representing the edge."""
    
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
        return (self.x2+renderer.labelPadDistance, self.y2-renderer.stringHeight(text)/2)

class SquareDendrogram(_RootedDendrogram):
    FormClass = SquareDendrogramForm
    use_lengths_default = True

class StraightDendrogram(_RootedDendrogram):
    FormClass = StraightDendrogramForm
    use_lengths_default = True

class ContemporaneousDendrogram(_RootedDendrogram):
    FormClass = ContemporaneousDendrogramForm
    use_lengths_default = False

class ShelvedDendrogram(_RootedDendrogram):
    FormClass = ShelvedDendrogramForm
    use_lengths_default = False

class AlignedShelvedDendrogram(_RootedDendrogram):
    FormClass = AlignedShelvedDendrogramForm
    use_lengths_default = False
    
    def update_y_coordinates(self, style, y1=None):
        """The second pass through the tree.  Y coordinates only
        depend on the shape of the tree and yscale"""
        for child in self.children:
            child.update_y_coordinates(style, None)
        (self.y1, self.y2) = style.yCoords(self, None)
    

class ContemporaneousStraightDendrogram(_RootedDendrogram):
    FormClass = ContemporaneousStraightDendrogramForm
    use_lengths_default = False

class UnrootedDendrogram(_Dendrogram):
    use_lengths_default = True
    
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
        ty = renderer.stringHeight(text) / 2
        tx = renderer.stringWidth(text) / 2
        (x, y) = (self.x2+(tx+pad)*dx, self.y2+(ty+pad)*dy)
        return (x-tx, y-ty)
    
