from cogent3.core.tree import TreeNode
import numpy as np
from cogent3.draw.drawable import Drawable

__author__ = "Rahul Ghangas, Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell","Gavin Huttley","Rahul Ghangas"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class _Dendrogram(TreeNode, Drawable):
    aspect_distorts_lengths = True

    def __init__(self, tree, use_lengths=True):
        children = [type(self)(child) for child in tree.children]
        TreeNode.__init__(self, params=tree.params.copy(), children=children,
                          name=("" if children else tree.name))
        self.length = tree.length
        self.original = tree  # for edge_color_callback
        self.collapsed = False
        self.use_lengths_default = use_lengths
        self._trace = None

    def __repr__(self):
        return f'{self.length} {self.height} {self.children}'

    # make semi-private, def _update_geometry()..
    # (method used internally by object, but not for users)
    def _update_geometry(self, use_lengths, depth=None, track_coordinates=None):
        """Calculate tree node attributes such as height and depth.
        Despite the name this first pass is ignorant of issues like
        scale and orientation"""

        if self.length is None or not use_lengths:
            if depth is None:
                self.length = 0
            else:
                self.length = 1
        else:
            self.length = self.length

        self.depth = (depth or 0) + self.length

        children = self.children
        if children:
            for c in children:
                c._update_geometry(use_lengths, self.depth, track_coordinates)
            self.height = max([c.height for c in children]) + self.length
            self.leafcount = sum([c.leafcount for c in children])
            self.edgecount = sum([c.edgecount for c in children]) + 1
            self.longest_label = max([c.longest_label for c in children],
                                     key=len)
        else:
            self.height = self.length
            self.leafcount = self.edgecount = 1
            self.longest_label = self.name or ''

        if track_coordinates is not None and self.name != "root":
            self.track_y = track_coordinates[self.name]
        else:
            self.track_y = 0

    # semi-private
    def depth_dict(self):
        self._update_geometry(use_lengths=self.use_lengths_default)
        depth_dict = dict()

        for node in self.get_edge_vector():
            depth_dict[node] = node.depth

        return depth_dict

    def get_trace(self):
        try:
            self._trace
        except NameError:
            print("No existing figure detected: Creating new")
            self._set_initial_layout()
        else:
            pass

        return self._trace['data']

    def get_layout(self):
        return

    def set_layout(self, layout_updates):
        self._trace['layout'] = layout_updates

    def update_layout(self, layout_updates):
        self._trace['layout'].update(layout_updates)

    def _set_initial_layout(self, width=None, height=None, title=None, use_lengths=None, **kw):
        super(_Dendrogram, self)._set_initial_layout(width, height)
        if use_lengths is None:
            use_lengths = self.use_lengths_default
        else:
            pass  # deprecate setting use_lengths here?
        if use_lengths and self.aspect_distorts_lengths:
            self._trace['layout'].update(scene=dict(aspectmode="data"))

        if title is not None:
            self._trace['layout']['title'] = title


class Dimensions(object):

    def __init__(self, xscale, yscale, total_tree_height):
        self.x = xscale
        self.y = yscale
        self.height = total_tree_height


class _RootedDendrogram(_Dendrogram):
    """_RootedDendrogram subclasses provide y_coords and x_coords, which examine
    attributes of a node (its length, coordinates of its children) and return
    a tuple for start/end of the line representing the edge."""

    # semi-private
    def width_required(self):
        return self.leafcount

    # semi-private
    def x_coords(self, scale, x1):
        raise NotImplementedError

    # semi-private
    def y_coords(self, scale, y1):
        raise NotImplementedError

    # semi-private
    def update_coordinates(self, width, height):
        xscale = width / self.height
        yscale = height / self.width_required()
        scale = Dimensions(xscale, yscale, self.height)

        # y coords done postorder, x preorder, y first.
        # so it has to be done in 2 passes.
        self.update_y_coordinates(scale)
        self.update_x_coordinates(scale)
        return xscale

    # semi-private
    def update_y_coordinates(self, scale, y1=None):
        """The second pass through the tree.  Y coordinates only
        depend on the shape of the tree and yscale"""
        if y1 is None:
            y1 = self.width_required() * scale.y
        child_y = y1
        for child in self.children:
            child.update_y_coordinates(scale, child_y)
            child_y -= child.width_required() * scale.y
        (self.y1, self.y2) = self.y_coords(scale, y1)

    # semi-private
    def update_x_coordinates(self, scale, x1=0):
        """For non 'square' styles the x coordinates will depend
        (a bit) on the y coodinates, so they should be done first"""
        (self.x1, self.x2) = self.x_coords(scale, x1)
        for child in self.children:
            child.update_x_coordinates(scale, self.x2)


def _calc_row(clade, ycoords):
    for subclade in clade.children:
        if subclade not in ycoords:
            ycoords = _calc_row(subclade, ycoords)
    ycoords[clade] = (ycoords[clade.children[0]] +
                      ycoords[clade.children[-1]]) / 2
    return ycoords


class SquareDendrogram(_RootedDendrogram):
    aspect_distorts_lengths = False

    # semi-private
    def y_coords(self, scale, y1):
        cys = [c.y1 for c in self.children]
        if cys:
            y2 = (cys[0] + cys[-1]) / 2.0
        else:
            y2 = y1 - 0.5 * scale.y
        return (y2, y2)

    # semi-private
    def x_coords(self, scale, x1):
        dx = scale.x * self.length
        x2 = x1 + dx
        return (x1, x2)

    # semi-private
    def _get_clade_lines(self, orientation='horizontal', y_curr=0, x_start=0,
                         x_curr=0, y_bot=0, y_top=0,
                         line_color='rgb(25,25,25)', line_width=0.5):
        # define a Plotly shape of type 'line', for each branch

        branch_line = dict(type='line',
                           layer='below',
                           line=dict(color=line_color,
                                     width=line_width)
                           )
        if orientation == 'horizontal':
            branch_line.update(x0=x_start,
                               y0=y_curr,
                               x1=x_curr,
                               y1=y_curr)
        elif orientation == 'vertical':
            branch_line.update(x0=x_curr,
                               y0=y_bot,
                               x1=x_curr,
                               y1=y_top)
        else:
            raise ValueError("Line type can be 'horizontal' or 'vertical'")

        return branch_line

    # semi-private
    def _draw_clade(self, clade, x_start, line_shapes,
                    line_color='rgb(15,15,15)', line_width=1):

        x_curr = self.xcoords[clade]
        y_curr = self.ycoords[clade]

        branch_line = self._get_clade_lines(orientation='horizontal',
                                            y_curr=y_curr, x_start=x_start,
                                            x_curr=x_curr,
                                            line_color=line_color,
                                            line_width=line_width)

        line_shapes.append(branch_line)

        if clade.children:
            y_top = self.ycoords[clade.children[0]]
            y_bot = self.ycoords[clade.children[-1]]

            line_shapes.append(
                self._get_clade_lines(orientation='vertical', x_curr=x_curr,
                                      y_bot=y_bot, y_top=y_top,
                                      line_color=line_color,
                                      line_width=line_width))

            for child in clade:
                self._draw_clade(child, x_curr, line_shapes)

    # semi-private
    def get_all_x_coordinates(self):
        xcoords = self.depth_dict()
        self.xcoords = xcoords

    # semi-private
    def get_all_y_coordinates(self, scale=1.3):
        maxheight = len(self.tips())

        ycoords = dict((leaf, maxheight - i * scale)
                       for i, leaf in enumerate(reversed(self.tips())))

        if self.children:
            ycoords = _calc_row(self.root(), ycoords)
        self.ycoords = ycoords

    def _set_initial_layout(self, width=800, height=800, title=None, use_lengths=None, **kw):
        super(SquareDendrogram, self)._set_initial_layout(width=width,
                                                          title=title,
                                                          height=height,
                                                          use_length=use_lengths,
                                                          **kw)

        self.get_all_x_coordinates()
        self.get_all_y_coordinates()

        my_tree_clades = list(self.xcoords.keys())
        X = []
        Y = []
        text = []
        intermediate_node_color = 'rgb(100,100,100)'
        color = [intermediate_node_color] * len(my_tree_clades)

        for index in range(len(my_tree_clades)):
            cl = my_tree_clades[index]
            X.append(self.xcoords[cl])
            Y.append(self.ycoords[cl])

            if cl.is_tip():
                color[index] = 'rgb(255, 0, 0)'
                text.append(cl.name)
            else:
                text.append(None)

        line_shapes = []
        self._draw_clade(self.root(), 0, line_shapes,
                         line_color='rgb(25,25,25)', line_width=1)

        self._trace['data'][0].update(type='scatter',
                                    x=X,
                                    y=Y,
                                    mode='markers',
                                    marker=dict(color=color, size=5),
                                    opacity=1.0,
                                    text=text,
                                    hoverinfo='text')

        self._trace['layout'].update(xaxis=dict(showline=True,
                                              zeroline=False,
                                              showgrid=False,
                                              ticklen=4,
                                              showticklabels=True,
                                              title='branch length'),
                                   shapes=line_shapes)


def _assign_ycoord(clade, node_ycoord):
    for subclade in clade.children:
        if subclade not in node_ycoord:
            node_ycoord = _assign_ycoord(subclade, node_ycoord)
    node_ycoord[clade] = 0.5 * (node_ycoord[clade.children[0]] +
                                node_ycoord[clade.children[-1]])
    return node_ycoord


def _ycoord2theta(y, start_angle, end_angle, ymin, ymax):
    return start_angle + (end_angle - start_angle) * (y - ymin) / float(
        ymax - ymin)


def _get_points_on_radial_line(x_left, x_right, y_right, start_angle, end_angle,
                               ymin, ymax):
    theta = _ycoord2theta(y_right, start_angle, end_angle, ymin,
                          ymax)
    X = [x_left * np.cos(theta), x_right * np.cos(theta), None]
    Y = [x_left * np.sin(theta), x_right * np.sin(theta), None]
    return X, Y


def _get_points_on_angular_line(x_right, y_bot, y_top, start_angle, end_angle,
                                ymin, ymax):
    theta_b = _ycoord2theta(y_bot, start_angle, end_angle, ymin,
                            ymax)
    theta_t = _ycoord2theta(y_top, start_angle, end_angle, ymin,
                            ymax)
    t = np.linspace(0, 1, 10)
    theta = (1 - t) * theta_b + t * theta_t
    X = list(x_right * np.cos(theta)) + [None]
    Y = list(x_right * np.sin(theta)) + [None]
    return X, Y


def _get_line_lists(clade, x_left, xlines, ylines, xarc, yarc, node_radius,
                    start_angle, end_angle, ymin, ymax, node_ycoord):
    """Recursively compute the lists of points that span the tree branches"""

    # xlines, ylines  - the lists of x-coords, resp y-coords of radial edge ends
    # xarc, yarc - the lists of points generating arc segments for tree branches

    x_right = node_radius[clade]
    y_right = node_ycoord[clade]

    X, Y = _get_points_on_radial_line(x_left, x_right, y_right,
                                      start_angle, end_angle,
                                      ymin, ymax)

    xlines.extend(X)
    ylines.extend(Y)

    if clade.children:

        y_top = node_ycoord[clade.children[0]]
        y_bot = node_ycoord[clade.children[-1]]

        X, Y = _get_points_on_angular_line(x_right, y_bot, y_top,
                                           start_angle, end_angle,
                                           ymin, ymax)
        xarc.extend(X)
        yarc.extend(Y)

        for child in clade.children:
            _get_line_lists(child, x_right, xlines, ylines, xarc, yarc,
                            node_radius, start_angle, end_angle, ymin, ymax,
                            node_ycoord)


class CircularDendrogram(_RootedDendrogram):
    aspect_distorts_lengths = False

    def _get_vertical_position(self, start_leaf):
        """
        returns a dict {clade: ycoord}, where y-coord is the cartesian y-coordinate
        of a  clade root in a rectangular phylogram

        """
        # n_leafs = self.tips()

        if start_leaf == 'first':
            node_ycoord = dict(
                (leaf, k) for k, leaf in enumerate(self.tips()))
        elif start_leaf == 'last':
            node_ycoord = dict(
                (leaf, k) for k, leaf in enumerate(reversed(self.tips())))
        else:
            raise ValueError("start leaf can be only 'first' or 'last'")

        if self.root().children:
            node_ycoord = _assign_ycoord(self.root(), node_ycoord)
        return node_ycoord

    def _get_circular_tree_data(self, order='level', dist=1, start_angle=0,
                                end_angle=360, start_leaf='first'):

        start_angle *= np.pi / 180
        end_angle *= np.pi / 180

        node_radius = self.depth_dict()
        node_ycoord = self._get_vertical_position(start_leaf)
        y_vals = node_ycoord.values()
        ymin, ymax = min(y_vals), max(y_vals)
        ymin -= dist

        xlines = []
        ylines = []
        xarc = []
        yarc = []
        _get_line_lists(self.root(), 0, xlines, ylines, xarc, yarc, node_radius,
                        start_angle, end_angle, ymin, ymax, node_ycoord)
        xnodes = []
        ynodes = []

        for clade in self.postorder():
            theta = _ycoord2theta(node_ycoord[clade], start_angle, end_angle,
                                  ymin, ymax)
            xnodes.append(node_radius[clade] * np.cos(theta))
            ynodes.append(node_radius[clade] * np.sin(theta))

        return xnodes, ynodes, xlines, ylines, xarc, yarc

    def _set_initial_layout(self, width=800, height=800, title=None, use_lengths=None, **kw):
        super(CircularDendrogram, self)._set_initial_layout(width=width,
                                                            height=height,
                                                            title=title,
                                                            use_length=use_lengths,
                                                            **kw)

        all_clades = list(self.postorder())
        for k in range(len(all_clades)):
            all_clades[k].id = k

        traverse_order = 'postorder'
        xnodes, ynodes, xlines, ylines, xarc, yarc = self._get_circular_tree_data(
            order=traverse_order,
            start_leaf='last')

        my_tree_clades = list(self.postorder())
        text = []  # list of text to be displayed on hover over nodes
        intermediate_node_color = 'rgb(100,100,100)'
        color = [intermediate_node_color] * len(my_tree_clades)

        for index in range(len(my_tree_clades)):
            cl = my_tree_clades[index]

            if cl.is_tip():
                color[index] = 'rgb(255, 0, 0)'
                text.append(cl.name)
            else:
                text.append(None)

        size = [9 if c != -1 else 7 for c in color]

        trace_nodes = dict(type='scatter',
                           x=xnodes,
                           y=ynodes,
                           mode='markers',
                           marker=dict(color=color,
                                       size=size,
                                       # colorscale=pl_colorscale
                                       ),
                           text=text,
                           hoverinfo='text',
                           opacity=1)

        trace_radial_lines = dict(type='scatter',
                                  x=xlines,
                                  y=ylines,
                                  mode='lines',
                                  line=dict(color='rgb(20,20,20)', width=1),
                                  hoverinfo='none')

        trace_arcs = dict(type='scatter',
                          x=xarc,
                          y=yarc,
                          mode='lines',
                          line=dict(color='rgb(20,20,20)', width=1,
                                    shape='spline'),
                          hoverinfo='none')

        self._trace['data'] = [trace_radial_lines, trace_arcs, trace_nodes]
