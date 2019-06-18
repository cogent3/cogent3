from cogent3.core.tree import TreeNode, PhyloNode
import numpy as np
from cogent3.draw.drawable import Drawable
from cogent3.util.misc import extend_docstring_from

__author__ = "Rahul Ghangas, Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rahul Ghangas"]
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
        Drawable.__init__(self, showlegend=False, visible_axes=False)
        self.length = tree.length
        self.height = None
        self.original = tree  # for edge_color_callback
        self.collapsed = False
        self.use_lengths_default = use_lengths

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
                c._update_geometry(
                    use_lengths, self.depth, track_coordinates)
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

    def _depth_dict(self):
        self._update_geometry(use_lengths=self.use_lengths_default)
        depth_dict = dict()

        for node in self.get_edge_vector():
            depth_dict[node] = node.depth

        return depth_dict


class _RootedDendrogram(_Dendrogram):
    def _width_required(self):
        return self.leafcount

    def x_coords(self, scale, x1):
        raise NotImplementedError

    def y_coords(self, scale, y1):
        raise NotImplementedError

    def _update_coordinates(self, width, height):
        xscale = width / self.height
        yscale = height / self._width_required()
        scale = Dimensions(xscale, yscale, self.height)

        # y coords done postorder, x preorder, y first.
        # so it has to be done in 2 passes.
        self._update_y_coordinates(scale)
        self._update_x_coordinates(scale)
        return xscale

    def _update_y_coordinates(self, scale, y1=None):
        """The second pass through the tree.  Y coordinates only
        depend on the shape of the tree and yscale"""
        if y1 is None:
            y1 = self._width_required() * scale.y
        child_y = y1
        for child in self.children:
            child._update_y_coordinates(scale, child_y)
            child_y -= child._width_required() * scale.y
        (self.y1, self.y2) = self.y_coords(scale, y1)

    def _update_x_coordinates(self, scale, x1=0):
        """For non 'square' styles the x coordinates will depend
        (a bit) on the y coordinates, so they should be done first"""
        (self.x1, self.x2) = self.x_coords(scale, x1)
        for child in self.children:
            child._update_x_coordinates(scale, self.x2)


class SquareDendrogram(_RootedDendrogram):
    aspect_distorts_lengths = False

    def y_coords(self, scale, y1):
        cys = [c.y1 for c in self.children]
        if cys:
            y2 = (cys[0] + cys[-1]) / 2.0
        else:
            y2 = y1 - 0.5 * scale.y
        return (y2, y2)

    def x_coords(self, scale, x1):
        dx = scale.x * self.length
        x2 = x1 + dx
        return (x1, x2)

    def _get_clade_lines(self, orientation='horizontal', y_curr=0, x_start=0,
                         x_curr=0, y_bot=0, y_top=0,
                         line_color='rgb(25,25,25)', line_width=0.5):
        # define a Plotly shape of type 'line', for each branch

        branch_line = dict(type='line',
                           layer='below',
                           line=dict(color=line_color,
                                     width=line_width))
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

    def _draw_clade(self, clade, x_start, line_shapes,
                    line_color='rgb(15,15,15)', line_width=2):

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

    # todo this is not a 'get' method, it's a set method
    def _get_all_x_coordinates(self):
        xcoords = self._depth_dict()
        self.xcoords = xcoords

    # todo this is not a 'get' method, it's a set method
    def _get_all_y_coordinates(self, scale=1.3):
        maxheight = len(self.tips())

        ycoords = dict((leaf, maxheight - i * scale)
                       for i, leaf in enumerate(reversed(self.tips())))

        if self.children:
            ycoords = _calc_row(self.root(), ycoords)
        self.ycoords = ycoords

    def _build_fig(self, width=800, height=800, title=None, use_lengths=None,
                   **kw):
        layout = {k: v for k, v in locals().items() if k != 'self' and v}
        self.layout.update(layout)
        self._get_all_x_coordinates()
        self._get_all_y_coordinates()

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
                text.append('')

        # todo refactor so modification of line_shapes is not hidden
        line_shapes = []
        self._draw_clade(self.root(), 0, line_shapes,
                         line_color='rgb(25,25,25)', line_width=2)

        trace = dict(type='scatter', x=X, y=Y, mode='markers+text',
                     marker=dict(color=color, size=5), opacity=1.0,
                     text=text, textposition='middle right', hoverinfo='skip')
        if not self.traces:
            self.traces.append(trace)
        else:
            self.traces[0].update(trace)

        layout = dict(font=dict(family='Courier New, monospace', size=10,
                                color='#000000'),
                      xaxis=dict(showline=True,
                                 zeroline=False,
                                 showgrid=False,
                                 ticklen=4,
                                 showticklabels=True,
                                 title='branch length'),
                      shapes=line_shapes)
        self.layout.update(layout)


class CircularDendrogram(_RootedDendrogram):
    aspect_distorts_lengths = False

    def _get_vertical_position(self, start_leaf):
        """
        returns a dict {clade: ycoord}, where y-coord is the cartesian y-coordinate
        of a  clade root in a rectangular phylogram

        """
        # n_leafs = self.tips()

        if start_leaf:
            node_ycoord = {leaf: k for k, leaf in enumerate(self.tips())}
        else:
            node_ycoord = {leaf: k for k, leaf in
                           enumerate(reversed(self.tips()))}

        if self.root().children:
            node_ycoord = _assign_ycoord(self.root(), node_ycoord)
        return node_ycoord

    def _get_circular_tree_data(self, dist=1, start_angle=0,
                                end_angle=360, start_leaf=True):

        start_angle *= np.pi / 180
        end_angle *= np.pi / 180

        node_radius = self._depth_dict()
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
        angles = []

        for clade in self.postorder():
            theta = _ycoord2theta(node_ycoord[clade], start_angle, end_angle,
                                  ymin, ymax)
            angles.append(theta)
            xnodes.append(node_radius[clade] * np.cos(theta))
            ynodes.append(node_radius[clade] * np.sin(theta))

        return xnodes, ynodes, xlines, ylines, xarc, yarc, angles

    def _build_fig(self, width=800, height=800, title=None,
                   use_lengths=None, **kw):
        layout = {k: v for k, v in locals().items() if k != 'self' and v}
        self.layout.update(layout)

        all_clades = list(self.postorder())
        for k in range(len(all_clades)):
            all_clades[k].id = k

        xnodes, ynodes, xlines, ylines, xarc, yarc, angles = self._get_circular_tree_data(
            start_leaf=False)

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
                text.append('')

        size = [9 if c != -1 else 7 for c in color]

        trace_nodes = dict(type='scatter',
                           x=xnodes,
                           y=ynodes,
                           mode='markers',
                           marker=dict(color=color,
                                       size=size,
                                       ),
                           hoverinfo='none',
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

        annotations = []
        for i in range(len(xnodes)):

            point = dict(x=xnodes[i] + 0.3 * np.cos(angles[i]),
                         y=ynodes[i] - 0.3 * np.sin(- angles[i]),
                         showarrow=False,
                         text=text[i],
                         textangle=-angles[i] * 180 / np.pi if np.cos(
                             angles[i]) > 0 else -180 - angles[i] * 180 / np.pi
                         )
            annotations += [point]

        layout = dict(font=dict(family='Courier New, monospace', size=10,
                                color='#000000'),
                      annotations=annotations)
        self.layout.update(layout)

        self.traces.extend([trace_radial_lines, trace_arcs, trace_nodes])


class Dimensions(object):

    def __init__(self, xscale, yscale, total_tree_height):
        self.x = xscale
        self.y = yscale
        self.height = total_tree_height


def _calc_row(clade, ycoords):
    for subclade in clade.children:
        if subclade not in ycoords:
            ycoords = _calc_row(subclade, ycoords)
    ycoords[clade] = (ycoords[clade.children[0]] +
                      ycoords[clade.children[-1]]) / 2
    return ycoords


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


class TreeGeometryBase(PhyloNode):
    """base class that computes geometric coordinates for display"""

    def __init__(self, tree, length_attr='length'):
        """
        Parameters
        ----------
        tree : either a PhyloNode or TreeNode instance
        length_attr : str
            name of the attribute to use for length, defaults to 'length'
        """
        children = [type(self)(child, length_attr=length_attr)
                    for child in tree.children]
        PhyloNode.__init__(self, params=tree.params.copy(), children=children,
                           name=tree.name)
        # todo do we need to validate the length_attr key exists?
        self._length = length_attr
        self._node_space = 1.3
        self._start = None
        self._end = None
        self._y = None
        self._x = None
        self._tip_rank = None

    def propagate_properties(self):
        self._init_tip_ranks()
        self._init_length_depth_attr()

    @property
    def x(self):
        if self.is_root():
            self._x = 0
        elif self._x is None:
            val = self.params['cum_length']
            self._x = val
        return self._x

    @property
    def tip_rank(self):
        return self._tip_rank

    def _max_child_depth(self):
        """computes the maximum number of nodes to the tip"""
        if 'max_child_depth' not in self.params:
            if self.is_tip():
                self.params['max_child_depth'] = self.params['depth']
            else:
                depths = [child._max_child_depth()
                          for child in self.children]
                self.params['max_child_depth'] = max(depths)
        return self.params['max_child_depth']

    def _init_length_depth_attr(self):
        """check it exists, if not, creates with default value of 1"""
        for edge in self.preorder():
            if edge.is_root():
                edge.params['depth'] = 0
                edge.params[self._length] = 0
            else:
                length = edge.params.get(self._length, None) or 1
                edge.params[self._length] = length
                edge.params['depth'] = edge.parent.params.get('depth', 0) + 1

        self._max_child_depth()
        for edge in self.preorder():
            if edge.is_root():
                continue

            parent_frac = edge.parent.params.get('cum_length', 0)
            if edge.is_tip():
                frac = 1 - parent_frac
            else:
                frac = 1 / edge.params['max_child_depth']
            edge.params['frac_pos'] = frac
            edge.params['cum_length'] = parent_frac + edge.params['frac_pos']

    @property
    def depth(self):
        return self.params['depth']

    def get_segment_to_parent(self, children):
        """returns (self.start, self.end)"""
        return self.start, self.end

    def value_and_coordinate(self, attr, padding=0.1):
        """
        Parameters
        ----------
        attr : str
            attribute of self, e.g. 'name', or key in self.params
        padding : float
            distance from self coordinate

        Returns
        -------
        (value of attr, (x, y)
        """
        # todo, possibly also return a rotation?
        raise NotImplementedError('implement in sub-class')


class SquareTreeGeometry(TreeGeometryBase):
    """represents Square dendrograms, contemporaneous or not"""

    def __init__(self, *args, **kwargs):
        super(SquareTreeGeometry, self).__init__(*args, **kwargs)

        if self.is_root():
            self._init_tip_ranks()

    def _init_tip_ranks(self):
        tips = self.tips()
        num_tips = len(tips)
        for index, tip in enumerate(tips):
            tip._tip_rank = index
            tip._y = ((num_tips - 1) / 2 - index) * self._node_space

    @property
    def y(self):
        if self.is_root():
            self._y = 0
        elif self._y is None:
            val = (self.children[0].y + self.children[-1].y) / 2
            self._y = val
        return self._y

    @property
    def start(self):
        """x, y coordinate for line connecting parent to this node"""
        # needs to ask parent, but has to do more than just get the parent
        # end, parent needs to know
        if self.is_root():
            val = 0, 0
        else:
            val = self.parent.x, self.y
        return val

    @property
    def end(self):
        """x, y coordinate for this node"""
        return self.x, self.y

    def get_segment_to_children(self):
        """returns coordinates connecting all children to self.end"""
        # if tip needs to
        a = self.children[0].start
        b = self.children[-1].start
        return a, b

    @extend_docstring_from(TreeGeometryBase.value_and_coordinate)
    def value_and_coordinate(self, attr, padding=0.1):
        # todo, possibly also return a rotation?
        x = self.x + padding
        y = self.y
        value = self.params.get(attr, None)
        if value is None:
            value = getattr(self, attr, None)
        return value, (x, y)


class Dendrogram(Drawable):
    def __init__(self, tree, style='square', label_pad=.1, contemporaneous=True,
                 *args, **kwargs):
        super(Dendrogram, self).__init__(
            visible_axes=False, showlegend=False, *args, **kwargs)
        klass = {'square': SquareTreeGeometry}[style]
        kwargs = dict(length_attr='frac_pos') if contemporaneous else {}
        self.tree = klass(tree, **kwargs)
        self.tree.propagate_properties()
        self.label_pad = label_pad
        self._tip_font = dict(size=16, family='sans serif')
        self._line_width = 2
        self._line_color = 'black'

    @property
    def tip_font(self):
        return self._tip_font

    def _get_tip_name_annotations(self):
        annotations = []
        for tip in self.tree.tips():
            name, (x, y) = tip.value_and_coordinate('name',
                                                    padding=self.label_pad)
            anote = dict(x=x,
                         y=y,
                         xref='x',
                         yref='y',
                         showarrow=False,
                         xanchor='left',
                         text=name,
                         font=self.tip_font)
            annotations.append(anote)
        return annotations

    def _build_fig(self, **kwargs):
        X = []
        Y = []
        tree = self.tree
        text = {'type': 'scatter', 'text': [], 'x': [], 'y': [],
                'hoverinfo': 'text', 'mode': 'markers',
                'marker': {'symbol': 'circle', 'color': 'black'}}
        for edge in tree.preorder():
            x0, y0 = edge.start
            x1, y1 = edge.end
            X.extend([x0, x1, None])
            Y.extend([y0, y1, None])
            if edge.is_tip():
                continue
            name, (x, y) = edge.value_and_coordinate('name', padding=0)
            text['x'].append(x)
            text['y'].append(y)
            text['text'].append(name)
            segment = []
            for x, y in edge.get_segment_to_children():
                segment += [(x, y)]
            xs, ys = list(zip(*segment))
            xs += (None,)
            ys += (None,)
            X.extend(xs)
            Y.extend(ys)
        trace = dict(type='scatter', x=X, y=Y, mode='lines',
                     line=dict(width=self._line_width,
                               color=self._line_color))
        self.traces.extend([trace, text])
        self.layout.annotations = tuple(self._get_tip_name_annotations())
