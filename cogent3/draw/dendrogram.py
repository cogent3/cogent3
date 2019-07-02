from collections import defaultdict
from math import floor

import numpy as np

from cogent3.core.tree import PhyloNode
from cogent3.draw.drawable import Drawable
from cogent3.util.misc import extend_docstring_from
from cogent3.util.union_dict import UnionDict


__author__ = "Rahul Ghangas, Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rahul Ghangas"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class TreeGeometryBase(PhyloNode):
    """base class that computes geometric coordinates for display"""

    def __init__(self, tree=None, length_attr="length", *args, **kwargs):
        """
        Parameters
        ----------
        tree : either a PhyloNode or TreeNode instance
        length_attr : str
            name of the attribute to use for length, defaults to 'length'
        """
        if tree is not None:
            children = [type(self)(child, *args, **kwargs) for child in tree.children]
            PhyloNode.__init__(
                self, params=tree.params.copy(), children=children, name=tree.name
            )
        else:
            PhyloNode.__init__(self, **kwargs)

        # todo do we need to validate the length_attr key exists?
        self._length = length_attr
        self._node_space = 1.3
        self._start = None
        self._end = None
        self._y = None
        self._x = None
        self._tip_rank = None
        self._max_x = 0
        self._min_x = 0
        self._max_y = 0
        self._min_y = 0
        self._theta = 0
        self._num_tips = 1

    def propagate_properties(self):
        self._init_length_depth_attr()
        self._init_tip_ranks()

    def _max_child_depth(self):
        """computes the maximum number of nodes to the tip"""
        if self.is_tip():
            self.params["max_child_depth"] = self.params["depth"]
        else:
            depths = [child._max_child_depth() for child in self.children]
            self.params["max_child_depth"] = max(depths)
        return self.params["max_child_depth"]

    def _init_tip_ranks(self):
        tips = self.tips()
        num_tips = len(tips)
        for index, tip in enumerate(tips):
            tip._tip_rank = index
            tip._y = ((num_tips - 1) / 2 - index) * self.node_space

    def _init_length_depth_attr(self):
        """check it exists, if not, creates with default value of 1"""
        # we compute cumulative lengths first with sorting, as that dictates the ordering of children
        # then determin
        for edge in self.preorder():
            # making sure children have the correct setting for ._length
            edge._length = self._length
            edge._num_tips = self._num_tips
            edge.node_space = self.node_space
            if len(edge.children) > 0:
                # we order children by their length
                children = sorted(
                    [(c.params.get(self._length, 0) or 0, c) for c in edge.children]
                )
                edge.children = [c for _, c in children]

            if edge.is_root():
                edge.params["depth"] = 0
                edge.params[self._length] = 0
                continue

            length = edge.params.get(self._length, None) or 1
            edge.params[self._length] = length
            edge.params["depth"] = edge.parent.params.get("depth", 0) + 1

        self._max_child_depth()
        for edge in self.preorder():
            if edge.is_root():
                edge.params["cum_length"] = 0
                continue

            parent_frac = edge.parent.params.get("cum_length", 0)
            if edge.is_tip():
                frac = 1 - parent_frac
            else:
                frac = 1 / edge.params["max_child_depth"]
            edge.params["frac_pos"] = frac
            edge.params["cum_length"] = parent_frac + edge.params[self._length]

    @property
    def x(self):
        if self.is_root():
            self._x = 0
        elif self._x is None:
            val = self.params["cum_length"]
            self._x = val
        return self._x

    @property
    def max_x(self):
        if not self._max_x:
            self._max_x = max(e.x for e in self.tips())
        return self._max_x

    @property
    def min_x(self):
        if not self._min_x:
            self._min_x = min([e.x for e in self.postorder()])
        return self._min_x

    @property
    def max_y(self):
        if not self._max_y:
            self._max_y = max(e.y for e in self.tips())
        return self._max_y

    @property
    def min_y(self):
        if not self._min_y:
            self._min_y = min(e.y for e in self.tips())
        return self._min_y

    @property
    def node_space(self):
        return self._node_space

    @node_space.setter
    def node_space(self, value):
        if value < 0:
            raise ValueError("node spacing must be > 0")
        self._node_space = value

    @property
    def tip_rank(self):
        return self._tip_rank

    @property
    def theta(self):
        return self._theta

    @property
    def start(self):
        """x, y coordinate for line connecting parent to this node"""
        # needs to ask parent, but has to do more than just get the parent
        # end, parent needs to know
        if self.is_root():
            val = 0, self.y
        else:
            val = self.parent.x, self.y
        return val

    @property
    def end(self):
        """x, y coordinate for this node"""
        return self.x, self.y

    @property
    def depth(self):
        return self.params["depth"]

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
        raise NotImplementedError("implement in sub-class")


class SquareTreeGeometry(TreeGeometryBase):
    """represents Square dendrograms, contemporaneous or not"""

    def __init__(self, *args, **kwargs):
        super(SquareTreeGeometry, self).__init__(*args, **kwargs)

    @property
    def y(self):
        if self._y is None:
            num_kids = len(self.children)
            even = num_kids % 2 == 0
            if even:
                i = floor(num_kids / 2)
                val = (self.children[i].y + self.children[i - 1].y) / 2
            else:
                i = floor(num_kids / 2)
                val = self.children[i].y
            self._y = val
        return self._y

    # todo should there should be a single line connecting each parent to child?
    def get_segment_to_children(self):
        """returns coordinates connecting all children to self.end"""
        # if tip needs to
        ordered = list(sorted((c.y, c) for c in self.children))
        a = ordered[0][1].start
        b = ordered[-1][1].start
        return a, b

    @extend_docstring_from(TreeGeometryBase.value_and_coordinate)
    def value_and_coordinate(self, attr, padding=0.05):
        # todo, possibly also return a rotation?
        x = self.x + padding
        y = self.y
        value = self.params.get(attr, None)
        if value is None:
            value = getattr(self, attr, None)
        data = UnionDict(
            x=x,
            y=y,
            textangle=self.theta,
            showarrow=False,
            text=self.name,
            xanchor="left",
        )

        return data


r_2_d = np.pi / 180


def polar_2_cartesian(θ, radius):
    radians = θ * r_2_d
    x = np.cos(radians) * radius
    y = np.sin(radians) * radius
    return x, y


class RadialTreeGeometry(TreeGeometryBase):
    def __init__(self, *args, **kwargs):
        super(RadialTreeGeometry, self).__init__(*args, **kwargs)
        self._num_tips = 1
        self._theta = None
        self._node_space = None

    def propagate_properties(self):
        self._num_tips = len(self.tips())
        self._init_tip_ranks()
        self._init_length_depth_attr()

    @property
    def node_space(self):
        if self._node_space is None:
            self._node_space = 360 / self._num_tips
        return self._node_space

    @node_space.setter
    def node_space(self, value):
        if value < 0:
            raise ValueError("node spacing must be > 0")
        if self._num_tips * value > 360:
            raise ValueError(f"{value} * {(self._num_tips + 1)} is > 360")
        self._node_space = value

    @property
    def theta(self):
        if self._theta is None:
            if self.is_root():
                self._theta = 0
            elif self.is_tip():
                self._theta = (self.tip_rank + 1) * self.node_space
            else:
                self._theta = sum(c.theta for c in self.children) / len(self.children)

        return self._theta

    @property
    def x(self):
        if self._x is None:
            y = self.y  # triggers populating values
        return self._x

    @property
    def y(self):
        radius = self.params["cum_length"]
        if self._y is None and self.is_root():
            self._x = self._y = 0
        elif self._x is None or self._y is None:
            x, y = polar_2_cartesian(self.theta, radius)
            self._x = x
            self._y = y

        return self._y

    @property
    def start(self):
        """x, y coordinate for line connecting parent to this node"""
        # needs to ask parent, but has to do more than just get the parent
        # end, parent needs to know
        if self.is_root():
            val = 0, 0
        else:
            # radius comes from parent
            radius = self.parent.params["cum_length"]
            val = polar_2_cartesian(self.theta, radius)
        return val

    def get_segment_to_children(self):
        """returns coordinates connecting all children to self.end"""
        # if tip needs to
        a = self.children[0].start
        b = self.end
        c = self.children[-1].start
        return a, b, c

    @extend_docstring_from(TreeGeometryBase.value_and_coordinate)
    def value_and_coordinate(self, attr, padding=0.05):
        if 90 < self.theta <= 270:
            radius = np.sqrt(self.x ** 2 + self.y ** 2) + 2 * padding
            textangle = 180 - self.theta
        else:
            radius = np.sqrt(self.x ** 2 + self.y ** 2) + padding
            textangle = 360 - self.theta

        x, y = polar_2_cartesian(self.theta, radius)

        data = UnionDict(
            x=x,
            y=y,
            textangle=textangle,
            showarrow=False,
            text=self.name,
            xanchor="center",
            yanchor="middle",
            align="left",
        )
        return data


class Dendrogram(Drawable):
    def __init__(
        self,
        tree,
        style="square",
        label_pad=0.05,
        contemporaneous=True,
        *args,
        **kwargs,
    ):
        super(Dendrogram, self).__init__(
            visible_axes=False, showlegend=False, *args, **kwargs
        )
        klass = {"square": SquareTreeGeometry, "radial": RadialTreeGeometry}[style]
        kwargs = UnionDict(length_attr="frac_pos") if contemporaneous else {}
        self.tree = klass(tree, **kwargs)
        self.tree.propagate_properties()
        assert 0 <= label_pad <= 1
        self._label_pad = label_pad
        self._tip_font = UnionDict(size=16, family="sans serif")
        self._line_width = 2
        self._line_color = "black"
        self._scale_bar = "bottom left"
        self._edge_sets = {}
        self._edge_mapping = {}
        self._contemporaneous = contemporaneous

    @property
    def label_pad(self):
        return self._label_pad

    @label_pad.setter
    def label_pad(self, value):
        if not 0 <= value <= 1:
            raise ValueError("label_pad must be in range [0, 1]")
        self._label_pad = value

    @property
    def contemporaneous(self):
        return self._contemporaneous

    @contemporaneous.setter
    def contemporaneous(self, value):
        if not type(value) == bool:
            raise TypeError
        if not self._contemporaneous == value:
            klass = self.tree.__class__
            self.tree = klass(self.tree, length_attr="frac_pos")
            self.tree.propagate_properties()
            self._traces = []
            if value:  # scale bar not needed
                self._scale_bar = False

        self._contemporaneous = value

    @property
    def tip_font(self):
        return self._tip_font

    def _scale_label_pad(self):
        """returns the label pad scaled by maximum dist to tip"""
        label_pad = self.tree.max_x * self.label_pad
        return label_pad

    def _get_tip_name_annotations(self):
        annotations = []
        label_pad = self._scale_label_pad()
        for tip in self.tree.tips():
            anote = tip.value_and_coordinate("name", padding=label_pad)
            anote |= UnionDict(xref="x", yref="y", font=self.tip_font)
            annotations.append(anote)
        return annotations

    def _get_scale_bar(self):
        if not self.scale_bar:
            return None, None

        if "left" in self.scale_bar:
            x = self.tree.min_x
        else:
            x = self.tree.max_x

        if "bottom" in self.scale_bar:
            y = self.tree.min_y
        else:
            y = self.tree.max_y

        scale = 0.1 * self.tree.max_x
        if scale < 1e-4:
            text = "{:.2e}".format(scale)
        else:
            text = "{:.2f}".format(scale)

        shape = {
            "type": "line",
            "x0": x,
            "y0": y,
            "x1": x + scale,
            "y1": y,
            "line": {"color": self._line_color, "width": self._line_width},
        }
        annotation = UnionDict(
            x=x + scale + (0.5 * scale),
            y=y,
            xref="x",
            yref="y",
            text=text,
            showarrow=False,
            ax=0,
            ay=0,
        )
        return shape, annotation

    def _build_fig(self, **kwargs):
        grouped = {}

        tree = self.tree
        text = {
            "type": "scatter",
            "text": [],
            "x": [],
            "y": [],
            "hoverinfo": "text",
            "mode": "markers",
            "marker": {"symbol": "circle", "color": "black"},
            "showlegend": False,
        }

        get_edge_group = self._edge_mapping.get
        for edge in tree.preorder():
            key = get_edge_group(edge.name, None)
            if key not in grouped:
                grouped[key] = defaultdict(list)
            group = grouped[key]
            x0, y0 = edge.start
            x1, y1 = edge.end
            group["x"].extend([x0, x1, None])
            group["y"].extend([y0, y1, None])
            if edge.is_tip():
                continue

            edge_label = edge.value_and_coordinate("name", padding=0)
            text["x"].append(edge_label.x)
            text["y"].append(edge_label.y)
            text["text"].append(edge_label.text)

            child_groups = set(get_edge_group(c.name, None) for c in edge.children)
            segment = []
            for x, y in edge.get_segment_to_children():
                segment += [(x, y)]
            xs, ys = list(zip(*segment))
            xs += (None,)
            ys += (None,)
            # todo this needs to be able to cope with children belonging to
            # 2 different groups
            if key not in child_groups:
                # different affiliation
                key = child_groups.pop()
                if key not in grouped:
                    grouped[key] = defaultdict(list)
                group = grouped[key]
            group["x"].extend(xs)
            group["y"].extend(ys)

        traces = []
        for key in grouped:
            group = grouped[key]
            style = self._edge_sets.get(
                key,
                UnionDict(
                    line=UnionDict(
                        width=self._line_width,
                        color=self._line_color,
                        shape="spline",
                        smoothing=1.3,
                    )
                ),
            )
            trace = UnionDict(type="scatter", x=group["x"], y=group["y"], mode="lines")
            trace |= style
            if "legendgroup" not in style:
                trace["showlegend"] = False
            else:
                trace["name"] = style["legendgroup"]
            traces.append(trace)

        scale_shape, scale_text = self._get_scale_bar()
        traces.extend([text])
        self.traces.extend(traces)
        self.layout.annotations = tuple(self._get_tip_name_annotations())
        if scale_shape:
            self.layout.shapes = self.layout.get("shape", []) + [scale_shape]
            self.layout.annotations += (scale_text,)
        else:
            self.layout.pop("shapes", None)

        if isinstance(self.tree, RadialTreeGeometry):
            max_x = max(self.tree.max_x, abs(self.tree.min_x)) * 1.1
            max_y = max(self.tree.max_y, abs(self.tree.min_y)) * 1.1
            # making sure the coordinates centered on the origin to avoid
            # distortion (over/under rotation of the text)
            max_x = max_y = max(max_y, max_x)

            # must draw this square
            if self.layout.width and self.layout.height:
                dim = max(self.layout.width, self.layout.height)
            elif self.layout.width:
                dim = self.layout.width
            elif self.layout.height:
                dim = self.layout.height
            else:
                dim = 800
            self.layout.width = self.layout.height = dim

            # I'm assuming we span the origin
            axes_range = dict(
                xaxis=dict(range=[-max_x, max_x]), yaxis=dict(range=[-max_y, max_y])
            )
            self.layout |= axes_range

    def style_edges(self, edges, line, legendgroup=None):
        """adjust display layout for the edges

        Parameters
        ----------
        edges : str or series
            names of edges
        line : dict
            with plotly line style to applied to these edges
        legendgroup : str or None
            if str, a legend will be presented
        """
        if type(edges) == str:
            edges = [edges]
        edges = frozenset(edges)
        style = UnionDict(width=self._line_width, color=self._line_color)
        style.update(line)
        self._edge_sets[edges] = UnionDict(legendgroup=legendgroup, line=style)
        mapping = {e: edges for e in edges}
        self._edge_mapping.update(mapping)
        if legendgroup:
            self.layout["showlegend"] = True

        # need to trigger recreation of figure
        self._traces = []

    def get_edge_names(self, tip1, tip2, outgroup=None, stem=False, clade=True):
        names = self.tree.get_edge_names(
            tip1, tip2, stem=stem, clade=clade, outgroup_name=outgroup
        )
        return names

    @property
    def scale_bar(self):
        """where to place a scale bar"""
        return self._scale_bar

    @scale_bar.setter
    def scale_bar(self, value):
        if value is True:
            value = "bottom left"

        valid = {"bottom left", "bottom right", "top left", "top right", False, None}

        assert value in valid
        if value != self._scale_bar:
            self._traces = []
        self._scale_bar = value
