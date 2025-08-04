import contextlib
from collections import defaultdict
from collections.abc import Iterable, Sequence
from math import floor
from typing import Any, Literal, cast

import numpy as np

from cogent3.core.tree import PhyloNode
from cogent3.draw.drawable import Drawable
from cogent3.util.misc import extend_docstring_from
from cogent3.util.union_dict import UnionDict


class TreeGeometryBase(PhyloNode):
    """base class that computes geometric coordinates for display"""

    def __init__(
        self,
        tree: PhyloNode | None = None,
        length_attr: str = "length",
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """
        Parameters
        ----------
        tree : PhyloNode
        length_attr : str
            name of the attribute to use for length, defaults to 'length'
        """
        if tree is not None:
            children = [type(self)(child, *args, **kwargs) for child in tree.children]
            PhyloNode.__init__(
                self,
                params=tree.params.copy(),
                children=children,
                name=tree.name,
                name_loaded=tree.name_loaded,
                length=tree.length,
                support=tree.support,
            )
        else:
            PhyloNode.__init__(self, **kwargs)

        if length_attr == "length":
            self.params["length"] = self.length

        # TODO do we need to validate the length_attr key exists?
        self._length = length_attr
        self._node_space: float | None = 1.3
        self._start = None
        self._end = None
        self._y: float | None = None
        self._x: float | None = None
        self._tip_rank: int | None = None
        self._max_x = 0.0
        self._min_x = 0.0
        self._max_y = 0.0
        self._min_y = 0.0
        self._theta: float | None = 0.0
        self._num_tips = 1

    def propagate_properties(self) -> None:
        self._init_length_depth_attr()
        self._init_tip_ranks()

    def _max_child_depth(self) -> int:
        """computes the maximum number of nodes to the tip"""
        if self.is_tip():
            self.params["max_child_depth"] = self.params["depth"]
        else:
            depths = [child._max_child_depth() for child in self.children]
            self.params["max_child_depth"] = max(depths)
        return self.params["max_child_depth"]

    def _init_tip_ranks(self) -> None:
        tips = self.tips()
        num_tips = len(tips)
        for index, tip in enumerate(tips):
            tip._tip_rank = index
            tip._y = ((num_tips - 1) / 2 - index) * self.node_space

    def _init_length_depth_attr(self) -> None:
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
                    [(c.params.get(self._length, 0) or 0, c) for c in edge.children],
                )
                edge.children = [c for _, c in children]

            if edge.is_root():
                edge.params["depth"] = 0
                edge.params[self._length] = 0
                continue

            parent = cast("Self", edge.parent)

            length = edge.params.get(self._length, None) or 1
            edge.params[self._length] = length
            edge.params["depth"] = parent.params.get("depth", 0) + 1

        self._max_child_depth()
        for edge in self.preorder():
            if edge.is_root():
                edge.params["cum_length"] = 0
                continue

            parent = cast("Self", edge.parent)

            parent_frac = parent.params.get("cum_length", 0)
            if edge.is_tip():
                frac = 1 - parent_frac
            else:
                frac = 1 / edge.params["max_child_depth"]
            edge.params["frac_pos"] = frac
            edge.params["cum_length"] = parent_frac + edge.params[self._length]

    @property
    def x(self) -> float:
        if self.is_root():
            self._x = 0
        elif self._x is None:
            val: float = self.params["cum_length"]
            self._x = val
        return self._x

    @property
    def y(self) -> float:
        if self._y is None:
            return 0
        return self._y

    @property
    def max_x(self) -> float:
        if not self._max_x:
            self._max_x = max(e.x for e in self.tips())
        return self._max_x

    @property
    def min_x(self) -> float:
        if not self._min_x:
            self._min_x = min([e.x for e in self.postorder()])
        return self._min_x

    @property
    def max_y(self) -> float:
        if not self._max_y:
            self._max_y = max(e.y for e in self.tips())
        return self._max_y

    @property
    def min_y(self) -> float:
        if not self._min_y:
            self._min_y = min(e.y for e in self.tips())
        return self._min_y

    @property
    def node_space(self) -> float:
        if self._node_space is None:
            msg = "Node space should not be None."
            raise ValueError(msg)

        return self._node_space

    @node_space.setter
    def node_space(self, value: float) -> None:
        if value < 0:
            msg = "node spacing must be > 0"
            raise ValueError(msg)
        self._node_space = value

    @property
    def tip_rank(self) -> int:
        if self._tip_rank is None:
            msg = "Tip rank should not be None"
            raise ValueError(msg)
        return self._tip_rank

    @property
    def theta(self) -> float:
        if self._theta is None:
            msg = "Theta should not be None."
            raise ValueError(msg)
        return self._theta

    @property
    def start(self) -> tuple[float, float]:
        """x, y coordinate for line connecting parent to this node"""
        # needs to ask parent, but has to do more than just get the parent
        # end, parent needs to know
        if self.parent is None:
            msg = "Node has no parent"
            raise ValueError(msg)
        return (0, self.y) if self.is_root() else (self.parent.x, self.y)

    @property
    def end(self) -> tuple[float, float]:
        """x, y coordinate for this node"""
        return self.x, self.y

    @property
    def depth(self) -> int:
        return self.params["depth"]

    def get_segment_to_parent(
        self,
    ) -> tuple[tuple[float, float] | tuple[None, None], ...]:
        """returns (self.start, self.end)"""
        if self.is_root():
            return ((None, None),)

        parent = cast("Self", self.parent)

        segment_start = parent.get_segment_to_child(self)
        result: Sequence[tuple[float, float] | tuple[None, None]]
        if isinstance(segment_start, list):
            result = [
                *segment_start,
                (None, None),
                self.start,
                self.end,
            ]
            return tuple(result)
        return segment_start, self.start, (None, None), self.start, self.end

    def get_segment_to_child(
        self, child: "TreeGeometryBase"
    ) -> tuple[float, float] | list[tuple[float, float]]:
        """returns coordinates connecting a child to self and descendants"""
        return self.end

    def value_and_coordinate(
        self, attr: str, padding: float = 0.1, max_attr_length: int | None = None
    ) -> UnionDict:
        """
        Parameters
        ----------
        attr
            attribute of self, e.g. 'name', or key in self.params
        padding
            distance from self coordinate
        max_attr_length
            maximum text length of the attribute
        Returns
        -------
        (value of attr, (x, y)
        """
        # TODO, possibly also return a rotation?
        msg = "implement in sub-class"
        raise NotImplementedError(msg)

    def support_text_coord(
        self,
        xshift: int | None,
        yshift: int | None,
        threshold: float = 100,
        max_attr_length: int | None = 4,
    ) -> UnionDict | None:
        """
        Parameters
        ----------
        xshift, yshift
            relative position (in pixels) of text
        threshold
            support values below this will be displayed
        max_attr_length
            maximum text length of the attribute

        Returns
        -------
        None if threshold not met, else params['support'] and coords
        """
        if xshift is None:
            xshift = -14

        if yshift is None:
            yshift = 7

        if self.is_tip():
            return None

        val = self.support
        if val is None or val > threshold or self.is_tip():
            return None

        if threshold > 1:
            # assuming we have support as a percentage
            text = f"{int(round(val, 0))}"
        else:
            text = f"{val:.2f}"

        x = self.x
        return UnionDict(
            x=x,
            y=self.y,
            xshift=xshift,
            yshift=yshift,
            textangle=self.theta,
            showarrow=False,
            text=text,
            xanchor="center",
        )


class SquareTreeGeometry(TreeGeometryBase):
    """represents Square dendrograms, contemporaneous or not"""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

    @property
    def y(self) -> float:
        if self._y is None:
            num_kids = len(self.children)
            even = num_kids % 2 == 0
            i = floor(num_kids / 2)
            if even:
                val = (self.children[i].y + self.children[i - 1].y) / 2
            else:
                val = self.children[i].y
            self._y = val
        return self._y

    def get_segment_to_child(self, child: TreeGeometryBase) -> tuple[float, float]:
        """returns coordinates connecting a child to self and descendants"""

        if not hasattr(self, "_ordered"):
            self._ordered = sorted(
                [(c.y, c.start) for c in self.children] + [(self.y, self.end)],
            )
        ordered = self._ordered
        dist = child.y - self.y
        if np.allclose(dist, 0):
            return self.end
        if dist < 0:
            result = ordered[ordered.index((child.y, child.start)) + 1][1]
        else:
            result = ordered[ordered.index((child.y, child.start)) - 1][1]

        return result

    @extend_docstring_from(TreeGeometryBase.value_and_coordinate)
    def value_and_coordinate(
        self,
        attr: str = "name",
        padding: float = 0.05,
        max_attr_length: int | None = None,
    ) -> UnionDict:
        # TODO, possibly also return a rotation?
        x = self.x + padding
        y = self.y
        value = self.params.get(attr, None)
        if value is None:
            value = getattr(self, attr, None)
        return UnionDict(
            x=x,
            y=y,
            textangle=self.theta,
            showarrow=False,
            text=value,
            xanchor="left",
        )


class _AngularGeometry(TreeGeometryBase):
    """directly connects child to parents"""

    @property
    def start(self) -> tuple[float, float]:
        """x, y coordinate for line connecting parent to this node"""
        if self.is_root():
            return (0, self.y)
        parent = cast("Self", self.parent)
        return parent.end


class AngularTreeGeometry(_AngularGeometry, SquareTreeGeometry):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)


r_2_d = np.pi / 180


def polar_2_cartesian(
    theta: float | np.floating, radius: float | np.floating
) -> tuple[float, float]:
    radians = theta * r_2_d
    x = float(np.cos(radians) * radius)
    y = float(np.sin(radians) * radius)
    return x, y


class CircularTreeGeometry(TreeGeometryBase):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._num_tips = 1
        self._theta: float | None = None
        self._node_space = None

    def propagate_properties(self) -> None:
        self._num_tips = len(self.tips())
        self._init_length_depth_attr()
        self._init_tip_ranks()

    @property
    def node_space(self) -> float:
        if self._node_space is None:
            self._node_space = 360 / self._num_tips
        return self._node_space

    @node_space.setter
    def node_space(self, value: float) -> None:
        if value < 0:
            msg = "node spacing must be > 0"
            raise ValueError(msg)
        if self._num_tips * value > 360:
            msg = f"{value} * {(self._num_tips + 1)} is > 360"
            raise ValueError(msg)
        self._node_space = value

    @property
    def theta(self) -> float:
        if self._theta is None:
            if self.is_root():
                self._theta = 0
            elif self.is_tip():
                self._theta = (self.tip_rank + 1) * self.node_space
            else:
                self._theta = sum(c.theta for c in self.children) / len(self.children)

        return self._theta

    @property
    def x(self) -> float:
        if self._x is None:
            _ = self.y  # triggers populating values
            self._x = cast("float", self._x)
        return self._x

    @property
    def y(self) -> float:
        radius: float = self.params["cum_length"]
        if self._y is None and self.is_root():
            self._x = self._y = 0
        elif self._x is None or self._y is None:
            x, y = polar_2_cartesian(self.theta, radius)
            self._x = x
            self._y = y

        return self._y

    @property
    def start(self) -> tuple[float, float]:
        """x, y coordinate for line connecting parent to this node"""
        # needs to ask parent, but has to do more than just get the parent
        # end, parent needs to know
        if self.is_root():
            return (0, 0)

        # radius comes from parent
        parent = cast("Self", self.parent)
        radius = parent.params["cum_length"]
        return polar_2_cartesian(self.theta, radius)

    @extend_docstring_from(TreeGeometryBase.value_and_coordinate)
    def value_and_coordinate(
        self,
        attr: str = "name",
        padding: float = 0.05,
        max_attr_length: int | None = None,
    ) -> UnionDict:
        value: str | None = self.params.get(attr, None)
        if value is None:
            value = cast("str", getattr(self, attr, None))

        max_attr_length = len(value) if max_attr_length is None else max_attr_length
        if 90 < self.theta <= 270:
            textangle = 180 - self.theta
            value = value.rjust(max_attr_length)
        else:
            textangle = 360 - self.theta
            value = value.ljust(max_attr_length)

        radius = np.sqrt(self.x**2 + self.y**2) + padding
        x, y = polar_2_cartesian(self.theta, radius)

        return UnionDict(
            x=x,
            y=y,
            textangle=textangle,
            showarrow=False,
            text=value,
            xanchor="center",
            yanchor="middle",
        )

    @extend_docstring_from(TreeGeometryBase.support_text_coord)
    def support_text_coord(
        self,
        xshift: float | None,
        yshift: float | None,
        threshold: float = 100,
        max_attr_length: int | None = 4,
    ) -> UnionDict | None:
        if xshift is None:
            xshift = -18

        if yshift is None:
            yshift = 3

        if self.is_tip():
            return None

        val = self.support
        if val is None or val > threshold or self.is_tip():
            return None

        if threshold > 1:
            # assuming we have support as a percentage
            text = f"{int(round(val, 0))}"
        else:
            text = f"{val:.2f}"

        if 90 < self.theta <= 270:
            textangle = 180 - self.theta
        else:
            textangle = 360 - self.theta
            xshift = -xshift

        c = np.cos(textangle * (np.pi / 180))
        s = np.sin(textangle * (np.pi / 180))
        m = np.array([[c, s], [-s, c]])
        d = np.dot(m, [xshift, yshift])

        new_xshift = float(d.T[0])
        new_yshift = float(d.T[1])

        x = self.x
        return UnionDict(
            x=x,
            y=self.y,
            xshift=new_xshift,
            yshift=new_yshift,
            textangle=textangle,
            showarrow=False,
            text=text,
            xanchor="center",
        )

    def get_segment_to_child(
        self, child: TreeGeometryBase
    ) -> tuple[float, float] | list[tuple[float, float]]:
        """returns coordinates connecting a child to self and descendants"""

        if not hasattr(self, "_ordered"):
            self._ordered = sorted(
                [(c.theta, c.start) for c in self.children] + [(self.theta, self.end)],
            )
        ordered = self._ordered
        dist = child.theta - self.theta
        if np.allclose(dist, 0):
            return self.end
        if dist < 0:
            neighbours = [
                ordered[ordered.index((child.theta, child.start)) + 1],
                (child.theta, child.start),
            ]
        else:
            neighbours = [
                ordered[ordered.index((child.theta, child.start)) - 1],
                (child.theta, child.start),
            ]

        neighbours = sorted(neighbours)
        result = [
            polar_2_cartesian(theta, self.params["cum_length"])
            for theta in np.arange(neighbours[0][0], neighbours[1][0], 5)
        ]
        result.append(neighbours[1][1])

        return result


class RadialTreeGeometry(_AngularGeometry, CircularTreeGeometry):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

    def get_segment_to_child(self, child: TreeGeometryBase) -> tuple[float, float]:
        """returns coordinates connecting a child to self and descendants"""
        return self.end


class Dendrogram(Drawable):
    def __init__(
        self,
        tree: PhyloNode,
        style: Literal["square", "circular", "angular", "radial"] = "square",
        label_pad: float | None = None,
        contemporaneous: bool = False,
        show_support: bool = True,
        threshold: float = 1.0,
        support_is_percent: bool = True,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        length_attr: str | None = kwargs.pop("length_attr", None)
        super().__init__(
            *args,
            visible_axes=False,  # type: ignore[misc]
            showlegend=False,  # type: ignore[misc]
            **kwargs,
        )
        style_to_geometry: dict[str, type[TreeGeometryBase]] = {
            "square": SquareTreeGeometry,
            "circular": CircularTreeGeometry,
            "angular": AngularTreeGeometry,
            "radial": RadialTreeGeometry,
        }
        klass = style_to_geometry[style]
        if length_attr is None and not contemporaneous:
            contemporaneous = tree.children[0].length is None

        length_attr = "frac_pos" if contemporaneous else length_attr or "length"
        kwargs = UnionDict(length_attr=length_attr) if contemporaneous else {}
        self.tree = klass(tree, **kwargs)
        self.tree.propagate_properties()
        self._label_pad = label_pad
        self._tip_font = UnionDict(size=12, family="Inconsolata, monospace")
        self._line_width = 1.25
        self._marker_size = 3.0
        self._line_color = "black"
        self._scale_bar: str | bool | None = "bottom left"
        self._edge_sets: dict[frozenset[str] | None, UnionDict] = {}
        self._edge_mapping: dict[str, frozenset[str]] = {}
        self._contemporaneous = contemporaneous
        self._tips_as_text = True
        self._length_attr = self.tree._length  # type: ignore[reportPrivateUsage]
        self._tip_names = tuple(e.name for e in self.tree.tips())
        self._max_label_length = max(map(len, self._tip_names))
        if all(child.support is None for child in self.tree.children):
            show_support = False
        self._show_support = show_support

        if support_is_percent and threshold <= 1.0:
            threshold *= 100

        self._threshold = threshold
        self._support_xshift: int | None = None
        self._support_yshift: int | None = None
        self._default_layout.autosize = True
        self.layout = UnionDict(self._default_layout)

    @property
    def label_pad(self) -> float:
        default = 0.15 if isinstance(self.tree, CircularTreeGeometry) else 0.025
        if self._label_pad is None:
            if not self.contemporaneous:
                max_x = max(self.tree.max_x, abs(self.tree.min_x))
                self._label_pad = max_x * default
            else:
                self._label_pad = default
        return self._label_pad

    @label_pad.setter
    def label_pad(self, value: float | None) -> None:
        self._label_pad = value
        self._traces = []

    @property
    def support_xshift(self) -> int | None:
        """relative x position (in pixels) of support text. Can be negative or positive."""
        return self._support_xshift

    @support_xshift.setter
    def support_xshift(self, value: int) -> None:
        if value == self._support_xshift:
            return
        self._support_xshift = value
        self._traces = []

    @property
    def support_yshift(self) -> int | None:
        """relative y position (in pixels) of support text. Can be negative or positive."""
        return self._support_yshift

    @support_yshift.setter
    def support_yshift(self, value: int) -> None:
        if value == self._support_yshift:
            return
        self._support_yshift = value
        self._traces = []

    @property
    def contemporaneous(self) -> bool:
        return self._contemporaneous

    @contemporaneous.setter
    def contemporaneous(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError
        if self._contemporaneous != value:
            klass = self.tree.__class__
            length_attr = "frac_pos" if value else self._length_attr
            self.tree = klass(self.tree, length_attr=length_attr)
            self.tree.propagate_properties()
            self._traces = []
            self.layout["xaxis"] |= {"range": None, "autorange": True}
            self.layout["yaxis"] |= {"range": None, "autorange": True}
            if value:  # scale bar not needed
                self._scale_bar = False

        self._contemporaneous = value

    def _update_tip_font(self, val: dict[Any, Any] | int) -> None:
        """update tip font settings"""
        if isinstance(val, dict):
            self._tip_font |= val
            return
        if isinstance(val, int):
            self._tip_font.size = val
            return

        msg = f"{type(val)} not a supported type for tip_font"
        raise TypeError(msg)

    @property
    def tip_font(self) -> dict[Any, Any]:
        return self._tip_font

    @tip_font.setter
    def tip_font(self, val: dict[Any, Any] | int) -> None:
        """update tip font settings"""
        self._update_tip_font(val)

    def _scale_label_pad(self) -> float:
        """returns the label pad scaled by maximum dist to tip"""
        return self.label_pad

    def _get_tip_name_annotations(self) -> list[UnionDict]:
        annotations: list[UnionDict] = []
        for tip in self.tree.tips():
            anote = tip.value_and_coordinate(
                "name",
                padding=self.label_pad,
                max_attr_length=self._max_label_length,
            )
            anote |= UnionDict(xref="x", yref="y", font=self.tip_font)
            annotations.append(anote)
        return annotations

    def _get_scale_bar(self) -> tuple[dict[str, Any], UnionDict] | tuple[None, None]:
        if not self.scale_bar or self.contemporaneous:
            return None, None
        self.scale_bar = cast("str", self.scale_bar)

        # place scale bar above / below dendrogram area
        y_shift = (self.tree.max_y - self.tree.min_y) / 11
        x = self.tree.min_x if "left" in self.scale_bar else self.tree.max_x
        y = (
            self.tree.min_y - y_shift
            if "bottom" in self.scale_bar
            else self.tree.max_y + y_shift
        )
        scale = 0.1 * self.tree.max_x
        text = f"{scale:.1e}" if scale < 1e-2 else f"{scale:.2f}"
        shape: dict[str, Any] = {
            "type": "line",
            "x0": x,
            "y0": y,
            "x1": x + scale,
            "y1": y,
            "line": {"color": self._line_color, "width": self._line_width},
            "name": "scale_bar",
        }
        annotation = UnionDict(
            x=x + (0.5 * scale),
            y=y,
            xref="x",
            yref="y",
            yshift=10,
            text=text,
            showarrow=False,
            ax=0,
            ay=0,
        )
        return shape, annotation

    def _build_fig(self, **kwargs: Any) -> None:
        grouped: dict[frozenset[str] | None, defaultdict[str, list[Any]]] = {}

        tree = self.tree
        text = UnionDict(
            {
                "type": "scatter",
                "text": [],
                "x": [],
                "y": [],
                "hoverinfo": "text",
                "mode": "markers",
                "marker": {
                    "symbol": "circle",
                    "color": "black",
                    "size": self._marker_size,
                },
                "showlegend": False,
            },
        )
        support_text: list[UnionDict] = []
        get_edge_group = self._edge_mapping.get
        for edge in tree.preorder():
            key = get_edge_group(edge.name, None)
            if key not in grouped:
                grouped[key] = defaultdict(list)
            group = grouped[key]
            coords = edge.get_segment_to_parent()
            xs, ys = list(zip(*coords, strict=False))
            group["x"].extend((*xs, None))
            group["y"].extend((*ys, None))

            edge_label = edge.value_and_coordinate("name", padding=0)
            text["x"].append(edge_label.x)  # type: ignore[reportUnknownMemberType]
            text["y"].append(edge_label.y)  # type: ignore[reportUnknownMemberType]
            text["text"].append(edge_label.text)  # type: ignore[reportUnknownMemberType]
            if self.show_support:
                support = edge.support_text_coord(
                    self.support_xshift,
                    self.support_yshift,
                    threshold=self.support_threshold,
                )
                if support is not None:
                    support |= UnionDict(xref="x", yref="y", font=self.tip_font)
                    support_text.append(support)

        traces: list[UnionDict] = []
        for key, group in grouped.items():
            style = self._edge_sets.get(
                key,
                UnionDict(
                    line=UnionDict(
                        width=self._line_width,
                        color=self._line_color,
                        shape="spline",
                        smoothing=1.3,
                    ),
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
        if self.tips_as_text:
            self.layout.annotations = tuple(self._get_tip_name_annotations())

        if self.show_support and support_text:
            self.layout.annotations = self.layout.annotations + tuple(support_text)  # type: ignore[reportUnknownMemberType]

        if scale_shape:
            self.layout.shapes = [*self.layout.get("shape", []), scale_shape]  # type: ignore[reportUnknownMemberType]
            self.layout.annotations += (scale_text,)  # type: ignore[reportUnknownMemberType]
        else:
            self.layout.pop("shapes", None)  # type: ignore[reportUnknownMemberType]

        # TODO: look into TypedDict as all the type: ignores are getting a bit silly
        if isinstance(self.tree, CircularTreeGeometry):
            # must draw this square
            if self.layout.width and self.layout.height:  # type: ignore[reportUnknownMemberType]
                dim = max(self.layout.width, self.layout.height)  # type: ignore[reportUnknownMemberType]
            elif self.layout.width:  # type: ignore[reportUnknownMemberType]
                dim = self.layout.width  # type: ignore[reportUnknownMemberType]
            elif self.layout.height:  # type: ignore[reportUnknownMemberType]
                dim = self.layout.height  # type: ignore[reportUnknownMemberType]
            else:
                dim = 800
            self.layout.width = self.layout.height = dim

            # Span of tree along x-axis and Span of tree along y-axis
            x_diff = self.tree.max_x - self.tree.min_x
            y_diff = self.tree.max_y - self.tree.min_y

            # Maximum span
            max_span = max(x_diff, y_diff)

            # Use maximum span along both axes and pad the smaller one accordingly
            axes_range = {
                "xaxis": {
                    "range": [
                        self.tree.min_x - (1.4 * max_span - x_diff) / 2,
                        self.tree.max_x + (1.4 * max_span - x_diff) / 2,
                    ],
                },
                "yaxis": {
                    "range": [
                        self.tree.min_y - (1.4 * max_span - y_diff) / 2,
                        self.tree.max_y + (1.4 * max_span - y_diff) / 2,
                    ],
                },
            }
            self.layout |= axes_range

    def style_edges(
        self,
        edges: str | Iterable[str],
        line: dict[Any, Any],
        legendgroup: str | None = None,
        tip2: str | None = None,
        **kwargs: Any,
    ) -> None:
        """adjust display layout for the edges

        Parameters
        ----------
        edges : str or series
            names of edges
        line : dict
            with plotly line style to applied to these edges
        legendgroup : str or None
            if str, a legend will be presented
        tip2 : str
            if provided, and edges is a str, passes edges (as tip1) and kwargs to get_edge_names
        kwargs
            keyword arguments passed onto get_edge_names
        """
        if tip2:
            if not isinstance(edges, str):
                msg = "cannot use a series of edges and tip2"
                raise TypeError(msg)
            edges = self.get_edge_names(edges, tip2, **kwargs)

        if isinstance(edges, str):
            edges = [edges]
        edges = frozenset(edges)
        if not edges.issubset({edge.name for edge in self.tree.preorder()}):
            msg = "edge not present in tree"
            raise ValueError(msg)
        style = UnionDict(width=self._line_width, color=self._line_color)
        style.update(line)  # type: ignore[reportUnknownMemberType]
        self._edge_sets[edges] = UnionDict(legendgroup=legendgroup, line=style)
        mapping = dict.fromkeys(edges, edges)
        self._edge_mapping.update(mapping)
        if legendgroup:
            self.layout["showlegend"] = True

        # need to trigger recreation of figure
        self._traces = []

    def reorient(self, name: str, tip2: str | None = None, **kwargs: Any) -> None:
        """change orientation of tree
        Parameters
        ----------
        name : str
            name of an edge in the tree. If name is a tip, its parent becomes
            the new root, otherwise the edge becomes the root.
        tip2 : str
            if provided, passes name (as tip1) and all other args to get_edge_names,
            but sets clade=False and stem=True
        kwargs
            keyword arguments passed onto get_edge_names
        """
        if tip2:
            kwargs.update({"stem": True, "clade": False})
            edges = self.get_edge_names(name, tip2, **kwargs)
            name = edges[0]

        if name in self._tip_names:
            self.tree = self.tree.rooted_with_tip(name)
        else:
            self.tree = self.tree.rooted_at(name)

        self.tree.propagate_properties()
        self._traces = []

    def get_edge_names(
        self,
        tip1: str,
        tip2: str,
        outgroup: str | None = None,
        stem: bool = False,
        clade: bool = True,
    ) -> list[str]:
        """

        Parameters
        ----------
        tip1 : str
            name of tip 1
        tip2 : str
            name of tip 1
        outgroup : str
            name of tip outside clade of interest
        stem : bool
            include name of stem to clade defined by tip1, tip2, outgroup
        clade : bool
            include names of edges within clade defined by tip1, tip2, outgroup

        Returns
        -------
        list of edge names
        """
        return self.tree.get_edge_names(
            tip1,
            tip2,
            stem=stem,
            clade=clade,
            outgroup_name=outgroup,
        )

    @property
    def scale_bar(self) -> str | bool | None:
        """where to place a scale bar"""
        return self._scale_bar

    @scale_bar.setter
    def scale_bar(self, value: str | bool | None) -> None:
        if value is True:
            value = "bottom left"

        valid: set[str | bool | None] = {
            "bottom left",
            "bottom right",
            "top left",
            "top right",
            False,
            None,
        }

        if value not in valid:
            msg = "Invalid scale bar."
            raise ValueError(msg)
        if value != self._scale_bar:
            self._traces = []
        self._scale_bar = value

    @property
    def tips_as_text(self) -> bool:
        """displays tips as text"""
        return self._tips_as_text

    @tips_as_text.setter
    def tips_as_text(self, value: bool) -> None:
        if value == self._tips_as_text:
            return

        self._tips_as_text = value
        self._traces = []
        self.layout.annotations = ()

    @property
    def line_width(self) -> float:
        """width of dendrogram lines"""
        return self._line_width

    @line_width.setter
    def line_width(self, width: float) -> None:
        self._line_width = width
        if self.traces:
            setting = {"width": width}
            for trace in self.traces:
                with contextlib.suppress(KeyError):
                    trace["line"] |= setting

    @property
    def marker(self) -> float:
        return self._marker_size

    @marker.setter
    def marker(self, size: float) -> None:
        self._marker_size = size
        if self.traces:
            setting = {"size": size}
            for trace in self.traces:
                if trace.get("mode", None) == "markers":  # type: ignore[reportUnknownMemberType]
                    trace["marker"] |= setting

    @property
    def show_support(self) -> bool:
        """whether tree edge support entries are displayed"""
        return self._show_support

    @show_support.setter
    def show_support(self, value: bool) -> None:
        """whether tree edge support entries are displayed"""
        if value == self._show_support:
            return

        self._show_support = value
        self._traces = []
        self.layout.annotations = ()

    @property
    def support_threshold(self) -> float:
        """cutoff for dislaying support"""
        return self._threshold

    @support_threshold.setter
    def support_threshold(self, value: float) -> None:
        self._threshold = value
        self._traces = []
        self.layout.annotations = ()
