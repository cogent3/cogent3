import numpy

from numpy.testing import assert_allclose

from cogent3.draw.drawable import Drawable
from cogent3.maths.geometry import SimplexTransform
from cogent3.util.union_dict import UnionDict


__author__ = "Rahul Ghangas and Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rahul Ghangas", "Helmut Simon"]
__license__ = "BSD-3"
__version__ = "2019.10.17a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class Simplex(Drawable):
    def __init__(self, vertex_labels=None, **kwargs):
        super(Simplex, self).__init__(**kwargs)
        self.vertex_labels = vertex_labels
        self._transformer = None
        self._vertices = None
        self._groups = set()
        self._used_groups = set()

        # points
        self._points = []
        self._point_hovertext = []
        self._point_legendgroups = []

        # segments
        self._segments = []
        self._segment_hovertext = []
        self._segment_legendgroups = []

    @property
    def transformer(self):
        if self._transformer is None:
            self._transformer = SimplexTransform()
        return self._transformer

    @property
    def vertices(self):
        # todo when we get autoscaled, we need to combine both points and
        # segments
        if self._vertices is None:
            self._vertices = self.transformer.q
        return self._vertices

    def add_point(self, probs, hovertext=None, legendgroup=None):
        """add single point in the simplex
        Parameters
        ----------
        probs : dict or DictArray
            a probability vector to be displayed as a point in the simplex
        hovertext
            a corresponding annotation to be associated as hovertext
        legendgroup
            group to which the point belongs
        """
        if self.vertex_labels is None:
            labels = list(sorted(probs.keys()))
            assert len(labels) == 4, "4-state systems only"
            self.vertex_labels = labels

        probs = [probs[k] for k in self.vertex_labels]
        assert_allclose(sum(probs), 1.0)
        self._points.append(probs)
        self._point_hovertext.append(hovertext)
        self._point_legendgroups.append(legendgroup)
        self._groups.add(legendgroup)

    def add_points(self, probs, hovertext=None, legendgroup=None):
        """adds series of points via successive calls to add_points
        Parameters
        ----------
        probs : a series of dict or DictArray objects
            each element is a probability vector for display
        hovertext : series or None
            an equal length series to probs
        legendgroup : series or None
            an equal length series to probs
        """
        for i in range(len(probs)):
            p = probs[i]
            h = None if hovertext is None else hovertext[i]
            l = None if legendgroup is None else legendgroup[i]
            self.add_point(p, hovertext=h, legendgroup=l)

    def add_segment(self, segment, hovertext=None, legendgroup=None):
        """add series of points connected by a line
        Parameters
        ----------
        segment : series of dict or DictArray
            each element of segment is a probability vector that can be
            displayed as a point in the simplex
        hovertext
            a corresponding annotation to be associated as hovertext
        legendgroup
            group to which the segment belongs
        """
        assert segment[0], "must provide valid data"
        if self.vertex_labels is None:
            labels = list(sorted(segment[0].keys()))
            assert len(labels) == 4, "4-state systems only"
            self.vertex_labels = labels

        probs = []
        for point in segment:
            point = [point[k] for k in self.vertex_labels]
            assert_allclose(sum(point), 1.0)
            probs.append(point)

        self._segments.append(probs)
        self._segment_hovertext.append(hovertext)
        self._segment_legendgroups.append(legendgroup)
        self._groups.add(legendgroup)

    def _get_frame_trace(self):
        from itertools import combinations

        combos = numpy.array(list(combinations(self.vertices, 2)))
        trace = UnionDict(
            type="scatter3d",
            # Draw the edges of the Tetrahedron
            name="frame",
            x=combos[:, :, 0].ravel(),
            y=combos[:, :, 1].ravel(),
            z=combos[:, :, 2].ravel(),
            marker=UnionDict(size=4, color="#1f77b4", colorscale="Viridis"),
            line=UnionDict(color="#1f77b4", width=5),
            mode="lines",
            hoverinfo="skip",
            showlegend=False,
        )
        return trace

    def _get_vertex_label_trace(self):
        trace = UnionDict(
            type="scatter3d",
            # Draw the vertex labels
            x=self.vertices[:, 0],
            y=self.vertices[:, 1],
            z=self.vertices[:, 2],
            marker=UnionDict(size=4, color="#1f77b4", colorscale="Viridis"),
            text=self.vertex_labels,
            textfont=UnionDict(size=16, family="sans serif"),
            mode="markers+text",
            hoverinfo="skip",
            showlegend=False,
            name="labels",
        )
        return trace

    def _get_3d_scatter(
        self,
        data,
        mode,
        name=None,
        legendgroup=None,
        text=None,
        showlegend=None,
        line=None,
    ):
        data = numpy.array(data, dtype=float)
        points = numpy.array([v @ self.transformer for v in data])
        scatter = UnionDict(
            type="scatter3d",
            x=points[:, 0],
            y=points[:, 1],
            z=points[:, 2],
            marker=UnionDict(size=4, colorscale="Viridis"),
            showlegend=showlegend,
            mode=mode,
            name=name,
            text=text,
            opacity=0.75,
            legendgroup=legendgroup,
            line=line,
        )
        return scatter

    def _get_point_traces(self):
        """returns scatter 3D for points"""
        data = numpy.array(self._points, dtype=float)
        if any(self._point_hovertext):
            hovertext = numpy.array(self._point_hovertext, dtype="O")
        else:
            hovertext = None
        groups = set(self._point_legendgroups)
        multigroup = len(groups) > 1
        legendgroups = numpy.array(self._point_legendgroups, dtype="O")
        traces = []
        for group in groups:
            name = None
            if multigroup and group is None:
                name = "Other"
                showlegend = True
            elif not multigroup:
                name = None
                showlegend = False

            self._used_groups.add(group)
            indices = legendgroups == group
            group_data = data[indices, :]
            if hovertext is not None:
                group_text = hovertext[indices]
            else:
                group_text = None
            trace = self._get_3d_scatter(
                group_data,
                mode="markers",
                name=name,
                legendgroup=group,
                text=group_text,
                showlegend=showlegend,
            )
            traces.append(trace)

        return traces

    def _get_segment_traces(self):
        """returns scatter 3D for segments"""
        multigroup = len(self._groups) > 1
        traces = []
        for i, segment in enumerate(self._segments):
            group = self._segment_legendgroups[i]
            name = None
            if multigroup and group is None:
                name = "Other"
                showlegend = True
            elif not multigroup:
                name = None
                showlegend = False

            self._used_groups.add(group)
            data = numpy.array(segment, dtype=float)
            traces.append(
                self._get_3d_scatter(
                    data,
                    "lines+markers",
                    name=name,
                    legendgroup=group,
                    showlegend=showlegend,
                    line=UnionDict(width=3),
                )
            )
        return traces

    def _build_fig(self, **kwargs):
        self.traces.extend([self._get_frame_trace(), self._get_vertex_label_trace()])
        if self._points:
            self.traces.extend(self._get_point_traces())

        if self._segments:
            self.traces.extend(self._get_segment_traces())

        # Layout attributes for plotly
        axis_range = [self.vertices.min(), self.vertices.max()]
        layout = UnionDict(
            scene=UnionDict(
                xaxis=UnionDict(title="x", visible=False, range=axis_range),
                yaxis=UnionDict(title="y", visible=False, range=axis_range),
                zaxis=UnionDict(title="z", visible=False, range=axis_range),
            ),
            autosize=False,
        )
        self.layout |= layout
