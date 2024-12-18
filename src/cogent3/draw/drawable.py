import os
import pathlib

import numpy

from cogent3.util.misc import extend_docstring_from
from cogent3.util.union_dict import UnionDict

# user specified environment variable for plotly renderer
PLOTLY_RENDERER = os.environ.get("PLOTLY_RENDERER", None)


def get_domain(total, element, is_y, space=0.01):
    """returns evenly spaced domain for an element in a grid plot

    Parameters
    ----------
    total : int
        the total number of elements on the axis
    element : int
        the element number to compute the domain for
    is_y : bool
        if True, this is for a y-coordinate domain. This is reversed
        so the result is in cartesian, not array, coordinates
    space : float
        the separation between elements
    """
    if total == 1:
        return [0, 1]

    if element > total - 1:
        msg = f"{element} index too big for {total}"
        raise ValueError(msg)

    per_element = 1 / total
    space = min(space / 2, per_element / 10)
    bounds = [per_element * i for i in range(total + 1)]
    domains = [
        (bounds[k] + space, bounds[k + 1] - space) for k in range(len(bounds) - 1)
    ]
    if is_y:
        element = total - element - 1

    return domains[element]


def _show_(cls, renderer=None, **kwargs) -> None:
    """display figure

    Parameters
    ----------
    renderer : str
        names of renderers that control ability for display. If not specified,
        looks for PLOTLY_RENDERER environment variable, otherwise defaults to
        'notebook_connected+plotly_mimetype'. This setting supports display in
        JupyterLab and Jupyter Notebook, while keeping notebook size small (relies
        on internet connection for getting the javascript). See
        help(plotly.io.renderer) for more options.
    kwargs
        other arguments for plotly.io.show
    """
    from plotly.io import show

    if renderer is None and PLOTLY_RENDERER is None:
        renderer = "notebook_connected+plotly_mimetype"
    elif renderer is None:
        renderer = PLOTLY_RENDERER

    kwargs["renderer"] = renderer
    drawable = getattr(cls, "drawable", None) or cls
    fig = getattr(drawable, "figure", None)
    if fig is None:
        msg = f"{cls} does not have a drawable or figure attribute"
        raise TypeError(msg)

    width = kwargs.get("width", fig.layout.width)
    height = kwargs.get("height", fig.layout.height)
    kwargs["width"] = fig.layout.width = width
    kwargs["height"] = fig.layout.height = height
    show(fig, **kwargs)


def _iplot_(cls, width=None, height=None) -> None:
    from plotly.offline import iplot as _iplot

    layout = {}
    if width:
        layout["width"] = width
    if height:
        layout["height"] = height
    if layout:
        cls.drawable.layout |= dict(layout)
    _iplot(cls.drawable.figure)


def bind_drawable(obj, drawable):
    """binds drawable"""
    from types import MethodType

    obj.drawable = drawable
    obj.iplot = MethodType(_iplot_, obj)
    obj.show = MethodType(_show_, obj)
    return obj


class Drawable:
    """container object for Plotly figures"""

    def __init__(
        self,
        title=None,
        traces=None,
        width=None,
        height=None,
        showlegend=True,
        visible_axes=True,
        layout=None,
        xtitle=None,
        ytitle=None,
    ) -> None:
        if traces is None:
            self._traces = []
        else:
            try:
                self._traces = [UnionDict(trace) for trace in traces]
            except ValueError:
                msg = f"expected a series of dicts, got {traces}"
                raise TypeError(msg)
        title = title if title is None else {"text": title}
        self._default_layout = UnionDict(
            font={"family": "Balto", "size": 14},
            autosize=False,
            hovermode="closest",
            template=None,
            plot_bgcolor=None,
            margin={"l": 50, "r": 50, "t": 50, "b": 50, "pad": 4},
            xaxis={"visible": visible_axes},
            yaxis={"visible": visible_axes},
            title=title,
            width=width,
            height=height,
            showlegend=showlegend,
        )
        layout = layout or {}
        self.layout = UnionDict(self._default_layout)
        self.layout |= layout
        # constructor layout value over-rides
        overrides = UnionDict(
            title=title,
            width=width,
            height=height,
            showlegend=showlegend,
            xaxis={"visible": visible_axes},
            yaxis={"visible": visible_axes},
        )
        self.layout |= overrides
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.title = title

    def _repr_html_(self) -> None:
        self.show()

    @property
    def traces(self):
        return self._traces

    def add_trace(self, trace) -> None:
        self.traces.append(UnionDict(trace))

    def bound_to(self, obj):
        """returns obj with self bound to it"""
        return bind_drawable(obj, self)

    @property
    def figure(self):
        if not self.traces and hasattr(self, "_build_fig"):
            self._build_fig()

        traces = self.traces or [{}]

        xtitle = self.xtitle or self.layout.xaxis.get("title", None)
        ytitle = self.ytitle or self.layout.yaxis.get("title", None)
        self.layout.xaxis.title = xtitle
        self.layout.yaxis.title = ytitle
        return UnionDict(data=traces, layout=self.layout)

    @property
    def plotly_figure(self):
        """returns a plotly graph object"""
        from plotly.graph_objects import Figure

        return Figure(**self.figure)

    @extend_docstring_from(_show_)
    def show(self, renderer=None, **kwargs) -> None:
        _show_(self, renderer, **kwargs)

    def write(self, path, **kwargs) -> None:
        """writes static image file, suffix dictates format"""
        from plotly.io import write_image

        fig = self.figure
        kwargs["width"] = kwargs.get("width", fig.layout.width)
        kwargs["height"] = kwargs.get("height", fig.layout.height)

        path = pathlib.Path(path).expanduser().absolute()
        write_image(fig, str(path), **kwargs)

    def to_image(self, format="png", **kwargs):
        """creates static image, suffix dictates format"""
        from plotly.io import to_image

        fig = self.figure
        kwargs["width"] = kwargs.get("width", fig.layout.width)
        kwargs["height"] = kwargs.get("height", fig.layout.height)

        return to_image(fig, format=format, **kwargs)

    @property
    def fig_width(self):
        """figure width, also settable via .layout.width"""
        return self.layout.width

    @fig_width.setter
    def fig_width(self, width) -> None:
        self.layout.width = width

    @property
    def fig_height(self):
        """figure height, also settable via .layout.height"""
        return self.layout.height

    @fig_height.setter
    def fig_height(self, height) -> None:
        self.layout.height = height


_ticks_off = (
    ("showticklabels", False),
    ("mirror", True),
    ("showgrid", False),
    ("showline", True),
    ("ticks", ""),
)

_ticks_on = (
    ("showticklabels", True),
    ("mirror", True),
    ("showgrid", False),
    ("showline", True),
)


class AnnotatedDrawable(Drawable):
    """supports drawing with left and bottom tracks of annotations"""

    def __init__(
        self,
        core,
        left_track=None,
        bottom_track=None,
        xtitle=None,
        ytitle=None,
        title=None,
        xrange=None,
        yrange=None,
        width=500,
        height=500,
        layout=None,
    ) -> None:
        super().__init__(
            visible_axes=True,
            showlegend=True,
            width=width,
            height=height,
            layout=layout,
            xtitle=xtitle,
            ytitle=ytitle,
        )
        self.yrange = yrange
        self.xrange = xrange
        self._overlaying = False

        core.title = title or core.title
        self.core = core
        self.left_track = left_track
        self.bottom_track = bottom_track

    def _build_fig(self, xaxis="x", yaxis="y"):
        f = self.core.figure
        try:
            if self.layout.yaxis2.overlaying != "free":
                self._overlaying = True
        except AttributeError:
            pass

        traces = f.data
        self.layout |= dict(f.layout)
        for trace in traces:
            trace.xaxis = xaxis
            if self._overlaying and "yaxis" in trace:
                trace.yaxis = "y3"
            else:
                trace.yaxis = yaxis

        self._traces = traces
        ticks_on = dict(_ticks_on)
        f.layout.xaxis.title = self.xtitle
        f.layout.yaxis.title = self.ytitle
        f.layout.xaxis |= ticks_on
        f.layout.yaxis |= ticks_on
        return f

    def _build_2x2_fig(self):
        if not self.traces:
            _ = self._build_fig(xaxis="x2", yaxis="y2")

        layout = UnionDict(
            {
                "xaxis": {"anchor": "y", "domain": [0.0, 0.099]},
                "xaxis2": {"anchor": "y2", "domain": [0.109, 1.0]},
                "xaxis3": {"anchor": "y3", "domain": [0.109, 1.0]},
                "yaxis": {"anchor": "x", "domain": [0.109, 1.0]},
                "yaxis2": {"anchor": "x2", "domain": [0.109, 1.0]},
                "yaxis3": {"anchor": "x3", "domain": [0.0, 0.099]},
            },
        )
        layout |= self.layout
        fig = UnionDict(data=[], layout=layout)

        # common settings
        ticks_off_kwargs = dict(_ticks_off)
        ticks_on_kwargs = dict(_ticks_on)

        # core traces and layout
        fig.data.extend(self.traces)

        fig.layout.xaxis2 |= dict(range=self.xrange, **ticks_off_kwargs)
        fig.layout.yaxis2 |= dict(range=self.yrange, **ticks_off_kwargs)

        # left_track traces
        seen_types = set()
        max_x = 0
        traces = []
        for trace in self.left_track.traces:
            traces.append(trace)
            # convert to numpy array to handle None's
            x = numpy.array(trace.x, dtype=float)
            indices = numpy.logical_not(numpy.isnan(x))
            max_x = max(x[indices].max(), max_x)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        left_range = [0, int(max_x) + 1]

        # bottom_track traces
        max_y = 0
        for trace in self.bottom_track.traces:
            trace.xaxis = "x3"
            trace.yaxis = "y3"
            traces.append(trace)
            # convert to numpy array to handle None's
            y = numpy.array(trace.y, dtype=float)
            indices = numpy.logical_not(numpy.isnan(y))
            max_y = max(y[indices].max(), max_y)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        bottom_range = [0, int(max_y) + 1]

        # add all traces
        fig.data.extend(traces)
        # configure axes for titles, limits, border and ticks
        fig.layout.yaxis |= dict(
            title={"text": self.ytitle},
            range=self.yrange,
            **ticks_on_kwargs,
        )

        fig.layout.xaxis3 |= dict(
            title={"text": self.xtitle},
            range=self.xrange,
            **ticks_on_kwargs,
        )

        # adjust row width of left plot for number of feature tracks
        min_range = min(left_range[1], bottom_range[1])
        left_prop = left_range[1] / min_range

        # first the top row
        xaxis_domain = list(layout.xaxis.domain)
        xaxis_domain[1] = left_prop * xaxis_domain[1]
        fig.layout.xaxis |= dict(
            title=None,
            range=left_range,
            domain=xaxis_domain,
            **ticks_off_kwargs,
        )
        fig.layout.xaxis |= dict(
            title={},
            range=left_range,
            domain=xaxis_domain,
            **ticks_off_kwargs,
        )

        space = 0.01
        fig.layout.xaxis2.domain = (xaxis_domain[1] + space, 1.0)
        fig.layout.xaxis3.domain = (xaxis_domain[1] + space, 1.0)

        # now the right column
        bottom_prop = bottom_range[1] / min_range
        yaxis_domain = list(layout.yaxis3.domain)
        yaxis_domain[1] = bottom_prop * yaxis_domain[1]
        fig.layout.yaxis3 |= dict(
            title={},
            range=bottom_range,
            domain=yaxis_domain,
            **ticks_off_kwargs,
        )

        # and bottom of the boxes above
        fig.layout.yaxis.domain = (yaxis_domain[1] + space, 1.0)
        fig.layout.yaxis2.domain = (yaxis_domain[1] + space, 1.0)

        return fig

    def _build_2x1_fig(self):
        """2 rows, one column, dotplot and seq1 annotated"""
        if not self.traces:
            _ = self._build_fig()

        layout = UnionDict(
            xaxis={"anchor": "y2", "domain": [0.0, 1.0]},
            yaxis={"anchor": "free", "domain": [0.1135, 1.0], "position": 0.0},
            yaxis2={"anchor": "x", "domain": [0.0, 0.0985]},
        )
        if self._overlaying:
            self.layout.yaxis3 = self.layout.yaxis2
            self.layout.yaxis2 = {}
            self.layout.legend.x = 1.3
        layout |= dict(self.layout)
        fig = UnionDict(data=[], layout=layout)

        # common settings
        ticks_off_kwargs = dict(_ticks_off)
        ticks_on_kwargs = dict(_ticks_on)

        # core traces and layout
        fig.data.extend(self.traces)

        fig.layout.xaxis |= dict(
            title={"text": self.xtitle},
            range=self.xrange,
            **ticks_on_kwargs,
        )
        fig.layout.yaxis |= dict(
            title={"text": self.ytitle},
            range=self.yrange,
            **ticks_on_kwargs,
        )

        # bottom traces
        seen_types = set()
        max_y = 0
        traces = []
        for trace in self.bottom_track.traces:
            trace.yaxis = "y2"
            trace.xaxis = "x"
            traces.append(trace)
            y = numpy.array(trace.y, dtype=float)
            indices = numpy.logical_not(numpy.isnan(y))
            max_y = max(y[indices].max(), max_y)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        fig.data.extend(traces)
        fig.layout.yaxis2 |= dict(
            title={},
            range=[0, int(max_y) + 1],
            **ticks_off_kwargs,
        )
        return fig

    def _build_1x2_fig(self):
        if not self.traces:
            self._build_fig(xaxis="x2")
        layout = UnionDict(
            xaxis={"anchor": "y", "domain": [0.0, 0.099]},
            xaxis2={"anchor": "free", "domain": [0.109, 1.0], "position": 0.0},
            yaxis={"anchor": "x", "domain": [0.0, 1.0]},
        )

        layout |= self.layout
        fig = UnionDict(data=[], layout=layout)

        # common settings
        ticks_off_kwargs = dict(_ticks_off)
        ticks_on_kwargs = dict(_ticks_on)

        # core traces and layout
        fig.data.extend(self.traces)

        fig.layout.xaxis2 |= dict(
            title=self.xtitle,
            range=self.xrange,
            **ticks_on_kwargs,
        )
        fig.layout.yaxis |= dict(
            title=self.ytitle,
            range=self.yrange,
            **ticks_on_kwargs,
        )

        # left track
        seen_types = set()
        max_x = 0
        traces = []
        for trace in self.left_track.traces:
            trace.yaxis = "y"
            traces.append(trace)
            x = numpy.array(trace.x, dtype=float)
            indices = numpy.logical_not(numpy.isnan(x))
            max_x = max(x[indices].max(), max_x)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        fig.data.extend(traces)
        fig.layout.xaxis |= dict(
            title=None,
            range=[0, int(max_x) + 1],
            **ticks_off_kwargs,
        )
        return fig

    @property
    def figure(self):
        if self.bottom_track and self.left_track:
            func = self._build_2x2_fig
        elif self.bottom_track:
            func = self._build_2x1_fig
        elif self.left_track:
            func = self._build_1x2_fig
        else:
            func = self._build_fig

        return func()

    def remove_track(self, left_track=False, bottom_track=False) -> None:
        """
        Parameters
        ----------
        left_track : bool
            the left track is removed
        bottom_track : bool
            the bottom track is removed
        """
        if left_track:
            self.left_track = None

        if bottom_track:
            self.bottom_track = None

        if left_track or bottom_track:
            self.core._traces = []
            self._traces = []


class Shape:
    _mode = "lines"

    def __init__(
        self,
        name=None,
        text=None,
        filled=True,
        legendgroup=None,
        showlegend=True,
        hoverinfo=None,
        fillcolor=None,
        **kwargs,
    ) -> None:
        self.filled = filled
        self.fillcolor = fillcolor
        self._legendgroup = legendgroup
        self._showlegend = showlegend
        self.name = name
        self.text = text
        self._hoverinfo = hoverinfo or name

    def shift(self, x=0, y=0):
        if not isinstance(self.x, numpy.ndarray):
            self.x += x
            self.y += y
        else:
            self.x[self.x != None] += x  # noqa
            self.y[self.y != None] += y  # noqa

        return self

    @property
    def height(self):
        return self.top - self.bottom

    @property
    def top(self):
        if not isinstance(self.y, numpy.ndarray):
            return numpy.max(self.y)
        return numpy.max(self.y[self.y != None])  # noqa

    @property
    def bottom(self):
        if not isinstance(self.y, numpy.ndarray):
            return numpy.min(self.y)
        return numpy.min(self.y[self.y != None])  # noqa

    @property
    def middle(self):
        return self.height / 2 + self.bottom

    @property
    def T(self):
        self.x, self.y = self.y, self.x
        return self

    def as_trace(self, name=None):
        """returns component for plotly display"""
        name = name or self.name
        return UnionDict(
            type="scatter",
            x=self.x,
            y=self.y,
            mode=self._mode,
            fill="toself",
            fillcolor=self.fillcolor,
            line={"color": self.fillcolor},
            text=self.text,
            name=name,
            legendgroup=self._legendgroup,
            showlegend=self._showlegend,
            hoverinfo="text",
        )


class Rectangle(Shape):
    def __init__(self, coords, y=0, height=0.25, **kwargs) -> None:
        super().__init__(**kwargs)
        width = abs(coords[0][0] - coords[0][1])
        x_coord = min(coords[0][0], coords[0][1])
        xs = [x_coord, x_coord, x_coord + width, x_coord + width, x_coord]
        ys = [y, y + height, y + height, y, y]
        for i in range(1, len(coords)):
            # Add coordinates for connecting line segment
            xs += [None, coords[i - 1][1], coords[i][0], None]
            ys += [None, y + height / 2, y + height / 2, None]
            # Add coordinates for individual rectangle
            width = abs(coords[i][0] - coords[i][1])
            x_coord = min(coords[i][0], coords[i][1])
            xs += [x_coord, x_coord, x_coord + width, x_coord + width, x_coord]
            ys += [y, y + height, y + height, y, y]
        self.x = numpy.array(xs)
        self.y = numpy.array(ys)


class Diamond(Shape):
    def __init__(self, coords, y=0, height=0.25, **kwargs) -> None:
        super().__init__(**kwargs)
        width = abs(coords[0][0] - coords[0][1])
        x_coord = min(coords[0][0], coords[0][1])
        hh = height / 2
        xs = [
            x_coord,
            x_coord + width / 2,
            x_coord + width,
            x_coord + width / 2,
            x_coord,
        ]
        ys = [y, y + hh, y, y - hh, y]
        for i in range(1, len(coords)):
            # Add coordinates for connecting line segment
            xs += [None, coords[i - 1][1], coords[i][0], None]
            ys += [None, y, y, None]
            # Add coordinates for individual diamond
            width = abs(coords[i][0] - coords[i][1])
            x_coord = min(coords[i][0], coords[i][1])
            xs += [
                x_coord,
                x_coord + width / 2,
                x_coord + width,
                x_coord + width / 2,
                x_coord,
            ]
            ys += [y, y + hh, y, y - hh, y]
        self.x = numpy.array(xs)
        self.y = numpy.array(ys)


def _calc_arrow_width(
    *,
    parent_length: int,
    feature_width: float,
    frac: float = 0.05,
) -> float:
    aw = parent_length * frac
    if aw > feature_width:
        aw = feature_width * frac
    return aw


def _make_rectangles(coords: list[list[int]], y=0, height=0.25) -> tuple[list[float]]:
    xs = []
    ys = []
    y_coord = [y, y + height, y + height, y, y]
    for coord in coords:
        # Add coordinates for individual rectangle
        width = abs(coord[0] - coord[1])
        x_coord = min(coord[0], coord[1])
        xs.append([x_coord, x_coord, x_coord + width, x_coord + width, x_coord])
        ys.append(y_coord)
    return xs, ys


def _make_arrow_head(
    coord: list[float],
    y: float,
    height: float,
    head_width_frac: float,
    reverse: bool,
    parent_length: int,
) -> tuple[list[float]]:
    width = abs(coord[0] - coord[1])
    x_coord = min(coord[0], coord[1])
    hh = height * head_width_frac * 2
    hw = _calc_arrow_width(
        parent_length=parent_length,
        feature_width=width,
        frac=head_width_frac,
    )

    # Coordinates for arrow head
    arrow_x = [
        x_coord,
        x_coord + width - hw,
        x_coord + width - hw,
        x_coord + width,
        x_coord + width - hw,
        x_coord + width - hw,
        x_coord,
        x_coord,
    ]
    arrow_y = [
        y,
        y,
        y - hh,
        y + height / 2,
        y + height + hh,
        y + height,
        y + height,
        y,
    ]
    if reverse:
        arrow_x = numpy.array(arrow_x)
        arrow_x = list(numpy.flip(arrow_x.max() - arrow_x + arrow_x.min()))
        arrow_y = list(numpy.flip(arrow_y))
    return arrow_x, arrow_y


def _connecting_lines(
    x_grp: list[list[float]],
    y_grp: list[list[float]],
    connect_y: float,
) -> tuple[list[float | None]]:
    new_xs, new_ys = x_grp[0][:], y_grp[0][:]
    y_connect_coord = [None, connect_y, connect_y, None]
    for i, (xs, ys) in enumerate(zip(x_grp[1:], y_grp[1:], strict=False), 1):
        new_xs.extend([None, max(x_grp[i - 1]), min(x_grp[i]), None, *xs])
        new_ys.extend(y_connect_coord + ys)
    return new_xs, new_ys


class Arrow(Shape):
    def __init__(
        self,
        coords,
        y=0,
        height=0.25,
        arrow_head_w=0.05,
        reverse=False,
        parent_length: int | None = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)

        if len(coords) == 1 or reverse:
            head_index = 0
            rec_coords = None if len(coords) == 1 else coords[1:]
        else:
            head_index = -1
            rec_coords = coords[:-1]

        arrow_x, arrow_y = _make_arrow_head(
            coords[head_index],
            y=y,
            height=height,
            head_width_frac=arrow_head_w,
            reverse=reverse,
            parent_length=parent_length,
        )
        if rec_coords:
            xs, ys = _make_rectangles(rec_coords, y=y, height=height)
        else:
            xs, ys = [], []

        if reverse:
            xs.insert(0, arrow_x)
            ys.insert(0, arrow_y)
        else:
            xs.append(arrow_x)
            ys.append(arrow_y)

        xs, ys = _connecting_lines(xs, ys, y + height / 2)
        self.x = numpy.array(xs)
        self.y = numpy.array(ys)


# https://plot.ly/python/marker-style/
# https://plot.ly/python/reference/#scatter-marker-symbol
class Point(Shape):
    _mode = "markers"

    def __init__(self, x, y, size=14, symbol="square", **kwargs) -> None:
        super().__init__(**kwargs)
        self.x = numpy.array([x], dtype="O")
        self.y = numpy.array([y], dtype="O")
        self._size = size
        self._symbol = symbol


class _MakeShape:
    """container class that builds annotation shapes"""

    _colors = {
        "cds": "rgba(0,0,150,0.75)",
        "rrna": "rgba(0,0,150,0.75)",
        "misc_rna": "rgba(0,0,150,0.75)",
        "trna": "rgba(0,0,150,0.75)",
        "exon": "rgba(0,0,100,0.75)",
        "gene": "rgba(161,0,0,0.75)",
        "transcript": "rgba(140,102,139,0.75)",
        "mrna": "rgba(140,102,139,0.75)",
        "snp": "rgba(200,0,0,0.75)",
        "snv": "rgba(200,0,0,0.75)",
    }
    _shapes = {
        "cds": Arrow,
        "rrna": Arrow,
        "misc_rna": Arrow,
        "trna": Arrow,
        "exon": Arrow,
        "transcript": Arrow,
        "mrna": Arrow,
        "gene": Arrow,
        "repeat": Rectangle,
        "snp": Point,
        "snv": Point,
        "variation": Diamond,
    }

    def __call__(self, type_=None, name=None, coords=None, **kwargs):
        if hasattr(type_, "map"):
            if not type_.map.useful:
                return None

            name = type_.name if type_.map.complete else f"{type_.name} (incomplete)"
            coords = type_.map.get_coordinates()
            reverse = type_.reversed
            type_ = type_.biotype
        else:
            reverse = coords[0][0] > coords[-1][1]
            if coords is None:
                msg = "No coordinates defined"
                raise ValueError(msg)
        kwargs |= {"reverse": reverse}

        klass = self._shapes.get(type_.lower(), Rectangle)
        color = self.get_colour(type_)
        if klass != Arrow:
            kwargs.pop("reverse", None)

        if klass != Point:
            result = klass(
                name=type_,
                text=name,
                legendgroup=type_,
                coords=coords,
                fillcolor=color,
                **kwargs,
            )
        else:
            result = Point(
                name=type_,
                text=name,
                legendgroup=type_,
                x=min(coords[0][0], coords[-1][1]),
                y=1,
                size=14,
                symbol="square",
                fillcolor=color,
                **kwargs,
            )
        return result

    def get_colour(self, label):
        from plotly.colors import sample_colorscale

        label = label.lower()

        if label not in self._colors:
            self._colors[label] = sample_colorscale("Viridis", numpy.random.rand())[0]

        return self._colors[label]


make_shape = _MakeShape()
