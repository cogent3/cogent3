import numpy

from cogent3.util.union_dict import UnionDict


__author__ = "Rahul Ghangas and Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Rahul Ghangas", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.08.06a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


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
    ):
        self._traces = traces or []
        title = title if title is None else dict(text=title)
        self._layout = UnionDict(
            title=title,
            font=dict(family="Balto", size=14),
            width=width,
            height=height,
            autosize=False,
            showlegend=showlegend,
            xaxis=dict(visible=visible_axes),
            yaxis=dict(visible=visible_axes),
            hovermode="closest",
            plot_bgcolor=None,
            margin=dict(l=50, r=50, t=50),
        )
        layout = layout or {}
        self._layout |= layout
        self.xtitle = xtitle
        self.ytitle = ytitle

    def _repr_html_(self):
        self.iplot()

    @property
    def layout(self):
        return self._layout

    @layout.setter
    def layout(self, value):
        self.layout.update(value)

    @property
    def traces(self):
        return self._traces

    def get_trace_titles(self):
        titles = [tr.name for tr in self.traces]
        return titles

    def pop_trace(self, title):
        """removes the trace with a matching title attribute"""
        try:
            index = self.get_trace_titles().index(title)
        except ValueError:
            UserWarning(f"no trace with name {title}")
            return

        return self.traces.pop(index)

    def remove_traces(self, names):
        """removes traces by name

        Parameters
        ----------
        names : str or iterable of str
            trace names

        """
        if not self.traces:
            self._build_fig()

        names = names if type(names) != str else [names]
        for name in names:
            _ = self.pop_trace(name)

    def add_trace(self, trace):
        self.traces.append(trace)

    def bound_to(self, obj):
        """returns obj with self bound to it"""
        return bind_drawable(obj, self)

    @property
    def figure(self):
        if not self.traces:
            self._build_fig()
        xtitle = self.xtitle if not self.xtitle else dict(text=self.xtitle)
        ytitle = self.ytitle if not self.ytitle else dict(text=self.ytitle)
        self.layout.xaxis.title = xtitle
        self.layout.yaxis.title = ytitle
        return UnionDict(data=self.traces, layout=self.layout)

    def iplot(self, *args, **kwargs):
        from plotly.offline import iplot as _iplot

        _iplot(self.figure, *args, **kwargs)

    def write(self, path, **kwargs):
        """writes static image file, suffix dictates format"""
        from plotly.io import write_image

        write_image(self.figure, path, **kwargs)

    def to_image(self, format="png", **kwargs):
        """creates static image, suffix dictates format"""
        from plotly.io import to_image

        return to_image(self.figure, format=format, **kwargs)


def _iplot_(cls, width=None, height=None):
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
    return obj


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
    ):
        super(AnnotatedDrawable, self).__init__(
            visible_axes=True,
            showlegend=True,
            width=width,
            height=height,
            layout=layout,
        )

        self.xtitle = xtitle
        self.ytitle = ytitle
        self.yrange = yrange
        self.xrange = xrange

        core.title = title or core.title
        self.core = core
        self.left_track = left_track
        self.bottom_track = bottom_track

    def _build_fig(self, xaxis="x", yaxis="y"):
        f = self.core.figure
        try:
            traces = f.traces
            self.layout |= dict(f.layout)
        except AttributeError:
            traces = f["data"]
            self.layout |= f["layout"]
        for trace in traces:
            trace.xaxis = xaxis
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
            }
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
            max_x = max(numpy.max(trace.x), max_x)
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
            max_y = max(numpy.max(trace.y), max_y)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        bottom_range = [0, int(max_y) + 1]

        # add all traces
        fig.data.extend(traces)
        # configure axes for titles, limits, border and ticks
        fig.layout.yaxis |= dict(
            title=dict(text=self.ytitle), range=self.yrange, **ticks_on_kwargs
        )

        fig.layout.xaxis3 |= dict(
            title=dict(text=self.xtitle), range=self.xrange, **ticks_on_kwargs
        )

        # adjust row width of left plot for number of feature tracks
        min_range = min(left_range[1], bottom_range[1])
        left_prop = left_range[1] / min_range

        # first the top row
        xaxis_domain = list(layout.xaxis.domain)
        xaxis_domain[1] = left_prop * xaxis_domain[1]
        fig.layout.xaxis |= dict(
            title=None, range=left_range, domain=xaxis_domain, **ticks_off_kwargs
        )
        fig.layout.xaxis |= dict(
            title={}, range=left_range, domain=xaxis_domain, **ticks_off_kwargs
        )

        space = 0.01
        fig.layout.xaxis2.domain = (xaxis_domain[1] + space, 1.0)
        fig.layout.xaxis3.domain = (xaxis_domain[1] + space, 1.0)

        # now the right column
        bottom_prop = bottom_range[1] / min_range
        yaxis_domain = list(layout.yaxis3.domain)
        yaxis_domain[1] = bottom_prop * yaxis_domain[1]
        fig.layout.yaxis3 |= dict(
            title={}, range=bottom_range, domain=yaxis_domain, **ticks_off_kwargs
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
        layout |= dict(self.layout)
        fig = UnionDict(data=[], layout=layout)

        # common settings
        ticks_off_kwargs = dict(_ticks_off)
        ticks_on_kwargs = dict(_ticks_on)

        # core traces and layout
        fig.data.extend(self.traces)

        fig.layout.xaxis |= dict(
            title=dict(text=self.xtitle), range=self.xrange, **ticks_on_kwargs
        )
        fig.layout.yaxis |= dict(
            title=dict(text=self.ytitle), range=self.yrange, **ticks_on_kwargs
        )

        # bottom traces
        seen_types = set()
        max_y = 0
        traces = []
        for trace in self.bottom_track.traces:
            trace.yaxis = "y2"
            trace.xaxis = "x"
            traces.append(trace)
            max_y = max(numpy.max(trace.y), max_y)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        fig.data.extend(traces)
        fig.layout.yaxis2 |= dict(
            title={}, range=[0, int(max_y) + 1], **ticks_off_kwargs
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
            title=self.xtitle, range=self.xrange, **ticks_on_kwargs
        )
        fig.layout.yaxis |= dict(
            title=self.ytitle, range=self.yrange, **ticks_on_kwargs
        )

        # left track
        seen_types = set()
        max_x = 0
        traces = []
        for trace in self.left_track.traces:
            trace.yaxis = "y"
            traces.append(trace)
            max_x = max(numpy.max(trace.x), max_x)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        fig.data.extend(traces)
        fig.layout.xaxis |= dict(
            title=None, range=[0, int(max_x) + 1], **ticks_off_kwargs
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

        result = func()

        return result

    def remove_track(self, left_track=False, bottom_track=False):
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
    ):
        self.filled = filled
        self.fillcolor = fillcolor
        self._legendgroup = legendgroup
        self._showlegend = showlegend
        self.name = name
        self.text = text
        self._hoverinfo = hoverinfo or name

    def shift(self, x=0, y=0):
        self.x += x
        self.y += y
        return self

    @property
    def height(self):
        return self.top - self.bottom

    @property
    def top(self):
        return numpy.max(self.y)

    @property
    def bottom(self):
        return numpy.min(self.y)

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
        data = UnionDict(
            type="scatter",
            x=self.x,
            y=self.y,
            mode=self._mode,
            fill="toself",
            fillcolor=self.fillcolor,
            line=dict(color=self.fillcolor),
            text=self.text,
            name=name,
            legendgroup=self._legendgroup,
            showlegend=self._showlegend,
            hoverinfo="text",
        )
        return data


class Rectangle(Shape):
    def __init__(self, x, width, y=0, height=0.25, **kwargs):
        super(Rectangle, self).__init__(**kwargs)
        self.x = numpy.array([x, x, x + width, x + width, x])
        self.y = numpy.array([y, y + height, y + height, y, y])


class Diamond(Shape):
    def __init__(self, x, width, y=0, height=0.25, **kwargs):
        super(Diamond, self).__init__(**kwargs)
        hw = width / 2
        hh = height / 2
        self.x = numpy.array([x - hw, x, x + hw, x, x - hw])
        self.y = numpy.array([y, y + hh, y, y - hh, y])


class Arrow(Shape):
    def __init__(
        self, x, width, y=0, height=0.25, arrow_head_w=0.1, reverse=False, **kwargs
    ):
        super(Arrow, self).__init__(**kwargs)
        hw = width * arrow_head_w * 2
        hh = height * arrow_head_w * 2
        self.x = numpy.array(
            [
                x,
                x + width - hw,
                x + width - hw,
                x + width,
                x + width - hw,
                x + width - hw,
                x,
                x,
            ]
        )
        self.y = numpy.array(
            [y, y, y - hh, y + height / 2, y + height + hh, y + height, y + height, y]
        )
        if reverse:
            self.x = numpy.flip(self.x.max() - self.x + self.x.min())
            self.y = numpy.flip(self.y)


# https://plot.ly/python/marker-style/
# https://plot.ly/python/reference/#scatter-marker-symbol
class Point(Shape):
    _mode = "markers"

    def __init__(self, x, y, size=14, symbol="square", **kwargs):
        super(Point, self).__init__(**kwargs)
        self.x = numpy.array([x], dtype="O")
        self.y = numpy.array([y], dtype="O")
        self._size = size
        self._symbol = symbol


class _MakeShape:
    """container class that builds annotation shapes"""

    _colors = dict(
        cds="rgba(0,0,150,0.5)",
        exon="rgba(0,0,100,0.5)",
        gene="rgba(0,0,150,0.5)",
        transcript="rgba(0,0,200,0.5)",
        snp="rgba(200,0,0,0.5)",
        snv="rgba(200,0,0,0.5)",
    )
    _shapes = dict(
        cds=Arrow,
        exon=Arrow,
        transcript=Arrow,
        gene=Arrow,
        repeat=Rectangle,
        snp=Point,
        snv=Point,
    )

    def __call__(self, type_=None, name=None, coords=None, width=None, **kwargs):
        from cogent3.core.annotation import _Annotatable

        if isinstance(type_, _Annotatable):
            name = type_.name
            width = len(type_)
            map = type_.map.get_covering_span()
            reverse = map.reverse
            start = min(map.spans[0].start, map.spans[0].end)
            type_ = type_.type
            kwargs.update(dict(reverse=reverse))

        klass = self._shapes.get(type_.lower(), Rectangle)
        color = self._colors.get(type_.lower(), None)
        if klass != Arrow:
            kwargs.pop("reverse", None)

        if klass != Point:
            result = klass(
                name=type_,
                text=name,
                legendgroup=type_,
                x=start,
                width=width,
                fillcolor=color,
                **kwargs,
            )
        else:
            result = Point(
                name=type_,
                text=name,
                legendgroup=type_,
                x=start,
                y=1,
                size=14,
                symbol="square",
                fillcolor=color,
                **kwargs,
            )
        return result


make_shape = _MakeShape()
