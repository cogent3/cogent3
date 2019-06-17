import numpy
from plotly.offline import iplot as _iplot
from plotly import graph_objs as go, tools
from plotly.io import write_image, to_image

__author__ = "Rahul Ghangas and Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Rahul Ghangas", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class Drawable:
    """container object for Plotly figures"""

    def __init__(self, title=None, traces=None, width=None, height=None,
                 showlegend=True, visible_axes=True, layout=None,
                 xtitle=None, ytitle=None):
        self._traces = traces or []
        self._layout = go.Layout(title=title,
                                 font=dict(family='Balto', size=14),
                                 width=width,
                                 height=height,
                                 autosize=False,
                                 showlegend=showlegend,
                                 xaxis=dict(visible=visible_axes),
                                 yaxis=dict(visible=visible_axes),
                                 hovermode='closest',
                                 plot_bgcolor='rgb(245,245,245)',
                                 margin=dict(l=50, r=50, t=50)
                                 )
        layout = layout or {}
        self._layout.update(**layout)
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
        self.layout.xaxis.title = self.xtitle
        self.layout.yaxis.title = self.ytitle
        return dict(data=self.traces, layout=self.layout)

    def iplot(self, *args, **kwargs):
        _iplot(self.figure, *args, **kwargs)

    def write(self, path, **kwargs):
        """writes static image file, suffix dictates format"""
        write_image(self.figure, path, **kwargs)

    def to_image(self, format='png', **kwargs):
        return to_image(self.figure, format=format, **kwargs)


def _iplot_(cls, width=None, height=None):
    layout = {}
    if width:
        layout['width'] = width
    if height:
        layout['height'] = height
    if layout:
        cls.drawable.layout.update(layout)
    _iplot(cls.drawable.figure)


def bind_drawable(obj, drawable):
    """binds drawable"""
    from types import MethodType
    obj.drawable = drawable
    obj.iplot = MethodType(_iplot_, obj)
    return obj


class AnnotatedDrawable(Drawable):
    """supports drawing with left and bottom tracks of annotations"""

    def __init__(self, core, left_track=None, bottom_track=None, xtitle=None,
                 ytitle=None, xrange=None, yrange=None, width=500, height=500,
                 layout=None):
        super(AnnotatedDrawable, self).__init__(visible_axes=True,
                                                showlegend=True,
                                                width=width,
                                                height=height,
                                                layout=layout)

        self.xtitle = xtitle
        self.ytitle = ytitle
        self.yrange = yrange
        self.xrange = xrange
        self.core = core
        self.left_track = left_track
        self.bottom_track = bottom_track

    def _build_fig(self, xaxis='x', yaxis='y'):
        f = self.core.figure
        try:
            traces = f.traces
            self.layout.update(f.layout)
        except AttributeError:
            traces = f['data']
            self.layout.update(f['layout'])
        for trace in traces:
            trace.xaxis = xaxis
            trace.yaxis = yaxis
        self._traces = traces

    def _build_2x2_fig(self):
        if not self.traces:
            self._build_fig(xaxis='x2', yaxis='y2')

        space = 0.01
        fig = tools.make_subplots(rows=2, cols=2,
                                  specs=[[{}, {}],
                                         [None, {}]],
                                  horizontal_spacing=space,
                                  vertical_spacing=space,
                                  print_grid=False,
                                  row_width=[0.1, 0.9],
                                  column_width=[0.1, 0.9])

        layout = fig.layout
        layout.update(self.layout)

        # common settings
        ticks_off_kwargs = dict(showticklabels=False, mirror=True,
                                showgrid=False,
                                showline=True,
                                zeroline=False,
                                ticks='')
        ticks_on_kwargs = dict(showticklabels=True, mirror=True, showgrid=False,
                               showline=True,
                               zeroline=False)

        # dotpolot traces and layout
        fig.add_traces(self.traces)

        layout.xaxis2.update(range=self.xrange, **ticks_off_kwargs)
        layout.yaxis2.update(range=self.yrange, **ticks_off_kwargs)

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
            trace.xaxis = 'x3'
            trace.yaxis = 'y3'
            traces.append(trace)
            max_y = max(numpy.max(trace.y), max_y)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        bottom_range = [0, int(max_y) + 1]

        # add all traces
        fig.add_traces(traces)
        # configure easy axes for titles, limits, border and ticks
        layout.yaxis.update(title=self.ytitle,
                            range=self.yrange,
                            **ticks_on_kwargs)

        layout.xaxis3.update(title=self.xtitle,
                             range=self.xrange,
                             **ticks_on_kwargs)

        # adjust row width of left plot for number of feature tracks
        min_range = min(left_range[1], bottom_range[1])
        left_prop = left_range[1] / min_range

        # first the top row
        xaxis_domain = list(layout.xaxis.domain)
        xaxis_domain[1] = left_prop * xaxis_domain[1]
        layout.xaxis.update(title=None,
                            range=left_range,
                            domain=xaxis_domain,
                            **ticks_off_kwargs)
        layout.xaxis.update(title=None,
                            range=left_range,
                            domain=xaxis_domain,
                            **ticks_off_kwargs)

        layout.xaxis2.domain = (xaxis_domain[1] + space, 1.0)

        # now the right column
        bottom_prop = bottom_range[1] / min_range
        yaxis_domain = list(layout.yaxis3.domain)
        yaxis_domain[1] = bottom_prop * yaxis_domain[1]
        layout.yaxis3.update(title=None,
                             range=bottom_range,
                             domain=yaxis_domain,
                             **ticks_off_kwargs)

        # and bottom of the boxes above
        layout.yaxis.domain = (yaxis_domain[1] + space, 1.0)
        layout.yaxis2.domain = (yaxis_domain[1] + space, 1.0)

        return fig

    def _build_2x1_fig(self):
        """2 rows, one column, dotplot and seq1 annotated"""
        if not self.traces:
            self._build_fig()

        fig = tools.make_subplots(rows=2, cols=1,
                                  vertical_spacing=0.015,
                                  print_grid=False,
                                  row_width=[0.1, 0.9],
                                  shared_xaxes=True)
        layout = fig.layout
        layout.update(self.layout)

        # common settings
        ticks_off_kwargs = dict(showticklabels=False,
                                mirror=True,
                                showgrid=False,
                                showline=True)
        ticks_on_kwargs = dict(showticklabels=True,
                               mirror=True,
                               showgrid=False,
                               showline=True)

        # core traces and layout
        fig.add_traces(self.traces)

        layout.xaxis.update(title=self.xtitle,
                            range=self.xrange,
                            **ticks_on_kwargs)
        layout.yaxis.update(title=self.ytitle,
                            range=self.yrange,
                            **ticks_on_kwargs)

        # bottom traces
        seen_types = set()
        max_y = 0
        traces = []
        for trace in self.bottom_track.traces:
            trace.yaxis = 'y2'
            traces.append(trace)
            max_y = max(numpy.max(trace.y), max_y)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        fig.add_traces(traces)
        layout.yaxis2.update(title=None,
                             range=[0, int(max_y) + 1],
                             **ticks_off_kwargs)
        return fig

    def _build_1x2_fig(self):
        if not self.traces:
            self._build_fig(xaxis='x2')
        fig = tools.make_subplots(rows=1, cols=2,
                                  horizontal_spacing=0.01,
                                  print_grid=False,
                                  column_width=[0.1, 0.9],
                                  shared_yaxes=True)
        layout = fig.layout
        layout.update(self.layout)

        # common settings
        ticks_off_kwargs = dict(showticklabels=False,
                                mirror=True,
                                showgrid=False,
                                showline=True)
        ticks_on_kwargs = dict(showticklabels=True,
                               mirror=True,
                               showgrid=False,
                               showline=True)

        # core traces and layout
        fig.add_traces(self.traces)

        layout.xaxis2.update(title=self.xtitle,
                             range=self.xrange,
                             **ticks_on_kwargs)
        layout.yaxis.update(title=self.ytitle,
                            range=self.yrange,
                            **ticks_on_kwargs)

        # left track
        seen_types = set()
        max_x = 0
        traces = []
        for trace in self.left_track.traces:
            traces.append(trace)
            max_x = max(numpy.max(trace.x), max_x)
            if trace.legendgroup in seen_types:
                trace.showlegend = False
            seen_types.add(trace.legendgroup)

        fig.add_traces(traces)
        layout.xaxis.update(title=None,
                            range=[0, int(max_x) + 1],
                            **ticks_off_kwargs)
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


class Shape:
    _mode = 'lines'

    def __init__(self, name=None, text=None, filled=True, legendgroup=None,
                 showlegend=True, hoverinfo=None,
                 fillcolor=None):
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
        data = go.Scatter(x=self.x, y=self.y,
                          mode=self._mode,
                          fill='toself',
                          fillcolor=self.fillcolor,
                          line=dict(color=self.fillcolor),
                          text=self.text,
                          name=name,
                          legendgroup=self._legendgroup,
                          showlegend=self._showlegend,
                          hoverinfo='text')
        return data


class Rectangle(Shape):
    def __init__(self, x, width, y=0, height=0.5, **kwargs):
        super(Rectangle, self).__init__(**kwargs)
        self.x = numpy.array([x, x, x + width, x + width, x])
        self.y = numpy.array([y, y + height, y + height, y, y])


class Diamond(Shape):
    def __init__(self, x, width, y=0, height=0.5, **kwargs):
        super(Diamond, self).__init__(**kwargs)
        hw = width / 2
        hh = height / 2
        self.x = numpy.array([x - hw, x, x + hw, x, x - hw])
        self.y = numpy.array([y, y + hh, y, y - hh, y])


class Arrow(Shape):
    def __init__(self, x, width, y=0, height=0.5, arrow_head_w=0.1,
                 reverse=False, **kwargs):
        super(Arrow, self).__init__(**kwargs)
        hw = width * arrow_head_w * 2
        hh = height * arrow_head_w * 2
        self.x = numpy.array([x, x + width - hw, x + width - hw, x + width,
                              x + width - hw, x + width - hw, x, x])
        self.y = numpy.array([y, y, y - hh, y + height / 2,
                              y + height + hh, y + height, y + height, y])
        if reverse:
            self.x = numpy.flip(self.x.max() - self.x + self.x.min())
            self.y = numpy.flip(self.y)


# https://plot.ly/python/marker-style/
# https://plot.ly/python/reference/#scatter-marker-symbol
class Point(Shape):
    _mode = 'markers'

    def __init__(self, x, y, size=14, symbol='square', **kwargs):
        super(Point, self).__init__(**kwargs)
        self.x = numpy.array([x], dtype='O')
        self.y = numpy.array([y], dtype='O')
        self._size = size
        self._symbol = symbol


class _MakeShape:
    """container class that builds annotation shapes"""
    _colors = dict(cds='rgba(0,0,150,0.5)',
                   exon='rgba(0,0,100,0.5)',
                   gene='rgba(0,0,150,0.5)',
                   transcript='rgba(0,0,200,0.5)',
                   snp='rgba(200,0,0,0.5)',
                   snv='rgba(200,0,0,0.5)')
    _shapes = dict(cds=Arrow,
                   exon=Arrow,
                   transcript=Arrow,
                   gene=Arrow,
                   repeat=Rectangle,
                   snp=Point,
                   snv=Point)

    def __call__(self, type_=None, name=None, coords=None, width=None,
                 **kwargs):
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
            kwargs.pop('reverse', None)

        if klass != Point:
            result = klass(name=type_, text=name, legendgroup=type_,
                           x=start, width=width,
                           fillcolor=color, **kwargs)
        else:
            result = Point(name=type_, text=name, legendgroup=type_,
                           x=start, y=1, size=14, symbol='square',
                           fillcolor=color, **kwargs)
        return result


make_shape = _MakeShape()
