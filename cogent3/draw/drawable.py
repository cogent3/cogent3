from plotly.offline import iplot as _iplot
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
                 showlegend=True, visible_axes=True):
        self._traces = traces or []
        self._layout = dict(title=title,
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
        titles = [tr.title for tr in self.traces]
        return titles

    def pop_trace(self, title):
        """removes the trace with a matching title attribute"""
        index = self.get_trace_titles().index(title)
        return self.traces.pop(index)

    def add_trace(self, trace):
        self.traces.append(trace)

    def bound_to(self, obj):
        """returns obj with self bound to it"""
        return bind_drawable(obj, self)

    @property
    def figure(self):
        if not self.traces:
            # try:
            self._build_fig()
            # except AttributeError as err:
            #     raise RuntimeError(err.args[0])
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
