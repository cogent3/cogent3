import numpy
from plotly.offline import iplot as _iplot


class Drawable:
    # Superclass for objects which can generate a figure, in order
    # to supply consistent and convenient iplot() and write()
    # methods.
    # Subclasses must provide .make_figure() which will make use of
    # _makeFigure() matplotlib.pyplot import done at runtime to give the
    # user every chance to change the matplotlib backend first

    def _set_initial_layout(self, width=None, height=None, visible_axes=False,
                    showlegend=False, **kw):
        data = dict()
        layout = dict(title='',
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

        if 'layout' in kw.keys():
            layout.update(kw['layout'])

        self._trace = dict(data=[data], layout=layout)

    def iplot(self, *args, **kwargs):
        if not self._trace:
            self._set_initial_layout(*args, **kwargs)
        _iplot(self._trace)
