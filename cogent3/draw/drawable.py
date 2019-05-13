from plotly.offline import iplot as _iplot

__author__ = "Rahul Ghangasand Gavin Huttley"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Rahul Ghangas", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class Drawable:
    # Superclass for objects which can generate a figure, in order
    # to supply consistent and convenient get_figure() and write()
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

    def _repr_html_(self):
        self.iplot()

    def iplot(self, *args, **kwargs):
        if not self._trace:
            self._set_initial_layout(*args, **kwargs)
        _iplot(self._trace)
