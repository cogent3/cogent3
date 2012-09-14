"""Biplot and triplot for ordination results.
"""
from itertools import imap
import pylab; from pylab import xlim, ylim, plot, scatter
from matplotlib import cm
from matplotlib.colors import rgb2hex, Normalize, Colormap
from numpy import asarray, isscalar, concatenate, any
from pdb import set_trace

__author__ = "Zongzhi Liu"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Zongzhi Liu"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Zongzhi Liu"
__email__ = "zongzhi.liu@gmail.com"
__status__ = "Development"

def plot_ordination(res, choices=[1,2],
        axis_names='PC', constrained_names=None,
        samples_kw={}, species_kw={}, centroids_kw={}, biplot_kw={},
        axline_kw={}):
    """plot the ordination result dict.

    - res: ordination result as a dict with
        'eigvals': increasing eigen values as a 1darray.
        'samples': sample points as a nx2 array.
        'species': optional species points as a nx2.
        'centroids': optional centroid points as a nx2.
        'biplot': optional biplot arrowheads as a nx2.
    - choices=[1,2]: a 2 item list, including axes(base from 1) to choose
    - axis_names='PC': a name prefix or a list of names for each axis.
    - constrained_names=None: The prefix or a list of constrained axis names.
    - samples_kw, centroids_kw: control the plot of sample points and
      centroid points respectively.  pass to scatter_points(label=, c=, s=, ...)
    - species_kw: control the plot of species points.  pass to
      plot_points(points, cml=, label=, ...)
        - cml: a str for color+marker+linestyle. eg. 'r+:'
    - biplot_kw: pass to arrows(points_from, points_to, label=, ...)
    - axline_kw: pass to axhline(...) and axvline(...).
    """
    choices = asarray(choices) -1
    samples, species, centroids, biplot, evals = [asarray(res.get(k, None))
            for k in ['samples', 'species', 'centroids', 'biplot', 'eigvals']]
    if isinstance(axis_names, str):
        axis_names = ['%s%i' % (axis_names, i+1) for i in range(len(evals))]
    # draw the axis lines
    axline_kw = dict({'color':'gray'}, **axline_kw)
    pylab.axvline(**axline_kw)
    pylab.axhline(**axline_kw)
    # calc percentages from evals and label them
    evals = asarray(evals)
    evals[evals<0] = 0
    percs = 100 * (evals / evals.sum())
    pylab.xlabel('%s - %.2f%%' % (axis_names[choices[0]], percs[choices[0]]))
    pylab.ylabel('%s - %.2f%%' % (axis_names[choices[1]], percs[choices[1]]))
    #plot the speceies points in red +
    if any(species):
        species_kw = dict({'cml': 'r+', 'label_kw':{'size':'smaller'}},
                **species_kw) #set default
        plot_points(species[:, choices], **species_kw)
    # scatter the sample points in black
    default = {'c': 'k', 's': 50, 'alpha': 0.5, 'label_kw': {'size': 'medium'}}
    scatter_points(samples[:, choices], **dict(default, **samples_kw))
    # scatter the centroids
    if any(centroids):
        default = {'c':'b', 's':0, 'label': 'X', 'label_kw':{'size':'larger',
            'color':'b', 'ha':'center', 'va':'center'}}
        scatter_points(centroids[:, choices], **dict(default, **centroids_kw))
    # arrow the biplot points
    if any(biplot):
        default = {'c':'b', 'label_kw':{'size':'larger', 'color':'b'}}
        arrows([[0,0]] * len(biplot), biplot[:, choices],
                **dict(default, **biplot_kw))
    # calc the constrained percentage and title it
    if constrained_names:
        if isinstance(constrained_names, str): #a prefix
            constrained_names = [n for n in axis_names if
                    n.startswith(constrained_names)]
            con_idxs = [i for i, n in enumerate(axis_names)
                    if n in constrained_names]
        con_perc = percs[con_idxs].sum()
        pylab.title('%.2f%% constrained' % con_perc)


####
# support functions
def map_colors(colors, cmap=None, lut=None, mode='hexs', **norm_kw):
    """return a list of rgb tuples/hexs from color numbers.

    - colors: a seq of color numbers.
    - cmap: a Colormap or a name like 'jet' (passto cm.get_cmap(cmap, lut)
    - mode: one of  ['hexs', 'tuples', 'arrays']

    Ref: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
    """
    modes = ['hexs', 'tuples', 'arrays']
    if mode not in modes:
        raise ValueError('mode must be one of %s, but got %s'
                % (modes, mode))
    if not isinstance(cmap, Colormap):
        cmap = cm.get_cmap(cmap, lut=lut)
    rgba_arrays = cmap(Normalize(**norm_kw)(colors))
    rgb_arrays = rgba_arrays[:, :-1] #without alpha
    if mode == 'arrays':
        return rgb_arrays
    elif mode == 'tuples':
        return list(imap(tuple, rgb_arrays))
    else: # mode == 'hexs':
        return list(imap(rgb2hex, rgb_arrays))

def text_points(points, strs, **kw):
    """print each str at each point.

    - points: a seq of xy pairs.
    - strs: a str or a seq of strings with the length of points.
    **kw: params pass to pylab.text(x, y, s, **kw)
        - horizontalalignment or ha: ['center', 'left', 'right']
        - verticalalignment or va: ['center', 'top', 'bottom']
        - size, rotation, ...
    """
    xs, ys = asarray(points, float).T
    if isinstance(strs, str): #vectorize strs
        strs = [strs] * len(xs)
    for x, y, s in zip(xs, ys, strs):
        if not s: continue
        pylab.text(x, y, s, **kw)

def plot_points(points, cml=None, label=None, label_kw={}, **kw):
    """plot at each point (x,y).

    - points: a seq of xy pairs.
    - cml: a str combining color, marker and linestyle, default to be '-'
    - label, label_kw: pass to pylab,text(x,y, lbl, **kw)
    **kw: pass to plot(x,y,**kw)
        -color, marker(+,), linestyle(-, --, :,)
        Note: label was overwritten to be text label for each point.
    """
    points = asarray(points, float)
    xs, ys = points.T
    if cml:
        pylab.plot(xs, ys, cml, **kw)
    else:
        pylab.plot(xs, ys, **kw)
    if label:#add label to each point
        text_points(points, label, **label_kw)

def scatter_points(points, s=10, c='b', marker='o', label=None,
        label_kw={} , **kw):
    """scatter enhanced with points, markers, and labels.

    - points: a seq of xy-pairs
    - s=10: size or sizes of bubbles.
    - c='b': a color (str or tuple) or a list of colors
    - marker: a marker (str or tuple) or a list of markers.
    - label: a label (str) or a list of labels.
    - label_kw: will be passed to text() to print the labels.
    **kw: params pass to pylab.scatter()
    """
    points = asarray(points, float)
    xs, ys = points.T
    num_points = len(xs)
    #vectorize size, color, and marker
    if isscalar(s):
        s = [s] * num_points
    if isinstance(c, (str, tuple)):
        c = [c] * num_points
    if isinstance(marker, (str, tuple)):
        marker = [marker] * num_points
    #numer list to hex list
    if isinstance(c[0], (int, float)):
        c = map_colors(c)
    #scatter each point; colormap() will not work properly now.
    for xi, yi, si, ci, mi in zip(xs, ys, s, c, marker):
        scatter([xi], [yi], si, ci, mi, **kw)
    if label: #add label to each point
        text_points(points, label, **label_kw)

def arrows(points_from, points_to, color='k',
        width=None, width_ratio=0.002,
        label=None, label_kw={}, margin_ratio=None, **kw):
    """draw arrows from each point_from to each point_to.

    - points_from, points_to: each is a seq of xy pairs.
    - color or c: arrow color(s).
    - width: arrow width(s), default to auto adjust with width_ratio.
    - width_ratio: set width to minimal_axisrange * width_ratio.
    - label=None: a list of strs to label the arrow heads.
    - label_kw: pass to pylab.text(x,y,txt,...)
    - margin_ratio: if provided as a number, canvas margins will be set to
        axisrange * margin_ratio.
    **kw: params pass to pylab.arrow(,...)
    """
    # convert and valid check of the points
    points_from = asarray(points_from, float)
    points_to = asarray(points_to, float)
    num_points = len(points_from)
    if len(points_to) != num_points:
        raise ValueError('Length of two group of points unmatched')
    if not width: #auto adjust arrow widt
        min_range = min(asarray(xlim()).ptp(), asarray(ylim()).ptp())
        width = min_range * width_ratio
    #vectorize colors and width if necessary
    color = kw.pop('c', color)
    if isscalar(color): color = [color] * num_points
    if isscalar(width): width = [width] * num_points
    #draw each arrow
    #?? length_include_head=True will break??
    default = {'length_includes_head':False, 'alpha':0.75}
    kw = dict(default, **kw)
    for (x0, y0), (x1, y1), w, c in zip(points_from, points_to, width, color):
        pylab.arrow(x0, y0, x1-x0, y1-y0,
            edgecolor=c, facecolor=c, width=w, **kw)
    if label: #label the arrow heads
        text_points(points_to, label, **label_kw)
    #hack fix of axis limits, otherwise some of the arrows will be invisible
    #not neccessary if the points were also scattered
    if margin_ratio is not None:
        x, y = concatenate((points_from, points_to)).T #all the points
        xmargin = x.ptp() * margin_ratio
        ymargin = y.ptp() * margin_ratio
        xlim(x.min()-xmargin, x.max()+xmargin)
        ylim(y.min()-ymargin, y.max()+xmargin)
    return
