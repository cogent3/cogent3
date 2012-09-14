#/usr/bin/env python
"""Provides different kinds of generally useful plots using matplotlib.

Some of these plots are enhancements of the matplotlib versions (e.g. 
hist() or copies of plot types that have been withdrawn from matplotlib
(e.g. scatter_classic).

Notable capabilities include automated series coloring and drawing of 
regression lines, the ability to plot scatterplots with correlated histograms,
etc.

See individual docstrings for more info.
"""
from __future__ import division
from matplotlib import use, rc, rcParams

__author__ = "Stephanie Wilson"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Stephanie Wilson"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

#use('Agg')  #suppress graphical rendering
#rc('text', usetex=True)
rc('font', family='serif')  #required to match latex text and equations
try:
    import Image
    import ImageFilter
except ImportError:
    Image = ImageFilter = None  #only used to smooth contours: skip if no PIL
from numpy import array, shape, fromstring, sqrt, zeros, pi
from cogent.core.usage import UnsafeCodonUsage as CodonUsage
from cogent.maths.stats.test import regress, correlation
from pylab import plot, cm, savefig, gca, gcf, arange, text, subplot, \
    asarray, iterable, searchsorted, sort, diff, concatenate, silent_list, \
    is_string_like, Circle, mean, std, normpdf, legend, contourf, \
    colorbar, ravel, imshow, contour
from matplotlib.font_manager import FontProperties
from os.path import split
#module-level constants
standard_series_colors=['k','r','g','b', 'm','c']

def hist(x, bins=10, normed='height', bottom=0, \
    align='edge', orientation='vertical', width=None, axes=None, **kwargs):
    """Just like the matplotlib hist, but normalizes bar heights to 1.
    
    axes uses gca() by default (built-in hist is a method of Axes).
    
    Original docs from matplotlib:
 
    HIST(x, bins=10, normed=0, bottom=0, orientiation='vertical', **kwargs)

    Compute the histogram of x.  bins is either an integer number of
    bins or a sequence giving the bins.  x are the data to be binned.

    The return values is (n, bins, patches)

    If normed is true, the first element of the return tuple will
    be the counts normalized to form a probability density, ie,
    n/(len(x)*dbin)


    orientation = 'horizontal' | 'vertical'.  If horizontal, barh
    will be used and the "bottom" kwarg will be the left.

    width: the width of the bars.  If None, automatically compute
    the width.

    kwargs are used to update the properties of the
    hist bars
    """
    if axes is None:
        axes = gca()
    if not axes._hold: axes.cla()
    n, bins = norm_hist_bins(x, bins, normed)
    if width is None: width = 0.9*(bins[1]-bins[0])
    if orientation=='horizontal':
        patches = axes.barh(bins, n, height=width, left=bottom, \
            align=align)
    else:
        patches = axes.bar(bins, n, width=width, bottom=bottom, \
            align=align)
    for p in patches:
        p.update(kwargs)
    return n, bins, silent_list('Patch', patches)

def norm_hist_bins(y, bins=10, normed='height'):
    """Just like the matplotlib mlab.hist, but can normalize by height.

    normed can be 'area' (produces matplotlib behavior, area is 1), 
    any False value (no normalization), or any True value (normalization).

    Original docs from matplotlib:

    Return the histogram of y with bins equally sized bins.  If bins
    is an array, use the bins.  Return value is
    (n,x) where n is the count for each bin in x

    If normed is False, return the counts in the first element of the
    return tuple.  If normed is True, return the probability density
    n/(len(y)*dbin)
    
    If y has rank>1, it will be raveled
    Credits: the Numeric 22 documentation
    """
    y = asarray(y)
    if len(y.shape)>1: y = ravel(y)

    if not iterable(bins):
        ymin, ymax = min(y), max(y)
        if ymin==ymax:
            ymin -= 0.5
            ymax += 0.5

        if bins==1: bins=ymax
        dy = (ymax-ymin)/bins
        bins = ymin + dy*arange(bins)
    n = searchsorted(sort(y), bins)
    n = diff(concatenate([n, [len(y)]]))
    if normed:
        if normed == 'area':
            db = bins[1]-bins[0]
        else:
            db = 1.0
        return 1/(len(y)*db)*n, bins
    else:
        return n, bins

def scatter_classic(x, y, s=None, c='b'):
    """
    SCATTER_CLASSIC(x, y, s=None, c='b')

    Make a scatter plot of x versus y.  s is a size (in data coords) and
    can be either a scalar or an array of the same length as x or y.  c is
    a color and can be a single color format string or an length(x) array
    of intensities which will be mapped by the colormap jet.

    If size is None a default size will be used

    Copied from older version of matplotlib -- removed in version 0.9.1
    for whatever reason.
    """
    self = gca()
    if not self._hold: self.cla()
    if is_string_like(c):
        c = [c]*len(x)
    elif not iterable(c):
        c = [c]*len(x)
    else:
        norm = normalize()
        norm(c)
        c = cm.jet(c)

    if s is None:
        s = [abs(0.015*(amax(y)-amin(y)))]*len(x)
    elif not iterable(s):
        s = [s]*len(x)

    if len(c)!=len(x):
        raise ValueError, 'c and x are not equal lengths'
    if len(s)!=len(x):
        raise ValueError, 's and x are not equal lengths'

    patches = []
    for thisX, thisY, thisS, thisC in zip(x,y,s,c):
        circ = Circle( (thisX, thisY),
                       radius=thisS,
                       )
        circ.set_facecolor(thisC)
        self.add_patch(circ)
        patches.append(circ)
    self.autoscale_view()
    return patches

def as_species(name, leave_path=False):
    """Cleans up a filename into a species name, italicizing it in latex."""
    #trim extension if present
    dot_location = name.rfind('.')
    if dot_location > -1:
        name = name[:dot_location]
    #get rid of _small if present -- used for debugging
    if name.endswith('_small'):
        name = name[:-len('_small')]
    if name.endswith('_codon_usage'):
        name = name[:-len('_codon_usage')]
    #get rid of path unless told to leave it
    name = split(name)[-1]
    #replace underscores with spaces
    name = name.replace('_', ' ')
    #make sure the first letter of the genus is caps, and not the first letter
    #of the species
    fields = name.split()
    fields[0] = fields[0].title()
    #assume second field is species name
    if len(fields) > 1:
        fields[1] = fields[1].lower()
    binomial = ' '.join(fields)
    if rcParams.get('text.usetex'):
        binomial = r'\emph{' + binomial + '}'
    return binomial

def frac_to_psq(frac, graph_size):
    """Converts diameter as fraction of graph to points squared for scatter.
    
    frac: fraction of graph (e.g. .01 is 1% of graph size)
    graph_size: graph size in inches
    """
    points = frac * graph_size * 72
    return pi * (points/2.0)**2
    
def init_graph_display(title=None, aux_title=None, size=4.0, \
    graph_shape='sqr', graph_grid=None, x_label='', y_label='', \
    dark=False, with_parens=True, prob_axes=True, axes=None, num_genes=None):
    """Initializes a range of graph settings for standard plots.

    These settings include:
        - font sizes based on the size of the graph
        - graph shape
        - grid, including lines for x=y or at x and y = 0.5
        - title, auxillary title, and x and y axis labels
                
    Parameters:
        title: displayed on left of graph, at the top, latex-format string
        
        aux_title: displayed on top right of graph, latex-format string.
        typically used for number of genes.

        size:   size of graph, in inches

        graph_shape: 'sqr' for square graphs, 'rect' for graphs that include
        a colorbar, 3to1: width 3 to height 1.

        graph_grid: background grid for the graph. Currently recognized grids
        are '/' (line at x=y) and 't' (cross at x=.5 and y=.5).

        x_label: label for x axis, latex-format string.

        y_label: label for y axis, latex-format string.

        dark: set to True if dark background, reverses text and tick colors.
        
        with_parens: if True (default), puts parens around auxillary title
        
    returns font, label_font_size (for use in producing additional labels in 
    calling function).
    """
    if dark:
        color='w'
    else:
        color='k'

    rect_scale_factor = 1.28    #need to allow for legend while keeping graph
                                #square; empirically determined at 1.28
    font_size = int(size*3-1)   #want 11pt font w/ default graph size 4" sqr
    label_scale_factor = 0.8
    label_font_size = font_size * label_scale_factor
    label_offset = label_font_size * 0.5
    axis_label_font={'fontsize':font_size}
    font={'fontsize':font_size, 'color':color}
    

    if graph_shape == 'sqr':
        gcf().set_size_inches(size,size)
    elif graph_shape == 'rect':
        #scaling for sqr graphs with colorbar
        gcf().set_size_inches(size*rect_scale_factor,size)
    elif graph_shape == '3to1':
        gcf().set_size_inches(3*size, size)
    elif graph_shape == '2to1':
        gcf().set_size_inches(2*size, size)
    else:
        raise ValueError, "Got unknown graph shape %s" % graph_shape
    
    #set or create axes
    if axes is None:
        axes = gca()

    min_x, max_x =axes.get_xlim()
    min_y, max_y = axes.get_ylim()
    x_range = abs(max_x - min_x)
    y_range = abs(max_y - min_y)

    
    min_offset = (x_range * 0.05) + min_x #minimum offset, e.g. for text
    max_offset = max_y - (y_range * 0.05) 

    #draw grid manually: these are in data coordinates. 
    if graph_grid == 't':
        #grid lines at 0.5 on each axis, horiz & vertic
        axes.axvline(x=.5, ymin=0, ymax=1, color=color, linestyle=':')
        axes.axhline(y=.5, xmin=0, xmax=1, color=color, linestyle=':')
    elif graph_grid == '/':
        #diagonal gridlines from 0,0 to 1,1.
        axes.plot([0,1], color=color, linestyle=':')
    else:
        pass    #ignore other choices
        
    #remove default grid
    axes.grid(False)

    #set x and y labels
    axes.set_ylabel(y_label, axis_label_font)
    axes.set_xlabel(x_label, axis_label_font)

    #add title/aux_title to graph directly. Note that we want 
    #the tops of these to be fixed, and we want the label to be 
    #left-justified and the number of genes to be right justified, 
    #so that it still works when we resize the graph.
    if title is not None:
        axes.text(min_offset, max_offset, str(title), font, \
            verticalalignment='top', horizontalalignment='left')
    #use num_genes as aux_title by default
    aux_title = num_genes or aux_title
    if aux_title is not None:
        if with_parens:
            aux_title='('+str(aux_title)+')'
        axes.text(max_offset, max_offset, str(aux_title), font,
             verticalalignment='top', horizontalalignment='right')
    if prob_axes:
        init_ticks(axes, label_font_size, dark)
    #set x and y label offsets -- currently though rcParams, but should be
    #able to do at instance level?
    #rc('xtick.major', pad=label_offset)
    #rc('ytick.major', pad=label_offset)
    return font, label_font_size

def init_ticks(axes=None, label_font_size=None, dark=False):
    """Initializes ticks for fingerprint plots or other plots ranging from 0-1.
    
    takes axis argument a from a = gca(), or a specified axis 
    sets the ticks to span from 0 to 1 with .1 intervals
    changes the size of the ticks and the corresponding number labels
    """
    if axes is None:
        axes = gca()
    axes.set_xticks(arange(0,1.01,.1),)
    axes.set_yticks(arange(0,1.01,.1))

    #reset sizes for x and y labels
    x = axes.get_xticklabels()
    y = axes.get_yticklabels()
    if label_font_size is not None:
        for l in axes.get_xticklabels() + axes.get_yticklabels():
            l.set_fontsize(label_font_size)
    #if dark, need to reset color of internal ticks to white
    if dark:
        for l in axes.get_xticklines() + axes.get_yticklines():
            l.set_markeredgecolor('white')

def set_axis_to_probs(axes=None):
    """sets the axes to span from 0 to 1.

    Useful for forcing axes to range over probabilities. Axes are
    sometimes reset by other calls.
    """
    #set axis for probabilities (range 0 to 1)
    if axes is None:
        axes = gca()
    axes.set_xlim([0,1])
    axes.set_ylim([0,1])

def plot_regression_line(x,y,line_color='r', axes=None, prob_axes=False, \
    axis_range=None):
    """Plots the regression line, and returns the equation.
    
    x and y are the x and y data for a single series
    line_color is a matplotlib color, will be used for the line
    axes is the name of the axes the regression will be plotted against
    prob_axes, if true, forces the axes to be between 0 and 1
    range, if not None, forces the axes to be between (xmin, xmax, ymin, ymax).
    """
    if axes is None:
        axes = gca()
    m, b = regress(x, y)
    r, significance = correlation(x,y)
    #set the a, b, and r values. a is the slope, b is the intercept.
    r_str = '%0.3g'% (r**2)
    m_str ='%0.3g' % m
    b_str = '%0.3g' % b

    #want to clip the line so it's contained entirely within the graph
    #coordinates. Basically, we need to find the values of y where x
    #is at x_min and x_max, and the values of x where y is at y_min and
    #y_max.

    #if we didn't set prob_axis or axis_range, just find empirical x and y
    if (not prob_axes) and (axis_range is None):
       x1, x2 = min(x), max(x)
       y1, y2 = m*x1 + b, m*x2 + b
       x_min, x_max = x1, x2
    else:
        if prob_axes:
            x_min, x_max = 0, 1
            y_min, y_max = 0, 1
        else: #axis range must have been set
            x_min, x_max, y_min, y_max = axis_range
        #figure out bounds for x_min and y_min
        y_at_x_min = m*x_min + b
        if y_at_x_min < y_min:  #too low: find x at y_min
            y1 = y_min
            x1 = (y_min-b)/m
        elif y_at_x_min > y_max: #too high: find x at y_max
            y1 = y_max
            x1 = (y_max-b)/m
        else:   #just right
            x1, y1 = x_min, y_at_x_min

        y_at_x_max = m*x_max + b
        if y_at_x_max < y_min:  #too low: find x at y_min
            y2 = y_min
            x2 = (y_min-b)/m
        elif y_at_x_max > y_max: #too high: find x at y_max
            y2 = y_max
            x2 = (y_max-b)/m
        else:   #just right
            x2, y2 = x_max, y_at_x_max

        #need to check that the series wasn't entirely in range
    if (x_min <= x1 <= x_max) and (x_min <= x2 <= x_max):
        axes.plot([x1,x2],[y1,y2], color=line_color, linewidth=0.5)

    if b >= 0:
        sign_str = ' + '
    else:
        sign_str = ' '
    
    equation=''.join(['y= ',m_str,'x',sign_str,b_str,'\nr$^2$=',r_str])
    return equation, line_color

def add_regression_equations(equations, axes=None, prob_axes=False, \
    horizontalalignment='right', verticalalignment='bottom'):
    """Writes list of regression equations to graph.

    equations: list of regression equations

    size: size of the graph in inches
    """
    if axes is None:
        axes = gca()
    if prob_axes:
        min_x, max_x = 0, 1
        min_y, max_y = 0, 1
    else:
        min_x, max_x = axes.get_xlim()
        min_y, max_y = axes.get_ylim()
    x_range = abs(max_x - min_x)
    y_range = abs(max_y - min_y)

    for i, (eq_text, eq_color) in enumerate(equations):
        axes.text((x_range * 0.98) + min_x, \
            (y_range * 0.02 + min_y +(y_range * .1 * i)), \
            str(eq_text), \
            horizontalalignment=horizontalalignment, \
            verticalalignment=verticalalignment, \
            color=eq_color)

def broadcast(i, n):
    """Broadcasts i to a vector of length n."""
    try:
        i = list(i)
    except:
        i = [i]
    reps, leftovers = divmod(n, len(i))
    return (i * reps) + i[:leftovers]
    
    
#scatterplot functions and helpers

def plot_scatter(data, series_names=None, \
    series_color=standard_series_colors, line_color=standard_series_colors,\
    alpha=0.25, marker_size=.015, scale_markers=True,
    show_legend=True,legend_loc='center right',
    show_regression=True, show_equation=True,
    prob_axes=False, size=8.0, axes=None,
    **kwargs):
    """helper plots one or more series of scatter data of specified color,
    calls the initializing functions, doesn't print graph
    
    takes: plotted_pairs, series_names, show_legend, legend_loc, and
        **kwargs passed on to init_graph_display (these include title,
        aux_title, size, graph_shape, graph_grid, x_label, y_label,
        dark, with_parens).
                 
    plotted_pairs = (first_pos, second_pos, dot_color, line_color,
    alpha, show_regression, show_equation)

    returns the regression str equation (list) if regression is set true

    suppresses legend if series not named, even if show_legend is True.
    """
    if not axes:
        axes = gca()
    #initialize fonts, shape and labels
    font,label_font_size=init_graph_display(prob_axes=prob_axes, \
        size=size, axes=axes, **kwargs)
    equations = []
    #figure out how many series there are, and scale vals accordingly
    num_series = int(len(data)/2)
    series_color = broadcast(series_color, num_series)
    line_color = broadcast(line_color, num_series)
    alpha = broadcast(alpha, num_series)
    marker_size = broadcast(marker_size, num_series)
    if scale_markers:
        marker_size = [frac_to_psq(m, size) for m in marker_size]
    
    series = []
    for i in range(num_series):
        x, y = data[2*i], data[2*i+1]
        series.append(axes.scatter(x,y,s=marker_size[i],c=series_color[i],\
        alpha=alpha[i]))
        #find the equation and plots the regression line if True
        if show_regression:
            equation = plot_regression_line(x,y,line_color[i], axes=axes, \
                prob_axes=prob_axes)
        if show_equation:
            equations.append(equation)  #will be (str, color) tuple
    #update graph size for new data
    axes.autoscale_view(tight=True)
    #print all the regression equations at once -- need to know how many
    if show_regression:
        add_regression_equations(equations, axes=axes, prob_axes=prob_axes)
    #clean up axes if necessary
    if show_legend and series_names: #suppress legend if series not named
        axes.legend(series, series_names, legend_loc)

    if prob_axes:
        set_axis_to_probs(axes)
    return equations, font

#Contour plots and related functions

def plot_filled_contour(plot_data, xy_data=None, show_regression=False, \
    show_equation=False, fill_cmap=cm.hot, graph_shape='rect', \
    num_contour_lines=10, prob_axes=False, **kwargs):
    """helper plots one or more series of contour data
    calls the initializing functions, doesn't output figure
    
    takes: plot_data, xy_data, show_regression, show_equation, fill_cmap, 
    and **kwargs passed on to init_graph_display.
                 
           plot_data = (x_bin, y_bin, data_matrix dot_colors)
    """
    if show_regression:
        equation = plot_regression_line(xy_data[:,0],xy_data[:,1], \
            prob_axes=prob_axes)
        if show_equation:
            add_regression_equations([equation])
    #init graph display, rectangular due to needed colorbar space
    init_graph_display(graph_shape=graph_shape, **kwargs)
    #plots the contour data
    for x_bin,y_bin,data_matrix in plot_data:
        contourf(x_bin,y_bin,data_matrix, num_contour_lines, cmap=fill_cmap)
    #add the colorbar legend to the side
    colorbar()

def plot_contour_lines(plot_data, xy_data=None, show_regression=False, \
        show_equation=False, smooth_steps=0, num_contour_lines=10, \
        label_contours=False, line_cmap=cm.hot, fill_cmap=cm.gray,dark=True,
        graph_shape='rect', prob_axes=False, **kwargs):
    """helper plots one or more series of contour line data
    calls the initializing functions, doesn't output figure
    
    takes: plot_data, xy_data, show_regression, show_equation, smooth,
        num_contour_lines, label_contours, line_cmap, fill_cmap, graph_shape,
        and **kwargs passed on to init_graph_display.
                 
           plot_data = (x_bin, y_bin, data_matrix dot_colors)
    """
    if prob_axes:
        extent = (0,1,0,1)
    else:
        a = gca()
        extent = a.get_xlim()+a.get_ylim()
    #init graph display, rectangular due to needed colorbar space
    init_graph_display(graph_shape=graph_shape,
        dark=dark, **kwargs)
    #plots the contour data
    for x_bin,y_bin,data in plot_data:
        orig_max = max(ravel(data))
        scaled_data = (data/orig_max*255).astype('b')
        if smooth_steps and (Image is not None):
            orig_shape = data.shape
            im = Image.fromstring('L', data.shape, scaled_data)
            for i in range(smooth_steps):
                im = im.filter(ImageFilter.BLUR)
            new_data = fromstring(im.tostring(), 'b')
            data = reshape(new_data.astype('i')/255.0 * orig_max, orig_shape)
        
        if fill_cmap is not None:
            im = imshow(data, interpolation='bicubic', extent=extent, \
                origin='lower', cmap=fill_cmap)
        result=contour(x_bin,y_bin,data, num_contour_lines,
                              origin='lower',linewidth=2,
                              extent=extent, cmap=line_cmap)
        if label_contours:
            clabel(result, fmt='%1.1g')

    #add the colorbar legend to the side
    cb = colorbar()
    cb.ax.axisbg = 'black'

    if show_regression:
        equation=plot_regression_line(xy_data[0],xy_data[1],prob_axes=prob_axes)
        if show_equation:
            add_regression_equations([equation])

def plot_histograms(data, graph_name='histogram.png', bins=20,\
        normal_fit=True, normed=True, colors=None, linecolors=None, \
        alpha=0.75, prob_axes=True, series_names=None, show_legend=False,\
        y_label=None, **kwargs):
    """Outputs a histogram with multiple series (must provide a list of series).
    
    takes:  data: list of arrays of values to plot (needs to be list of arrays
            so you can pass in arrays with different numbers of elements)

            graph_name: filename to write graph to
            bins: number of bins to use
            normal_fit: whether to show the normal curve best fitting the data
            normed: whether to normalize the histogram (e.g. so bars sum to 1)
            colors: list of colors to use for bars
            linecolors: list of colors to use for fit lines

            **kwargs are pssed on to init_graph_display.

    """
    rc('patch', linewidth=.2)
    if y_label is None:
        if normed:
            y_label='Frequency'
        else:
            y_label='Count'
    num_series = len(data)
    if colors is None:
        if num_series == 1:
            colors = ['white']
        else:
            colors = standard_series_colors
    if linecolors is None:
        if num_series == 1:
            linecolors = ['red']
        else:
            linecolors = standard_series_colors
    
    init_graph_display(prob_axes=prob_axes, y_label=y_label, **kwargs)
    all_patches = []
    for i, d in enumerate(data):
        fc = colors[i % len(colors)]
        lc = linecolors[i % len(linecolors)]
        
        counts, x_bins, patches = hist(d, bins=bins, normed=normed, \
            alpha=alpha, facecolor=fc)

        all_patches.append(patches[0])

        if normal_fit and len(d) > 1:
            maxv, minv = max(d), min(d)
            mu = mean(d)
            sigma = std(d)
            bin_width = x_bins[-1] - x_bins[-2]
            #want normpdf to extend over the range
            normpdf_bins = arange(minv,maxv,(maxv - minv)*.01)
            y = normpdf(normpdf_bins, mu, sigma)
            orig_area = sum(counts) * bin_width
            y = y * orig_area   #normpdf area is 1 by default
            plot(normpdf_bins, y, linestyle='--', color=lc, linewidth=1)

    if show_legend and series_names:
        fp = FontProperties()
        fp.set_size('x-small')
        legend(all_patches, series_names, prop = fp) 
    
   #output figure if graph name set -- otherwise, leave for further changes
    if graph_name is not None:
        savefig(graph_name)

def plot_monte_histograms(data, graph_name='gene_histogram.png', bins=20,\
        normal_fit=True, normed=True, colors=None, linecolors=None, \
        alpha=0.75, prob_axes=True, series_names=None, show_legend=False,\
        y_label=None, x_label=None, **kwargs):
    """Outputs a histogram with multiple series (must provide a list of series).

    Differs from regular histogram in that p-value works w/exactly two 
    datasets, where the first dataset is the reference set. Calculates the
    mean of the reference set, and compares this to the second set (which is
    assumed to contain the means of many runs producing data comparable to the
    data in the reference set).

    takes:  data: list of arrays of values to plot (needs to be list of arrays
            so you can pass in arrays with different numbers of elements)

            graph_name: filename to write graph to
            bins: number of bins to use
            normal_fit: whether to show the normal curve best fitting the data
            normed: whether to normalize the histogram (e.g. so bars sum to 1)
            colors: list of colors to use for bars
            linecolors: list of colors to use for fit lines

            **kwargs are passed on to init_graph_display.
    """
    rc('patch', linewidth=.2)
    rc('font', size='x-small')

    rc('axes', linewidth=.2)
    rc('axes', labelsize=7)
    rc('xtick', labelsize=7)
    rc('ytick', labelsize=7)

    if y_label is None:
        if normed:
            y_label='Frequency'
        else:
            y_label='Count'
    num_series = len(data)
    if colors is None:
        if num_series == 1:
            colors = ['white']
        else:
            colors = standard_series_colors
    if linecolors is None:
        if num_series == 1:
            linecolors = ['red']
        else:
            linecolors = standard_series_colors
   
    init_graph_display(prob_axes=prob_axes, y_label=y_label, **kwargs)
    all_patches = []
    
    for i, d in enumerate(data):
        fc = colors[i % len(colors)]
        lc = linecolors[i % len(linecolors)]
        
        counts, x_bins, patches = hist(d, bins=bins, normed=normed, \
            alpha=alpha, facecolor=fc)

        all_patches.append(patches[0])

        if normal_fit and len(d) > 1:
            mu = mean(d)
            sigma = std(d)
            minv = min(d)
            maxv = max(d)
            bin_width = x_bins[-1] - x_bins[-2]
            
            #set range for normpdf
            normpdf_bins = arange(minv,maxv,0.01*(maxv-minv))
            y = normpdf(normpdf_bins, mu, sigma)
            orig_area = sum(counts) * bin_width
            y = y * orig_area   #normpdf area is 1 by default
            plot(normpdf_bins, y, linestyle='--', color=lc, linewidth=1)
            font = { 'color': lc, 
                     'fontsize': 11}
            text(mu, 0.0 , "*", font, verticalalignment='center',
                             horizontalalignment='center')

    xlabel(x_label)
    if show_legend and series_names:

        fp = FontProperties()
        fp.set_size('x-small')

        legend(all_patches, series_names, prop = fp) 
    
   #output figure if graph name set -- otherwise, leave for further changes
    if graph_name is not None:
        savefig(graph_name)

def plot_scatter_with_histograms(data, graph_name='histo_scatter.png', \
    graph_grid='/', prob_axes=False, bins=20, frac=0.9, scatter_alpha=0.5, \
    hist_alpha=0.8, colors=standard_series_colors, normed='height', **kwargs):
    """Plots a scatter plot with histograms showing distribution of x and y.

    Data should be list of [x1, y1, x2, y2, ...].
    """

    #set up subplot coords
    tl=subplot(2,2,1)
    br=subplot(2,2,4)
    bl=subplot(2,2,3, sharex=tl, sharey=br)

    #get_position returns a Bbox relative to figure
    tl_coords = tl.get_position()
    bl_coords = bl.get_position()
    br_coords = br.get_position()

    left = tl_coords.xmin
    bottom = bl_coords.ymin

    width = br_coords.xmax - left
    height = tl_coords.ymax - bottom

    bl.set_position([left, bottom, frac*width, frac*height])
    tl.set_position([left, bottom+(frac*height), frac*width, (1-frac)*height])
    br.set_position([left+(frac*width), bottom, (1-frac)*width, frac*height])

    #suppress frame and axis for histograms
    for i in [tl,br]:
        i.set_frame_on(False)
        i.xaxis.set_visible(False)
        i.yaxis.set_visible(False)
    
    plot_scatter(data=data, alpha=scatter_alpha, axes=bl, **kwargs)
    
    for i in range(0, len(data), 2):
        x, y = data[i], data[i+1]
        color = colors[int((i/2))%len(colors)]

        hist(x, facecolor=color, bins=bins, alpha=hist_alpha, normed=normed, axes=tl)
        hist(y, facecolor=color, bins=bins, alpha=hist_alpha, normed=normed, \
        axes=br, orientation='horizontal')
    if prob_axes:
        bl.set_xlim(0,1)
        bl.set_ylim(0,1)
        br.set_ylim(0,1)
        tl.set_xlim(0,1)
  
   #output figure if graph name set -- otherwise, leave for further changes
    if graph_name is not None:
        savefig(graph_name)

def format_contour_array(data, points_per_cell=20, bulk=0.8):
    """Formats [x,y] series of data into x_bins, y_bins and data for contour().

    data: 2 x n array of float representing x,y coordinates
    
    points_per_cell: average points per unit cell in the bulk of the data,
    default 3
    
    bulk: fraction containing the 'bulk' of the data in x and y, default
    0.8 (i.e. 80% of the data will be used in the calculation).
    
    returns: x-bin, y-bin, and a square matrix of frequencies to be plotted

    WARNING: Assumes x and y are in the range 0-1.
    """
    #bind x and y data
    data_x = sort(data[0]) #note: numpy sort returns a sorted copy
    data_y = sort(data[1])
    num_points = len(data_x)

    #calculate the x and y bounds holding the bulk of the data
    low_prob = (1-bulk)/2.0
    low_tail = int(num_points*low_prob)
    high_tail = int(num_points*(1-low_prob))
    
    x_low = data_x[low_tail]
    x_high = data_x[high_tail]
    y_low = data_y[low_tail]
    y_high = data_y[high_tail]
    
    #calculate the side length in the bulk that holds the right number of
    #points
    delta_x = x_high - x_low
    delta_y = y_high - y_low
    
    points_in_bulk = num_points * bulk #approximate: assumes no correlation
    area_of_bulk = delta_x * delta_y
    points_per_area = points_in_bulk/area_of_bulk
    side_length = sqrt(points_per_cell / points_per_area)
    #correct the side length so we get an integer number of bins.
    num_bins = int(1/side_length)
    corrected_side_length = 1.0/num_bins
    
    #figure out how many items are in each grid square in x and y
    #
    #this is the tricky part, because contour() takes as its data matrix
    #the points at the vertices of each cell, rather than the points at
    #the centers of each cell. this means that if we were going to make
    #a 3 x 3 grid, we actually have to estimate a 4 x 4 matrix that's offset
    #by half a unit cell in both x and y.
    #
    #if the data are between 0 and 1, the first and last bin in our range are
    #superfluous because searchsorted will put items before the first
    #bin into bin 0, and items after the last bin into bin n+1, where
    #n is the maximum index in the original array. for example, if we
    #have 3 bins, the values .33 and .66 would suffice to find the centers, 
    #because anything below .33 gets index 0 and anything above .66 gets index 
    #2 (anything between them gets index 1). incidentally, this prevents 
    #issues with floating-point error and values slightly below 0 or above 1 
    #that might otherwise arise.
    #
    #however, for our 3 x 3 case, we actually want to start estimating at the
    #cell centered at 0, i.e. starting at -.33/2, so that we get the four
    #estimates centered at (rather than starting at) 0, .33, .66, and 1.
    #because the data are constrained to be between 0 and 1, we will need to
    #double the counts at the edges (and quadruple them at the corners) to get
    #a fair estimate of the density.
    csl = corrected_side_length #save typing below 
    eps = csl/10 #don't ever want max value to be in the list precisely
    half_csl = .5*csl
    bins = arange(half_csl, 1+half_csl-eps, csl)
    x_coords = searchsorted(bins, data[0])
    y_coords = searchsorted(bins, data[1])
    #matrix has dimension 1 more than num bins, b/c can be above largest
    matrix = zeros((num_bins+1, num_bins+1))
    #for some reason, need to swap x and y to match up with normal
    #scatter plots
    for coord in zip(y_coords, x_coords):
        matrix[coord] += 1
    
    #we now have estimates of the densities at the edge of each of the
    #n x n cells in the grid. for example, if we have a 3 x 3 grid, we have
    #16 densities, one at the center of each grid cell (0, .33, .66, 1 in each
    #dimension). need to double the counts at edges to reflect places where
    #we can't observe data because of range restrictions.
    matrix[0]*=2
    matrix[:,0]*=2
    matrix[-1]*=2
    matrix[:,-1]*=2
    #return adjusted_bins as centers, rather than boundaries, of the range
    x_bins = csl*arange(num_bins+1)
    return x_bins, x_bins, matrix


if __name__ == '__main__':
    from numpy.random import normal
    x = normal(0.3, 0.05, 1000)
    y = normal(0.5, 0.1, 1000)
    plot_scatter_with_histograms([x,x+y, y, (x+y)/2], prob_axes=True)
