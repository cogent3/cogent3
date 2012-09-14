#/usr/bin/env python
"""Provides different kinds of codon usage plots.

See individual docstrings for more info.
"""
from matplotlib import use, rc
use('Agg')  #suppress graphical rendering
rc('text', usetex=True)
rc('font', family='serif')  #required to match latex text and equations
from cogent.core.usage import UnsafeCodonUsage as CodonUsage
from cogent.draw.util import scatter_classic, \
   init_graph_display, init_ticks, set_axis_to_probs, \
    broadcast, plot_scatter, plot_filled_contour, \
    plot_contour_lines, standard_series_colors 
from pylab import plot, savefig, gca, text, figlegend

__author__ = "Stephanie Wilson"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Stephanie Wilson"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

#module-level constants

#historical doublet order for fingerprint plot; not currently used, but
#same order that the colors were entered in. Matches Sueoka 2002.
doublet_order = ['GC','CG','GG','CU','CC','UC','AC','GU','UU','CA','AU',\
                 'AA','AG','GA','UA','UG']
color_order = ["#000000","#FF0000","#00FF00","#FFFF00",
          "#CC99FF","#FFCC99","#CCFFFF","#C0C0C0",
          "#6D6D6D","#2353FF","#00FFFF","#FF8800",
          "#238853","#882353","#EC008C","#000099"]
#map doublets to colors so we can make sure the same doublet always
#gets the same colors
doublets_to_colors = dict(zip(doublet_order, color_order))
#creates a dictionary for the amino acid labels, less to input
aa_labels={'ALANINE':'GCN', 'ARGININE4':'CGN', 'GLYCINE':'GGN',
          'LEUCINE4':'CTN', 'PROLINE':'CCN', 'SERINE4':'TCN',
          'THREONINE':'ACN', 'VALINE':'GTN'}
standard_series_colors=['k','r','g','b', 'm','c']

#scatterplot functions and helpers
def plot_cai_p3_scatter(data, graph_name='cai_p3_scat.png', **kwargs):
    """Outputs a CAI vs P3 scatter plot.

    expects data as ([P3s_1, CAIs_1, P3s_2, CAIs_2, ...])
    """
    plot_scatter(data, graph_shape='sqr', graph_grid=None,\
        x_label="$P_3$",y_label="CAI", prob_axes=True,**kwargs)
    savefig(graph_name)

def plot_p12_p3(data, graph_name='p12_p3.png', **kwargs):
    """Outputs a P12 versus P3 scatter graph, optionally including regression.

    expects data as [P3_1, P12_1, P3_2, P12_2, ...n ].
    """
    plot_scatter(data, graph_shape='sqr', graph_grid='/',\
        x_label="$P_3$",y_label="$P_{12}$", prob_axes=True, **kwargs)
    savefig(graph_name)

def plot_p123_gc(data, graph_name='p123_gc.png', use_p3_as_x=False, **kwargs):
    """Output a scatter plot of p1,p2,p3 vs gc content
    
    Expects data as array with rows as GC, P1, P2, P3
    p1=blue, p2=green, p3=red

    """
    #unpack common x axis, and decide on series names
    if use_p3_as_x:
        series_names = ['$P_1$', '$P_2$']
        colors=['b','g']
        x_label='$P_3$'
        y_label='$P_{12}$'
        xy_pairs = [data[3], data[1], data[3], data[2]]
    else:
        series_names = ['$P_1$', '$P_2$', '$P_3$']
        colors=['b','g','r']
        x_label='GC'
        y_label='$P_{123}$'
        xy_pairs = [data[0], data[1], data[0], data[2], data[0], data[3]]
    
    #plot points and write graph
    plot_scatter(xy_pairs, graph_grid='/',x_label=x_label,y_label=y_label,
        series_names=series_names, prob_axes=True, **kwargs)
    savefig(graph_name)

def plot_fingerprint(data, alpha=0.7, \
    show_legend=True, graph_name='fingerprint.png', has_mean=True,
    which_blocks='quartets', multiple=False, graph_grid='t', prob_axes=True, \
    edge_colors='k', **kwargs):
    """Outputs a bubble plot of four-codon amino acid blocks
    labeled with the colors from Sueoka 2002.

    takes: data:  array-elements in the col order x, y, r of
           each of the four codon Amino Acids in the row order:
           ALA, ARG4, GLY, LEU4, PRO, SER, THR, VAL
           (for traditional fingerprint), or:
           UU -> GG (for 16-block fingerprint).
           last row is the mean (if has_mean is set True)

        **kwargs passed on to init_graph_display (these include 
        graph_shape, graph_grid, x_label, y_label, dark, with_parens).
                 
           title: will be printed on graph (default: 'Unknown Species')
           
           num_genes (number of genes contributing to graph: default None)
           NOTE: will not print if None.)
        
           size: of graph in inches (default = 8.0)

           alpha: transparency of bubbles
           (ranges from 0, transparent, to 1, opaque; default 0.7)
           
           show_legend: bool, default True, whether to print legend

           graph_name: name of file to write (default 'fingerprint.png')

           has_mean: whether the data contain the mean (default: True)

           which_blocks: which codon blocks to print (default is 'quartets'
           for the 4-codon amino acid blocks, but can also use 'all' for all 
           quartets or 'split' for just the split quartets.)

           multiple: if False (the default), assumes it got a single block
           of data. Otherwise, assumes multiple blocks of data in a list or
           array.

           edge_colors: if multiple is True (ignored otherwise), uses this
           sequence of edge color strings to hand out edge colors to successive
           series. Will iterate over this, so can be a string of 1-letter
           color codes or a list of color names.

    note: that the data are always expected to be in the range (0,1)
    since we're plotting frequencies. axes, gid, etc. are hard-coded
    to these values. 
    """
    #figure out which type of fingerprint plot we're doing, and get the
    #right colors
    if which_blocks == 'quartets':
        blocks = CodonUsage.SingleAABlocks
    elif which_blocks == 'split':
        blocks = CodonUsage.SplitBlocks
    else:
        blocks = CodonUsage.Blocks

    colors = [doublets_to_colors[i] for i in blocks]
      
    #formatting the labels in latex
    x_label="$G_3/(G_3+C_3)$"
    y_label="$A_3/(A_3+T_3)$"

    #initializing components of the graph
    font,label_font_size=init_graph_display(graph_shape='sqr', \
        graph_grid=graph_grid, x_label=x_label, \
        y_label=y_label, prob_axes=prob_axes, **kwargs)

    if not multiple:
        data = [data]
 
    alpha = broadcast(alpha, len(data))
    edge_colors = broadcast(edge_colors, len(data))
  
    for al, d, edge_color in zip(alpha, data, edge_colors):
        #skip this series if no data
        if d is None or not d.any():
            continue
        for i, color in enumerate(colors):
            j = i+1
            #note: doing these as slices because scatter_classic needs the
            #extra level of nesting
            patches = scatter_classic(d[i:j,0], d[i:j,1],
                        s=(d[i:j,2]/2), c=color)
            #set alpha for the patches manually
            for p in patches:
                p.set_alpha(al)
                p.set_edgecolor(edge_color)
        
        #plot mean as its own point -- can't do cross with scatter
        if has_mean:
            mean_index = len(blocks)    #next index after the blocks
            plot([d[mean_index,0]], [d[mean_index,1]],
                 '-k+',markersize=label_font_size, alpha=al)
               

    abbrev = CodonUsage.BlockAbbreviations

    a = gca()
    #if show_legend is True prints a legend in the right center area
    if show_legend:
        legend_key = [abbrev[b] for b in blocks]
        #copy legend font properties from the x axis tick labels
        legend_font_props = \
            a.xaxis.get_label().get_fontproperties().copy()
        legend_font_scale_factor = 0.7
        curr_size = legend_font_props.get_size()
        legend_font_props.set_size(curr_size*legend_font_scale_factor)
        l = figlegend(a.patches[:len(blocks)],
                  legend_key,
                  prop=legend_font_props,
                  loc='center right',borderpad=0.1,labelspacing=0.5,
                  handlelength=1.0,handletextpad=0.5, borderaxespad=0.0)
        #fix transparency of patches
        for p in l.get_patches():
            p.set_alpha(1)

    #initialize the ticks
    set_axis_to_probs()
    init_ticks(a, label_font_size)
    a.set_xticks([0, 0.5, 1])
    a.set_yticks([0,0.5,1])
    
    #output the figure
    if graph_name is not None:
        savefig(graph_name)

#Contour plots and related functions

def plot_cai_p3_contour(x_bin,y_bin,data,xy_data,
                        graph_name='cai_contour.png',
                        prob_axes=True, **kwargs):
    """Output a contour plot of cai vs p3 with colorbar on side

    takes: x_bin, y_bin, data (data matrix)
    
           label (default 'Unknown Species')

           num_genes (default 0 will not print, other numbers will)

           size: of graph in inches (default = 8.0)

           graph_name: default 'cai_contour.png'
    """
    plot_data =[(x_bin,y_bin,data)]
    plot_filled_contour(plot_data, graph_grid='/',x_label="$P_3$", \
        y_label="CAI", prob_axes=prob_axes, **kwargs)
    set_axis_to_probs()
    if graph_name is not None:
        savefig(graph_name)

def plot_cai_p3_contourlines(x_bin,y_bin,data,xy_data,
                             graph_name='cai_contourlines.png',
                             prob_axes=True, **kwargs):
    """Output a contour plot of cai
    
    takes: x_bin, y_bin, data (data matrix)
    
           label (default 'Unknown Species')

           num_genes (default 0 will not print, other numbers will)

           size: of graph in inches (default = 8.0)

           graph_name: default 'cai_contourlines.png'
    """
    plot_data =[(x_bin,y_bin,data)]
    plot_contour_lines(plot_data, graph_grid='/', x_label="$P_3$", \
        y_label="CAI", prob_axes=prob_axes,**kwargs)
    if graph_name is not None:
        savefig(graph_name)

def plot_p12_p3_contour(x_bin,y_bin,data,xy_data,
                        graph_name='p12_p3_contour.png',
                        prob_axes=True, **kwargs):
    """Outputs a P12 versus P3 contour graph
    and the mean equation of the plot

    takes: x_bin, y_bin, data (data matrix)
    
           label (default 'Unknown Species')

           num_genes (default 0 will not print, other numbers will)

           size: of graph in inches (default = 8.0)

           graph_name: default 'p12_p3_contourlines.png'
    """
    plot_data =[(x_bin,y_bin,data)]
    plot_filled_contour(plot_data, graph_grid='/', x_label="$P_3$", \
        y_label="$P_{12}$", prob_axes=prob_axes,**kwargs)
    set_axis_to_probs()
    if graph_name is not None:
        savefig(graph_name)

def plot_p12_p3_contourlines(x_bin,y_bin,data,xy_data, prob_axes=True,\
    graph_name='p12_p3_contourlines.png', **kwargs):
    """Outputs a P12 versus P3 contourline graph
    and the mean equation of the plot

    takes: x_bin, y_bin, data (data matrix)
    
           label (default 'Unknown Species')

           num_genes (default 0 will not print, other numbers will)

           size: of graph in inches (default = 8.0)

           graph_name: default 'p12_p3_contourlines.png
    """
    plot_data =[(x_bin,y_bin,data)]
    plot_contour_lines(plot_data, graph_grid='/', x_label="$P_3$",\
        y_label="$P_{12}$", prob_axes=prob_axes, **kwargs)
    set_axis_to_probs()
    if graph_name is not None:
        savefig(graph_name)

#Other graphs

def plot_pr2_bias(data, title='ALANINE', graph_name='pr2_bias.png', \
    num_genes='ignored', **kwargs):
    """Outputs a PR2-Bias plot of:
    -isotypic transversions (base swapping)
    with G3/(G3+C3) and A3/(A3+T3)
    -Transitions (deaminations)
    with G3/(G3+A3) and C3/(C3+T3)
    -Allotypic transversions (G- oxidations)
    with G3/(G3+T3) and C3/(C3+A3)

    takes: an array in the order: x,G3/(G3+C3),A3/(A3+T3),
    G3/(G3/A3),C3/(C3+T3),G3/(G3+T3),C3/(C3+A3)

    label: default 'ALANINE'
    one amino acid written out in caps:
    ALANINE, ARGININE4, GLYCINE, LEUCINE4,
    PROLINE, SERINE4, THREONINE, VALINE
       from one of the amino acids program will add acronym
       C2 type: ala(GCN), pro(CCN), ser4(TCN), thr(ACN)
       G2 type: arg4 (CGN), an gly(GGN)
       T2 type: leu4(CTN), val (GTN)

    size: of graph in inches (default = 8.0)

    graph_name: default 'pr2_bias.png'
    
    num_genes: number of genes contributing to graph, currently ignored.
    """
    #we can't put anything in the top right, so print num_genes after the title
    #if it was supplied
    #initializes the graph display and font
    font,label_font_size=init_graph_display(graph_shape='sqr', \
        graph_grid='/', x_label="$P_3$", y_label="Y axis", prob_axes=True, \
        title=title, **kwargs)
    #sets the marker_size relative to the font and thus the graph size
    marker_size = (label_font_size-1)
    
    #plots the pr2bias in order G3/(G3+C3),A3/(A3+T3),
    #                           G3/(G3/A3),C3/(C3+T3),
    #                           G3/(G3+T3),C3/(C3+A3)
    #colors and symbols coded from Sueoka 2002
    plot(data[:,0], data[:,1], '-ko', c='k',
         markersize=marker_size)
    plot(data[:,0], data[:,2], '-kv', c='k',
         markersize=marker_size)
    plot(data[:,0], data[:,3], '-ro', c='r',
         markersize=marker_size)
    plot(data[:,0], data[:,4], '-rv', c='r',
         markersize=marker_size)
    plot(data[:,0], data[:,5], '-wo', c='k', mfc='w',
         markersize=marker_size)
    plot(data[:,0], data[:,6], '-wv', c='k', mfc='w',
         markersize=marker_size)

    #aaLabel based on the amino acid that is graphed
    #C2 type: ala(GCN), pro(CCN), ser4(TCN), thr(ACN)
    #G2 type: arg4 (CGN), an gly(GGN)
    #T2 type: leu4(CTN), val (GTN) (Sueoka 2002)
    text(.95, .05, aa_labels[title], font, verticalalignment='bottom',
         horizontalalignment='right')

    #output the figure
    set_axis_to_probs()
    if graph_name is not None:
        savefig(graph_name)
