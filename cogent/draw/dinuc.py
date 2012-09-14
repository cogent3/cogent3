#/usr/bin/env python
"""Provides a plot for visualizing dinucleotide frequencies..
"""
from matplotlib import use, rc
from matplotlib.font_manager import FontProperties
use('Agg')  #suppress graphical rendering
from matplotlib.ticker import FixedFormatter
from pylab import arange, axvline, array, plot, scatter, legend, gca, xlim, \
    savefig

__author__ = "Stephanie Wilson"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Stephanie Wilson"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

rc('text', usetex=True)
rc('font', family='serif')  #required to match latex text and equations

colors = ["#000000","#FF0000","#00FF00","#FFFF00",
              "#CC99FF","#FFCC99","#CCFFFF","#C0C0C0",
              "#6D6D6D","#2353FF","#00FFFF","#FF8800",
              "#238853","#882353","#EC008C","#000099"]
axis_names= [a+b for a in 'UCAG' for b in 'UCAG']

#place to put the averages, and the genes on x axis
line_positions = arange(16)

#standard_gene_type_formats provides arguments to plot() that depend on the type
#of gene. All contents of this dict must be valid arguments and values for
#plot()!
standard_gene_type_formats = {
    'hgt':{'marker':'d'},
    'ribosomal':{'marker':'s'},
    None:{'marker':'o'},
}

#offsets maps labels onto horizontal positions within each column. Offset of
#0.5 means it will appear in the middle of the column.
offsets = {
    'hgt': 0.3,
    'ribosomal':0.7,
    None: 0.5
}

light_gray = '#CCCCCC' #define standard color for background lines

def dinuc_plot(data, graph_name="ADinucTester.png", title = "", \
        gene_type_formats=standard_gene_type_formats, \
        background_line_color=light_gray, avg_formats={},
        point_formats={}):
    """Returns dinucleotide usage plot from given species data.

    Data is a dict with the following structure:

    {name: {gene_type: [values]}}

    where name is e.g. the name of the species (will be displayed on graph),
    and gene_type will be the name of the type of gene (e.g. ribosomal,
    hgt, non-hgt, all).

    Values should be the 16-element array obtained from calling normalize() on
    the DinucUsage objects giving the DinucUsage frequencies, as fractions,
    for each of the 16 dinucleotides in alphabet order (TT -> GG).

    Calculates the average of each gene type on the fly and connects these
    averages with lines. If you have precalculated the average, just pass it
    in as a single 'gene' (the average of an array with one item will be 
    itself...).

    avg_formats should be passable to plot(), e.g. markersize.
    point_formats should be passable to scatter(), e.g. 'c'.
    """
    #constructing the vertical lines
    for a in line_positions:
        axvline(a,ymin=0,ymax=1, color=background_line_color)
    
    curr_color_index = 0
    lines = []
    labels = []
    for label, gene_types in sorted(data.items()):
        gene_types = data[label]
        curr_color = colors[curr_color_index]
        for type_label, usages in sorted(gene_types.items()):
            labels.append(':'.join(map(str, [label,type_label])))
            x_positions = line_positions + \
                offsets.get(type_label, offsets[None])
            curr_format = gene_type_formats.get(type_label, \
                gene_type_formats[None])
            #plot the averages...connected with a line for easier visual
            avg_format = curr_format.copy()
            avg_format.update(avg_formats)
            avg = sum(array(usages))/len(usages)
            curr_line = plot(x_positions, avg, color=curr_color,
                 markerfacecolor=curr_color, label=label,
                 **avg_format)
            lines.append(curr_line)
            point_format = curr_format.copy()
            if len(usages) > 1: #got more than one gene: plot points
                point_format.update(point_formats)
                for u in usages:
                    scatter(x_positions, u, c=curr_color, **point_format)
        #get the next color 
        try:
            curr_color_index += 1
        except IndexError:      #fell off end of list -- wrap around
            curr_color_index = 0
            
    legend(lines, labels, loc='best', \
            prop=FontProperties(size='smaller'))
    a = gca()
    a.set_yticks(arange(0,1.01,.1)) #set the yaxis labels
    a.set_xticks(arange(0.5,16.5,1))
    a.xaxis.set_major_formatter(FixedFormatter(axis_names))
    xlim(0,16)
    if graph_name is not None:
        savefig(graph_name)
