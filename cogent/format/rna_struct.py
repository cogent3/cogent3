#!/usr/bin/env python
"""
Writer functions for RNA 2D structures.

    NOTE: Still in beta testing.
"""
from matplotlib import use
use('Agg')  #suppress graphical rendering

from cogent.app.vienna_package import plot_from_seq_and_struct
from cogent.parse.rna_plot import RnaPlotParser
from cogent.parse.record_finder import LabeledRecordFinder
from matplotlib.patches import Circle, Rectangle, Polygon
from matplotlib.text import Text
from pylab import gcf,gca, draw, savefig, clf
from math import ceil

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"



################################################################################
##### Code for postscript RNA 2d Struct ########################################
################################################################################

COLOR_FUNCTION_STRING = \
"""/Colormark { % i color Colormark   draw circle around base i
   % setcolor
   setrgbcolor
   newpath 1 sub coor exch get aload pop
   fsize 2 div 0 360 arc fill stroke
} bind def

/Backcolor { % i color Backcolor   draw circle around base i
   % setcolor
   setrgbcolor
   newpath 1 sub coor exch get aload pop
   fsize 2.6 div 0 360 arc fill stroke
} bind def

/Backgroundbases {
  indices {
    WHITE Backcolor
  } forall
} bind def"""

#Draw structure outline, base pairs, and base background.
INIT_STRUCT_COMMAND = \
"""drawoutline
drawpairs
Backgroundbases
"""

def get_indices_string(sequence):
    """Returns postscript index definition string.
        
        - sequence: string, list, etc representing the sequence of the RNA.
            Calling len(sequence) must give the length of the sequence.
    """
    index_string = '/indices [%s] def'%(' '.join(map(str,\
                                [i+1 for i in range(len(sequence))])))
    return index_string

def get_color_map_string(color_map):
    """Returns postscript color map string given color map.
        
        - color_map: mapping from name to RGB color. {name:[R,G,B]}
            - Name must be a string with no spaces.
            - [R, G, B] must be a list of fractional RGB values from 0 to 1.
                eg {"red": [1.0, 0.0, 0.0]
    """
    color_list = []
    for name, color in color_map.items():
        r,g,b = color
        color_list.append('/%s { %s %s %s } def'%(name,r,g,b))
    return '\n'.join(color_list)

def get_color_commands(indices, colors):
    """Returns postscript coloring string given indices and colors.
    
        - indices: base 1 index of color.  NOT BASE 0.
        - colors: color name corresponding to color name used to generate
            color_map.
        - indices and colors must be lists of the same length.
    """
    color_commands = []
    for index, color in zip(indices,colors):
        color_commands.append('%s %s Colormark'%(str(index),color))
    return ' '.join(color_commands)

def get_circle_commands(indices, color='seqcolor'):
    """Returns postscript circling string given indices.
    
        - indices: base 1 index of color.  NOT BASE 0.
        - color: color name of circle. 
    """
    circle_commands = []
    if indices:
        circle_commands.append(color)
        for index in indices:
            circle_commands.append(str(index)+' cmark')
    return ' '.join(circle_commands)

def get_rnaplot_postscript(sequence, struct):
    """Returns postscript string for seq and struct.
    """
    #Params for RNAplot
    params = {'-t':'0',\
              '--pre':'%PreTextHere'}

    #Get the postscript list
    ps_list = plot_from_seq_and_struct(sequence,\
        struct,params=params).split('\n')
    
    #parse it into prefix and suffix lists
    pre_finder = LabeledRecordFinder(\
        is_label_line=lambda x: x.startswith('%PreTextHere'))
    prefix,suffix = list(pre_finder(ps_list))
    
    #Remove drawoutline and drawpairs commands form suffix
    new_suffix = []
    for s in suffix:
        if not (s.startswith('drawpairs') or s.startswith('drawoutline')):
            new_suffix.append(s)
    
    return '\n'.join(prefix), '\n'.join(new_suffix)
    
def color_on_structure(sequence,struct,color_map,indices=None,colors=None,\
    circle_indices=None):
    """Returns a postscript string colored at indices.
    """
    if indices is None:
        indices = []
    if colors is None:
        colors = []
    
    if len(indices) != len(colors):
        raise ValueError, 'indices and colors must be equal sized lists'
    #Get indices string
    indices_string = get_indices_string(str(sequence))
    
    #Get color map string
    color_map_string = get_color_map_string(color_map)
    
    #Get color commands
    color_commands = get_color_commands(indices,colors)
    
    #Get circle commands
    circle_commands = get_circle_commands(circle_indices)
    
    #get RNAplot postscript
    prefix, suffix = get_rnaplot_postscript(sequence, struct)
    
    return '\n'.join([prefix,COLOR_FUNCTION_STRING,\
            indices_string,INIT_STRUCT_COMMAND,\
            color_map_string,color_commands,circle_commands,suffix])

################################################################################
##### END Code for postscript RNA 2d Struct ####################################
################################################################################

################################################################################
##### Code for matplotlib RNA 2d Struct ########################################
################################################################################

def scale_coords(coords):
    """Returns coordinate list scaled to matplotlib coordintates.
    
        - coords: list of lists, of x,y coordinates:
            [[x1,y1],[x2,y2],[x3,y3],...]
    """
    new_coords = []
    #get min and max x coordinates
    max_x = max([c[0] for c in coords])
    min_x = min([c[0] for c in coords])
    #get min and max y coordinates
    max_y = max([c[1] for c in coords])
    min_y = min([c[1] for c in coords])
    
    #Get scaled max values for x and y
    scaled_max_x = max_x - min_x
    scaled_max_y = max_y - min_y
    
    #max scale value
    max_scale = max(scaled_max_x, scaled_max_y)
    
    scale_x = min_x
    scale_y = min_y
    
    for x,y in coords:
    
        new_coords.append([(x-scale_x)/max_scale,\
            (y-scale_y)/max_scale])
    
    return new_coords

def set_axis_limits(axis, all_coords):
    """Sets axis limits based on all coordinates preventing clipping.
        
        axis: from calling gca()
        all_coords: all coordinates for shapes and labels.
    """
    all_x = [x[0] for x in all_coords]
    all_y = [x[1] for x in all_coords]
    
    offset = abs(all_coords[0][0]-all_coords[1][0]) + \
        abs(all_coords[0][1]-all_coords[1][1])
    
    min_x = min(all_x)
    max_x = max(all_x)
    
    min_y = min(all_y)
    max_y = max(all_y)
    
    axis.set_xlim((-offset,1.0+offset))
    axis.set_ylim((-offset,1.0+offset))
    

def make_circles(coords,facecolors,edgecolors,radius=0.02,alpha=1.0,fill=False):
    """Returns list of Circle objects, given list of coordinates.
    
        - coords: list of [x,y] coordinates, already scaled to matplotlib axes.
        - facecolor: color of circle face.
        - edgecolor: color of circle edge.
        - radius: radius of circle.
    """
    recenter_divide = radius*100.
    recenter = radius/recenter_divide
    circles = []
    for coord, facecolor,edgecolor in zip(coords,facecolors,edgecolors):
        x_coord,y_coord = coord
        curr_circle = Circle([x_coord+recenter,y_coord+recenter],\
            radius=radius,facecolor=facecolor,edgecolor=edgecolor,alpha=alpha,\
            fill=fill)
        circles.append(curr_circle)
    
    return circles

def make_boxes(coords, facecolor='white', edgecolor='black', edge_size=0.03,\
    alpha=1.0,fill=False):
    """Returns list of Rectangle objects, given list of coordinates.
    
        - coords: list of [x,y] coordinates, already scaled to matplotlib axes.
        - facecolor: color of box face.
        - edgecolor: color of box edge.
        - edge_size: length of box edges.
    """
    boxes = []
    recenter_divide = edge_size*200.
    recenter = edge_size/recenter_divide
    for x_coord,y_coord in coords:
        curr_box = Rectangle([x_coord-recenter,y_coord-recenter],\
            edge_size,edge_size,\
            facecolor=facecolor,edgecolor=edgecolor,alpha=alpha,fill=fill)
        boxes.append(curr_box)
    
    return boxes

def make_letters(coords,letters,color='black'):
    """Returns list of Text objects, given list of coordinates and letteres.
    
        - coords: list of [x,y] coordinates, already scaled to matplotlib axes.
        - letters: list of letters to be drawn at given coordinates.  Must be
            same size list as coords.
        - color: color of the letters.
    """
    letter_list = []
    for coord, letter in zip(coords, letters):
        x_coord, y_coord = coord
        curr_letter = Text(x_coord, y_coord, letter)
        letter_list.append(curr_letter)
    
    return letter_list

def make_pairs(coords,pair_indices,offset=.01):
    """Returns list of Polygon objects, given a list of coordinates and indices.
    
        - coords: list of [x,y] coordinates, already scaled to matplotlib axes.
        - pair_indices: indices in the coordinate list that are paired.
    """
    pairs = []
    for first, second in pair_indices:
        fx, fy = coords[first]
        sx, sy = coords[second]
        pairs.append(Polygon([[fx+offset,fy+offset],[sx+offset,sy+offset]],\
            alpha=0.2,linewidth=2))
    return pairs

def make_outline(coords,offset=.01):
    """Returns Polygon object given coords. 
    """
    outline_coords = [[x+offset,y+offset] for x,y in coords]
    outline = Polygon(outline_coords,alpha=0.2,linewidth=2,facecolor='white')
    finish = Polygon([outline_coords[0],outline_coords[-1]],\
            alpha=1,linewidth=3,edgecolor='white',facecolor='white')
    return [outline,finish]

def make_labels(coords):
    """Returns Text objects with 5' and 3' labels.
    """
    #get five prime coordinates
    fp_x = coords[0][0] - (coords[1][0] - coords[0][0])
    fp_y = coords[0][1] - (coords[1][1] - coords[0][1])
    fp_label = Text(fp_x,fp_y,"5'")
    
    #get three prime coordinates
    tp_x = coords[-1][0] - (coords[-2][0] - coords[-1][0])
    tp_y = coords[-1][1] - (coords[-2][1] - coords[-1][1])
    tp_label = Text(tp_x,tp_y,"3'")

    return [fp_label,tp_label], [[fp_x,fp_y],[tp_x,tp_y]]

def draw_structure(sequence,struct,indices=None,colors=None,\
    circle_indices=None, square_indices=None,radial=True):
    """Returns a postscript string colored at indices.
        
        sequence: string of sequence characters.
        struct: string of ViennaStructure for sequence.  Must be valid structure
            same length as sequence.
        indices: list of indices in sequence that will be colored as a solid
            circle.
        colors: list of colors, same length as list of indices.
        circle_indices: list of indices in sequence to draw an empty circle 
            around.
        square_indices: list of indices in sequence to draw an empty square
            around.
        radial: draw structue in radial format (default=True).
    """
    seq_len_scale = int(len(sequence)/50.)
    
    if seq_len_scale < 1:
        circle_scale_size = 0.
        square_scale_size = 0.
    else:
        #Get circle radius.  Proportional to sequence length
        circle_scale_size = (.02/seq_len_scale)/4.0
        #Get edge size.  Proportional to sequence length
        square_scale_size = (.03/seq_len_scale)/4.0
        
    circle_radius = .02 - circle_scale_size
    square_edge_size = .03 - square_scale_size
    
    if indices is None:
        indices = []
    if colors is None:
        colors = []
    if circle_indices is None:
        circle_indices = []
    if square_indices is None:
        square_indices = []
    
    if len(indices) != len(colors):
        raise ValueError, 'indices and colors must be equal sized lists'
    
    if radial:
        params = {'-t':'0'}
    else:
        params = {'-t':'1'}
    #Get the postscript list
    ps_list = plot_from_seq_and_struct(sequence,\
        struct,params=params).split('\n')
    #Parse out seq, coords, and pairs
    seq, coords, pair_list = RnaPlotParser(ps_list)
    coords = scale_coords(coords)
    
    #get letters
    letters = make_letters(coords, list(seq))
    #get pairs
    pairs = make_pairs(coords, pair_list)
    #get outline
    outline = make_outline(coords)
    #get labels
    labels,label_coords = make_labels(coords)


    #get plain circle coords
    circle_coords = [coords[i] for i in circle_indices]
    circle_faces = ['white']*len(circle_coords)
    circle_edges = ['black']*len(circle_coords)
    plain_circles = make_circles(circle_coords,circle_faces, circle_edges,\
        radius=circle_radius)
    
    #get motif circles
    motif_coords = [coords[i] for i in indices]
    motif_circles = make_circles(motif_coords,colors,colors,fill=True,\
        radius=circle_radius)
    
    #Get square coords
    square_coords = [coords[i] for i in square_indices]
    plain_squares = make_boxes(square_coords,edge_size=square_edge_size)
    
    axis = gca()
    axis.set_axis_off()
    
    all_coords = coords + label_coords
    
    set_axis_limits(axis, all_coords)
    
    for l in [letters, pairs, outline, motif_circles, plain_circles, \
        plain_squares,labels]:
        for shape in l:
            axis.add_artist(shape)
    
    fig = gcf()
    
    largest = max([fig.get_figheight(), fig.get_figwidth()])
    fig.set_figheight(largest)
    fig.set_figwidth(largest)

def save_structure(sequence,struct,out_filename,indices=None,colors=None,\
    circle_indices=None, square_indices=None,radial=True,format='png',\
    dpi=75):
    """Saves figure of 2D structure generated by draw_structure().
    
        sequence: string of sequence characters.
        struct: string of ViennaStructure for sequence.  Must be valid structure
            same length as sequence.
        indices: list of indices in sequence that will be colored as a solid
            circle.
        colors: list of colors, same length as list of indices.
        circle_indices: list of indices in sequence to draw an empty circle 
            around.
        square_indices: list of indices in sequence to draw an empty square
            around.
        radial: draw structue in radial format (default=True).
        format: format of structure figure (default=png).  Must be a valid
            matplotlib format.
        dpi: resolution of figure (default=75).

    """
    draw_structure(sequence,struct,indices=indices,colors=colors,\
        circle_indices=circle_indices,\
        square_indices=square_indices,radial=radial)
    savefig(out_filename,format=format,dpi=dpi)
    clf()
        
################################################################################
##### End Code for matplotlib RNA 2d Struct ####################################
################################################################################


    
