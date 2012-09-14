#!/usr/bin/env python

"""Visualisation of phylogenetic compatibility within an alignment.
Jakobsen & Easteal, CABIOS 12(4), 1996
Jakobsen, Wilson & Easteal, Mol. Biol. Evol. 14(5), 1997 
"""

from __future__ import division
import sys
import math
import numpy
import numpy.random
import operator
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.colors

from cogent.draw.linear import Display

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

def order_to_cluster_similar(S, elts=None, start=None):
    """Order so as to keep the most similar parts adjacent to each other
    S is expected to be a square matrix, and the returned list of 
    ordinals is len(S) long."""
    position = {}
    unavailable = set()
    if elts is not None:
        if len(elts) < 2:
            return elts
        elts = set(elts)
        for x in range(len(S)):
            if x != start and x not in elts:
                unavailable.add(x)
    if start is not None:
        position[start] = (False, [start])
    similarity = [numpy.unravel_index(p, S.shape) 
                for p in numpy.argsort(S, axis=None)]
    for (x, y) in similarity[::-1]:
        if x==y or x in unavailable or y in unavailable:
            continue
        (x_end, x_run) = position.get(x, (None, [x]))
        (y_end, y_run) = position.get(y, (None, [y]))
        if x_run is y_run:
            continue
        if x_end is not None:
            if not x_end:
                x_run.reverse()
            unavailable.add(x)
        if y_end is not None:
            if y_end:
                y_run.reverse()
            unavailable.add(y)
        run = x_run + y_run
        position[run[0]] = (False, run)
        position[run[-1]] = (True, run)
    if start is not None:
        if run[-1] == start:
            run.reverse()
        run = run[1:]
    return run

def tied_segments(scores):
    """(start, end) of each run of equal values in scores
    
    >>> tied_segments([1,1,1,2])
    [(0, 3), (3, 4)]
    """
    pos = numpy.flatnonzero(numpy.diff(scores))
    pos = numpy.concatenate(([0],pos+1,[len(scores)]))
    return zip(pos[:-1],pos[1:])

def order_tied_to_cluster_similar(S, scores):
    """Use similarity measure S to make similar elements
    adjacent, but only to break ties in the primary order
    defined by the list of scores"""
    assert S.shape == (len(scores), len(scores))
    new_order = []
    start = None    
    for (a,b) in tied_segments(scores):
        useful = range(a,b)
        if start is not None:
            useful.append(start)
            start = len(useful)-1
        useful = numpy.array(useful)
        S2 = S[useful,:]
        S2 = S2[:,useful]
        sub_order = order_to_cluster_similar(S2, range(b-a), start)
        new_order.extend([useful[i] for i in sub_order])
        start = new_order[-1]
    assert set(new_order) == set(range(len(scores))) 
    return new_order

def bit_encode(x, _bool2num=numpy.array(["0","1"]).take):
    """Convert a boolean array into an integer"""
    return int(_bool2num(x).tostring(), 2)

def bit_decode(x, numseqs):
    """Convert an integer into a boolean array"""
    result = numpy.empty([numseqs], bool)
    bit = 1 << (numseqs-1)
    for i in range(numseqs):
        result[i] = bit & x
        bit >>= 1
    return result


def binary_partitions(alignment):
    """Returns (sites, columns, partitions) 
    sites[informative column number] = alignment position number
    columns[informative column number] = distinct partition number
    partitions[distinct partition number] = (partition, mask) as ints
    """
    sites = []
    columns = []
    partitions = []
    partition_index = {}
    for (site, column) in enumerate(alignment.Positions):
        column = numpy.array(column)
        (A, T, C, G, R, Y, W, S, U) = [
            bit_encode(column == N) for N in "ATCGRYWSU"]
        T |= U
        for split in [([A, G, R], [C, T, Y]), ([A, T, W], [G, C, S])]:
            halves = []
            for char_group in split:
                X = reduce(operator.or_, char_group)
                if not (X & (X - 1)):
                    break  # fewer than 2 bits set in X
                halves.append(X)
            else:
                (X, Z) = sorted(halves)
                partition = (X,X|Z)
                if partition not in partition_index:
                    partition_index[partition] = len(partitions)
                    partitions.append(partition)
                sites.append(site)
                columns.append(partition_index[partition])
                break  # if R/Y split OK no need to consider W/S split. 
    return (sites, columns, partitions)

def min_edges(columns):
    """Given two boolean arrays each representing an informative alignment 
    position, there are 4 possible combinations for each sequence: 
    TT, TF, FT and FF.
    If N of these 4 possibilities are found then there must be at least 
    N-1 tree edges on which mutations occured
    As a special case, the diagonal values are set to 0 rather than, 
    as theory suggests, 1.  This is simply a convenience for later 
    drawing code"""  
    N = len(columns)
    result = numpy.zeros([N, N], int)
    for i in range(0, N-1):
        (a, mask_a) = columns[i]
        for j in range(i+1, N):
            (b, mask_b) = columns[j]
            mask = mask_a & mask_b
            (na, nb) = (~a, ~b)
            combos = [c & mask for c in [a&b, a&nb, na&b, na&nb]]
            combos = [c for c in combos if c]
            result[i,j] = result[j,i] = len(combos) - 1
    return result

def neighbour_similarity_score(matrix):
    left = matrix[:-1]
    right = matrix[1:]
    upper = matrix[:,:-1]
    lower = matrix[:,1:]
    same = (lower == upper).sum() + (left == right).sum()
    neighbours = numpy.product(left.shape)+numpy.product(upper.shape)
    return same / neighbours

def shuffled(matrix):
    assert matrix.shape == (len(matrix), len(matrix)), matrix.shape
    index = numpy.random.permutation(numpy.arange(len(matrix)))
    return matrix[index,:][:,index]
    
def nss_significance(matrix, samples=10000):
    score = neighbour_similarity_score(matrix)
    scores = numpy.empty([samples])
    for i in range(samples):
        s = neighbour_similarity_score(shuffled(matrix))
        scores[i] = s
    scores.sort()
    p = (samples-scores.searchsorted(score)+1) / samples
    return (score, sum(scores)/samples, p)
    
def inter_region_average(a):
    return a.sum()/numpy.product(a.shape)
    
def intra_region_average(a):
    d = numpy.diag(a)    # ignore the diagonal
    return (a.sum()-d.sum())/(numpy.product(a.shape)-len(d))
    
def integer_tick_label(sites):
    def _formatfunc(x, pos, _sites=sites, _n=len(sites)):
        if 0 < x < _n:
            return str(_sites[int(x)])
        else:
            return ""
    return _formatfunc    

def boolean_similarity(matrix):
    # same as numpy.equal.outer(matrix, matrix).trace(axis1=1, axis2=3)
    # but that would use much memory
    true = matrix.T.astype(int)
    false = (~matrix).T.astype(int)
    both_true = numpy.inner(true, true)
    both_false = numpy.inner(false, false)
    return both_true + both_false
    
def partimatrix(alignment, display=False, samples=0, s_limit=0, title="",
        include_incomplete=False, print_stats=True, max_site_labels=50):
    if print_stats:
        print "%s sequences in %s bp alignment" % (
                alignment.getNumSeqs(), len(alignment))
    (sites, columns, partitions) = binary_partitions(alignment)
    if print_stats:
        print "%s unique binary partitions from %s informative sites" % (
                len(partitions), len(sites))
    partpart = min_edges(partitions)      # [partition,partition]
    partimatrix = partpart[columns,:]     # [site, partition]
    sitematrix = partimatrix[:,columns]   # [site, site]
    
    # RETICULATE, JE 1996
    
    compatiblity = sitematrix <= 2
    if print_stats:
        print "Overall compatibility %.6f" % intra_region_average(compatiblity)
        if samples == 0:
            print "Neighbour similarity score = %.6f" % \
                    neighbour_similarity_score(compatiblity)
        else:
            print "Neighbour similarity = %.6f, avg random = %.6f, p < %s" % \
                    nss_significance(compatiblity, samples=samples)
        
    # PARTIMATRIX, JWE 1997
    
    # Remove the incomplete partitions with gaps or other ambiguities
    mask = 2**alignment.getNumSeqs()-1
    complete = [i for (i,(x, xz)) in enumerate(partitions) if xz==mask]
    if not include_incomplete:
        partimatrix = partimatrix[:,complete]
        partitions = [partitions[i] for i in complete]
    # For scoring/ordering purposes, also remove the incomplete sequences
    complete_columns = [i for (i,c) in enumerate(columns) if c in complete]
    scoreable_partimatrix = partimatrix[complete_columns, :]
    
    # Order partitions by increasing conflict score
    conflict = (scoreable_partimatrix > 2).sum(axis=0)
    conflict_order = numpy.argsort(conflict)
    partimatrix = partimatrix[:, conflict_order]
    partitions = [partitions[i] for i in conflict_order]
    scoreable_partimatrix = partimatrix[complete_columns, :]
    support = (scoreable_partimatrix == 0).sum(axis=0)
    consist = (scoreable_partimatrix <= 2).sum(axis=0)
    conflict = (scoreable_partimatrix > 2).sum(axis=0)
    
    # Similarity measure between partitions
    O = boolean_similarity(scoreable_partimatrix <= 2)
    s = 1.0*len(complete_columns)
    O = O.astype(float) / s
    p,q = consist/s, conflict/s
    E = numpy.outer(p,p) + numpy.outer(q,q)
    S = (O-E)/numpy.sqrt(E*(1-E)/s)
    
    # Order partitions for better visual grouping
    if "order_by_conflict":
        order = order_tied_to_cluster_similar(S, conflict)
    else:
        order = order_to_cluster_similar(S)
        half = len(order) // 2
        if sum(conflict[order[:half]]) > sum(conflict[order[half:]]):
            order.reverse()
    
    partimatrix = partimatrix[:, order]
    conflict = conflict[order]
    support = support[order]
    partitions = [partitions[i] for i in order]
    
    if display:
        figwidth = 8.0

        (c_size, p_size) = partimatrix.shape
        s_size = num_seqs = alignment.getNumSeqs()
        
        # Layout (including figure height) chosen to get aspect ratio of
        # 1.0 for the compatibility matrix, and if possible the other
        # matrices.
        
        if s_size > s_limit:
            # too many species to show
            s_size = 0
        else:
            # distort squares to give enough space for species names
            extra = max(1.0, (12/80)/(figwidth/(c_size + p_size)))
            p_size *= numpy.sqrt(extra)
            s_size *= extra
        
        genemap = Display(alignment, recursive=s_size>0, 
                colour_sequences=False, draw_bases=False)
        annot_width = max(genemap.height / 80, 0.1)
        figwidth = max(figwidth, figwidth/2 + annot_width)
        
        bar_height = 0.5
        link_width = 0.3
        x_margin = 0.60
        y_margin = 0.35
        xpad = 0.05
        ypad = 0.2
        (x, y) = (c_size + p_size, c_size + s_size)
        x_scale = y_scale = (figwidth-2*x_margin-xpad-link_width-annot_width)/x
        figheight = y_scale * y + 2*y_margin + 2*ypad + bar_height
        x_scale /= figwidth
        y_scale /= figheight
        x_margin /= figwidth
        y_margin /= figheight
        xpad /= figwidth
        ypad /= figheight
        bar_height /= figheight
        link_width /= figwidth
        annot_width /= figwidth
        (c_width, c_height) = (c_size*x_scale, c_size*y_scale)
        (p_width, s_height) = (p_size*x_scale, s_size*y_scale)
        vert = (x_margin + xpad + c_width)
        top = (y_margin + c_height + ypad)
        fig = plt.figure(figsize=(figwidth,figheight))
        kw = dict(axisbg=fig.get_facecolor())
        axC = fig.add_axes([x_margin, y_margin, c_width, c_height], **kw)
        axP = fig.add_axes([vert, y_margin, p_width, c_height], 
                sharey=axC, **kw)
        axS = fig.add_axes([vert, top, p_width, s_height or .001], 
                sharex=axP, **kw)
        axB = fig.add_axes([vert, top+ypad+s_height, p_width, bar_height], 
                sharex=axP, **kw)
        axZ = fig.add_axes([vert+p_width, y_margin, link_width, c_height], 
            frameon=False)
            
        axA = genemap.asAxes(
            fig, [vert+p_width+link_width, y_margin, annot_width, c_height], 
            vertical=True, labeled=True)
            
        axP.yaxis.set_visible(False)
        #for ax in [axC, axP, axS]:
            #ax.set_aspect(adjustable='box', aspect='equal')
        
        fig.text(x_margin+c_width/2, .995, title, ha='center', va='top')
        
        if not s_size:
            axS.set_visible(False)
        # No ticks for these non-float dimensions
        for axes in [axB, axC, axS, axP]:
            for axis in [axes.xaxis, axes.yaxis]:
                for tick in axis.get_major_ticks():
                    tick.gridOn = False
                    tick.tick1On = False
                    tick.tick2On = False
                    tick.label1.set_size(8)
                    tick.label2.set_size(8)
                    if axis is axes.xaxis:
                        tick.label1.set_rotation('vertical')

        # Partition dimension
        for axis in [axS.xaxis, axP.xaxis, axB.xaxis, axB.yaxis]:
            axis.set_major_formatter(matplotlib.ticker.NullFormatter())
            axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        
        # Site dimension
        if c_size > max_site_labels:
            for axis in [axC.yaxis, axC.xaxis]:
                axis.set_visible(False)
        else:
            isl = integer_tick_label(sites)
            for axis in [axC.yaxis, axC.xaxis]:            
                axis.set_minor_locator(matplotlib.ticker.IndexLocator(1,0))
                axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                axis.set_major_locator(matplotlib.ticker.IndexLocator(1,0.5))
                axis.set_major_formatter(matplotlib.ticker.FuncFormatter(isl))
        
        # Species dimension
        if s_size:
            seq_names = [name.split('  ')[0] 
                    for name in alignment.getSeqNames()]
            axS.yaxis.set_minor_locator(matplotlib.ticker.IndexLocator(1,0))
            axS.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            axS.yaxis.set_major_locator(matplotlib.ticker.IndexLocator(1,0.5))
            axS.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(seq_names))
            #axS.yaxis.grid(False) #, 'minor')
    
        # Display the main matrices: compatibility and partimatrix
        axC.pcolorfast(compatiblity, cmap=plt.cm.gray)
        partishow = partimatrix <= 2
        axP.pcolorfast(partishow, cmap=plt.cm.gray)
        axP.set_autoscale_on(False)
        axC.plot([0,c_size], [0, c_size], color='lightgreen')
        (sx, sy) = numpy.nonzero(partimatrix.T==0)
        axP.scatter(sx+0.5, sy+0.5, color='lightgreen', marker='^',
            s=15)
        
        # Make [partition, sequence] matrix
        # Not a good idea with too many sequences
        if s_size:
            partseq1 = numpy.empty([len(partitions), num_seqs], bool)
            partseq2 = numpy.empty([len(partitions), num_seqs], bool)
            for (i, (x, xz)) in enumerate(partitions):
                partseq1[i] = bit_decode(x, num_seqs)
                partseq2[i] = bit_decode(xz^x, num_seqs)
            
            # Order sequqnces so as to place similar sequences adjacent
            O = boolean_similarity(partseq1)
            order = order_to_cluster_similar(O)
            partseq1 = partseq1[:,order]
            partseq2 = partseq2[:,order]
            seq_names = [seq_names[i] for i in order]
            axS.set_ylim(0, len(seq_names))
            axS.set_autoscale_on(False)
    
            for (halfpart,color) in [(partseq1, 'red'),(partseq2, 'blue')]:
                (sx, sy) = numpy.nonzero(halfpart)
                axS.scatter(sx+0.5, sy+0.5, color=color, marker='o')
            axS.grid(False)
            #axS.yaxis.tick_right()
            #axS.yaxis.set_label_position('right')
        
        # Bar chart of partition support and conflict scores
        #axB.set_autoscalex_on(False)
        if conflict.sum():
            axB.bar(numpy.arange(len(partitions)), -conflict/conflict.sum(), 
                1.0, color='black', align='edge')
        if support.sum():
            axB.bar(numpy.arange(len(partitions)), +support/support.sum(), 
                1.0, color='lightgreen', align='edge')
        axB.set_xlim(0.0, len(partitions))
        
        # Alignment features
        axA.set_ylim(0, len(alignment))
        axA.set_autoscale_on(False)
        axA.yaxis.set_major_formatter(
                matplotlib.ticker.FuncFormatter(lambda y,pos:str(int(y))))
        axA.yaxis.tick_right()
        axA.yaxis.set_label_position('right')
        axA.xaxis.tick_top()
        axA.xaxis.set_label_position('top')
        #axA.xaxis.set_visible(False)
        
        # "Zoom lines" linking informative-site coords to alignment coords 
        from matplotlib.patches import PathPatch
        from matplotlib.path import Path
        axZ.set_xlim(0.0,1.0)
        axZ.set_xticks([])
        axZ.set_ylim(0, len(alignment))
        axZ.set_yticks([])
        zoom = len(alignment) / len(sites)
        vertices = []
        for (i,p) in enumerate(sites):
            vertices.extend([(.1, (i+0.5)*zoom), (.9,p+0.5)])
            axA.axhspan(p, p+1, facecolor='green', edgecolor='green', alpha=0.3)
        ops = [Path.MOVETO, Path.LINETO] * (len(vertices)//2)
        path = Path(vertices, ops)
        axZ.add_patch(PathPatch(path, fill=False, linewidth=0.25))
        
        # interactive navigation messes up axZ.  Could use callbacks but
        # probably not worth the extra complexity.
        for ax in [axC, axP, axS, axB, axZ, axA]:
            ax.set_navigate(False)

        return fig

    
if __name__ == '__main__':
    from cogent import LoadSeqs, DNA
    import sys, optparse, os.path
    parser = optparse.OptionParser("usage: %prog [options] alignment")
    parser.add_option("-p", "--print", action="store_true", 
            default=True, dest="print_stats", 
            help="print neighbour similarity score etc.")
    parser.add_option("-d", "--display", action="store_true", 
            default=False, dest="display", 
            help="show matrices via matplotlib")
    parser.add_option("-i", "--incomplete", action="store_true", 
            default=False, dest="include_incomplete", 
            help="include partitions containing ambiguities")
    parser.add_option("-t", "--taxalimit", 
                dest="s_limit", default=20, type="int", 
                help="maximum number of species that can be displayed")
    parser.add_option("-s", "--samples",
                dest="samples", default=10000, type="int",
                help="samples for significance test")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit(1)        
    alignment = LoadSeqs(args[0], moltype=DNA)
    kw = vars(options)
    kw['title'] = os.path.splitext(os.path.basename(args[0]))[0]
    fig = partimatrix(alignment, **kw)
    if fig:
        plt.show()

