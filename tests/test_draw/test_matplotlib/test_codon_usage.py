#/usr/bin/env python
"""Tests of the codon usage graphs.

Note: currently, this must be invoked manually because all the output is
graphical. Expects to be invoked from the test directory.
"""
from cogent.maths.stats.cai.adaptor import data_from_file, adapt_p12, \
    adapt_p123gc, adapt_cai_p3, adapt_cai_histogram, adapt_fingerprint, \
    adapt_pr2_bias, file_to_codon_list, \
    bin_by_p3
from cogent.draw.codon_usage import plot_cai_p3_scatter, \
    plot_p12_p3, plot_p123_gc, plot_fingerprint, plot_pr2_bias, \
    plot_p12_p3_contour, \
    plot_p12_p3_contourlines, aa_labels
from cogent.draw.util import as_species, \
    plot_scatter_with_histograms, plot_histograms, \
    plot_filled_contour, format_contour_array
from pylab import gca, clf
from numpy import transpose
from sys import argv
from os import getcwd
from os.path import sep, join
test_path = getcwd().split(sep)
index = test_path.index('tests')
fields = test_path[:index+1] + ["data"]
test_path = sep + join(*fields)
test_file_name = join(test_path, 'Homo_sapiens_codon_usage.pri')

__author__ = "Stephanie Wilson"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Stephanie Wilson"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def insert_before_extension(name, item):
    """Inserts item before extension in name."""
    last_dot = name.rfind('.')
    return name[:last_dot+1]+str(item)+'.'+ name[last_dot+1:]
    
#TESTS USING ADAPTORS
def make_generic_adaptor_test(adaptor_f, plot_f, default_outfilename):
    """Makes adaptor test for generic graphs."""
    def adaptor_test(codons, infilename, name, outfilename=default_outfilename):
        print "=> outfile:", outfilename
        print "from file:", infilename
        graph_data = adaptor_f(codons)
        plot_f(graph_data, num_genes=len(codons), graph_name=outfilename,\
            title=name)
        gca().clear()
        clf()
    return adaptor_test

def make_contour_adaptor_test(adaptor_f, plot_f, default_outfilename):
    """Makes adaptor test for contour graphs."""
    def adaptor_test(codons, infilename, name, outfilename=default_outfilename):
        print "=> outfile:", outfilename
        print "from file:", infilename
        xy_data = adaptor_f(codons)
        x, y, data = format_contour_array(xy_data)
        plot_f(x, y, data, xy_data, num_genes=len(codons), \
            graph_name=outfilename, title=name)
        gca().clear()
        clf()
    return adaptor_test

def make_pr2bias_adaptor_test(adaptor_f, plot_f, default_outfilename):
    """Makes adaptor test for pr2bias graphs."""
    def adaptor_test(codons, infilename, name, outfilename=default_outfilename):
        print "=> base outfile:", outfilename
        print "from file:", infilename
        for aa, triplet in aa_labels.items():
            triplet = triplet.replace('T','U')
            graph_data = adaptor_f(codons, block=triplet[:2])
            curr_outfilename = insert_before_extension(outfilename, triplet)
            plot_f(graph_data, num_genes=len(codons), \
                graph_name=curr_outfilename, title=aa)
            gca().clear()
            clf()
    return adaptor_test

def make_gc_gradient_adaptor_test(adaptor_f, plot_f, default_outfilename, \
    one_graph=False, one_series=False):
    """Makes adaptor test for replicated graphs over e.g. a GC gradient"""
    def adaptor_test(codons, infilename, name, outfilename=default_outfilename):
        min_gene_threshold=10   #suppress bins with few genes
        print "=>base outfile:", outfilename
        print "from file:", infilename
        print "one graph:", one_graph
        print "one series:", one_series
        gc_bins = bin_by_p3(codons)
        if one_series:  #assume we want to adapt the list of codon usages
            data = adaptor_f(gc_bins)
            plot_f(data, num_genes=len(codons), graph_name=outfilename,\
                title=name)
        else:
            data = []
            for b in gc_bins:
                try:
                    data.append(adaptor_f(b))
                except:
                    data.append([])
            if one_graph:   #assume downstream f copes with list of data
                total_genes = len(codons)
                alpha = [float(len(b))/total_genes for b in gc_bins]
                plot_f(data, num_genes=total_genes, graph_name=outfilename, \
                    alpha=alpha, title=name, multiple=True)
            else:   #multiple graphs: feed one at a time
                for i, c in enumerate(gc_bins):
                    curr_p3 = i*0.1
                    if len(c) > min_gene_threshold:
                        curr_outfilename = insert_before_extension(\
                            outfilename, curr_p3)
                        graph_data = adaptor_f(c)
                        plot_f(graph_data, num_genes=len(c), \
                            graph_name=curr_outfilename, \
                            title=name+' '+str(curr_p3)+'-'+str(curr_p3+0.1))
                    gca().clear()
                    clf()
        gca().clear()
        clf()
    return adaptor_test
    
mgat = make_generic_adaptor_test    #save typing in what follows
mcat = make_contour_adaptor_test    #ditto
mpat = make_pr2bias_adaptor_test
mggat = make_gc_gradient_adaptor_test

#scatterplot adaptors
test_p12_p3_adaptor = mgat(adapt_p12, plot_p12_p3, 'test_p12_p3_A.png')
test_p123_gc_adaptor = mgat(adapt_p123gc, plot_p123_gc, \
    'test_p123_gc_A.png')

def plot_p12_p3_from_gc(*args, **kwargs):
    return plot_p123_gc(use_p3_as_x=True, graph_shape='sqr',*args, **kwargs)
test_p12_p3gc_adaptor = mgat(adapt_p123gc, plot_p12_p3_from_gc, \
    'test_p12_from_gc_A.png')
test_cai_p3_adaptor = mgat(adapt_cai_p3, plot_cai_p3_scatter, \
    'test_p3_cai_A.png')
def adapt_cai_p3_twoseries(*args, **kwargs):
    return adapt_cai_p3(both_series=True)
def adapt_cai_p3_twoseries(*args, **kwargs): \
    return adapt_cai_p3(both_series=True, *args, **kwargs)
test_cai_p3_twoseries_adaptor = mgat(adapt_cai_p3_twoseries, \
    plot_cai_p3_scatter, 'test_p3_cai_twoseries_A.png')
def scat_hist_cai_p3(data, *args, **kwargs):
    return plot_scatter_with_histograms(data, x_label='$P_3$', y_label='CAI',\
        *args, **kwargs)
test_cai_p3_twoseries_adaptor_hist = mgat(adapt_cai_p3_twoseries, \
    scat_hist_cai_p3, 'test_p3_cai_twoseries_hist_A.png')

#hist adaptors
def cai_histogram(data, *args, **kwargs):
    return plot_histograms(data, x_label='CAI', \
        series_names=['others', 'ribosomal'], show_legend=True, \
        colors=['white','red'], linecolors=['black','red'], alpha=0.7, \
        *args, **kwargs)
test_cai_hist_adaptor = mgat(adapt_cai_histogram, cai_histogram,\
    'test_cai_hist.png')

#fingerprint adaptors
test_fingerprint_adaptor = mgat(adapt_fingerprint, plot_fingerprint, \
    'test_fingerprint_A.png')
test_fingerprint_gradient_adaptor = mggat(adapt_fingerprint, plot_fingerprint,\
    'test_fingerprint_gradient_A.png')
test_fingerprint_gradient_adaptor_one_graph = mggat(adapt_fingerprint, \
    plot_fingerprint, 'test_fingerprint_gradient_onegraph_A.png', one_graph=True)

#pr2 bias adaptors
test_pr2bias_adaptor = mpat(adapt_pr2_bias, plot_pr2_bias, \
    'test_pr2_bias_A.png')

#contour adaptors
test_p12_p3_contour_adaptor = mcat(adapt_p12, plot_p12_p3_contour, \
    'test_p12_p3_contour_A.png')
test_p12_p3_contourlines_adaptor = mcat(adapt_p12, plot_p12_p3_contourlines, \
    'test_p12_p3_contourlines_A.png')

scatter_adaptor = [test_p12_p3_adaptor, test_p123_gc_adaptor, \
    test_p12_p3gc_adaptor, test_cai_p3_adaptor, test_cai_p3_twoseries_adaptor,\
    test_cai_p3_twoseries_adaptor_hist]
hist_adaptor=[test_cai_hist_adaptor]
fingerprint_adaptor = [test_fingerprint_adaptor,
    test_fingerprint_gradient_adaptor, \
    test_fingerprint_gradient_adaptor_one_graph]
pr2bias_adaptor = [test_pr2bias_adaptor]
contour_adaptor = [test_p12_p3_contour_adaptor, \
    test_p12_p3_contourlines_adaptor]

all_adaptor = scatter_adaptor + fingerprint_adaptor + pr2bias_adaptor \
    + hist_adaptor + contour_adaptor

#take in pre-constructed codon usage objects, output requested graphs
def codons_to_graph(codons, as_file, species, which_tests):
    """Function for directly passing in codons instead of 
    reading them in from a file"""
    
    adaptor_tests= all_adaptor
    codon_data_fname=as_file
    print "Running adaptor tests..."
    codon_data = codons
    if which_tests:
        for i in which_tests:
            print "doing test %s" % i
            a = adaptor_tests[i]
            a(codon_data, codon_data_fname, species)
    else:
        for i, a in enumerate(adaptor_tests):
            print "doing test %s" % i
            a(codon_data,codon_data_fname, species)
            
if __name__ == '__main__':
    """Tests if the graphs will all compile
    and outputs the graph from the current
    version of code:
    test_fingerprint.png
    """
    adaptor_tests = all_adaptor
    if len(argv) > 1:
        codon_data_fname = argv[1]
    else:
        codon_data_fname = test_file_name

    if len(argv) > 2:
        which_tests = map(int, argv[2].split(','))
    else:
        which_tests = None
        
    print "Running adaptor tests..."
    if codon_data_fname.endswith('.nuc'):   #assume FASTA from KEGG
        codon_data = kegg_fasta_to_codon_list(open(codon_data_fname))
    else:
        codon_data = file_to_codon_list(codon_data_fname)

    if which_tests:
        for i in which_tests:
            print "doing test %s" % i
            a = adaptor_tests[i]
            a(codon_data, codon_data_fname, as_species(codon_data_fname))
    else:
        for i, a in enumerate(adaptor_tests):
            print "doing test %s" % i
            a(codon_data, codon_data_fname, as_species(codon_data_fname))
