#/usr/bin/env python
from cogent.draw.multivariate_plot import (plot_ordination,
        map_colors, text_points, scatter_points, plot_points, arrows)
from cogent.util.unit_test import TestCase, main

import os, pylab
from tempfile import mkstemp
from numpy import asarray, c_
from pdb import set_trace

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

######
# util class
class TestCasePlot(TestCase):
    Debug = True

    def fig(self, fname=None):
        if self.Debug:
            pylab.show()
        else:
            if not fname:
                fd, fname = mkstemp(prefix='PlotTest_', suffix='.png')
            pylab.savefig(fname)
            pylab.clf()
            os.remove(fname)
    
    def p(self, obj):
        if self.Debug:
            print obj

class TestI(TestCasePlot):
    def setUp(self):
        self.points = [(0, 0), (1,1), (1,2), (2, 2)]
        #adjust canvas
        #pylab.xlim(-0.5, 2.5)
        #pylab.ylim(-0.5, 2.5)

class FunctionsTests(TestI):
    def test_map_colors(self):
        #default
        self.assertEqual(map_colors(range(3)),
                ['#000080', '#7dff7a', '#800000'])
        #alternative cmap
        self.assertEqual(map_colors(range(3), cmap='hot'),
                ['#0b0000', '#ff5c00', '#ffffff'])
        #return tuples
        self.assertFloatEqual(map_colors(range(3), mode='tuples'),
                [(0.0, 0.0, 0.5), (0.490196, 1.0, 0.477546), (0.5, 0.0, 0.0)])

    def test_text_points(self):
        #same text for all the points
        text_points(self.points, 'X')
        self.fig()

    def test_text_points_diff_texts(self):
        #diff text for all the points
        text_points(self.points, ['A', 'B', 'C', 'A'])
        self.fig()

    def test_plot_points(self):
        plot_points(self.points, label='X')
        self.fig()

class scatter_points_tests(TestI):
    def test_basic(self):
        scatter_points(self.points, label='X')
        self.fig()

    def test_color_list(self):
        scatter_points(self.points, c=['k', 'red', '#00FF00', (0, 0, 1)],
                s=[100, 200, 300, 400], marker=['o', 's', 'd', 'h'])
        pylab.xlim(-0.5, 2.5)
        pylab.ylim(-0.5, 2.5)
        self.fig()

    def test_color_shades(self):
        #self.Debug = True
        scatter_points(self.points, c=[1, 2, 3, 4],
                s=[100, 200, 300, 400], marker=['o', 's', 'd', 'h'])
        pylab.xlim(-0.5, 2.5)
        pylab.ylim(-0.5, 2.5)
        self.fig()

class arrows_tests(TestI):
    def test_arrows(self):
        points = self.points
        arrows([[0,0]]*len(points), points)
        self.fig()

class plot_ordination_tests(TestCasePlot):
    def setUp(self):
        self.points = points_3d
        self.keys = ['eigvals', 'samples', 'species', 'centroids', 'biplot']
        self.values = [
                [0.3, 0.1, 0, -0.5], #eigvals
                self.points, #samples
                self.points + 0.2, #species
                [(0, 1, 0), (2, -2, 0)], #centroids
                [(0.5, 0.5, 1), (-1, 1.5, 1)], #biplot
                ]

    def test_basic(self):
        res = dict(zip(self.keys[:2], self.values))
        plot_ordination(res)
        self.fig()

    def test_species(self):
        res = dict(zip(self.keys[:3], self.values))
        plot_ordination(res)
        self.fig()

    def test_centroids(self):
        res = dict(zip(self.keys[:4], self.values))
        plot_ordination(res)
        self.fig()

    def test_biplot(self):
        res = dict(zip(self.keys, self.values))
        plot_ordination(res,
                species_kw={'label': 'sp'}, biplot_kw={'label':['b1', 'b2']},
                samples_kw={'label': 'sa'})
        self.fig()

    def test_choices(self):
        res = dict(zip(self.keys, self.values))
        plot_ordination(res, choices=[2,3])
        self.fig()

    def test_axis_names(self):
        #self.Debug = True
        res = dict(zip(self.keys[:2], self.values))
        plot_ordination(res, axis_names=['CCA1', 'CA1', 'CA2'],
                constrained_names='CCA')
        self.fig()



#####
# test data
points = asarray([(0,0), (-1, 1), (1, -2), (2, 3)], float)
points_3d = c_[points, [[3]]*len(points)]
if __name__ == '__main__':
    main()
