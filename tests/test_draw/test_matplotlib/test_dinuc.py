#/usr/bin/env python
from cogent.draw.dinuc import dinuc_plot
from numpy import array, clip
from numpy.random import random
from cogent.util.unit_test import TestCase, main
from os import remove

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def add_random(a):
    """Adds a small random component to gene a, while maintaining the sum."""
    r = (random(a.shape)-.5)/10
    return clip(a + r, 0, 1)
    
class dinuc_tests(TestCase):
    """Tests of top-level function, primarily checking that it writes the file.
    
    WARNING: must visually inspect output to check correctness!
    """
    def test_dinuc(self):
        """dinuc_plot should write correct file and not raise exception"""

        spec_a_ave = array([0.25, 0.20, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
                               0.45, 0.40, 0.55, 0.60, 0.65, 0.70, 0.60, 0.70])
        spec_b_ave = array([0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45,
                               0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.20, 0.15])
        spec_c_ave = array([0.40, 0.50, 0.55, 0.50, 0.45, 0.50, 0.45, 0.50,
                               0.45, 0.55, 0.45, 0.45, 0.45, 0.55, 0.50, 0.55])
        
        spec_a_ht_gene = array([0.41, 0.52, 0.53, 0.53, 0.45, 0.51, 0.46, 0.52,
                               0.43, 0.55, 0.45, 0.45, 0.45, 0.55, 0.53, 0.55])
        spec_b_ht_gene = array([0.26, 0.22, 0.17, 0.19, 0.23, 0.32, 0.36, 0.41,
                               0.44, 0.41, 0.53, 0.62, 0.63, 0.71, 0.64, 0.72])

        spec_b_rb_gene = array([0.82, 0.74, 0.73, 0.63, 0.62, 0.54, 0.51, 0.43,
                               0.41, 0.34, 0.32, 0.24, 0.23, 0.17, 0.28, 0.19])
        spec_c_rb_gene = array([0.43, 0.54, 0.56, 0.51, 0.44, 0.51, 0.42, 0.53,
                                0.44, 0.53, 0.43, 0.47, 0.46, 0.57, 0.53, 0.52])

        
        a_data = {'hgt':[spec_a_ht_gene], \
            None: map(add_random, [spec_a_ave.copy() for i in range(10)])}

        b_data = {'hgt':[spec_b_ht_gene], 'ribosomal':[spec_b_rb_gene], \
            None:[spec_b_ave]}

        c_data = {'ribosomal':[spec_b_rb_gene], None:[spec_c_ave]}
        
        data = {'Species A': a_data, 'Species B': b_data, 'Species C':c_data}

        dinuc_plot(data, avg_formats={'markersize':5}, \
            point_formats={'s':2, 'alpha':0.2}, graph_name='test.png')
        #note:comment out the next line to see the test file
        #remove('test.png')

if __name__ == '__main__':
    main()
    
