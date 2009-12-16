#/usr/bin/env python
from cogent.draw.arrow_rates import make_arrow_plot, sample_data
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

class arrow_rates_tests(TestCase):
    """Tests of top-level function, primarily checking that it writes the file.
    
    WARNING: must visually inspect output to check correctness!
    """
    def test_make_arrow_plot(self):
        """arrow_plot should write correct file and not raise exception"""
        make_arrow_plot(sample_data, graph_name='arrows.png')
        #comment out line below to see the result
        #remove('arrows.png')

if __name__ == '__main__':
    main()
    
