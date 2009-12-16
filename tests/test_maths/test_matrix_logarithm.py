#!/usr/bin/env python
"""Unit tests for matrix logarithm."""
from numpy import array
from cogent.util.unit_test import TestCase, main
from cogent.maths.matrix_logarithm import logm, logm_taylor

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class logarithm_tests(TestCase):
    """Tests of top-level matrix logarithm functions."""
    def test_logm(self):
        """logm results should match scipy's"""
        p = array([[ 0.86758487,  0.05575623,  0.0196798 ,  0.0569791 ],
       [ 0.01827347,  0.93312148,  0.02109664,  0.02750842],
       [ 0.04782582,  0.1375742 ,  0.80046869,  0.01413129],
       [ 0.23022035,  0.22306947,  0.06995306,  0.47675713]])

        q = logm(p)
        self.assertFloatEqual(q, \
        array([[-0.15572053,  0.04947485,  0.01918653,  0.08705915],
       [ 0.01405019, -0.07652296,  0.02252941,  0.03994336],
       [ 0.05365208,  0.15569116, -0.22588966,  0.01654642],
       [ 0.35144866,  0.31279003,  0.10478999, -0.76902868]]))
    
    def test_logm_taylor(self):
        """logm_taylor should return same result as logm"""
        q_eig = logm([[ 0.86758487,  0.05575623,  0.0196798 ,  0.0569791 ],
                       [ 0.01827347,  0.93312148,  0.02109664,  0.02750842],
                       [ 0.04782582,  0.1375742 ,  0.80046869,  0.01413129],
                       [ 0.23022035,  0.22306947,  0.06995306,  0.47675713]])
        q_taylor = logm_taylor([[0.86758487, 0.05575623, 0.0196798, 0.0569791],
                  [ 0.01827347,  0.93312148,  0.02109664,  0.02750842],
                  [ 0.04782582,  0.1375742 ,  0.80046869,  0.01413129],
                  [ 0.23022035,  0.22306947,  0.06995306,  0.47675713]])
        self.assertFloatEqual(q_taylor, q_eig)
    

if __name__ == '__main__':
    main()
