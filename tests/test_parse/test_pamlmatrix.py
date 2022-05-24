#!/usr/bin/env python
from io import StringIO
from unittest import TestCase, main

from cogent3.evolve.models import DSO78_freqs, DSO78_matrix
from cogent3.parse.paml_matrix import PamlMatrixParser


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"

from numpy.testing import assert_equal


data = """
       27									    
       98  32									    
      120   0 905								    
       36  23   0   0								    
       89 246 103 134   0							    
      198   1 148 1153  0 716							    
      240   9 139 125  11  28  81						    
       23 240 535  86  28 606  43  10						    
       65  64  77  24  44  18  61   0   7					    
       41  15  34   0   0  73  11   7  44 257					    
       26 464 318  71   0 153  83  27  26  46  18				    
       72  90   1   0   0 114  30  17   0 336 527 243				    
       18  14  14   0   0   0   0  15  48 196 157   0  92			    
      250 103  42  13  19 153  51  34  94  12  32  33  17  11			    
      409 154 495  95 161  56  79 234  35  24  17  96  62  46 245		    
      371  26 229  66  16  53  34  30  22 192  33 136 104  13  78 550		    
        0 201  23   0   0   0   0   0  27   0  46   0   0  76   0  75   0	    
       24   8  95   0  96   0  22   0 127  37  28  13   0 698   0  34  42  61	    
      208  24  15  18  49  35  37  54  44 889 175  10 258  12  48  30 157   0  28 

      0.087127 0.040904 0.040432 0.046872 0.033474 0.038255 0.049530
      0.088612 0.033618 0.036886 0.085357 0.080482 0.014753 0.039772
      0.050680 0.069577 0.058542 0.010494 0.029916 0.064718

      Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val

      S_ij = S_ji and PI_i for the Dayhoff model, with the rate Q_ij=S_ij*PI_j
      The rest of the file is not used.
      Prepared by Z. Yang, March 1995.


      See the following reference for notation used here:

      Yang, Z., R. Nielsen and M. Hasegawa. 1998. Models of amino acid substitution and
      applications to mitochondrial protein evolution. Mol. Biol. Evol. 15:1600-1611.
"""


class TestParsePamlMatrix(TestCase):
    def test_parse(self):
        matrix, freqs = PamlMatrixParser(StringIO(data))
        assert_equal(DSO78_matrix, matrix)
        assert_equal(DSO78_freqs, freqs)


if __name__ == "__main__":
    main()
