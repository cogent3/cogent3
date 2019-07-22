from unittest import TestCase, main

from numpy import array

from cogent3.parse import cisbp, jaspar


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.07.10a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class TestPwmParsers(TestCase):
    def test_jaspar(self):
        """correctly load jaspar formatted counts matrix"""
        path = "data/sample.jaspar"
        mid, bases, mat = jaspar.read(path)
        assert mid == ["PSSMid", "HGNCsymbol"], "ID line wrong"
        assert bases == list("ACGT"), "base order wrong"
        assert mat == [
            [352, 3, 354, 268, 360, 222, 155],
            [0, 10, 0, 0, 3, 2, 44],
            [2, 2, 5, 0, 10, 44, 157],
            [35, 374, 30, 121, 6, 121, 33],
        ], "Matrix counts wrong"

    def test_cisbp(self):
        """correctly read a wights matrix"""
        path = "data/M0926_1.02.txt"
        bases, mat = cisbp.read(path)
        assert bases == list("ACGT"), "base order wrong"
        mat = array(mat)
        mat = mat.round(2).tolist()
        expect = [
            [0.2, 0.17, 0.53, 0.78, 0.11, 0.16, 0.26, 0.18, 0.19],
            [0.24, 0.19, 0.07, 0.09, 0.15, 0.28, 0.08, 0.29, 0.28],
            [0.18, 0.09, 0.24, 0.06, 0.1, 0.28, 0.48, 0.27, 0.18],
            [0.37, 0.55, 0.17, 0.06, 0.64, 0.28, 0.18, 0.26, 0.35],
        ]
        assert mat == expect, "Matrix weights wrong"


if __name__ == "__main__":
    main()
