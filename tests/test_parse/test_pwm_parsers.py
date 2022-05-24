from unittest import TestCase, main

from numpy import array
from numpy.testing import assert_allclose, assert_array_equal

from cogent3.core.moltype import get_moltype
from cogent3.parse import cisbp, jaspar


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class TestPwmParsers(TestCase):
    def test_jaspar(self):
        """correctly load jaspar formatted counts matrix"""
        path = "data/sample.jaspar"
        mid, pwm = jaspar.read(path)
        assert mid == ["PSSMid", "HGNCsymbol"], "ID line wrong"
        # note state indices are ordered by moltype
        list(get_moltype("dna"))
        expect = [
            [35, 374, 30, 121, 6, 121, 33],
            [0, 10, 0, 0, 3, 2, 44],
            [352, 3, 354, 268, 360, 222, 155],
            [2, 2, 5, 0, 10, 44, 157],
        ]
        assert_array_equal(pwm.array, array(expect).T)
        self.assertEqual(pwm[0, "A"], 352)
        self.assertEqual(pwm[3, "T"], 121)

    def test_cisbp(self):
        """correctly read a wights matrix"""
        path = "data/M0926_1.02.txt"
        pfm = cisbp.read(path)
        expect = [
            [0.37, 0.55, 0.17, 0.06, 0.64, 0.28, 0.18, 0.26, 0.35],
            [0.24, 0.19, 0.07, 0.09, 0.15, 0.28, 0.08, 0.29, 0.28],
            [0.2, 0.17, 0.53, 0.78, 0.11, 0.16, 0.26, 0.18, 0.19],
            [0.18, 0.09, 0.24, 0.06, 0.1, 0.28, 0.48, 0.27, 0.18],
        ]
        assert_allclose(pfm.array, array(expect).T, atol=1e-2)
        assert_allclose(pfm[0, "A"], 0.199862209150251)
        self.assertEqual(pfm[6, "C"], 0.0787969447816471)


if __name__ == "__main__":
    main()
