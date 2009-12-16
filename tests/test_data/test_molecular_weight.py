#!/usr/bin/env python
"""Tests for molecular weight.
"""
from cogent.util.unit_test import TestCase, main
from cogent.data.molecular_weight import WeightCalculator, DnaMW, RnaMW, \
    ProteinMW

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class WeightCalculatorTests(TestCase):
    """Tests for WeightCalculator, which should calculate molecular weights.
    """
    def test_call(self):
        """WeightCalculator should return correct molecular weight"""
        r = RnaMW
        p = ProteinMW
        self.assertEqual(p(''), 0)
        self.assertEqual(r(''), 0)
        self.assertFloatEqual(p('A'), 107.09)
        self.assertFloatEqual(r('A'), 375.17)
        self.assertFloatEqual(p('AAA'), 285.27)
        self.assertFloatEqual(r('AAA'), 1001.59)
        self.assertFloatEqual(r('AAACCCA'), 2182.37)
     
#run if called from command-line
if __name__ == "__main__":
    main()
