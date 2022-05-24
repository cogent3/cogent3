#!/usr/bin/env python
"""Tests for molecular weight.
"""
from unittest import TestCase, main

from cogent3.data.molecular_weight import ProteinMW, RnaMW


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

from numpy.testing import assert_allclose


class WeightCalculatorTests(TestCase):
    """Tests for WeightCalculator, which should calculate molecular weights."""

    def test_call(self):
        """WeightCalculator should return correct molecular weight"""
        r = RnaMW
        p = ProteinMW
        self.assertEqual(p(""), 0)
        self.assertEqual(r(""), 0)
        assert_allclose(p("A"), 89.09)
        assert_allclose(r("A"), 375.17)
        assert_allclose(p("AAA"), 231.27)
        assert_allclose(r("AAA"), 1001.59)
        assert_allclose(r("AAACCCA"), 2182.37)
        assert_allclose(
            p(
                "MVQQAESLEAESNLPREALDTEEGEFMACSPVALDESDPDWCKTASGHIKRPMNAFMVWSKIERRKIMEQSPDMHNAEISKRLGKR\
                                 WKMLKDSEKIPFIREAERLRLKHMADYPDYKYRPRKKPKMDPSAKPSASQSPEKSAAGGGGGSAGGGAGGAKTSKGSSKKCGKLKA\
                                 PAAAGAKAGAGKAAQSGDYGGAGDDYVLGSLRVSGSGGGGAGKTVKCVFLDEDDDDDDDDDELQLQIKQEPDEEDEEPPHQQLLQP\
                                 PGQQPSQLLRRYNVAKVPASPTLSSSAESPEGASLYDEVRAGATSGAGGGSRLYYSFKNITKQHPPPLAQPALSPASSRSVSTSSS\
                                 SSSGSSSGSSGEDADDLMFDLSLNFSQSAHSASEQQLGGGAAAGNLSLSLVDKDLDSFSEGSLGSHFEFPDYCTPELSEMIAGDWL\
                                 EANFSDLVFTY"
            ),
            46685.97,
        )


# run if called from command-line
if __name__ == "__main__":
    main()
