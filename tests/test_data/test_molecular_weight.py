#!/usr/bin/env python
"""Tests for molecular weight.
"""
from cogent3.data.molecular_weight import (
    DnaMW,
    ProteinMW,
    RnaMW,
    WeightCalculator,
)
from cogent3.util.unit_test import TestCase, main


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2020, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
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
        self.assertEqual(p(""), 0)
        self.assertEqual(r(""), 0)
        self.assertFloatEqual(p("A"), 89.09)
        self.assertFloatEqual(r("A"), 375.17)
        self.assertFloatEqual(p("AAA"), 231.27)
        self.assertFloatEqual(r("AAA"), 1001.59)
        self.assertFloatEqual(r("AAACCCA"), 2182.37)
        self.assertFloatEqual(
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
