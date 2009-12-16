#!/usr/bin/env python
"""Unit tests of the get_by_cai filter classes."""
from cogent.util.unit_test import TestCase, main
import cogent.maths.stats.cai.get_by_cai

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"


class get_by_cai_tests(TestCase):
    """Tests of top-level functionality.
    
    NOTE: The adaptors are currently tested in an integration test with the
    drawing modules in test_draw/test_matplotlib/test_codon_usage. There are
    not individual unit tests at present, although these should possibly be
    added later.
    """

if __name__ == '__main__':
    main()
