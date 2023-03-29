#!/usr/bin/env python

"""Unit tests for union_dict.
"""
from unittest import TestCase, main

from cogent3.util.warning import deprecate


__author__ = "Richard Morris"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Richard Morris"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class WarningTests(TestCase):
    """Tests of functions in warnings"""

    def test_deprecate(self):
        """test that deprecated decorator functions"""

        def new_function():
            return True

        @deprecate(new_function, "2023.3", "test deprecation")
        def old_function():
            return False

        self.assertTrue(new_function())
        self.assertTrue(old_function())


if __name__ == "__main__":
    main()
