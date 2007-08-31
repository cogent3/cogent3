#!/usr/bin/env python
"""Tests of data retrieval from PDB."""
from cogent.util.unit_test import TestCase, main
from cogent.db.pdb import Pdb

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class PdbTests(TestCase):
    """Tests of the Pdb class."""
    def test_simple_get(self):
        """Simple access of an item should work."""
        rec = Pdb()
        result = rec['1RMN'].read()
        assert result.startswith('HEADER')
        assert result.endswith('1RMN1101\n')
        assert 'HAMMERHEAD RIBOZYME' in result

if __name__ == '__main__':
    main()
