#!/usr/bin/env python

"""Unit tests for table.
"""
from cogent3.util.table import Table
from cogent3.util.unit_test import TestCase, main


__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2016, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La"]
__license__ = "GPL"
__version__ = "3.0a2"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class TableTests(TestCase):
    """Tests of individual functions in table"""

    # test table 1
    t1_header = ['id', 'foo', 'bar']
    t1_rows = [
        [1, 'abc', 11],
        [2, 'bca', 22],
        [3, 'cab', 33],
        [4, 'abc', 44],
        [5, 'bca', 55],
    ]
    t1 = Table(header=t1_header, rows=t1_rows)

    # test table 2
    t2_header = ['chrom', 'stableid', 'length']
    t2_rows = [
        ['X', 'ENSG00000005893', 1353],
        ['A', 'ENSG00000019485', 1827],
        ['A', 'ENSG00000019102', 999],
        ['X', 'ENSG00000012174', 1599],
        ['X', 'ENSG00000010671', 1977],
        ['A', 'ENSG00000019186', 1554],
        ['A', 'ENSG00000019144', 4185],
        ['X', 'ENSG00000008056', 2307],
        ['A', 'ENSG00000018408', 1383],
        ['A', 'ENSG00000019169', 1698],
    ]
    t2 = Table(header=t2_header, rows=t2_rows)

    def test_count(self):
        """test the table count methods"""
        t1 = self.t1
        self.assertEqual(t1.count('foo == "abc"'), 2)
        self.assertEqual(t1.count('foo == "cab"'), 1)
        self.assertEqual(t1.count('bar % 2 == 0'), 2)
        self.assertEqual(t1.count('id == 0'), 0)

        t2 = self.t2
        self.assertEqual(t2.count('chrom == "X"'), 4)
        self.assertEqual(t2.count('stableid.endswith("6")'), 2)
        self.assertEqual(t2.count('length % 2 == 0'), 2)
        self.assertEqual(t2.count('chrom == "Y"'), 0)
        self.assertEqual(t2.count('length % 2 == 0 and chrom == "A"'), 2)
        self.assertEqual(t2.count('length % 2 == 0 or chrom == "X"'), 6)


if __name__ == "__main__":
    main()
