#!/usr/bin/env python

"""Unit tests for table.
"""
from cogent3.util.table import Table
from cogent3.util.unit_test import TestCase, main


__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2019.8.20a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class TableTests(TestCase):
    """Tests of individual functions in table"""

    # test table 1
    t1_header = ["chrom", "stableid", "length"]
    t1_rows = [
        ["X", "ENSG00000005893", 1353],
        ["A", "ENSG00000019485", 1827],
        ["A", "ENSG00000019102", 999],
        ["X", "ENSG00000012174", 1599],
        ["X", "ENSG00000010671", 1977],
        ["A", "ENSG00000019186", 1554],
        ["A", "ENSG00000019144", 4185],
        ["X", "ENSG00000008056", 2307],
        ["A", "ENSG00000018408", 1383],
        ["A", "ENSG00000019169", 1698],
    ]

    # test table 2
    t2_header = ["id", "foo", "bar"]
    t2_rows = [
        [1, "abc", 11],
        [2, "bca", 22],
        [3, "cab", 33],
        [4, "abc", 44],
        [5, "bca", 55],
    ]

    # test table 3
    t3_header = ["id", "foo", "bar"]
    t3_rows = [[6, "abc", 66], [7, "bca", 77]]

    # test table 4
    t4_header = ["id", "foo", "bar"]
    t4_rows = [[8, "abc", 88], [9, "bca", 99]]

    # test table 5
    t5_header = ["a", "b", "c", "d"]
    t5_rows = [[1, 1, 1, 1], [2, 0, 1, 1], [1, 3, 2, 2]]

    def test_appended(self):
        """test the table appended method"""
        t2 = Table(header=self.t2_header, rows=self.t2_rows)
        t3 = Table(header=self.t3_header, rows=self.t3_rows)
        t4 = Table(header=self.t4_header, rows=self.t4_rows)

        append_0 = t2.appended("foo", [], title="self")
        self.assertEqual(append_0.shape[0], t2.shape[0])
        # test the title feature
        self.assertEqual(append_0.title, "self")

        append_1 = t2.appended("foo", [t3])
        self.assertEqual(append_1.shape[0], t2.shape[0] + t3.shape[0])
        # test the new_column feature
        self.assertEqual(append_1.shape[1], 4)
        self.assertEqual(append_1.header[0], "foo")

        append_2 = t2.appended("foo", [t3, t4])
        self.assertEqual(append_2.shape[0], t2.shape[0] + t3.shape[0] + t4.shape[0])

    def test_count(self):
        """test the table count method"""
        t1 = Table(header=self.t1_header, rows=self.t1_rows)
        self.assertEqual(t1.count('chrom == "X"'), 4)
        self.assertEqual(t1.count('stableid.endswith("6")'), 2)
        self.assertEqual(t1.count("length % 2 == 0"), 2)
        self.assertEqual(t1.count('chrom == "Y"'), 0)
        self.assertEqual(t1.count('length % 2 == 0 and chrom == "A"'), 2)
        self.assertEqual(t1.count('length % 2 == 0 or chrom == "X"'), 6)

        t2 = Table(header=self.t2_header, rows=self.t2_rows)
        self.assertEqual(t2.count('foo == "abc"'), 2)
        self.assertEqual(t2.count('foo == "cab"'), 1)
        self.assertEqual(t2.count("bar % 2 == 0"), 2)
        self.assertEqual(t2.count("id == 0"), 0)

    def test_distinct_values(self):
        """test the table distinct_values method"""
        t1 = Table(header=self.t1_header, rows=self.t1_rows)
        self.assertEqual(len(t1.distinct_values("chrom")), 2)
        self.assertEqual(len(t1.distinct_values("stableid")), 10)
        self.assertEqual(len(t1.distinct_values("length")), 10)

        t2 = Table(header=self.t2_header, rows=self.t2_rows)
        self.assertEqual(len(t2.distinct_values("id")), 5)
        self.assertEqual(len(t2.distinct_values("foo")), 3)
        self.assertEqual(len(t2.distinct_values("bar")), 5)

    def test_filtered(self):
        """test the table filtered method"""
        t1 = Table(header=self.t1_header, rows=self.t1_rows)
        self.assertEqual(t1.filtered('chrom == "X"').shape[0], 4)
        self.assertEqual(t1.filtered('stableid.endswith("6")').shape[0], 2)
        self.assertEqual(t1.filtered("length % 2 == 0").shape[0], 2)
        self.assertEqual(t1.filtered('chrom == "Y"').shape[0], 0)
        self.assertEqual(t1.filtered('length % 2 == 0 and chrom == "A"').shape[0], 2)
        self.assertEqual(t1.filtered('length % 2 == 0 or chrom == "X"').shape[0], 6)

        t2 = Table(header=self.t2_header, rows=self.t2_rows)
        self.assertEqual(t2.filtered('foo == "abc"').shape[0], 2)
        self.assertEqual(t2.filtered('foo == "cab"').shape[0], 1)
        self.assertEqual(t2.filtered("bar % 2 == 0").shape[0], 2)
        self.assertEqual(t2.filtered("id == 0").shape[0], 0)

    def test_filtered_by_column(self):
        """test the table filtered_by_column method"""
        t1 = Table(header=self.t1_header, rows=self.t1_rows)
        t2 = Table(header=self.t2_header, rows=self.t2_rows)

        def is_numeric(values):
            try:
                sum(values)
            except TypeError:
                return False
            return True

        self.assertEqual(t1.filtered_by_column(is_numeric).shape[1], 1)
        self.assertEqual(t2.filtered_by_column(is_numeric).shape[1], 2)

    def test_get_columns(self):
        """test the table get_columns method"""
        t1 = Table(header=self.t1_header, rows=self.t1_rows)
        self.assertEqual(t1.get_columns("chrom").shape[0], t1.shape[0])
        self.assertEqual(t1.get_columns("chrom").shape[1], 1)

        self.assertEqual(t1.get_columns(["chrom", "length"]).shape[0], t1.shape[0])
        self.assertEqual(t1.get_columns(["chrom", "length"]).shape[1], 2)

    def test_joined(self):
        """test the table joined method"""
        t2 = Table(header=self.t2_header, rows=self.t2_rows)
        t3 = Table(header=self.t3_header, rows=self.t3_rows)

        # inner join test
        self.assertEqual(
            t2.joined(t3, columns_self="foo", columns_other="foo").shape[0], 4
        )
        # merged 'foo' column, so (6-1) columns in join
        self.assertEqual(
            t2.joined(t3, columns_self="foo", columns_other="foo").shape[1], 5
        )

        # non-inner join test (cartesian product of rows)
        self.assertEqual(
            t2.joined(t3, inner_join=False).shape[0], t2.shape[0] * t3.shape[0]
        )
        self.assertEqual(
            t2.joined(t3, inner_join=False).shape[1], t2.shape[1] + t3.shape[1]
        )

    def test_normalized(self):
        """test the table normalized method"""
        t5 = Table(header=self.t5_header, rows=self.t5_rows)
        self.assertEqual(
            t5.normalized().tolist(t5.header),
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.5, 0.0, 0.25, 0.25],
                [0.125, 0.375, 0.25, 0.25],
            ],
        )
        self.assertEqual(
            t5.normalized(by_row=False).tolist(t5.header),
            [[0.25, 0.25, 0.25, 0.25], [0.5, 0.0, 0.25, 0.25], [0.25, 0.75, 0.5, 0.5]],
        )

    def test_sorted(self):
        """test the table sorted method"""
        t1 = Table(header=self.t1_header, rows=self.t1_rows)
        self.assertEqual(
            t1.sorted("length").tolist("length"),
            [999, 1353, 1383, 1554, 1599, 1698, 1827, 1977, 2307, 4185],
        )

        t5 = Table(header=self.t5_header, rows=self.t5_rows)
        self.assertEqual(t5.sorted("b").tolist("b"), [0, 1, 3])

    def test_summed(self):
        """test the table summed method"""
        t5 = Table(header=self.t5_header, rows=self.t5_rows)
        self.assertEqual(t5.summed(), [4, 4, 4, 4])
        self.assertEqual(t5.summed(col_sum=False), [4, 4, 8])

    def test_tolist(self):
        """test the table tolist method"""
        t3 = Table(header=self.t3_header, rows=self.t3_rows)
        self.assertEqual(t3.tolist("id"), [6, 7])
        self.assertEqual(t3.tolist("foo"), ["abc", "bca"])

    def test_transposed(self):
        """test the table transposed method"""
        t1 = Table(header=self.t1_header, rows=self.t1_rows)
        # note the transpose has one less row
        self.assertEqual(t1.transposed("").shape[0], t1.shape[1] - 1)
        # note the transpose has an extra column
        self.assertEqual(t1.transposed("").shape[1], t1.shape[0] + 1)

    def test_with_new_column(self):
        """test the table with_new_column method"""
        t5 = Table(header=self.t5_header, rows=self.t5_rows)
        t5_row_sum = t5.with_new_column("sum", sum, t5.header)
        self.assertEqual(t5_row_sum.get_columns("sum").tolist(), [[4], [4], [8]])

    def test_with_new_header(self):
        """test the table with_new_header method"""
        t2 = Table(header=self.t2_header, rows=self.t2_rows)
        t2 = t2.with_new_header("id", "no")
        self.assertEqual(t2.header[0], "no")
        t2 = t2.with_new_header("foo", "moo")
        self.assertEqual(t2.header[1], "moo")
        t2 = t2.with_new_header("moo", "foo")
        self.assertEqual(t2.header[1], "foo")


if __name__ == "__main__":
    main()
