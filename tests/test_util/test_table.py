#!/usr/bin/env python

"""Unit tests for table.
"""
import os

from pandas import DataFrame

from cogent3 import load_table, make_table
from cogent3.util.table import Table
from cogent3.util.unit_test import TestCase, main


__author__ = "Thomas La"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Thomas La", "Christopher Bradley"]
__license__ = "BSD-3"
__version__ = "2019.11.15.a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class TrapOutput:
    def __call__(self, data):
        self.data, _ = data._get_repr_()
        self.output = repr(data)


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

    # test table 6
    t6_header = ["id", "foo", "bar"]
    t6_rows = [["60", " | ", "666"], ["70", "bca", "777"]]

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

    def test_grid_table_format(self):
        """test the table grid_table_format method"""
        from cogent3.format.table import grid_table_format

        formatted_grid = grid_table_format(
            self.t6_header, self.t6_rows, title="Test", legend="Units"
        )
        self.assertEqual(len(formatted_grid.split("\n")), len(self.t6_rows) * 2 + 7)

        formatted_grid = grid_table_format(
            self.t6_header,
            self.t6_rows,
            title="Really Long Title",
            legend="Extra Long Legend",
        )
        self.assertEqual(len(formatted_grid.split("\n")), len(self.t6_rows) * 2 + 7 + 2)

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

    def test_markdown(self):
        """Exercising the table markdown method"""
        from cogent3.format.table import markdown

        markdown_table = markdown(self.t6_header, self.t6_rows, justify="crl")
        markdown_list = markdown_table.split("\n")
        self.assertEqual(markdown_list[2].count(r"|"), 5)
        # the pipe symbol should have been escaped
        self.assertEqual(markdown_list[2].count(r"\|"), 1)

        with self.assertRaises(ValueError):
            _ = markdown(self.t6_header, self.t6_rows, justify="cr1")

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
        self.assertEqual(t5.sorted().tolist("a"), [1, 1, 2])
        self.assertEqual(t5.sorted(reverse="a").tolist("a"), [2, 1, 1])

        path = os.path.dirname(os.path.dirname(__file__))
        path = os.path.join(path, "data/sample.tsv")
        table = load_table(path)

        table = table.sorted(columns=["chrom", "stableid"])
        last_index = len(table) - 1
        self.assertEqual(table[0, "stableid"], "ENSG00000018408")
        self.assertEqual(table[last_index, "stableid"], "ENSG00000012174")

        table = table.sorted(reverse="stableid")
        self.assertEqual(table[0, "stableid"], "ENSG00000019485")
        self.assertEqual(table[last_index, "stableid"], "ENSG00000005893")

        table = table.sorted(reverse="chrom", columns="length")
        self.assertEqual(table[0, "stableid"], "ENSG00000019102")
        self.assertEqual(table[last_index, "stableid"], "ENSG00000019144")

    def test_summed(self):
        """test the table summed method"""
        t5 = Table(header=self.t5_header, rows=self.t5_rows)
        self.assertEqual(t5.summed(), [4, 4, 4, 4])
        self.assertEqual(t5.summed(col_sum=False), [4, 4, 8])
        t2 = Table(header=self.t2_header, rows=self.t2_rows)
        self.assertEqual(t2.summed(indices=2), 165)

        mix = make_table(header=["A", "B"], rows=[[0, ""], [1, 2], [3, 4]])
        self.assertEqual(mix.summed("B", strict=False), 6)
        self.assertEqual(mix.summed(0, col_sum=False, strict=False), 0)
        self.assertEqual(mix.summed(1, col_sum=False), 3)
        self.assertEqual(mix.summed(strict=False), [4, 6])
        self.assertEqual(mix.summed(col_sum=False, strict=False), [0, 3, 7])
        with self.assertRaises(RuntimeError):
            _ = mix.summed([0, 2], col_sum=False, strict=False)
        with self.assertRaises(TypeError):
            _ = mix.summed(strict=True)

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

    def test_formats(self):
        """exercising the different supported formats"""
        last = ""
        for format, startwith in (
            ("md", "|"),
            ("rst", "+"),
            ("latex", r"\begin"),
            ("markdown", "|"),
            ("csv", "id,"),
            ("tsv", "id\t"),
            ("simple", "="),
        ):
            t3 = Table(self.t3_header, rows=self.t3_rows, format=format)
            got = str(t3)
            self.assertIsInstance(got, str)
            self.assertNotEqual(got, last)
            self.assertTrue(got.startswith(startwith))
            last = got

    def test_separator_format(self):
        """testing separator_format with title and legend, and contents that match the separator"""
        from cogent3.format.table import separator_format

        with self.assertRaises(RuntimeError):
            _ = separator_format(self.t6_header, self.t6_rows)
        separated_table = separator_format(
            self.t6_header, self.t6_rows, sep=" | ", title="Test", legend="Units"
        )
        self.assertEqual(len(separated_table.split("\n")), len(self.t6_rows) + 3)

    def test_separator_format_writer(self):
        """exercising separator_format_writer"""
        from cogent3.format.table import separator_formatter

        t3 = Table(self.t3_header, rows=self.t3_rows)
        comma_sep = t3.to_string(sep=",").splitlines()
        writer = separator_formatter(sep=" | ")
        formatted = [
            f for f in writer([l.split(",") for l in comma_sep], has_header=True)
        ]
        expected_format = ["id | foo | bar", " 6 | abc |  66", " 7 | bca |  77"]
        self.assertEqual(formatted, expected_format)

    def test_set_repr_policy(self):
        """exercising setting repr policy"""
        t = Table(self.t2_header, rows=self.t2_rows)
        t.set_repr_policy(random=2)
        r = repr(t)
        self.assertIsInstance(r, str)
        r, _ = t._get_repr_()
        self.assertEqual(r.shape[0], 2)
        t.set_repr_policy(head=1)
        r, _ = t._get_repr_()
        self.assertEqual(r.shape[0], 1)
        t.set_repr_policy(tail=3)
        r, _ = t._get_repr_()
        self.assertEqual(r.shape[0], 3)

    def test_head(self):
        """returns the head of the table!"""
        from cogent3.util import table

        display = table.display
        head = TrapOutput()
        table.display = head
        t = Table(self.t1_header, rows=self.t1_rows)
        t.head(nrows=3)
        self.assertEqual(head.data.shape[0], 3)
        self.assertEqual(len(head.output.splitlines()), 9)
        self.assertEqual(head.data.tolist(), self.t1_rows[:3])
        table.display = display

    def test_tail(self):
        """returns the tail of the table!"""
        from cogent3.util import table

        display = table.display
        tail = TrapOutput()
        table.display = tail
        t = Table(self.t1_header, rows=self.t1_rows)
        t.tail(nrows=3)
        self.assertEqual(tail.data.shape[0], 3)
        self.assertEqual(len(tail.output.splitlines()), 9)
        self.assertEqual(tail.data.tolist(), self.t1_rows[-3:])
        table.display = display

    def test_to_dataframe(self):
        """produces a dataframe"""
        t = Table(self.t1_header, rows=self.t1_rows)
        df = t.to_dataframe()
        self.assertIsInstance(df, DataFrame)
        data = df.to_numpy()
        self.assertEqual(data.tolist(), self.t1_rows)

    def test_load_table(self):
        """exercising load table"""
        path = os.path.dirname(os.path.dirname(__file__))
        path = os.path.join(path, "data/sample.tsv")
        table = load_table(path)
        self.assertEqual(table.shape, (10, 3))


if __name__ == "__main__":
    main()
