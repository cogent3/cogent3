from unittest import TestCase

from cogent3.core.table import Table


class FormatBedgraph(TestCase):
    def test_only_required_columns(self):
        """generate bedgraph from minimal data"""
        table = Table(
            header=["chrom", "start", "end", "value"],
            data=[["1", 100, i, 0] for i in range(101, 111)]
            + [["1", 150, i, 10] for i in range(151, 161)],
        )

        bgraph = table.to_string(
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
        )
        assert bgraph, "\n".join(
            [
                'track type=bedGraph name="test track" '
                'description="test of bedgraph" color=255,0,0',
                "1\t100\t110\t0",
                "1\t150\t160\t10",
            ],
        )

    def test_merged_overlapping_spans(self):
        """bedgraph merged overlapping spans, one chrom"""
        rows = [["1", i, i + 1, 0] for i in range(100, 121)] + [
            ["1", i, i + 1, 10] for i in range(150, 161)
        ]
        table = Table(header=["chrom", "start", "end", "value"], data=rows)

        bgraph = table.to_string(
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
        )
        assert bgraph, "\n".join(
            [
                'track type=bedGraph name="test track" '
                'description="test of bedgraph" color=255,0,0',
                "1\t100\t120\t0",
                "1\t150\t160\t10",
            ],
        )

    def test_merged_overlapping_spans_multichrom(self):
        """bedgraph merged overlapping spans, two crhoms"""
        rows = [["1", i, i + 1, 0] for i in range(100, 121)] + [
            ["1", i, i + 1, 10] for i in range(150, 161)
        ]
        rows += [["2", i, i + 1, 0] for i in range(100, 121)]
        table = Table(header=["chrom", "start", "end", "value"], data=rows)
        bgraph = table.to_string(
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
        )

        assert bgraph, "\n".join(
            [
                'track type=bedGraph name="test track" '
                'description="test of bedgraph" color=255,0,0',
                "1\t100\t120\t1",
                "1\t150\t160\t10",
                "2\t105\t120\t1",
            ],
        )

    def test_invalid_args_fail(self):
        """incorrect bedgraph args causes RuntimeError"""
        rows = [["1", i, i + 1, 0] for i in range(100, 121)] + [
            ["1", i, i + 1, 10] for i in range(150, 161)
        ]
        table = Table(header=["chrom", "start", "end", "value"], data=rows)

        self.assertRaises(
            RuntimeError,
            table.to_string,
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
            abc=None,
        )

    def test_invalid_table_fails(self):
        """assertion error if table has > 4 columns"""
        rows = [["1", i, i + 1, 0, 1] for i in range(100, 121)] + [
            ["1", i, i + 1, 10, 1] for i in range(150, 161)
        ]
        table = Table(header=["chrom", "start", "end", "value", "blah"], data=rows)

        self.assertRaises(
            AssertionError,
            table.to_string,
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
            abc=None,
        )

    def test_boolean_correctly_formatted(self):
        """boolean setting correctly formatted"""
        rows = [["1", i, i + 1, 0] for i in range(100, 121)] + [
            ["1", i, i + 1, 10] for i in range(150, 161)
        ]
        table = Table(header=["chrom", "start", "end", "value"], data=rows)

        bgraph = table.to_string(
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
            autoScale=True,
        )

        assert bgraph, "\n".join(
            [
                'track type=bedGraph name="test track" '
                'description="test of bedgraph" color=255,0,0 autoScale=on',
                "1\t100\t110\t1",
                "1\t150\t160\t10",
            ],
        )

    def test_int_correctly_formatted(self):
        """int should be correctly formatted"""
        rows = [["1", i, i + 1, 0] for i in range(100, 121)] + [
            ["1", i, i + 1, 10] for i in range(150, 161)
        ]
        table = Table(header=["chrom", "start", "end", "value"], data=rows)

        bgraph = table.to_string(
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
            smoothingWindow=10,
        )

        assert bgraph, "\n".join(
            [
                'track type=bedGraph name="test track" '
                'description="test of bedgraph" color=255,0,0 smoothingWindow=10',
                "1\t100\t110\t1",
                "1\t150\t160\t10",
            ],
        )

    def test_raises_on_incorrect_format_val(self):
        """raise AssertionError when provide incorrect format value"""
        rows = [["1", i, i + 1, 0] for i in range(100, 121)] + [
            ["1", i, i + 1, 10] for i in range(150, 161)
        ]
        table = Table(header=["chrom", "start", "end", "value"], data=rows)

        self.assertRaises(
            AssertionError,
            table.to_string,
            format="bedgraph",
            name="test track",
            description="test of bedgraph",
            color=(255, 0, 0),
            windowingFunction="sqrt",
        )
