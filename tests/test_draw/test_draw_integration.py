import pathlib
import unittest

from numpy.testing import assert_allclose

from cogent3 import load_aligned_seqs, make_aligned_seqs, make_table
from cogent3.draw.drawable import AnnotatedDrawable, Drawable, get_domain
from cogent3.util.union_dict import UnionDict


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2020.2.7a"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


def load_alignment(annotate1=False, annotate2=False):
    """creates an alignment with None, one or two sequences annotated"""
    path = str(pathlib.Path(__file__).parent.parent / "data/brca1_5.paml")
    aln = load_aligned_seqs(path, array_align=False, moltype="dna")
    aln = aln.omit_gap_pos()
    if annotate1:
        x1 = aln.get_seq(aln.names[0]).add_feature("gene", "abcde1", [(20, 50)])
        x2 = aln.get_seq(aln.names[0]).add_feature("variation", "one", [(11, 12)])
    if annotate2:
        y1 = aln.get_seq(aln.names[1]).add_feature("gene", "abcde2", [(20, 50)])
        y2 = aln.get_seq(aln.names[1]).add_feature("domain", "abcde2", [(10, 15)])
    return aln


class UtilDrawablesTests(unittest.TestCase):
    """testing utility functions"""

    def test_get_domain_1(self):
        """handles single domain"""
        x = get_domain(1, 0, is_y=False)
        y = get_domain(1, 0, is_y=True)
        assert_allclose(x, [0, 1])
        assert_allclose(y, [0, 1])

    def test_get_domain_multi(self):
        """produces reversed coordinates for y"""
        # plotly coords are cartesian, index coords are array, so they need to
        # be converted
        y = get_domain(2, 0, is_y=True, space=0)
        assert_allclose(y, [0.5, 1])
        y = get_domain(2, 1, is_y=True, space=0)
        assert_allclose(y, [0, 0.5])
        # x-coords unaffected
        x = get_domain(2, 0, is_y=False, space=0)
        assert_allclose(x, [0, 0.5])
        x = get_domain(2, 1, is_y=False, space=0)
        assert_allclose(x, [0.5, 1])

    def test_get_domain_sep(self):
        """space argument shifts boundaries"""
        sep = 0.1
        x = get_domain(2, 0, is_y=False, space=sep)
        assert_allclose(x, [0 + sep / 2, 0.5 - sep / 2])
        x = get_domain(2, 1, is_y=False, space=sep)
        assert_allclose(x, [0.5 + sep / 2, 1 - sep / 2])
        # if space too big relative to the span of each domain, it's reduced
        # to 1/10th the domain span
        sep = 0.6
        x = get_domain(2, 0, is_y=False, space=sep)
        exp_sep = 0.5 / 10
        assert_allclose(x, [0 + exp_sep, 0.5 - exp_sep])

    def test_domain_element_size(self):
        """domain element value must not exceed num domains - 1"""
        with self.assertRaises(ValueError):
            x = get_domain(2, 2, is_y=False)


class BaseDrawablesTests(unittest.TestCase):
    """methods for checking drawables"""

    def _check_drawable_attrs(self, fig, type_):
        self.assertIsInstance(fig, UnionDict)
        # should have a layout and a data key
        self.assertTrue("layout" in fig)
        self.assertTrue("data" in fig)
        # data traces should be of type "scatter"
        self.assertEqual({tr.type for tr in fig.data}, {type_})

    def _check_drawable_styles(self, method, styles, **kwargs):
        for style in styles:
            obj = method(drawable=style, **kwargs)
            self._check_drawable_attrs(obj.drawable.figure, style)


class CustomDrawable(BaseDrawablesTests):
    """custom applications of Drawable"""

    def test_no_trace(self):
        """should still produce a valid figure"""
        d = Drawable()
        f = d.figure
        self.assertEqual(f.data, [{}])

    def test_one_trace(self):
        """should still produce a valid figure"""
        d = Drawable()
        trace = dict(type="scatter", x=[0, 1], y=[0, 1])
        d.add_trace(trace)
        f = d.figure
        self.assertEqual(f.data, [trace])

    def test_layout_with_titles(self):
        """if provided layout has axis titles, keep them"""
        layout = dict(xaxis=dict(title="X"), yaxis=dict(title="Y"))
        d = Drawable(layout=layout)
        fig = d.figure
        self.assertEqual(fig.layout.xaxis.title, "X")
        self.assertEqual(fig.layout.yaxis.title, "Y")
        d = Drawable()
        fig = d.figure
        self.assertEqual(fig.layout.xaxis.title, None)
        self.assertEqual(fig.layout.yaxis.title, None)


class AlignmentDrawablesTests(BaseDrawablesTests):
    """methods on SequenceCollection produce Drawables"""

    def _check_drawable_attrs(self, fig, type_):
        self.assertIsInstance(fig, UnionDict)
        # should have a layout and a data key
        self.assertTrue("layout" in fig)
        self.assertTrue("data" in fig)
        # data traces should be of type "scatter"
        self.assertEqual({tr.type for tr in fig.data}, {type_})

    def _check_drawable_styles(self, method, styles, **kwargs):
        for style in styles:
            obj = method(drawable=style, **kwargs)
            self._check_drawable_attrs(obj.drawable.figure, style)

    def test_dotplot_regression(self):
        """Tests whether dotplot produces traces and in correct ordering. Also tests if pop_trace() works"""
        aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
        aln = aln.take_seqs(["Human", "Chimpanzee"])
        aln = aln[:200]
        dp = aln.dotplot()
        _ = dp.figure
        trace_names = dp.get_trace_titles()

        self.assertTrue(
            dp.get_trace_titles() != [] and len(trace_names) == len(dp.traces),
            "No traces found for dotplot",
        )
        self.assertTrue(
            [trace_names[i] == dp.traces[i]["name"] for i in range(len(trace_names))],
            "Order of traces don't match with get_trace_titles()",
        )

        for trace_name in trace_names:
            dp.pop_trace(trace_name)
            self.assertFalse(
                trace_name in dp.get_trace_titles(),
                "Trace name still present in get_trace_titles() even after popping off trace",
            )

    def test_dotplot_annotated(self):
        """alignment with / without annotated seqs"""
        from cogent3.draw.dotplot import Dotplot

        aln = load_alignment()
        # no annotation
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        self._check_drawable_attrs(dp.figure, "scatter")
        self.assertIsInstance(dp, Dotplot)

        # seq1 annotated
        aln = load_alignment(annotate1=True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        self.assertIsInstance(dp, AnnotatedDrawable)
        fig = dp.figure
        self._check_drawable_attrs(fig, "scatter")
        self.assertIs(dp.left_track, None)
        self.assertIsInstance(dp.bottom_track, Drawable)

        # seq2 annotated
        aln = load_alignment(annotate2=True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        self.assertIsInstance(dp, AnnotatedDrawable)
        fig = dp.figure
        self._check_drawable_attrs(fig, "scatter")
        self.assertIs(dp.bottom_track, None)
        self.assertIsInstance(dp.left_track, Drawable)

        # both annotated
        aln = load_alignment(True, True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        self.assertIsInstance(dp, AnnotatedDrawable)
        fig = dp.figure
        self._check_drawable_attrs(fig, "scatter")
        self.assertIsInstance(dp.bottom_track, Drawable)
        self.assertIsInstance(dp.left_track, Drawable)

    def test_annotated_dotplot_remove_tracks(self):
        """removing annotation tracks from dotplot should work"""
        # both annotated
        aln = load_alignment(True, True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        orig_fig = dp.figure
        # remove left
        dp.remove_track(left_track=True)
        self.assertIs(dp.left_track, None)
        self.assertEqual(dp._traces, [])
        self.assertIsNot(dp.figure, orig_fig)
        self.assertIsInstance(dp.figure, UnionDict)

        # both annotated, remove right
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        dp.remove_track(bottom_track=True)
        self.assertIs(dp.bottom_track, None)
        self.assertEqual(dp._traces, [])
        self.assertIsNot(dp.figure, orig_fig)
        self.assertIsInstance(dp.figure, UnionDict)

        # both annotated, remove both
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        dp.remove_track(True, True)
        self.assertIs(dp.bottom_track, None)
        self.assertIs(dp.left_track, None)
        self.assertEqual(dp._traces, [])
        self.assertIsNot(dp.figure, orig_fig)
        self.assertIsInstance(dp.figure, UnionDict)

    def test_count_gaps_per_seq(self):
        """creation of drawables works"""
        styles = "bar", "box", "violin"
        aln = load_alignment(False, False)
        # no drawable
        counts = aln.count_gaps_per_seq()
        self.assertFalse(hasattr(counts, "drawable"))
        self._check_drawable_styles(aln.count_gaps_per_seq, styles)
        # now ArrayAlignment
        aln = aln.to_type(array_align=True)
        self._check_drawable_styles(aln.count_gaps_per_seq, styles)

    def test_coevo_drawables(self):
        """coevolution produces drawables"""
        styles = "box", "heatmap", "violin"
        aln = load_alignment(False, False)
        aln = aln[:10]
        coevo = aln.coevolution(show_progress=False)
        self.assertFalse(hasattr(coevo, "drawable"))
        self._check_drawable_styles(aln.coevolution, styles, show_progress=False)
        # now ArrayAlignment
        aln = aln.to_type(array_align=True)
        self._check_drawable_styles(aln.coevolution, styles, show_progress=False)

    def test_coevo_annotated(self):
        """coevolution on alignment with annotated seqs should add to heatmap plot"""
        aln = load_alignment(True, False)
        aln = aln[30:]
        coevo = aln.coevolution(show_progress=False, drawable="heatmap")
        drawable = coevo.drawable
        self.assertIsInstance(drawable, AnnotatedDrawable)
        self.assertIsInstance(drawable.left_track, Drawable)
        self.assertIsInstance(drawable.bottom_track, Drawable)

    def test_information_plot(self):
        """infoprmation plot makes a drawable"""
        aln = load_alignment(False, False)
        drawable = aln.information_plot()
        self.assertIsInstance(drawable, Drawable)
        self._check_drawable_attrs(drawable.figure, "scatter")
        # ArrayAlignment
        aln = aln.to_type(array_align=True)
        drawable = aln.information_plot()
        self._check_drawable_attrs(drawable.figure, "scatter")

    def test_get_drawable(self):
        """sliced alignment with features returns a drawable"""
        aln = make_aligned_seqs(
            data=dict(a="AAACGGTTT", b="CAA--GTAA"), array_align=False
        )
        _ = aln.get_seq("b").add_feature("domain", "1", [(1, 5)])
        _ = aln.get_seq("b").add_feature("variation", "1", [(1, 5)])
        _ = aln.get_seq("b").add_feature("gene", "1", [(1, 5)])
        _ = aln.get_seq("b").add_feature("gene", "1", [(5, 1)])
        drawable = aln.get_drawable()


class TableDrawablesTest(BaseDrawablesTests):
    """single Table method produces a Drawable"""

    def test_to_plotly(self):
        """exercise producing a plotly table"""
        table = make_table(header=["a", "b"], data=[[0, 1]], row_ids="a")
        drawable = table.to_plotly()
        self.assertIsInstance(drawable, Drawable)
        self._check_drawable_attrs(drawable.figure, "table")


if __name__ == "__main__":
    unittest.main()
