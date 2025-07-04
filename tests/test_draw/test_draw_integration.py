import pathlib
import unittest
import uuid

import pytest
from numpy.testing import assert_allclose

from cogent3 import load_aligned_seqs, make_aligned_seqs, make_table
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.draw.drawable import (
    AnnotatedDrawable,
    Arrow,
    Drawable,
    _calc_arrow_width,
    get_domain,
    make_shape,
)
from cogent3.util.union_dict import UnionDict


def load_alignment(annotate1=False, annotate2=False):
    """creates an alignment with None, one or two sequences annotated"""
    db = GffAnnotationDb()

    path = str(pathlib.Path(__file__).parent.parent / "data/brca1_5.paml")
    aln = load_aligned_seqs(path, moltype="dna")
    aln = aln.omit_gap_pos()
    if annotate1:
        db.add_feature(
            seqid=aln.names[0],
            biotype="gene",
            name="abcde1",
            spans=[(20, 50)],
        )
        db.add_feature(
            seqid=aln.names[0],
            biotype="variation",
            name="one",
            spans=[(11, 12)],
        )
    if annotate2:
        db.add_feature(
            seqid=aln.names[1],
            biotype="gene",
            name="abcde2",
            spans=[(20, 50)],
        )

        db.add_feature(
            seqid=aln.names[1],
            biotype="domain",
            name="abcde2",
            spans=[(10, 15)],
        )
    aln.annotation_db = db
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
        with pytest.raises(ValueError):
            get_domain(2, 2, is_y=False)


class DrawableObjectTests(unittest.TestCase):
    """testing Drawable object methods and properties"""

    def test_traces(self):
        """test trace initialisation"""
        d = Drawable()
        assert d.traces == []
        trace = {"type": "scatter", "x": [0, 1], "y": [0, 1]}
        d = Drawable(traces=[trace])
        assert d.traces == [trace]
        assert isinstance(d.traces[0], UnionDict)
        with pytest.raises(TypeError):
            trace = {"type": "scatter", "x": [0, 1], "y": [0, 1]}
            _ = Drawable(traces=trace)

    def test_add_traces(self):
        """test trace add method"""
        d = Drawable()
        assert d.traces == []
        trace = {"type": "scatter", "x": [0, 1], "y": [0, 1]}
        d.add_trace(trace)
        assert d.traces == [trace]
        assert isinstance(d.traces[0], UnionDict)

    def test_bound_to(self):
        """bound object should have the drawable object and show methods"""

        class TestObj:
            pass

        o = TestObj()
        d = Drawable()
        b = d.bound_to(o)
        assert b.drawable == d
        assert hasattr(b, "iplot")
        assert hasattr(b, "show")

    def test_figure(self):
        """figure should contain the same data and layout as Drawable"""
        trace = {"type": "scatter", "x": [0, 1], "y": [0, 1]}
        layout = {"title": "layout", "width": 20}
        d = Drawable(traces=[trace], layout=layout)
        f = d.figure
        assert f.data == d.traces
        assert f.layout == d.layout

    def test_plotly_figure(self):
        """is a plotly graph object Figure instance"""
        from plotly.graph_objects import Figure

        trace = {"type": "scatter", "x": [0, 1], "y": [0, 1]}
        layout = {"title": "layout", "width": 20}
        d = Drawable(traces=[trace], layout=layout)
        assert isinstance(d.plotly_figure, Figure)


class AnnotatedDrawableObjectTests(unittest.TestCase):
    """testing AnnotatedDrawable object methods and properties"""

    def test__build_fig(self):
        """figure built should have the same traces as core and set yaxis to y3 if yaxis2 is overlaying"""
        trace = {
            "type": "scatter",
            "x": [0, 1],
            "y": [0, 1],
            "xaxis": "x",
            "yaxis": "y",
        }
        layout = {"title": "layout", "width": 20, "yaxis2": {"overlaying": "free"}}
        cd = Drawable(traces=[trace])

        ad = AnnotatedDrawable(cd, layout=layout)
        f = ad._build_fig()
        assert ad._traces == f["data"]
        assert f["data"][0]["yaxis"] != "y3"

        layout = {"title": "layout", "width": 20, "yaxis2": {"overlaying": "y"}}
        ad = AnnotatedDrawable(cd, layout=layout)
        f = ad._build_fig()
        assert f["data"][0]["yaxis"] == "y3"

    def test_plotly_figure(self):
        """is a plotly graph object Figure instance"""
        from plotly.graph_objects import Figure

        trace = {
            "type": "scatter",
            "x": [0, 1],
            "y": [0, 1],
            "xaxis": "x",
            "yaxis": "y",
        }
        layout = {"title": "layout", "width": 20, "yaxis2": {"overlaying": "free"}}
        cd = Drawable(traces=[trace])

        ad = AnnotatedDrawable(cd, layout=layout)
        assert isinstance(ad.plotly_figure, Figure)


class BaseDrawablesTests(unittest.TestCase):
    """methods for checking drawables"""

    def _check_drawable_attrs(self, fig, type_):
        assert isinstance(fig, UnionDict)
        # should have a layout and a data key
        assert "layout" in fig
        assert "data" in fig
        # data traces should be of type "scatter"
        assert {tr.type for tr in fig.data} == {type_}

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
        assert f.data == [{}]

    def test_one_trace(self):
        """should still produce a valid figure"""
        d = Drawable()
        trace = {"type": "scatter", "x": [0, 1], "y": [0, 1]}
        d.add_trace(trace)
        f = d.figure
        assert f.data == [trace]

    def test_layout_with_titles(self):
        """if provided layout has axis titles, keep them"""
        layout = {"xaxis": {"title": "X"}, "yaxis": {"title": "Y"}}
        d = Drawable(layout=layout)
        fig = d.figure
        assert fig.layout.xaxis.title == "X"
        assert fig.layout.yaxis.title == "Y"
        d = Drawable()
        fig = d.figure
        assert fig.layout.xaxis.title is None
        assert fig.layout.yaxis.title is None

    def test_layout_with_overrides(self):
        """provided layout attributes should be overridden with provided parameters"""
        layout = {"title": "layout", "width": 20}
        d = Drawable(layout=layout)
        fig = d.figure
        assert fig.layout.title is None
        assert fig.layout.width is None
        d = Drawable(layout=layout, title="parameter", width=50)
        fig = d.figure
        assert fig.layout.title.text == "parameter"
        assert fig.layout.width == 50


class AlignmentDrawablesTests(BaseDrawablesTests):
    """methods on SequenceCollection produce Drawables"""

    def _check_drawable_attrs(self, fig, type_):
        assert isinstance(fig, UnionDict)
        # should have a layout and a data key
        assert "layout" in fig
        assert "data" in fig
        # data traces should be of type "scatter"
        assert {tr.type for tr in fig.data} == {type_}

    def _check_drawable_styles(self, method, styles, **kwargs):
        for style in styles:
            obj = method(drawable=style, **kwargs)
            self._check_drawable_attrs(obj.drawable.figure, style)

    def test_dotplot_regression(self):  # ported
        """Tests whether dotplot produces traces and in correct ordering."""
        aln = load_aligned_seqs("data/brca1.fasta", moltype="dna")
        aln = aln.take_seqs(["Human", "Chimpanzee"])
        aln = aln[:200]
        dp = aln.dotplot()
        _ = dp.figure
        trace_names = [tr.name for tr in dp.traces]

        assert [tr.name for tr in dp.traces] != [] and len(trace_names) == len(
            dp.traces,
        ), "No traces found for dotplot"
        assert [
            trace_names[i] == dp.traces[i]["name"] for i in range(len(trace_names))
        ], "Order of traces don't match with dp traces"

        for trace_name in trace_names:
            index = [tr.name for tr in dp.traces].index(trace_name)
            dp.traces.pop(index)

            assert trace_name not in [tr.name for tr in dp.traces], (
                "Trace name still present in dp traces even after popping off trace"
            )

    def test_dotplot_annotated(self):  # ported
        """alignment with / without annotated seqs"""
        from cogent3.draw.dotplot import Dotplot

        aln = load_alignment()
        # no annotation
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        self._check_drawable_attrs(dp.figure, "scatter")
        assert isinstance(dp, Dotplot)

        # seq1 annotated
        aln = load_alignment(annotate1=True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        assert isinstance(dp, AnnotatedDrawable)
        fig = dp.figure
        self._check_drawable_attrs(fig, "scatter")
        assert dp.left_track is None
        assert isinstance(dp.bottom_track, Drawable)

        # seq2 annotated
        aln = load_alignment(annotate2=True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        assert isinstance(dp, AnnotatedDrawable)
        fig = dp.figure
        self._check_drawable_attrs(fig, "scatter")
        assert dp.bottom_track is None
        assert isinstance(dp.left_track, Drawable)

        # both annotated
        aln = load_alignment(True, True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        assert isinstance(dp, AnnotatedDrawable)
        fig = dp.figure
        self._check_drawable_attrs(fig, "scatter")
        assert isinstance(dp.bottom_track, Drawable)
        assert isinstance(dp.left_track, Drawable)

    def test_annotated_dotplot_remove_tracks(self):  # ported
        """removing annotation tracks from dotplot should work"""
        # both annotated
        aln = load_alignment(True, True)
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        orig_fig = dp.figure
        # remove left
        dp.remove_track(left_track=True)
        assert dp.left_track is None
        assert dp._traces == []
        assert dp.figure is not orig_fig
        assert isinstance(dp.figure, UnionDict)

        # both annotated, remove right
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        dp.remove_track(bottom_track=True)
        assert dp.bottom_track is None
        assert dp._traces == []
        assert dp.figure is not orig_fig
        assert isinstance(dp.figure, UnionDict)

        # both annotated, remove both
        dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
        dp.remove_track(True, True)
        assert dp.bottom_track is None
        assert dp.left_track is None
        assert dp._traces == []
        assert dp.figure is not orig_fig
        assert isinstance(dp.figure, UnionDict)

    def test_count_gaps_per_seq(self):  # ported
        """creation of drawables works"""
        styles = "bar", "box", "violin"
        aln = load_alignment(False, False)
        # no drawable
        counts = aln.count_gaps_per_seq()
        assert not hasattr(counts, "drawable")
        self._check_drawable_styles(aln.count_gaps_per_seq, styles)
        self._check_drawable_styles(aln.count_gaps_per_seq, styles)

    def test_coevo_drawables(self):  # ported
        """coevolution produces drawables"""
        styles = "box", "heatmap", "violin"
        aln = load_alignment(False, False)
        aln = aln[:10]
        coevo = aln.coevolution(show_progress=False)
        assert not hasattr(coevo, "drawable")
        self._check_drawable_styles(aln.coevolution, styles, show_progress=False)
        self._check_drawable_styles(aln.coevolution, styles, show_progress=False)

    def test_coevo_annotated(self):  # ported
        """coevolution on alignment with annotated seqs should add to heatmap plot"""
        aln = load_alignment(True, False)
        aln = aln[:30]
        coevo = aln.coevolution(show_progress=False, drawable="heatmap")
        drawable = coevo.drawable
        assert isinstance(drawable, AnnotatedDrawable)
        assert isinstance(drawable.left_track, Drawable)
        assert isinstance(drawable.bottom_track, Drawable)

    def test_information_plot(self):  # ported
        """infoprmation plot makes a drawable"""
        aln = load_alignment(False, False)
        drawable = aln.information_plot()
        assert isinstance(drawable, Drawable)
        self._check_drawable_attrs(drawable.figure, "scatter")
        drawable = aln.information_plot()
        self._check_drawable_attrs(drawable.figure, "scatter")

    def test_get_drawable(self):  # ported
        """sliced alignment with features returns a drawable"""
        aln = make_aligned_seqs(
            {"a": "AAACGGTTT", "b": "CAA--GTAA"},
            moltype="dna",
        )
        db = GffAnnotationDb()
        db.add_feature(seqid="b", biotype="domain", name="1", spans=[(1, 5)])
        db.add_feature(seqid="b", biotype="variation", name="1", spans=[(1, 5)])
        db.add_feature(seqid="b", biotype="gene", name="1", spans=[(1, 5)])
        db.add_feature(seqid="b", biotype="gene", name="1", spans=[(5, 1)])
        aln.annotation_db = db
        g = aln.get_drawable()
        assert isinstance(g, Drawable)


class TableDrawablesTest(BaseDrawablesTests):
    """single Table method produces a Drawable"""

    def test_to_plotly(self):
        """exercise producing a plotly table"""
        table = make_table(header=["a", "b"], data=[[0, 1]], index_name="a")
        drawable = table.to_plotly()
        assert isinstance(drawable, Drawable)
        self._check_drawable_attrs(drawable.figure, "table")


def test_calculating_arrow_width_adjusted():
    aw = _calc_arrow_width(parent_length=100, feature_width=10, frac=1.0)
    assert aw == 10
    aw = _calc_arrow_width(parent_length=100, feature_width=10, frac=0.05)
    assert aw == 5


def test_colour_choice():
    label = str(uuid.uuid4())
    assert label not in make_shape._colors
    shape1 = make_shape(type_=label, coords=[(10, 20)])
    assert label in make_shape._colors
    shape2 = make_shape(type_=label, coords=[(30, 60)])
    assert shape1.fillcolor == shape2.fillcolor == make_shape._colors[label]


@pytest.mark.parametrize("biotype", ["mRNA", "MRNA", "gene", "misc_RNA", "tRNA", "CDS"])
def test_arrow_shapes(biotype):
    shape = make_shape(type_=biotype, coords=[(10, 20)], parent_length=20)
    assert isinstance(shape, Arrow)
