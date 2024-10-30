import pathlib

import pytest

from cogent3 import load_aligned_seqs, load_seq
from cogent3.core import new_alignment
from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.draw.dotplot import Dotplot
from cogent3.draw.drawable import AnnotatedDrawable, Drawable
from cogent3.util.union_dict import UnionDict


@pytest.fixture(scope="function")
def annotated_seq(DATA_DIR):
    return load_seq(
        DATA_DIR / "c_elegans_WS199_dna_shortened.fasta",
        annotation_path=DATA_DIR / "c_elegans_WS199_shortened_gff.gff3",
        moltype="dna",
        new_type=True,
    )


@pytest.fixture
def load_alignment():
    """Fixture to create an alignment with None, one, or two sequences annotated."""

    def _load_alignment(annotate1=False, annotate2=False):
        db = GffAnnotationDb()

        path = str(pathlib.Path(__file__).parent.parent / "data/brca1_5.paml")
        aln = load_aligned_seqs(path, moltype="dna", new_type=True)
        aln = aln.omit_gap_pos()
        if annotate1:
            db.add_feature(
                seqid=aln.names[0], biotype="gene", name="abcde1", spans=[(20, 50)]
            )
            db.add_feature(
                seqid=aln.names[0], biotype="variation", name="one", spans=[(11, 12)]
            )
        if annotate2:
            db.add_feature(
                seqid=aln.names[1], biotype="gene", name="abcde2", spans=[(20, 50)]
            )
            db.add_feature(
                seqid=aln.names[1], biotype="domain", name="abcde2", spans=[(10, 15)]
            )
        aln.annotation_db = db
        return aln

    return _load_alignment


def check_drawable_attrs(fig, type_):
    assert isinstance(fig, UnionDict), f"Expected fig to be UnionDict, got {type(fig)}"
    # should have a layout and a data key
    assert "layout" in fig, "Expected 'layout' key in fig"
    assert "data" in fig, "Expected 'data' key in fig"
    # data traces should be of type "scatter"
    assert {tr.type for tr in fig.data} == {
        type_
    }, f"Expected data traces of type {type_}"


def check_drawable_styles(method, styles, **kwargs):
    for style in styles:
        obj = method(drawable=style, **kwargs)
        check_drawable_attrs(obj.drawable.figure, style)


@pytest.fixture(scope="function")
def dotplot_seqs():
    data = {
        "Human": "CAGATTTGGCAGTT-",
        "Mouse": "CAGATTCAGCAGGTG",
        "Rat": "CAGATTCAGCAGG--G",
    }
    return new_alignment.make_unaligned_seqs(data, moltype="dna")


def test_dotplot_base_cases(dotplot_seqs):
    """exercising dotplot method"""
    # with provided names
    plot = dotplot_seqs.dotplot(name1="Human", name2="Mouse")
    assert str(plot.seq1) != str(plot.seq2)

    # without providing names
    plot = dotplot_seqs.dotplot()
    assert str(plot.seq1) != str(plot.seq2)

    # providing only one name
    plot = dotplot_seqs.dotplot(name1="Human")
    assert str(plot.seq1) != str(plot.seq2)

    # a collection of one sequence should make dotplot with itself
    less_seqs = dotplot_seqs.take_seqs("Human")
    plot = less_seqs.dotplot()
    assert str(plot.seq1) == str(plot.seq2)

    # k larger than window should raise an error
    with pytest.raises(AssertionError):
        dotplot_seqs.dotplot(window=5, k=11)

    # names not in the collection should raise an error
    with pytest.raises(ValueError):
        dotplot_seqs.dotplot(name1="Human", name2="Dog")


@pytest.mark.parametrize("with_annotations", [True, False])
def test_dotplot_annotated(annotated_seq, with_annotations):
    if not with_annotations:
        annotated_seq.replace_annotation_db(None)  # this drops all annotations

    coll = new_alignment.make_unaligned_seqs(
        {"c_elegans": annotated_seq}, moltype="dna"
    )
    dp = coll.dotplot()
    if with_annotations:
        assert len(dp.figure.data) > 2
    else:
        assert len(dp.figure.data) == 2


def test_dotplot_regression():
    """Tests whether dotplot produces traces and in correct ordering."""
    aln = load_aligned_seqs("data/brca1.fasta", moltype="dna", new_type=True)
    aln = aln.take_seqs(["Human", "Chimpanzee"])
    aln = aln[:200]
    dp = aln.dotplot()
    _ = dp.figure
    trace_names = [tr.name for tr in dp.traces]

    assert [tr.name for tr in dp.traces] != [] and len(trace_names) == len(
        dp.traces
    ), "No traces found for dotplot"
    assert all(
        trace_names[i] == dp.traces[i]["name"] for i in range(len(trace_names))
    ), "Order of traces doesn't match with dp traces"

    for trace_name in trace_names:
        index = [tr.name for tr in dp.traces].index(trace_name)
        dp.traces.pop(index)
        assert trace_name not in [
            tr.name for tr in dp.traces
        ], "Trace name still present in dp traces even after popping off trace"


def test_dotplot_with_diff_annotation_permutations(load_alignment):
    """alignment with / without annotated seqs"""

    # no annotation
    aln = load_alignment()
    dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
    check_drawable_attrs(dp.figure, "scatter")
    assert isinstance(dp, Dotplot)

    # seq1 annotated
    aln = load_alignment(annotate1=True)
    dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
    assert isinstance(dp, AnnotatedDrawable)
    fig = dp.figure
    check_drawable_attrs(fig, "scatter")
    assert dp.left_track is None
    assert isinstance(dp.bottom_track, Drawable)

    # seq2 annotated
    aln = load_alignment(annotate2=True)
    dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
    assert isinstance(dp, AnnotatedDrawable)
    fig = dp.figure
    check_drawable_attrs(fig, "scatter")
    assert dp.bottom_track is None
    assert isinstance(dp.left_track, Drawable)

    # both annotated
    aln = load_alignment(annotate1=True, annotate2=True)
    dp = aln.dotplot(name1=aln.names[0], name2=aln.names[1])
    assert isinstance(dp, AnnotatedDrawable)
    fig = dp.figure
    check_drawable_attrs(fig, "scatter")
    assert isinstance(dp.bottom_track, Drawable)
    assert isinstance(dp.left_track, Drawable)


def test_annotated_dotplot_remove_tracks(load_alignment):
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


def test_count_gaps_per_seq(load_alignment):
    """creation of drawables works"""
    styles = "bar", "box", "violin"
    aln = load_alignment(False, False)
    # no drawable
    counts = aln.count_gaps_per_seq()
    assert not hasattr(counts, "drawable")
    check_drawable_styles(aln.count_gaps_per_seq, styles)


def test_coevo_drawables(load_alignment):
    """coevolution produces drawables"""
    styles = "box", "heatmap", "violin"
    aln = load_alignment(False, False)
    aln = aln[:10]
    coevo = aln.coevolution(show_progress=False)
    assert not (hasattr(coevo, "drawable"))
    check_drawable_styles(aln.coevolution, styles, show_progress=False)


def test_coevo_annotated(load_alignment):
    """coevolution on alignment with annotated seqs should add to heatmap plot"""
    aln = load_alignment(True, False)
    aln = aln[:30]
    coevo = aln.coevolution(show_progress=False, drawable="heatmap")
    drawable = coevo.drawable
    assert isinstance(drawable, AnnotatedDrawable)
    assert isinstance(drawable.left_track, Drawable)
    assert isinstance(drawable.bottom_track, Drawable)


def test_information_plot(load_alignment):
    """infoprmation plot makes a drawable"""
    aln = load_alignment(False, False)
    drawable = aln.information_plot()
    assert isinstance(drawable, Drawable)
    check_drawable_attrs(drawable.figure, "scatter")


def test_get_drawable():
    """sliced alignment with features returns a drawable"""
    aln = new_alignment.make_aligned_seqs(
        dict(a="AAACGGTTT", b="CAA--GTAA"), moltype="dna"
    )
    db = GffAnnotationDb()
    db.add_feature(seqid="b", biotype="domain", name="1", spans=[(1, 5)])
    db.add_feature(seqid="b", biotype="variation", name="1", spans=[(1, 5)])
    db.add_feature(seqid="b", biotype="gene", name="1", spans=[(1, 5)])
    db.add_feature(seqid="b", biotype="gene", name="1", spans=[(5, 1)])
    aln.annotation_db = db
    got = aln.get_drawable()
    assert got is not None


def test_seqlogo():
    """exercise producing a seq logo"""
    data = {
        "seq1": "CAGGTCGACCTCGGC---------CACGAC",
        "seq2": "CAGATCGACCTCGGC---------CACGAC",
        "seq3": "CAGATCGACCTCGGT---------CACGAT",
        "seq4": "CAGATCGACCTCGGCGAACACGGCCATGAT",
        "seq5": "CCGATCGACATGGGC---------CACGAT",
        "seq6": "GCC---------------------------",
    }
    aln = new_alignment.make_aligned_seqs(data, moltype="dna")
    logo = aln.seqlogo()
    # using wrap argument
    logo = aln.seqlogo(wrap=20)

    # for a sliced alignment
    aln = aln[:20]
    logo = aln.seqlogo()


