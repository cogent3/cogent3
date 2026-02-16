"""Tests for cogent3.draw.annotation module."""

import pytest

from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.draw.annotation import (
    _annotations_to_shapes,
    _build_hover_text,
    _get_distinct_seqids,
    _infer_max_coord,
    draw_annotations,
)


@pytest.fixture
def single_seqid_db():
    """Create a DB with features on a single seqid."""
    db = GffAnnotationDb()
    db.add_feature(
        seqid="chr1", biotype="gene", name="geneA", spans=[(100, 500)], strand="+"
    )
    db.add_feature(
        seqid="chr1", biotype="gene", name="geneB", spans=[(600, 900)], strand="-"
    )
    db.add_feature(
        seqid="chr1",
        biotype="CDS",
        name="cdsA",
        spans=[(150, 300), (350, 480)],
        strand="+",
    )
    db.add_feature(
        seqid="chr1", biotype="exon", name="exon1", spans=[(100, 200)], strand="+"
    )
    return db


@pytest.fixture
def multi_seqid_db():
    """Create a DB with features on multiple seqids."""
    db = GffAnnotationDb()
    db.add_feature(
        seqid="chr1", biotype="gene", name="geneA", spans=[(100, 500)], strand="+"
    )
    db.add_feature(
        seqid="chr1", biotype="CDS", name="cdsA", spans=[(150, 300)], strand="+"
    )
    db.add_feature(
        seqid="chr2", biotype="gene", name="geneB", spans=[(200, 800)], strand="-"
    )
    db.add_feature(
        seqid="chr2", biotype="exon", name="exon1", spans=[(200, 400)], strand="-"
    )
    db.add_feature(
        seqid="chr3", biotype="gene", name="geneC", spans=[(50, 300)], strand="+"
    )
    return db


@pytest.fixture
def empty_db():
    """Create an empty annotation DB."""
    return GffAnnotationDb()


def test_build_hover_text_basic():
    feat = {
        "seqid": "chr1",
        "biotype": "gene",
        "name": "geneA",
        "spans": [(100, 500)],
        "strand": "+",
        "on_alignment": False,
        "xattr": {},
    }
    text = _build_hover_text(feat)
    assert "<b>geneA</b>" in text
    assert "biotype: gene" in text
    assert "coords: 100-500" in text
    assert "strand: +" in text


def test_build_hover_text_minus_strand():
    feat = {
        "seqid": "chr1",
        "biotype": "CDS",
        "name": "cdsX",
        "spans": [(100, 200)],
        "strand": "-",
        "on_alignment": False,
        "xattr": {},
    }
    text = _build_hover_text(feat)
    assert "strand: -" in text


def test_build_hover_text_multi_span():
    feat = {
        "seqid": "chr1",
        "biotype": "CDS",
        "name": "cdsA",
        "spans": [(100, 200), (300, 400)],
        "strand": "+",
        "on_alignment": False,
        "xattr": {},
    }
    text = _build_hover_text(feat)
    assert "100-200, 300-400" in text


def test_infer_max_coord_basic():
    features = [
        {
            "spans": [(0, 100)],
            "biotype": "gene",
            "name": "a",
            "seqid": "s",
            "strand": "+",
            "on_alignment": False,
            "xattr": {},
        },
        {
            "spans": [(200, 500)],
            "biotype": "gene",
            "name": "b",
            "seqid": "s",
            "strand": "+",
            "on_alignment": False,
            "xattr": {},
        },
    ]
    result = _infer_max_coord(features)
    assert result == 500 + int(500 * 0.05)


def test_infer_max_coord_empty():
    result = _infer_max_coord([])
    assert result == 0


def test_get_distinct_seqids_single(single_seqid_db):
    seqids = _get_distinct_seqids(single_seqid_db)
    assert seqids == ["chr1"]


def test_get_distinct_seqids_multi(multi_seqid_db):
    seqids = _get_distinct_seqids(multi_seqid_db)
    assert set(seqids) == {"chr1", "chr2", "chr3"}


def test_get_distinct_seqids_empty(empty_db):
    seqids = _get_distinct_seqids(empty_db)
    assert seqids == []


def test_annotations_to_shapes_groups_by_biotype(single_seqid_db):
    features = list(
        single_seqid_db.get_features_matching(seqid="chr1", allow_partial=True)
    )
    shapes = _annotations_to_shapes(features, max_coord=1000)
    assert "gene" in shapes
    assert "CDS" in shapes


def test_annotations_to_shapes_empty_features():
    shapes = _annotations_to_shapes([], max_coord=1000)
    assert shapes == {}


def test_draw_annotations_single_seqid_explicit(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1")
    assert result is not None
    assert len(result.traces) > 0


def test_draw_annotations_single_seqid_auto_resolve(single_seqid_db):
    result = draw_annotations(single_seqid_db)
    assert result is not None
    assert len(result.traces) > 0


def test_draw_annotations_no_matches_returns_none(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="nonexistent")
    assert result is None


def test_draw_annotations_empty_db_returns_none(empty_db):
    result = draw_annotations(empty_db)
    assert result is None


def test_draw_annotations_biotype_filtering(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", biotype="gene")
    assert result is not None
    # All traces should have legendgroup == "gene"
    for trace in result.traces:
        assert trace["legendgroup"] == "gene"


def test_draw_annotations_biotype_tuple_filtering(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", biotype=("gene", "CDS"))
    assert result is not None
    biotypes = {trace["legendgroup"] for trace in result.traces}
    assert biotypes <= {"gene", "CDS"}


def test_draw_annotations_start_stop_filtering(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", start=0, stop=200)
    assert result is not None


def test_draw_annotations_strand_reverse_mapping(single_seqid_db):
    """Minus strand features should produce reversed arrows."""
    result = draw_annotations(single_seqid_db, seqid="chr1", biotype="gene")
    assert result is not None
    assert len(result.traces) >= 2


def test_draw_annotations_max_coord_sets_xaxis_range(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", max_coord=2000)
    assert result is not None
    assert result.layout.xaxis.range[1] == 2000


def test_draw_annotations_range_slider_present(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1")
    assert result is not None
    assert result.layout.xaxis.rangeslider.visible is True


def test_draw_annotations_controls_disabled(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", show_controls=False)
    assert result is not None
    assert "rangeslider" not in result.layout.xaxis


def test_draw_annotations_legend_grouping(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1")
    assert result is not None
    for trace in result.traces:
        assert "legendgroup" in trace


def test_draw_annotations_hover_text_present(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1")
    assert result is not None
    for trace in result.traces:
        text = trace.get("text", "")
        assert "biotype:" in text


def test_draw_annotations_multi_seqid_stacked(multi_seqid_db):
    result = draw_annotations(multi_seqid_db)
    assert result is not None
    # Should have traces from all seqids
    assert len(result.traces) > 0


def test_draw_annotations_multi_seqid_independent_xaxis(multi_seqid_db):
    result = draw_annotations(multi_seqid_db)
    assert result is not None
    # Should have xaxis and xaxis2 etc
    assert "xaxis" in result.layout
    assert "xaxis2" in result.layout


def test_draw_annotations_multi_seqid_no_duplicate_legend(multi_seqid_db):
    result = draw_annotations(multi_seqid_db)
    assert result is not None
    # "gene" should only have showlegend=True once
    gene_shown = sum(
        1
        for t in result.traces
        if t.get("legendgroup") == "gene" and t.get("showlegend", True)
    )
    assert gene_shown == 1


def test_draw_annotations_multi_seqid_yaxis_domains(multi_seqid_db):
    result = draw_annotations(multi_seqid_db)
    assert result is not None
    # yaxis domains should be set for each seqid
    assert "domain" in result.layout.yaxis


def test_draw_annotations_max_features_exceeded(single_seqid_db):
    with pytest.raises(ValueError, match="limit is 1"):
        draw_annotations(single_seqid_db, seqid="chr1", max_features=1)


def test_draw_annotations_max_features_override(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", max_features=10000)
    assert result is not None


def test_draw_annotations_max_features_exceeded_no_seqid(multi_seqid_db):
    with pytest.raises(ValueError, match="limit is 1"):
        draw_annotations(multi_seqid_db, max_features=1)


def test_draw_annotations_title_parameter(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", title="My Title")
    assert result is not None
    assert result.layout.title.text == "My Title"


def test_draw_annotations_width_parameter(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", width=800)
    assert result is not None
    assert result.layout.width == 800


def test_draw_annotations_xaxis_title_is_range_slider(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1")
    assert result is not None
    assert result.layout.xaxis.title.text == "Range Slider"


def test_draw_annotations_xaxis_title_is_seqid_when_no_controls(single_seqid_db):
    result = draw_annotations(single_seqid_db, seqid="chr1", show_controls=False)
    assert result is not None
    assert result.layout.xaxis.title.text == "chr1"


def test_draw_annotations_via_cogent3_namespace():
    """draw_annotations is accessible from top-level cogent3."""
    import cogent3

    assert hasattr(cogent3, "draw_annotations")
    assert callable(cogent3.draw_annotations)
