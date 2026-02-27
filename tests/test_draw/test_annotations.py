"""Tests for cogent3.draw.annotation module."""

import warnings

import pytest

from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.draw.annotation import (
    _annotations_to_shapes,
    _build_hover_text,
    _find_anchor_features,
    _get_distinct_seqids,
    _get_span_extent,
    _infer_coord_range,
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


def test_infer_coord_range_basic():
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
    lo, hi = _infer_coord_range(features)
    padding = int(500 * 0.05)
    assert lo == 0
    assert hi == 500 + padding


def test_infer_coord_range_offset():
    """Features far from origin should produce min > 0."""
    features = [
        {
            "spans": [(1000, 2000)],
            "biotype": "gene",
            "name": "a",
            "seqid": "s",
            "strand": "+",
            "on_alignment": False,
            "xattr": {},
        },
    ]
    lo, hi = _infer_coord_range(features)
    padding = int(1000 * 0.05)
    assert lo == 1000 - padding
    assert hi == 2000 + padding


def test_infer_coord_range_empty():
    lo, hi = _infer_coord_range([])
    assert (lo, hi) == (0, 0)


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
    assert result.layout.xaxis.range[0] == 0
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


# --- Fixtures and tests for center_on functionality ---


@pytest.fixture
def shared_feature_db():
    """Create a DB with a shared feature across multiple seqids."""
    db = GffAnnotationDb()
    # Shared gene "dnaA" on three seqids at different positions
    db.add_feature(
        seqid="genome1",
        biotype="gene",
        name="dnaA",
        spans=[(10000, 11000)],
        strand="+",
    )
    db.add_feature(
        seqid="genome1",
        biotype="gene",
        name="otherGene",
        spans=[(12000, 13000)],
        strand="+",
    )
    db.add_feature(
        seqid="genome2",
        biotype="gene",
        name="dnaA",
        spans=[(20000, 21200)],
        strand="-",
    )
    db.add_feature(
        seqid="genome2",
        biotype="CDS",
        name="someCDS",
        spans=[(19000, 19500)],
        strand="+",
    )
    db.add_feature(
        seqid="genome3",
        biotype="gene",
        name="dnaA",
        spans=[(5000, 5800)],
        strand="+",
    )
    db.add_feature(
        seqid="genome3",
        biotype="exon",
        name="exon1",
        spans=[(5000, 5400)],
        strand="+",
    )
    # genome4 does NOT have dnaA
    db.add_feature(
        seqid="genome4",
        biotype="gene",
        name="unrelated",
        spans=[(1000, 2000)],
        strand="+",
    )
    return db


def test_find_anchor_features(shared_feature_db):
    anchors = _find_anchor_features(shared_feature_db, name="dnaA")
    assert set(anchors.keys()) == {"genome1", "genome2", "genome3"}
    for sid, feat in anchors.items():
        assert feat["name"] == "dnaA"
        assert feat["seqid"] == sid


def test_find_anchor_features_with_biotype(shared_feature_db):
    anchors = _find_anchor_features(
        shared_feature_db,
        name="dnaA",
        biotype="gene",
    )
    assert set(anchors.keys()) == {"genome1", "genome2", "genome3"}


def test_find_anchor_features_no_match(shared_feature_db):
    anchors = _find_anchor_features(shared_feature_db, name="nonexistent")
    assert anchors == {}


def test_center_on_produces_stacked_subplots(shared_feature_db):
    result = draw_annotations(shared_feature_db, center_on="dnaA", biotype="gene")
    assert result is not None
    assert len(result.traces) > 0
    # Should have multiple x-axes for multiple seqids
    assert "xaxis2" in result.layout


def test_center_on_only_includes_matching_seqids(shared_feature_db):
    result = draw_annotations(shared_feature_db, center_on="dnaA", biotype="gene")
    assert result is not None
    # genome4 should not appear (no dnaA)
    # Check x-axis titles contain only genome1, genome2, genome3
    axis_titles = set()
    for key in result.layout:
        if key.startswith("xaxis"):
            title_text = result.layout[key].get("title", {})
            if isinstance(title_text, dict):
                axis_titles.add(title_text.get("text", ""))
            elif hasattr(title_text, "text"):
                axis_titles.add(title_text.text)
    assert "genome4" not in axis_titles


def test_center_on_flank_controls_window(shared_feature_db):
    result_small = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        biotype="gene",
        flank=100,
    )
    result_large = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        biotype="gene",
        flank=50000,
    )
    assert result_small is not None
    assert result_large is not None
    # Larger flank should produce wider x-axis range
    small_range = result_small.layout.xaxis.range
    large_range = result_large.layout.xaxis.range
    small_width = small_range[1] - small_range[0]
    large_width = large_range[1] - large_range[0]
    assert large_width > small_width
    # x-axis range should NOT start at 0 (features are far from origin)
    assert small_range[0] > 0


def test_center_on_xaxis_range_matches_window(shared_feature_db):
    """x-axis range should be inferred from features in the window."""
    flank = 2000
    result = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        biotype="gene",
        flank=flank,
        max_seqids=1,
    )
    assert result is not None
    xrange = result.layout.xaxis.range
    # max_seqids=1 selects genome1 only; dnaA at (10000, 11000)
    # Window queries [8000, 13000], finding dnaA (10000-11000) and
    # otherGene (12000-13000). _infer_coord_range computes bounds from
    # feature spans with padding.
    # The range should NOT start at 0 (features are far from origin)
    assert xrange[0] > 0
    # The range should cover the anchor feature
    assert xrange[0] <= 10000
    assert xrange[1] >= 11000


def test_center_on_max_seqids_limits_output(shared_feature_db):
    result = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        biotype="gene",
        max_seqids=2,
    )
    assert result is not None
    # Count distinct x-axes
    n_axes = sum(1 for k in result.layout if k.startswith("xaxis"))
    assert n_axes == 2


def test_center_on_max_seqids_one(shared_feature_db):
    result = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        biotype="gene",
        max_seqids=1,
    )
    assert result is not None
    # Single seqid should produce single subplot with range slider
    n_axes = sum(1 for k in result.layout if k.startswith("xaxis"))
    assert n_axes == 1


def test_center_on_returns_none_when_not_found(shared_feature_db):
    result = draw_annotations(shared_feature_db, center_on="nonexistent")
    assert result is None


def test_center_on_biotype_filters_displayed_features(shared_feature_db):
    result = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        biotype="gene",
    )
    assert result is not None
    # All traces should be gene biotype
    for trace in result.traces:
        assert trace.get("legendgroup") == "gene"


def test_center_on_honors_seqid(shared_feature_db):
    """When seqid is set with center_on, only that seqid is shown."""
    result = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        seqid="genome1",
    )
    assert result is not None
    # Single seqid should produce one subplot
    n_axes = sum(1 for k in result.layout if k.startswith("xaxis"))
    assert n_axes == 1
    # The x-axis title should be genome1
    assert result.layout.xaxis.title.text in ("genome1", "Range Slider")


def test_center_on_warns_when_max_coord_set(shared_feature_db):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        draw_annotations(
            shared_feature_db,
            center_on="dnaA",
            max_coord=50000,
        )
        assert any("max_coord is ignored" in str(warning.message) for warning in w)


def test_center_on_warns_when_start_stop_set(shared_feature_db):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        draw_annotations(
            shared_feature_db,
            center_on="dnaA",
            start=0,
            stop=100,
        )
        assert any("start/stop are ignored" in str(warning.message) for warning in w)


def test_center_on_default_title(shared_feature_db):
    result = draw_annotations(shared_feature_db, center_on="dnaA", biotype="gene")
    assert result is not None
    assert "dnaA" in result.layout.title.text


@pytest.mark.parametrize("spans", [[], None], ids=["empty_list", "missing_key"])
def test_get_span_extent_no_spans(spans):
    """_get_span_extent raises ValueError when feature has no spans."""
    feature = {"spans": spans} if spans is not None else {}
    with pytest.raises(ValueError, match="empty"):
        _get_span_extent(feature)


def test_center_on_custom_title(shared_feature_db):
    result = draw_annotations(
        shared_feature_db,
        center_on="dnaA",
        biotype="gene",
        title="Custom",
    )
    assert result is not None
    assert result.layout.title.text == "Custom"
