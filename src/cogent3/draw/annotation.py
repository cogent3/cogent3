"""Visualise genomic annotations directly from an AnnotationDb."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from cogent3.core.location import Strand
from cogent3.draw.drawable import Drawable, Shape, get_domain, make_shape, stack_shapes

if TYPE_CHECKING:
    from cogent3.core.annotation_db import AnnotationDbABC, FeatureDataType


def _build_hover_text(feature: FeatureDataType) -> str:
    """Build HTML hover string from a feature record."""
    spans = feature["spans"]
    if hasattr(spans, "tolist"):
        spans = spans.tolist()
    coords_str = ", ".join(f"{s}-{e}" for s, e in spans)
    strand_val = Strand.from_value(feature["strand"])
    strand_map = {Strand.PLUS: "+", Strand.MINUS: "-"}
    strand_str = strand_map.get(strand_val, ".")
    return (
        f"<b>{feature['name']}</b><br>"
        f"biotype: {feature['biotype']}<br>"
        f"coords: {coords_str}<br>"
        f"strand: {strand_str}"
    )


def _infer_max_coord(
    features: list[FeatureDataType],
    padding_frac: float = 0.05,
) -> int:
    """Determine x-axis range from feature spans with padding."""
    max_val = 0
    for feat in features:
        for span in feat["spans"]:
            max_val = max(max_val, span[0], span[1])
    padding = int(max_val * padding_frac)
    return max_val + padding


def _get_distinct_seqids(db: AnnotationDbABC) -> list[str]:
    """Query DB for available seqids."""
    table = db.count_distinct(seqid=True)
    if table is None:
        return []
    return [str(v) for v in table.columns["seqid"]]


def _annotations_to_shapes(
    features: list[FeatureDataType],
    max_coord: int,
) -> dict[str, list[Shape]]:
    """Convert feature records to Shape objects grouped by biotype."""
    result: dict[str, list[Shape]] = {}
    for feat in features:
        biotype = feat["biotype"]
        hover_text = _build_hover_text(feat)
        spans = feat["spans"]
        if hasattr(spans, "tolist"):
            spans = spans.tolist()
        spans = [tuple(s) for s in spans]

        strand_val = Strand.from_value(feat["strand"])
        reverse = strand_val == Strand.MINUS

        shape = make_shape(
            type_=biotype,
            name=hover_text,
            coords=spans,
            reverse=reverse,
            parent_length=max_coord,
        )
        if shape is not None:
            result.setdefault(biotype, []).append(shape)
    return result


def _draw_single_seqid(
    features: list[FeatureDataType],
    *,
    seqid: str,
    max_coord: int | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any], float]:
    """Render features for one seqid.

    Returns (traces, xaxis_config, yaxis_config, top_y).
    """
    if max_coord is None:
        max_coord = _infer_max_coord(features)

    drawables = _annotations_to_shapes(features, max_coord)
    if not drawables:
        return [], {}, {}, 0.0

    annotes, top = stack_shapes(drawables)
    all_traces = [t.as_trace() for t in annotes]

    xaxis: dict[str, Any] = {
        "range": [0, max_coord],
        "zeroline": False,
        "showline": True,
        "title": {"text": seqid},
    }
    yaxis: dict[str, Any] = {
        "range": [0, top],
        "visible": False,
        "zeroline": True,
        "showline": True,
    }

    return all_traces, xaxis, yaxis, top


def _query_features(
    db: AnnotationDbABC,
    *,
    seqid: str | None = None,
    biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
    name: str | None = None,
    start: int | None = None,
    stop: int | None = None,
    strand: str | None = None,
) -> list[FeatureDataType]:
    """Query features from DB, building kwargs from non-None params."""
    kwargs: dict[str, Any] = {}
    if seqid is not None:
        kwargs["seqid"] = seqid
    if biotype is not None:
        kwargs["biotype"] = biotype
    if name is not None:
        kwargs["name"] = name
    if start is not None:
        kwargs["start"] = start
    if stop is not None:
        kwargs["stop"] = stop
    if strand is not None:
        kwargs["strand"] = strand
    return list(db.get_features_matching(allow_partial=True, **kwargs))


def _check_max_features(n_features: int, max_features: int) -> None:
    """Raise ValueError if feature count exceeds limit."""
    if n_features > max_features:
        msg = (
            f"Query returned {n_features} features (limit is {max_features}). "
            f"Narrow your query with biotype, start/stop, or name filters, "
            f"or increase max_features."
        )
        raise ValueError(msg)


def _make_single_seqid_drawable(
    features: list[FeatureDataType],
    *,
    seqid: str,
    max_coord: int | None,
    width: float,
    title: str | None,
    show_controls: bool = True,
) -> Drawable | None:
    """Build a Drawable for a single seqid with optional range slider."""
    traces, xaxis, yaxis, top = _draw_single_seqid(
        features,
        seqid=seqid,
        max_coord=max_coord,
    )
    if not traces:
        return None

    effective_max = max_coord or _infer_max_coord(features)
    height = max((top / max(1, effective_max)) * width, 300)
    if show_controls:
        xaxis["rangeslider"] = {"visible": True}
        xaxis["title"] = {"text": "Range Slider"}

    drawer = Drawable(
        title=title or f"Annotations: {seqid}",
        traces=traces,
        width=width,
        height=height,
    )
    drawer.layout.update(xaxis=xaxis, yaxis=yaxis)
    return drawer


def _make_multi_seqid_drawable(
    by_seqid: dict[str, list[FeatureDataType]],
    *,
    max_coord: int | None,
    width: float,
    title: str | None,
) -> Drawable | None:
    """Build a stacked subplot Drawable for multiple seqids."""
    seqids = sorted(by_seqid)
    n_seqids = len(seqids)
    height_per_seqid = 200

    all_traces: list[dict[str, Any]] = []
    layout_updates: dict[str, Any] = {}
    seen_legendgroups: set[str] = set()

    for idx, sid in enumerate(seqids):
        axis_num = idx + 1
        xaxis_key = f"xaxis{axis_num}" if axis_num > 1 else "xaxis"
        yaxis_key = f"yaxis{axis_num}" if axis_num > 1 else "yaxis"
        xref = f"x{axis_num}" if axis_num > 1 else "x"
        yref = f"y{axis_num}" if axis_num > 1 else "y"

        traces, xaxis_cfg, yaxis_cfg, _ = _draw_single_seqid(
            by_seqid[sid],
            seqid=sid,
            max_coord=max_coord,
        )

        domain = get_domain(n_seqids, idx, is_y=True)
        yaxis_cfg["domain"] = list(domain)
        xaxis_cfg["anchor"] = yref

        for trace in traces:
            trace["xaxis"] = xref
            trace["yaxis"] = yref
            lg = trace.get("legendgroup", "")
            if lg in seen_legendgroups:
                trace["showlegend"] = False
            else:
                seen_legendgroups.add(lg)
            all_traces.append(trace)

        layout_updates[xaxis_key] = xaxis_cfg
        layout_updates[yaxis_key] = yaxis_cfg

    if not all_traces:
        return None

    drawer = Drawable(
        title=title or "Annotations",
        traces=all_traces,
        width=width,
        height=height_per_seqid * n_seqids,
    )
    drawer.layout.update(**layout_updates)
    return drawer


def draw_annotations(
    db: AnnotationDbABC,
    *,
    seqid: str | None = None,
    biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
    name: str | None = None,
    start: int | None = None,
    stop: int | None = None,
    strand: str | None = None,
    max_coord: int | None = None,
    max_features: int = 2000,
    width: float = 600,
    title: str | None = None,
    show_controls: bool = True,
) -> Drawable | None:
    """Visualise annotations from an AnnotationDb

    Parameters
    ----------
    db
        An annotation database instance.
    seqid
        Sequence identifier to display. If None, features for all seqids
        are shown as stacked subplots.
    biotype
        Feature type(s) to include. If None, all biotypes are shown.
    name
        Feature name filter.
    start
        Start coordinate filter.
    stop
        Stop coordinate filter.
    strand
        Strand filter.
    max_coord
        Maximum x-axis coordinate. Inferred from features if not provided.
    max_features
        Maximum number of features to render. Raise ValueError if exceeded.
    width
        Figure width in pixels.
    title
        Figure title.
    show_controls
        If True (default), display a range slider for navigating coordinates.
        Set to False for clean static image export.

    Returns
    -------
    Drawable or None
        A Drawable instance, or None if no features match.
    """
    features = _query_features(
        db,
        seqid=seqid,
        biotype=biotype,
        name=name,
        start=start,
        stop=stop,
        strand=strand,
    )
    if not features:
        return None

    _check_max_features(len(features), max_features)

    if seqid is not None:
        return _make_single_seqid_drawable(
            features,
            seqid=seqid,
            max_coord=max_coord,
            width=width,
            title=title,
            show_controls=show_controls,
        )

    # Group features by seqid
    by_seqid: dict[str, list[FeatureDataType]] = {}
    for feat in features:
        by_seqid.setdefault(feat["seqid"], []).append(feat)

    seqids = sorted(by_seqid)

    if len(seqids) == 1:
        return _make_single_seqid_drawable(
            by_seqid[seqids[0]],
            seqid=seqids[0],
            max_coord=max_coord,
            width=width,
            title=title,
            show_controls=show_controls,
        )

    return _make_multi_seqid_drawable(
        by_seqid,
        max_coord=max_coord,
        width=width,
        title=title,
    )
