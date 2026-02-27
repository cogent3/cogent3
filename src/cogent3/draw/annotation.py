"""Visualise genomic annotations directly from an AnnotationDb."""

from __future__ import annotations

import warnings
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


def _infer_coord_range(
    features: list[FeatureDataType],
    padding_frac: float = 0.05,
) -> tuple[int, int]:
    """Determine x-axis range from feature spans with padding."""
    min_val = float("inf")
    max_val = 0
    for feat in features:
        for span in feat["spans"]:
            min_val = min(min_val, span[0], span[1])
            max_val = max(max_val, span[0], span[1])
    if min_val == float("inf"):
        return 0, 0
    span_width = max_val - min_val
    padding = int(span_width * padding_frac) if span_width > 0 else 0
    return int(max(0, min_val - padding)), max_val + padding


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
        min_coord, max_coord = _infer_coord_range(features)
        axis_range = [min_coord, max_coord]
        shape_max_coord = max_coord - min_coord
    else:
        axis_range = [0, max_coord]
        shape_max_coord = max_coord

    drawables = _annotations_to_shapes(features, shape_max_coord)
    if not drawables:
        return [], {}, {}, 0.0

    annotes, top = stack_shapes(drawables)
    all_traces = [t.as_trace() for t in annotes]

    xaxis: dict[str, Any] = {
        "range": axis_range,
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

    if max_coord is None:
        inferred_min, inferred_max = _infer_coord_range(features)
        effective_max = inferred_max - inferred_min
    else:
        effective_max = max_coord
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
    height_per_seqid = 300

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

        domain = get_domain(n_seqids, idx, is_y=True, space=0.12)
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


def _find_anchor_features(
    db: AnnotationDbABC,
    *,
    name: str,
    biotype: str | tuple[str, ...] | list[str] | set[str] | None = None,
    seqid: str | None = None,
    max_seqids: int | None = None,
) -> dict[str, FeatureDataType]:
    """Find the first feature matching ``name`` on each seqid.

    Parameters
    ----------
    db
        An annotation database instance.
    name
        Feature name to search for.
    biotype
        Optional biotype filter.
    seqid
        If provided, restrict search to this seqid only.
    max_seqids
        If provided, stop collecting after this many distinct seqids.

    Returns
    -------
    dict mapping seqid to the first matching feature record.
    """
    kwargs: dict[str, Any] = {"name": name, "allow_partial": True}
    if biotype is not None:
        kwargs["biotype"] = biotype
    if seqid is not None:
        kwargs["seqid"] = seqid
    anchors: dict[str, FeatureDataType] = {}
    for feat in db.get_features_matching(**kwargs):
        sid = feat["seqid"]
        if sid not in anchors:
            anchors[sid] = feat
            if max_seqids is not None and len(anchors) >= max_seqids:
                break
    return anchors


def _get_span_extent(feature: FeatureDataType) -> tuple[int, int]:
    """Return (min_start, max_end) across all spans in a feature."""
    if not (spans := feature.get("spans")):
        msg = "feature['spans'] is empty, cannot compute span extent"
        raise ValueError(msg)

    if hasattr(spans, "tolist"):
        spans = spans.tolist()
    all_starts = [s[0] for s in spans]
    all_ends = [s[1] for s in spans]
    return min(all_starts), max(all_ends)


def _draw_centered(
    db: AnnotationDbABC,
    *,
    center_on: str,
    seqid: str | None,
    biotype: str | tuple[str, ...] | list[str] | set[str] | None,
    flank: int,
    max_seqids: int | None,
    max_features: int,
    width: float,
    title: str | None,
    show_controls: bool,
) -> Drawable | None:
    """Build a drawable centered on a named feature across seqids."""
    anchors = _find_anchor_features(
        db,
        name=center_on,
        biotype=biotype,
        seqid=seqid,
        max_seqids=max_seqids,
    )
    if not anchors:
        return None

    selected_seqids = sorted(anchors)

    # Compute uniform window width from the maximum anchor span
    span_extents = {sid: _get_span_extent(anchors[sid]) for sid in selected_seqids}
    max_span = max(end - start for (start, end) in span_extents.values())
    half = (max_span + 2 * flank) // 2

    by_seqid: dict[str, list[FeatureDataType]] = {}
    for sid in selected_seqids:
        start, end = span_extents[sid]
        mid = (start + end) // 2
        w_start = max(0, mid - half)
        w_stop = mid + half
        if feats := _query_features(
            db,
            seqid=sid,
            biotype=biotype,
            start=w_start,
            stop=w_stop,
        ):
            by_seqid[sid] = feats

    if not by_seqid:
        return None

    total = sum(len(v) for v in by_seqid.values())
    _check_max_features(total, max_features)

    default_title = title or f"Centered on: {center_on}"
    seqids = sorted(by_seqid)
    if len(seqids) == 1:
        return _make_single_seqid_drawable(
            by_seqid[seqids[0]],
            seqid=seqids[0],
            max_coord=None,
            width=width,
            title=default_title,
            show_controls=show_controls,
        )

    return _make_multi_seqid_drawable(
        by_seqid,
        max_coord=None,
        width=width,
        title=default_title,
    )


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
    center_on: str | None = None,
    flank: int = 5000,
    max_seqids: int | None = None,
) -> Drawable | None:
    """Visualise annotations from an AnnotationDb

    Parameters
    ----------
    db
        An annotation database instance.
    seqid
        Sequence identifier to display. If None, features for all seqids
        are shown as stacked subplots. When used with ``center_on``,
        restricts the anchor search to this seqid.
    biotype
        Feature type(s) to include. If None, all biotypes are shown.
    name
        Feature name filter. Ignored when ``center_on`` is set.
    start
        Start coordinate filter. Ignored when ``center_on`` is set.
    stop
        Stop coordinate filter. Ignored when ``center_on`` is set.
    strand
        Strand filter. Ignored when ``center_on`` is set.
    max_coord
        Maximum x-axis coordinate. Inferred from features if not provided.
        Ignored when ``center_on`` is set.
    max_features
        Maximum number of features to render. Raise ValueError if exceeded.
    width
        Figure width in pixels.
    title
        Figure title.
    show_controls
        If True (default), display a range slider for navigating coordinates.
        Set to False for clean static image export.
    center_on
        Feature name to center on. When provided, finds features with this
        name across all seqids and displays a window around each match.
        Only seqids containing the feature are shown.
    flank
        Number of bases on each side of the anchor feature's midpoint.
        Only used when ``center_on`` is set.
    max_seqids
        Maximum number of seqids to display. Applied after filtering to
        those containing the anchor. Sorted alphabetically, first N taken.

    Returns
    -------
    Drawable or None
        A Drawable instance, or None if no features match.
    """
    if center_on is not None:
        if name is not None:
            warnings.warn(
                "name is redundant when center_on is set and thus ignored",
                UserWarning,
                stacklevel=2,
            )
        if start is not None or stop is not None or strand is not None:
            warnings.warn(
                "start/stop/strand are ignored when center_on is set",
                UserWarning,
                stacklevel=2,
            )
        if max_coord is not None:
            warnings.warn(
                "max_coord is ignored when center_on is set",
                UserWarning,
                stacklevel=2,
            )

        return _draw_centered(
            db,
            center_on=center_on,
            seqid=seqid,
            biotype=biotype,
            flank=flank,
            max_seqids=max_seqids,
            max_features=max_features,
            width=width,
            title=title,
            show_controls=show_controls,
        )

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
