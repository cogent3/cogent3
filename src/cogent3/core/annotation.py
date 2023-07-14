import copy
import json

from typing import Iterable, Optional

from cogent3._version import __version__
from cogent3.util import warning as c3warn
from cogent3.util.misc import get_object_provenance

from .location import Map


class _Serialisable:
    def to_rich_dict(self):
        """returns {'name': name, 'seq': sequence, 'moltype': moltype.label}"""
        data = copy.deepcopy(self._serialisable)
        # the first constructor argument will be the instance recreating
        # so we pop out the two possible keys
        data.pop("parent", None)
        data.pop("seq", None)
        if "original" in data:
            data.pop("original")
        # convert the map to coordinates
        data["map"] = data.pop("map").to_rich_dict()
        data = dict(annotation_construction=data)
        data["type"] = get_object_provenance(self)
        data["version"] = __version__

        try:
            annotations = [a.to_rich_dict() for a in self.annotations]
            if annotations:
                data["annotations"] = annotations
        except AttributeError:
            pass

        return data

    def to_json(self):
        return json.dumps(self.to_rich_dict())


# todo gah implement __repr__ and __str__ methods
# todo gah write docstrings!
class Feature:
    """new style annotation, created on demand"""

    __slots__ = (
        "parent",
        "seqid",
        "map",
        "biotype",
        "name",
        "_serialisable",
    )

    # todo gah implement a __new__ to trap args for serialisation purposes?
    def __init__(self, *, parent, seqid: str, map: Map, biotype: str, name: str):
        # _serialisable is used for creating derivative instances
        d = locals()
        exclude = ("self", "__class__", "kw")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}
        self.biotype = biotype
        self.name = name
        self.parent = parent
        self.seqid = seqid
        assert map.parent_length == len(parent), (map, len(parent))
        self.map = map

    def get_slice(self, complete: bool = False, allow_gaps: bool = False):
        """
        The corresponding sequence fragment.

        Parameters
        ----------
        complete
            if feature not complete on parent,causes an exception to be
            raised. If False, gaps are removed.
        allow_gaps
            if on an alignment, includes the gap positions

        Returns
        -------
        a slice of self.parent

        Notes
        -----
        If 'complete' is true and the full length of this feature is not
        present in the sequence, then this method will fail.
        """
        # todo gah set allow_gaps=True as the default
        map = self.map
        if not (complete or map.complete):
            map = map.without_gaps()
        if not allow_gaps:
            return self.parent[map]
        return self.parent[map.start : map.end]

    def without_lost_spans(self):
        """Keeps only the parts which are actually present in the underlying sequence"""
        if self.map.complete:
            return self
        keep = self.map.nongap()
        kwargs = {**self._serialisable, **dict(map=self.map[keep])}
        return self.__class__(**kwargs)

    def as_one_span(self):
        """returns a feature that preserves any gaps"""
        kwargs = {
            **self._serialisable,
            **dict(map=self.map.get_covering_span(), name=f"one-span {self.name}"),
        }
        return self.__class__(**kwargs)

    def shadow(self):
        """returns new instance corresponding to disjoint of self coordinates"""
        kwargs = {
            **self._serialisable,
            **{
                "map": self.map.shadow(),
                "biotype": f"not {self.biotype}",
                "name": f"not {self.name}",
            },
        }
        return self.__class__(**kwargs)

    @c3warn.deprecated_callable(
        "2023.7", reason="new method", new="<instance>.shadow()"
    )
    def get_shadow(self):  # pragma: no cover
        return self.shadow()

    def __len__(self):
        return len(self.map)

    def __repr__(self):
        name = self.__class__.__name__
        txt = ", ".join(
            f"{attr}={getattr(self, attr)!r}"
            for attr in ("seqid", "biotype", "name", "map", "parent")
        )
        return f"{name}({txt})"

    def remapped_to(self, grandparent, gmap):
        seqid = grandparent.name or f"from {self.seqid!r}"
        kwargs = {
            **self._serialisable,
            **{"map": gmap[self.map], "parent": grandparent, "seqid": seqid},
        }
        return self.__class__(**kwargs)

    def get_coordinates(self):
        """returns sequence coordinates of this Feature as
        [(start1, end1), ...]"""
        return self.map.get_coordinates()

    def get_children(self, biotype: Optional[str] = None, **kwargs):
        """generator returns sub-features of self optionally matching biotype"""
        make_feature = self.parent.make_feature
        db = self.parent.annotation_db
        for child in db.get_feature_children(
            biotype=biotype,
            name=self.name,
            start=self.map.start,
            end=self.map.end,
            **kwargs,
        ):
            yield make_feature(feature=child)

    def get_parent(self, **kwargs):
        """generator returns parent features of self optionally matching biotype"""
        make_feature = self.parent.make_feature
        db = self.parent.annotation_db
        for child in db.get_feature_parent(
            name=self.name, start=self.map.start, end=self.map.end, **kwargs
        ):
            yield make_feature(feature=child)

    def union(self, features: Iterable):
        """return as a single Feature

        Notes
        -----
        Overlapping spans are merged
        """
        combined = self.map.spans[:]
        feat_names = [self.name] if self.name else set()
        biotypes = {self.biotype} if self.biotype else set()
        seqids = {self.seqid} if self.seqid else set()
        for feature in features:
            if feature.parent is not self.parent:
                raise ValueError(f"cannot merge annotations from different objects")

            combined.extend(feature.map.spans)
            if feature.name:
                feat_names.append(feature.name)
            if feature.seqid:
                seqids.add(feature.seqid)
            if feature.biotype:
                biotypes.add(feature.biotype)
        name = ", ".join(feat_names)
        map = Map(spans=combined, parent_length=len(self.parent))
        map = map.covered()  # No overlaps
        # the covered method drops reversed status so we need to
        # resurrect that, but noting we've not checked consistency
        # across the features
        if self.map.reverse != map.reverse:
            map = map.reversed()
        seqid = ", ".join(seqids) if seqids else None
        biotype = ", ".join(biotypes)
        kwargs = {
            **self._serialisable,
            **{"map": map, "seqid": seqid, "biotype": biotype, "name": name},
        }

        return self.__class__(**kwargs)

    def get_drawable(self):
        """returns plotly trace"""
        from cogent3.draw.drawable import make_shape

        return make_shape(type_=self)

    def to_dict(self):
        """returns"""
        result = {
            **self._serialisable,
            **dict(
                spans=self.map.get_coordinates(),
                strand="-" if self.map.reverse else "+",
            ),
        }
        for key in ("map", "parent"):
            result.pop(key, None)
        return result
