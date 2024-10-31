from typing import Iterable, Optional

from numpy import array

from .location import FeatureMap


# todo gah write docstrings!
class Feature:
    """new style annotation, created on demand"""

    # we make the object immutable by making public attributes a property
    __slots__ = (
        "_parent",
        "_seqid",
        "_map",
        "_biotype",
        "_name",
        "_serialisable",
        "_id",
        "_strand",
    )

    # todo gah implement a __new__ to trap args for serialisation purposes?
    def __init__(
        self,
        *,
        parent,
        seqid: str,
        map: FeatureMap,
        biotype: str,
        name: str,
        strand: str,
    ):
        # _serialisable is used for creating derivative instances
        d = locals()
        exclude = ("self", "__class__", "kw")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}
        self._parent = parent
        self._seqid = seqid
        assert map.parent_length == len(parent), (map, len(parent))
        self._map = map
        self._biotype = biotype
        self._name = name
        data = [id(self.parent), tuple(self.map.get_coordinates())]
        data.extend((self.seqid, self.biotype, self.name))
        self._id = hash(tuple(data))
        self._strand = strand

    def __eq__(self, other):
        return self._id == other._id

    def __hash__(self):
        """Features can be used in a dictionary!"""
        return self._id

    @property
    def parent(self):
        return self._parent

    @property
    def seqid(self):
        return self._seqid

    @property
    def map(self):
        return self._map

    @property
    def biotype(self):
        return self._biotype

    @property
    def name(self):
        return self._name

    def get_slice(self, complete: bool = False, allow_gaps: bool = False):
        """
        The corresponding sequence fragment.

        Parameters
        ----------
        complete
            if feature not complete on parent, causes an exception to be
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
        fmap = self.map

        if not complete and not fmap.complete:
            fmap = fmap.without_gaps()
        if not allow_gaps:
            result = self.parent[fmap]
            return self._do_seq_slice(result)
        # all slicing now requires start < end
        result = self.parent[fmap.start : fmap.end]
        return self._do_seq_slice(result)

    def _do_seq_slice(self, result):
        if self.reversed:
            result = result.rc()
        if self.map.num_spans > 1:
            # db querying will be incorrect so make sure it can't be done
            result.annotation_db = None
        return result

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
        # grandparent can be either a Sequence or an Alignment
        if not isinstance(gmap, FeatureMap):
            # due to separation of IndelMap and Map, change class
            gmap = gmap.to_feature_map()

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
        offset = getattr(self.parent, "annotation_offset", 0)
        start = self.map.start + offset
        stop = self.map.end + offset

        make_feature = self.parent.make_feature
        db = self.parent.annotation_db
        for record in db.get_feature_children(
            biotype=biotype,
            name=self.name,
            start=start,
            stop=stop,
            **kwargs,
        ):
            record["spans"] = array(record["spans"]) - offset
            feature = make_feature(feature=record)
            # For GenBank, the names can be shared between parent and child,
            # but we don't want to return self. The cleanest way to not
            # return self is to just check.
            if feature == self:
                continue
            yield feature

    def get_parent(self, **kwargs):
        """generator returns parent features of self optionally matching biotype"""
        offset = getattr(self.parent, "annotation_offset", 0)
        start = self.map.start + offset
        stop = self.map.end + offset

        make_feature = self.parent.make_feature
        db = self.parent.annotation_db
        for record in db.get_feature_parent(
            name=self.name,
            start=start,
            stop=stop,
            **kwargs,
        ):
            record["spans"] = array(record["spans"]) - offset
            feature = make_feature(feature=record)
            # For GenBank, the names can be shared between parent and child,
            # but we don't want to return self. The cleanest way to not
            # return self is to just check.
            if feature == self:
                continue
            yield feature

    def union(self, features: Iterable):
        """return as a single Feature

        Notes
        -----
        Overlapping spans are merged
        """
        # spans always on the plus strand, irrespective of whether
        # a feature is reversed
        combined = list(self.map.spans)
        feat_names = [self.name] if self.name else set()
        biotypes = {self.biotype} if self.biotype else set()
        seqids = {self.seqid} if self.seqid else set()

        same_orientation = True
        for feature in features:
            if feature.parent is not self.parent:
                raise ValueError("cannot merge annotations from different objects")

            if same_orientation and feature.reversed != self.reversed:
                same_orientation = False

            combined.extend(feature.map.spans)
            if feature.name:
                feat_names.append(feature.name)
            if feature.seqid:
                seqids.add(feature.seqid)
            if feature.biotype:
                biotypes.add(feature.biotype)

        name = ", ".join(feat_names)
        fmap = FeatureMap(spans=combined, parent_length=len(self.parent))
        fmap = fmap.covered()  # No overlaps
        # the covered method drops reversed status so we need to
        # resurrect that, but noting we've not checked consistency
        # across the features
        strand = self._strand if same_orientation else "+"
        seqid = ", ".join(seqids) if seqids else None
        biotype = ", ".join(biotypes)
        kwargs = {
            **self._serialisable,
            **{
                "map": fmap,
                "seqid": seqid,
                "biotype": biotype,
                "name": name,
                "strand": strand,
            },
        }

        return self.__class__(**kwargs)

    def get_drawable(self):
        """returns plotly trace"""
        from cogent3.draw.drawable import make_shape

        return make_shape(type_=self, parent_length=len(self.parent))

    def to_dict(self):
        """returns"""
        result = {
            **self._serialisable,
            **dict(
                spans=self.map.get_coordinates(),
            ),
        }
        for key in ("map", "parent"):
            result.pop(key, None)
        return result

    @property
    def reversed(self):
        """whether Feature is on the reverse strand relative to bound object"""
        return self._strand == "-"
