import typing
from collections.abc import Iterable

import typing_extensions
from numpy import array

from .location import FeatureMap, SeqCoordTypes, Strand

if typing.TYPE_CHECKING:
    from cogent3.core.new_alignment import Alignment
    from cogent3.core.new_sequence import Sequence
    from cogent3.draw.drawable import Shape

SeqORAlign = typing.Union["Sequence", "Alignment"]


# TODO gah write docstrings!
class Feature:
    """new style annotation, created on demand"""

    # we make the object immutable by making public attributes a property
    __slots__ = (
        "_biotype",
        "_id",
        "_map",
        "_name",
        "_parent",
        "_seqid",
        "_serialisable",
        "_strand",
        "_xattr",
    )

    # TODO gah implement a __new__ to trap args for serialisation purposes?
    def __init__(
        self,
        *,
        parent: SeqORAlign,
        seqid: str,
        map: FeatureMap,
        biotype: str,
        name: str,
        strand: int | str,
        xattr: dict[str, typing.Any] | None = None,
    ) -> None:
        # _serialisable is used for creating derivative instances
        d = locals()
        self._serialisable = {k: v for k, v in d.items() if k not in ("self", "d")}
        self._parent = parent
        self._seqid = seqid
        assert map.parent_length == len(parent), (map, len(parent))
        self._map = map
        self._biotype = biotype
        self._name = name
        data = [id(self.parent), tuple(self.map.get_coordinates())]
        data.extend((self.seqid, self.biotype, self.name))
        self._id = hash(tuple(data))
        self._strand = Strand.from_value(strand)
        self._xattr = xattr

    def __eq__(self, other: object) -> bool:
        return self._id == getattr(other, "_id", None)

    def __hash__(self) -> int:
        """Features can be used in a dictionary!"""
        return self._id

    @property
    def parent(self) -> SeqORAlign:
        """sequence or aligned or alignment"""
        return self._parent

    @property
    def seqid(self) -> str:
        """the sequence id of the parent sequence"""
        return self._seqid

    @property
    def map(self) -> FeatureMap:
        """coordinate map and properties of this feature"""
        return self._map

    @property
    def biotype(self) -> str:
        """type of biological feature"""
        return self._biotype

    @property
    def name(self) -> str:
        """name of the feature"""
        return self._name

    def get_slice(
        self,
        complete: bool = False,
        allow_gaps: bool = False,
        apply_name: bool = True,
    ) -> SeqORAlign:
        """
        The corresponding sequence fragment.

        Parameters
        ----------
        complete
            if feature not complete on parent, causes an exception to be
            raised. If False, gaps are removed.
        allow_gaps
            if on an alignment, includes the gap positions
        apply_name
            assigns self.name to the resulting seq.name

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
            return self._do_seq_slice(result, apply_name)
        # all slicing now requires start < end
        result = self.parent[fmap.start : fmap.end]
        return self._do_seq_slice(result, apply_name)

    def _do_seq_slice(self, result: SeqORAlign, apply_name: bool) -> SeqORAlign:
        if self.reversed:
            result = result.rc()
        if self.map.num_spans > 1:
            # db querying will be incorrect so make sure it can't be done
            result.annotation_db = None
        if apply_name:
            result.name = self.name
        return result

    def without_lost_spans(self) -> typing_extensions.Self:
        """Keeps only the parts which are actually present in the underlying sequence"""
        if self.map.complete:
            return self
        keep = self.map.nongap()
        kwargs = {**self._serialisable, "map": self.map[keep]}
        return self.__class__(**kwargs)

    def as_one_span(self, name: str | None = None) -> typing_extensions.Self:
        """returns a feature that preserves any gaps in the underlying sequence

        Parameters
        ----------
        name
            The name of the one-span feature, by default 'one-span <self name>'
        """
        kwargs = {
            **self._serialisable,
            "map": self.map.get_covering_span(),
            "name": name or f"one-span {self.name}",
        }
        return self.__class__(**kwargs)

    def shadow(self, name: str | None = None) -> typing_extensions.Self:
        """returns new instance corresponding to disjoint of self coordinates

        Parameters
        ----------
        name
            The name of the shadow feature, by default 'not <self name>'
        """
        kwargs = {
            **self._serialisable,
            "map": self.map.shadow(),
            "biotype": f"not {self.biotype}",
            "name": name or f"not {self.name}",
        }
        return self.__class__(**kwargs)

    def __len__(self) -> int:
        return len(self.map)

    def __repr__(self) -> str:
        name = self.__class__.__name__
        txt = ", ".join(
            f"{attr}={getattr(self, attr)!r}"
            for attr in ("seqid", "biotype", "name", "map", "parent")
        )
        return f"{name}({txt})"

    def remapped_to(self, grandparent, gmap) -> typing_extensions.Self:
        # grandparent can be either a Sequence or an Alignment
        if not isinstance(gmap, FeatureMap):
            # due to separation of IndelMap and Map, change class
            gmap = gmap.to_feature_map()

        seqid = getattr(grandparent, "name", None) or f"from {self.seqid!r}"
        kwargs = {
            **self._serialisable,
            "map": gmap[self.map],
            "parent": grandparent,
            "seqid": seqid,
        }
        return self.__class__(**kwargs)

    def get_coordinates(self) -> SeqCoordTypes:
        """returns sequence coordinates of this Feature as
        [(start1, end1), ...]"""
        return self.map.get_coordinates()

    def get_children(
        self,
        biotype: str | None = None,
        **kwargs,
    ) -> typing.Iterable[typing_extensions.Self]:
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

    def get_parent(self, **kwargs) -> typing.Iterable[typing_extensions.Self]:
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

    def union(self, features: Iterable) -> typing_extensions.Self:
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
                msg = "cannot merge annotations from different objects"
                raise ValueError(msg)

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
        strand = self._strand.value if same_orientation else Strand.PLUS.value
        seqid = ", ".join(seqids) if seqids else None
        biotype = ", ".join(biotypes)
        kwargs = {
            **self._serialisable,
            "map": fmap,
            "seqid": seqid,
            "biotype": biotype,
            "name": name,
            "strand": strand,
        }

        return self.__class__(**kwargs)

    def get_drawable(self) -> "Shape":
        """returns plotly trace"""
        from cogent3.draw.drawable import make_shape

        return make_shape(type_=self, parent_length=len(self.parent))

    def to_dict(self) -> dict[str, typing.Any]:
        """returns"""
        result = {
            **self._serialisable,
            "spans": self.map.get_coordinates(),
        }
        for key in ("map", "parent", "xattr"):
            result.pop(key, None)
        return result

    @property
    def reversed(self) -> bool:
        """whether Feature is on the reverse strand relative to bound object"""
        return self._strand is Strand.MINUS

    @property
    def xattr(self) -> dict[str, typing.Any] | None:
        """extra attributes for this feature"""
        return self._xattr

    @xattr.setter
    def xattr(self, val: typing.Any) -> None:
        msg = "setting xattr is not supported"
        raise TypeError(msg)
