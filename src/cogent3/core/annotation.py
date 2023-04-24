import copy
import json

from collections import defaultdict
from fnmatch import fnmatch
from typing import Optional

import numpy

from cogent3.util.misc import get_object_provenance

from .location import Map, as_map


__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Katherine Caley"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


class _AnnotationMixin:
    def add_feature(self, biotype, name, spans):
        feat = Feature(self, biotype, name, spans)

        if feat.parent is not self:
            raise ValueError("doesn't belong here")

        if feat.attached:
            raise ValueError("already attached")

        if self.annotations is self.__class__.annotations:
            self.annotations = []

        self.annotations.extend([feat])

        feat.attached = True

        return feat

    def get_features_matching(self, feature_type, name=None, extend_query=False):
        """

        Parameters
        ----------
        feature_type : string
            name of the feature type. Wild-cards allowed.
        name : string
            name of the instance. Wild-cards allowed.
        extend_query : boolean
            queries sub-annotations if True
        Returns
        -------
        list of AnnotatableFeatures
        """
        result = []
        if len(self.annotations) == 0:
            return result
        for annotation in self.annotations:
            if fnmatch(annotation.type, feature_type) and (
                name is None or fnmatch(annotation.name, name)
            ):
                result.append(annotation)
            if extend_query:
                result.extend(
                    annotation.get_features_matching(
                        feature_type, name, extend_query=extend_query
                    )
                )
        return result

    def get_drawables(self):
        """returns a dict of drawables, keyed by type"""
        result = defaultdict(list)
        for a in self.annotations:
            result[a.type].append(a.get_drawable())
        return result

    def get_drawable(self, width=600, vertical=False):
        """returns Drawable instance"""
        from cogent3.draw.drawable import Drawable

        drawables = self.get_drawables()
        if not drawables:
            return None
        # we order by tracks
        top = 0
        space = 0.25
        annotes = []
        for feature_type in drawables:
            new_bottom = top + space
            for i, annott in enumerate(drawables[feature_type]):
                annott.shift(y=new_bottom - annott.bottom)
                if i > 0:
                    annott._showlegend = False
                annotes.append(annott)

            top = annott.top

        top += space
        height = max((top / len(self)) * width, 300)
        xaxis = dict(range=[0, len(self)], zeroline=False, showline=True)
        yaxis = dict(range=[0, top], visible=False, zeroline=True, showline=True)

        if vertical:
            all_traces = [t.T.as_trace() for t in annotes]
            width, height = height, width
            xaxis, yaxis = yaxis, xaxis
        else:
            all_traces = [t.as_trace() for t in annotes]

        drawer = Drawable(
            title=self.name, traces=all_traces, width=width, height=height
        )
        drawer.layout.update(xaxis=xaxis, yaxis=yaxis)
        return drawer


class _AnnotationCore:
    """core methods for annotation handling, used by both new and old
    style features"""

    def _sliced_annotations(self, new, slice):
        result = []
        if self.annotations:
            slicemap = self._as_map(slice)
            # try:
            newmap = slicemap.inverse()
            # except ValueError, detail:
            #    print "Annotations dropped because %s" % detail
            #    return []
            if slicemap.useful:
                for annot in self.annotations:
                    if not annot.map.useful:
                        continue
                    if (
                        annot.map.start < slicemap.end
                        and annot.map.end > slicemap.start
                    ):
                        annot = annot.remapped_to(new, newmap)
                        if annot.map.useful:
                            result.append(annot)
        return result

    def _as_map(self, index):
        """Can take a slice, integer, map or feature, or even a list/tuple of those"""
        if type(index) in [list, tuple]:
            spans = []
            for i in index:
                spans.extend(self._as_map(i).spans)
            map = Map(spans=spans, parent_length=len(self))
        elif isinstance(index, _Feature):
            feature = index
            map = feature.map
            base = feature.parent
            containers = []
            while feature and base is not self and hasattr(base, "parent"):
                containers.append(base)
                base = base.parent
            if base is not self:
                raise ValueError(
                    f"Can't map {index} onto {repr(self)} via {containers}"
                )
            for base in containers:
                feature = feature.remapped_to(base, base.map)
        else:
            map = as_map(index, len(self))
        return map

    def __getitem__(self, index):
        map = self._as_map(index)
        new = self._mapped(map)
        sliced_annots = self._sliced_annotations(new, map)
        new.attach_annotations(sliced_annots)
        if hasattr(self, "_repr_policy"):
            new._repr_policy.update(self._repr_policy)
        return new

    def _mapped(self, map):
        raise NotImplementedError

    def get_region_covering_all(
        self, annotations, feature_class=None, extend_query=False
    ):
        from cogent3.core.sequence import Sequence

        if isinstance(self, Sequence):
            from cogent3.util.warning import deprecated

            deprecated(
                "method",
                "Sequence.get_region_covering_all",
                "_Annotatable.get_region_covering_all",
                " 2023.3",
                "method .get_region_covering_all will be discontinued for Sequence objects",
            )

        if extend_query:
            annotations = [annot._projected_to_base(self) for annot in annotations]
        spans = []
        annotation_types = []
        for annot in annotations:
            spans.extend(annot.map.spans)
            if annot.type not in annotation_types:
                annotation_types.append(annot.type)
        map = Map(spans=spans, parent_length=len(self))
        map = map.covered()  # No overlaps
        name = ",".join(annotation_types)

        if feature_class is None:
            feature_class = _Feature

        return feature_class(self, map, type="region", name=name)

    def _annotations_nucleic_reversed_on(self, new):
        """applies self.annotations to new with coordinates adjusted for
        reverse complement."""
        assert len(new) == len(self)
        annotations = []
        for annot in self.annotations:
            new_map = annot.map.nucleic_reversed()
            annotations.append(annot.__class__(new, new_map, annot))
        new.attach_annotations(annotations)


class _Annotatable(_AnnotationCore, _AnnotationMixin):
    # default
    annotations = ()

    # Subclasses should provide __init__, getOwnTracks, and a _mapped for use by
    # __getitem__

    def _shifted_annotations(self, new, shift):
        result = []
        if self.annotations:
            newmap = Map([(shift, shift + len(self))], parent_length=len(new))
            for annot in self.annotations:
                annot = annot.remapped_to(new, newmap)
                result.append(annot)
        return result

    def add_annotation(self, klass, *args, **kw):
        from cogent3.core.sequence import Sequence

        if isinstance(self, Sequence):
            from cogent3.util.warning import deprecated

            deprecated(
                "method",
                "Sequence.add_annotation",
                "_Annotatable.add_annotation",
                " 2023.3",
                "method .add_annotation will be discontinued for Sequence objects",
            )

        annot = klass(self, *args, **kw)
        self.attach_annotations([annot])
        return annot

    def clear_annotations(self):
        from cogent3.core.sequence import Sequence

        if isinstance(self, Sequence):
            from cogent3.util.warning import deprecated

            deprecated(
                "method",
                "Sequence.clear_annotations",
                "_Annotatable.clear_annotations",
                " 2023.3",
                "method .clear_annotations will be discontinued for Sequence objects",
            )

        self.annotations = []

    def attach_annotations(self, annots):
        from cogent3.core.sequence import Sequence

        if isinstance(self, Sequence):
            from cogent3.util.warning import deprecated

            deprecated(
                "method",
                "Sequence.attach_annotations",
                "_Annotatable.attach_annotations",
                " 2023.3",
                "method .attach_annotations will be discontinued for Sequence objects",
            )

        for annot in annots:
            if annot.parent is not self:
                raise ValueError("doesn't belong here")
            if annot.attached:
                raise ValueError("already attached")
        if self.annotations is self.__class__.annotations:
            self.annotations = []
        self.annotations.extend(annots)
        for annot in annots:
            annot.attached = True

    def detach_annotations(self, annots):
        from cogent3.core.sequence import Sequence

        if isinstance(self, Sequence):
            from cogent3.util.warning import deprecated

            deprecated(
                "method",
                "Sequence.detach_annotations",
                "_Annotatable.detach_annotations",
                " 2023.3",
                "method .detach_annotations will be discontinued for Sequence objects",
            )
        for annot in annots:
            if annot.parent is not self:
                raise ValueError("doesn't live here")
        for annot in annots:
            if annot.attached:
                self.annotations.remove(annot)
                annot.attached = False

    def get_annotations_matching(self, annotation_type, name=None, extend_query=False):
        """

        Parameters
        ----------
        annotation_type : string
            name of the annotation type. Wild-cards allowed.
        name : string
            name of the instance. Wild-cards allowed.
        extend_query : boolean
            queries sub-annotations if True
        Returns
        -------
        list of AnnotatableFeatures
        """
        from cogent3.util.warning import deprecated

        deprecated(
            "method",
            "get_annotations_matching",
            "get_features_matching",
            " 2023.3",
        )

        result = []
        if len(self.annotations) == 0:
            return result
        for annotation in self.annotations:
            if fnmatch(annotation.type, annotation_type) and (
                name is None or fnmatch(annotation.name, name)
            ):
                result.append(annotation)
            if extend_query:
                result.extend(
                    annotation.get_annotations_matching(
                        annotation_type, name, extend_query=extend_query
                    )
                )
        return result

    def get_by_annotation(self, annotation_type, name=None, ignore_partial=False):
        """yields the sequence segments corresponding to the specified
        annotation_type and name one at a time.

        Parameters
        ----------
        ignore_partial
            if True, annotations that extend beyond the
            current sequence are ignored.

        """
        from cogent3.core.sequence import Sequence

        if isinstance(self, Sequence):

            from cogent3.util.warning import deprecated

            deprecated(
                "method",
                "Sequence.get_by_annotation",
                "_Annotatable.get_by_annotation",
                " 2023.3",
                "method .get_by_annotation will be discontinued for Sequence objects",
            )
        for annotation in self.get_features_matching(annotation_type, name):
            try:
                seq = self[annotation.map]
            except ValueError as msg:
                if ignore_partial:
                    continue
                raise msg
            seq.info["name"] = annotation.name
            yield seq

    def _projected_to_base(self, base):
        raise NotImplementedError


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
class Annotation(_AnnotationCore, _Serialisable):
    """new style annotation, created on demand"""

    __slots__ = (
        "parent",
        "seqid",
        "map",
        "biotype",
        "name",
        "_serialisable",
        "base",  # todo gah do we need this?
        "base_map",  # todo gah do we need this?
    )

    # todo gah implement a __new__ to trap args for serialisation purposes?
    def __init__(self, *, parent, seqid: str, map: Map, biotype: str, name: str):
        d = locals()
        exclude = ("self", "__class__", "kw")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}
        self.biotype = biotype
        self.name = name
        self.parent = parent
        if isinstance(map, Map):
            assert map.parent_length == len(parent), (map, len(parent))
        else:
            map = Map(locations=map, parent_length=len(parent))

        self.map = map

    def get_drawable(self):
        """returns plotly trace"""
        from cogent3.draw.drawable import make_shape

        return make_shape(type_=self)

    def _mapped(self, slicemap):
        # the self._serialisable dict is used for serialisation, so we need to
        # use a copy to create the new instance
        kwargs = {
            **self._serialisable,
            **{"map": slicemap, "name": f"{repr(slicemap)} of {self.name}"},
        }
        return self.__class__(**kwargs)

    def get_slice(self, complete=True):
        """The corresponding sequence fragment.  If 'complete' is true
        and the full length of this feature is not present in the sequence
        then this method will fail."""
        map = self.map
        if not (complete or map.complete):
            map = map.without_gaps()
        return self.parent[map]

    def without_lost_spans(self):
        """Keeps only the parts which are actually present in the underlying sequence"""
        if self.map.complete:
            return self
        keep = self.map.nongap()
        return self.__class__(
            self.parent, self.map[keep], biotype=self.biotype, name=self.name
        )

    def as_one_span(self):
        new_map = self.map.get_covering_span()
        return self.__class__(
            self.parent, new_map, biotype=self.biotype, name=f"one-span {self.name}"
        )

    def get_shadow(self):
        return self.__class__(
            self.parent, self.map.shadow(), type="region", name=f"not {self.name}"
        )

    def __len__(self):
        return len(self.map)

    def __repr__(self):
        name = f' "{self.name}"'
        return f"{self.biotype}{name} at {self.map}"

    def _projected_to_base(self, base):
        if self.parent == base:
            return self.__class__(base, self.map)
        return self.remapped_to(base, self.parent._projected_to_base(base).map)

    def remapped_to(self, grandparent, gmap):
        map = gmap[self.map]
        return self.__class__(grandparent, map, biotype=self.biotype, name=self.name)

    def get_coordinates(self):
        """returns sequence coordinates of this Feature as
        [(start1, end1), ...]"""
        return self.map.get_coordinates()

    def get_children(self, biotype: Optional[str] = None):
        """generator returns sub-features of self optionally matching biotype"""
        make_feature = self.parent.make_feature
        db = self.parent.annotation_db
        for child in db.get_feature_children(biotype=biotype, name=self.name):
            yield make_feature(child)


class FeatureNew(_AnnotationMixin, _Serialisable):
    def __init__(self, parent, biotype, name, spans, reversed):
        self.parent = parent
        self.biotype = biotype
        self.name = name
        self.reversed = reversed
        self.map = Map(locations=spans, parent_length=len(parent))

    # remapped_to takes an alignment and the "inverted" Aligned map and returns a new instance of self
    def from_seq_to_align(self, aln, aligned):
        """transfer instance into aln through alignment coords

        Parameters
        ----------
        aln : Alignment
            self should be bound to an Aligned member of obj
        """
        assert len(self.parent) == len(aln)
        # the Map.inverse() method transforms the seq coords into alignment
        # coords, returning a new Map which can itself be sliced
        mapp = aligned.map.inverse()
        new_map = mapp[self.map]
        return self.__class__(
            aln,
            biotype=self.biotype,
            name=self.name,
            spans=new_map.get_coordinates(),
            reversed=self.reversed,
        )

    def get_slice(self):
        # feature maps and spans are always oriented with respect to the + stand sequence, this means
        # that when we really do need to be grabbing out the - strand sequence, that the features need to be
        # be adjusted so that the features spans need to be made relative w.r.t. the end.

        # todo: need to cope with cases where span is a lostSpan...

        # separation of concerns issue, this prevents slicing using an alignment
        parent_reversed = self.parent._seq.reversed
        feature_reversed = self.reversed

        if not parent_reversed and not feature_reversed:
            segments = [
                str(self.parent[span.start : span.end]) for span in self.map.spans
            ]

        elif not parent_reversed and feature_reversed:
            # because the underlying sequence is unchanged, just grab the segments out as strings
            # join them, and then do the reverse complement
            segments = [
                str(self.parent[span.start : span.end]) for span in self.map.spans
            ]
            segments = self.parent.moltype.complement("".join(segments)[::-1])

        elif parent_reversed and not feature_reversed:
            parent = self.parent[::-1]
            segments = [str(parent[span.start : span.end]) for span in self.map.spans]

        elif parent_reversed and feature_reversed:
            segments = [
                str(self.parent[span.start : span.end])
                for span in self.map.nucleic_reversed().reversed()
            ]

        return self.parent.moltype.make_seq(seq="".join(segments), name=self.name)

    def get_children(self, biotype=None):
        for child in self.parent.annotation_db.get_feature_children(
            biotype=biotype, name=self.name
        ):
            yield Annotation(self.parent, **child)

    def __len__(self):
        return len(self.map)

    def __repr__(self):
        return f"{self.name} at {self.map}"


# https://pythonspeed.com/products/filmemoryprofiler/


class _Feature(_Annotatable, _Serialisable):
    qualifier_names = ["type", "name"]
    __slots__ = ["parent", "map", "original", "_serialisable", "base", "base_map"]

    def __init__(self, parent, map, original=None, **kw):
        d = locals()
        exclude = ("self", "__class__", "kw")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}
        self._serialisable.update(kw)

        assert isinstance(parent, _Annotatable), parent
        self.parent = parent
        self.attached = False
        if isinstance(map, Map):
            assert map.parent_length == len(parent), (map, len(parent))
        else:
            map = Map(locations=map, parent_length=len(parent))

        self.map = map
        if hasattr(parent, "base"):
            self.base = parent.base
            self.base_map = parent.base_map[self.map]
        else:
            self.base = parent
            self.base_map = map

        for n in self.qualifier_names:
            if n in kw:
                setattr(self, n, kw.pop(n))
            else:
                val = getattr(original, n)
                setattr(self, n, val)
                self._serialisable[n] = val
        assert not kw, kw

    def get_drawable(self):
        """returns plotly trace"""
        from cogent3.draw.drawable import make_shape

        return make_shape(type_=self)

    def attach(self):
        self.parent.attach_annotations([self])

    def detach(self):
        self.parent.detach_annotations([self])

    def _mapped(self, slicemap):
        name = f"{repr(slicemap)} of {self.name}"
        return self.__class__(self, slicemap, type="slice", name=name)

    def get_slice(self, complete=True):
        """The corresponding sequence fragment.  If 'complete' is true
        and the full length of this feature is not present in the sequence
        then this method will fail."""
        map = self.base_map
        if not (complete or map.complete):
            map = map.without_gaps()
        return self.base[map]

    def without_lost_spans(self):
        """Keeps only the parts which are actually present in the underlying sequence"""
        if self.map.complete:
            return self
        keep = self.map.nongap()
        new = self.__class__(self.parent, self.map[keep], original=self)
        if self.annotations:
            sliced_annots = self._sliced_annotations(new, keep)
            new.attach_annotations(sliced_annots)
        return new

    def as_one_span(self):
        new_map = self.map.get_covering_span()
        return self.__class__(self.parent, new_map, type="span", name=self.name)

    def get_shadow(self):
        return self.__class__(
            self.parent, self.map.shadow(), type="region", name=f"not {self.name}"
        )

    def __len__(self):
        return len(self.map)

    def __repr__(self):
        name = getattr(self, "name", "")
        if name:
            name = f' "{name}"'
        return f"{self.type}{name} at {self.map}"

    def _projected_to_base(self, base):
        if self.parent == base:
            return self.__class__(base, self.map, original=self)
        return self.remapped_to(base, self.parent._projected_to_base(base).map)

    def remapped_to(self, grandparent, gmap):
        map = gmap[self.map]
        return self.__class__(grandparent, map, original=self)

    def get_coordinates(self):
        """returns sequence coordinates of this Feature as
        [(start1, end1), ...]"""
        return self.map.get_coordinates()

    def copy_annotations_to(self, annotatable):
        """copies annotations to another annotatable object

        Parameters
        ----------
        annotatable : _Annotatable
            another annotatable object

        Returns
        -------
        the processed annotatable object
        """
        assert len(annotatable) == len(self.parent)
        serialisable = copy.deepcopy(self._serialisable)
        serialisable["parent"] = annotatable
        new = self.__class__(**serialisable)
        annotatable.attach_annotations([new])
        return annotatable


class AnnotatableFeature(_Feature):
    """These features can themselves be annotated."""

    def _mapped(self, slicemap):
        new_map = self.map[slicemap]
        return self.__class__(self.parent, new_map, type="slice", name="")

    def remapped_to(self, grandparent, gmap):
        new = _Feature.remapped_to(self, grandparent, gmap)
        new.annotations = [annot for annot in self.annotations if annot.map.useful]
        return new


class Source(_Feature):
    # Has two maps - where it is on the sequence it annotates, and
    # where it is on the original sequence.
    type = "source"

    def __init__(self, seq, map, accession, basemap):
        d = locals()
        exclude = ("self", "__class__")
        self._serialisable = {k: v for k, v in d.items() if k not in exclude}

        self.accession = accession
        self.name = repr(basemap) + " of " + accession
        self.parent = seq
        self.attached = False
        self.map = map
        self.basemap = basemap

    def remapped_to(self, grandparent, gmap):
        new_map = gmap[self.map]
        # unlike other annotations, sources are divisible, so throw
        # away gaps.  since they don't have annotations it's simple.
        ng = new_map.nongap()
        new_map = new_map[ng]
        basemap = self.basemap[ng]
        return self.__class__(grandparent, new_map, self.accession, basemap)

    def without_lost_spans(self):
        return self


def Feature(parent, type, name, spans, value=None):
    if isinstance(spans, Map):
        map = spans
        assert map.parent_length == len(parent), (map, len(parent))
    else:
        map = Map(locations=spans, parent_length=len(parent))
    return AnnotatableFeature(parent, map, type=type, name=name)


class _Variable(_Feature):
    qualifier_names = _Feature.qualifier_names + ["xxy_list"]

    def without_lost_spans(self):
        if self.map.complete:
            return self
        raise NotImplementedError


def Variable(parent, type, name, xxy_list):
    """A variable that has 2 x-components (start, end) and a single y component.
    Currently used by Vestige - BMC Bioinformatics, 6:130, 2005."""
    start = min(min(x1, x2) for ((x1, x2), y) in xxy_list)
    end = max(max(x1, x2) for ((x1, x2), y) in xxy_list)
    if start != 0:
        xxy_list = [((x1 - start, x2 - start), y) for ((x1, x2), y) in xxy_list]
        end -= start
    # values = [location.Span(x1-start, x2-start, True, True, y) for ((x1, x2), y) in xxy]
    map = Map([(start, end)], parent_length=len(parent))
    return _Variable(parent, map, type=type, name=name, xxy_list=xxy_list)


class _SimpleVariable(_Feature):
    qualifier_names = _Feature.qualifier_names + ["data"]

    def without_lost_spans(self):
        if self.map.complete:
            return self
        keep = self.map.nongap()
        indices = numpy.concatenate([list(span) for span in keep.spans])
        data = numpy.asarray(self.data)[indices]
        return self.__class__(self.parent, self.map[keep], data=data, original=self)


def SimpleVariable(parent, type, name, data):
    """A simple variable type of annotation, such as a computed property of
    a sequence that varies spatially."""
    assert len(data) == len(parent), (len(data), len(parent))
    map = Map([(0, len(data))], parent_length=len(parent))
    return _SimpleVariable(parent, map, type=type, name=name, data=data)


def make_annotation_for_seq(seq, feature):
    ...
