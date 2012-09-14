from location import as_map, Map
import numpy

__author__ = "Peter Maxwell and Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class _Annotatable(object):
    # default
    annotations = ()
    
    # Subclasses should provide __init__, getOwnTracks, and a _mapped for use by
    # __getitem__
    
    def _slicedAnnotations(self, new, slice):
        result = []
        if self.annotations:
            slicemap = self._as_map(slice)
            #try:
            newmap = slicemap.inverse()
            #except ValueError, detail:
            #    print "Annotations dropped because %s" % detail
            #    return []
            if slicemap.useful:
                for annot in self.annotations:
                    if not annot.map.useful:
                        continue
                    if annot.map.Start < slicemap.End and \
                            annot.map.End > slicemap.Start:
                        annot = annot.remappedTo(new, newmap)
                        if annot.map.useful:
                            result.append(annot)
        return result
    
    def _shiftedAnnotations(self, new, shift):
        result = []
        if self.annotations:
            newmap = Map([(shift, shift+len(self))], parent_length=len(new))
            for annot in self.annotations:
                annot = annot.remappedTo(new, newmap)
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
            while feature and base is not self and hasattr(base, 'parent'):
                containers.append(base)
                base = base.parent
            if base is not self:
                raise ValueError("Can't map %s onto %s via %s" % (index, repr(self), containers))
            for base in containers:
                feature = feature.remappedTo(base, base.map)
            index = map
        else:
            map = as_map(index, len(self))
        return map
    
    def __getitem__(self, index):
        map = self._as_map(index)
        new = self._mapped(map)
        sliced_annots = self._slicedAnnotations(new, map)
        new.attachAnnotations(sliced_annots)
        return new
    
    def _mapped(self, map):
        raise NotImplementedError
    
    def getAnnotationTracks(self, policy):
        result = []
        for annot in self.annotations:
            result.extend(annot.getTracks(policy))
        return result
    
    def addAnnotation(self, klass, *args, **kw):
        annot = klass(self, *args, **kw)
        self.attachAnnotations([annot])
        return annot
    
    def attachAnnotations(self, annots):
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
    
    def detachAnnotations(self, annots):
        for annot in annots:
            if annot.parent is not self:
                raise ValueError("doesn't live here")
        for annot in annots:
            if annot.attached:
                self.annotations.remove(annot)
                annot.attached = False
    
    def addFeature(self, type, Name, spans):
        return self.addAnnotation(Feature, type, Name, spans)
    
    def getAnnotationsMatching(self, annotation_type, Name=None):
        result = []
        for annotation in self.annotations:
            if annotation_type == annotation.type and (
                    Name is None or Name == annotation.Name):
                result.append(annotation)
        return result
    
    def getRegionCoveringAll(self, annotations):
        spans = []
        annotation_types = []
        for annot in annotations:
            spans.extend(annot.map.spans)
            if annot.type not in annotation_types:
                annotation_types.append(annot.type)
        map = Map(spans=spans, parent_length=len(self))
        map = map.covered() # No overlaps
        Name = ','.join(annotation_types)
        return _Feature(self, map, type='region', Name=Name)
    
    def getByAnnotation(self, annotation_type, Name=None, ignore_partial=False):
        """yields the sequence segments corresponding to the specified
        annotation_type and Name one at a time.
        
        Arguments:
            - ignore_partial: if True, annotations that extend beyond the
            current sequence are ignored."""
        for annotation in self.getAnnotationsMatching(annotation_type, Name):
            try:
                seq = self[annotation.map]
            except ValueError, msg:
                if ignore_partial:
                    continue
                raise msg
            seq.Info['Name'] = annotation.Name
            yield seq
    
    def _annotations_nucleic_reversed_on(self, new):
        """applies self.annotations to new with coordinates adjusted for
        reverse complement."""
        assert len(new) == len(self)
        annotations = []
        for annot in self.annotations:
            new_map = annot.map.nucleicReversed()
            annotations.append(annot.__class__(new, new_map, annot))
        new.attachAnnotations(annotations)


class _Feature(_Annotatable):
    qualifier_names = ['type', 'Name']
    
    def __init__(self, parent, map, original=None, **kw):
        assert isinstance(parent, _Annotatable), parent
        self.parent = parent
        self.attached = False
        self.map = map
        if hasattr(parent, 'base'):
            self.base = parent.base
            self.base_map = parent.base_map[self.map]
        else:
            self.base = parent
            self.base_map = map
        
        for n in self.qualifier_names:
            if n in kw:
                setattr(self, n, kw.pop(n))
            else:
                setattr(self, n, getattr(original, n))
        assert not kw, kw
    
    def attach(self):
        self.parent.attachAnnotations([self])
    
    def detach(self):
        self.parent.detachAnnotations([self])
    
    def _mapped(self, slicemap):
        Name = "%s of %s" % (repr(slicemap), self.Name)
        return _Feature(self, slicemap, type="slice", Name=Name)
    
    def getSlice(self, complete=True):
        """The corresponding sequence fragment.  If 'complete' is true
        and the full length of this feature is not present in the sequence
        then this method will fail."""
        map = self.base_map
        if not (complete or map.complete):
            map = map.withoutGaps()
        return self.base[map]
    
    def withoutLostSpans(self):
        """Keeps only the parts which are actually present in the underlying sequence"""
        if self.map.complete:
            return self
        keep = self.map.nongap()
        new = type(self)(self.parent, self.map[keep], original=self)
        if self.annotations:
            sliced_annots = self._slicedAnnotations(new, keep)
            new.attachAnnotations(sliced_annots)
        return new
    
    def asOneSpan(self):
        new_map = self.map.getCoveringSpan()
        return _Feature(self.parent, new_map, type="span", Name=self.Name)
    
    def getShadow(self):
        return _Feature(self.parent, self.map.shadow(), type='region',
                Name='not '+ self.Name)
    
    def __len__(self):
        return len(self.map)
    
    def __repr__(self):
        Name = getattr(self, 'Name', '')
        if Name: Name = ' "%s"' % Name
        return '%s%s at %s' % (self.type, Name, self.map)
    
    def remappedTo(self, grandparent, gmap):
        map = gmap[self.map]
        return self.__class__(grandparent, map, original=self)
    
    def getCoordinates(self):
        """returns sequence coordinates of this Feature as
        [(start1, end1), ...]"""
        coordinates = [(span.Start, span.End) for span in self.map.spans]
        return coordinates
    

class AnnotatableFeature(_Feature):
    """These features can themselves be annotated."""
    def _mapped(self, slicemap):
        new_map = self.map[slicemap]
        return _Feature(self.parent, new_map, type='slice', Name='')
    
    def remappedTo(self, grandparent, gmap):
        new = _Feature.remappedTo(self, grandparent, gmap)
        new.annotations = [annot for annot in self.annotations if annot.map.useful]
        return new
    
    def getTracks(self, policy):
        return policy.at(self.map).tracksForFeature(self)
    

class Source(_Feature):
    # Has two maps - where it is on the sequence it annotates, and
    # where it is on the original sequence.
    type = 'source'
    
    def __init__(self, seq, map, accession, basemap):
        self.accession = accession
        self.Name = repr(basemap) + ' of ' + accession
        self.parent = seq
        self.attached = False
        self.map = map
        self.basemap = basemap
    
    def remappedTo(self, grandparent, gmap):
        new_map = gmap[self.map]
        # unlike other annotations, sources are divisible, so throw
        # away gaps.  since they don't have annotations it's simple.
        ng = new_map.nongap()
        new_map = new_map[ng]
        basemap = self.basemap[ng]
        return self.__class__(grandparent, new_map, self.accession, basemap)
    
    def withoutLostSpans(self):
        return self

def Feature(parent, type, Name, spans, value=None):
    if isinstance(spans, Map):
        map = spans
        assert map.parent_length == len(parent), (map, len(parent))
    else:
        map = Map(locations=spans, parent_length=len(parent))
    return AnnotatableFeature(parent, map, type=type, Name=Name)

class _Variable(_Feature):
    qualifier_names = _Feature.qualifier_names + ['xxy_list']
    
    def getTracks(self, policy):
        return policy.tracksForVariable(self)
    
    def withoutLostSpans(self):
        if self.map.complete:
            return self
        raise NotImplementedError
        

def Variable(parent, type, Name, xxy_list):
    """A variable that has 2 x-components (start, end) and a single y component.
    Currently used by Vestige - BMC Bioinformatics, 6:130, 2005."""
    start = min([min(x1, x2) for ((x1, x2), y) in xxy_list])
    end = max([max(x1, x2) for ((x1, x2), y) in xxy_list])
    if start != 0:
        xxy_list = [((x1-start, x2-start), y) for ((x1, x2), y) in xxy_list]
        end -= start
    #values = [location.Span(x1-start, x2-start, True, True, y) for ((x1, x2), y) in xxy]
    map = Map([(start, end)], parent_length=len(parent))
    return _Variable(parent, map, type=type, Name=Name, xxy_list=xxy_list)

class _SimpleVariable(_Feature):
    qualifier_names = _Feature.qualifier_names + ['data']
    
    def getTracks(self, policy):
        return policy.tracks_for_value(self)

    def withoutLostSpans(self):
        if self.map.complete:
            return self
        keep = self.map.nongap()
        indicies = numpy.concatenate([list(span) for span in keep.Spans])
        data = numpy.asarray(data)[indicies]
        new = type(self)(self.parent, self.map[keep], data=data, original=self)
        return new
        
def SimpleVariable(parent, type, Name, data):
    """A simple variable type of annotation, such as a computed property of
    a sequence that varies spatially."""
    assert len(data) == len(parent), (len(data), len(parent))
    map = Map([(0, len(data))], parent_length=len(parent))
    return _SimpleVariable(parent, map, type=type, Name=Name, data=data)

