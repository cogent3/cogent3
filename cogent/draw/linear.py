#!/usr/bin/env python
from reportlab.graphics import shapes, renderPDF
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4

import copy
import logging
from cogent.core.location import Map, Span, _norm_slice

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.3.0.dev"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

log = logging.getLogger('cogent')

# map - a 1D space.  must support len() and in some cases [i]
# track - a panel within a sequence display holding annotations
# display - a collection of stacked tracks with a shared base sequence
# base - the map that provides the X scale
# annotation - an annotated disjoint location on a map
# span - a contiguous part of an annotation
# feature, code, variable, shade - types of anotation

# TODO: subtracks, spectra/width variables, strides, windows
#       broken ends, vscales, zero lines, sub-pixel features.
#       fuzzy ends.  more ids/qualifiers.

def nice_round_step(pixels_per_base):
    """Pick a scale step given the available resolution"""
    # xxx redo with log10?
    assert pixels_per_base > 0
    step = 1
    while step * pixels_per_base < 5:
        step *= 10
    return step

class TrackDefn(object):
    def __init__(self, tag, features, label=None):
        assert tag
        self.tag = tag
        self.label = label
        self.feature_types = features
    
    def makeTrack(self, features, **kw):
        """level=0, label=None, needs_border=False, max_y=None, height=None"""
        return Track(self.tag, features, **kw)
    
    def __iter__(self):
        return iter(self.feature_types)
    

class Track(object):
    def __init__(self, tag, features, level=0, label=None,
            needs_border=False, max_y=None, height=None):
        assert tag
        if label is None:
            label = tag
        self.tag = tag
        self.label = label
        assert isinstance(label, str), (tag, label)
        self.features = features
        self.height = max(height, max([f.height for f in self.features]))
        self.range = max_y or max(
            # xxx this works for zero-based only
            [getattr(f, 'value', None) for f in self.features]) or 0
        self.level = level
        self.needs_border = needs_border
    
    def getShapes(self, width, length, height, yrange=None, done_border=False):
        shape_list = [feature.shape(1.0*width/length, height,
                yrange or self.range) for feature in self.features]
        if self.needs_border and not done_border:
            border = shapes.Group(
                    shapes.Line(0, 0, width, 0,
                            strokeWidth=.5, strokeColor=colors.black),
                    shapes.Line(0, height, width, height,
                            strokeWidth=.5, strokeColor=colors.black)
                    )
            shape_list = [border] + shape_list
        return shape_list
    
    def __repr__(self):
        return "Track(%(tag)s,%(label)s)" % vars(self)
    

class CompositeTrack(Track):
    """Overlayed tracks"""
    def __init__(self, tag, tracks, label=None):
        if label is None:
            labels = dict([(track.label, True)
                    for track in tracks if track.label]).keys()
            if len(labels) == 1:
                label = labels[0]
            else:
                label = ''
        self.tag = tag
        self.label = label
        self.tracks = tracks
        self.height = max([track.height for track in tracks])
        self.level = max([track.level for track in tracks])
        self.range = max([track.range for track in tracks])
    
    def getShapes(self, width, length, height, yrange=None, done_border=False):
        if yrange is None:
            yrange = self.range
        shape_list = []
        for track in self.tracks:
            if track.needs_border and not done_border:
                border = shapes.Group(
                        shapes.Line(0, 0, width, 0,
                                strokeWidth=.5, strokeColor=colors.black),
                        shapes.Line(0, height, width, height,
                                strokeWidth=.5, strokeColor=colors.black)
                        )
                shape_list.append(border)
                done_border = True
            shape_list.extend(track.getShapes(width, length, height,
                    yrange=yrange, done_border=True))
        return shape_list
    

class Annotation(object):
    """A map, a style, and some values"""
    
    def __init__(self, map, *args, **kw):
        self.map = map
        self.values = self._make_values(*args, **kw)
    
    def _make_values(self, *args, **kw):
        # override for variables etc.
        return []
    
    def __repr__(self):
        return "%s at %s" % (type(self), getattr(self, 'map', '?'))
    
    # xxx styles do this.  what about maps/ others
    def shape(self, pixels_per_base, height, yrange):
        g = shapes.Group()
        posn = 0
        for span in self.map.spans:
            if not span.lost:
                g.add(self._item_shape(
                        span, self.values[posn:posn+span.length],
                        pixels_per_base, height, yrange))
            posn += span.length
        return g
    

class Code(Annotation):
    height = 20
    def _make_values(self, sequence, colours=None):
        self.characters = []
        for c in sequence:
            if c not in self.characters:
                self.characters.append(c)
        w = (0, '?')
        for c in self.characters:
            w = max(w, (shapes.String(0,0,c, fontSize=20).getEast(), c))
        self.widest_character = w[1]
        self.colours = colours or {}
        for c in self.colours.keys():
            self.colours[c.upper()] = self.colours[c.lower()] = self.colours[c]
        return sequence
    
    def _item_shape(self, span, values, pixels_per_base, height, yrange):
        w = pixels_per_base + 1
        font_size = height + (int(height)/2)
        while font_size > 6 and w > pixels_per_base:
            font_size -= 1
            w = shapes.String(0, 0, self.widest_character,
                    fontSize=font_size).getEast()
        g = shapes.Group()
        if font_size > 6:
            g.translate(0.0, 1.0+font_size/5) # desenders, rough guess
            g.scale(1.0, min(float(height-2)/font_size, 2.0))
            text_transform = shapes.scale([1.0, -1.0][span.Reverse], 1.0)
            for (i,c) in enumerate(values):
                offset = [span.Start+0.5, -1*span.End-0.5][span.Reverse]
                s = shapes.String(
                    (i+offset)*pixels_per_base, 0, c,
                    fontSize=font_size, textAnchor='middle',
                    transform = text_transform,
                    fillColor = self.colours.get(c, colors.black))
                g.add(s)
        else:
            g.add(shapes.Line(span.Start*pixels_per_base, height/2,
                    span.End*pixels_per_base, height/2))
        return g
    

class CodeLine(Annotation):
    height = 10
    def _item_shape(self, span, values, pixels_per_base, height, yrange):
        g = shapes.Group()
        g.add(shapes.Line(span.Start*pixels_per_base, height/2,
                span.End*pixels_per_base, height/2))
        return g
    

class DensityGraph(Annotation):
    height = 20
    def _make_values(self, data):
        return data
    
    def _item_shape(self, span, values, pixels_per_base, height, yrange):
        #max_density = int(1.0/pixels_per_base) + 1
        #colours = [colors.linearlyInterpolatedColor(
        #        colors.white, colors.red, 0, height, i)
        #        for i in range(max_density)]
        mid = colors.linearlyInterpolatedColor(colors.pink, colors.red,
            0, 10, 5)
        
        g = shapes.Group()
        yscale = height / yrange
        last_x = -1
        ys = None
        for (p, q) in enumerate(values):
            if q is not None:
                g.add(shapes.Circle(
                    p*pixels_per_base, q*yscale, pixels_per_base/2,
                    strokeWidth=0.0, strokeColor=colors.red,
                    fillColor=colors.red))
            
            #x = int(p * pixels_per_base)
            #y = int(q * yscale)
            #if x > last_x:
                #if ys:
                #    ys.sort()
                #    c = len(ys)
                    #g.add(shapes.Line(x, ys[0], x, ys[-1],
                    #    strokeColor=colors.pink))
                    #g.add(shapes.Line(x, ys[c/3], x, ys[2*c/3],
                    #    strokeColor=mid))
                    #g.add(shapes.Line(x, ys[c/2]-.5, x, ys[c/2]+.5,
                    #    strokeColor=colors.red))
            #    ys = []
            #    last_x = x
            #ys.append(y)
        return g
    

class Scale(Annotation):
    height = 15
    
    def _make_values(self, length, orientation=1):
        assert orientation in [-1, 1]
        self.orientation = orientation
        # Yes, a span being used as a value.  It acts like a list but
        # is much more efficiently stored.
        return Span(0, length)
    
    def _item_shape(self, span, values, pixels_per_base, height, yrange):
        step = nice_round_step(pixels_per_base)
        g = shapes.Group()
        skip_zero = 1 #pixels_per_base > 1
        if self.orientation == -1:
            (textbase, tick_lo, tick_hi) = (0, height-5, height-1)
        else:
            (textbase, tick_lo, tick_hi) = (6, 5, 1)
        
        for posn in range(int(values.Start/step+skip_zero)*step,
                values.End+1, step):
            x = posn - values.Start
            x = [span.Start+x, span.End-x][span.Reverse]
            x = (x-0.5)*pixels_per_base
            label = posn % (step*10) == 0
            g.add(shapes.Line(x, tick_lo, x, tick_hi))
            if label:
                g.add(shapes.String(x, textbase, str(posn),
                    textAnchor="middle", fontSize=height-5))
        return g
    

class Feature(Annotation):
    """An Annotation with a style and location rather than values"""
    def __init__(self, map, style, label=None, value=None):
        self.map = map
        self.style = style
        self.label = label
        self.value = value
        self.height = style.height
        #self.values = self._make_values(*args, **kw)
    
    def shape(self, pixels_per_base, height, yrange):
        return self.style(pixels_per_base, height, self.label, self.map,
                self.value, yrange)
    

class _FeatureStyle(object):
    range_required = False
    def __init__(self, fill=True, color=colors.black, min_width=0.5,
            showLabel=False, height=1, thickness=0.6, closed=True,
            one_span=False, **kw):
        opts = {}
        if fill:
            opts['fillColor'] = color
            opts['strokeColor'] = None
            opts['strokeWidth'] = 0
        else:
            opts['fillColor'] = None
            opts['strokeColor'] = color
            opts['strokeWidth'] = 1
        self.filled = fill
        self.closed = closed
        opts.update(kw)
        self.opts = opts
        self.min_width = min_width
        self.showLabel = showLabel
        self.height = height
        self.proportion_of_track = thickness
        self.one_span = one_span
    
    def __call__(self, pixels_per_base, height, label, map, value, yrange):
        #return self.FeatureClass(label, map)
        g = shapes.Group()
        last = first = None
        if self.range_required and not yrange:
            log.warning("'%s' graph values are all zero" % label)
            yrange = 1.0
        if map.useful and self.one_span:
            map = map.getCoveringSpan()
        for (i, span) in enumerate(map.spans):
            #if last is not None:
            #    g.add(shapes.Line(last*pixels_per_base, height,
            #        part.Start*pixels_per_base, height))
            if span.lost or (value is None and self.range_required):
                continue
            if span.Reverse:
                (start, end) = (span.End, span.Start)
                (tidy_start, tidy_end) = (span.tidy_end, span.tidy_start)
            else:
                (start, end) = (span.Start, span.End)
                (tidy_start, tidy_end) = (span.tidy_start, span.tidy_end)
            shape = self._item_shape(
                start*pixels_per_base, end*pixels_per_base,
                tidy_start, tidy_end, height, value, yrange,
                last=i==len(map.spans)-1)
            g.add(shape)
            last = end
            if first is None:
                first = start
        if self.showLabel and label and last is not None and height > 7:
            label_shape = shapes.String(
                    (first+last)/2*pixels_per_base, height/4+2, label,
                    textAnchor="middle", fontSize=height/2-2)
            # bug in getEast() - it thinks textAnchor=="left"
            if ((label_shape.getEast()-label_shape.x) <
                    abs(first-last)*pixels_per_base):
                g.add(label_shape)
            else:
                log.info("couldn't fit feature label '%s'" % label)
        return g
    

class _VariableThicknessFeatureStyle(_FeatureStyle):
    def _item_shape(self, start, end, tidy_start, tidy_end, height, value,
            yrange, last=False):
        if yrange:
            thickness = 1.0*value/yrange*height
        else:
            thickness = height*self.proportion_of_track
        return self._item_shape_scaled(start, end, tidy_start, tidy_end,
                height/2, max(2, thickness), last)
    

class _End(object):
    def __init__(self, x_near, x_far, y_first, y_second, **kw):
        self.x_near = x_near
        self.x_far = x_far
        self.y_first = y_first
        self.y_second = y_second
        for (n, v) in kw.items():
            setattr(self, n, v)
    
    def moveToStart(self, path):
        path.moveTo(*self.startPoint())
    
    def drawToStart(self, path):
        path.lineTo(*self.startPoint())
    
    def finish(self, path):
        path.closePath()
    
    def startPoint(self):
        return (self.x_near, self.y_first)
    

class Open(_End):
    def finish(self, path):
        self.drawToStart(path)
    
    def drawEnd(self, path):
        path.moveTo(self.x_near, self.y_second)
    

class Square(_End):
    def drawEnd(self, path):
        path.lineTo(self.x_near, self.y_second)
    

class Rounded(_End):
    def startPoint(self):
        return (self.x_near + self.dx, self.y_first)
    
    def drawEnd(self, path):
        path.curveTo(self.x_near, self.y_first, self.x_near, self.y_first,
                self.x_near, self.y_first + self.dy)
        path.lineTo(self.x_near, self.y_second - self.dy)
        path.curveTo(self.x_near, self.y_second, self.x_near, self.y_second,
                self.x_near + self.dx, self.y_second)
    

class Pointy(_End):
    def _effective_dx(self):
        return max(abs(self.dx), abs(self.dy))*self.dx/abs(self.dx)
    
    def startPoint(self):
        return (self.x_near + self._effective_dx(), self.y_first)
    
    def drawEnd(self, path):
        head_start = self.x_near + self._effective_dx()
        middle = (self.y_first + self.y_second) / 2
        if self.blunt:
            for (x, y) in [
                    (head_start, self.y_first + self.dy),
                    (self.x_near, self.y_first),
                    (self.x_near, self.y_second),
                    (head_start, self.y_second - self.dy),
                    (head_start, self.y_second)]:
                path.lineTo(x, y)
        else:
            for (x, y) in [
                    (head_start, self.y_first + self.dy),
                    (self.x_near, middle),
                    (head_start, self.y_second - self.dy),
                    (head_start, self.y_second)]:
                path.lineTo(x, y)
    

def sign(x):
    return x and x/abs(x)

class Box(_VariableThicknessFeatureStyle):
    arrow = False
    
    def _item_shape_scaled(self, start, end, tidy_start, tidy_end, middle,
            thickness, last):
        p = shapes.Path(**self.opts)
        ends = []
        (top, bottom) = (middle+thickness/2, middle-thickness/2)
        for (x1, x2, y1, y2, tidy, pointy) in [
                (start, end, bottom, top, tidy_start, False),
                (end, start, top, bottom, tidy_end, last and self.arrow)]:
            inwards = sign(x2 - x1)
            span = max(abs(x2 - x1), self.min_width)
            if pointy:
                head_size = min(thickness, span)
                head_size = max(head_size, self.min_width/2)
                height = thickness / self.proportion_of_track
                spare = (height - thickness) / 2
                end = Pointy(x1, x2, y1, y2,
                        dx=inwards*head_size/2, dy=sign(y1-y2)*spare,
                        blunt=self.blunt)
            elif not self.closed:
                end = Open(x1, x2, y1, y2)
            elif tidy:
                ry = thickness/4
                rx = min(span, ry)
                end = Rounded(x1, x2, y1, y2, dx=rx*inwards, dy=ry*sign(y2-y1))
            else:
                end = Square(x1, x2, y1, y2)
            ends.append(end)
        ends[0].moveToStart(p)
        ends[0].drawEnd(p)
        ends[1].drawToStart(p)
        ends[1].drawEnd(p)
        ends[0].finish(p)
        return p
    

class Arrow(Box):
    arrow = True
    blunt = False

class BluntArrow(Box):
    arrow = True
    blunt = True

class Diamond(_VariableThicknessFeatureStyle):
    """diamond"""
    def _item_shape_scaled(self, start, end, tidy_start, tidy_end, middle,
            thickness, last):
        x = (start+end)/2
        spread = max(abs(start-end), self.min_width) / 2
        return shapes.Polygon(
            [x-spread, middle, x, middle+thickness/2, x+spread, middle, x,
                    middle-thickness/2], **self.opts)
    

class Line(_FeatureStyle):
    """For a line segment graph"""
    range_required = True
    def _item_shape(self, start, end, tidy_start, tidy_end, height, value,
            yrange, last=False):
        altitude = value * (height-1) / yrange
        #if self.orientation < 0:
        #    altitude = height - altitude
        return shapes.Line(start, altitude, end, altitude, **self.opts)
    

class Area(_FeatureStyle):
    """For a line segment graph"""
    range_required = True
    def _item_shape(self, start, end, tidy_start, tidy_end, height, value,
            yrange, last=False):
        altitude = value * (height-1) / yrange
        #if self.orientation < 0:
        #    altitude = height - altitude
        if end < start:
            start, end = end, start
            tidy_start, tidy_end = tidy_end, tidy_start
        return shapes.Rect(start, 0, end-start, altitude, **self.opts)
    

class DisplayPolicy(object):
    def _makeFeatureStyles(self):
        return {
            #gene structure
            'misc_RNA': Box(True, colors.paleturquoise),
            'precursor_RNA': Box(True, colors.lightcyan),
            'prim_transcript': Box(True, colors.lightcyan),
            "3'clip": Box(True, colors.lightcyan),
            "5'clip": Box(True, colors.lightcyan),
            'mRNA': Box(True, colors.cyan),
            'exon': Box(True, colors.cyan),
            'intron': Box(False, colors.cyan, closed = False),
            "3'UTR": Box(True, colors.cyan),
            "5'UTR": Box(True, colors.cyan),
            'CDS': Box(True, colors.blue),
            'mat_peptide': Box(True, colors.blue),
            'sig_peptide': Box(True, colors.navy),
            'transit_peptide': Box(True, colors.navy),
            'polyA_signal': Box(True, colors.lightgreen),
            'polyA_site': Diamond(True, colors.lightgreen),
            'gene': BluntArrow(False, colors.blue,
                    showLabel=True, closed = False),
            'operon': BluntArrow(False, colors.royalblue,
                    showLabel=True, closed = False),
            #regulation
            'attenuator': Box(False, colors.red),
            'enhancer': Box(True, colors.green),
            'CAAT_signal': Diamond(True, colors.blue),
            'TATA_signal': Diamond(True, colors.teal),
            'promoter': Box(False, colors.seagreen),
            'GC_signal': Box(True, colors.purple),
            'protein_bind': Box(True, colors.orange),
            'misc_binding': Box(False, colors.black),
            '-10_signal': Diamond(True, colors.blue),
            '-35_signal': Diamond(True, colors.teal),
            'terminator': Diamond(True, colors.red),
            'misc_signal': Box(False, colors.maroon),
            'rep_origin': Box(True, colors.linen),
            'RBS': Diamond(True, colors.navy),
            #repeats
            'repeat_region': Box(True, colors.brown),
            'repeat_unit': Arrow(True, colors.brown),
            'LTR': Box(False, colors.black),
            'satellite': Box(False, colors.brown),
            'stem_loop': Box(False, colors.dimgray),
            'misc_structure': Box(False, colors.darkslategray),
            #rna genes
            'rRNA': Arrow(False, colors.darkorchid, showLabel=True),
            'scRNA': Arrow(False, colors.darkslateblue, showLabel=True),
            'snRNA': Arrow(False, colors.darkviolet, showLabel=True),
            'snoRNA': Arrow(False, colors.darkviolet, showLabel=True),
            'tRNA': Arrow(False, colors.darkturquoise, showLabel=True),
            #sequence
            'source': Box(False, colors.black, showLabel=True),
            'misc_recomb': Box(False, colors.black, showLabel=True),
            'variation': Diamond(True, colors.violet, showLabel=True),
            'domain': Box(False, colors.darkorange, showLabel=True),
            'bluediamond': Diamond(True, colors.blue),
            'reddiamond': Diamond(True, colors.red),
            'misc_feature': Box(True, colors.darkorange, showLabel=True),
            'old_sequence': Box(False, colors.darkslategray),
            'unsure': Diamond(False, colors.crimson, min_width=2,),
            'misc_difference': Diamond(False, colors.darkorange),
            'conflict': Box(False, colors.darkorange),
            'modified_base': Diamond(True, colors.black),
            'primer_bind': Arrow(False, colors.green, showLabel=True),
            'STS': Box(False, colors.black),
            'gap': Box(True, colors.gray),
            #graphs
            'blueline': Line(False, colors.blue),
            'redline': Line(False, colors.red),
            #other
            ##immune system specific
            #'C_region': Diamond(True, colors.mediumblue),
            #'N_region': Box(False, colors.linen),
            #'S_region': Box(False, colors.linen),
            #'V_region': Box(False, colors.linen),
            #'D_segment': Diamond(True, colors.mediumpurple),
            #'J_segment': Box(False, colors.linen),
            #'V_segment': Box(False, colors.linen),
            #'iDNA': Box(False, colors.grey),
            ##Mitocondria specific
            #'D-loop': Diamond(True, colors.linen),
            ##Bacterial element specific
            #'oriT': Box(False, colors.linen),
        }
    
    def _makeTrackDefns(self):
        return [TrackDefn(*args) for args in [
            ('Gene Structure',[
                'misc_RNA',
                'precursor_RNA',
                'prim_transcript',
                "3'clip",
                "5'clip",
                'mRNA',
                'exon',
                'intron',
                "3'UTR",
                "5'UTR",
                'CDS',
                'mat_peptide',
                'sig_peptide',
                'transit_peptide',
                'polyA_signal',
                'polyA_site',
                'gene',
                'operon',
                ]),
            ('Regulation',[
                'attenuator',
                'enhancer',
                'CAAT_signal',
                'TATA_signal',
                'promoter',
                'GC_signal',
                'protein_bind',
                'misc_binding',
                '-10_signal',
                '-35_signal',
                'terminator',
                'misc_signal',
                'rep_origin',
                'RBS',
                ]),
            ('Repeats',[
                'repeat_region',
                'repeat_unit',
                'LTR',
                'satellite',
                'stem_loop',
                'misc_structure',
                ]),
            ('Rna Genes',[
                'rRNA',
                'scRNA',
                'snRNA',
                'snoRNA',
                'tRNA',
                ]),
            ('Sequence',[
                'source',
                'misc_recomb',
                'domain',
                'variation',
                'bluediamond',
                'reddiamond',
                'misc_feature',
                'old_sequence',
                'unsure',
                'misc_difference',
                'conflict',
                'modified_base',
                'primer_bind',
                'STS',
                'gap',
                ]),
            ('Graphs',[
                'blueline',
                'redline',
                ]),
        ]]
    
    _default_ignored_features =    ['C_region','N_region','S_region','V_region',
        'D_segment','J_segment','V_segment','iDNA','D-loop','oriT',]
    _default_keep_unexpected_tracks = True
    
    dont_merge = []
    show_text = None # auto
    colour_sequences = False
    
    def __init__(self,
            min_feature_height = 20,
            min_graph_height = None,
            ignored_features=None,
            keep_unexpected_tracks=None,
            **kw):
        
        if min_graph_height is None:
            min_graph_height = min_feature_height * 2
        
        feature_styles = self._makeFeatureStyles()
        # yuk
        for style in feature_styles.values():
            if style.range_required:
                style.height = max(style.height, min_graph_height)
            else:
                style.height = max(style.height, min_feature_height)
        
        self._track_defns = self._makeTrackDefns()
        
        if ignored_features is None:
            ignored_features = self._default_ignored_features
        self._ignored_features = ignored_features
        if keep_unexpected_tracks is None:
            keep_unexpected_tracks = self._default_keep_unexpected_tracks
        self.keep_unexpected_tracks = keep_unexpected_tracks
        
        if not hasattr(self, '_track_map'):
            self._track_map = {}
            for track_defn in self._track_defns:
                for (level, feature_tag) in enumerate(track_defn):
                    feature_style = feature_styles[feature_tag]
                    self._track_map[feature_tag] = (
                            track_defn, level, feature_style)
            for ft in self._ignored_features:
                self._track_map[ft] = (None, 0, None)
            for ft in feature_styles:
                if ft not in self._track_map:
                    self._track_map[ft] = (None, level, feature_style)
        self.map = None
        self.depth = 0
        self.orientation = -1
        self.show_code = True
        self.show_alignment_scale = True
        self.show_scale = True
        self._logged_drops = []
        
        if kw:
            self._setattrs(kw)
    
    def _setattrs(self, kw):
        for (n,v) in kw.items():
            setattr(self, n, v)
    
    def included(self, **kw):
        new = copy.copy(self)
        new._setattrs(kw)
        return new
    
    def at(self, map):
        if map is None:
            return self
        else:
            return self.included(map=self.map[map], depth=self.depth+1)
    
    def mergeTracks(self, orig_tracks, keep_unexpected=None):
        # merge tracks with same names
        # order features within a track by level  # xxx remerge
        tracks = {}
        orig_track_tags = []
        for track in orig_tracks:
            if not track.tag in tracks:
                tracks[track.tag] = {}
                orig_track_tags.append(track.tag) # ordered list
            if not track.level in tracks[track.tag]:
                tracks[track.tag][track.level] = []
            tracks[track.tag][track.level].append(track)
        
        track_order = [track.tag for track in self._track_defns
                if track.tag in tracks]
        unexpected = [tag for tag in orig_track_tags if tag not in track_order]
        if keep_unexpected is None:
            keep_unexpected = self.keep_unexpected_tracks
        if keep_unexpected:
            track_order += unexpected
        elif unexpected:
            log.warning('dropped tracks ' + ','.join(unexpected), stacklevel=2)
        
        sorted_tracks = []
        for track_tag in track_order:
            annots = []
            levels = tracks[track_tag].keys()
            levels.sort()
            for level in levels:
                annots.extend(tracks[track_tag][level])
            if len(annots)> 1 and track_tag not in self.dont_merge:
                sorted_tracks.append(CompositeTrack(track_tag, annots))
            else:
                sorted_tracks.extend(annots)
        return sorted_tracks
    
    def tracksForAlignment(self, alignment):
        annot_tracks = alignment.getAnnotationTracks(self)
        seq_tracks = alignment.getChildTracks(
                self.included(show_scale = not self.show_alignment_scale))
        if self.show_alignment_scale:
            own_tracks = [Track('scale',
                    [Scale(self.map, len(alignment), self.orientation)],
                    level=2, label='')]
        else:
            own_tracks = []
        annot_tracks = self.mergeTracks(annot_tracks)
        return seq_tracks + annot_tracks + own_tracks
    
    def tracksForSequence(self, sequence=None):
        result = []
        length = None
        if length is None and sequence is not None:
            length = len(sequence)
        
        label = getattr(self, 'seqname', '')
        
        if self.show_code and sequence is not None:
            # this should be based on resolution, not rowlen, but that's all
            # we have at this point
            # also these 2 tracks should be able to share a label.
            show_text = self.show_text
            if show_text is None:
                show_text = self.rowlen <= 100
            if show_text and self.rowlen <= 200:
                if self.colour_sequences:
                    colours = sequence.getColourScheme(colors)
                else:
                    colours = {}
                feature = Code(self.map, str(sequence), colours=colours)
                label1 = label
                label2 = ''
            else:
                feature = CodeLine(self.map, str(sequence))
                label1 = label
                label2 = ''
            result.append(Track('seq', [feature], level=2, label=label1))
        else:
            label2 = label
        if self.show_scale and length is not None:
            result.append(Track('scale',
                    [Scale(self.map, length, self.orientation)],
                    level=2, label=label2))
        
        annot_tracks = sequence.getAnnotationTracks(self)
        annot_tracks = self.mergeTracks(annot_tracks)
        
        return annot_tracks + result
    
    def getStyleDefnForFeature(self, feature):
        if feature.type in self._track_map:
            (track_defn, level, style) = self._track_map[feature.type]
        elif self.keep_unexpected_tracks:
            (track_defn, level, style) = self._track_map['misc_feature']
        else:
            if feature.type not in self._logged_drops:
                log.warning('dropped feature ' + repr(feature.type))
                self._logged_drops.append(feature.type)
            return (None, None, None)
        
        if track_defn is None:
            log.info('dropped feature ' + repr(feature.type))
            return (None, None, None)
        else:
            track_tag = track_defn.tag or feature.type
        return (track_tag, style, level)
    
    def tracksForFeature(self, feature):
        (track_tag, style, level) = self.getStyleDefnForFeature(feature)
        if style is None:
            return []
        annot_tracks = feature.getAnnotationTracks(self)
        return annot_tracks + [Track(track_tag,
                [Feature(self.map, style, feature.Name)], level=level)]
    
    def tracksForVariable(self, variable):
        (track_tag, style, level) = self.getStyleDefnForFeature(variable)
        if style is None:
            return []
        segments = []
        max_y = 0.0
        for ((x1, x2), y) in variable.xxy_list:
            map = self.map[x1:x2]
            segments.append(Feature(map, style, variable.Name, value=y))
            if type(y) is tuple: y = max(y)
            if y > max_y: max_y = y
        return [Track(track_tag, segments, max_y=max_y, needs_border=True,
                label=variable.Name, level=level)]
    

class Display(object):
    """Holds a list of tracks and displays them all aligned"""
    
    def __init__(self, base, policy=DisplayPolicy, _policy=None, pad=1,
            yrange=None, label_width=None, **kw):
        self.pad = pad
        self.base = base
        self.yrange = yrange
        self.label_width = label_width
        assert len(base) > 0, len(base)
        
        if _policy is None:
            policy = policy(**kw).included(
                map=Map([(0, len(base))], parent_length=len(base)),
                        depth=0, rowlen=len(base))
        else:
            policy = _policy
        self.policy = policy
        self._calc_tracks()
    
    def __len__(self):
        return len(self.policy.map.inverse())
    
    def _calc_tracks(self):
        y = 0
        self._tracks = []
        for p in self.base.getTracks(self.policy)[::-1]:
            if not isinstance(p, Track):
                if not isinstance(p, list):
                    p = [p]
                p = Track('', p)
            y2 = y + p.height + self.pad
            self._tracks.append((y+self.pad/2, (y+y2)/2, p))
            y = y2
        self.height = y
        
        if self.yrange is None:
            self.yrange = {}
            for (y, ym, p) in self._tracks:
                self.yrange[p.tag] = max(self.yrange.get(p.tag, 0), p.range)
        
        if self.label_width is None:
            self.label_width = 0
            for (y, ym, p) in self._tracks:
                self.label_width = max(self.label_width,
                        shapes.String(0, 0, p.label).getEast()+2)
    
    def included(self, **kw):
        import copy
        new = copy.copy(self)
        new.policy = self.policy.included(**kw)
        new._calc_tracks()
        return new
    
    def __getitem__(self, slice):
        return self.included(map=self.policy.map.inverse()[slice].inverse())
    
    def shape(self, label_width, data_width, border=True):
        offset = 0
        g = shapes.Group()
        
        if border:
            g.add(shapes.Rect(label_width, offset, data_width, self.height,
                fillColor=None, strokeWidth=.5, strokeColor=colors.black))
        
        for (y, ym, p) in self._tracks:
            if p.height>8 and label_width:
                g.add(shapes.String(0, offset+ym-5, p.label))
            for s in p.getShapes(data_width, self.policy.rowlen,
                    float(p.height), yrange=self.yrange[p.tag]):
                s.translate(label_width, float(y+offset))
                g.add(s)
        return g
    
    def asShapes(self, total_width, rowlen=None, wraps=None,
                border=True, withTrackLabelColumn=True):
        # Yields ReportLab group objects, one per wrap of the alignment.
        # Wrapping is caused by setting 'wraps' or 'rowlen'.
        if withTrackLabelColumn:
            label_width = self.label_width
        else:
            label_width = 0
        
        data_width = float(total_width) - label_width
        if data_width <= 0:
            raise ValueError('margins plus labels too wide')
        seq_length = len(self)
        assert seq_length
        if wraps:
            assert not rowlen
            rowlen = seq_length / wraps
            pixels_per_base = 1.0 * data_width / rowlen
            step = nice_round_step(pixels_per_base/wraps)
            rowlen = (rowlen-1 - (rowlen-1) % step) + step
        if not rowlen:
            rowlen = seq_length
        wraps = seq_length / rowlen
        if seq_length % rowlen:
            wraps += 1
        
        self = self.included(rowlen=rowlen)
        
        if wraps > 1:
            for i in range(wraps):
                s = self[i*rowlen:(i+1)*rowlen]
                g2 = s.shape(label_width, data_width, border=border)
                yield g2
        else:
            yield self.shape(label_width, data_width, border=border)
    
    def asShape(self, total_width, border=True, withTrackLabelColumn=True):
        result = list(self.asShapes(total_width,
                border=border, withTrackLabelColumn=withTrackLabelColumn))
        assert len(result) == 1
        return result[0]
    
    def asDrawing(self, total_width, rowlen=None, wraps=None, margin=10,
                border=True, withTrackLabelColumn=True):
        panels = self.asShapes(total_width-2*margin, rowlen=rowlen, wraps=wraps,
                border=border, withTrackLabelColumn=withTrackLabelColumn)
        panels = list(panels)
        wraps = len(panels)
        if wraps == 1:
            g = panels[0]
        else:
            g = shapes.Group()
            for (i,d) in enumerate(panels):
                offset = (wraps-i-1) * (self.height+margin)
                d.translate(0.0, offset)
                g.add(d)
        g.translate(margin, margin)
        D = shapes.Drawing(total_width, self.height*wraps+margin*(wraps+2))
        D.add(g)
        return D
    
    def asDrawings(self, total_width, rowlen=None, wraps=None, margin=0,
                border=True, withTrackLabelColumn=True):
        panels = self.asShapes(total_width-2*margin, rowlen=rowlen, wraps=wraps,
                border=border, withTrackLabelColumn=withTrackLabelColumn)
        for g in panels:
            g.translate(margin, margin)
            D = shapes.Drawing(total_width, self.height+2*margin)
            D.add(g)
            yield D
    
    def drawToPDF(self, filename, total_width, *args, **kw):
        """rowlen=None, wraps=None, margin=10, border=True,
        withTrackLabelColumn=True"""
        D = self.asDrawing(total_width, *args, **kw)
        renderPDF.drawToFile(D, filename)
    
