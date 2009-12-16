#!/usr/bin/env python

import unittest
from cogent.core.location import Map, Span
from cogent.core.annotation import _Annotatable, Feature, _Feature
from cogent import LoadSeqs, DNA

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Rob Knight",
                    "Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def SimpleAnnotation(parent, locations, Name):
    return Feature(parent, '', Name, locations)

def annotate(parent, start, end, Name):
    annot = parent.addAnnotation(SimpleAnnotation, locations=[(start, end)], Name=Name)
    return annot

def structure(a, depth=0):
    annots = [structure(annot, depth+1) for annot in a.annotations]
    if not isinstance(a, _Feature):
        return ('seq', len(a), annots)
    elif annots:
        return (a.Name, repr(a.map), annots)
    else:
        return (a.Name, repr(a.map))

class MapTest(unittest.TestCase):
    """Testing annotation of and by maps"""
    def test_spans(self):
        # a simple two part map of length 10
        map = Map([(0,5), (5,10)], parent_length=10)
        # try different spans on the above map
        for ((start, end), expected) in [
                ((0, 4), "[0:4]"),
                ((0, 5), "[0:5]"),
                ((0, 6), "[0:5, 5:6]"),
                ((5, 10), "[5:10]"),
                ((-1, 10), "[-1-, 0:5, 5:10]"),
                ((5, 11), "[5:10, -1-]"),
                ((0, 10), "[0:5, 5:10]"),
                ((10, 0), "[10:5, 5:0]"),
                ]:
            r = repr(Span(start, end, Reverse=start>end).remapWith(map))
            #print (start, end), r,
            if r != expected:
                self.fail(repr((r, expected)))
            
    def test_maps_on_maps(self):
        seq = DNA.makeSequence('ATCGATCGAT' * 5, Name='base')
        feat1 = annotate(seq, 10, 20, 'fake')
        feat2 = annotate(feat1, 3, 5, 'fake2')
        feat3 = annotate(seq, 1, 3, 'left')
        
        seq2 = seq[5:]
        self.assertEqual(structure(seq), ('seq', 50,
            [('fake', '[10:20]/50',
                [('fake2', '[3:5]/10')]),
            ('left', '[1:3]/50')])
            )
        self.assertEqual(structure(seq2), ('seq', 45,
            [('fake', '[5:15]/45',
                [('fake2', '[3:5]/10')])])
            )
            
    
    def test_getByAnnotation(self):
        seq = DNA.makeSequence('ATCGATCGAT' * 5, Name='base')
        seq.addAnnotation(Feature, 'test_type', 'test_label', [(5,10)])
        seq.addAnnotation(Feature, 'test_type', 'test_label2', [(15,18)])
        
        answer = list(seq.getByAnnotation('test_type'))
        self.assertEqual( len(answer), 2)
        self.assertEqual( str(answer[0]), 'TCGAT')
        self.assertEqual( str(answer[1]), 'TCG')
        
        answer = list(seq.getByAnnotation('test_type', 'test_label'))
        self.assertEqual( len(answer), 1)
        self.assertEqual( str(answer[0]), 'TCGAT')
        
        # test ignoring of a partial annotation
        sliced_seq = seq[:17]
        answer = list(sliced_seq.getByAnnotation('test_type', ignore_partial=True))
        self.assertEqual(len(answer), 1)
        self.assertEqual( str(answer[0]), 'TCGAT')
    
    def test_getBySequenceAnnotation(self):
        aln = LoadSeqs(data={
                'a': 'ATCGAAATCGAT',
                'b': 'ATCGA--TCGAT'})
        b = aln.getSeq('b')
        b.addAnnotation(Feature, 'test_type', 'test_label', [(4,6)])
        
        answer = aln.getBySequenceAnnotation('b', 'test_type')[0].todict()
        self.assertEqual(answer, {'b':'A--T', 'a':'AAAT'})
    

if 0:  # old, needs fixes
    # Maps
    a = Map([(10,20)], parent_length=100)
    
    for (desc, map, expected) in [
        ('a ', a,                                "Map([10:20] on base)"),
        ('i ', a.inverse(),                     "Map([-10-, 0:10, -80-] on Map([10:20] on base))"),
        ('1 ', a[5:],                            "Map([5:10] on Map([10:20] on base))"),
        ('1r', a[5:].relative_to(b),            "Map([15:20] on base)"),
        ('2 ', a[:5],                            "Map([0:5] on Map([10:20] on base))"),
        ('2r', a[:5].relative_to(b),            "Map([10:15] on base)"),
        ('r ', a.relative_to(a[5:]),            "Map([-5-, 0:5] on Map([5:10] on Map([10:20] on base)))"),
        ('r ', a[2:4].relative_to(a[2:6]),      "Map([0:2] on Map([2:6] on Map([10:20] on base)))"),
        ('r ', a[2:4].relative_to(a[2:6][0:3]), "Map([0:2] on Map([0:3] on Map([2:6] on Map([10:20] on base))))")]:
        print desc, repr(map),
        if repr(map) == expected:
            print
        else:
            print ' <--- ', expected
            bad = True

if __name__ == '__main__':
        unittest.main()
        
