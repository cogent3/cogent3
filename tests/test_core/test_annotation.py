#!/usr/bin/env python

import unittest

from cogent import DNA, LoadSeqs
from cogent.core.annotation import Feature, Variable
from cogent.core.location import Map, Span

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def makeSampleSequence(mid_gaps=False):
    raw_seq = 'AACCCAAAATTTTTTGGGGGGGGGGCCCC'
    cds = (15, 25)
    utr = (12, 15)
    if mid_gaps:
        rev_seq = raw_seq[:5] + '-----' +raw_seq[10:]
        raw_seq = rev_seq
        # annotations only make sense when they're on the raw sequence
        cds = (10, 20)
        utr = (5, 8)
    seq = DNA.makeSequence(raw_seq)
    seq.addAnnotation(Feature, 'CDS', 'CDS', [cds])
    seq.addAnnotation(Feature, "5'UTR", "5' UTR", [utr])
    return seq

def makeSampleAlignment():
    seq1 = makeSampleSequence()
    seq2 = makeSampleSequence(mid_gaps=True)
    seqs = {'FAKE01': seq1, 'FAKE02': seq2}
    aln = LoadSeqs(data = seqs)
    aln.addAnnotation(Feature, 'misc_feature', 'misc', [(12,25)])
    aln.addAnnotation(Feature, 'CDS', 'blue', [(15, 25)])
    aln.addAnnotation(Feature, "5'UTR", 'red', [(2, 4)])
    aln.addAnnotation(Feature, "LTR", "fake", [(2,15)])
    return aln

class TestAnnotations(unittest.TestCase):
    def setUp(self):
        self.seq = makeSampleSequence()
        self.aln = makeSampleAlignment()
    
    def test_slice_seq_with_annotations(self):
        newseq = self.seq[:5] + self.seq[10:]
        for annot_type in ["CDS", "5'UTR"]:
            orig = str(list(self.seq.getByAnnotation(annot_type))[0])
            new = str(list(newseq.getByAnnotation(annot_type))[0])
            assert orig == new, (annot_type, orig, new)
    
    def test_aln_annotations(self):
        """test that annotations to alignment and its' sequences"""
        aln_expecteds = {"misc_feature":{'FAKE01': 'TTTGGGGGGGGGG',
                                     'FAKE02': 'TTTGGGGGGGGGG'},
                     "CDS": {'FAKE01': 'GGGGGGGGGG', 'FAKE02': 'GGGGGGGGGG'},
                     "5'UTR": {'FAKE01': 'CC', 'FAKE02': 'CC'},
                     "LTR" : {"FAKE01": "CCCAAAATTTTTT",
                              "FAKE02": "CCC-----TTTTT"}
                    }
        seq_expecteds = {"CDS": {"FAKE01": "GGGGGGGGGG",
                                "FAKE02": "GGGGGGGGGG"},
                        "5'UTR": {"FAKE01": "TTT",
                                "FAKE02": "TTT"}}
        for annot_type in ["misc_feature", "CDS", "5'UTR", "LTR"]:
            observed = list(self.aln.getByAnnotation(annot_type))[0].todict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, (annot_type, expected, observed)
            if annot_type in ["misc_feature", "LTR"]:
                continue # because seqs haven't been annotated with it
            for name in self.aln.Names:
                observed = list(self.aln.NamedSeqs[name].data.\
                                        getByAnnotation(annot_type))[0]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, (annot_type, name, expected,
                                                    observed)
    
    def test_slice_aln_with_annotations(self):
        """test that annotations of sequences and alignments survive alignment
        slicing."""
        aln_expecteds = {"misc_feature":{'FAKE01': 'TTTGGGGGGGGGG',
                                         'FAKE02': 'TTTGGGGGGGGGG'},
                     "CDS": {'FAKE01': 'GGGGGGGGGG', 'FAKE02': 'GGGGGGGGGG'},
                     "5'UTR": {'FAKE01': 'CC', 'FAKE02': 'CC'},
                     "LTR" : {"FAKE01": "CCCTTTTT",
                              "FAKE02": "CCCTTTTT"}}
        newaln = self.aln[:5]+self.aln[10:]
        feature_list = newaln.getAnnotationsMatching("LTR")
        for annot_type in ["LTR", "misc_feature", "CDS", "5'UTR"]:
            feature_list = newaln.getAnnotationsMatching(annot_type)
            new = newaln.getRegionCoveringAll(feature_list).getSlice().todict()
            expected = aln_expecteds[annot_type]
            assert expected == new, (annot_type, expected, new)
            if annot_type in ["misc_feature", "LTR"]:
                continue # because seqs haven't been annotated with it
            for name in self.aln.Names:
                orig = str(list(self.aln.getAnnotationsFromSequence(name,
                                                    annot_type))[0].getSlice())
                new = str(list(newaln.getAnnotationsFromSequence(name,
                                                    annot_type))[0].getSlice())
                assert orig == new, (name, annot_type, orig, new)
    
    def test_reversecomplement(self):
        """test correct translation of annotations on reverse complement."""
        aln_expecteds = {"misc_feature":{'FAKE01': 'TTTGGGGGGGGGG',
                                     'FAKE02': 'TTTGGGGGGGGGG'},
                     "CDS": {'FAKE01': 'GGGGGGGGGG', 'FAKE02': 'GGGGGGGGGG'},
                     "5'UTR": {'FAKE01': 'CC', 'FAKE02': 'CC'},
                     "LTR" : {"FAKE01": "CCCAAAATTTTTT",
                              "FAKE02": "CCC-----TTTTT"}
                    }
        
        seq_expecteds = {"CDS": {"FAKE01": "GGGGGGGGGG",
                                "FAKE02": "GGGGGGGGGG"},
                        "5'UTR": {"FAKE01": "TTT",
                                "FAKE02": "TTT"}}
        
        rc = self.aln.rc()
        # rc'ing an Alignment or Sequence rc's their annotations too. This means
        # slicing returns the same sequence as the non-rc'd alignment/seq
        for annot_type in ["misc_feature", "CDS", "5'UTR", "LTR"]:
            observed = list(self.aln.getByAnnotation(annot_type))[0].todict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, ("+", annot_type, expected, observed)
            observed = list(rc.getByAnnotation(annot_type))[0].todict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, ("-", annot_type, expected, observed)
            
            if annot_type in ["misc_feature", "LTR"]:
                continue # because seqs haven't been annotated with it
            for name in self.aln.Names:
                observed = list(self.aln.NamedSeqs[name].data.\
                                        getByAnnotation(annot_type))[0]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, ("+", annot_type, name, expected,
                                                    observed)
                observed = list(rc.NamedSeqs[name].data.\
                                        getByAnnotation(annot_type))[0]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, ("-", annot_type, name, expected,
                                                    observed)

class TestMapSpans(unittest.TestCase):
    """Test attributes of Map & Spans classes critical to annotation
    manipulation."""
    def test_span(self):
        length = 100
        forward = Span(20, 30)
        reverse = Span(70, 80, Reverse=True)
        assert forward.reversedRelativeTo(100) == reverse
        assert reverse.reversedRelativeTo(100) == forward
    
    def test_map(self):
        """reversing a map with multiple spans should preserve span relative
        order"""
        forward = [Span(20,30), Span(40,50)]
        fmap = Map(spans=forward, parent_length=100)
        fmap_reversed = fmap.nucleicReversed()
        reverse = [Span(70,80, Reverse=True), Span(50,60, Reverse=True)]
        rmap = Map(spans=reverse, parent_length=100)
        for i in range(2):
            self.assertEquals(fmap_reversed.spans[i], rmap.spans[i])
    
    

if __name__ == '__main__':
    unittest.main()
