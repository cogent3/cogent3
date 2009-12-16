#!/usr/bin/env python
import unittest, sys, os
from cogent import DNA, LoadSeqs
from cogent.parse.cigar import map_to_cigar, cigar_to_map, aligned_from_cigar, \
                                slice_cigar, CigarParser

__author__ = "Hua Ying"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Hua Ying", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Hua Ying"
__email__ = "hua.ying@anu.edu.au"
__status__ = "Production"

class TestCigar(unittest.TestCase):
    def setUp(self):
        self.cigar_text = '3D2M3D6MDM2D3MD'
        self.aln_seq = DNA.makeSequence('---AA---GCTTAG-A--CCT-')
        self.aln_seq1 = DNA.makeSequence('CCAAAAAA---TAGT-GGC--G')
        self.map, self.seq = self.aln_seq.parseOutGaps()
        self.map1, self.seq1 = self.aln_seq1.parseOutGaps()
        self.slices = [(1, 4), (0, 8), (7, 12), (0, 1), (3, 5)]
        self.aln = LoadSeqs(data = {"FAKE01": self.aln_seq, "FAKE02": self.aln_seq1})
        self.cigars = {"FAKE01": self.cigar_text, "FAKE02": map_to_cigar(self.map1)}
        self.seqs = {"FAKE01": str(self.seq), "FAKE02": str(self.seq1)}
    
    def test_map_to_cigar(self):
        """convert a Map to cigar string"""
        assert map_to_cigar(self.map) == self.cigar_text
    
    def test_cigar_to_map(self):
        """test generating a Map from cigar"""
        map = cigar_to_map(self.cigar_text)
        assert str(map) == str(self.map)
    
    def test_aligned_from_cigar(self):
        """test generating aligned seq from cigar"""
        aligned_seq = aligned_from_cigar(self.cigar_text, self.seq)
        assert aligned_seq == self.aln_seq
    
    def test_slice_cigar(self):
        """test slicing cigars"""
        for start, end in self.slices:
            # test by_align = True
            map1, loc1 = slice_cigar(self.cigar_text, start, end)
            ori1 = self.aln_seq[start:end]
            if loc1:
                slicealn1 = self.seq[loc1[0]:loc1[1]].gappedByMap(map1)
                assert ori1 == slicealn1
            else:
                assert map1.length == len(ori1)
            
            # test by_align = False
            map2, loc2 = slice_cigar(self.cigar_text, start, end, by_align = False)
            slicealn2 = self.seq[start:end].gappedByMap(map2)
            ori2 = self.aln_seq[loc2[0]:loc2[1]]
            assert slicealn2 == ori2
    
    def test_CigarParser(self):
        """test without slice"""
        aln = CigarParser(self.seqs, self.cigars)
        assert aln == self.aln
        # test slice
        i = 1
        for start, end in self.slices:
            self.aln.getSeq("FAKE01").addFeature("annot%d"%i, "annot", [(start, end)])
            annot = self.aln.getAnnotationsFromAnySequence("annot%d"%i)
            slice_aln = aln.getRegionCoveringAll(annot).asOneSpan().getSlice()
            i += 1
            
            cmp_aln = CigarParser(self.seqs, self.cigars, sliced = True,
                                  ref_seqname = "FAKE01", start = start, end = end)
            assert cmp_aln == slice_aln 
    


if __name__ == '__main__':
    unittest.main()