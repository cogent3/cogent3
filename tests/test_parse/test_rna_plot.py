#!/usr/bin/env python
"""Tests for the RNAplot parser
"""
from __future__ import division
from cogent.util.unit_test import TestCase, main
import string
import re
from cogent.parse.rna_plot import get_sequence, get_coordinates, get_pairs,\
    RnaPlotParser

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.2"
__maintainer__ = "Jeremy Widmann"
__email__ = "jermey.widmann@colorado.edu"
__status__ = "Production"

class RnaPlotParserTests(TestCase):
    """Tests for RnaPlotParser.
    """

    def setUp(self):
        """Setup function for RnaPlotParser tests.
        """
        self.sequence_lines = SEQUENCE_LINES.split('\n')
        self.expected_seq = 'AAACCGCCUUU'
        
        self.coordinate_lines = COORDINATE_LINES.split('\n')
        self.expected_coords = [[92.500,92.500],\
                                [92.500,77.500],\
                                [92.500,62.500],\
                                [82.218,50.185],\
                                [85.577,34.497],\
                                [100.000,27.472],\
                                [114.423,34.497],\
                                [117.782,50.185],\
                                [107.500,62.500],\
                                [107.500,77.500],\
                                [107.500,92.500]]
                                
        self.pairs_lines = PAIRS_LINES.split('\n')
        self.expected_pairs = [[0,10],\
                               [1,9],\
                               [2,8]]
        
        self.rna_plot_lines = RNA_PLOT_FILE.split('\n')
    
    def test_get_sequence(self):
        """get_sequence should properly parse out sequence.
        """
        #test real data
        obs_seq = get_sequence(self.sequence_lines)
        self.assertEqual(obs_seq, self.expected_seq)
        
        #test empty list
        self.assertEqual(get_sequence([]),'')
        
    
    def test_get_coordinates(self):
        """get_coordinates should proplerly parse out coordinates.
        """
        obs_coords = get_coordinates(self.coordinate_lines)
        for (obs1, obs2), (exp1, exp2) in zip(obs_coords,self.expected_coords):
            self.assertFloatEqual(obs1,exp1)
            self.assertFloatEqual(obs2,exp2)
        
        #test empty list
        self.assertEqual(get_coordinates([]),[])

                            
    def test_get_pairs(self):
        """get_pairs should proplerly parse out pairs.
        """
        obs_pairs = get_pairs(self.pairs_lines)
        self.assertEqual(obs_pairs, self.expected_pairs)
        
        #test empty list
        self.assertEqual(get_pairs([]),[])

    
    def test_RnaPlotParser(self):
        """RnaPlotParser should properly parse full RNAplot postscript file.
        """
        obs_seq, obs_coords, obs_pairs = RnaPlotParser(self.rna_plot_lines)
        #test seq is correctly parsed
        self.assertEqual(obs_seq, self.expected_seq)
        
        #test coords are correctly parsed
        for (obs1, obs2), (exp1, exp2) in zip(obs_coords,self.expected_coords):
            self.assertFloatEqual(obs1,exp1)
            self.assertFloatEqual(obs2,exp2)
        
        #test pairs are correctly parsed
        self.assertEqual(obs_pairs, self.expected_pairs)
        
        #test empty list
        self.assertEqual(RnaPlotParser([]),('',[],[]))
        
        
        
SEQUENCE_LINES = """
/sequence (\
AAACCGCCUUU\
) def
/coor [
[92.500 92.500]
[92.500 77.500]
[92.500 62.500]
[82.218 50.185]
[85.577 34.497]
[100.000 27.472]
[114.423 34.497]
[117.782 50.185]
[107.500 62.500]
[107.500 77.500]
[107.500 92.500]
] def
/pairs [
[1 11]
[2 10]
[3 9]
] def

init

% switch off outline pairs or bases by removing these lines
drawoutline
drawpairs
drawbases
% show it
showpage
end
%%EOF
"""

COORDINATE_LINES = """
/coor [
[92.500 92.500]
[92.500 77.500]
[92.500 62.500]
[82.218 50.185]
[85.577 34.497]
[100.000 27.472]
[114.423 34.497]
[117.782 50.185]
[107.500 62.500]
[107.500 77.500]
[107.500 92.500]
] def
/pairs [
[1 11]
[2 10]
[3 9]
] def

init

% switch off outline pairs or bases by removing these lines
drawoutline
drawpairs
drawbases
% show it
showpage
end
%%EOF
"""

PAIRS_LINES = """
/pairs [
[1 11]
[2 10]
[3 9]
] def

init

% switch off outline pairs or bases by removing these lines
drawoutline
drawpairs
drawbases
% show it
showpage
end
%%EOF
"""

RNA_PLOT_FILE = """
%!PS-Adobe-3.0 EPSF-3.0
%%Creator: PS_dot.c,v 1.38 2007/02/02 15:18:13 ivo Exp $, ViennaRNA-1.8.2
%%CreationDate: Wed Apr 14 12:08:23 2010
%%Title: RNA Secondary Structure Plot
%%BoundingBox: 66 210 518 662
%%DocumentFonts: Helvetica
%%Pages: 1
%%EndComments

%Options: 
% to switch off outline pairs of sequence comment or
% delete the appropriate line near the end of the file

%%BeginProlog
/RNAplot 100 dict def
RNAplot begin
/fsize  14 def
/outlinecolor {0.2 setgray} bind def
/paircolor    {0.2 setgray} bind def
/seqcolor     {0   setgray} bind def
/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def
/min { 2 copy gt { exch } if pop } bind def
/max { 2 copy lt { exch } if pop } bind def
/drawoutline {
  gsave outlinecolor newpath
  coor 0 get aload pop 0.8 0 360 arc % draw 5' circle of 1st sequence
  currentdict /cutpoint known        % check if cutpoint is defined
  {coor 0 cutpoint getinterval
   {aload pop lineto} forall         % draw outline of 1st sequence
   coor cutpoint get aload pop
   2 copy moveto 0.8 0 360 arc       % draw 5' circle of 2nd sequence
   coor cutpoint coor length cutpoint sub getinterval
   {aload pop lineto} forall}        % draw outline of 2nd sequence
  {coor {aload pop lineto} forall}   % draw outline as a whole
  ifelse
  stroke grestore
} bind def
/drawpairs {
  paircolor
  0.7 setlinewidth
  [9 3.01] 9 setdash
  newpath
  pairs {aload pop
     coor exch 1 sub get aload pop moveto
     coor exch 1 sub get aload pop lineto
  } forall
  stroke
} bind def
% draw bases
/drawbases {
  [] 0 setdash
  seqcolor
  0
  coor {
    aload pop moveto
    dup sequence exch 1 getinterval cshow
    1 add
  } forall
  pop
} bind def

/init {
  /Helvetica findfont fsize scalefont setfont
  1 setlinejoin
  1 setlinecap
  0.8 setlinewidth
  72 216 translate
  % find the coordinate range
  /xmax -1000 def /xmin 10000 def
  /ymax -1000 def /ymin 10000 def
  coor {
      aload pop
      dup ymin lt {dup /ymin exch def} if
      dup ymax gt {/ymax exch def} {pop} ifelse
      dup xmin lt {dup /xmin exch def} if
      dup xmax gt {/xmax exch def} {pop} ifelse
  } forall
  /size {xmax xmin sub ymax ymin sub max} bind def
  72 6 mul size div dup scale
  size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div
  translate
} bind def
end
%%EndProlog
RNAplot begin
% data start here
/sequence (\
AAACCGCCUUU\
) def
/coor [
[92.500 92.500]
[92.500 77.500]
[92.500 62.500]
[82.218 50.185]
[85.577 34.497]
[100.000 27.472]
[114.423 34.497]
[117.782 50.185]
[107.500 62.500]
[107.500 77.500]
[107.500 92.500]
] def
/pairs [
[1 11]
[2 10]
[3 9]
] def

init

% switch off outline pairs or bases by removing these lines
drawoutline
drawpairs
drawbases
% show it
showpage
end
%%EOF
"""

#run if called from command-line
if __name__ == "__main__":
    main()
