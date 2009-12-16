#!/usr/bin/env python
"""Tests for the RNAfold dot plot parser
"""
from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.parse.rna_fold import *

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

class RnaFoldParserTests(TestCase):
    """Tests for RnaFoldParser.
    """

    def setUp(self):
        """Setup function for RnaFoldParser tests.
        """
        self.rna_fold_parser_results = ('ACGUGCUAG',
            [(1, 7, float(0.01462)),
            (2, 9, float(0.11118)),
            (3, 7, float(0.00985)),
            (4, 8, float(0.01005)),
            (4, 9, float(0.01586))])
        
        self.sequence_lines = ['/sequence line before where sequences start\n',
            '  ACCUGUCUAUCGCUGC&*$#@(*\n',
            'ACGGUUAUAUUAUCUCUG\\\n',
            ') end of sequence\n']
        self.sequence_lines_empty = ['/sequence \n',
            '\n',
            ')\n']
        self.index_lines = ['unimportant line',
            '1 3 0.332 ubox',
            '1 4 0.003 ubox']
        self.index_lines_no_ubox = ['unimportant line',
            '1 2 0.432 u box',
            '1 4 0.32 ubo x']
    
    def test_getSequence(self):
        self.assertEqual(getSequence(self.sequence_lines),
            'ACCUGUCUAUCGCUGC&*$#@(*ACGGUUAUAUUAUCUCUG')
        self.assertEqual(getSequence(self.sequence_lines_empty),'')
    
    def test_getIndices(self):
        self.assertEqual(getIndices(self.index_lines),[(1,3,float(0.332)),
                                                        (1,4,float(0.003))])
        self.assertEqual(getIndices(self.index_lines_no_ubox),[])
    
    def test_RnaFoldParser(self):
        self.assertEqual(RnaFoldParser([]), ('',[]))
        self.assertEqual(RnaFoldParser(RNA_FOLD_RESULTS),
            self.rna_fold_parser_results)


RNA_FOLD_RESULTS = ['%!PS-Adobe-3.0 EPSF-3.0\n',
 '%%Title: RNA DotPlot\n',
 '%%Creator: PS_dot.c,v 1.24 2003/08/07 09:01:00 ivo Exp $, ViennaRNA-1.5\n',
 '%%CreationDate: Fri Oct  8 13:15:01 2004\n',
 '%%BoundingBox: 66 211 518 662\n',
 '%%DocumentFonts: Helvetica\n',
 '%%Pages: 1\n',
 '%%EndComments\n',
 '\n',
 '%Options: \n',
 '%This file contains the square roots of the base pair probabilities in the form\n',
 '% i  j  sqrt(p(i,j)) ubox\n',
 '100 dict begin\n',
 '\n',
 '/logscale false def\n',
 '\n',
 '%delete next line to get rid of title\n',
 '270 665 moveto /Helvetica findfont 14 scalefont setfont (dot.ps) show\n',
 '\n',
 '/lpmin {\n',
 '   1e-05 log  % log(pmin) only probs>pmin will be shown\n',
 '} bind def\n',
 '\n',
 '/box { %size x y box - draws box centered on x,y\n',
 '   2 index 0.5 mul add            % x += 0.5\n',
 '   exch 2 index 0.5 mul add exch  % x += 0.5\n',
 '   newpath\n',
 '   moveto\n',
 '   dup neg   0 rlineto\n',
 '   dup neg   0 exch rlineto\n',
 '             0 rlineto\n',
 '   closepath\n',
 '   fill\n',
 '} bind def\n',
 '\n',
 '/sequence { (\\\n',
 'ACGUGCUAG\\\n',
 ') } def\n',
 '/len { sequence length } bind def\n',
 '\n',
 '/ubox {\n',
 '   logscale {\n',
 '      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n',
 '   } if\n',
 '   3 1 roll\n',
 '   exch len exch sub 1 add box\n',
 '} bind def\n',
 '\n',
 '/lbox {\n',
 '   3 1 roll\n',
 '   len exch sub 1 add box\n',
 '} bind def\n',
 '\n',
 '72 216 translate\n',
 '72 6 mul len 1 add div dup scale\n',
 '/Helvetica findfont 0.95 scalefont setfont\n',
 '\n',
 '% print sequence along all 4 sides\n',
 '[ [0.7 -0.3 0 ]\n',
 '  [0.7 0.7 len add 0]\n',
 '  [0.7 -0.2 90]\n',
 '  [-0.3 len sub 0.7 len add -90]\n',
 '] {\n',
 '  gsave\n',
 '    aload pop rotate translate\n',
 '    0 1 len 1 sub {\n',
 '     dup 0 moveto\n',
 '     sequence exch 1 getinterval\n',
 '     show\n',
 '    } for\n',
 '  grestore\n',
 '} forall\n',
 '\n',
 '0.5 dup translate\n',
 '% draw diagonal\n',
 '0.04 setlinewidth\n',
 '0 len moveto len 0 lineto stroke \n',
 '\n',
 '%draw grid\n',
 '0.01 setlinewidth\n',
 'len log 0.9 sub cvi 10 exch exp  % grid spacing\n',
 'dup 1 gt {\n',
 '   dup dup 20 div dup 2 array astore exch 40 div setdash\n',
 '} { [0.3 0.7] 0.1 setdash } ifelse\n',
 '0 exch len {\n',
 '   dup dup\n',
 '   0 moveto\n',
 '   len lineto \n',
 '   dup\n',
 '   len exch sub 0 exch moveto\n',
 '   len exch len exch sub lineto\n',
 '   stroke\n',
 '} for\n',
 '0.5 neg dup translate\n',
 '\n',
 '1 7 0.01462 ubox\n',
 '2 9 0.11118 ubox\n',
 '3 7 0.00985 ubox\n',
 '4 8 0.01005 ubox\n',
 '4 9 0.01586 ubox\n',
 'showpage\n',
 'end\n',
 '%%EOF\n']
     
#run if called from command-line
if __name__ == "__main__":
    main()
