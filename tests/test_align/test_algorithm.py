#!/usr/bin/env python

from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.align.algorithm import ScoreCell, MatchScorer, equality_scorer,\
    default_gap, default_gap_symbol, ScoreMatrix, NeedlemanWunschMatrix, \
    SmithWatermanMatrix, nw_align, sw_align
from copy import copy, deepcopy

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

"""Tests of the cogent.align.algorithm module.
"""

class ScoreCellTests(TestCase):
    """Tests for ScoreCell class.
    """
    
    def setUp(self):
        """Setup for ScoreCell tests."""
        self.score_cell = ScoreCell()
    
    def test_init(self):
        """Tests for ScoreCell __init__ function."""
        #Test for empty score cell
        self.assertEqual(self.score_cell.Score,0)
        self.assertEqual(self.score_cell.Pointer,None)
    
    def test_update(self):
        """Tests for ScoreCell update function."""
        #Test tie up wins
        self.score_cell.update(1,1,1)
        self.assertEqual(self.score_cell.Score,1)
        self.assertEqual(self.score_cell.Pointer,"up")
        #Test tie diag wins
        self.score_cell.update(3,1,3)
        self.assertEqual(self.score_cell.Score,3)
        self.assertEqual(self.score_cell.Pointer,"diag")
        #Test up
        self.score_cell.update(1,1,1)
        self.assertEqual(self.score_cell.Score,1)
        self.assertEqual(self.score_cell.Pointer,"up")
        #Test diag
        self.score_cell.update(3,1,2)
        self.assertEqual(self.score_cell.Score,3)
        self.assertEqual(self.score_cell.Pointer,"diag")
        #Test left
        self.score_cell.update(1,2,3)
        self.assertEqual(self.score_cell.Score,3)
        self.assertEqual(self.score_cell.Pointer,"left")

class MatchScorerTests(TestCase):
    """Tests for MatchScorer function.
    """
    def setUp(self):
        """Setup for MatchScorer tests."""
        self.equality_scorer = MatchScorer(1,-1)
        self.mismatch_worse_scorer = MatchScorer(1,-2)
        self.match_worse_scorer = MatchScorer(-1,1)
    
    def test_scorer(self):
        """Tests for MatchScorer function."""
        #Test equality_scorer
        self.assertEqual(self.equality_scorer('A','A'),1)
        self.assertEqual(self.equality_scorer('A','B'),-1)
        #Test mismatch_worse_scorer
        self.assertEqual(self.mismatch_worse_scorer('A','A'),1)
        self.assertEqual(self.mismatch_worse_scorer('A','B'),-2)
        #Test match_worse_scorer
        self.assertEqual(self.match_worse_scorer('A','A'),-1)
        self.assertEqual(self.match_worse_scorer('A','B'),1)
        

class ScoreMatrixTests(TestCase):
    """Tests for ScoreCell class.
    """
    
    def setUp(self):
        """Setup for ScoreMatrix tests."""
        self.score_matrix = ScoreMatrix('AACGU','ACAGU')
        def fill():
            pass
        def traceback():
            pass
        self.score_matrix.fill = fill
        self.score_matrix.traceback = traceback
        
        self.empty_score_matrix = ScoreMatrix('','')
    def test_init(self):
        """Tests for ScoreMatrix __init__ function."""
        #Test empty ScoreMatrix
        self.empty_score_matrix = ScoreMatrix('','')
        self.assertEqual(self.empty_score_matrix.First,'')
        self.assertEqual(self.empty_score_matrix.Second,'')
        self.assertEqual(self.empty_score_matrix.Cols,1)
        self.assertEqual(self.empty_score_matrix.Rows,1)
        self.assertEqual(self.empty_score_matrix.GapScore,-1)
        self.assertEqual(self.empty_score_matrix.GapSymbol,'-')
        self.assertEqual(self.empty_score_matrix.FirstAlign,[])
        self.assertEqual(self.empty_score_matrix.SecondAlign,[])
    
    def test_str(self):
        """Tests for ScoreMatrix __str__ function."""
        #Test empty ScoreMatrix
        self.assertEqual(self.empty_score_matrix.__str__(),"Empty Score Matrix")
        
        #Test full ScoreMatrix
        self.assertEqual(self.score_matrix.__str__(),\
"""\t\tA\tA\tC\tG\tU\n\t0\t0\t0\t0\t0\t0\nA\t0\t0\t0\t0\t0\t0\nC\t0\t0\t0\t0\t0\t0\nA\t0\t0\t0\t0\t0\t0\nG\t0\t0\t0\t0\t0\t0\nU\t0\t0\t0\t0\t0\t0""")

    def test_alignment(self):
        """Tests for ScoreMatrix alignment function."""
        #Should not align since ScoreMatrix base object does not have fill()
        #or traceback() methods.  For testing purposes fill() and traceback()
        #for self.score_matrix do nothing.
        self.assertEqual(self.score_matrix.alignment(),('',''))

class NeedlemanWunschMatrixTests(TestCase):
    """Tests for NeedlemanWunschMatrix class.
    """
    
    def setUp(self):
        """Setup for NeedlemanWunschMatrix tests."""
        self.nw_matrix = NeedlemanWunschMatrix('ACGU','CAGU')
        #Since _init_first_row, _init_first_col and _init_first_cell are
        #automatically called when NeedlemanWunschMatrix is initialized,
        #need to change all elements in nw_matrix_empty to 0 to test
        #that each _init works.
        self.nw_matrix_empty = \
            NeedlemanWunschMatrix('ACGU','CAGU')
        for i in range(len(self.nw_matrix_empty)):
            for j in range(len(self.nw_matrix_empty[i])):
                self.nw_matrix_empty[i][j]=ScoreCell(0)
    
    def test_init_first_row(self):
        """Tests for NeedlemanWunschMatrix _init_first_row function."""
        nw_matrix_empty = copy(self.nw_matrix_empty)
        #Init first row
        nw_matrix_empty._init_first_row()
        #matrix after first row init
        matrix_scores_first_row_init = [[0,-1,-2,-3,-4],
                                        [0,0,0,0,0],
                                        [0,0,0,0,0],
                                        [0,0,0,0,0],
                                        [0,0,0,0,0]]
        for i in range(len(nw_matrix_empty)):
            for j in range(len(nw_matrix_empty[i])):
                self.assertEqual(nw_matrix_empty[i][j].Score,\
                                matrix_scores_first_row_init[i][j])
    
    def test_init_first_col(self):
        """Tests for NeedlemanWunschMatrix _init_first_col function."""
        nw_matrix_empty = copy(self.nw_matrix_empty)
        #Init first row
        nw_matrix_empty._init_first_col()
        #matrix after first row init
        matrix_scores_first_row_init = [[0,0,0,0,0],
                                        [-1,0,0,0,0],
                                        [-2,0,0,0,0],
                                        [-3,0,0,0,0],
                                        [-4,0,0,0,0]]
        for i in range(len(nw_matrix_empty)):
            for j in range(len(nw_matrix_empty[i])):
                self.assertEqual(nw_matrix_empty[i][j].Score,\
                                matrix_scores_first_row_init[i][j])
    
    def test_init_first_cell(self):
        """Tests for NeedlemanWunschMatrix _init_first_cell function."""
        nw_matrix_empty = copy(self.nw_matrix_empty)
        #Init first row
        nw_matrix_empty._init_first_cell()
        self.assertEqual(nw_matrix_empty[0][0].Score,0)
    
    def test_initialized_matrix(self):
        """Tests for NeedlemanWunschMatrix after full initialization."""
        matrix_scores_first_row_init = [[0,-1,-2,-3,-4],
                                        [-1,0,0,0,0],
                                        [-2,0,0,0,0],
                                        [-3,0,0,0,0],
                                        [-4,0,0,0,0]]
        for i in range(len(self.nw_matrix)):
            for j in range(len(self.nw_matrix[i])):
                self.assertEqual(self.nw_matrix[i][j].Score,\
                                matrix_scores_first_row_init[i][j])
    
    def test_fill(self):
        """Tests for NeedlemanWunschMatrix fill function."""
        filled_nw_matrix = copy(self.nw_matrix)
        filled_nw_matrix.fill()
        matrix_scores_filled = [[0,-1,-2,-3,-4],
                                [-1,-1,0,-1,-2],
                                [-2,0,-1,-1,-2],
                                [-3,-1,-1,0,-1],
                                [-4,-2,-2,-1,1]]
        matrix_pointers_filled = [[None,'left','left','left','left'],
                                  ['up','diag','diag','left','left'],
                                  ['up','diag','up','diag','diag'],
                                  ['up','up','diag','diag','left'],
                                  ['up','up','up','up','diag']]
        
        self.assertEqual(filled_nw_matrix.Filled,True)
        
        for i in range(len(filled_nw_matrix)):
            for j in range(len(filled_nw_matrix[i])):
                self.assertEqual(filled_nw_matrix[i][j].Score,\
                                    matrix_scores_filled[i][j])
                self.assertEqual(filled_nw_matrix[i][j].Pointer,\
                                    matrix_pointers_filled[i][j])
        self.assertEqual(filled_nw_matrix.MaxScore,(1,4,4))


    def test_traceback(self):
        """Tests for NeedlemanWunschMatrix traceback function."""
        self.nw_matrix.traceback()
        self.assertEqual(self.nw_matrix.FirstAlign,['A','C','-','G','U'])
        self.assertEqual(self.nw_matrix.SecondAlign,['-','C','A','G','U'])

class SmithWatermanMatrixTests(TestCase):
    """Tests for SmithWatermanMatrix class.
    """
    
    def setUp(self):
        """Setup for SmithWatermanMatrix tests."""
        self.sw_matrix = SmithWatermanMatrix('ACGU','CAGU')
            
    def test_fill(self):
        """Tests for SmithWatermanMatrix fill function."""
        filled_sw_matrix = copy(self.sw_matrix)
        filled_sw_matrix.fill()
        matrix_scores_filled = [[0,0,0,0,0],
                                [0,0,1,0,0],
                                [0,1,0,0,0],
                                [0,0,0,1,0],
                                [0,0,0,0,2]]
        matrix_pointers_filled = [[None,None,None,None,None],
                                  [None,None,'diag',None,None],
                                  [None,'diag',None,None,None],
                                  [None,None,None,'diag',None],
                                  [None,None,None,None,'diag']]
        
        
        for i in range(len(filled_sw_matrix)):
            for j in range(len(filled_sw_matrix[i])):
                self.assertEqual(filled_sw_matrix[i][j].Score,\
                                    matrix_scores_filled[i][j])
                self.assertEqual(filled_sw_matrix[i][j].Pointer,\
                                    matrix_pointers_filled[i][j])
        self.assertEqual(filled_sw_matrix.MaxScore,(2,4,4))

    def test_traceback(self):
        """Tests for SmithWatermanMatrix traceback function."""
        self.sw_matrix.fill()
        self.sw_matrix.traceback()
        self.assertEqual(self.sw_matrix.FirstAlign,['G','U'])
        self.assertEqual(self.sw_matrix.SecondAlign,['G','U'])

class NwAlignTests(TestCase):
    """Tests for nw_align fuction.
    """
    def test_nw_align_empty(self):
        """Tests for nw_align function."""
        (first,second),score = nw_align('','',return_score=True)
        self.assertEqual(first,'')
        self.assertEqual(second,'')
        self.assertEqual(score,0)

    def test_nw_align(self):
        """Tests for nw_align function."""
        (first,second),score = nw_align('ACGU','CAGU',return_score=True)
        self.assertEqual(first,'AC-GU')
        self.assertEqual(second,'-CAGU')
        self.assertEqual(score,1)

class SwAlignTests(TestCase):
    """Tests for sw_align function.
    """
    def test_sw_align_empty(self):
        """Tests for sw_align function."""
        (first,second),score = sw_align('','',return_score=True)
        self.assertEqual(first,'')
        self.assertEqual(second,'')
        self.assertEqual(score,0)

    def test_sw_align(self):
        """Tests for sw_align function."""
        (first,second),score = sw_align('ACGU','CAGU',return_score=True)
        self.assertEqual(first,'GU')
        self.assertEqual(second,'GU')
        self.assertEqual(score,2)
#run if called from command-line
if __name__ == "__main__":
    main()
