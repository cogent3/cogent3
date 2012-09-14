#!/usr/bin/env python
"""Code for performing alignments by Needleman-Wunsch and Smith-Waterman.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

class ScoreCell(object):
    """Cell in a ScoreMatrix object. Contains score and pointer."""
    __slots__ = ['Score', 'Pointer']
    def __init__(self, Score=0, Pointer=None):
        """Returns new ScoreCell object.

        Score should be a number.
        Pointer should be a direction (typically None, "up", "diag", or "left")
        """
        self.Score = Score
        self.Pointer = Pointer

    def update(self, diag, up, left):
        """Updates the ScoreCell object with scores for its three neighbors.

        Finds max of (diag, up, and left), and sets the score and pointer
        accordingly. In case of tie, picks up (i.e. tries to minimize
        gaps since the first sequence is the shortest), then diag, then left.
        """
        best = max(diag, up, left)
        if best == up:
            self.Score = up
            self.Pointer = "up"
        elif best == diag:
            self.Score = diag
            self.Pointer = "diag"
        else:
            self.Score = left
            self.Pointer = "left"

def MatchScorer(match, mismatch):
    """Factory function that returns a score function set to match and mismatch.

    match and mismatch should both be numbers. Typically, match should be 
    positive and mismatch should be negative.

    Resulting function has signature f(x,y) -> number.
    """
    def scorer(x, y):
        if x == y:
            return match
        else:
            return mismatch
    return scorer

equality_scorer = MatchScorer(1, -1)
default_gap = -1
default_gap_symbol = '-'

class ScoreMatrix(list):
    """Matrix that contains (score, pointer) pairs for sequence alignment."""
    
    def __init__(self, First, Second, Score=equality_scorer, \
        GapScore=default_gap, GapSymbol=default_gap_symbol):
        """Returns new ScoreMatrix object, initialized but not filled.

        Calls internal methods to initialize first row, column and cell.

        First and Second are the two sequences that will be aligned. Columns
        correspond to positions in First; rows correspond to positions in 
        Second.
        
        Score is a scoring function that takes two corresponding items in
        First and Second, and returns a number that can be compared and added.

        Gap is the score for inserting a gap.

        Note that the matrix is one row and one column bigger than the lengths
        of the sequences, since the first row and first column contain 
        initialization data.
        """
        self.First = First          #first sequence
        self.Second = Second        #second sequence
        self.Cols = len(First) + 1  #need to add 1 for the start col
        self.Rows = len(Second) + 1 #need to add 1 for the start row
        self.Score = Score          #scoring function: f(x,y) = num
        self.GapScore = GapScore    #gap penalty
        self.GapSymbol = GapSymbol  #symbol for gaps
        self.Filled = False         #whether the matrix has been filled
        self.FirstAlign = []        #first sequence, aligned, as list
        self.SecondAlign = []       #second sequence, aligned, as list
        for r in range(self.Rows):
            self.append([ScoreCell() for c in range(self.Cols)])
        self._init_first_row()
        self._init_first_col()
        self._init_first_cell()
        #print "INITIALIZED WITH:"
        #print self

    def __str__(self):
        """Prints score matrix, including labels for row and column."""
        first = self.First
        second = self.Second
        if not first or not second:
            return "Empty Score Matrix"
        #add first sequence as first row: skip 2 positions for row label and
        #for the first column, which is initialized to default values
        rows = ['\t\t' + '\t'.join(self.First)]
        for index, row in enumerate(self):
            #if it's not the first row, add the appropriate char from the seq
            if index:
                curr = second[index-1]
            else:
                curr = ''
            #add the row containing the char and the data
            rows.append('%s\t'%curr + '\t'.join([str(i.Score) for i in row]))
        return '\n'.join(rows)

    def _init_first_row(self):
        """Hook for matrices that need to initialize the first row."""
        pass

    def _init_first_col(self):
        """Hook for matrices that need to initialize the first column."""
        pass

    def _init_first_cell(self):
        """Hook for matrices that need to initialize the first cell."""
        pass

    def alignment(self):
        """Returns alignment of first and second sequences.

        Calls fill() and traceback() if necessary.
        
        Converts the aligned versions to the same class as the originals.
        """
        seq1, seq2 = self.First, self.Second
        aln1, aln2 = self.FirstAlign, self.SecondAlign
        if not aln1 or not aln2:
            self.fill()
            self.traceback()
            aln1, aln2 = self.FirstAlign, self.SecondAlign
        #Remember to return sequences that are the correct class. If the class
        #subclasses str, you probably want ''.join rather than str to feed into
        #the constructor, since str() on a list prints the brackets and commas.
        if isinstance(seq1, str):
            first_result = seq1.__class__(''.join(aln1))
        else:
            first_result = seq1.__class__(aln1)
            
        if isinstance(seq2, str):
            second_result = seq2.__class__(''.join(aln2))
        else:
            second_result = seq2.__class__(aln2)

        return first_result, second_result
    

class NeedlemanWunschMatrix(ScoreMatrix):
    """Score matrix for global alignment using Needleman-Wunsch."""

    def _init_first_row(self):
        """First row initialized with index * gap penalty."""
        gap = self.GapScore
        self[0] = [ScoreCell(gap * i, 'left') for i in range(self.Cols)]
        #print "AFTER FIRST ROW:"
        #print self

    def _init_first_col(self):
        """First column initialized with index * gap penalty."""
        gap = self.GapScore
        for index, row in enumerate(self):
            row[0] = ScoreCell(gap * index, 'up')
        #print "AFTER FIRST COL:"
        #print self

    def _init_first_cell(self):
        """First cell initialized to 0 with no pointer."""
        self[0][0] = ScoreCell(0)
        #print "AFTER FIRST CELL:"
        #print self

    def fill(self):
        """Fills each cell with its best score and the direction of the next.

        For moving up or left, calculates the gap penalty plus the score of the
        cell. For moving diagonally, calculates the match/mismatch score for
        the corresponding positions in the sequence plus the score of the cell.

        Performs all calculations in place.
        """
        score = self.Score
        gap = self.GapScore
        seq_1, seq_2 = self.First, self.Second
        curr_row = self[0]
        for row_i in range(1, self.Rows):
            prev_row = curr_row
            curr_row = self[row_i]
            curr_cell = curr_row[0]
            for col_i in range(1, self.Cols):
                prev_cell = curr_cell
                curr_cell = curr_row[col_i]
                #remember to subtract 1 to find corresponding pos in seq
                diag_score = score(seq_1[col_i-1], seq_2[row_i-1]) + \
                    prev_row[col_i-1].Score
                up_score = gap + prev_row[col_i].Score
                left_score = gap + prev_cell.Score
                #print 'row: %s col: %s '%(row_i,col_i), \
                #    diag_score,up_score,left_score
                curr_cell.update(diag_score, up_score, left_score)
        self.Filled = True
        self.MaxScore = (self[-1][-1].Score, len(self)-1, len(self[0])-1)
        #print "FINISHED FILL"
        #print self

    def traceback(self):
        """Returns the optimal alignment as a (first, second) tuple w/ gaps.

        Follows the pointers assigned in fill(), inserting gaps whenever the
        movement is not diagonal.
        """
        if not self.Filled:
            self.fill()
        align_1 = []
        align_2 = []
        seq_1 = self.First
        seq_2 = self.Second
        curr_row = self.Rows - 1    #subtract 1 to use length as index
        curr_col = self.Cols - 1    #subtract 1 to use length as index
        p = self[curr_row][curr_col].Pointer
        while 1:
            #print "ROW: %s, COL: %s" % (curr_row, curr_col),p
            if p == 'diag':
                align_1.append(seq_1[curr_col-1])
                align_2.append(seq_2[curr_row-1])
                curr_row -= 1
                curr_col -= 1
            elif p == 'left':
                align_1.append(seq_1[curr_col-1])
                align_2.append('-')
                curr_col -= 1
            elif p == 'up':
                align_1.append('-')
                align_2.append(seq_2[curr_row-1])
                curr_row -= 1
            else:
                break
            p = self[curr_row][curr_col].Pointer
        align_1.reverse()
        align_2.reverse()
        self.FirstAlign, self.SecondAlign = align_1, align_2

class SmithWatermanMatrix(ScoreMatrix):

    def fill(self):
        max_row, max_col, max_score = 0, 0, 0
        score = self.Score
        gap = self.GapScore
        seq_1, seq_2 = self.First, self.Second
        curr_row = self[0]
        for row_i in range(1, self.Rows):
            prev_row = curr_row
            curr_row = self[row_i]
            curr_cell = curr_row[0]
            for col_i in range(1, self.Cols):
                prev_cell = curr_cell
                curr_cell = curr_row[col_i]
                #remember to subtract 1 to find corresponding pos in seq
                diag_score = score(seq_1[col_i-1], seq_2[row_i-1]) + \
                    prev_row[col_i-1].Score
                up_score = gap + prev_row[col_i].Score
                left_score = gap + prev_cell.Score
                if max(up_score, left_score, diag_score) <= 0:
                    continue #leave uninitialized if scores all below threshold
                curr_cell.update(diag_score, up_score, left_score)
                if curr_cell.Score > max_score:
                    max_score = curr_cell.Score
                    max_row = row_i
                    max_col = col_i
        self.MaxScore = (max_score, max_row, max_col)
                    
        #print "FINISHED FILL"
        #print self


    def traceback(self):
        align_1 = []
        align_2 = []
        seq_1 = self.First
        seq_2 = self.Second
        max_score, curr_row, curr_col = self.MaxScore
        p = self[curr_row][curr_col].Pointer
        while 1:
            #print "ROW: %s, COL: %s" % (curr_row, curr_col)
            if p == 'diag':
                align_1.append(seq_1[curr_col-1])
                align_2.append(seq_2[curr_row-1])
                curr_row -= 1
                curr_col -= 1
            elif p == 'left':
                align_1.append(seq_1[curr_col-1])
                align_2.append('-')
                curr_col -= 1
            elif p == 'up':
                align_1.append('-')
                align_2.append(seq_2[curr_row-1])
                curr_row -= 1
            else:
                break
            p = self[curr_row][curr_col].Pointer
        align_1.reverse()
        align_2.reverse()
        self.FirstAlign, self.SecondAlign = align_1, align_2

def nw_align(seq1, seq2, scorer=equality_scorer, gap=default_gap, return_score=False):
    """Returns globally optimal alignment of seq1 and seq2."""
    N = NeedlemanWunschMatrix(seq1, seq2, scorer, gap)
    if return_score:
        return N.alignment(), N.MaxScore[0]
    else:
        return N.alignment()

def sw_align(seq1, seq2, scorer=equality_scorer, gap=default_gap, return_score=False):
    """Returns locally optimal alignment of seq1 and seq2."""
    S = SmithWatermanMatrix(seq1, seq2, scorer, gap)
    if return_score:
        return S.alignment(), S.MaxScore[0]
    else:
        return S.alignment()

def demo(seq1, seq2):
    result = []
    result.append("Global alignment:")
    result.extend(nw_align(seq1, seq2))
    result.append("Local alignment:")
    result.extend(sw_align(seq1, seq2))
    return "\n".join(result)

#initialization

if __name__ == '__main__':
    from sys import argv
    print demo(argv[1], argv[2])
