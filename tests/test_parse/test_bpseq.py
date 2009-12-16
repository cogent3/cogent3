#!/usr/bin/env python
"""Provides Tests for BpseqParser and related functions.
"""
from cogent.util.unit_test import TestCase, main
from cogent.core.info import Info
from cogent.struct.knots import inc_order
from cogent.parse.bpseq import BpseqParseError, construct_sequence,\
    parse_header, parse_residues, MinimalBpseqParser, BpseqParser,\
    bpseq_specify_output

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class BpseqParserTests(TestCase):
    """Provides tests for BpseqParser and related functions"""

    def test_parse_header(self):
        """parse_header: should work on standard header"""
        
        h1 = ['Filename: d.16.b.E.coli.bpseq','Organism: Escherichia coli',\
            'Accession Number: J01695', 'Citation and related information'+\
            ' available at http://www.rna.icmb.utexas.edu']
        self.assertEqual(parse_header(h1),{'Filename':'d.16.b.E.coli.bpseq',\
            'Accession Number': 'J01695', 'Organism': 'Escherichia coli',\
            'Refs': {},'Citation':'http://www.rna.icmb.utexas.edu'})
        assert isinstance(parse_header(h1), Info)
        
        # lines without ':' are skipped
        h2 = ['Filename: d.16.b.E.coli.bpseq','Organism: Escherichia coli',\
            'Accession Number: J01695', 'Remark this is an interesting seq']
        exp = {'Filename':'d.16.b.E.coli.bpseq', 'Refs': {},\
            'Organism': 'Escherichia coli', 'Accession Number':'J01695'}
        self.assertEqual(parse_header(h2),exp)

    def test_construct_sequence(self):
        """construct_sequence: should return correct sequence or raise error
        """
        d = {0:'A',1:'C',2:'G',3:'U'}
        self.assertEqual(construct_sequence(d),'ACGU')
        # doesn't check residue identity
        d = {0:'A',1:'-',2:'R',3:'U'}
        self.assertEqual(construct_sequence(d),'A-RU')
        # error when sequence isn't continuous
        d = {0:'A',1:'C',2:'G',5:'U'}
        self.assertRaises(BpseqParseError, construct_sequence, d)
        # error when first index is not zero
        d = {1:'C',2:'G',3:'U',4:'A'}
        self.assertRaises(BpseqParseError, construct_sequence, d)

    def test_parse_residues(self):
        """parse_residues: should work on valid data
        """
        lines = RES_LINES.split('\n')
        exp_seq = 'UGGUAAUACGUUGCGAAGCC'
        exp_pairs = [(2,8),(3,7),(4,11),(5,10),(6,9),(12,18),(13,17)]
        self.assertEqual(parse_residues(lines, num_base=1,\
            unpaired_symbol='0'), (exp_seq, exp_pairs))
        
    def test_parse_residues_errors(self):
        """parse_residues: should raise BpseqParseErrors in several cases
        """
        not_all_lines = RES_LINES_NOT_ALL.split('\n')
        wrong_lines = RES_LINES_WRONG.split('\n')
        conflict_lines = RES_LINES_CONFLICT.split('\n')
        bp_conflict = RES_LINES_BP_CONFLICT.split('\n')
        self.assertRaises(BpseqParseError, parse_residues, not_all_lines,\
            num_base=1, unpaired_symbol='0')
        self.assertRaises(BpseqParseError, parse_residues, wrong_lines,\
            num_base=1, unpaired_symbol='0')
        self.assertRaises(BpseqParseError, parse_residues, conflict_lines,\
            num_base=1, unpaired_symbol='0')
        self.assertRaises(BpseqParseError, parse_residues, bp_conflict,\
            num_base=1, unpaired_symbol='0')

    def test_parse_residues_diff_base(self):
        """parse_residues: should work with diff base and unpaired_symbol"""
        lines = RES_LINES_DIFF_BASE.split('\n')
        exp_seq = 'CAGACU'
        exp_pairs = [(1,5),(2,4)]
        obs = parse_residues(lines, num_base=3, unpaired_symbol='xxx')
        self.assertEqual(obs, (exp_seq, exp_pairs))

    def test_MinimalBpseqParser(self):
        """MinimalBpseqParser: should separate lines correctly"""
        lines = ['Accesion: J01234', 'LABEL : label', '1 U 4', '2 A 10', 'xx',\
            'A B C D E']
        exp = {'HEADER': ['Accesion: J01234', 'LABEL : label'],\
            'SEQ_STRUCT': ['1 U 4', '2 A 10']}
        self.assertEqual(MinimalBpseqParser(lines), exp)

    def test_BpseqParser(self):
        """BpseqParser: should work on valid data, returning Vienna or Pairs
        """
        lines = RES_LINES_W_HEADER.split('\n')
        exp_seq = 'UGGUAAUACGUUGCGAAGCC'
        exp_pairs = [(2,8),(3,7),(4,11),(5,10),(6,9),(12,18),(13,17)]
        self.assertEqual(BpseqParser(lines),(exp_seq, exp_pairs))
        self.assertEqual(BpseqParser(lines)[0].Info,\
            {'Filename':'d.16.b.E.coli.bpseq',\
            'Accession Number': 'J01695', 'Organism': 'Escherichia coli',\
            'Refs': {},'Citation':'http://www.rna.icmb.utexas.edu'})

        # should work with different base
        lines = RES_LINES_DIFF_BASE.split('\n')
        exp_seq = 'CAGACU'
        exp_pairs = [(1,5),(2,4)]
        obs_seq, obs_pairs = BpseqParser(lines, num_base=3,\
            unpaired_symbol='xxx')
        self.assertEqual(obs_seq, exp_seq)
        self.assertEqual(obs_seq.Info, {'Refs':{}})
        self.assertEqual(obs_pairs, exp_pairs)

    def test_BpseqParser_errors(self):
        """BpseqParser: should skip lines in unknown format"""
        exp_seq = 'UGGUAAUACGUUGCGAAGCC'
        exp_vienna_m = '....(((..)))((...)).'
        exp_pairs = [(2,8),(3,7),(4,11),(5,10),(6,9),(12,18),(13,17)]
        
        #skips lines in unknown format
        lines = RES_LINES_UNKNOWN.split('\n')
        obs_seq, obs_pairs = BpseqParser(lines)

        self.assertEqual(obs_seq, exp_seq)
        self.assertEqual(obs_pairs, exp_pairs)
        self.assertEqual(obs_seq.Info,\
            {'Filename':'d.16.b.E.coli.bpseq',\
            'Accession Number': 'J01695', 'Organism': 'Escherichia coli',\
            'Refs': {},'Citation':'http://www.rna.icmb.utexas.edu'})

class ConvenienceFunctionTests(TestCase):
    """Tests for convenience functions"""

    def test_bpseq_specify_output(self):
        """bpseq_specify_output: different return values"""
        f = bpseq_specify_output
        lines = RES_LINES_W_HEADER.split('\n')
        exp_seq = 'UGGUAAUACGUUGCGAAGCC'
        exp_pairs = [(2,8),(3,7),(4,11),(5,10),(6,9),(12,18),(13,17)]
        exp_pairs_majority = [(4,11),(5,10),(6,9),(12,18),(13,17)]
        exp_pairs_first = [(2,8),(3,7),(12,18),(13,17)]
        exp_vienna_majority = '....(((..)))((...)).'

        self.assertEqual(f(lines),(exp_seq, exp_pairs))
        self.assertEqual(f(lines, remove_pseudo=True),\
            (exp_seq, exp_pairs_majority))
        self.assertEqual(f(lines, remove_pseudo=True, pseudoknot_function=inc_order),\
            (exp_seq, exp_pairs_first))
        self.assertEqual(f(lines, return_vienna=True),\
            (exp_seq, exp_vienna_majority))


RES_LINES=\
"""1 U 0
2 G 0
3 G 9
4 U 8
5 A 12
6 A 11
7 U 10
8 A 4
9 C 3
10 G 7
11 U 6
12 U 5
13 G 19
14 C 18
15 G 0
16 A 0
17 A 0
18 G 14
19 C 13
20 C 0"""


RES_LINES_NOT_ALL=\
"""1 U 0
2 G 0
3 G 0
6 A 0"""

RES_LINES_WRONG=\
"""1 U0
2 G 0
3 G 0
6 A 0"""

RES_LINES_CONFLICT=\
"""1 U 4
2 G 3
3 G 2
4 A 1
4 A 0"""

RES_LINES_BP_CONFLICT=\
"""1 U 0
2 G 4
3 G 4
4 C 2
5 A 0"""

RES_LINES_DIFF_BASE=\
"""3 C xxx
4 A 8
5 G 7
6 A xxx
7 C 5
8 U 4"""

RES_LINES_W_HEADER=\
"""Filename: d.16.b.E.coli.bpseq
Organism: Escherichia coli
Accession Number: J01695
Citation and related information available at http://www.rna.icmb.utexas.edu
1 U 0
2 G 0
3 G 9
4 U 8
5 A 12
6 A 11
7 U 10
8 A 4
9 C 3
10 G 7
11 U 6
12 U 5
13 G 19
14 C 18
15 G 0
16 A 0
17 A 0
18 G 14
19 C 13
20 C 0"""

RES_LINES_UNKNOWN=\
"""Filename: d.16.b.E.coli.bpseq
Organism: Escherichia coli
Accession Number: J01695
Citation and related information available at http://www.rna.icmb.utexas.edu
1 U 0
2 G 0
3 G 9
UNKNOWN LINE
4 U 8
5 A 12
6 A 11
7 U 10
8 A 4
9 C 3
10 G 7
11 U 6
12 U 5
13 G 19
14 C 18
15 G 0
16 A 0
17 A 0
18 G 14
19 C 13
20 C 0"""

#run if called from command-line
if __name__ == "__main__":
    main()
