#!/usr/bin/env python
from __future__ import division
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import app_path
from cogent.format.pdb_color import get_aligned_muscle, make_color_list, \
    ungapped_to_pdb_numbers, get_matching_chains, get_chains, \
    get_best_muscle_hits, chains_to_seqs, align_subject_to_pdb

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

"""Tests of the pdb_color module.

Owner: Jeremy Widmann jeremy.widmann@colorado.edu

Revision History:

October 2006 Jeremy Widmann: File created
"""

MUSCLE_PATH = app_path('muscle')

class PdbColorTests(TestCase):
    """Tests for pdb_color functions.
    """
    def setUp(self):
        """Setup for pdb_color tests."""
        #Nucleotide test data results
        self.test_pdb_chains_1 = {'A': [(1, 'G'), (2, 'C'), (3, 'C'), (4, 'A'),
                                      (5, 'C'), (6, 'C'), (7, 'C'), (8, 'U'),
                                      (9, 'G')],
                                'B': [(10, 'C'), (11, 'A'), (12, 'G'),
                                      (13, 'G'), (14, 'G'), (15, 'U'),
                                      (16, 'C'), (17, 'G'), (18, 'G'),
                                      (19, 'C')]}
        self.ungapped_to_pdb_1 = {'A': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6,
                                        6: 7, 7: 8, 8: 9},
                                  'B': {0: 10, 1: 11, 2: 12, 3: 13, 4: 14,
                                        5: 15, 6: 16, 7: 17, 8: 18, 9: 19}}
        
        self.test_pdb_seqs_1 = {'A': 'GCCACCCUG', 'B': 'CAGGGUCGGC'}
        self.test_pdb_types_1 = {'A': 'Nucleotide', 'B': 'Nucleotide'}
        #Protein test data results
        self.test_pdb_chains_2 = {'A': [(1, 'ALA'), (2, 'PRO'), (3, 'ILE'),
                                        (4, 'LYS'), (5, 'VAL'), (6, 'GLY'),
                                        (7, 'ASP'), (8, 'ALA')]}
        self.ungapped_to_pdb_2 = {'A': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6,
                                        6: 7, 7: 8}}
        self.test_pdb_seqs_2 = {'A': 'APIKVGDA'}
        self.test_pdb_types_2 = {'A': 'Protein'}
    
    def test_get_aligned_muscle(self):
        """Tests for get_aligned_muscle function.
        """
        if not MUSCLE_PATH:
            return 'skipping test'
        seq1 = 'ACCUG'
        seq2 = 'ACGGUG'
        seq1_aligned_known = 'AC-CUG'
        seq2_aligned_known = 'ACGGUG'
        frac_same_known = 4/5.0
        
        seq1_aln, seq2_aln, frac_same = get_aligned_muscle(seq1,seq2)

        self.assertEqual(seq1_aln,seq1_aligned_known)
        self.assertEqual(seq2_aln,seq2_aligned_known)
        self.assertEqual(frac_same,frac_same_known)
            
    def test_get_chains_nucleotide(self):
        """Tests for get_chains function using nucleotide pdb lines.
        """
        chains_nuc = get_chains(TEST_PDB_STRING_1.split('\n'))
        self.assertEqual(chains_nuc, self.test_pdb_chains_1)
    
    def test_get_chains_protein(self):
        """Tests for get_chains function using protein pdb lines.
        """    
        chains_prot = get_chains(TEST_PDB_STRING_2.split('\n'))
        self.assertEqual(chains_prot, self.test_pdb_chains_2)
    
    def test_ungapped_to_pdb_nucleotide(self):
        """Tests for ungapped_to_pdb function using nucleotide pdb chains.
        """
        for k,v in self.test_pdb_chains_1.items():
            self.assertEqual(ungapped_to_pdb_numbers(v),\
                            self.ungapped_to_pdb_1[k])                        
    
    def test_ungapped_to_pdb_protein(self):
        """Tests for ungapped_to_pdb function using protein pdb chains.
        """
        for k,v in self.test_pdb_chains_2.items():
            self.assertEqual(ungapped_to_pdb_numbers(v),\
                            self.ungapped_to_pdb_2[k])
    
    def test_chains_to_seqs_nucleotide(self):
        """Tests for chains_to_seqs function using nucleotide pdb chains.
        """
        seqs, seqtypes = chains_to_seqs(self.test_pdb_chains_1)
        self.assertEqual(seqs, self.test_pdb_seqs_1)
        self.assertEqual(seqtypes, self.test_pdb_types_1)
        
    
    def test_chains_to_seqs_protein(self):
        """Tests for chains_to_seqs function using protein pdb chains.
        """
        seqs, seqtypes = chains_to_seqs(self.test_pdb_chains_2)
        self.assertEqual(seqs, self.test_pdb_seqs_2)
        self.assertEqual(seqtypes, self.test_pdb_types_2)
    
    def test_get_best_muscle_hits(self):
        """Tests for get_best_muscle_hits function.
        """
        if not MUSCLE_PATH:
            return 'skipping test'
        subject_seq =  'AACCGGUU'
        query_aln = {1:'CCCCCCCC',
                     2:'GGGGGGGG',
                     3:'AAGGGGUU',
                     4:'AACCGGGU'}
        res_20 = {1:'CCCCCCCC',
                  2:'GGGGGGGG',
                  3:'AAGGGGUU',
                  4:'AACCGGGU'}
        res_50 = {3:'AAGGGGUU',
                  4:'AACCGGGU'}
        res_80 = {4:'AACCGGGU'}
        res_100 = {}
        
        self.assertEqual(get_best_muscle_hits(subject_seq,query_aln,.2),res_20)
        self.assertEqual(get_best_muscle_hits(subject_seq,query_aln,.5),res_50)
        self.assertEqual(get_best_muscle_hits(subject_seq,query_aln,.8),res_80)
        self.assertEqual(get_best_muscle_hits(subject_seq,query_aln,1),res_100)
    
    def test_get_matching_chains(self):
        """Tests for get_matching_chains function.
        """
        if not MUSCLE_PATH:
            return 'skipping test'
        subject_seq = 'GCGACCCUG'
        res_30 = {'A': 'GCCACCCUG', 'B': 'CAGGGUCGGC'}
        res_80 = {'A': 'GCCACCCUG'}
        res_100 = {}
        #Threshold of .3
        test_30, ungapped_to_pdb = get_matching_chains(subject_seq, \
                                        TEST_PDB_STRING_1.split('\n'),\
                                        subject_type='Nucleotide',\
                                        threshold=.3)
        #Threshold of .8
        test_80, ungapped_to_pdb = get_matching_chains(subject_seq, \
                                        TEST_PDB_STRING_1.split('\n'),\
                                        subject_type='Nucleotide',\
                                        threshold=.8)
        #Threshold of 1
        test_100, ungapped_to_pdb = get_matching_chains(subject_seq, \
                                        TEST_PDB_STRING_1.split('\n'),\
                                        subject_type='Nucleotide',\
                                        threshold=1)
        #Incorrect subject_type
        #Threshold of .3
        test_wrong_subject, ungapped_to_pdb = get_matching_chains(subject_seq, \
                                        TEST_PDB_STRING_1.split('\n'),\
                                        subject_type='Protein',\
                                        threshold=.3)
        self.assertEqual(test_30,res_30)
        self.assertEqual(test_80,res_80)
        self.assertEqual(test_100,res_100)
        self.assertEqual(test_wrong_subject,{})
    
    def test_align_subject_to_pdb(self):
        """Tests for align_subject_to_pdb function.
        """
        if not MUSCLE_PATH:
            return 'skipping test'
        subject_seq = 'GCGACCCUG'
        pdb_matching = {'A': 'GCCACCCUG', 'B': 'CAGGGUCGGC'}
        result = {'A':('GCGACCCUG', 'GCCACCCUG'), \
                  'B':('GCGACCCUG-', 'CAGGGUCGGC')}
        self.assertEqual(align_subject_to_pdb(subject_seq,pdb_matching),result)
    
    def test_make_color_list(self):
        """Tests for make_color_list function.
        """
        colors = [(1.0,1.0,1.0),(1.0,0.0,1.0),(.5,.5,.5)]
        res = [('color_1',(1.0,1.0,1.0)),\
               ('color_2',(1.0,0.0,1.0)),\
               ('color_3',(.5,.5,.5))]
        
        self.assertEqual(make_color_list(colors),res)
    
    
    
TEST_PDB_STRING_1 = """
HEADER    RIBONUCLEIC ACID                        04-JAN-00   1DQH              
TITLE     CRYSTAL STRUCTURE OF HELIX II OF THE X. LAEVIS SOMATIC 5S             
TITLE    2 RRNA WITH A CYTOSINE BULGE IN TWO CONFORMATIONS                               
CRYST1   32.780   32.780  102.500  90.00  90.00  90.00 P 43 21 2     8                                 
ATOM      1  O5*   G A   1      38.612  13.536  39.204  1.00 37.41           O  
ATOM      2  C5*   G A   1      39.496  12.419  39.356  1.00 34.43           C  
ATOM      3  C4*   G A   1      38.750  11.165  39.729  1.00 33.33           C  
ATOM      4  O4*   G A   1      38.129  11.341  41.036  1.00 33.01           O  
ATOM      5  C3*   G A   1      37.614  10.780  38.799  1.00 33.57           C  
ATOM      6  O3*   G A   1      38.110   9.968  37.740  1.00 34.48           O  
ATOM      7  C2*   G A   1      36.722   9.952  39.715  1.00 32.60           C  
ATOM      8  O2*   G A   1      37.155   8.606  39.837  1.00 31.59           O  
ATOM      9  C1*   G A   1      36.871  10.693  41.048  1.00 31.18           C  
ATOM     10  N9    G A   1      35.863  11.709  41.312  1.00 29.41           N  
ATOM     11  C8    G A   1      36.087  13.049  41.490  1.00 28.46           C  
ATOM     12  N7    G A   1      34.996  13.725  41.725  1.00 28.77           N  
ATOM     13  C5    G A   1      33.986  12.775  41.711  1.00 28.42           C  
ATOM     14  C6    G A   1      32.586  12.918  41.929  1.00 28.30           C  
ATOM     15  O6    G A   1      31.954  13.943  42.190  1.00 27.86           O  
ATOM     16  N1    G A   1      31.921  11.699  41.818  1.00 27.81           N  
ATOM     17  C2    G A   1      32.533  10.496  41.534  1.00 27.95           C  
ATOM     18  N2    G A   1      31.731   9.428  41.426  1.00 27.06           N  
ATOM     19  N3    G A   1      33.844  10.357  41.352  1.00 28.40           N  
ATOM     20  C4    G A   1      34.498  11.525  41.450  1.00 28.58           C  
ATOM     21  P     C A   2      37.584  10.166  36.231  1.00 35.40           P  
ATOM     22  O1P   C A   2      38.502   9.371  35.370  1.00 37.17           O  
ATOM     23  O2P   C A   2      37.371  11.599  35.951  1.00 34.58           O  
ATOM     24  O5*   C A   2      36.188   9.406  36.193  1.00 33.54           O  
ATOM     25  C5*   C A   2      36.118   7.997  36.285  1.00 32.07           C  
ATOM     26  C4*   C A   2      34.683   7.567  36.457  1.00 30.69           C  
ATOM     27  O4*   C A   2      34.156   8.055  37.727  1.00 29.94           O  
ATOM     28  C3*   C A   2      33.728   8.137  35.431  1.00 30.66           C  
ATOM     29  O3*   C A   2      33.789   7.366  34.249  1.00 31.59           O  
ATOM     30  C2*   C A   2      32.391   7.983  36.148  1.00 28.79           C  
ATOM     31  O2*   C A   2      31.937   6.639  36.197  1.00 29.01           O  
ATOM     32  C1*   C A   2      32.779   8.394  37.563  1.00 27.55           C  
ATOM     33  N1    C A   2      32.625   9.836  37.823  1.00 25.64           N  
ATOM     34  C2    C A   2      31.353  10.331  38.113  1.00 24.78           C  
ATOM     35  O2    C A   2      30.382   9.563  38.052  1.00 23.25           O  
ATOM     36  N3    C A   2      31.210  11.641  38.434  1.00 24.08           N  
ATOM     37  C4    C A   2      32.265  12.446  38.446  1.00 23.49           C  
ATOM     38  N4    C A   2      32.070  13.718  38.781  1.00 22.55           N  
ATOM     39  C5    C A   2      33.573  11.981  38.110  1.00 24.21           C  
ATOM     40  C6    C A   2      33.702  10.676  37.806  1.00 24.76           C  
ATOM     41  P     C A   3      33.724   8.080  32.829  1.00 33.04           P  
ATOM     42  O1P   C A   3      34.035   6.997  31.842  1.00 33.90           O  
ATOM     43  O2P   C A   3      34.513   9.331  32.832  1.00 32.17           O  
ATOM     44  O5*   C A   3      32.184   8.423  32.636  1.00 31.52           O  
ATOM     45  C5*   C A   3      31.244   7.380  32.516  1.00 29.58           C  
ATOM     46  C4*   C A   3      29.856   7.885  32.841  1.00 27.91           C  
ATOM     47  O4*   C A   3      29.809   8.390  34.201  1.00 26.77           O  
ATOM     48  C3*   C A   3      29.352   9.042  32.005  1.00 28.23           C  
ATOM     49  O3*   C A   3      28.917   8.548  30.744  1.00 30.13           O  
ATOM     50  C2*   C A   3      28.212   9.549  32.887  1.00 26.35           C  
ATOM     51  O2*   C A   3      27.079   8.690  32.854  1.00 23.86           O  
ATOM     52  C1*   C A   3      28.855   9.454  34.273  1.00 25.79           C  
ATOM     53  N1    C A   3      29.575  10.688  34.648  1.00 24.80           N  
ATOM     54  C2    C A   3      28.832  11.781  35.132  1.00 24.14           C  
ATOM     55  O2    C A   3      27.577  11.677  35.209  1.00 23.22           O  
ATOM     56  N3    C A   3      29.482  12.906  35.489  1.00 23.57           N  
ATOM     57  C4    C A   3      30.820  12.983  35.357  1.00 23.98           C  
ATOM     58  N4    C A   3      31.421  14.109  35.709  1.00 23.48           N  
ATOM     59  C5    C A   3      31.596  11.888  34.855  1.00 24.07           C  
ATOM     60  C6    C A   3      30.938  10.773  34.525  1.00 24.08           C  
ATOM     61  P     A A   4      28.825   9.524  29.474  1.00 32.73           P  
ATOM     62  O1P   A A   4      28.402   8.666  28.346  1.00 34.23           O  
ATOM     63  O2P   A A   4      30.032  10.358  29.329  1.00 32.79           O  
ATOM     64  O5*   A A   4      27.664  10.569  29.796  1.00 32.19           O  
ATOM     65  C5*   A A   4      26.305  10.158  29.955  1.00 30.34           C  
ATOM     66  C4*   A A   4      25.469  11.306  30.483  1.00 29.48           C  
ATOM     67  O4*   A A   4      25.962  11.694  31.789  1.00 27.46           O  
ATOM     68  C3*   A A   4      25.495  12.610  29.694  1.00 29.92           C  
ATOM     69  O3*   A A   4      24.603  12.575  28.598  1.00 32.61           O  
ATOM     70  C2*   A A   4      25.037  13.609  30.748  1.00 28.81           C  
ATOM     71  O2*   A A   4      23.640  13.574  30.987  1.00 29.76           O  
ATOM     72  C1*   A A   4      25.761  13.085  31.989  1.00 26.51           C  
ATOM     73  N9    A A   4      27.073  13.686  32.195  1.00 25.35           N  
ATOM     74  C8    A A   4      28.277  13.137  31.857  1.00 24.29           C  
ATOM     75  N7    A A   4      29.306  13.870  32.184  1.00 23.88           N  
ATOM     76  C5    A A   4      28.741  14.980  32.793  1.00 24.49           C  
ATOM     77  C6    A A   4      29.305  16.106  33.385  1.00 23.89           C  
ATOM     78  N6    A A   4      30.621  16.294  33.487  1.00 23.21           N  
ATOM     79  N1    A A   4      28.468  17.036  33.883  1.00 23.83           N  
ATOM     80  C2    A A   4      27.145  16.832  33.779  1.00 24.11           C  
ATOM     81  N3    A A   4      26.497  15.805  33.254  1.00 23.77           N  
ATOM     82  C4    A A   4      27.365  14.897  32.777  1.00 24.53           C  
ATOM     83  P     C A   5      24.910  13.447  27.286  1.00 35.53           P  
ATOM     84  O1P   C A   5      23.998  12.939  26.228  1.00 37.30           O  
ATOM     85  O2P   C A   5      26.366  13.448  27.050  1.00 36.57           O  
ATOM     86  O5*   C A   5      24.512  14.941  27.678  1.00 34.43           O  
ATOM     87  C5*   C A   5      23.191  15.257  28.063  1.00 33.62           C  
ATOM     88  C4*   C A   5      23.144  16.642  28.633  1.00 33.30           C  
ATOM     89  O4*   C A   5      23.865  16.674  29.889  1.00 32.25           O  
ATOM     90  C3*   C A   5      23.834  17.712  27.809  1.00 33.29           C  
ATOM     91  O3*   C A   5      23.016  18.137  26.731  1.00 35.25           O  
ATOM     92  C2*   C A   5      24.027  18.781  28.860  1.00 32.42           C  
ATOM     93  O2*   C A   5      22.811  19.399  29.229  1.00 32.98           O  
ATOM     94  C1*   C A   5      24.494  17.930  30.040  1.00 31.25           C  
ATOM     95  N1    C A   5      25.956  17.721  30.055  1.00 29.62           N  
ATOM     96  C2    C A   5      26.739  18.692  30.649  1.00 28.72           C  
ATOM     97  O2    C A   5      26.182  19.698  31.113  1.00 29.29           O  
ATOM     98  N3    C A   5      28.074  18.527  30.707  1.00 28.00           N  
ATOM     99  C4    C A   5      28.632  17.446  30.193  1.00 27.58           C  
ATOM    100  N4    C A   5      29.955  17.322  30.311  1.00 27.48           N  
ATOM    101  C5    C A   5      27.862  16.439  29.546  1.00 27.81           C  
ATOM    102  C6    C A   5      26.532  16.614  29.508  1.00 28.20           C  
ATOM    103  P     C A   6      23.692  18.693  25.380  1.00 36.34           P  
ATOM    104  O1P   C A   6      22.530  18.995  24.495  1.00 37.64           O  
ATOM    105  O2P   C A   6      24.762  17.807  24.897  1.00 35.91           O  
ATOM    106  O5*   C A   6      24.353  20.071  25.823  1.00 34.59           O  
ATOM    107  C5*   C A   6      23.528  21.138  26.251  1.00 33.84           C  
ATOM    108  C4*   C A   6      24.359  22.254  26.818  1.00 33.26           C  
ATOM    109  O4*   C A   6      25.110  21.734  27.937  1.00 32.57           O  
ATOM    110  C3*   C A   6      25.439  22.842  25.933  1.00 33.36           C  
ATOM    111  O3*   C A   6      24.918  23.787  25.011  1.00 34.78           O  
ATOM    112  C2*   C A   6      26.314  23.522  26.974  1.00 32.29           C  
ATOM    113  O2*   C A   6      25.695  24.687  27.482  1.00 33.47           O  
ATOM    114  C1*   C A   6      26.331  22.447  28.060  1.00 31.49           C  
ATOM    115  N1    C A   6      27.434  21.495  27.881  1.00 29.26           N  
ATOM    116  C2    C A   6      28.720  21.876  28.301  1.00 29.08           C  
ATOM    117  O2    C A   6      28.868  22.999  28.797  1.00 29.07           O  
ATOM    118  N3    C A   6      29.752  21.012  28.149  1.00 28.04           N  
ATOM    119  C4    C A   6      29.534  19.814  27.609  1.00 28.05           C  
ATOM    120  N4    C A   6      30.564  18.975  27.495  1.00 28.92           N  
ATOM    121  C5    C A   6      28.230  19.407  27.162  1.00 28.53           C  
ATOM    122  C6    C A   6      27.226  20.270  27.318  1.00 28.47           C  
ATOM    123  P     C A   7      25.633  23.980  23.581  1.00 34.67           P  
ATOM    124  O1P   C A   7      24.754  24.860  22.776  1.00 34.93           O  
ATOM    125  O2P   C A   7      26.052  22.667  23.033  1.00 33.74           O  
ATOM    126  O5*   C A   7      26.947  24.804  23.919  1.00 33.62           O  
ATOM    127  C5*   C A   7      26.834  26.089  24.487  1.00 32.39           C  
ATOM    128  C4*   C A   7      28.195  26.619  24.841  1.00 32.43           C  
ATOM    129  O4*   C A   7      28.796  25.780  25.861  1.00 31.40           O  
ATOM    130  C3*   C A   7      29.216  26.582  23.727  1.00 32.59           C  
ATOM    131  O3*   C A   7      29.039  27.661  22.832  1.00 33.34           O  
ATOM    132  C2*   C A   7      30.510  26.706  24.503  1.00 30.81           C  
ATOM    133  O2*   C A   7      30.730  28.021  24.964  1.00 31.51           O  
ATOM    134  C1*   C A   7      30.207  25.796  25.698  1.00 30.43           C  
ATOM    135  N1    C A   7      30.670  24.427  25.465  1.00 28.25           N  
ATOM    136  C2    C A   7      32.021  24.155  25.622  1.00 27.81           C  
ATOM    137  O2    C A   7      32.798  25.101  25.915  1.00 26.72           O  
ATOM    138  N3    C A   7      32.458  22.883  25.457  1.00 26.97           N  
ATOM    139  C4    C A   7      31.591  21.912  25.152  1.00 26.73           C  
ATOM    140  N4    C A   7      32.056  20.670  25.049  1.00 25.90           N  
ATOM    141  C5    C A   7      30.206  22.176  24.955  1.00 26.90           C  
ATOM    142  C6    C A   7      29.796  23.434  25.113  1.00 28.38           C  
ATOM    143  P     U A   8      29.356  27.448  21.284  1.00 34.70           P  
ATOM    144  O1P   U A   8      28.951  28.708  20.571  1.00 35.62           O  
ATOM    145  O2P   U A   8      28.821  26.145  20.836  1.00 34.46           O  
ATOM    146  O5*   U A   8      30.942  27.354  21.220  1.00 32.95           O  
ATOM    147  C5*   U A   8      31.739  28.423  21.710  1.00 32.04           C  
ATOM    148  C4*   U A   8      33.182  28.009  21.784  1.00 32.28           C  
ATOM    149  O4*   U A   8      33.347  26.969  22.780  1.00 31.32           O  
ATOM    150  C3*   U A   8      33.778  27.374  20.544  1.00 32.35           C  
ATOM    151  O3*   U A   8      34.101  28.349  19.563  1.00 34.03           O  
ATOM    152  C2*   U A   8      35.026  26.728  21.124  1.00 31.88           C  
ATOM    153  O2*   U A   8      36.029  27.682  21.404  1.00 31.73           O  
ATOM    154  C1*   U A   8      34.481  26.181  22.445  1.00 31.32           C  
ATOM    155  N1    U A   8      34.022  24.800  22.296  1.00 30.13           N  
ATOM    156  C2    U A   8      34.954  23.806  22.406  1.00 29.56           C  
ATOM    157  O2    U A   8      36.151  24.027  22.587  1.00 30.22           O  
ATOM    158  N3    U A   8      34.454  22.540  22.298  1.00 28.90           N  
ATOM    159  C4    U A   8      33.147  22.177  22.079  1.00 28.56           C  
ATOM    160  O4    U A   8      32.873  20.991  21.974  1.00 28.38           O  
ATOM    161  C5    U A   8      32.230  23.269  21.949  1.00 29.20           C  
ATOM    162  C6    U A   8      32.691  24.519  22.059  1.00 29.91           C  
ATOM    163  P     G A   9      34.056  27.954  18.027  1.00 35.48           P  
ATOM    164  O1P   G A   9      34.253  29.197  17.239  1.00 36.23           O  
ATOM    165  O2P   G A   9      32.850  27.108  17.775  1.00 36.21           O  
ATOM    166  O5*   G A   9      35.348  27.048  17.835  1.00 32.95           O  
ATOM    167  C5*   G A   9      36.637  27.612  18.017  1.00 31.54           C  
ATOM    168  C4*   G A   9      37.696  26.545  17.935  1.00 29.88           C  
ATOM    169  O4*   G A   9      37.533  25.628  19.041  1.00 29.65           O  
ATOM    170  C3*   G A   9      37.696  25.633  16.719  1.00 29.94           C  
ATOM    171  O3*   G A   9      38.321  26.196  15.566  1.00 30.72           O  
ATOM    172  C2*   G A   9      38.513  24.452  17.217  1.00 29.09           C  
ATOM    173  O2*   G A   9      39.906  24.698  17.166  1.00 30.51           O  
ATOM    174  C1*   G A   9      38.035  24.354  18.672  1.00 28.52           C  
ATOM    175  N9    G A   9      36.951  23.383  18.827  1.00 27.70           N  
ATOM    176  C8    G A   9      35.600  23.640  18.864  1.00 27.66           C  
ATOM    177  N7    G A   9      34.873  22.551  18.981  1.00 27.66           N  
ATOM    178  C5    G A   9      35.807  21.518  19.032  1.00 27.08           C  
ATOM    179  C6    G A   9      35.622  20.108  19.171  1.00 26.99           C  
ATOM    180  O6    G A   9      34.557  19.479  19.295  1.00 27.32           O  
ATOM    181  N1    G A   9      36.834  19.430  19.166  1.00 25.88           N  
ATOM    182  C2    G A   9      38.067  20.024  19.072  1.00 26.62           C  
ATOM    183  N2    G A   9      39.121  19.203  19.094  1.00 26.35           N  
ATOM    184  N3    G A   9      38.253  21.334  18.960  1.00 26.19           N  
ATOM    185  C4    G A   9      37.092  22.013  18.942  1.00 26.94           C  
TER     186        G A   9                                                      
ATOM    187  O5*   C B  10      37.876  10.866  21.876  1.00 38.18           O  
ATOM    188  C5*   C B  10      39.087  10.527  21.197  1.00 34.41           C  
ATOM    189  C4*   C B  10      39.746  11.780  20.669  1.00 34.90           C  
ATOM    190  O4*   C B  10      38.931  12.392  19.627  1.00 33.24           O  
ATOM    191  C3*   C B  10      39.904  12.927  21.657  1.00 34.06           C  
ATOM    192  O3*   C B  10      40.989  12.675  22.550  1.00 35.88           O  
ATOM    193  C2*   C B  10      40.214  14.067  20.695  1.00 32.72           C  
ATOM    194  O2*   C B  10      41.499  13.919  20.131  1.00 33.32           O  
ATOM    195  C1*   C B  10      39.202  13.792  19.582  1.00 31.82           C  
ATOM    196  N1    C B  10      37.945  14.536  19.750  1.00 30.40           N  
ATOM    197  C2    C B  10      37.944  15.911  19.492  1.00 29.71           C  
ATOM    198  O2    C B  10      39.014  16.463  19.179  1.00 29.35           O  
ATOM    199  N3    C B  10      36.795  16.602  19.605  1.00 29.04           N  
ATOM    200  C4    C B  10      35.677  15.980  19.976  1.00 29.04           C  
ATOM    201  N4    C B  10      34.555  16.703  20.064  1.00 28.12           N  
ATOM    202  C5    C B  10      35.656  14.584  20.270  1.00 29.17           C  
ATOM    203  C6    C B  10      36.802  13.911  20.146  1.00 29.79           C  
ATOM    204  P     A B  11      40.901  13.173  24.071  1.00 38.04           P  
ATOM    205  O1P   A B  11      42.040  12.552  24.800  1.00 37.87           O  
ATOM    206  O2P   A B  11      39.521  13.016  24.582  1.00 36.68           O  
ATOM    207  O5*   A B  11      41.192  14.729  23.921  1.00 35.30           O  
ATOM    208  C5*   A B  11      42.473  15.172  23.503  1.00 33.64           C  
ATOM    209  C4*   A B  11      42.479  16.676  23.341  1.00 31.60           C  
ATOM    210  O4*   A B  11      41.587  17.051  22.255  1.00 30.20           O  
ATOM    211  C3*   A B  11      41.947  17.462  24.525  1.00 30.78           C  
ATOM    212  O3*   A B  11      42.937  17.594  25.538  1.00 30.99           O  
ATOM    213  C2*   A B  11      41.593  18.790  23.867  1.00 29.71           C  
ATOM    214  O2*   A B  11      42.736  19.561  23.571  1.00 29.00           O  
ATOM    215  C1*   A B  11      40.990  18.307  22.547  1.00 29.66           C  
ATOM    216  N9    A B  11      39.547  18.127  22.638  1.00 28.32           N  
ATOM    217  C8    A B  11      38.844  16.967  22.817  1.00 28.32           C  
ATOM    218  N7    A B  11      37.548  17.134  22.849  1.00 28.39           N  
ATOM    219  C5    A B  11      37.382  18.501  22.678  1.00 27.65           C  
ATOM    220  C6    A B  11      36.236  19.324  22.601  1.00 27.42           C  
ATOM    221  N6    A B  11      34.973  18.877  22.692  1.00 27.87           N  
ATOM    222  N1    A B  11      36.433  20.651  22.423  1.00 26.83           N  
ATOM    223  C2    A B  11      37.687  21.105  22.330  1.00 27.23           C  
ATOM    224  N3    A B  11      38.837  20.433  22.378  1.00 27.13           N  
ATOM    225  C4    A B  11      38.610  19.123  22.553  1.00 27.90           C  
ATOM    226  P     G B  12      42.493  17.715  27.072  1.00 32.74           P  
ATOM    227  O1P   G B  12      43.728  17.509  27.884  1.00 34.58           O  
ATOM    228  O2P   G B  12      41.277  16.922  27.381  1.00 32.04           O  
ATOM    229  O5*   G B  12      42.027  19.234  27.207  1.00 30.99           O  
ATOM    230  C5*   G B  12      42.962  20.289  27.056  1.00 30.37           C  
ATOM    231  C4*   G B  12      42.254  21.614  27.115  1.00 29.14           C  
ATOM    232  O4*   G B  12      41.457  21.794  25.913  1.00 28.24           O  
ATOM    233  C3*   G B  12      41.237  21.787  28.229  1.00 28.75           C  
ATOM    234  O3*   G B  12      41.830  22.073  29.500  1.00 30.19           O  
ATOM    235  C2*   G B  12      40.436  22.958  27.690  1.00 27.60           C  
ATOM    236  O2*   G B  12      41.166  24.177  27.802  1.00 27.01           O  
ATOM    237  C1*   G B  12      40.309  22.557  26.218  1.00 27.19           C  
ATOM    238  N9    G B  12      39.139  21.706  26.021  1.00 26.10           N  
ATOM    239  C8    G B  12      39.095  20.333  25.971  1.00 25.61           C  
ATOM    240  N7    G B  12      37.878  19.869  25.854  1.00 25.39           N  
ATOM    241  C5    G B  12      37.080  21.003  25.805  1.00 24.51           C  
ATOM    242  C6    G B  12      35.677  21.130  25.709  1.00 25.12           C  
ATOM    243  O6    G B  12      34.819  20.224  25.646  1.00 24.39           O  
ATOM    244  N1    G B  12      35.282  22.463  25.710  1.00 23.14           N  
ATOM    245  C2    G B  12      36.136  23.544  25.801  1.00 25.12           C  
ATOM    246  N2    G B  12      35.566  24.781  25.809  1.00 23.32           N  
ATOM    247  N3    G B  12      37.447  23.430  25.888  1.00 23.89           N  
ATOM    248  C4    G B  12      37.847  22.143  25.892  1.00 24.92           C  
ATOM    249  P     G B  13      41.050  21.646  30.847  1.00 31.79           P  
ATOM    250  O1P   G B  13      41.913  21.970  32.013  1.00 33.13           O  
ATOM    251  O2P   G B  13      40.502  20.266  30.697  1.00 30.96           O  
ATOM    252  O5*   G B  13      39.754  22.570  30.921  1.00 30.06           O  
ATOM    253  C5*   G B  13      39.862  23.958  31.171  1.00 29.46           C  
ATOM    254  C4*   G B  13      38.543  24.649  30.914  1.00 28.17           C  
ATOM    255  O4*   G B  13      38.049  24.334  29.586  1.00 27.88           O  
ATOM    256  C3*   G B  13      37.356  24.322  31.810  1.00 27.17           C  
ATOM    257  O3*   G B  13      37.506  24.978  33.065  1.00 27.57           O  
ATOM    258  C2*   G B  13      36.240  24.941  30.983  1.00 27.19           C  
ATOM    259  O2*   G B  13      36.251  26.358  31.018  1.00 26.53           O  
ATOM    260  C1*   G B  13      36.629  24.474  29.573  1.00 26.87           C  
ATOM    261  N9    G B  13      36.031  23.162  29.330  1.00 26.99           N  
ATOM    262  C8    G B  13      36.660  21.949  29.286  1.00 25.89           C  
ATOM    263  N7    G B  13      35.829  20.951  29.099  1.00 26.13           N  
ATOM    264  C5    G B  13      34.583  21.549  28.997  1.00 25.76           C  
ATOM    265  C6    G B  13      33.304  20.970  28.786  1.00 25.44           C  
ATOM    266  O6    G B  13      33.020  19.783  28.643  1.00 26.51           O  
ATOM    267  N1    G B  13      32.308  21.927  28.756  1.00 26.07           N  
ATOM    268  C2    G B  13      32.506  23.276  28.917  1.00 26.27           C  
ATOM    269  N2    G B  13      31.394  24.054  28.887  1.00 26.49           N  
ATOM    270  N3    G B  13      33.702  23.829  29.104  1.00 25.76           N  
ATOM    271  C4    G B  13      34.683  22.911  29.131  1.00 25.78           C  
ATOM    272  P     G B  14      36.688  24.480  34.364  1.00 28.04           P  
ATOM    273  O1P   G B  14      37.232  25.265  35.498  1.00 28.24           O  
ATOM    274  O2P   G B  14      36.700  23.014  34.418  1.00 27.78           O  
ATOM    275  O5*   G B  14      35.189  24.948  34.089  1.00 26.64           O  
ATOM    276  C5*   G B  14      34.821  26.309  34.234  1.00 26.48           C  
ATOM    277  C4*   G B  14      33.329  26.470  34.021  1.00 26.73           C  
ATOM    278  O4*   G B  14      32.973  26.023  32.691  1.00 27.34           O  
ATOM    279  C3*   G B  14      32.455  25.631  34.927  1.00 26.32           C  
ATOM    280  O3*   G B  14      32.289  26.283  36.163  1.00 25.18           O  
ATOM    281  C2*   G B  14      31.150  25.594  34.155  1.00 26.28           C  
ATOM    282  O2*   G B  14      30.462  26.832  34.292  1.00 28.24           O  
ATOM    283  C1*   G B  14      31.685  25.398  32.734  1.00 27.44           C  
ATOM    284  N9    G B  14      31.902  23.979  32.496  1.00 26.19           N  
ATOM    285  C8    G B  14      33.102  23.325  32.481  1.00 25.54           C  
ATOM    286  N7    G B  14      32.989  22.052  32.226  1.00 25.02           N  
ATOM    287  C5    G B  14      31.622  21.853  32.065  1.00 24.97           C  
ATOM    288  C6    G B  14      30.885  20.664  31.768  1.00 25.65           C  
ATOM    289  O6    G B  14      31.314  19.514  31.556  1.00 24.40           O  
ATOM    290  N1    G B  14      29.512  20.911  31.724  1.00 23.81           N  
ATOM    291  C2    G B  14      28.927  22.133  31.942  1.00 25.97           C  
ATOM    292  N2    G B  14      27.586  22.167  31.875  1.00 26.21           N  
ATOM    293  N3    G B  14      29.603  23.248  32.211  1.00 25.92           N  
ATOM    294  C4    G B  14      30.936  23.030  32.250  1.00 25.92           C  
ATOM    295  P     U B  15      32.149  25.416  37.480  1.00 25.35           P  
ATOM    296  O1P   U B  15      32.057  26.334  38.616  1.00 27.54           O  
ATOM    297  O2P   U B  15      33.208  24.374  37.459  1.00 26.82           O  
ATOM    298  O5*   U B  15      30.740  24.685  37.321  1.00 24.63           O  
ATOM    299  C5*   U B  15      29.549  25.454  37.213  1.00 25.46           C  
ATOM    300  C4*   U B  15      28.349  24.546  37.085  1.00 24.23           C  
ATOM    301  O4*   U B  15      28.377  23.938  35.774  1.00 23.55           O  
ATOM    302  C3*   U B  15      28.321  23.355  38.026  1.00 24.30           C  
ATOM    303  O3*   U B  15      27.841  23.737  39.309  1.00 25.05           O  
ATOM    304  C2*   U B  15      27.338  22.453  37.289  1.00 23.86           C  
ATOM    305  O2*   U B  15      25.979  22.899  37.359  1.00 23.23           O  
ATOM    306  C1*   U B  15      27.826  22.628  35.851  1.00 24.23           C  
ATOM    307  N1    U B  15      28.871  21.650  35.545  1.00 23.38           N  
ATOM    308  C2    U B  15      28.448  20.398  35.155  1.00 24.01           C  
ATOM    309  O2    U B  15      27.270  20.112  35.029  1.00 22.67           O  
ATOM    310  N3    U B  15      29.458  19.491  34.928  1.00 24.00           N  
ATOM    311  C4    U B  15      30.813  19.710  35.044  1.00 24.09           C  
ATOM    312  O4    U B  15      31.590  18.793  34.769  1.00 24.56           O  
ATOM    313  C5    U B  15      31.178  21.033  35.445  1.00 24.19           C  
ATOM    314  C6    U B  15      30.217  21.944  35.673  1.00 23.69           C  
ATOM    315  P     C B  16      28.396  23.004  40.622  1.00 24.91           P  
ATOM    316  O1P   C B  16      29.881  23.085  40.622  1.00 24.84           O  
ATOM    317  O2P   C B  16      27.730  21.661  40.668  1.00 24.89           O  
ATOM    318  O5*   C B  16      27.946  23.943  41.822  1.00 25.30           O  
ATOM    319  C5*   C B  16      26.884  23.591  42.712  1.00 25.79           C  
ATOM    320  C4*   C B  16      25.915  24.750  42.836  1.00 25.83           C  
ATOM    321  O4*   C B  16      26.638  25.986  43.088  1.00 25.96           O  
ATOM    322  C3*   C B  16      25.135  25.023  41.571  1.00 25.48           C  
ATOM    323  O3*   C B  16      23.978  24.188  41.595  1.00 24.15           O  
ATOM    324  C2*   C B  16      24.736  26.489  41.738  1.00 25.72           C  
ATOM    325  O2*   C B  16      23.575  26.576  42.539  1.00 25.74           O  
ATOM    326  C1*   C B  16      25.947  27.061  42.477  1.00 25.87           C  
ATOM    327  N1    C B  16      26.880  27.811  41.611  1.00 26.24           N  
ATOM    328  C2    C B  16      26.585  29.152  41.350  1.00 26.56           C  
ATOM    329  O2    C B  16      25.575  29.642  41.883  1.00 25.85           O  
ATOM    330  N3    C B  16      27.407  29.874  40.549  1.00 25.78           N  
ATOM    331  C4    C B  16      28.502  29.306  40.033  1.00 27.00           C  
ATOM    332  N4    C B  16      29.310  30.069  39.261  1.00 26.12           N  
ATOM    333  C5    C B  16      28.832  27.943  40.289  1.00 26.61           C  
ATOM    334  C6    C B  16      27.998  27.237  41.079  1.00 26.28           C  
ATOM    335  P     G B  17      23.384  23.642  40.245  1.00 23.95           P  
ATOM    336  O1P   G B  17      23.436  24.574  39.098  1.00 25.14           O  
ATOM    337  O2P   G B  17      22.082  23.017  40.615  1.00 27.22           O  
ATOM    338  O5*   G B  17      24.421  22.475  39.840  1.00 26.25           O  
ATOM    339  C5*   G B  17      24.386  21.204  40.495  1.00 26.12           C  
ATOM    340  C4*   G B  17      23.742  20.170  39.592  1.00 24.76           C  
ATOM    341  O4*   G B  17      24.528  20.039  38.380  1.00 24.76           O  
ATOM    342  C3*   G B  17      23.732  18.759  40.174  1.00 25.57           C  
ATOM    343  O3*   G B  17      22.581  18.586  40.985  1.00 25.04           O  
ATOM    344  C2*   G B  17      23.674  17.891  38.931  1.00 23.99           C  
ATOM    345  O2*   G B  17      22.367  17.810  38.389  1.00 26.46           O  
ATOM    346  C1*   G B  17      24.591  18.671  37.983  1.00 24.16           C  
ATOM    347  N9    G B  17      26.001  18.255  38.013  1.00 22.61           N  
ATOM    348  C8    G B  17      27.064  18.970  38.508  1.00 22.62           C  
ATOM    349  N7    G B  17      28.216  18.388  38.321  1.00 22.00           N  
ATOM    350  C5    G B  17      27.898  17.197  37.680  1.00 22.60           C  
ATOM    351  C6    G B  17      28.740  16.163  37.189  1.00 22.09           C  
ATOM    352  O6    G B  17      29.987  16.125  37.154  1.00 22.31           O  
ATOM    353  N1    G B  17      28.003  15.117  36.663  1.00 20.99           N  
ATOM    354  C2    G B  17      26.634  15.088  36.560  1.00 21.61           C  
ATOM    355  N2    G B  17      26.130  13.961  36.055  1.00 21.54           N  
ATOM    356  N3    G B  17      25.834  16.080  36.938  1.00 20.35           N  
ATOM    357  C4    G B  17      26.525  17.085  37.502  1.00 22.36           C  
ATOM    358  P     G B  18      22.687  17.714  42.302  1.00 26.50           P  
ATOM    359  O1P   G B  18      21.410  17.993  43.015  1.00 27.27           O  
ATOM    360  O2P   G B  18      23.973  17.912  42.983  1.00 23.81           O  
ATOM    361  O5*   G B  18      22.672  16.216  41.755  1.00 23.89           O  
ATOM    362  C5*   G B  18      21.485  15.647  41.200  1.00 25.00           C  
ATOM    363  C4*   G B  18      21.782  14.282  40.588  1.00 24.30           C  
ATOM    364  O4*   G B  18      22.684  14.431  39.468  1.00 22.31           O  
ATOM    365  C3*   G B  18      22.476  13.275  41.487  1.00 24.94           C  
ATOM    366  O3*   G B  18      21.499  12.593  42.250  1.00 26.00           O  
ATOM    367  C2*   G B  18      23.118  12.339  40.468  1.00 22.26           C  
ATOM    368  O2*   G B  18      22.152  11.549  39.777  1.00 22.69           O  
ATOM    369  C1*   G B  18      23.616  13.361  39.455  1.00 22.64           C  
ATOM    370  N9    G B  18      24.946  13.869  39.792  1.00 20.90           N  
ATOM    371  C8    G B  18      25.300  15.055  40.391  1.00 21.43           C  
ATOM    372  N7    G B  18      26.600  15.220  40.446  1.00 20.70           N  
ATOM    373  C5    G B  18      27.116  14.070  39.870  1.00 20.39           C  
ATOM    374  C6    G B  18      28.475  13.671  39.631  1.00 21.70           C  
ATOM    375  O6    G B  18      29.522  14.295  39.877  1.00 21.10           O  
ATOM    376  N1    G B  18      28.545  12.408  39.051  1.00 22.16           N  
ATOM    377  C2    G B  18      27.475  11.625  38.743  1.00 21.99           C  
ATOM    378  N2    G B  18      27.770  10.418  38.218  1.00 23.07           N  
ATOM    379  N3    G B  18      26.205  11.986  38.939  1.00 21.12           N  
ATOM    380  C4    G B  18      26.110  13.214  39.499  1.00 21.28           C  
ATOM    381  P     C B  19      21.883  12.026  43.669  1.00 27.44           P  
ATOM    382  O1P   C B  19      20.671  11.313  44.167  1.00 28.17           O  
ATOM    383  O2P   C B  19      22.528  13.044  44.526  1.00 26.74           O  
ATOM    384  O5*   C B  19      23.028  10.969  43.369  1.00 27.50           O  
ATOM    385  C5*   C B  19      22.755   9.790  42.650  1.00 28.82           C  
ATOM    386  C4*   C B  19      24.029   8.991  42.502  1.00 29.58           C  
ATOM    387  O4*   C B  19      24.980   9.694  41.662  1.00 29.75           O  
ATOM    388  C3*   C B  19      24.790   8.765  43.794  1.00 29.78           C  
ATOM    389  O3*   C B  19      24.156   7.716  44.580  1.00 32.33           O  
ATOM    390  C2*   C B  19      26.182   8.418  43.280  1.00 29.77           C  
ATOM    391  O2*   C B  19      26.195   7.085  42.797  1.00 31.42           O  
ATOM    392  C1*   C B  19      26.302   9.364  42.073  1.00 28.44           C  
ATOM    393  N1    C B  19      27.024  10.602  42.393  1.00 27.04           N  
ATOM    394  C2    C B  19      28.387  10.641  42.164  1.00 26.22           C  
ATOM    395  O2    C B  19      28.935   9.617  41.731  1.00 25.20           O  
ATOM    396  N3    C B  19      29.078  11.788  42.418  1.00 25.55           N  
ATOM    397  C4    C B  19      28.436  12.863  42.885  1.00 25.98           C  
ATOM    398  N4    C B  19      29.141  14.002  43.089  1.00 25.03           N  
ATOM    399  C5    C B  19      27.039  12.836  43.162  1.00 26.00           C  
ATOM    400  C6    C B  19      26.375  11.697  42.903  1.00 26.97           C  
TER     401        C B  19                                                      
MASTER      238    0    0    0    0    0    0    6  456    2    0    2          
END                                                                             
"""

TEST_PDB_STRING_2 = """
HEADER    ANTIOXIDANT ENZYME                      06-NOV-00   1HD2              
TITLE     HUMAN PEROXIREDOXIN 5                                                 
ATOM      1  N   ALA A   1      -7.101  53.135  16.957  1.00 88.42           N  
ANISOU    1  N   ALA A   1    12714   7605  13277    523  -2633   3491       N  
ATOM      2  CA  ALA A   1      -8.014  52.075  17.450  1.00 63.39           C  
ANISOU    2  CA  ALA A   1     8225   7477   8383   2990   -789   1435       C  
ATOM      3  C   ALA A   1      -7.241  50.817  17.757  1.00 46.53           C  
ANISOU    3  C   ALA A   1     5793   6074   5811   2042    989     16       C  
ATOM      4  O   ALA A   1      -6.073  50.698  17.346  1.00 53.54           O  
ANISOU    4  O   ALA A   1     5678   8269   6398   1327   1320   1080       O  
ATOM      5  CB  ALA A   1      -9.119  51.791  16.443  1.00 77.54           C  
ANISOU    5  CB  ALA A   1     9616  10505   9342   2125  -2228   2932       C  
ATOM      6  N   PRO A   2      -7.796  49.873  18.488  1.00 34.26           N  
ANISOU    6  N   PRO A   2     4045   4756   4215   1116    109  -2030       N  
ATOM      7  CA  PRO A   2      -6.966  48.670  18.750  1.00 29.30           C  
ANISOU    7  CA  PRO A   2     3413   4041   3677    360   -116  -2345       C  
ATOM      8  C   PRO A   2      -6.707  47.922  17.451  1.00 23.66           C  
ANISOU    8  C   PRO A   2     1982   4034   2972    261   -662  -1762       C  
ATOM      9  O   PRO A   2      -7.549  47.657  16.601  1.00 23.73           O  
ANISOU    9  O   PRO A   2     1855   3960   3202    209   -836  -1320       O  
ATOM     10  CB  PRO A   2      -7.774  47.860  19.732  1.00 30.97           C  
ANISOU   10  CB  PRO A   2     3871   4866   3032    146   -165  -2287       C  
ATOM     11  CG  PRO A   2      -8.779  48.807  20.281  1.00 35.77           C  
ANISOU   11  CG  PRO A   2     3448   6110   4032    639     23  -2099       C  
ATOM     12  CD  PRO A   2      -9.080  49.777  19.155  1.00 36.34           C  
ANISOU   12  CD  PRO A   2     3188   6479   4142    982   -759  -2220       C  
ATOM     13  N   ILE A   3      -5.409  47.566  17.357  1.00 19.45           N  
ANISOU   13  N   ILE A   3     1760   3495   2134    -75   -484  -1027       N  
ATOM     14  CA  ILE A   3      -5.040  46.822  16.152  1.00 18.02           C  
ANISOU   14  CA  ILE A   3     1703   3171   1973     10   -766   -936       C  
ATOM     15  C   ILE A   3      -5.689  45.452  16.192  1.00 18.71           C  
ANISOU   15  C   ILE A   3     1903   3233   1974   -124   -474   -787       C  
ATOM     16  O   ILE A   3      -5.915  44.901  17.260  1.00 20.34           O  
ANISOU   16  O   ILE A   3     2262   3409   2057    315   -565   -507       O  
ATOM     17  CB  ILE A   3      -3.513  46.734  16.025  1.00 17.09           C  
ANISOU   17  CB  ILE A   3     1740   2849   1903    104   -657   -844       C  
ATOM     18  CG1 ILE A   3      -3.034  46.332  14.628  1.00 20.89           C  
ANISOU   18  CG1 ILE A   3     2079   3804   2056   -100   -461  -1173       C  
ATOM     19  CG2 ILE A   3      -2.939  45.866  17.110  1.00 18.41           C  
ANISOU   19  CG2 ILE A   3     1664   3016   2316   -146   -718   -494       C  
ATOM     20  CD1 ILE A   3      -1.553  46.546  14.371  1.00 22.47           C  
ANISOU   20  CD1 ILE A   3     2492   3346   2701   -814    225   -681       C  
ATOM     21  N   LYS A   4      -6.016  44.945  14.979  1.00 18.15           N  
ANISOU   21  N   LYS A   4     2109   2759   2028   -146   -425   -709       N  
ATOM     22  CA  LYS A   4      -6.669  43.672  14.871  1.00 19.21           C  
ANISOU   22  CA  LYS A   4     1734   2945   2619   -166   -240  -1004       C  
ATOM     23  C   LYS A   4      -6.136  42.905  13.666  1.00 18.35           C  
ANISOU   23  C   LYS A   4     1621   2920   2432   -469   -242   -951       C  
ATOM     24  O   LYS A   4      -5.523  43.516  12.787  1.00 18.36           O  
ANISOU   24  O   LYS A   4     1727   2866   2383   -290   -380   -717       O  
ATOM     25  CB  LYS A   4      -8.179  43.880  14.717  1.00 22.82           C  
ANISOU   25  CB  LYS A   4     1797   3115   3760     78   -321  -1132       C  
ATOM     26  CG  LYS A   4      -8.515  44.601  13.433  1.00 36.33           C  
ANISOU   26  CG  LYS A   4     3420   5560   4825    996  -1411   -178       C  
ATOM     27  CD  LYS A   4      -9.973  44.994  13.306  1.00 47.94           C  
ANISOU   27  CD  LYS A   4     3683   7877   6655   1611  -1898     19       C  
ATOM     28  CE  LYS A   4     -10.268  45.683  11.970  1.00 55.24           C  
ANISOU   28  CE  LYS A   4     4630   8886   7472   2301  -2486    712       C
ATOM     29  NZ  LYS A   4      -9.697  47.057  11.880  1.00 64.74           N  
ANISOU   29  NZ  LYS A   4     6567   9962   8070    776  -3006   2044       N  
ATOM     30  N   VAL A   5      -6.390  41.610  13.638  1.00 17.46           N  
ANISOU   30  N   VAL A   5     1965   2729   1939    -44   -514   -643       N  
ATOM     31  CA  VAL A   5      -6.062  40.803  12.447  1.00 17.65           C  
ANISOU   31  CA  VAL A   5     2357   2607   1741     11   -466   -478       C  
ATOM     32  C   VAL A   5      -6.706  41.431  11.221  1.00 16.38           C  
ANISOU   32  C   VAL A   5     1801   2477   1947     27   -578   -599       C  
ATOM     33  O   VAL A   5      -7.842  41.860  11.225  1.00 20.73           O  
ANISOU   33  O   VAL A   5     1855   3251   2769    339   -647   -999       O  
ATOM     34  CB  VAL A   5      -6.540  39.355  12.621  1.00 18.16           C  
ANISOU   34  CB  VAL A   5     2548   2687   1666   -199   -700   -506       C  
ATOM     35  CG1 VAL A   5      -6.490  38.556  11.331  1.00 20.68           C  
ANISOU   35  CG1 VAL A   5     3654   2669   1532   -286   -788   -390       C  
ATOM     36  CG2 VAL A   5      -5.643  38.711  13.693  1.00 21.32           C  
ANISOU   36  CG2 VAL A   5     3412   3031   1657   -320   -951   -136       C  
ATOM     37  N   GLY A   6      -5.884  41.470  10.169  1.00 17.57           N  
ANISOU   37  N   GLY A   6     1818   3021   1838     -5   -684   -127       N  
ATOM     38  CA  GLY A   6      -6.293  42.101   8.926  1.00 17.93           C  
ANISOU   38  CA  GLY A   6     2091   2771   1951    266  -1009   -287       C  
ATOM     39  C   GLY A   6      -5.787  43.509   8.756  1.00 19.23           C  
ANISOU   39  C   GLY A   6     2774   2651   1880    339   -804   -321       C  
ATOM     40  O   GLY A   6      -5.730  44.041   7.631  1.00 18.94           O  
ANISOU   40  O   GLY A   6     2552   2787   1858    397   -657   -317       O  
ATOM     41  N   ASP A   7      -5.412  44.192   9.821  1.00 16.89           N  
ANISOU   41  N   ASP A   7     2075   2493   1851    359   -716   -218       N  
ATOM     42  CA  ASP A   7      -4.884  45.527   9.687  1.00 17.92           C  
ANISOU   42  CA  ASP A   7     2197   2589   2023    267   -381   -346       C  
ATOM     43  C   ASP A   7      -3.441  45.489   9.193  1.00 15.85           C  
ANISOU   43  C   ASP A   7     2048   2223   1750    430   -734   -281       C  
ATOM     44  O   ASP A   7      -2.691  44.572   9.409  1.00 18.38           O  
ANISOU   44  O   ASP A   7     2451   2196   2337    656  -1042   -611       O  
ATOM     45  CB  ASP A   7      -4.884  46.242  11.037  1.00 19.75           C  
ANISOU   45  CB  ASP A   7     2401   2759   2345    469   -423   -702       C  
ATOM     46  CG  ASP A   7      -6.256  46.558  11.580  1.00 19.89           C  
ANISOU   46  CG  ASP A   7     2495   3124   1940    591   -278   -147       C  
ATOM     47  OD1 ASP A   7      -7.246  46.486  10.836  1.00 22.27           O  
ANISOU   47  OD1 ASP A   7     2413   3699   2348    528   -394   -258       O  
ATOM     48  OD2 ASP A   7      -6.286  46.895  12.814  1.00 23.63           O  
ANISOU   48  OD2 ASP A   7     3176   3849   1952    665   -161   -245       O  
ATOM     49  N   ALA A   8      -3.094  46.594   8.534  1.00 16.73           N  
ANISOU   49  N   ALA A   8     2285   2228   1845    255   -469   -376       N  
ATOM     50  CA  ALA A   8      -1.686  46.796   8.224  1.00 19.26           C  
ANISOU   50  CA  ALA A   8     2282   3100   1936    119   -555   -155       C  
ATOM     51  C   ALA A   8      -0.940  47.209   9.477  1.00 18.08           C  
ANISOU   51  C   ALA A   8     2249   2719   1900    181   -453   -242       C  
ATOM     52  O   ALA A   8      -1.418  47.960  10.308  1.00 20.61           O  
ANISOU   52  O   ALA A   8     2470   3061   2299    280   -308   -530       O  
ATOM     53  CB  ALA A   8      -1.558  47.881   7.175  1.00 22.45           C  
ANISOU   53  CB  ALA A   8     2906   3904   1718     16    -78     87       C  
MASTER      245    0    6    6    7    0    3    6 1429    1    9   13          
END                                                                             
"""

#run if called from command-line
if __name__ == "__main__":
    main()
