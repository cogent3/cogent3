#!/usr/bin/env python
"""Tests for the Gibbs parser
"""
from __future__ import division
from cogent.util.unit_test import TestCase, main
import string
import re
from cogent.motif.util import Motif,Module
from cogent.core.moltype import DNA,RNA,PROTEIN
from cogent.parse.record import DelimitedSplitter
from cogent.parse.record_finder import LabeledRecordFinder
from cogent.parse.gibbs import get_sequence_and_motif_blocks, get_sequence_map,\
    get_motif_blocks, get_motif_sequences, get_motif_p_value, guess_alphabet,\
    build_module_objects, module_ids_to_int, GibbsParser
from math import exp

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Production"

class GibbsTests(TestCase):
    """Tests for gibbs parser.
    """

    def setUp(self):
        """Setup function for gibbs tests.
        """
        self.gibbs_lines = GIBBS_FILE.split('\n')
        self.sequence_map = {'1':'1091044',\
                            '10':'135765',\
                            '11':'1388082',\
                            '12':'140543',\
                            '13':'14286173',\
                            '14':'14578634',\
                            '15':'14600438',\
                            '16':'15218394',\
                            '17':'15597673',\
                            '18':'15599256',\
                            '19':'15602312',\
                            '2':'11467494',\
                            '20':'15605725',\
                            '21':'15605963',\
                            '22':'15609375',\
                            '23':'15609658',\
                            '24':'15613511',\
                            '25':'15614085',\
                            '26':'15614140',\
                            '27':'15615431',\
                            '28':'15643152',\
                            '29':'15672286',\
                            '3':'11499727',\
                            '30':'15790738',\
                            '31':'15791337',\
                            '32':'15801846',\
                            '33':'15805225',\
                            '34':'15805374',\
                            '35':'15807234',\
                            '36':'15826629',\
                            '37':'15899007',\
                            '38':'15899339',\
                            '39':'15964668',\
                            '4':'1174686',\
                            '40':'15966937',\
                            '41':'15988313',\
                            '42':'16078864',\
                            '43':'16123427',\
                            '44':'16125919',\
                            '45':'16330420',\
                            '46':'1633495',\
                            '47':'16501671',\
                            '48':'1651717',\
                            '49':'16759994',\
                            '5':'12044976',\
                            '50':'16761507',\
                            '51':'16803644',\
                            '52':'16804867',\
                            '53':'17229033',\
                            '54':'17229859',\
                            '55':'1729944',\
                            '56':'17531233',\
                            '57':'17537401',\
                            '58':'17547503',\
                            '59':'18309723',\
                            '6':'13186328',\
                            '60':'18313548',\
                            '61':'18406743',\
                            '62':'19173077',\
                            '63':'19554157',\
                            '64':'19705357',\
                            '65':'19746502',\
                            '66':'20092028',\
                            '67':'20151112',\
                            '68':'21112072',\
                            '69':'21222859',\
                            '7':'13358154',\
                            '70':'21223405',\
                            '71':'21227878',\
                            '72':'21283385',\
                            '73':'21674812',\
                            '74':'23098307',\
                            '75':'2649838',\
                            '76':'267116',\
                            '77':'27375582',\
                            '78':'2822332',\
                            '79':'30021713',\
                            '8':'13541053',\
                            '80':'3261501',\
                            '81':'3318841',\
                            '82':'3323237',\
                            '83':'4155972',\
                            '84':'4200327',\
                            '85':'4433065',\
                            '86':'4704732',\
                            '87':'4996210',\
                            '88':'5326864',\
                            '89':'6322180',\
                            '9':'13541117',\
                            '90':'6323138',\
                            '91':'6687568',\
                            '92':'6850955',\
                            '93':'7109697',\
                            '94':'7290567',\
                            '95':'9955016',\
                            '96':'15677788',\
                            }
        self.motif_a_lines = """
10 columns
Num Motifs: 27
   2,  1      72 klstq ILAISVDSPFSH lqyll      83   1.00 F 11467494
   6,  1      66 nlntk IYAISNDSHFVQ knwie      77   1.00 F 13186328
   8,  1      68 kknte VISVSEDTVYVH kawvq      79   1.00 F 13541053
   9,  1      66 kfkak VIGISVDSPFSL aefak      77   1.00 F 13541117
""".split('\n')
        self.motif_b_lines = """                          MOTIF b

15 columns
Num Motifs: 6
   2,  1     161 riles IQYVKENPGYACPVNWNFG dqvfy     179   1.00 F 11467494
  47,  1     160 lrmvd ALQFHEEHGDVCPAQWEKG kegmn     178   1.00 F 16501671
  67,  1     154 rkika AQYVAAHPGEVCPAKWKEG eatla     172   1.00 F 20151112
  81,  1     166 lrvvi SLQLTAEKRVATPVDWKDG dsvmv     184   1.00 F 3318841
  87,  1     163 lrvlk SLQLTNTHPVATPVNWKEG dkcci     181   1.00 F 4996210
  95,  1     160 lrlvq AFQYTDEHGEVCPAGWKPG sdtik     178   1.00 F 9955016
                       **** * ******* ** *

Log Motif portion of MAP for motif b = -187.76179
Log Fragmentation portion of MAP for motif b = -7.77486

-------------------------------------------------------------------------
""".split('\n')
    
    def test_get_sequence_and_motif_blocks(self):
        """get_sequence_and_motif_blocks tests."""
        seq_motif_lines = ['before line',\
                           '=====MAP MAXIMIZATION RESULTS=====',\
                           'after line'\
                           ]
        exp_seq_block=['before line']
        exp_motif_block=['=====MAP MAXIMIZATION RESULTS=====','after line']
        seq_block,motif_block = get_sequence_and_motif_blocks(seq_motif_lines)
        self.assertEqual(seq_block,exp_seq_block)
        self.assertEqual(motif_block,exp_motif_block)
    
    def test_get_sequence_map(self):
        """get_sequence_map tests."""
        sequence_map = get_sequence_map(self.gibbs_lines)
        self.assertEqual(sequence_map,self.sequence_map)
    
    def test_get_motif_blocks(self):
        """get_motif_blocks tests."""
        motif_lines = ['before motifs',\
                       'first MOTIF a',\
                       ' motif a data',\
                       'second MOTIF b',\
                       ' motif b data',
                       'after motifs'
                       ]
        exp_motif_blocks = [['first MOTIF a', 'motif a data'],\
                            ['second MOTIF b', 'motif b data', 'after motifs']]
        motif_blocks = get_motif_blocks(motif_lines)
        self.assertEqual(motif_blocks,exp_motif_blocks)
    
    def test_get_motif_sequences(self):
        """get_motif_sequences tests."""
        motif_list = get_motif_sequences(self.motif_a_lines)
        exp_motif_list = [('2', 71, 'ILAISVDSPFSH', 1.0, '1'),\
                          ('6', 65, 'IYAISNDSHFVQ', 1.0, '1'),\
                          ('8', 67, 'VISVSEDTVYVH', 1.0, '1'),\
                          ('9', 65, 'VIGISVDSPFSL', 1.0, '1')]
        self.assertEqual(motif_list,exp_motif_list)
    
    def test_get_motif_p_value(self):
        """get_motif_p_value tests."""
        log_list = ['Column 8 :  Sequence Description from Fast A input',\
                    'Log Motif portion of MAP for motif a = -469.15170',\
                    'Log Fragmentation portion of MAP for motif a = -3.80666',\
                    ]
        exp_p_val = exp(-469.15170)
        self.assertEqual(get_motif_p_value(log_list),exp_p_val)

    def test_guess_alphabet(self):
        """guess_alphabet tests."""
        motif_list = [('2', 71, 'ILAISVDSPFSH', 1.0, '1'),\
                      ('6', 65, 'IYAISNDSHFVQ', 1.0, '1'),\
                      ('8', 67, 'VISVSEDTVYVH', 1.0, '1'),\
                      ('9', 65, 'VIGISVDSPFSL', 1.0, '1')]
        alphabet = guess_alphabet(motif_list)
        self.assertEqual(alphabet,PROTEIN)
    
    def test_build_module_objects(self):
        """build_module_objects tests."""
        module = list(build_module_objects(self.motif_b_lines,\
            self.sequence_map))[0]
        exp_module_dict = {('20151112', 153): 'AQYVAAHPGEVCPAKWKEG',\
                           ('9955016', 159): 'AFQYTDEHGEVCPAGWKPG',\
                           ('16501671', 159): 'ALQFHEEHGDVCPAQWEKG',\
                           ('11467494', 160): 'IQYVKENPGYACPVNWNFG',\
                           ('4996210', 162): 'SLQLTNTHPVATPVNWKEG',\
                           ('3318841', 165): 'SLQLTAEKRVATPVDWKDG',\
                           }
        #module.AlignedSeqs.items() == exp_module_dict.items()
        for k1,k2 in zip(module.AlignedSeqs.keys(),\
                               exp_module_dict.keys()):
            self.assertEqual(k1,k2)
            v1 = str(module.AlignedSeqs[k1])
            v2 = exp_module_dict[k2]
            self.assertEqual(v1,v2)
    
    def test_module_ids_to_int(self):
        """module_ids_to_int tests."""
        module = list(build_module_objects(self.motif_b_lines,\
            self.sequence_map))[0]
        module_ids_to_int([module])
        self.assertEqual(module.ID,'0')
       

GIBBS_FILE = """
Gibbs.linux superfamily_aln_gis.fasta 10,15,20,25 5,5,5,5 
i = 20 range = 20 high = 11 low = -9

Gibbs 2.06.024  Jul 21 2005
Data file: superfamily_aln_gis.fasta
Current directory: /home/widmannj/superfamilies
The following options are set:
Concentrated Region          False    Sequence type        False
Collapsed Alphabet           False    Pseudocount weight   False
Use Expectation/Maximization False    Don't Xnu sequence   False
Help flag                    False    Near optimal cutoff  False
Number of iterations         False    Don't fragment       False
Don't use map maximization   False    Repeat regions       False
Output file                  False    Informed priors file False
Plateau periods              False    palindromic sequence False
Don't Reverse complement     False    Number of seeds      False
Seed Value                   False    Pseudosite weight    False
Suboptimal sampler output    False    Overlap              False
Allow width to vary          False    Wilcoxon signed rank False
Sample along length          False    Output Scan File     False
Output prior file            False    Modular Sampler      False
Ignore Spacing Model         False    Sample Background    False
Bkgnd Comp Model             False    Init from prior      False
Homologous Seq pairs         False    Parallel Tempering   False
Group Sampler                False    No progress info     False
Fragment from middle         False    Verify Mode          False
Alternate sample on k        False    No freq. soln.       False
Calc. def. pseudo wt.        False    Motif/Recur smpl     False
Phylogenetic Sampling        False    Supress Near Opt.    False
Nearopt display cutoff       False

site_samp            =            0
nMotifLen            =           10, 15, 20, 25
nAlphaLen            =           20
nNumMotifs           =            5 ,5 ,5 ,5
dPseudoCntWt         =          0.1
dPseudoSiteWt        =          0.8
nMaxIterations       =          500
lSeedVal             =   1149743202
nPlateauPeriods      =           20
nSeeds               =           10
nNumMotifTypes       =            4
dCutoff              =         0.01
dNearOptDispCutoff   =          0.5
RevComplement        =            0
glOverlapParam       =            0
Rcutoff factor       =        0.001
Post Plateau Samples =            0
Frag/Shft Per.       =            5
Frag width           =           15,22,30,37


Sequences to be Searched:
_________________________
#1   1091044
#2   11467494
#3   11499727
#4   1174686
#5   12044976
#6   13186328
#7   13358154
#8   13541053
#9   13541117
#10  135765
#11  1388082
#12  140543
#13  14286173
#14  14578634
#15  14600438
#16  15218394
#17  15597673
#18  15599256
#19  15602312
#20  15605725
#21  15605963
#22  15609375
#23  15609658
#24  15613511
#25  15614085
#26  15614140
#27  15615431
#28  15643152
#29  15672286
#30  15790738
#31  15791337
#32  15801846
#33  15805225
#34  15805374
#35  15807234
#36  15826629
#37  15899007
#38  15899339
#39  15964668
#40  15966937
#41  15988313
#42  16078864
#43  16123427
#44  16125919
#45  16330420
#46  1633495
#47  16501671
#48  1651717
#49  16759994
#50  16761507
#51  16803644
#52  16804867
#53  17229033
#54  17229859
#55  1729944
#56  17531233
#57  17537401
#58  17547503
#59  18309723
#60  18313548
#61  18406743
#62  19173077
#63  19554157
#64  19705357
#65  19746502
#66  20092028
#67  20151112
#68  21112072
#69  21222859
#70  21223405
#71  21227878
#72  21283385
#73  21674812
#74  23098307
#75  2649838
#76  267116
#77  27375582
#78  2822332
#79  30021713
#80  3261501
#81  3318841
#82  3323237
#83  4155972
#84  4200327
#85  4433065
#86  4704732
#87  4996210
#88  5326864
#89  6322180
#90  6323138
#91  6687568
#92  6850955
#93  7109697
#94  7290567
#95  9955016
#96  15677788
Processed Sequence Length: 16216 Total sequence length: 16307

Seed = 1149743202






motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 1 **

1
2
3
4[] motif A cycle 4 AP 0.0 (0 sites)
[] motif B cycle 4 AP -567.7 (18 sites)
[] motif C cycle 4 AP -1161.7 (31 sites)
[] motif D cycle 4 AP -245.0 (4 sites)
Total Map : 412.899 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 53

5[] motif A cycle 5 AP -26.6 (1 sites)
[] motif B cycle 5 AP -426.5 (17 sites)
[] motif C cycle 5 AP -1245.5 (33 sites)
[------] motif D cycle 5 AP -210.4 (4 sites)
Total Map : 499.315 Prev: 412.899 Diff: 86.4157 Motifs: 55

6[] motif A cycle 6 AP -178.0 (10 sites)
[] motif B cycle 6 AP -691.2 (26 sites)
[] motif C cycle 6 AP -1189.4 (32 sites)
[------] motif D cycle 6 AP -315.4 (6 sites)
Total Map : 605.7 Prev: 499.315 Diff: 106.385 Motifs: 74

7[] motif A cycle 7 AP -260.3 (16 sites)
[] motif B cycle 7 AP -664.9 (25 sites)
[] motif C cycle 7 AP -1189.4 (32 sites)
[------] motif D cycle 7 AP -366.4 (7 sites)
Total Map : 646.637 Prev: 605.7 Diff: 40.937 Motifs: 80

8[] motif A cycle 8 AP -331.5 (20 sites)
[] motif B cycle 8 AP -725.3 (27 sites)
[] motif C cycle 8 AP -1141.6 (31 sites)
[------] motif D cycle 8 AP -366.4 (7 sites)
Total Map : 660.38 Prev: 646.637 Diff: 13.7434 Motifs: 85

9
10[] motif A cycle 10 AP -371.3 (22 sites)
[] motif B cycle 10 AP -727.7 (27 sites)
[] motif C cycle 10 AP -1245.5 (33 sites)
[] motif D cycle 10 AP -346.6 (7 sites)
Total Map : 671.561 Prev: 660.38 Diff: 11.181 Motifs: 89

11[] motif A cycle 11 AP -365.7 (22 sites)
[] motif B cycle 11 AP -760.1 (28 sites)
[] motif C cycle 11 AP -1189.4 (32 sites)
[] motif D cycle 11 AP -346.6 (7 sites)
Total Map : 685.28 Prev: 671.561 Diff: 13.719 Motifs: 89

12[] motif A cycle 12 AP -446.5 (26 sites)
[] motif B cycle 12 AP -725.3 (27 sites)
[] motif C cycle 12 AP -1141.6 (31 sites)
[] motif D cycle 12 AP -346.6 (7 sites)
Total Map : 689.5 Prev: 685.28 Diff: 4.21965 Motifs: 91

13
14
15[] motif A cycle 15 AP -422.8 (26 sites)
[] motif B cycle 15 AP -726.1 (27 sites)
[] motif C cycle 15 AP -1041.5 (29 sites)
[] motif D cycle 15 AP -348.1 (7 sites)
Total Map : 714.246 Prev: 689.5 Diff: 24.7462 Motifs: 89

16[] motif A cycle 16 AP -440.3 (27 sites)
[] motif B cycle 16 AP -725.3 (27 sites)
[] motif C cycle 16 AP -1041.5 (29 sites)
[] motif D cycle 16 AP -348.1 (7 sites)
Total Map : 717.139 Prev: 714.246 Diff: 2.89264 Motifs: 90

17[] motif A cycle 17 AP -416.4 (26 sites)
[] motif B cycle 17 AP -725.3 (27 sites)
[] motif C cycle 17 AP -1041.5 (29 sites)
[] motif D cycle 17 AP -348.1 (7 sites)
Total Map : 723.199 Prev: 717.139 Diff: 6.05975 Motifs: 89

18[] motif A cycle 18 AP -476.5 (29 sites)
[] motif B cycle 18 AP -660.6 (25 sites)
[] motif C cycle 18 AP -1041.5 (29 sites)
[] motif D cycle 18 AP -348.1 (7 sites)
Total Map : 725.194 Prev: 723.199 Diff: 1.99565 Motifs: 90

19[] motif A cycle 19 AP -387.7 (25 sites)
[] motif B cycle 19 AP -725.3 (27 sites)
[] motif C cycle 19 AP -1041.5 (29 sites)
[] motif D cycle 19 AP -348.1 (7 sites)
Total Map : 730.839 Prev: 725.194 Diff: 5.64428 Motifs: 88

20[] motif A cycle 20 AP -454.8 (28 sites)
[] motif B cycle 20 AP -691.2 (26 sites)
[] motif C cycle 20 AP -1090.2 (30 sites)
[] motif D cycle 20 AP -348.1 (7 sites)
Total Map : 740.315 Prev: 730.839 Diff: 9.47605 Motifs: 91

21
22[] motif A cycle 22 AP -347.4 (23 sites)
[] motif B cycle 22 AP -691.2 (26 sites)
[] motif C cycle 22 AP -1041.0 (29 sites)
[] motif D cycle 22 AP -348.1 (7 sites)
Total Map : 742.668 Prev: 740.315 Diff: 2.35327 Motifs: 85

23
24[] motif A cycle 24 AP -368.0 (24 sites)
[] motif B cycle 24 AP -728.8 (27 sites)
[] motif C cycle 24 AP -1041.0 (29 sites)
[] motif D cycle 24 AP -348.1 (7 sites)
Total Map : 742.863 Prev: 742.668 Diff: 0.19554 Motifs: 87

25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
MAX :: 742.863379 (Seed = 1149743202, Iteration = 24   Motif A = 24 Motif B = 27 Motif C = 29 Motif D = 7 )
motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 2 **

1
2
3
4[] motif A cycle 4 AP -53.8 (3 sites)
[] motif B cycle 4 AP -751.2 (33 sites)
[] motif C cycle 4 AP -1873.4 (54 sites)
[] motif D cycle 4 AP -667.4 (12 sites)
Total Map : 1679.61 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 102

5[] motif A cycle 5 AP -51.0 (3 sites)
[] motif B cycle 5 AP -684.4 (33 sites)
[] motif C cycle 5 AP -1599.1 (54 sites)
[------] motif D cycle 5 AP -534.8 (12 sites)
Total Map : 2005.21 Prev: 1679.61 Diff: 325.602 Motifs: 102

6[] motif A cycle 6 AP -71.4 (4 sites)
[] motif B cycle 6 AP -741.8 (35 sites)
[] motif C cycle 6 AP -1599.1 (54 sites)
[------] motif D cycle 6 AP -468.9 (11 sites)
Total Map : 2020.82 Prev: 2005.21 Diff: 15.6117 Motifs: 104

7[] motif A cycle 7 AP -51.0 (3 sites)
[] motif B cycle 7 AP -741.8 (35 sites)
[] motif C cycle 7 AP -1599.1 (54 sites)
[------] motif D cycle 7 AP -468.9 (11 sites)
Total Map : 2022.05 Prev: 2020.82 Diff: 1.23593 Motifs: 103

8
9
10[] motif A cycle 10 AP -48.5 (3 sites)
[] motif B cycle 10 AP -741.8 (35 sites)
[] motif C cycle 10 AP -1599.1 (54 sites)
[] motif D cycle 10 AP -534.8 (12 sites)
Total Map : 2026.1 Prev: 2022.05 Diff: 4.04393 Motifs: 104

11
12
13
14
15[] motif A cycle 15 AP -47.4 (3 sites)
[] motif B cycle 15 AP -741.8 (35 sites)
[] motif C cycle 15 AP -1599.1 (54 sites)
[] motif D cycle 15 AP -534.8 (12 sites)
Total Map : 2026.45 Prev: 2026.1 Diff: 0.353957 Motifs: 104

16
17
18[] motif A cycle 18 AP -64.6 (4 sites)
[] motif B cycle 18 AP -741.8 (35 sites)
[] motif C cycle 18 AP -1599.1 (54 sites)
[] motif D cycle 18 AP -468.9 (11 sites)
Total Map : 2028.63 Prev: 2026.45 Diff: 2.1753 Motifs: 104

19
20
21
22[] motif A cycle 22 AP -64.6 (4 sites)
[] motif B cycle 22 AP -850.9 (38 sites)
[] motif C cycle 22 AP -1599.1 (54 sites)
[] motif D cycle 22 AP -468.9 (11 sites)
Total Map : 2029.4 Prev: 2028.63 Diff: 0.769306 Motifs: 107

23
24
25[] motif A cycle 25 AP -64.6 (4 sites)
[] motif B cycle 25 AP -850.9 (38 sites)
[] motif C cycle 25 AP -1599.1 (54 sites)
[] motif D cycle 25 AP -534.8 (12 sites)
Total Map : 2029.8 Prev: 2029.4 Diff: 0.400632 Motifs: 108

26
27
28
29
30[] motif A cycle 30 AP -47.5 (3 sites)
[] motif B cycle 30 AP -850.9 (38 sites)
[] motif C cycle 30 AP -1658.1 (55 sites)
[] motif D cycle 30 AP -467.4 (11 sites)
Total Map : 2035.48 Prev: 2029.8 Diff: 5.68003 Motifs: 107

31
32
33
34[] motif A cycle 34 AP -47.5 (3 sites)
[] motif B cycle 34 AP -816.5 (37 sites)
[] motif C cycle 34 AP -1599.1 (54 sites)
[] motif D cycle 34 AP -467.4 (11 sites)
Total Map : 2036.53 Prev: 2035.48 Diff: 1.0477 Motifs: 105

35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
MAX :: 2036.525370 (Seed = 1149743202, Iteration = 34   Motif A = 3 Motif B = 37 Motif C = 54 Motif D = 11 )
motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 3 **

1
2
3
4[] motif A cycle 4 AP -408.4 (26 sites)
[] motif B cycle 4 AP -71.4 (2 sites)
[] motif C cycle 4 AP -1871.1 (53 sites)
[] motif D cycle 4 AP -174.1 (3 sites)
Total Map : 1146.38 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 84

5[] motif A cycle 5 AP -296.5 (26 sites)
[] motif B cycle 5 AP -71.4 (2 sites)
[] motif C cycle 5 AP -1665.3 (53 sites)
[] motif D cycle 5 AP -174.1 (3 sites)
Total Map : 1455.53 Prev: 1146.38 Diff: 309.149 Motifs: 84

6[] motif A cycle 6 AP -343.8 (29 sites)
[] motif B cycle 6 AP -71.4 (2 sites)
[] motif C cycle 6 AP -1696.2 (54 sites)
[] motif D cycle 6 AP -174.1 (3 sites)
Total Map : 1498.45 Prev: 1455.53 Diff: 42.9205 Motifs: 88

7
8[] motif A cycle 8 AP -384.0 (31 sites)
[] motif B cycle 8 AP -71.4 (2 sites)
[] motif C cycle 8 AP -1696.2 (54 sites)
[] motif D cycle 8 AP -174.1 (3 sites)
Total Map : 1502.55 Prev: 1498.45 Diff: 4.0994 Motifs: 90

9[] motif A cycle 9 AP -422.5 (33 sites)
[] motif B cycle 9 AP -71.4 (2 sites)
[] motif C cycle 9 AP -1696.2 (54 sites)
[] motif D cycle 9 AP -174.1 (3 sites)
Total Map : 1503.54 Prev: 1502.55 Diff: 0.99024 Motifs: 92

10
11[] motif A cycle 11 AP -419.9 (34 sites)
[] motif B cycle 11 AP -71.4 (2 sites)
[] motif C cycle 11 AP -1753.3 (55 sites)
[] motif D cycle 11 AP -174.1 (3 sites)
Total Map : 1516.48 Prev: 1503.54 Diff: 12.938 Motifs: 94

12[] motif A cycle 12 AP -419.9 (34 sites)
[] motif B cycle 12 AP -71.4 (2 sites)
[] motif C cycle 12 AP -1696.2 (54 sites)
[] motif D cycle 12 AP -174.1 (3 sites)
Total Map : 1518.67 Prev: 1516.48 Diff: 2.18498 Motifs: 93

13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 4 **

1
2
3
4[] motif A cycle 4 AP 0.0 (0 sites)
[] motif B cycle 4 AP -1240.2 (53 sites)
[] motif C cycle 4 AP -192.8 (5 sites)
[] motif D cycle 4 AP -1139.0 (22 sites)
Total Map : 1140.99 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 80

5[] motif A cycle 5 AP -27.5 (1 sites)
[] motif B cycle 5 AP -1108.9 (54 sites)
[] motif C cycle 5 AP -107.7 (4 sites)
[+++] motif D cycle 5 AP -1112.4 (24 sites)
Total Map : 1548.32 Prev: 1140.99 Diff: 407.325 Motifs: 83

6
7
8
9[] motif A cycle 9 AP -370.2 (22 sites)
[] motif B cycle 9 AP -1103.5 (54 sites)
[] motif C cycle 9 AP -1257.6 (33 sites)
[+++] motif D cycle 9 AP 0.0 (0 sites)
Total Map : 1561.07 Prev: 1548.32 Diff: 12.75 Motifs: 109

10[] motif A cycle 10 AP -434.4 (26 sites)
[] motif B cycle 10 AP -1071.7 (53 sites)
[] motif C cycle 10 AP -1241.9 (34 sites)
[] motif D cycle 10 AP -75.4 (1 sites)
Total Map : 1606.41 Prev: 1561.07 Diff: 45.3433 Motifs: 114

11
12[] motif A cycle 12 AP -551.0 (33 sites)
[] motif B cycle 12 AP -1144.9 (55 sites)
[] motif C cycle 12 AP -1232.3 (34 sites)
[] motif D cycle 12 AP 0.0 (0 sites)
Total Map : 1660.34 Prev: 1606.41 Diff: 53.9249 Motifs: 122

13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 5 **

1
2
3
4[] motif A cycle 4 AP -44.7 (2 sites)
[] motif B cycle 4 AP -1352.6 (53 sites)
[] motif C cycle 4 AP -98.3 (2 sites)
[] motif D cycle 4 AP -1454.8 (28 sites)
Total Map : 989.762 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 85

5[] motif A cycle 5 AP -61.7 (3 sites)
[] motif B cycle 5 AP -1093.0 (53 sites)
[] motif C cycle 5 AP -82.2 (2 sites)
[-----] motif D cycle 5 AP -1250.5 (29 sites)
Total Map : 1543.73 Prev: 989.762 Diff: 553.967 Motifs: 87

6[] motif A cycle 6 AP -117.1 (6 sites)
[] motif B cycle 6 AP -1118.4 (54 sites)
[] motif C cycle 6 AP -82.2 (2 sites)
[-----] motif D cycle 6 AP -1228.9 (29 sites)
Total Map : 1579.26 Prev: 1543.73 Diff: 35.536 Motifs: 91

7
8
9
10[] motif A cycle 10 AP -119.2 (6 sites)
[] motif B cycle 10 AP -1156.4 (55 sites)
[] motif C cycle 10 AP -81.4 (2 sites)
[] motif D cycle 10 AP -1227.7 (29 sites)
Total Map : 1584.6 Prev: 1579.26 Diff: 5.33507 Motifs: 92

11[] motif A cycle 11 AP -117.1 (6 sites)
[] motif B cycle 11 AP -1117.5 (54 sites)
[] motif C cycle 11 AP -81.4 (2 sites)
[] motif D cycle 11 AP -1227.7 (29 sites)
Total Map : 1584.81 Prev: 1584.6 Diff: 0.213828 Motifs: 91

12[] motif A cycle 12 AP -117.1 (6 sites)
[] motif B cycle 12 AP -1156.4 (55 sites)
[] motif C cycle 12 AP -81.4 (2 sites)
[] motif D cycle 12 AP -1227.7 (29 sites)
Total Map : 1587.59 Prev: 1584.81 Diff: 2.77506 Motifs: 92

13
14
15
16
17
18
19
20[] motif A cycle 20 AP -97.4 (5 sites)
[] motif B cycle 20 AP -1156.8 (55 sites)
[] motif C cycle 20 AP -79.9 (2 sites)
[] motif D cycle 20 AP -1227.7 (29 sites)
Total Map : 1587.99 Prev: 1587.59 Diff: 0.397059 Motifs: 91

21
22
23
24
25
26
27
28
29
30[] motif A cycle 30 AP -39.0 (2 sites)
[] motif B cycle 30 AP -1156.8 (55 sites)
[] motif C cycle 30 AP -80.5 (2 sites)
[] motif D cycle 30 AP -1227.7 (29 sites)
Total Map : 1588.16 Prev: 1587.99 Diff: 0.175261 Motifs: 88

31
32[] motif A cycle 32 AP -58.3 (3 sites)
[] motif B cycle 32 AP -1156.8 (55 sites)
[] motif C cycle 32 AP -80.5 (2 sites)
[] motif D cycle 32 AP -1227.7 (29 sites)
Total Map : 1588.59 Prev: 1588.16 Diff: 0.429544 Motifs: 89

33
34
35
36[] motif A cycle 36 AP -232.0 (13 sites)
[] motif B cycle 36 AP -1117.1 (54 sites)
[] motif C cycle 36 AP 0.0 (0 sites)
[] motif D cycle 36 AP -1227.7 (29 sites)
Total Map : 1601.78 Prev: 1588.59 Diff: 13.1946 Motifs: 96

37[] motif A cycle 37 AP -229.0 (13 sites)
[] motif B cycle 37 AP -1156.8 (55 sites)
[] motif C cycle 37 AP -234.4 (5 sites)
[] motif D cycle 37 AP -1227.7 (29 sites)
Total Map : 1611.86 Prev: 1601.78 Diff: 10.0767 Motifs: 102

38[] motif A cycle 38 AP -247.4 (14 sites)
[] motif B cycle 38 AP -1156.8 (55 sites)
[] motif C cycle 38 AP -544.7 (13 sites)
[] motif D cycle 38 AP -1227.7 (29 sites)
Total Map : 1688.66 Prev: 1611.86 Diff: 76.798 Motifs: 111

39
40[] motif A cycle 40 AP -199.9 (12 sites)
[] motif B cycle 40 AP -1156.8 (55 sites)
[] motif C cycle 40 AP -692.6 (17 sites)
[] motif D cycle 40 AP -1227.7 (29 sites)
Total Map : 1767.8 Prev: 1688.66 Diff: 79.137 Motifs: 113

41[] motif A cycle 41 AP -311.4 (18 sites)
[] motif B cycle 41 AP -1156.8 (55 sites)
[] motif C cycle 41 AP -784.5 (19 sites)
[] motif D cycle 41 AP -1227.7 (29 sites)
Total Map : 1785.93 Prev: 1767.8 Diff: 18.1288 Motifs: 121

42
43
44
45
46[] motif A cycle 46 AP -483.4 (27 sites)
[] motif B cycle 46 AP -1120.7 (54 sites)
[] motif C cycle 46 AP -928.7 (22 sites)
[] motif D cycle 46 AP -1227.7 (29 sites)
Total Map : 1799.6 Prev: 1785.93 Diff: 13.6793 Motifs: 132

47[] motif A cycle 47 AP -529.5 (29 sites)
[] motif B cycle 47 AP -1156.8 (55 sites)
[] motif C cycle 47 AP -880.8 (21 sites)
[] motif D cycle 47 AP -1227.7 (29 sites)
Total Map : 1804.67 Prev: 1799.6 Diff: 5.06939 Motifs: 134

48[] motif A cycle 48 AP -459.6 (26 sites)
[] motif B cycle 48 AP -1118.4 (54 sites)
[] motif C cycle 48 AP -876.6 (21 sites)
[] motif D cycle 48 AP -1227.7 (29 sites)
Total Map : 1808.9 Prev: 1804.67 Diff: 4.22417 Motifs: 130

49
50
51
52[] motif A cycle 52 AP -487.1 (27 sites)
[] motif B cycle 52 AP -1156.8 (55 sites)
[] motif C cycle 52 AP -914.8 (22 sites)
[] motif D cycle 52 AP -1227.7 (29 sites)
Total Map : 1817.39 Prev: 1808.9 Diff: 8.49225 Motifs: 133

53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 6 **

1
2
3
4[] motif A cycle 4 AP -67.0 (3 sites)
[] motif B cycle 4 AP -1223.6 (54 sites)
[] motif C cycle 4 AP -145.7 (3 sites)
[] motif D cycle 4 AP -224.8 (4 sites)
Total Map : 981.16 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 64

5[] motif A cycle 5 AP -37.3 (2 sites)
[] motif B cycle 5 AP -989.0 (54 sites)
[] motif C cycle 5 AP -304.6 (8 sites)
[----] motif D cycle 5 AP -187.9 (4 sites)
Total Map : 1275.57 Prev: 981.16 Diff: 294.41 Motifs: 68

6[] motif A cycle 6 AP -81.8 (5 sites)
[] motif B cycle 6 AP -989.0 (54 sites)
[] motif C cycle 6 AP -590.1 (16 sites)
[----] motif D cycle 6 AP -358.3 (7 sites)
Total Map : 1416.85 Prev: 1275.57 Diff: 141.283 Motifs: 82

7[] motif A cycle 7 AP -96.8 (6 sites)
[] motif B cycle 7 AP -989.0 (54 sites)
[] motif C cycle 7 AP -633.2 (17 sites)
[----] motif D cycle 7 AP -358.3 (7 sites)
Total Map : 1424.61 Prev: 1416.85 Diff: 7.75982 Motifs: 84

8
9
10[] motif A cycle 10 AP -91.9 (6 sites)
[] motif B cycle 10 AP -989.0 (54 sites)
[] motif C cycle 10 AP -674.5 (18 sites)
[-] motif D cycle 10 AP -337.3 (7 sites)
Total Map : 1448.23 Prev: 1424.61 Diff: 23.6156 Motifs: 85

11[] motif A cycle 11 AP -91.9 (6 sites)
[] motif B cycle 11 AP -989.0 (54 sites)
[] motif C cycle 11 AP -667.6 (18 sites)
[-] motif D cycle 11 AP -337.3 (7 sites)
Total Map : 1456.26 Prev: 1448.23 Diff: 8.03531 Motifs: 85

12
13[] motif A cycle 13 AP -114.3 (7 sites)
[] motif B cycle 13 AP -989.0 (54 sites)
[] motif C cycle 13 AP -854.9 (22 sites)
[-] motif D cycle 13 AP -337.3 (7 sites)
Total Map : 1463.23 Prev: 1456.26 Diff: 6.96803 Motifs: 90

14[] motif A cycle 14 AP -137.2 (8 sites)
[] motif B cycle 14 AP -989.0 (54 sites)
[] motif C cycle 14 AP -803.3 (21 sites)
[-] motif D cycle 14 AP -337.3 (7 sites)
Total Map : 1465.75 Prev: 1463.23 Diff: 2.52048 Motifs: 90

15[] motif A cycle 15 AP -134.8 (8 sites)
[] motif B cycle 15 AP -989.0 (54 sites)
[] motif C cycle 15 AP -803.3 (21 sites)
[+] motif D cycle 15 AP -336.6 (7 sites)
Total Map : 1470.32 Prev: 1465.75 Diff: 4.56864 Motifs: 90

16[] motif A cycle 16 AP -108.3 (7 sites)
[] motif B cycle 16 AP -989.0 (54 sites)
[] motif C cycle 16 AP -850.5 (22 sites)
[+] motif D cycle 16 AP -411.5 (8 sites)
Total Map : 1480.76 Prev: 1470.32 Diff: 10.4373 Motifs: 91

17
18
19
20[] motif A cycle 20 AP -128.6 (8 sites)
[] motif B cycle 20 AP -989.0 (54 sites)
[] motif C cycle 20 AP -940.1 (24 sites)
[--] motif D cycle 20 AP -405.0 (8 sites)
Total Map : 1500.43 Prev: 1480.76 Diff: 19.6736 Motifs: 94

21[] motif A cycle 21 AP -108.3 (7 sites)
[] motif B cycle 21 AP -989.0 (54 sites)
[] motif C cycle 21 AP -940.1 (24 sites)
[--] motif D cycle 21 AP -405.0 (8 sites)
Total Map : 1501.82 Prev: 1500.43 Diff: 1.38679 Motifs: 93

22[] motif A cycle 22 AP -108.3 (7 sites)
[] motif B cycle 22 AP -989.0 (54 sites)
[] motif C cycle 22 AP -986.8 (25 sites)
[--] motif D cycle 22 AP -405.0 (8 sites)
Total Map : 1503.53 Prev: 1501.82 Diff: 1.70712 Motifs: 94

23
24
25
26
27[] motif A cycle 27 AP -128.6 (8 sites)
[] motif B cycle 27 AP -989.0 (54 sites)
[] motif C cycle 27 AP -939.0 (24 sites)
[] motif D cycle 27 AP -404.6 (8 sites)
Total Map : 1504 Prev: 1503.53 Diff: 0.475394 Motifs: 94

28[] motif A cycle 28 AP -108.3 (7 sites)
[] motif B cycle 28 AP -989.0 (54 sites)
[] motif C cycle 28 AP -939.0 (24 sites)
[] motif D cycle 28 AP -404.6 (8 sites)
Total Map : 1505.38 Prev: 1504 Diff: 1.37798 Motifs: 93

29
30
31
32
33
34
35[] motif A cycle 35 AP -142.7 (9 sites)
[] motif B cycle 35 AP -989.0 (54 sites)
[] motif C cycle 35 AP -889.8 (23 sites)
[] motif D cycle 35 AP -471.7 (9 sites)
Total Map : 1508.79 Prev: 1505.38 Diff: 3.40696 Motifs: 95

36
37[] motif A cycle 37 AP -159.8 (10 sites)
[] motif B cycle 37 AP -989.0 (54 sites)
[] motif C cycle 37 AP -939.0 (24 sites)
[] motif D cycle 37 AP -404.6 (8 sites)
Total Map : 1514.12 Prev: 1508.79 Diff: 5.3377 Motifs: 96

38[] motif A cycle 38 AP -159.8 (10 sites)
[] motif B cycle 38 AP -989.0 (54 sites)
[] motif C cycle 38 AP -986.4 (25 sites)
[] motif D cycle 38 AP -404.6 (8 sites)
Total Map : 1514.79 Prev: 1514.12 Diff: 0.670946 Motifs: 97

39
40[] motif A cycle 40 AP -175.3 (11 sites)
[] motif B cycle 40 AP -989.0 (54 sites)
[] motif C cycle 40 AP -939.0 (24 sites)
[] motif D cycle 40 AP -404.6 (8 sites)
Total Map : 1518.92 Prev: 1514.79 Diff: 4.12608 Motifs: 97

41[] motif A cycle 41 AP -151.0 (10 sites)
[] motif B cycle 41 AP -989.0 (54 sites)
[] motif C cycle 41 AP -937.3 (24 sites)
[] motif D cycle 41 AP -404.6 (8 sites)
Total Map : 1522.39 Prev: 1518.92 Diff: 3.46783 Motifs: 96

42[] motif A cycle 42 AP -132.3 (9 sites)
[] motif B cycle 42 AP -989.0 (54 sites)
[] motif C cycle 42 AP -986.4 (25 sites)
[] motif D cycle 42 AP -404.6 (8 sites)
Total Map : 1523.36 Prev: 1522.39 Diff: 0.973725 Motifs: 96

43
44
45
46
47
48
49
50[] motif A cycle 50 AP -151.0 (10 sites)
[] motif B cycle 50 AP -989.0 (54 sites)
[] motif C cycle 50 AP -986.4 (25 sites)
[] motif D cycle 50 AP -404.6 (8 sites)
Total Map : 1524.12 Prev: 1523.36 Diff: 0.756895 Motifs: 97

51
52
53
54[] motif A cycle 54 AP -169.9 (11 sites)
[] motif B cycle 54 AP -989.0 (54 sites)
[] motif C cycle 54 AP -939.0 (24 sites)
[] motif D cycle 54 AP -404.6 (8 sites)
Total Map : 1524.49 Prev: 1524.12 Diff: 0.371963 Motifs: 97

55
56
57
58
59
60
61
62
63
64
65[] motif A cycle 65 AP -151.0 (10 sites)
[] motif B cycle 65 AP -989.0 (54 sites)
[] motif C cycle 65 AP -986.8 (25 sites)
[] motif D cycle 65 AP -402.9 (8 sites)
Total Map : 1529.11 Prev: 1524.49 Diff: 4.61554 Motifs: 97

66[] motif A cycle 66 AP -169.9 (11 sites)
[] motif B cycle 66 AP -989.0 (54 sites)
[] motif C cycle 66 AP -986.8 (25 sites)
[] motif D cycle 66 AP -402.9 (8 sites)
Total Map : 1530.16 Prev: 1529.11 Diff: 1.05073 Motifs: 98

67
68
69
70[] motif A cycle 70 AP -169.9 (11 sites)
[] motif B cycle 70 AP -989.0 (54 sites)
[] motif C cycle 70 AP -986.4 (25 sites)
[] motif D cycle 70 AP -402.9 (8 sites)
Total Map : 1532.06 Prev: 1530.16 Diff: 1.89983 Motifs: 98

71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 7 **

1
2
3
4[] motif A cycle 4 AP -201.0 (12 sites)
[] motif B cycle 4 AP -113.6 (3 sites)
[] motif C cycle 4 AP -930.1 (23 sites)
[] motif D cycle 4 AP -2959.1 (54 sites)
Total Map : 903.055 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 92

5[] motif A cycle 5 AP -150.4 (11 sites)
[] motif B cycle 5 AP -117.6 (4 sites)
[] motif C cycle 5 AP -871.9 (23 sites)
[] motif D cycle 5 AP -2452.5 (52 sites)
Total Map : 1437.59 Prev: 903.055 Diff: 534.536 Motifs: 90

6[] motif A cycle 6 AP -185.8 (13 sites)
[] motif B cycle 6 AP -63.1 (2 sites)
[] motif C cycle 6 AP -817.6 (22 sites)
[] motif D cycle 6 AP -2556.7 (54 sites)
Total Map : 1469.43 Prev: 1437.59 Diff: 31.8425 Motifs: 91

7[] motif A cycle 7 AP -185.8 (13 sites)
[] motif B cycle 7 AP -63.1 (2 sites)
[] motif C cycle 7 AP -817.6 (22 sites)
[] motif D cycle 7 AP -2469.2 (53 sites)
Total Map : 1490.24 Prev: 1469.43 Diff: 20.8077 Motifs: 90

8
9
10[] motif A cycle 10 AP -185.8 (13 sites)
[] motif B cycle 10 AP -60.9 (2 sites)
[] motif C cycle 10 AP -923.4 (24 sites)
[--] motif D cycle 10 AP -2454.7 (53 sites)
Total Map : 1505.14 Prev: 1490.24 Diff: 14.9023 Motifs: 92

11
12
13[] motif A cycle 13 AP -185.8 (13 sites)
[] motif B cycle 13 AP -60.9 (2 sites)
[] motif C cycle 13 AP -817.6 (22 sites)
[--] motif D cycle 13 AP -2454.7 (53 sites)
Total Map : 1505.52 Prev: 1505.14 Diff: 0.372181 Motifs: 90

14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 8 **

1
2
3
4[] motif A cycle 4 AP -56.8 (3 sites)
[] motif B cycle 4 AP -594.9 (23 sites)
[] motif C cycle 4 AP -194.1 (4 sites)
[] motif D cycle 4 AP -367.4 (6 sites)
Total Map : 271.959 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 36

5[] motif A cycle 5 AP -46.0 (3 sites)
[] motif B cycle 5 AP -486.2 (21 sites)
[] motif C cycle 5 AP -979.1 (29 sites)
[-----] motif D cycle 5 AP -276.6 (6 sites)
Total Map : 840.305 Prev: 271.959 Diff: 568.345 Motifs: 59

6[] motif A cycle 6 AP -46.0 (3 sites)
[] motif B cycle 6 AP -481.1 (20 sites)
[] motif C cycle 6 AP -1396.9 (39 sites)
[-----] motif D cycle 6 AP -486.7 (10 sites)
Total Map : 950.671 Prev: 840.305 Diff: 110.367 Motifs: 72

7
8[] motif A cycle 8 AP -46.0 (3 sites)
[] motif B cycle 8 AP -316.3 (10 sites)
[] motif C cycle 8 AP -3004.0 (74 sites)
[-----] motif D cycle 8 AP -552.2 (11 sites)
Total Map : 954.552 Prev: 950.671 Diff: 3.88022 Motifs: 98

9
10[] motif A cycle 10 AP -46.5 (3 sites)
[] motif B cycle 10 AP -243.9 (8 sites)
[] motif C cycle 10 AP -2992.1 (74 sites)
[-] motif D cycle 10 AP -530.1 (11 sites)
Total Map : 1001.93 Prev: 954.552 Diff: 47.3832 Motifs: 96

11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 9 **

1
2
3
4[] motif A cycle 4 AP -86.2 (4 sites)
[] motif B cycle 4 AP -122.1 (3 sites)
[] motif C cycle 4 AP -2149.3 (51 sites)
[] motif D cycle 4 AP -234.1 (4 sites)
Total Map : 509.201 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 62

5[] motif A cycle 5 AP -66.9 (4 sites)
[] motif B cycle 5 AP -75.7 (3 sites)
[] motif C cycle 5 AP -1697.1 (50 sites)
[] motif D cycle 5 AP -461.9 (9 sites)
Total Map : 1092.03 Prev: 509.201 Diff: 582.825 Motifs: 66

6[] motif A cycle 6 AP -311.5 (21 sites)
[] motif B cycle 6 AP -96.6 (3 sites)
[] motif C cycle 6 AP -1790.1 (52 sites)
[] motif D cycle 6 AP -742.9 (14 sites)
Total Map : 1213.73 Prev: 1092.03 Diff: 121.707 Motifs: 90

7[] motif A cycle 7 AP -368.8 (24 sites)
[] motif B cycle 7 AP -112.7 (4 sites)
[] motif C cycle 7 AP -1737.2 (51 sites)
[] motif D cycle 7 AP -1041.8 (19 sites)
Total Map : 1253.72 Prev: 1213.73 Diff: 39.9908 Motifs: 98

8
9
10[] motif A cycle 10 AP -481.8 (31 sites)
[] motif B cycle 10 AP -58.0 (2 sites)
[] motif C cycle 10 AP -1790.1 (52 sites)
[] motif D cycle 10 AP -1237.8 (23 sites)
Total Map : 1285.2 Prev: 1253.72 Diff: 31.4769 Motifs: 108

11[] motif A cycle 11 AP -489.8 (32 sites)
[] motif B cycle 11 AP -99.2 (3 sites)
[] motif C cycle 11 AP -1790.1 (52 sites)
[] motif D cycle 11 AP -1417.6 (26 sites)
Total Map : 1303.74 Prev: 1285.2 Diff: 18.5344 Motifs: 113

12[] motif A cycle 12 AP -492.2 (32 sites)
[] motif B cycle 12 AP -58.0 (2 sites)
[] motif C cycle 12 AP -1842.5 (53 sites)
[] motif D cycle 12 AP -1419.1 (26 sites)
Total Map : 1309.87 Prev: 1303.74 Diff: 6.13088 Motifs: 113

13
14
15
16
17[] motif A cycle 17 AP -469.5 (31 sites)
[] motif B cycle 17 AP -58.7 (2 sites)
[] motif C cycle 17 AP -1790.1 (52 sites)
[] motif D cycle 17 AP -1419.1 (26 sites)
Total Map : 1311.38 Prev: 1309.87 Diff: 1.50905 Motifs: 111

18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37motif A: 5 (+/- 7.88) out of 15393   a = 20; b = 61552; p = 0.000323771
motif B: 5 (+/- 7.87) out of 14888   a = 20; b = 59532; p = 0.000334158
motif C: 5 (+/- 7.86) out of 14383   a = 20; b = 57512; p = 0.000345232
motif D: 5 (+/- 7.85) out of 13878   a = 20; b = 55492; p = 0.000357066

** 10 **

1
2
3
4[] motif A cycle 4 AP -374.5 (22 sites)
[] motif B cycle 4 AP -74.4 (2 sites)
[] motif C cycle 4 AP -1774.4 (54 sites)
[] motif D cycle 4 AP -1794.7 (36 sites)
Total Map : 1821.24 Prev: -1.79769e+308 Diff: 1.79769e+308 Motifs: 114

5[] motif A cycle 5 AP -364.1 (22 sites)
[] motif B cycle 5 AP -74.4 (2 sites)
[] motif C cycle 5 AP -1607.6 (54 sites)
[++] motif D cycle 5 AP -1738.3 (36 sites)
Total Map : 2005.44 Prev: 1821.24 Diff: 184.204 Motifs: 114

6[] motif A cycle 6 AP -355.0 (22 sites)
[] motif B cycle 6 AP -74.4 (2 sites)
[] motif C cycle 6 AP -1607.6 (54 sites)
[++] motif D cycle 6 AP -1800.3 (37 sites)
Total Map : 2026.99 Prev: 2005.44 Diff: 21.5427 Motifs: 115

7[] motif A cycle 7 AP -305.7 (20 sites)
[] motif B cycle 7 AP -74.4 (2 sites)
[] motif C cycle 7 AP -1607.6 (54 sites)
[++] motif D cycle 7 AP -1800.3 (37 sites)
Total Map : 2042.47 Prev: 2026.99 Diff: 15.4851 Motifs: 113

8[] motif A cycle 8 AP -360.2 (23 sites)
[] motif B cycle 8 AP -74.4 (2 sites)
[] motif C cycle 8 AP -1607.6 (54 sites)
[++] motif D cycle 8 AP -1734.6 (36 sites)
Total Map : 2049.89 Prev: 2042.47 Diff: 7.41789 Motifs: 115

9
10[] motif A cycle 10 AP -321.7 (22 sites)
[] motif B cycle 10 AP -74.4 (2 sites)
[] motif C cycle 10 AP -1607.6 (54 sites)
[] motif D cycle 10 AP -1800.3 (37 sites)
Total Map : 2065.41 Prev: 2049.89 Diff: 15.5161 Motifs: 115

11[] motif A cycle 11 AP -345.0 (23 sites)
[] motif B cycle 11 AP -74.4 (2 sites)
[] motif C cycle 11 AP -1607.6 (54 sites)
[] motif D cycle 11 AP -1800.3 (37 sites)
Total Map : 2068.72 Prev: 2065.41 Diff: 3.31307 Motifs: 116

12
13
14
15
16[] motif A cycle 16 AP -467.4 (29 sites)
[] motif B cycle 16 AP -74.4 (2 sites)
[] motif C cycle 16 AP -1607.6 (54 sites)
[] motif D cycle 16 AP -1800.3 (37 sites)
Total Map : 2074.61 Prev: 2068.72 Diff: 5.89359 Motifs: 122

17[] motif A cycle 17 AP -469.2 (29 sites)
[] motif B cycle 17 AP -74.4 (2 sites)
[] motif C cycle 17 AP -1607.6 (54 sites)
[] motif D cycle 17 AP -1800.3 (37 sites)
Total Map : 2075.76 Prev: 2074.61 Diff: 1.14881 Motifs: 122

18
19
20
21
22
23
24
25
26
27
28
29
30
31
32[] motif A cycle 32 AP -377.5 (25 sites)
[] motif B cycle 32 AP -107.2 (3 sites)
[] motif C cycle 32 AP -1607.6 (54 sites)
[] motif D cycle 32 AP -1800.3 (37 sites)
Total Map : 2090.08 Prev: 2075.76 Diff: 14.3209 Motifs: 119

33
34[] motif A cycle 34 AP -398.0 (26 sites)
[] motif B cycle 34 AP -107.2 (3 sites)
[] motif C cycle 34 AP -1607.6 (54 sites)
[] motif D cycle 34 AP -1800.3 (37 sites)
Total Map : 2090.33 Prev: 2090.08 Diff: 0.249073 Motifs: 120

35[] motif A cycle 35 AP -377.0 (25 sites)
[] motif B cycle 35 AP -89.9 (3 sites)
[] motif C cycle 35 AP -1607.6 (54 sites)
[] motif D cycle 35 AP -1800.3 (37 sites)
Total Map : 2095.44 Prev: 2090.33 Diff: 5.10984 Motifs: 119

36[] motif A cycle 36 AP -377.0 (25 sites)
[] motif B cycle 36 AP -159.5 (5 sites)
[] motif C cycle 36 AP -1607.6 (54 sites)
[] motif D cycle 36 AP -1800.3 (37 sites)
Total Map : 2099.24 Prev: 2095.44 Diff: 3.79902 Motifs: 121

37[] motif A cycle 37 AP -377.0 (25 sites)
[] motif B cycle 37 AP -160.6 (5 sites)
[] motif C cycle 37 AP -1607.6 (54 sites)
[] motif D cycle 37 AP -1800.3 (37 sites)
Total Map : 2102.31 Prev: 2099.24 Diff: 3.07348 Motifs: 121

38
39
40
41[] motif A cycle 41 AP -356.4 (24 sites)
[] motif B cycle 41 AP -160.7 (5 sites)
[] motif C cycle 41 AP -1607.6 (54 sites)
[] motif D cycle 41 AP -1800.3 (37 sites)
Total Map : 2103.18 Prev: 2102.31 Diff: 0.864532 Motifs: 120

42
43[] motif A cycle 43 AP -377.0 (25 sites)
[] motif B cycle 43 AP -160.7 (5 sites)
[] motif C cycle 43 AP -1607.6 (54 sites)
[] motif D cycle 43 AP -1800.3 (37 sites)
Total Map : 2103.36 Prev: 2103.18 Diff: 0.182404 Motifs: 121

44
45[] motif A cycle 45 AP -422.3 (27 sites)
[] motif B cycle 45 AP -229.0 (7 sites)
[] motif C cycle 45 AP -1607.6 (54 sites)
[] motif D cycle 45 AP -1800.3 (37 sites)
Total Map : 2104.94 Prev: 2103.36 Diff: 1.57855 Motifs: 125

46[] motif A cycle 46 AP -422.3 (27 sites)
[] motif B cycle 46 AP -187.8 (6 sites)
[] motif C cycle 46 AP -1607.6 (54 sites)
[] motif D cycle 46 AP -1668.3 (35 sites)
Total Map : 2109.23 Prev: 2104.94 Diff: 4.29029 Motifs: 122

47[] motif A cycle 47 AP -446.7 (28 sites)
[] motif B cycle 47 AP -187.8 (6 sites)
[] motif C cycle 47 AP -1607.6 (54 sites)
[] motif D cycle 47 AP -1668.3 (35 sites)
Total Map : 2109.87 Prev: 2109.23 Diff: 0.644012 Motifs: 123

48
49[] motif A cycle 49 AP -490.2 (30 sites)
[] motif B cycle 49 AP -187.8 (6 sites)
[] motif C cycle 49 AP -1607.6 (54 sites)
[] motif D cycle 49 AP -1668.3 (35 sites)
Total Map : 2110.99 Prev: 2109.87 Diff: 1.11211 Motifs: 125

50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
MAX :: 2110.985589 (Seed = 1149743202, Iteration = 49   Motif A = 30 Motif B = 6 Motif C = 54 Motif D = 35 )
Max subopt MAP found on seed 10





======================== NEAR OPTIMAL RESULTS ========================
======================================================================
MAP = 584 maybe = 588 discard = 64640
Max set 2111.157969 at 4

5
10
15
20
25
30
35
40
45
50
55
60
65
70
75
80
85
90
95
100
105
110
115
120
125
130
135
140
145
150
155
160
165
170
175
180
185
190
195
200
205
210
215
220
225
230
235
240
245
250
255
260
265
270
275
280
285
290
295
300
305
310
315
320
325
330
335
340
345
350
355
360
365
370
375
380
385
390
395
400
405
410
415
420
425
430
435
440
445
450
455
460
465
470
475
480
485
490
495
500

=============================================================
======                Results by Sequence               =====
=============================================================


   1,  1,  3      37 ldfdk EFRDKTVVIVAIPGAFTPTCTANHIPPF vekft      64   1.00 F 1091044
   1,  2,  0      79 agvda VIVLSANDPFVQ              safgk      90   1.00 F 1091044
   2,  1,  3      32 irlsd YRGKKYVILFFYPANFTAISPTELMLLS drise      59   1.00 F 11467494
   2,  2,  0      72 klstq ILAISVDSPFSH              lqyll      83   1.00 F 11467494
   2,  3,  1     161 riles IQYVKENPGYACPVNWNFG       dqvfy     179   1.00 F 11467494
   3,  1,  2      17 viqsd KLVVVDFYADWCMPCRYISPILEKL skeyn      41   1.00 F 11499727
   4,  1,  2      22 llntt QYVVADFYADWCGPCKAIAPMYAQF aktfs      46   1.00 F 1174686
   5,  1,  2      19 ifsak KNVIVDFWAAWCGPCKLTSPEFQKA adefs      43   1.00 F 12044976
   6,  1,  3      26 eikei DLKSNWNVFFFYPYSYSFICPLELKNIS nkike      53   1.00 F 13186328
   6,  2,  0      66 nlntk IYAISNDSHFVQ              knwie      77   1.00 F 13186328
   7,  1,  2      21 kenhs KPILIDFYADWCPPCRMLIPVLDSI ekkhg      45   1.00 F 13358154
   8,  1,  3      28 kirls SYRGKWVVLFFYPADFTFVCPTEVEGFA edyek      55   1.00 F 13541053
   8,  2,  0      68 kknte VISVSEDTVYVH              kawvq      79   1.00 F 13541053
   9,  1,  3      26 mrkls EFRGQNVVLAFFPGAFTSVCTKEMCTFR dsman      53   1.00 F 13541117
   9,  2,  0      66 kfkak VIGISVDSPFSL              aefak      77   1.00 F 13541117
  10,  1,  2      17 aetse GVVLADFWAPWCGPCKMIAPVLEEL dqemg      41   1.00 F 135765
  11,  1,  2      29 akesn KLIVIDFTASWCPPCRMIAPIFNDL akkfm      53   1.00 F 1388082
  12,  1,  2      44 likqn DKLVIDFYATWCGPCKMMQPHLTKL iqayp      68   1.00 F 140543
  13,  1,  3      25 melpd EFEGKWFILFSHPADFTPVCTTEFVAFQ evype      52   1.00 F 14286173
  13,  2,  0      65 eldce LVGLSVDQVFSH              ikwie      76   1.00 F 14286173
  14,  1,  2      80 selrg KVVMLQFTASWCGVCRKEMPFIEKD iwlkh     104   1.00 F 14578634
  15,  1,  3      25 kirls DFRGRIVVLYFYPRAMTPGCTREGVRFN ellde      52   1.00 F 14600438
  15,  2,  0      65 klgav VIGVSTDSVEKN              rkfae      76   1.00 F 14600438
  16,  1,  2      23 lqnsd KPVLVDFYATWCGPCQLMVPILNEV setlk      47   1.00 F 15218394
  17,  1,  2     157 adfrg RPLVINLWASWCPPCRREMPVLQQA qaenp     181   1.00 F 15597673
  18,  1,  2      26 ensfh KPVLVDFWADWCAPCKALMPLLAQI aesyq      50   1.00 F 15599256
  19,  1,  2      67 negkg KTILLNFWSETCGVCIAELKTFEQL lqsyp      91   1.00 F 15602312
  20,  1,  2      61 eefkg KVLLINFWATWCPPCKEEIPMFKEI yekyr      85   1.00 F 15605725
  21,  1,  0      80 megvd VTVVSMDLPFAQ              krfce      91   1.00 F 15605963
  22,  1,  3      26 vtlrg YRGAKNVLLVFFPLAFTGICQGELDQLR dhlpe      53   1.00 F 15609375
  23,  1,  3      30 nvsla DYRGRRVIVYFYPAASTPGCTKQACDFR dnlgd      57   1.00 F 15609658
  23,  2,  0      70 tagln VVGISPDKPEKL              atfrd      81   1.00 F 15609658
  24,  1,  3      24 tvsls DFKGKNIVLYFYPKDMTPGCTTEACDFR drved      51   1.00 F 15613511
  24,  2,  0      64 glntv ILGVSPDPVERH              kkfie      75   1.00 F 15613511
  25,  1,  2      60 sdyrg DVVILNVWASWCEPCRKEMPALMEL qsdye      84   1.00 F 15614085
  26,  1,  2      63 releg KGVFLNFWGTYCPPCEREMPHMEKL ygeyk      87   1.00 F 15614140
  27,  1,  2      72 sslrg QPVILHFFATWCPVCQDEMPSLVKL dkeyr      96   1.00 F 15615431
  28,  1,  3      20 tfthv DLYGKYTILFFFPKAGTSGCTREAVEFS renfe      47   1.00 F 15643152
  28,  2,  0      56 fekaq VVGISRDSVEAL              krfke      67   1.00 F 15643152
  30,  1,  3      61 gltda LADNRAVVLFFYPFDFSPVCATELCAIQ narwf      88   1.00 F 15790738
  30,  2,  0     101 tpgla VWGISPDSTYAH              eafad     112   1.00 F 15790738
  31,  1,  2       2     m TVTLKDFYADWCGPCKTQDPILEEL eadyd      26   1.00 F 15791337
  32,  1,  0      80 idntv VLCISADLPFAQ              srfcg      91   1.00 F 15801846
  33,  1,  2      72 adyrg RPVVLNFWASWCGPCREEAPLFAKL aahpg      96   1.00 F 15805225
  34,  1,  2      78 taaqg KPVVINFWASWCVPCRQEAPLFSKL sqeta     102   1.00 F 15805374
  35,  1,  3      26 itlss YRGQSHVVLVFYPLDFSPVCSMQLPEYS gsqdd      53   1.00 F 15807234
  35,  2,  0      66 eagav VLGINRDSVYAH              rawaa      77   1.00 F 15807234
  36,  1,  3      28 vnlae LFKGKKGVLFGVPGAFTPGCSKTHLPGF veqae      55   1.00 F 15826629
  37,  1,  2      49 fitkn KIVVVDFWAEWCAPCLILAPVIEEL andyp      73   1.00 F 15899007
  38,  1,  3      26 vkips DFKGKVVVLAFYPAAFTSVCTKEMCTFR dsmak      53   1.00 F 15899339
  38,  2,  0      66 evnav VIGISVDPPFSN              kafke      77   1.00 F 15899339
  39,  1,  3      30 vttel LFKGKRVVLFAVPGAFTPTCSLNHLPGY lenrd      57   1.00 F 15964668
  40,  1,  2      61 easrq QPVLVDFWAPWCGPCKQLTPVIEKV vreaa      85   1.00 F 15966937
  41,  1,  2      61 sdfrg KTLLVNLWATWCVPCRKEMPALDEL qgkls      85   1.00 F 15988313
  42,  1,  2      60 qdakg KKVLLNFWATWCKPCRQEMPAMEKL qkeya      84   1.00 F 16078864
  43,  1,  2      53 llqdd LPMVIDFWAPWCGPCRSFAPIFAET aaera      77   1.00 F 16123427
  44,  1,  3      50 fnlak ALKKGPVVLYFFPAAYTAGCTAEAREFA eatpe      77   1.00 F 16125919
  46,  1,  2      21 vlkad GAILVDFWAEWCGPCKMIAPILDEI adeyq      45   1.00 F 1633495
  47,  1,  3      31 fnfkq HTNGKTTVLFFWPMDFTFVCPSELIAFD kryee      58   1.00 F 16501671
  47,  2,  0      71 krgve VVGVSFDSEFVH              nawrn      82   1.00 F 16501671
  47,  3,  1     160 lrmvd ALQFHEEHGDVCPAQWEKG       kegmn     178   1.00 F 16501671
  48,  1,  2      34 vlqcp KPILVYFGAPWCGLCHFVKPLLNHL hgewq      58   1.00 F 1651717
  49,  1,  2      60 tlsee RPVLLYFWASWCGVCRFTTPAVAHL aaege      84   1.00 F 16759994
  50,  1,  2      53 llkdd LPVVIDFWAPWCGPCRNFAPIFEDV aeers      77   1.00 F 16761507
  51,  1,  3      33 slekn IEDDKWTILFFYPMDFTFVCPTEIVAIS arsde      60   1.00 F 16803644
  52,  1,  2      19 iissh PKILLNFWAEWCAPCRCFWPTLEQF aemee      43   1.00 F 16804867
  53,  1,  3      31 vttdd LFAGKTVAVFSLPGAFTPTCSSTHLPGY nelak      58   1.00 F 17229033
  53,  2,  0      73 ngvde IVCISVNDAFVM              newak      84   1.00 F 17229033
  54,  1,  2      22 vlsed KVVVVDFTATWCGPCRLVSPLMDQL adeyk      46   1.00 F 17229859
  55,  1,  2      18 vlegt GYVLVDYFSDGCVPCKALMPAVEEL skkye      42   1.00 F 1729944
  56,  1,  2      28 rqhpe KIIILDFYATWCGPCKAIAPLYKEL atthk      52   1.00 F 17531233
  57,  1,  2      27 ehlkg KIIGLYFSASWCPPCRAFTPKLKEF feeik      51   1.00 F 17537401
  58,  1,  2      63 safrg QPVVINFWAPWCGPCVEEMPELSAL aqeqk      87   1.00 F 17547503
  59,  1,  2     286 seykg KTIFLNFWATWCPPCRGEMPYIDEL ykeyn     310   1.00 F 18309723
  60,  1,  3      28 rlsev LKRGRPVVLLFFPGAFTSVCTKELCTFR dkmal      55   1.00 F 18313548
  60,  2,  0      68 kanae VLAISVDSPFAL              kafkd      79   1.00 F 18313548
  61,  1,  2      44 dsllg KKIGLYFSAAWCGPCQRFTPQLVEV ynels      68   1.00 F 18406743
  61,  2,  2     364 sdlvg KTILMYFSAHWCPPCRAFTPKLVEV ykqik     388   1.00 F 18406743
  62,  1,  3      26 eislq DYIGKYVVLAFYPLDFTFVCPTEINRFS dlkga      53   1.00 F 19173077
  62,  2,  0      66 rrnav VLLISCDSVYTH              kawas      77   1.00 F 19173077
  63,  1,  2      15 sdfeg EVVVLNAWGQWCAPCRAEVDDLQLV qetld      39   1.00 F 19554157
  64,  1,  2      39 eeykg KVVVINFWATWCGYCVEEMPGFEKV ykefg      63   1.00 F 19705357
  66,  1,  2       7 agdfm KPMLLDFSATWCGPCRMQKPILEEL ekkyg      31   1.00 F 20092028
  67,  1,  3      27 evtek DTEGRWSVFFFYPADFTFVCPTELGDVA dhyee      54   1.00 F 20151112
  67,  2,  0      67 klgvd VYSVSTDTHFTH              kawhs      78   1.00 F 20151112
  67,  3,  1     154 rkika AQYVAAHPGEVCPAKWKEG       eatla     172   1.00 F 20151112
  68,  1,  3      29 vdtht LFTGRKVVLFAVPGAFTPTCSAKHLPGY veqfe      56   1.00 F 21112072
  69,  1,  2     103 adykg KVVVLNVWGSWCPPCRAEAKNFEKV yqdvk     127   1.00 F 21222859
  70,  1,  3      32 qinhk TYEGQWKVVFAWPKDFTFVCPTEIAAFG klnde      59   1.00 F 21223405
  70,  2,  0      72 drdaq ILGFSGDSEFVH              hawrk      83   1.00 F 21223405
  71,  1,  3      28 eihly DLKGKKVLLSFHPLAWTQVCAQQMKSLE enyel      55   1.00 F 21227878
  72,  1,  0      78 keegi VLTISADLPFAQ              krwca      89   1.00 F 21283385
  73,  1,  3      25 mvsls EFKGRKVLLIFYPGDDTPVCTAQLCDYR nnvaa      52   1.00 F 21674812
  73,  2,  0      65 srgit VIGISGDSPESH              kqfae      76   1.00 F 21674812
  74,  1,  2      53 sdfkg ERVLINFWTTWCPPCRQEMPDMQRF yqdlq      77   1.00 F 23098307
  76,  1,  2      20 kylqh QRVVVDFSAEWCGPCRAIAPVFDKL sneft      44   1.00 F 267116
  77,  1,  2      81 aafkg KVSLVNVWASWCVPCHDEAPLLTEL gkdkr     105   1.00 F 27375582
  78,  1,  2      34 vtsdn DVVLADFYADWCGPCQMLEPVVETL aeqtd      58   1.00 F 2822332
  79,  1,  2      77 sdlkg KKVILNFWATWCGPCQQEMPDMEAF ykehk     101   1.00 F 30021713
  80,  1,  2      19 tisan SNVLVYFWAPLCAPCDLFTPTYEAS srkhf      43   1.00 F 3261501
  81,  1,  3      28 irfhd FLGDSWGILFSHPRDFTPVCTTELGRAA klape      55   1.00 F 3318841
  81,  2,  0      68 krnvk LIALSIDSVEDH              lawsk      79   1.00 F 3318841
  81,  3,  1     166 lrvvi SLQLTAEKRVATPVDWKDG       dsvmv     184   1.00 F 3318841
  82,  1,  2      19 tietn PLVIVDFWAPWCGSCKMLGPVLEEV esevg      43   1.00 F 3323237
  83,  1,  2      17 ektah QAVVVNVGASWCPDCRKIEPIMENL aktyk      41   1.00 F 4155972
  84,  1,  2      79 vvnse TPVVVDFHAQWCGPCKILGPRLEKM vakqh     103   1.00 F 4200327
  85,  1,  3      10 eidin EYKGKYVVLLFYPLDWTFVCPTEMIGYS evagq      37   1.00 F 4433065
  85,  2,  0      50 eince VIGVSVDSVYCH              qawce      61   1.00 F 4433065
  86,  1,  3      32 vsvhs IAAGKKVILFGVPGAFTPTCSMSHVPGF igkae      59   1.00 F 4704732
  86,  2,  0      74 kgide IICFSVNDPFVM              kawgk      85   1.00 F 4704732
  87,  1,  3      28 fdfyk YVGDNWAILFSHPHDFTPVCTTELAEFG kmhee      55   1.00 F 4996210
  87,  2,  0      68 klnck LIGFSCNSKESH              dqwie      79   1.00 F 4996210
  87,  3,  1     163 lrvlk SLQLTNTHPVATPVNWKEG       dkcci     181   1.00 F 4996210
  88,  1,  3      41 ynask EFANKKVVLFALPGAFTPVCSANHVPEY iqklp      68   1.00 F 5326864
  89,  1,  3      88 slkki TENNRVVVFFVYPRASTPGCTRQACGFR dnyqe     115   1.00 F 6322180
  90,  1,  3      43 ewskl ISENKKVIITGAPAAFSPTCTVSHIPGY inyld      70   1.00 F 6323138
  91,  1,  2      20 nenkg RLIVVDFFAQWCGPCRNIAPKVEAL akeip      44   1.00 F 6687568
  92,  1,  0      68 klgve VLSVSVDSVFVH              kmwnd      79   1.00 F 6850955
  93,  1,  2      18 llttn KKVVVDFYANWCGPCKILGPIFEEV aqdkk      42   1.00 F 7109697
  94,  1,  2      21 ilaed KLVVIDFYADWCGPCKIIAPKLDEL aqqys      45   1.00 F 7290567
  95,  1,  3      31 evkls DYKGKYVVLFFYPLDFTFVCPTEIIAFS nraed      58   1.00 F 9955016
  95,  2,  0      71 klgce VLGVSVDSQFTH              lawin      82   1.00 F 9955016
  95,  3,  1     160 lrlvq AFQYTDEHGEVCPAGWKPG       sdtik     178   1.00 F 9955016
  96,  1,  2      49 adlqg KVTLINFWFPSCPGCVSEMPKIIKT andyk      73   1.00 F 15677788
124 motifs




Column 1 :  Sequence Number, Site Number
Column 2 :  Motif type
Column 3 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

-------------------------------------------------------------------------
                          MOTIF a





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    .   68  .   .   .   .   .   .   20  .   .   10  .   .   .   .   .   .   .   .   2.4
   2 |    .   17  .   .   .   .   .   .   34   3  .   34  .   .    6  .   .   .   .    3  1.7
   3 |    13   6  10  .   .   .   51  .   .   .   .    3  .   .   .   .   .   .   10   3  1.9
   4 |    .   31  .   .   .   10  .   .   48  .   .   10  .   .   .   .   .   .   .   .   2.1
   5 |    .   .   .   .   .   .   .   .   .   .   .   .   .    3  .   .   .   .   96  .   3.8

   7 |    .   .   .   86  .   .   .   .   .   .   .   .   .   13  .   .   .   .   .   .   3.2
   8 |    .   .   .   10  .   .   .   .   .   .    3  10  .   .   .    6   3  .   58   6  2.0
   9 |     3  34  .   .    6  .   .    6  .   .    3  .   .   .   .   37   3  .   .    3  1.8
  10 |    .   .   .   .   24  58  .   .   .   .   .   .   .   .   17  .   .   .   .   .   2.9

  12 |    .   .   .   .   .   .   .   55  .   .   .   13   6   6  .   .   17  .   .   .   3.4

nonsite    8   8  .    6   7   4   7   1   6  .    7   9   2   4   2   4   3   4   5   5
site       1  15   1   9   3   6   5   6  10  .   .    8  .    2   2   4   2  .   16   1

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.001 0.679 0.000 0.001 0.001 0.001 0.001 0.000 0.204 0.000 0.001 0.103 0.000 0.001 0.000 0.001 0.001 0.001 0.001 0.001 
   2 |  0.001 0.171 0.000 0.001 0.001 0.001 0.001 0.000 0.340 0.034 0.001 0.341 0.000 0.001 0.068 0.001 0.001 0.001 0.001 0.035 
   3 |  0.137 0.069 0.102 0.001 0.001 0.001 0.510 0.000 0.001 0.000 0.001 0.035 0.000 0.001 0.000 0.001 0.001 0.001 0.103 0.035 
   4 |  0.001 0.306 0.000 0.001 0.001 0.103 0.001 0.000 0.476 0.000 0.001 0.103 0.000 0.001 0.000 0.001 0.001 0.001 0.001 0.001 
   5 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.035 0.000 0.001 0.001 0.001 0.950 0.001 

   7 |  0.001 0.001 0.000 0.849 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.136 0.000 0.001 0.001 0.001 0.001 0.001 
   8 |  0.001 0.001 0.000 0.103 0.001 0.001 0.001 0.000 0.001 0.000 0.035 0.103 0.000 0.001 0.000 0.069 0.034 0.001 0.577 0.069 
   9 |  0.035 0.340 0.000 0.001 0.069 0.001 0.001 0.068 0.001 0.000 0.035 0.002 0.000 0.001 0.000 0.374 0.034 0.001 0.001 0.035 
  10 |  0.001 0.001 0.000 0.001 0.238 0.577 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.001 0.170 0.001 0.001 0.001 0.001 0.001 

  12 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.543 0.001 0.000 0.001 0.137 0.068 0.068 0.000 0.001 0.170 0.001 0.001 0.001 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



10 columns
Num Motifs: 29
   1,  1      79 agvda VIVLSANDPFVQ safgk      90   1.00 F 1091044
   2,  1      72 klstq ILAISVDSPFSH lqyll      83   1.00 F 11467494
   6,  1      66 nlntk IYAISNDSHFVQ knwie      77   1.00 F 13186328
   8,  1      68 kknte VISVSEDTVYVH kawvq      79   1.00 F 13541053
   9,  1      66 kfkak VIGISVDSPFSL aefak      77   1.00 F 13541117
  13,  1      65 eldce LVGLSVDQVFSH ikwie      76   1.00 F 14286173
  15,  1      65 klgav VIGVSTDSVEKN rkfae      76   1.00 F 14600438
  21,  1      80 megvd VTVVSMDLPFAQ krfce      91   1.00 F 15605963
  23,  1      70 tagln VVGISPDKPEKL atfrd      81   1.00 F 15609658
  24,  1      64 glntv ILGVSPDPVERH kkfie      75   1.00 F 15613511
  28,  1      56 fekaq VVGISRDSVEAL krfke      67   1.00 F 15643152
  30,  1     101 tpgla VWGISPDSTYAH eafad     112   1.00 F 15790738
  32,  1      80 idntv VLCISADLPFAQ srfcg      91   1.00 F 15801846
  35,  1      66 eagav VLGINRDSVYAH rawaa      77   1.00 F 15807234
  38,  1      66 evnav VIGISVDPPFSN kafke      77   1.00 F 15899339
  47,  1      71 krgve VVGVSFDSEFVH nawrn      82   1.00 F 16501671
  53,  1      73 ngvde IVCISVNDAFVM newak      84   1.00 F 17229033
  60,  1      68 kanae VLAISVDSPFAL kafkd      79   1.00 F 18313548
  62,  1      66 rrnav VLLISCDSVYTH kawas      77   1.00 F 19173077
  67,  1      67 klgvd VYSVSTDTHFTH kawhs      78   1.00 F 20151112
  70,  1      72 drdaq ILGFSGDSEFVH hawrk      83   1.00 F 21223405
  72,  1      78 keegi VLTISADLPFAQ krwca      89   1.00 F 21283385
  73,  1      65 srgit VIGISGDSPESH kqfae      76   1.00 F 21674812
  81,  1      68 krnvk LIALSIDSVEDH lawsk      79   1.00 F 3318841
  85,  1      50 eince VIGVSVDSVYCH qawce      61   1.00 F 4433065
  86,  1      74 kgide IICFSVNDPFVM kawgk      85   1.00 F 4704732
  87,  1      68 klnck LIGFSCNSKESH dqwie      79   1.00 F 4996210
  92,  1      68 klgve VLSVSVDSVFVH kmwnd      79   1.00 F 6850955
  95,  1      71 klgce VLGVSVDSQFTH lawin      82   1.00 F 9955016
                       ***** **** *


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif a = -469.15170
Log Fragmentation portion of MAP for motif a = -3.80666



=============================================================
====== ELEMENTS OCCURRING GREATER THAN  50% OF THE TIME =====
======                    Motif a                       =====
=============================================================


Listing of those elements occurring greater than 50% of the time
in near optimal sampling using 500 iterations



Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    .   74  .   .   .   .   .   .   14  .   .   11  .   .   .   .   .   .   .   .   2.5
   2 |    .   14  .   .   .    3  .   .   29   3  .   37  .   .    7  .   .   .   .    3  1.6
   3 |    14   3   3  .   .   .   59  .   .   .   .    3  .   .   .   .   .   .   11   3  1.9
   4 |    .   33  .   .   .    7  .   .   48  .   .   11  .   .   .   .   .   .   .   .   2.1
   5 |    .   .   .   .   .   .   .   .   .   .   .   .   .    3  .   .   .   .   96  .   3.8

   7 |    .   .   .   96  .   .   .   .   .   .   .   .   .    3  .   .   .   .   .   .   3.5
   8 |    .   .   .   .   .   .   .   .   .   .    3  11  .   .   .    7   3  .   66   7  2.4
   9 |    .   40  .   .    7  .   .    7  .   .    3  .   .   .   .   33   3  .   .    3  1.9
  10 |    .   .   .   .   25  51  .   .   .   .   .   .   .   .   18  .   .   .   .    3  2.7

  12 |    .   .   .   .   .   .   .   59  .   .   .   14  .    7  .   .   18  .   .   .   3.6

nonsite    8   8  .    6   7   4   7   1   6  .    7   9   2   4   2   4   3   4   5   5
site       1  16  .    9   3   6   5   6   9  .   .    8  .    1   2   4   2  .   17   2

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.002 0.729 0.000 0.001 0.001 0.001 0.001 0.000 0.147 0.000 0.001 0.111 0.000 0.001 0.000 0.001 0.001 0.001 0.001 0.001 
   2 |  0.002 0.147 0.000 0.001 0.001 0.037 0.001 0.000 0.292 0.037 0.001 0.365 0.000 0.001 0.073 0.001 0.001 0.001 0.001 0.037 
   3 |  0.147 0.038 0.037 0.001 0.001 0.001 0.583 0.000 0.001 0.000 0.001 0.038 0.000 0.001 0.000 0.001 0.001 0.001 0.110 0.037 
   4 |  0.002 0.329 0.000 0.001 0.001 0.074 0.001 0.000 0.474 0.000 0.001 0.111 0.000 0.001 0.000 0.001 0.001 0.001 0.001 0.001 
   5 |  0.002 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.037 0.000 0.001 0.001 0.001 0.946 0.001 

   7 |  0.002 0.001 0.000 0.947 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.037 0.000 0.001 0.001 0.001 0.001 0.001 
   8 |  0.002 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.038 0.111 0.000 0.001 0.000 0.074 0.037 0.001 0.655 0.074 
   9 |  0.002 0.401 0.000 0.001 0.074 0.001 0.001 0.073 0.001 0.000 0.038 0.002 0.000 0.001 0.000 0.328 0.037 0.001 0.001 0.037 
  10 |  0.002 0.001 0.000 0.001 0.256 0.510 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.001 0.182 0.001 0.001 0.001 0.001 0.037 

  12 |  0.002 0.001 0.000 0.001 0.001 0.001 0.001 0.582 0.001 0.000 0.001 0.147 0.000 0.073 0.000 0.001 0.182 0.001 0.001 0.001 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



10 columns
Num Motifs: 27
   2,  1      72 klstq ILAISVDSPFSH lqyll      83   1.00 F 11467494
   6,  1      66 nlntk IYAISNDSHFVQ knwie      77   1.00 F 13186328
   8,  1      68 kknte VISVSEDTVYVH kawvq      79   1.00 F 13541053
   9,  1      66 kfkak VIGISVDSPFSL aefak      77   1.00 F 13541117
  13,  1      65 eldce LVGLSVDQVFSH ikwie      76   0.98 F 14286173
  15,  1      65 klgav VIGVSTDSVEKN rkfae      76   1.00 F 14600438
  21,  1      80 megvd VTVVSMDLPFAQ krfce      91   0.65 F 15605963
  23,  1      70 tagln VVGISPDKPEKL atfrd      81   0.99 F 15609658
  24,  1      64 glntv ILGVSPDPVERH kkfie      75   1.00 F 15613511
  28,  1      56 fekaq VVGISRDSVEAL krfke      67   1.00 F 15643152
  30,  1     101 tpgla VWGISPDSTYAH eafad     112   0.98 F 15790738
  32,  1      80 idntv VLCISADLPFAQ srfcg      91   0.99 F 15801846
  35,  1      66 eagav VLGINRDSVYAH rawaa      77   1.00 F 15807234
  38,  1      66 evnav VIGISVDPPFSN kafke      77   1.00 F 15899339
  47,  1      71 krgve VVGVSFDSEFVH nawrn      82   1.00 F 16501671
  60,  1      68 kanae VLAISVDSPFAL kafkd      79   1.00 F 18313548
  62,  1      66 rrnav VLLISCDSVYTH kawas      77   1.00 F 19173077
  67,  1      67 klgvd VYSVSTDTHFTH kawhs      78   1.00 F 20151112
  70,  1      72 drdaq ILGFSGDSEFVH hawrk      83   1.00 F 21223405
  72,  1      78 keegi VLTISADLPFAQ krwca      89   0.99 F 21283385
  73,  1      65 srgit VIGISGDSPESH kqfae      76   1.00 F 21674812
  81,  1      68 krnvk LIALSIDSVEDH lawsk      79   1.00 F 3318841
  85,  1      50 eince VIGVSVDSVYCH qawce      61   1.00 F 4433065
  87,  1      68 klnck LIGFSCNSKESH dqwie      79   0.56 F 4996210
  89,  1     127 kkyaa VFGLSADSVTSQ kkfqs     138   0.51 F 6322180
  92,  1      68 klgve VLSVSVDSVFVH kmwnd      79   1.00 F 6850955
  95,  1      71 klgce VLGVSVDSQFTH lawin      82   1.00 F 9955016
                       ***** **** *

-------------------------------------------------------------------------
                          MOTIF b





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    50  .   .   .   .   .   .   .   16  .   .   .   .   .   .   .   .   .   33  .   1.9
   2 |    .   .   .   .   .   16  .   .   .   .   .   50  .   .   .   .   33  .   .   .   2.1
   3 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   33  .   66  .   .   .   3.4
   4 |    .   33  .   .   .   16  .   .   .   .   .   33  .   .   16  .   .   .   .   .   1.7

   6 |    33  .   .   16  33  .   .   .   .   .   .   .   .   16  .   .   .   .   .   .   1.5

   8 |    .   .   .   .   .   .   .   50  .   .   16  .   .   .   .   33  .   .   .   .   3.2
   9 |    .   .   .   .   .   .   66  .   .   .   .   .   .   .   .   16  .   16  .   .   2.3
  10 |    .   33  .   16  33  .   .   .   .   .   .   .   .   .   16  .   .   .   .   .   1.7
  11 |    50  50  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   2.1
  12 |    .   .   66  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   33  4.4
  13 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  100  .   .   .   .   3.8
  14 |    50  50  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   2.1

  16 |    .   .   .   .   .   .   .   .   .  100  .   .   .   .   .   .   .   .   .   .   5.8
  17 |    .   .   .   .   16  .   .   .   .   .   66  .   .   16  .   .   .   .   .   .   2.1

  19 |    .   .   .   .   .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   3.2

nonsite    8   7  .    6   7   4   7   1   6  .    7   9   2   4   2   4   3   4   5   5
site      12  11   4   2   5   2  11   3   1   6   5   5  .    2   4  10   6   1   2   2

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.468 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.158 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.312 0.004 
   2 |  0.006 0.006 0.001 0.005 0.005 0.158 0.005 0.001 0.004 0.001 0.006 0.469 0.002 0.003 0.002 0.004 0.310 0.003 0.004 0.004 
   3 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.310 0.004 0.618 0.003 0.004 0.004 
   4 |  0.006 0.314 0.001 0.005 0.005 0.158 0.005 0.001 0.004 0.001 0.006 0.315 0.002 0.003 0.156 0.004 0.002 0.003 0.004 0.004 

   6 |  0.314 0.006 0.001 0.159 0.313 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.157 0.002 0.004 0.002 0.003 0.004 0.004 

   8 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.463 0.004 0.001 0.159 0.007 0.002 0.003 0.002 0.312 0.002 0.003 0.004 0.004 
   9 |  0.006 0.006 0.001 0.005 0.005 0.004 0.621 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.158 0.002 0.157 0.004 0.004 
  10 |  0.006 0.314 0.001 0.159 0.313 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.156 0.004 0.002 0.003 0.004 0.004 
  11 |  0.468 0.468 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 
  12 |  0.006 0.006 0.617 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.312 
  13 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.927 0.002 0.003 0.004 0.004 
  14 |  0.468 0.468 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 

  16 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.924 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 
  17 |  0.006 0.006 0.001 0.005 0.159 0.004 0.005 0.001 0.004 0.001 0.621 0.007 0.002 0.157 0.002 0.004 0.002 0.003 0.004 0.004 

  19 |  0.006 0.006 0.001 0.005 0.005 0.004 0.928 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



15 columns
Num Motifs: 6
   2,  1     161 riles IQYVKENPGYACPVNWNFG dqvfy     179   1.00 F 11467494
  47,  1     160 lrmvd ALQFHEEHGDVCPAQWEKG kegmn     178   1.00 F 16501671
  67,  1     154 rkika AQYVAAHPGEVCPAKWKEG eatla     172   1.00 F 20151112
  81,  1     166 lrvvi SLQLTAEKRVATPVDWKDG dsvmv     184   1.00 F 3318841
  87,  1     163 lrvlk SLQLTNTHPVATPVNWKEG dkcci     181   1.00 F 4996210
  95,  1     160 lrlvq AFQYTDEHGEVCPAGWKPG sdtik     178   1.00 F 9955016
                       **** * ******* ** *


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif b = -187.76179
Log Fragmentation portion of MAP for motif b = -7.77486



=============================================================
====== ELEMENTS OCCURRING GREATER THAN  50% OF THE TIME =====
======                    Motif b                       =====
=============================================================


Listing of those elements occurring greater than 50% of the time
in near optimal sampling using 500 iterations



Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    50  .   .   .   .   .   .   .   16  .   .   .   .   .   .   .   .   .   33  .   1.9
   2 |    .   .   .   .   .   16  .   .   .   .   .   50  .   .   .   .   33  .   .   .   2.1
   3 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   33  .   66  .   .   .   3.4
   4 |    .   33  .   .   .   16  .   .   .   .   .   33  .   .   16  .   .   .   .   .   1.7

   6 |    33  .   .   16  33  .   .   .   .   .   .   .   .   16  .   .   .   .   .   .   1.5

   8 |    .   .   .   .   .   .   .   50  .   .   16  .   .   .   .   33  .   .   .   .   3.2
   9 |    .   .   .   .   .   .   66  .   .   .   .   .   .   .   .   16  .   16  .   .   2.3
  10 |    .   33  .   16  33  .   .   .   .   .   .   .   .   .   16  .   .   .   .   .   1.7
  11 |    50  50  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   2.1
  12 |    .   .   66  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   33  4.4
  13 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  100  .   .   .   .   3.8
  14 |    50  50  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   2.1

  16 |    .   .   .   .   .   .   .   .   .  100  .   .   .   .   .   .   .   .   .   .   5.8
  17 |    .   .   .   .   16  .   .   .   .   .   66  .   .   16  .   .   .   .   .   .   2.1

  19 |    .   .   .   .   .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   3.2

nonsite    8   7  .    6   7   4   7   1   6  .    7   9   2   4   2   4   3   4   5   5
site      12  11   4   2   5   2  11   3   1   6   5   5  .    2   4  10   6   1   2   2

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.468 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.158 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.312 0.004 
   2 |  0.006 0.006 0.001 0.005 0.005 0.158 0.005 0.001 0.004 0.001 0.006 0.469 0.002 0.003 0.002 0.004 0.310 0.003 0.004 0.004 
   3 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.310 0.004 0.618 0.003 0.004 0.004 
   4 |  0.006 0.314 0.001 0.005 0.005 0.158 0.005 0.001 0.004 0.001 0.006 0.315 0.002 0.003 0.156 0.004 0.002 0.003 0.004 0.004 

   6 |  0.314 0.006 0.001 0.159 0.313 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.157 0.002 0.004 0.002 0.003 0.004 0.004 

   8 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.463 0.004 0.001 0.159 0.007 0.002 0.003 0.002 0.312 0.002 0.003 0.004 0.004 
   9 |  0.006 0.006 0.001 0.005 0.005 0.004 0.621 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.158 0.002 0.157 0.004 0.004 
  10 |  0.006 0.314 0.001 0.159 0.313 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.156 0.004 0.002 0.003 0.004 0.004 
  11 |  0.468 0.468 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 
  12 |  0.006 0.006 0.617 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.312 
  13 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.927 0.002 0.003 0.004 0.004 
  14 |  0.468 0.468 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 

  16 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.924 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 
  17 |  0.006 0.006 0.001 0.005 0.159 0.004 0.005 0.001 0.004 0.001 0.621 0.007 0.002 0.157 0.002 0.004 0.002 0.003 0.004 0.004 

  19 |  0.006 0.006 0.001 0.005 0.005 0.004 0.928 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



15 columns
Num Motifs: 6
   2,  1     161 riles IQYVKENPGYACPVNWNFG dqvfy     179   1.00 F 11467494
  47,  1     160 lrmvd ALQFHEEHGDVCPAQWEKG kegmn     178   1.00 F 16501671
  67,  1     154 rkika AQYVAAHPGEVCPAKWKEG eatla     172   1.00 F 20151112
  81,  1     166 lrvvi SLQLTAEKRVATPVDWKDG dsvmv     184   1.00 F 3318841
  87,  1     163 lrvlk SLQLTNTHPVATPVNWKEG dkcci     181   1.00 F 4996210
  95,  1     160 lrlvq AFQYTDEHGEVCPAGWKPG sdtik     178   1.00 F 9955016
                       **** * ******* ** *

-------------------------------------------------------------------------
                          MOTIF c





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    .   .   .    5   3  .    5  .   .   .   53   3  .   .   .    3  11   7   1   3  1.6

   3 |    .   61  .   .   .   .   .   .   22  .   .    7   3  .   .   .   .   .    1   3  2.1
   4 |    .   38  .   .   .    3   3  .   11  .   .   40   1  .   .   .   .   .   .   .   1.7
   5 |     5  35  .   .   .   .   .   .   24  .    1  31   1  .   .   .   .   .   .   .   1.7
   6 |    .   .   .   48  .   .   .    1  .   .   .   .   .   37  11  .    1  .   .   .   2.7
   7 |     1   7  .   .   .   85  .   .   .   .   .    3  .   .    1  .   .   .   .   .   3.4
   8 |    .   .   .   .   .    5   3   1  .   55  .   .   .   .   18  .   .   .    9   5  3.5
   9 |    87  .   .   .   .    1   5  .   .   .   .   .   .   .   .   .   .   .    3   1  2.8
  10 |     3  .   .   14   9  .   .    1  .   .   .   .   .    1  .   16   5  .   20  25  1.5
  11 |    .   .   .   .   .   .    1  .   .   90  .    1  .   .    1  .   .   .    1   1  5.2
  12 |    .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   6.0
  13 |     9   7  .   .    1  .   53  .   .   .    1  .    1  .   .   24  .   .   .   .   2.0
  14 |    .    7  .    1  .   .    1  .   .   .   .    1  .   .    1  83  .   .    1  .   3.2
  15 |    .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   6.0
  16 |    .    5  .    1   1  .   .    3   1  .   27   1  .   .   .   .    9  46  .   .   2.1

  18 |    .    3  .   .   37  12  .   .   18  .   .   16   3  .   .   .    3  .   .    3  1.4

  20 |    .   .   .    1  .   .   .   .   .   .    3  .   .   .   .   94  .   .   .   .   3.9

  22 |    .    7  .   .   .   22  .   .    9  .   .   44  11  .    5  .   .   .   .   .   1.8

  24 |     7  .   .    3  37  .   .    3  .   .   27   1  .    1  .   .   11   1   1   1  1.4
  25 |     3  18  .    1  .    9  .   .    7  .   .   51   1  .   .   .   .   .    1   3  1.4

nonsite    8   7   1   6   7   4   6   1   5   1   7   9   2   4   2   4   3   4   4   4
site       5   9  10   3   4   7   3  .    4   7   5  10   1   2   2  11   2   2   2   2

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.001 0.001 0.000 0.056 0.037 0.000 0.056 0.000 0.001 0.000 0.533 0.038 0.000 0.000 0.000 0.037 0.110 0.074 0.019 0.037 

   3 |  0.001 0.606 0.000 0.001 0.001 0.000 0.001 0.000 0.221 0.000 0.001 0.074 0.037 0.000 0.000 0.000 0.000 0.000 0.019 0.037 
   4 |  0.001 0.386 0.000 0.001 0.001 0.037 0.037 0.000 0.111 0.000 0.001 0.405 0.019 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
   5 |  0.056 0.349 0.000 0.001 0.001 0.000 0.001 0.000 0.239 0.000 0.019 0.313 0.019 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
   6 |  0.001 0.001 0.000 0.478 0.001 0.000 0.001 0.018 0.001 0.000 0.001 0.001 0.000 0.367 0.110 0.000 0.019 0.000 0.000 0.000 
   7 |  0.019 0.074 0.000 0.001 0.001 0.845 0.001 0.000 0.001 0.000 0.001 0.038 0.000 0.000 0.019 0.000 0.000 0.000 0.000 0.000 
   8 |  0.001 0.001 0.000 0.001 0.001 0.056 0.037 0.018 0.001 0.551 0.001 0.001 0.000 0.000 0.184 0.000 0.000 0.000 0.092 0.056 
   9 |  0.863 0.001 0.000 0.001 0.001 0.019 0.056 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.037 0.019 
  10 |  0.037 0.001 0.000 0.147 0.092 0.000 0.001 0.018 0.001 0.000 0.001 0.001 0.000 0.019 0.000 0.166 0.055 0.000 0.202 0.257 
  11 |  0.001 0.001 0.000 0.001 0.001 0.000 0.019 0.000 0.001 0.899 0.001 0.019 0.000 0.000 0.019 0.000 0.000 0.000 0.019 0.019 
  12 |  0.001 0.001 0.991 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
  13 |  0.093 0.074 0.000 0.001 0.019 0.000 0.533 0.000 0.001 0.000 0.019 0.001 0.019 0.000 0.000 0.239 0.000 0.000 0.000 0.000 
  14 |  0.001 0.074 0.000 0.019 0.001 0.000 0.019 0.000 0.001 0.000 0.001 0.019 0.000 0.000 0.019 0.826 0.000 0.000 0.019 0.000 
  15 |  0.001 0.001 0.991 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
  16 |  0.001 0.056 0.000 0.019 0.019 0.000 0.001 0.037 0.019 0.000 0.276 0.019 0.000 0.000 0.000 0.000 0.092 0.459 0.000 0.000 

  18 |  0.001 0.037 0.000 0.001 0.368 0.129 0.001 0.000 0.184 0.000 0.001 0.166 0.037 0.000 0.000 0.000 0.037 0.000 0.000 0.037 

  20 |  0.001 0.001 0.000 0.019 0.001 0.000 0.001 0.000 0.001 0.000 0.037 0.001 0.000 0.000 0.000 0.936 0.000 0.000 0.000 0.000 

  22 |  0.001 0.074 0.000 0.001 0.001 0.221 0.001 0.000 0.092 0.000 0.001 0.441 0.110 0.000 0.055 0.000 0.000 0.000 0.000 0.000 

  24 |  0.074 0.001 0.000 0.037 0.368 0.000 0.001 0.037 0.001 0.000 0.276 0.019 0.000 0.019 0.000 0.000 0.110 0.019 0.019 0.019 
  25 |  0.037 0.184 0.000 0.019 0.001 0.092 0.001 0.000 0.074 0.000 0.001 0.515 0.019 0.000 0.000 0.000 0.000 0.000 0.019 0.037 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



20 columns
Num Motifs: 54
   3,  1      17 viqsd KLVVVDFYADWCMPCRYISPILEKL skeyn      41   1.00 F 11499727
   4,  1      22 llntt QYVVADFYADWCGPCKAIAPMYAQF aktfs      46   1.00 F 1174686
   5,  1      19 ifsak KNVIVDFWAAWCGPCKLTSPEFQKA adefs      43   1.00 F 12044976
   7,  1      21 kenhs KPILIDFYADWCPPCRMLIPVLDSI ekkhg      45   1.00 F 13358154
  10,  1      17 aetse GVVLADFWAPWCGPCKMIAPVLEEL dqemg      41   1.00 F 135765
  11,  1      29 akesn KLIVIDFTASWCPPCRMIAPIFNDL akkfm      53   1.00 F 1388082
  12,  1      44 likqn DKLVIDFYATWCGPCKMMQPHLTKL iqayp      68   1.00 F 140543
  14,  1      80 selrg KVVMLQFTASWCGVCRKEMPFIEKD iwlkh     104   1.00 F 14578634
  16,  1      23 lqnsd KPVLVDFYATWCGPCQLMVPILNEV setlk      47   1.00 F 15218394
  17,  1     157 adfrg RPLVINLWASWCPPCRREMPVLQQA qaenp     181   1.00 F 15597673
  18,  1      26 ensfh KPVLVDFWADWCAPCKALMPLLAQI aesyq      50   1.00 F 15599256
  19,  1      67 negkg KTILLNFWSETCGVCIAELKTFEQL lqsyp      91   1.00 F 15602312
  20,  1      61 eefkg KVLLINFWATWCPPCKEEIPMFKEI yekyr      85   1.00 F 15605725
  25,  1      60 sdyrg DVVILNVWASWCEPCRKEMPALMEL qsdye      84   1.00 F 15614085
  26,  1      63 releg KGVFLNFWGTYCPPCEREMPHMEKL ygeyk      87   1.00 F 15614140
  27,  1      72 sslrg QPVILHFFATWCPVCQDEMPSLVKL dkeyr      96   1.00 F 15615431
  31,  1       2     m TVTLKDFYADWCGPCKTQDPILEEL eadyd      26   1.00 F 15791337
  33,  1      72 adyrg RPVVLNFWASWCGPCREEAPLFAKL aahpg      96   1.00 F 15805225
  34,  1      78 taaqg KPVVINFWASWCVPCRQEAPLFSKL sqeta     102   1.00 F 15805374
  37,  1      49 fitkn KIVVVDFWAEWCAPCLILAPVIEEL andyp      73   1.00 F 15899007
  40,  1      61 easrq QPVLVDFWAPWCGPCKQLTPVIEKV vreaa      85   1.00 F 15966937
  41,  1      61 sdfrg KTLLVNLWATWCVPCRKEMPALDEL qgkls      85   1.00 F 15988313
  42,  1      60 qdakg KKVLLNFWATWCKPCRQEMPAMEKL qkeya      84   1.00 F 16078864
  43,  1      53 llqdd LPMVIDFWAPWCGPCRSFAPIFAET aaera      77   1.00 F 16123427
  46,  1      21 vlkad GAILVDFWAEWCGPCKMIAPILDEI adeyq      45   1.00 F 1633495
  48,  1      34 vlqcp KPILVYFGAPWCGLCHFVKPLLNHL hgewq      58   1.00 F 1651717
  49,  1      60 tlsee RPVLLYFWASWCGVCRFTTPAVAHL aaege      84   1.00 F 16759994
  50,  1      53 llkdd LPVVIDFWAPWCGPCRNFAPIFEDV aeers      77   1.00 F 16761507
  52,  1      19 iissh PKILLNFWAEWCAPCRCFWPTLEQF aemee      43   1.00 F 16804867
  54,  1      22 vlsed KVVVVDFTATWCGPCRLVSPLMDQL adeyk      46   1.00 F 17229859
  55,  1      18 vlegt GYVLVDYFSDGCVPCKALMPAVEEL skkye      42   1.00 F 1729944
  56,  1      28 rqhpe KIIILDFYATWCGPCKAIAPLYKEL atthk      52   1.00 F 17531233
  57,  1      27 ehlkg KIIGLYFSASWCPPCRAFTPKLKEF feeik      51   1.00 F 17537401
  58,  1      63 safrg QPVVINFWAPWCGPCVEEMPELSAL aqeqk      87   1.00 F 17547503
  59,  1     286 seykg KTIFLNFWATWCPPCRGEMPYIDEL ykeyn     310   1.00 F 18309723
  61,  1      44 dsllg KKIGLYFSAAWCGPCQRFTPQLVEV ynels      68   1.00 F 18406743
  61,  2     364 sdlvg KTILMYFSAHWCPPCRAFTPKLVEV ykqik     388   1.00 F 18406743
  63,  1      15 sdfeg EVVVLNAWGQWCAPCRAEVDDLQLV qetld      39   1.00 F 19554157
  64,  1      39 eeykg KVVVINFWATWCGYCVEEMPGFEKV ykefg      63   1.00 F 19705357
  66,  1       7 agdfm KPMLLDFSATWCGPCRMQKPILEEL ekkyg      31   1.00 F 20092028
  69,  1     103 adykg KVVVLNVWGSWCPPCRAEAKNFEKV yqdvk     127   1.00 F 21222859
  74,  1      53 sdfkg ERVLINFWTTWCPPCRQEMPDMQRF yqdlq      77   1.00 F 23098307
  76,  1      20 kylqh QRVVVDFSAEWCGPCRAIAPVFDKL sneft      44   1.00 F 267116
  77,  1      81 aafkg KVSLVNVWASWCVPCHDEAPLLTEL gkdkr     105   1.00 F 27375582
  78,  1      34 vtsdn DVVLADFYADWCGPCQMLEPVVETL aeqtd      58   1.00 F 2822332
  79,  1      77 sdlkg KKVILNFWATWCGPCQQEMPDMEAF ykehk     101   1.00 F 30021713
  80,  1      19 tisan SNVLVYFWAPLCAPCDLFTPTYEAS srkhf      43   1.00 F 3261501
  82,  1      19 tietn PLVIVDFWAPWCGSCKMLGPVLEEV esevg      43   1.00 F 3323237
  83,  1      17 ektah QAVVVNVGASWCPDCRKIEPIMENL aktyk      41   1.00 F 4155972
  84,  1      79 vvnse TPVVVDFHAQWCGPCKILGPRLEKM vakqh     103   1.00 F 4200327
  91,  1      20 nenkg RLIVVDFFAQWCGPCRNIAPKVEAL akeip      44   1.00 F 6687568
  93,  1      18 llttn KKVVVDFYANWCGPCKILGPIFEEV aqdkk      42   1.00 F 7109697
  94,  1      21 ilaed KLVVIDFYADWCGPCKIIAPKLDEL aqqys      45   1.00 F 7290567
  96,  1      49 adlqg KVTLINFWFPSCPGCVSEMPKIIKT andyk      73   1.00 F 15677788
                       * ************** * * * **


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif c = -1607.59351
Log Fragmentation portion of MAP for motif c = -10.42374



=============================================================
====== ELEMENTS OCCURRING GREATER THAN  50% OF THE TIME =====
======                    Motif c                       =====
=============================================================


Listing of those elements occurring greater than 50% of the time
in near optimal sampling using 500 iterations



Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    .   .   .    5   3  .    5  .   .   .   53   3  .   .   .    3  11   7   1   3  1.6

   3 |    .   61  .   .   .   .   .   .   22  .   .    7   3  .   .   .   .   .    1   3  2.1
   4 |    .   38  .   .   .    3   3  .   11  .   .   40   1  .   .   .   .   .   .   .   1.7
   5 |     5  35  .   .   .   .   .   .   24  .    1  31   1  .   .   .   .   .   .   .   1.7
   6 |    .   .   .   48  .   .   .    1  .   .   .   .   .   37  11  .    1  .   .   .   2.7
   7 |     1   7  .   .   .   85  .   .   .   .   .    3  .   .    1  .   .   .   .   .   3.4
   8 |    .   .   .   .   .    5   3   1  .   55  .   .   .   .   18  .   .   .    9   5  3.5
   9 |    87  .   .   .   .    1   5  .   .   .   .   .   .   .   .   .   .   .    3   1  2.8
  10 |     3  .   .   14   9  .   .    1  .   .   .   .   .    1  .   16   5  .   20  25  1.5
  11 |    .   .   .   .   .   .    1  .   .   90  .    1  .   .    1  .   .   .    1   1  5.2
  12 |    .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   6.0
  13 |     9   7  .   .    1  .   53  .   .   .    1  .    1  .   .   24  .   .   .   .   2.0
  14 |    .    7  .    1  .   .    1  .   .   .   .    1  .   .    1  83  .   .    1  .   3.2
  15 |    .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   6.0
  16 |    .    5  .    1   1  .   .    3   1  .   27   1  .   .   .   .    9  46  .   .   2.1

  18 |    .    3  .   .   37  12  .   .   18  .   .   16   3  .   .   .    3  .   .    3  1.4

  20 |    .   .   .    1  .   .   .   .   .   .    3  .   .   .   .   94  .   .   .   .   3.9

  22 |    .    7  .   .   .   22  .   .    9  .   .   44  11  .    5  .   .   .   .   .   1.8

  24 |     7  .   .    3  37  .   .    3  .   .   27   1  .    1  .   .   11   1   1   1  1.4
  25 |     3  18  .    1  .    9  .   .    7  .   .   51   1  .   .   .   .   .    1   3  1.4

nonsite    8   7   1   6   7   4   6   1   5   1   7   9   2   4   2   4   3   4   4   4
site       5   9  10   3   4   7   3  .    4   7   5  10   1   2   2  11   2   2   2   2

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.001 0.001 0.000 0.056 0.037 0.000 0.056 0.000 0.001 0.000 0.533 0.038 0.000 0.000 0.000 0.037 0.110 0.074 0.019 0.037 

   3 |  0.001 0.606 0.000 0.001 0.001 0.000 0.001 0.000 0.221 0.000 0.001 0.074 0.037 0.000 0.000 0.000 0.000 0.000 0.019 0.037 
   4 |  0.001 0.386 0.000 0.001 0.001 0.037 0.037 0.000 0.111 0.000 0.001 0.405 0.019 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
   5 |  0.056 0.349 0.000 0.001 0.001 0.000 0.001 0.000 0.239 0.000 0.019 0.313 0.019 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
   6 |  0.001 0.001 0.000 0.478 0.001 0.000 0.001 0.018 0.001 0.000 0.001 0.001 0.000 0.367 0.110 0.000 0.019 0.000 0.000 0.000 
   7 |  0.019 0.074 0.000 0.001 0.001 0.845 0.001 0.000 0.001 0.000 0.001 0.038 0.000 0.000 0.019 0.000 0.000 0.000 0.000 0.000 
   8 |  0.001 0.001 0.000 0.001 0.001 0.056 0.037 0.018 0.001 0.551 0.001 0.001 0.000 0.000 0.184 0.000 0.000 0.000 0.092 0.056 
   9 |  0.863 0.001 0.000 0.001 0.001 0.019 0.056 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.037 0.019 
  10 |  0.037 0.001 0.000 0.147 0.092 0.000 0.001 0.018 0.001 0.000 0.001 0.001 0.000 0.019 0.000 0.166 0.055 0.000 0.202 0.257 
  11 |  0.001 0.001 0.000 0.001 0.001 0.000 0.019 0.000 0.001 0.899 0.001 0.019 0.000 0.000 0.019 0.000 0.000 0.000 0.019 0.019 
  12 |  0.001 0.001 0.991 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
  13 |  0.093 0.074 0.000 0.001 0.019 0.000 0.533 0.000 0.001 0.000 0.019 0.001 0.019 0.000 0.000 0.239 0.000 0.000 0.000 0.000 
  14 |  0.001 0.074 0.000 0.019 0.001 0.000 0.019 0.000 0.001 0.000 0.001 0.019 0.000 0.000 0.019 0.826 0.000 0.000 0.019 0.000 
  15 |  0.001 0.001 0.991 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
  16 |  0.001 0.056 0.000 0.019 0.019 0.000 0.001 0.037 0.019 0.000 0.276 0.019 0.000 0.000 0.000 0.000 0.092 0.459 0.000 0.000 

  18 |  0.001 0.037 0.000 0.001 0.368 0.129 0.001 0.000 0.184 0.000 0.001 0.166 0.037 0.000 0.000 0.000 0.037 0.000 0.000 0.037 

  20 |  0.001 0.001 0.000 0.019 0.001 0.000 0.001 0.000 0.001 0.000 0.037 0.001 0.000 0.000 0.000 0.936 0.000 0.000 0.000 0.000 

  22 |  0.001 0.074 0.000 0.001 0.001 0.221 0.001 0.000 0.092 0.000 0.001 0.441 0.110 0.000 0.055 0.000 0.000 0.000 0.000 0.000 

  24 |  0.074 0.001 0.000 0.037 0.368 0.000 0.001 0.037 0.001 0.000 0.276 0.019 0.000 0.019 0.000 0.000 0.110 0.019 0.019 0.019 
  25 |  0.037 0.184 0.000 0.019 0.001 0.092 0.001 0.000 0.074 0.000 0.001 0.515 0.019 0.000 0.000 0.000 0.000 0.000 0.019 0.037 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



20 columns
Num Motifs: 54
   3,  1      17 viqsd KLVVVDFYADWCMPCRYISPILEKL skeyn      41   1.00 F 11499727
   4,  1      22 llntt QYVVADFYADWCGPCKAIAPMYAQF aktfs      46   1.00 F 1174686
   5,  1      19 ifsak KNVIVDFWAAWCGPCKLTSPEFQKA adefs      43   1.00 F 12044976
   7,  1      21 kenhs KPILIDFYADWCPPCRMLIPVLDSI ekkhg      45   1.00 F 13358154
  10,  1      17 aetse GVVLADFWAPWCGPCKMIAPVLEEL dqemg      41   1.00 F 135765
  11,  1      29 akesn KLIVIDFTASWCPPCRMIAPIFNDL akkfm      53   1.00 F 1388082
  12,  1      44 likqn DKLVIDFYATWCGPCKMMQPHLTKL iqayp      68   1.00 F 140543
  14,  1      80 selrg KVVMLQFTASWCGVCRKEMPFIEKD iwlkh     104   1.00 F 14578634
  16,  1      23 lqnsd KPVLVDFYATWCGPCQLMVPILNEV setlk      47   1.00 F 15218394
  17,  1     157 adfrg RPLVINLWASWCPPCRREMPVLQQA qaenp     181   1.00 F 15597673
  18,  1      26 ensfh KPVLVDFWADWCAPCKALMPLLAQI aesyq      50   1.00 F 15599256
  19,  1      67 negkg KTILLNFWSETCGVCIAELKTFEQL lqsyp      91   1.00 F 15602312
  20,  1      61 eefkg KVLLINFWATWCPPCKEEIPMFKEI yekyr      85   1.00 F 15605725
  25,  1      60 sdyrg DVVILNVWASWCEPCRKEMPALMEL qsdye      84   1.00 F 15614085
  26,  1      63 releg KGVFLNFWGTYCPPCEREMPHMEKL ygeyk      87   1.00 F 15614140
  27,  1      72 sslrg QPVILHFFATWCPVCQDEMPSLVKL dkeyr      96   1.00 F 15615431
  31,  1       2     m TVTLKDFYADWCGPCKTQDPILEEL eadyd      26   1.00 F 15791337
  33,  1      72 adyrg RPVVLNFWASWCGPCREEAPLFAKL aahpg      96   1.00 F 15805225
  34,  1      78 taaqg KPVVINFWASWCVPCRQEAPLFSKL sqeta     102   1.00 F 15805374
  37,  1      49 fitkn KIVVVDFWAEWCAPCLILAPVIEEL andyp      73   1.00 F 15899007
  40,  1      61 easrq QPVLVDFWAPWCGPCKQLTPVIEKV vreaa      85   1.00 F 15966937
  41,  1      61 sdfrg KTLLVNLWATWCVPCRKEMPALDEL qgkls      85   1.00 F 15988313
  42,  1      60 qdakg KKVLLNFWATWCKPCRQEMPAMEKL qkeya      84   1.00 F 16078864
  43,  1      53 llqdd LPMVIDFWAPWCGPCRSFAPIFAET aaera      77   1.00 F 16123427
  46,  1      21 vlkad GAILVDFWAEWCGPCKMIAPILDEI adeyq      45   1.00 F 1633495
  48,  1      34 vlqcp KPILVYFGAPWCGLCHFVKPLLNHL hgewq      58   1.00 F 1651717
  49,  1      60 tlsee RPVLLYFWASWCGVCRFTTPAVAHL aaege      84   1.00 F 16759994
  50,  1      53 llkdd LPVVIDFWAPWCGPCRNFAPIFEDV aeers      77   1.00 F 16761507
  52,  1      19 iissh PKILLNFWAEWCAPCRCFWPTLEQF aemee      43   1.00 F 16804867
  54,  1      22 vlsed KVVVVDFTATWCGPCRLVSPLMDQL adeyk      46   1.00 F 17229859
  55,  1      18 vlegt GYVLVDYFSDGCVPCKALMPAVEEL skkye      42   1.00 F 1729944
  56,  1      28 rqhpe KIIILDFYATWCGPCKAIAPLYKEL atthk      52   1.00 F 17531233
  57,  1      27 ehlkg KIIGLYFSASWCPPCRAFTPKLKEF feeik      51   1.00 F 17537401
  58,  1      63 safrg QPVVINFWAPWCGPCVEEMPELSAL aqeqk      87   1.00 F 17547503
  59,  1     286 seykg KTIFLNFWATWCPPCRGEMPYIDEL ykeyn     310   1.00 F 18309723
  61,  1      44 dsllg KKIGLYFSAAWCGPCQRFTPQLVEV ynels      68   1.00 F 18406743
  61,  2     364 sdlvg KTILMYFSAHWCPPCRAFTPKLVEV ykqik     388   1.00 F 18406743
  63,  1      15 sdfeg EVVVLNAWGQWCAPCRAEVDDLQLV qetld      39   1.00 F 19554157
  64,  1      39 eeykg KVVVINFWATWCGYCVEEMPGFEKV ykefg      63   1.00 F 19705357
  66,  1       7 agdfm KPMLLDFSATWCGPCRMQKPILEEL ekkyg      31   1.00 F 20092028
  69,  1     103 adykg KVVVLNVWGSWCPPCRAEAKNFEKV yqdvk     127   1.00 F 21222859
  74,  1      53 sdfkg ERVLINFWTTWCPPCRQEMPDMQRF yqdlq      77   1.00 F 23098307
  76,  1      20 kylqh QRVVVDFSAEWCGPCRAIAPVFDKL sneft      44   1.00 F 267116
  77,  1      81 aafkg KVSLVNVWASWCVPCHDEAPLLTEL gkdkr     105   1.00 F 27375582
  78,  1      34 vtsdn DVVLADFYADWCGPCQMLEPVVETL aeqtd      58   1.00 F 2822332
  79,  1      77 sdlkg KKVILNFWATWCGPCQQEMPDMEAF ykehk     101   1.00 F 30021713
  80,  1      19 tisan SNVLVYFWAPLCAPCDLFTPTYEAS srkhf      43   1.00 F 3261501
  82,  1      19 tietn PLVIVDFWAPWCGSCKMLGPVLEEV esevg      43   1.00 F 3323237
  83,  1      17 ektah QAVVVNVGASWCPDCRKIEPIMENL aktyk      41   1.00 F 4155972
  84,  1      79 vvnse TPVVVDFHAQWCGPCKILGPRLEKM vakqh     103   1.00 F 4200327
  91,  1      20 nenkg RLIVVDFFAQWCGPCRNIAPKVEAL akeip      44   1.00 F 6687568
  93,  1      18 llttn KKVVVDFYANWCGPCKILGPIFEEV aqdkk      42   1.00 F 7109697
  94,  1      21 ilaed KLVVIDFYADWCGPCKIIAPKLDEL aqqys      45   1.00 F 7290567
  96,  1      49 adlqg KVTLINFWFPSCPGCVSEMPKIIKT andyk      73   1.00 F 15677788
                       * ************** * * * **

-------------------------------------------------------------------------
                          MOTIF d





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |     2  .   .   28  17   2  .    2   8  .   .   17  .   .   11  .   .   .    2   5  1.1
   2 |     5   2  .   .    5  34  .   .   .   .    2  14  .   .   17  .   .    8   2   5  1.4
   3 |     8  .   .    5  11  .   14  .    2  .   28  .   .    5   2  .   .   17  .    2  1.0
   4 |     2  .   .   11  .   .   62  .   .   .    5  .   .   11  .   .    2  .    2  .   2.1
   5 |    .   .   .   .   .   .    2  .   .   .   57  .   .    5  .   .    5  22   5  .   2.2

   7 |     2  68  .   .   .    2   5  .    2  .    2  .   .    2  .   .   .   .    2   8  1.9
   8 |     2  62  .   .   .   .   .   .   25  .   .    8  .   .   .   .   .   .   .   .   2.3
   9 |    .    8  .   .   .    8  .   .    5  .   .   77  .   .   .   .   .   .   .   .   2.3
  10 |     8   8  .   .   .   57  .   .    2  .   .    5  .   .   11  .   .   .    2   2  2.1
  11 |    14   2  .   .   .   62   8  .   .   .   .   .   .   .   .   .   .   .   11  .   2.4
  12 |     2  11  .   .   .   14  .   11   2   5  .    5  .   .   45  .   .   .   .   .   2.4
  13 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  100  .   .   .   .   4.3
  14 |    22  .   .   .   .    2  28   2  .   .    8  17   5  .    2  .   .    8  .   .   1.2
  15 |    51  .   .   42  .   .   .   .   .   .   .   .   .    2  .   .   .   .    2  .   2.3
  16 |    .   .   .    2  .   71   2  .   .    5  .   .    5  .    5  .   .   .    5  .   2.9
  17 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   11  88  3.6
  18 |     5  .   .   .   .   25   2  .   .   .   .   .   .   .   .   51   2  .   11  .   2.4
  19 |    .   54  .   .   .   .   20  .    8  .   .   .   .   .   .   .   .   .   .   17  2.0
  20 |    .   .   97  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .    2  .   6.2
  21 |     5  .   .   .   .   .   .   .   .   .   .   .   .   .   .   28   2  .   20  42  2.3
  22 |    14   2  .   .   .   .    2  .   .   .   14   5   5  .   .   .    2   8   5  37  1.3
  23 |    .   .   .   .   62  .   .   .   .   .    2  .   .    8  .   .   14  .    5   5  2.2
  24 |    14   2  .   .   .    2   2  22  11  .   .   31  11  .   .   .   .   .   .   .   1.8

  27 |     2   2  .   .    2  45  17  .    8  .   .    8  .   .    8   2  .   .   .   .   1.6
  28 |    11  .   .    2   2   8   5  .   .   .   .   .   .    2  14  .    5  22  22  .   1.4

nonsite    8   7  .    6   7   4   7   1   5  .    7   9   2   4   2   4   3   4   5   5
site       7   9   3   3   4  13   7   1   3  .    4   7   1   1   4   7   1   3   4   8

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.029 0.001 0.000 0.283 0.170 0.029 0.001 0.028 0.085 0.000 0.001 0.170 0.000 0.001 0.113 0.001 0.000 0.001 0.029 0.057 
   2 |  0.058 0.029 0.000 0.001 0.057 0.339 0.001 0.000 0.001 0.000 0.029 0.142 0.000 0.001 0.169 0.001 0.000 0.085 0.029 0.057 
   3 |  0.086 0.001 0.000 0.057 0.114 0.001 0.142 0.000 0.029 0.000 0.283 0.001 0.000 0.057 0.029 0.001 0.000 0.170 0.001 0.029 
   4 |  0.029 0.001 0.000 0.114 0.001 0.001 0.621 0.000 0.001 0.000 0.057 0.001 0.000 0.113 0.000 0.001 0.029 0.001 0.029 0.001 
   5 |  0.001 0.001 0.000 0.001 0.001 0.001 0.029 0.000 0.001 0.000 0.564 0.001 0.000 0.057 0.000 0.001 0.057 0.226 0.057 0.001 

   7 |  0.029 0.677 0.000 0.001 0.001 0.029 0.057 0.000 0.029 0.000 0.029 0.001 0.000 0.029 0.000 0.001 0.000 0.001 0.029 0.085 
   8 |  0.029 0.621 0.000 0.001 0.001 0.001 0.001 0.000 0.254 0.000 0.001 0.086 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.001 
   9 |  0.001 0.086 0.000 0.001 0.001 0.085 0.001 0.000 0.057 0.000 0.001 0.762 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.001 
  10 |  0.086 0.086 0.000 0.001 0.001 0.564 0.001 0.000 0.029 0.000 0.001 0.058 0.000 0.001 0.113 0.001 0.000 0.001 0.029 0.029 
  11 |  0.142 0.029 0.000 0.001 0.001 0.620 0.085 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.113 0.001 
  12 |  0.029 0.114 0.000 0.001 0.001 0.142 0.001 0.113 0.029 0.057 0.001 0.058 0.000 0.001 0.451 0.001 0.000 0.001 0.001 0.001 
  13 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.987 0.000 0.001 0.001 0.001 
  14 |  0.227 0.001 0.000 0.001 0.001 0.029 0.283 0.028 0.001 0.000 0.086 0.170 0.057 0.001 0.029 0.001 0.000 0.085 0.001 0.001 
  15 |  0.508 0.001 0.000 0.423 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.029 0.000 0.001 0.000 0.001 0.029 0.001 
  16 |  0.001 0.001 0.000 0.029 0.001 0.705 0.029 0.000 0.001 0.057 0.001 0.001 0.057 0.001 0.057 0.001 0.000 0.001 0.057 0.001 
  17 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.113 0.874 
  18 |  0.058 0.001 0.000 0.001 0.001 0.254 0.029 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.508 0.029 0.001 0.113 0.001 
  19 |  0.001 0.536 0.000 0.001 0.001 0.001 0.198 0.000 0.085 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.170 
  20 |  0.001 0.001 0.958 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.029 0.001 
  21 |  0.058 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.282 0.029 0.001 0.198 0.423 
  22 |  0.142 0.029 0.000 0.001 0.001 0.001 0.029 0.000 0.001 0.000 0.142 0.058 0.057 0.001 0.000 0.001 0.029 0.085 0.057 0.367 
  23 |  0.001 0.001 0.000 0.001 0.621 0.001 0.001 0.000 0.001 0.000 0.029 0.001 0.000 0.085 0.000 0.001 0.141 0.001 0.057 0.057 
  24 |  0.142 0.029 0.000 0.001 0.001 0.029 0.029 0.226 0.113 0.000 0.001 0.311 0.113 0.001 0.000 0.001 0.000 0.001 0.001 0.001 

  27 |  0.029 0.029 0.000 0.001 0.029 0.451 0.170 0.000 0.085 0.000 0.001 0.086 0.000 0.001 0.085 0.029 0.000 0.001 0.001 0.001 
  28 |  0.114 0.001 0.000 0.029 0.029 0.085 0.057 0.000 0.001 0.000 0.001 0.001 0.000 0.029 0.141 0.001 0.057 0.226 0.226 0.001 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



25 columns
Num Motifs: 35
   1,  1      37 ldfdk EFRDKTVVIVAIPGAFTPTCTANHIPPF vekft      64   1.00 F 1091044
   2,  1      32 irlsd YRGKKYVILFFYPANFTAISPTELMLLS drise      59   1.00 F 11467494
   6,  1      26 eikei DLKSNWNVFFFYPYSYSFICPLELKNIS nkike      53   1.00 F 13186328
   8,  1      28 kirls SYRGKWVVLFFYPADFTFVCPTEVEGFA edyek      55   1.00 F 13541053
   9,  1      26 mrkls EFRGQNVVLAFFPGAFTSVCTKEMCTFR dsman      53   1.00 F 13541117
  13,  1      25 melpd EFEGKWFILFSHPADFTPVCTTEFVAFQ evype      52   1.00 F 14286173
  15,  1      25 kirls DFRGRIVVLYFYPRAMTPGCTREGVRFN ellde      52   1.00 F 14600438
  22,  1      26 vtlrg YRGAKNVLLVFFPLAFTGICQGELDQLR dhlpe      53   1.00 F 15609375
  23,  1      30 nvsla DYRGRRVIVYFYPAASTPGCTKQACDFR dnlgd      57   1.00 F 15609658
  24,  1      24 tvsls DFKGKNIVLYFYPKDMTPGCTTEACDFR drved      51   1.00 F 15613511
  28,  1      20 tfthv DLYGKYTILFFFPKAGTSGCTREAVEFS renfe      47   1.00 F 15643152
  30,  1      61 gltda LADNRAVVLFFYPFDFSPVCATELCAIQ narwf      88   1.00 F 15790738
  35,  1      26 itlss YRGQSHVVLVFYPLDFSPVCSMQLPEYS gsqdd      53   1.00 F 15807234
  36,  1      28 vnlae LFKGKKGVLFGVPGAFTPGCSKTHLPGF veqae      55   1.00 F 15826629
  38,  1      26 vkips DFKGKVVVLAFYPAAFTSVCTKEMCTFR dsmak      53   1.00 F 15899339
  39,  1      30 vttel LFKGKRVVLFAVPGAFTPTCSLNHLPGY lenrd      57   1.00 F 15964668
  44,  1      50 fnlak ALKKGPVVLYFFPAAYTAGCTAEAREFA eatpe      77   1.00 F 16125919
  47,  1      31 fnfkq HTNGKTTVLFFWPMDFTFVCPSELIAFD kryee      58   1.00 F 16501671
  51,  1      33 slekn IEDDKWTILFFYPMDFTFVCPTEIVAIS arsde      60   1.00 F 16803644
  53,  1      31 vttdd LFAGKTVAVFSLPGAFTPTCSSTHLPGY nelak      58   1.00 F 17229033
  60,  1      28 rlsev LKRGRPVVLLFFPGAFTSVCTKELCTFR dkmal      55   1.00 F 18313548
  62,  1      26 eislq DYIGKYVVLAFYPLDFTFVCPTEINRFS dlkga      53   1.00 F 19173077
  67,  1      27 evtek DTEGRWSVFFFYPADFTFVCPTELGDVA dhyee      54   1.00 F 20151112
  68,  1      29 vdtht LFTGRKVVLFAVPGAFTPTCSAKHLPGY veqfe      56   1.00 F 21112072
  70,  1      32 qinhk TYEGQWKVVFAWPKDFTFVCPTEIAAFG klnde      59   1.00 F 21223405
  71,  1      28 eihly DLKGKKVLLSFHPLAWTQVCAQQMKSLE enyel      55   1.00 F 21227878
  73,  1      25 mvsls EFKGRKVLLIFYPGDDTPVCTAQLCDYR nnvaa      52   1.00 F 21674812
  81,  1      28 irfhd FLGDSWGILFSHPRDFTPVCTTELGRAA klape      55   1.00 F 3318841
  85,  1      10 eidin EYKGKYVVLLFYPLDWTFVCPTEMIGYS evagq      37   1.00 F 4433065
  86,  1      32 vsvhs IAAGKKVILFGVPGAFTPTCSMSHVPGF igkae      59   1.00 F 4704732
  87,  1      28 fdfyk YVGDNWAILFSHPHDFTPVCTTELAEFG kmhee      55   1.00 F 4996210
  88,  1      41 ynask EFANKKVVLFALPGAFTPVCSANHVPEY iqklp      68   1.00 F 5326864
  89,  1      88 slkki TENNRVVVFFVYPRASTPGCTRQACGFR dnyqe     115   1.00 F 6322180
  90,  1      43 ewskl ISENKKVIITGAPAAFSPTCTVSHIPGY inyld      70   1.00 F 6323138
  95,  1      31 evkls DYKGKYVVLFFYPLDFTFVCPTEIIAFS nraed      58   1.00 F 9955016
                       ***** ******************  **


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif d = -1668.31468
Log Fragmentation portion of MAP for motif d = -7.86327



=============================================================
====== ELEMENTS OCCURRING GREATER THAN  50% OF THE TIME =====
======                    Motif d                       =====
=============================================================


Listing of those elements occurring greater than 50% of the time
in near optimal sampling using 500 iterations



Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |     2  .   .   28  17   2  .    2   8  .   .   17  .   .   11  .   .   .    2   5  1.1
   2 |     5   2  .   .    5  34  .   .   .   .    2  14  .   .   17  .   .    8   2   5  1.4
   3 |     8  .   .    5  11  .   14  .    2  .   28  .   .    5   2  .   .   17  .    2  1.0
   4 |     2  .   .   11  .   .   62  .   .   .    5  .   .   11  .   .    2  .    2  .   2.1
   5 |    .   .   .   .   .   .    2  .   .   .   57  .   .    5  .   .    5  22   5  .   2.2

   7 |     2  68  .   .   .    2   5  .    2  .    2  .   .    2  .   .   .   .    2   8  1.9
   8 |     2  62  .   .   .   .   .   .   25  .   .    8  .   .   .   .   .   .   .   .   2.3
   9 |    .    8  .   .   .    8  .   .    5  .   .   77  .   .   .   .   .   .   .   .   2.3
  10 |     8   8  .   .   .   57  .   .    2  .   .    5  .   .   11  .   .   .    2   2  2.1
  11 |    14   2  .   .   .   62   8  .   .   .   .   .   .   .   .   .   .   .   11  .   2.4
  12 |     2  11  .   .   .   14  .   11   2   5  .    5  .   .   45  .   .   .   .   .   2.4
  13 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  100  .   .   .   .   4.3
  14 |    22  .   .   .   .    2  28   2  .   .    8  17   5  .    2  .   .    8  .   .   1.2
  15 |    51  .   .   42  .   .   .   .   .   .   .   .   .    2  .   .   .   .    2  .   2.3
  16 |    .   .   .    2  .   71   2  .   .    5  .   .    5  .    5  .   .   .    5  .   2.9
  17 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   11  88  3.6
  18 |     5  .   .   .   .   25   2  .   .   .   .   .   .   .   .   51   2  .   11  .   2.4
  19 |    .   54  .   .   .   .   20  .    8  .   .   .   .   .   .   .   .   .   .   17  2.0
  20 |    .   .   97  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .    2  .   6.2
  21 |     5  .   .   .   .   .   .   .   .   .   .   .   .   .   .   28   2  .   20  42  2.3
  22 |    14   2  .   .   .   .    2  .   .   .   14   5   5  .   .   .    2   8   5  37  1.3
  23 |    .   .   .   .   62  .   .   .   .   .    2  .   .    8  .   .   14  .    5   5  2.2
  24 |    14   2  .   .   .    2   2  22  11  .   .   31  11  .   .   .   .   .   .   .   1.8

  27 |     2   2  .   .    2  45  17  .    8  .   .    8  .   .    8   2  .   .   .   .   1.6
  28 |    11  .   .    2   2   8   5  .   .   .   .   .   .    2  14  .    5  22  22  .   1.4

nonsite    8   7  .    6   7   4   7   1   5  .    7   9   2   4   2   4   3   4   5   5
site       7   9   3   3   4  13   7   1   3  .    4   7   1   1   4   7   1   3   4   8

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.029 0.001 0.000 0.283 0.170 0.029 0.001 0.028 0.085 0.000 0.001 0.170 0.000 0.001 0.113 0.001 0.000 0.001 0.029 0.057 
   2 |  0.058 0.029 0.000 0.001 0.057 0.339 0.001 0.000 0.001 0.000 0.029 0.142 0.000 0.001 0.169 0.001 0.000 0.085 0.029 0.057 
   3 |  0.086 0.001 0.000 0.057 0.114 0.001 0.142 0.000 0.029 0.000 0.283 0.001 0.000 0.057 0.029 0.001 0.000 0.170 0.001 0.029 
   4 |  0.029 0.001 0.000 0.114 0.001 0.001 0.621 0.000 0.001 0.000 0.057 0.001 0.000 0.113 0.000 0.001 0.029 0.001 0.029 0.001 
   5 |  0.001 0.001 0.000 0.001 0.001 0.001 0.029 0.000 0.001 0.000 0.564 0.001 0.000 0.057 0.000 0.001 0.057 0.226 0.057 0.001 

   7 |  0.029 0.677 0.000 0.001 0.001 0.029 0.057 0.000 0.029 0.000 0.029 0.001 0.000 0.029 0.000 0.001 0.000 0.001 0.029 0.085 
   8 |  0.029 0.621 0.000 0.001 0.001 0.001 0.001 0.000 0.254 0.000 0.001 0.086 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.001 
   9 |  0.001 0.086 0.000 0.001 0.001 0.085 0.001 0.000 0.057 0.000 0.001 0.762 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.001 
  10 |  0.086 0.086 0.000 0.001 0.001 0.564 0.001 0.000 0.029 0.000 0.001 0.058 0.000 0.001 0.113 0.001 0.000 0.001 0.029 0.029 
  11 |  0.142 0.029 0.000 0.001 0.001 0.620 0.085 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.113 0.001 
  12 |  0.029 0.114 0.000 0.001 0.001 0.142 0.001 0.113 0.029 0.057 0.001 0.058 0.000 0.001 0.451 0.001 0.000 0.001 0.001 0.001 
  13 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.987 0.000 0.001 0.001 0.001 
  14 |  0.227 0.001 0.000 0.001 0.001 0.029 0.283 0.028 0.001 0.000 0.086 0.170 0.057 0.001 0.029 0.001 0.000 0.085 0.001 0.001 
  15 |  0.508 0.001 0.000 0.423 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.029 0.000 0.001 0.000 0.001 0.029 0.001 
  16 |  0.001 0.001 0.000 0.029 0.001 0.705 0.029 0.000 0.001 0.057 0.001 0.001 0.057 0.001 0.057 0.001 0.000 0.001 0.057 0.001 
  17 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.113 0.874 
  18 |  0.058 0.001 0.000 0.001 0.001 0.254 0.029 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.508 0.029 0.001 0.113 0.001 
  19 |  0.001 0.536 0.000 0.001 0.001 0.001 0.198 0.000 0.085 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.170 
  20 |  0.001 0.001 0.958 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.029 0.001 
  21 |  0.058 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.282 0.029 0.001 0.198 0.423 
  22 |  0.142 0.029 0.000 0.001 0.001 0.001 0.029 0.000 0.001 0.000 0.142 0.058 0.057 0.001 0.000 0.001 0.029 0.085 0.057 0.367 
  23 |  0.001 0.001 0.000 0.001 0.621 0.001 0.001 0.000 0.001 0.000 0.029 0.001 0.000 0.085 0.000 0.001 0.141 0.001 0.057 0.057 
  24 |  0.142 0.029 0.000 0.001 0.001 0.029 0.029 0.226 0.113 0.000 0.001 0.311 0.113 0.001 0.000 0.001 0.000 0.001 0.001 0.001 

  27 |  0.029 0.029 0.000 0.001 0.029 0.451 0.170 0.000 0.085 0.000 0.001 0.086 0.000 0.001 0.085 0.029 0.000 0.001 0.001 0.001 
  28 |  0.114 0.001 0.000 0.029 0.029 0.085 0.057 0.000 0.001 0.000 0.001 0.001 0.000 0.029 0.141 0.001 0.057 0.226 0.226 0.001 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



25 columns
Num Motifs: 35
   1,  1      37 ldfdk EFRDKTVVIVAIPGAFTPTCTANHIPPF vekft      64   1.00 F 1091044
   2,  1      32 irlsd YRGKKYVILFFYPANFTAISPTELMLLS drise      59   1.00 F 11467494
   6,  1      26 eikei DLKSNWNVFFFYPYSYSFICPLELKNIS nkike      53   0.98 F 13186328
   8,  1      28 kirls SYRGKWVVLFFYPADFTFVCPTEVEGFA edyek      55   1.00 F 13541053
   9,  1      26 mrkls EFRGQNVVLAFFPGAFTSVCTKEMCTFR dsman      53   1.00 F 13541117
  13,  1      25 melpd EFEGKWFILFSHPADFTPVCTTEFVAFQ evype      52   1.00 F 14286173
  15,  1      25 kirls DFRGRIVVLYFYPRAMTPGCTREGVRFN ellde      52   1.00 F 14600438
  22,  1      26 vtlrg YRGAKNVLLVFFPLAFTGICQGELDQLR dhlpe      53   1.00 F 15609375
  23,  1      30 nvsla DYRGRRVIVYFYPAASTPGCTKQACDFR dnlgd      57   1.00 F 15609658
  24,  1      24 tvsls DFKGKNIVLYFYPKDMTPGCTTEACDFR drved      51   1.00 F 15613511
  28,  1      20 tfthv DLYGKYTILFFFPKAGTSGCTREAVEFS renfe      47   1.00 F 15643152
  30,  1      61 gltda LADNRAVVLFFYPFDFSPVCATELCAIQ narwf      88   1.00 F 15790738
  35,  1      26 itlss YRGQSHVVLVFYPLDFSPVCSMQLPEYS gsqdd      53   1.00 F 15807234
  36,  1      28 vnlae LFKGKKGVLFGVPGAFTPGCSKTHLPGF veqae      55   1.00 F 15826629
  38,  1      26 vkips DFKGKVVVLAFYPAAFTSVCTKEMCTFR dsmak      53   1.00 F 15899339
  39,  1      30 vttel LFKGKRVVLFAVPGAFTPTCSLNHLPGY lenrd      57   1.00 F 15964668
  44,  1      50 fnlak ALKKGPVVLYFFPAAYTAGCTAEAREFA eatpe      77   1.00 F 16125919
  47,  1      31 fnfkq HTNGKTTVLFFWPMDFTFVCPSELIAFD kryee      58   1.00 F 16501671
  51,  1      33 slekn IEDDKWTILFFYPMDFTFVCPTEIVAIS arsde      60   1.00 F 16803644
  53,  1      31 vttdd LFAGKTVAVFSLPGAFTPTCSSTHLPGY nelak      58   1.00 F 17229033
  60,  1      28 rlsev LKRGRPVVLLFFPGAFTSVCTKELCTFR dkmal      55   1.00 F 18313548
  62,  1      26 eislq DYIGKYVVLAFYPLDFTFVCPTEINRFS dlkga      53   1.00 F 19173077
  67,  1      27 evtek DTEGRWSVFFFYPADFTFVCPTELGDVA dhyee      54   1.00 F 20151112
  68,  1      29 vdtht LFTGRKVVLFAVPGAFTPTCSAKHLPGY veqfe      56   1.00 F 21112072
  70,  1      32 qinhk TYEGQWKVVFAWPKDFTFVCPTEIAAFG klnde      59   1.00 F 21223405
  71,  1      28 eihly DLKGKKVLLSFHPLAWTQVCAQQMKSLE enyel      55   1.00 F 21227878
  73,  1      25 mvsls EFKGRKVLLIFYPGDDTPVCTAQLCDYR nnvaa      52   1.00 F 21674812
  81,  1      28 irfhd FLGDSWGILFSHPRDFTPVCTTELGRAA klape      55   1.00 F 3318841
  85,  1      10 eidin EYKGKYVVLLFYPLDWTFVCPTEMIGYS evagq      37   1.00 F 4433065
  86,  1      32 vsvhs IAAGKKVILFGVPGAFTPTCSMSHVPGF igkae      59   1.00 F 4704732
  87,  1      28 fdfyk YVGDNWAILFSHPHDFTPVCTTELAEFG kmhee      55   1.00 F 4996210
  88,  1      41 ynask EFANKKVVLFALPGAFTPVCSANHVPEY iqklp      68   1.00 F 5326864
  89,  1      88 slkki TENNRVVVFFVYPRASTPGCTRQACGFR dnyqe     115   1.00 F 6322180
  90,  1      43 ewskl ISENKKVIITGAPAAFSPTCTVSHIPGY inyld      70   0.99 F 6323138
  95,  1      31 evkls DYKGKYVVLFFYPLDFTFVCPTEIIAFS nraed      58   1.00 F 9955016
                       ***** ******************  **


Log Background portion of Map = -39912.17887
Log Alignment portion of Map = -957.33606
Log Site/seq portion of Map = 0.00000
Log Null Map = -46943.36311
Log Map = 2111.15797


log MAP = sum of motif and fragmentation parts of MAP + background + alignment + sites/seq - Null



=============================================================
======                Results by Sequence               =====
====== ELEMENTS OCCURRING GREATER THAN  50% OF THE TIME =====
=============================================================


   1,  1,  3      37 ldfdk EFRDKTVVIVAIPGAFTPTCTANHIPPF vekft      64   1.00 F 1091044
   2,  1,  3      32 irlsd YRGKKYVILFFYPANFTAISPTELMLLS drise      59   1.00 F 11467494
   2,  2,  0      72 klstq ILAISVDSPFSH              lqyll      83   1.00 F 11467494
   2,  3,  1     161 riles IQYVKENPGYACPVNWNFG       dqvfy     179   1.00 F 11467494
   3,  1,  2      17 viqsd KLVVVDFYADWCMPCRYISPILEKL skeyn      41   1.00 F 11499727
   4,  1,  2      22 llntt QYVVADFYADWCGPCKAIAPMYAQF aktfs      46   1.00 F 1174686
   5,  1,  2      19 ifsak KNVIVDFWAAWCGPCKLTSPEFQKA adefs      43   1.00 F 12044976
   6,  1,  3      26 eikei DLKSNWNVFFFYPYSYSFICPLELKNIS nkike      53   0.98 F 13186328
   6,  2,  0      66 nlntk IYAISNDSHFVQ              knwie      77   1.00 F 13186328
   7,  1,  2      21 kenhs KPILIDFYADWCPPCRMLIPVLDSI ekkhg      45   1.00 F 13358154
   8,  1,  3      28 kirls SYRGKWVVLFFYPADFTFVCPTEVEGFA edyek      55   1.00 F 13541053
   8,  2,  0      68 kknte VISVSEDTVYVH              kawvq      79   1.00 F 13541053
   9,  1,  3      26 mrkls EFRGQNVVLAFFPGAFTSVCTKEMCTFR dsman      53   1.00 F 13541117
   9,  2,  0      66 kfkak VIGISVDSPFSL              aefak      77   1.00 F 13541117
  10,  1,  2      17 aetse GVVLADFWAPWCGPCKMIAPVLEEL dqemg      41   1.00 F 135765
  11,  1,  2      29 akesn KLIVIDFTASWCPPCRMIAPIFNDL akkfm      53   1.00 F 1388082
  12,  1,  2      44 likqn DKLVIDFYATWCGPCKMMQPHLTKL iqayp      68   1.00 F 140543
  13,  1,  3      25 melpd EFEGKWFILFSHPADFTPVCTTEFVAFQ evype      52   1.00 F 14286173
  13,  2,  0      65 eldce LVGLSVDQVFSH              ikwie      76   0.98 F 14286173
  14,  1,  2      80 selrg KVVMLQFTASWCGVCRKEMPFIEKD iwlkh     104   1.00 F 14578634
  15,  1,  3      25 kirls DFRGRIVVLYFYPRAMTPGCTREGVRFN ellde      52   1.00 F 14600438
  15,  2,  0      65 klgav VIGVSTDSVEKN              rkfae      76   1.00 F 14600438
  16,  1,  2      23 lqnsd KPVLVDFYATWCGPCQLMVPILNEV setlk      47   1.00 F 15218394
  17,  1,  2     157 adfrg RPLVINLWASWCPPCRREMPVLQQA qaenp     181   1.00 F 15597673
  18,  1,  2      26 ensfh KPVLVDFWADWCAPCKALMPLLAQI aesyq      50   1.00 F 15599256
  19,  1,  2      67 negkg KTILLNFWSETCGVCIAELKTFEQL lqsyp      91   1.00 F 15602312
  20,  1,  2      61 eefkg KVLLINFWATWCPPCKEEIPMFKEI yekyr      85   1.00 F 15605725
  21,  1,  0      80 megvd VTVVSMDLPFAQ              krfce      91   0.65 F 15605963
  22,  1,  3      26 vtlrg YRGAKNVLLVFFPLAFTGICQGELDQLR dhlpe      53   1.00 F 15609375
  23,  1,  3      30 nvsla DYRGRRVIVYFYPAASTPGCTKQACDFR dnlgd      57   1.00 F 15609658
  23,  2,  0      70 tagln VVGISPDKPEKL              atfrd      81   0.99 F 15609658
  24,  1,  3      24 tvsls DFKGKNIVLYFYPKDMTPGCTTEACDFR drved      51   1.00 F 15613511
  24,  2,  0      64 glntv ILGVSPDPVERH              kkfie      75   1.00 F 15613511
  25,  1,  2      60 sdyrg DVVILNVWASWCEPCRKEMPALMEL qsdye      84   1.00 F 15614085
  26,  1,  2      63 releg KGVFLNFWGTYCPPCEREMPHMEKL ygeyk      87   1.00 F 15614140
  27,  1,  2      72 sslrg QPVILHFFATWCPVCQDEMPSLVKL dkeyr      96   1.00 F 15615431
  28,  1,  3      20 tfthv DLYGKYTILFFFPKAGTSGCTREAVEFS renfe      47   1.00 F 15643152
  28,  2,  0      56 fekaq VVGISRDSVEAL              krfke      67   1.00 F 15643152
  30,  1,  3      61 gltda LADNRAVVLFFYPFDFSPVCATELCAIQ narwf      88   1.00 F 15790738
  30,  2,  0     101 tpgla VWGISPDSTYAH              eafad     112   0.98 F 15790738
  31,  1,  2       2     m TVTLKDFYADWCGPCKTQDPILEEL eadyd      26   1.00 F 15791337
  32,  1,  0      80 idntv VLCISADLPFAQ              srfcg      91   0.99 F 15801846
  33,  1,  2      72 adyrg RPVVLNFWASWCGPCREEAPLFAKL aahpg      96   1.00 F 15805225
  34,  1,  2      78 taaqg KPVVINFWASWCVPCRQEAPLFSKL sqeta     102   1.00 F 15805374
  35,  1,  3      26 itlss YRGQSHVVLVFYPLDFSPVCSMQLPEYS gsqdd      53   1.00 F 15807234
  35,  2,  0      66 eagav VLGINRDSVYAH              rawaa      77   1.00 F 15807234
  36,  1,  3      28 vnlae LFKGKKGVLFGVPGAFTPGCSKTHLPGF veqae      55   1.00 F 15826629
  37,  1,  2      49 fitkn KIVVVDFWAEWCAPCLILAPVIEEL andyp      73   1.00 F 15899007
  38,  1,  3      26 vkips DFKGKVVVLAFYPAAFTSVCTKEMCTFR dsmak      53   1.00 F 15899339
  38,  2,  0      66 evnav VIGISVDPPFSN              kafke      77   1.00 F 15899339
  39,  1,  3      30 vttel LFKGKRVVLFAVPGAFTPTCSLNHLPGY lenrd      57   1.00 F 15964668
  40,  1,  2      61 easrq QPVLVDFWAPWCGPCKQLTPVIEKV vreaa      85   1.00 F 15966937
  41,  1,  2      61 sdfrg KTLLVNLWATWCVPCRKEMPALDEL qgkls      85   1.00 F 15988313
  42,  1,  2      60 qdakg KKVLLNFWATWCKPCRQEMPAMEKL qkeya      84   1.00 F 16078864
  43,  1,  2      53 llqdd LPMVIDFWAPWCGPCRSFAPIFAET aaera      77   1.00 F 16123427
  44,  1,  3      50 fnlak ALKKGPVVLYFFPAAYTAGCTAEAREFA eatpe      77   1.00 F 16125919
  46,  1,  2      21 vlkad GAILVDFWAEWCGPCKMIAPILDEI adeyq      45   1.00 F 1633495
  47,  1,  3      31 fnfkq HTNGKTTVLFFWPMDFTFVCPSELIAFD kryee      58   1.00 F 16501671
  47,  2,  0      71 krgve VVGVSFDSEFVH              nawrn      82   1.00 F 16501671
  47,  3,  1     160 lrmvd ALQFHEEHGDVCPAQWEKG       kegmn     178   1.00 F 16501671
  48,  1,  2      34 vlqcp KPILVYFGAPWCGLCHFVKPLLNHL hgewq      58   1.00 F 1651717
  49,  1,  2      60 tlsee RPVLLYFWASWCGVCRFTTPAVAHL aaege      84   1.00 F 16759994
  50,  1,  2      53 llkdd LPVVIDFWAPWCGPCRNFAPIFEDV aeers      77   1.00 F 16761507
  51,  1,  3      33 slekn IEDDKWTILFFYPMDFTFVCPTEIVAIS arsde      60   1.00 F 16803644
  52,  1,  2      19 iissh PKILLNFWAEWCAPCRCFWPTLEQF aemee      43   1.00 F 16804867
  53,  1,  3      31 vttdd LFAGKTVAVFSLPGAFTPTCSSTHLPGY nelak      58   1.00 F 17229033
  54,  1,  2      22 vlsed KVVVVDFTATWCGPCRLVSPLMDQL adeyk      46   1.00 F 17229859
  55,  1,  2      18 vlegt GYVLVDYFSDGCVPCKALMPAVEEL skkye      42   1.00 F 1729944
  56,  1,  2      28 rqhpe KIIILDFYATWCGPCKAIAPLYKEL atthk      52   1.00 F 17531233
  57,  1,  2      27 ehlkg KIIGLYFSASWCPPCRAFTPKLKEF feeik      51   1.00 F 17537401
  58,  1,  2      63 safrg QPVVINFWAPWCGPCVEEMPELSAL aqeqk      87   1.00 F 17547503
  59,  1,  2     286 seykg KTIFLNFWATWCPPCRGEMPYIDEL ykeyn     310   1.00 F 18309723
  60,  1,  3      28 rlsev LKRGRPVVLLFFPGAFTSVCTKELCTFR dkmal      55   1.00 F 18313548
  60,  2,  0      68 kanae VLAISVDSPFAL              kafkd      79   1.00 F 18313548
  61,  1,  2      44 dsllg KKIGLYFSAAWCGPCQRFTPQLVEV ynels      68   1.00 F 18406743
  61,  2,  2     364 sdlvg KTILMYFSAHWCPPCRAFTPKLVEV ykqik     388   1.00 F 18406743
  62,  1,  3      26 eislq DYIGKYVVLAFYPLDFTFVCPTEINRFS dlkga      53   1.00 F 19173077
  62,  2,  0      66 rrnav VLLISCDSVYTH              kawas      77   1.00 F 19173077
  63,  1,  2      15 sdfeg EVVVLNAWGQWCAPCRAEVDDLQLV qetld      39   1.00 F 19554157
  64,  1,  2      39 eeykg KVVVINFWATWCGYCVEEMPGFEKV ykefg      63   1.00 F 19705357
  66,  1,  2       7 agdfm KPMLLDFSATWCGPCRMQKPILEEL ekkyg      31   1.00 F 20092028
  67,  1,  3      27 evtek DTEGRWSVFFFYPADFTFVCPTELGDVA dhyee      54   1.00 F 20151112
  67,  2,  0      67 klgvd VYSVSTDTHFTH              kawhs      78   1.00 F 20151112
  67,  3,  1     154 rkika AQYVAAHPGEVCPAKWKEG       eatla     172   1.00 F 20151112
  68,  1,  3      29 vdtht LFTGRKVVLFAVPGAFTPTCSAKHLPGY veqfe      56   1.00 F 21112072
  69,  1,  2     103 adykg KVVVLNVWGSWCPPCRAEAKNFEKV yqdvk     127   1.00 F 21222859
  70,  1,  3      32 qinhk TYEGQWKVVFAWPKDFTFVCPTEIAAFG klnde      59   1.00 F 21223405
  70,  2,  0      72 drdaq ILGFSGDSEFVH              hawrk      83   1.00 F 21223405
  71,  1,  3      28 eihly DLKGKKVLLSFHPLAWTQVCAQQMKSLE enyel      55   1.00 F 21227878
  72,  1,  0      78 keegi VLTISADLPFAQ              krwca      89   0.99 F 21283385
  73,  1,  3      25 mvsls EFKGRKVLLIFYPGDDTPVCTAQLCDYR nnvaa      52   1.00 F 21674812
  73,  2,  0      65 srgit VIGISGDSPESH              kqfae      76   1.00 F 21674812
  74,  1,  2      53 sdfkg ERVLINFWTTWCPPCRQEMPDMQRF yqdlq      77   1.00 F 23098307
  76,  1,  2      20 kylqh QRVVVDFSAEWCGPCRAIAPVFDKL sneft      44   1.00 F 267116
  77,  1,  2      81 aafkg KVSLVNVWASWCVPCHDEAPLLTEL gkdkr     105   1.00 F 27375582
  78,  1,  2      34 vtsdn DVVLADFYADWCGPCQMLEPVVETL aeqtd      58   1.00 F 2822332
  79,  1,  2      77 sdlkg KKVILNFWATWCGPCQQEMPDMEAF ykehk     101   1.00 F 30021713
  80,  1,  2      19 tisan SNVLVYFWAPLCAPCDLFTPTYEAS srkhf      43   1.00 F 3261501
  81,  1,  3      28 irfhd FLGDSWGILFSHPRDFTPVCTTELGRAA klape      55   1.00 F 3318841
  81,  2,  0      68 krnvk LIALSIDSVEDH              lawsk      79   1.00 F 3318841
  81,  3,  1     166 lrvvi SLQLTAEKRVATPVDWKDG       dsvmv     184   1.00 F 3318841
  82,  1,  2      19 tietn PLVIVDFWAPWCGSCKMLGPVLEEV esevg      43   1.00 F 3323237
  83,  1,  2      17 ektah QAVVVNVGASWCPDCRKIEPIMENL aktyk      41   1.00 F 4155972
  84,  1,  2      79 vvnse TPVVVDFHAQWCGPCKILGPRLEKM vakqh     103   1.00 F 4200327
  85,  1,  3      10 eidin EYKGKYVVLLFYPLDWTFVCPTEMIGYS evagq      37   1.00 F 4433065
  85,  2,  0      50 eince VIGVSVDSVYCH              qawce      61   1.00 F 4433065
  86,  1,  3      32 vsvhs IAAGKKVILFGVPGAFTPTCSMSHVPGF igkae      59   1.00 F 4704732
  87,  1,  3      28 fdfyk YVGDNWAILFSHPHDFTPVCTTELAEFG kmhee      55   1.00 F 4996210
  87,  2,  0      68 klnck LIGFSCNSKESH              dqwie      79   0.56 F 4996210
  87,  3,  1     163 lrvlk SLQLTNTHPVATPVNWKEG       dkcci     181   1.00 F 4996210
  88,  1,  3      41 ynask EFANKKVVLFALPGAFTPVCSANHVPEY iqklp      68   1.00 F 5326864
  89,  1,  3      88 slkki TENNRVVVFFVYPRASTPGCTRQACGFR dnyqe     115   1.00 F 6322180
  89,  2,  0     127 kkyaa VFGLSADSVTSQ              kkfqs     138   0.51 F 6322180
  90,  1,  3      43 ewskl ISENKKVIITGAPAAFSPTCTVSHIPGY inyld      70   0.99 F 6323138
  91,  1,  2      20 nenkg RLIVVDFFAQWCGPCRNIAPKVEAL akeip      44   1.00 F 6687568
  92,  1,  0      68 klgve VLSVSVDSVFVH              kmwnd      79   1.00 F 6850955
  93,  1,  2      18 llttn KKVVVDFYANWCGPCKILGPIFEEV aqdkk      42   1.00 F 7109697
  94,  1,  2      21 ilaed KLVVIDFYADWCGPCKIIAPKLDEL aqqys      45   1.00 F 7290567
  95,  1,  3      31 evkls DYKGKYVVLFFYPLDFTFVCPTEIIAFS nraed      58   1.00 F 9955016
  95,  2,  0      71 klgce VLGVSVDSQFTH              lawin      82   1.00 F 9955016
  95,  3,  1     160 lrlvq AFQYTDEHGEVCPAGWKPG       sdtik     178   1.00 F 9955016
  96,  1,  2      49 adlqg KVTLINFWFPSCPGCVSEMPKIIKT andyk      73   1.00 F 15677788
122 motifs


Column 1 :  Sequence Number, Site Number
Column 2 :  Motif type
Column 3 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input






======================== MAP MAXIMIZATION RESULTS ====================
======================================================================


=============================================================
======                Results by Sequence               =====
=============================================================


   1,  1,  3      37 ldfdk EFRDKTVVIVAIPGAFTPTCTANHIPPF vekft      64   1.00 F 1091044
   1,  2,  0      79 agvda VIVLSANDPFVQ              safgk      90   0.48 F 1091044
   2,  1,  3      32 irlsd YRGKKYVILFFYPANFTAISPTELMLLS drise      59   1.00 F 11467494
   2,  2,  0      72 klstq ILAISVDSPFSH              lqyll      83   1.00 F 11467494
   2,  3,  1     161 riles IQYVKENPGYACPVNWNFG       dqvfy     179   1.00 F 11467494
   3,  1,  2      17 viqsd KLVVVDFYADWCMPCRYISPILEKL skeyn      41   1.00 F 11499727
   4,  1,  2      22 llntt QYVVADFYADWCGPCKAIAPMYAQF aktfs      46   1.00 F 1174686
   5,  1,  2      19 ifsak KNVIVDFWAAWCGPCKLTSPEFQKA adefs      43   1.00 F 12044976
   6,  1,  3      26 eikei DLKSNWNVFFFYPYSYSFICPLELKNIS nkike      53   0.98 F 13186328
   6,  2,  0      66 nlntk IYAISNDSHFVQ              knwie      77   1.00 F 13186328
   7,  1,  2      21 kenhs KPILIDFYADWCPPCRMLIPVLDSI ekkhg      45   1.00 F 13358154
   8,  1,  3      28 kirls SYRGKWVVLFFYPADFTFVCPTEVEGFA edyek      55   1.00 F 13541053
   8,  2,  0      68 kknte VISVSEDTVYVH              kawvq      79   1.00 F 13541053
   9,  1,  3      26 mrkls EFRGQNVVLAFFPGAFTSVCTKEMCTFR dsman      53   1.00 F 13541117
   9,  2,  0      66 kfkak VIGISVDSPFSL              aefak      77   1.00 F 13541117
  10,  1,  2      17 aetse GVVLADFWAPWCGPCKMIAPVLEEL dqemg      41   1.00 F 135765
  11,  1,  2      29 akesn KLIVIDFTASWCPPCRMIAPIFNDL akkfm      53   1.00 F 1388082
  12,  1,  2      44 likqn DKLVIDFYATWCGPCKMMQPHLTKL iqayp      68   1.00 F 140543
  13,  1,  3      25 melpd EFEGKWFILFSHPADFTPVCTTEFVAFQ evype      52   1.00 F 14286173
  13,  2,  0      65 eldce LVGLSVDQVFSH              ikwie      76   0.98 F 14286173
  14,  1,  2      80 selrg KVVMLQFTASWCGVCRKEMPFIEKD iwlkh     104   1.00 F 14578634
  15,  1,  3      25 kirls DFRGRIVVLYFYPRAMTPGCTREGVRFN ellde      52   1.00 F 14600438
  15,  2,  0      65 klgav VIGVSTDSVEKN              rkfae      76   1.00 F 14600438
  16,  1,  2      23 lqnsd KPVLVDFYATWCGPCQLMVPILNEV setlk      47   1.00 F 15218394
  17,  1,  2     157 adfrg RPLVINLWASWCPPCRREMPVLQQA qaenp     181   1.00 F 15597673
  18,  1,  2      26 ensfh KPVLVDFWADWCAPCKALMPLLAQI aesyq      50   1.00 F 15599256
  19,  1,  2      67 negkg KTILLNFWSETCGVCIAELKTFEQL lqsyp      91   1.00 F 15602312
  20,  1,  2      61 eefkg KVLLINFWATWCPPCKEEIPMFKEI yekyr      85   1.00 F 15605725
  21,  1,  0      80 megvd VTVVSMDLPFAQ              krfce      91   0.65 F 15605963
  22,  1,  3      26 vtlrg YRGAKNVLLVFFPLAFTGICQGELDQLR dhlpe      53   1.00 F 15609375
  23,  1,  3      30 nvsla DYRGRRVIVYFYPAASTPGCTKQACDFR dnlgd      57   1.00 F 15609658
  23,  2,  0      70 tagln VVGISPDKPEKL              atfrd      81   0.99 F 15609658
  24,  1,  3      24 tvsls DFKGKNIVLYFYPKDMTPGCTTEACDFR drved      51   1.00 F 15613511
  24,  2,  0      64 glntv ILGVSPDPVERH              kkfie      75   1.00 F 15613511
  25,  1,  2      60 sdyrg DVVILNVWASWCEPCRKEMPALMEL qsdye      84   1.00 F 15614085
  26,  1,  2      63 releg KGVFLNFWGTYCPPCEREMPHMEKL ygeyk      87   1.00 F 15614140
  27,  1,  2      72 sslrg QPVILHFFATWCPVCQDEMPSLVKL dkeyr      96   1.00 F 15615431
  28,  1,  3      20 tfthv DLYGKYTILFFFPKAGTSGCTREAVEFS renfe      47   1.00 F 15643152
  28,  2,  0      56 fekaq VVGISRDSVEAL              krfke      67   1.00 F 15643152
  30,  1,  3      61 gltda LADNRAVVLFFYPFDFSPVCATELCAIQ narwf      88   1.00 F 15790738
  30,  2,  0     101 tpgla VWGISPDSTYAH              eafad     112   0.98 F 15790738
  31,  1,  2       2     m TVTLKDFYADWCGPCKTQDPILEEL eadyd      26   1.00 F 15791337
  32,  1,  0      80 idntv VLCISADLPFAQ              srfcg      91   0.99 F 15801846
  33,  1,  2      72 adyrg RPVVLNFWASWCGPCREEAPLFAKL aahpg      96   1.00 F 15805225
  34,  1,  2      78 taaqg KPVVINFWASWCVPCRQEAPLFSKL sqeta     102   1.00 F 15805374
  35,  1,  3      26 itlss YRGQSHVVLVFYPLDFSPVCSMQLPEYS gsqdd      53   1.00 F 15807234
  35,  2,  0      66 eagav VLGINRDSVYAH              rawaa      77   1.00 F 15807234
  36,  1,  3      28 vnlae LFKGKKGVLFGVPGAFTPGCSKTHLPGF veqae      55   1.00 F 15826629
  37,  1,  2      49 fitkn KIVVVDFWAEWCAPCLILAPVIEEL andyp      73   1.00 F 15899007
  38,  1,  3      26 vkips DFKGKVVVLAFYPAAFTSVCTKEMCTFR dsmak      53   1.00 F 15899339
  38,  2,  0      66 evnav VIGISVDPPFSN              kafke      77   1.00 F 15899339
  39,  1,  3      30 vttel LFKGKRVVLFAVPGAFTPTCSLNHLPGY lenrd      57   1.00 F 15964668
  40,  1,  2      61 easrq QPVLVDFWAPWCGPCKQLTPVIEKV vreaa      85   1.00 F 15966937
  41,  1,  2      61 sdfrg KTLLVNLWATWCVPCRKEMPALDEL qgkls      85   1.00 F 15988313
  42,  1,  2      60 qdakg KKVLLNFWATWCKPCRQEMPAMEKL qkeya      84   1.00 F 16078864
  43,  1,  2      53 llqdd LPMVIDFWAPWCGPCRSFAPIFAET aaera      77   1.00 F 16123427
  44,  1,  3      50 fnlak ALKKGPVVLYFFPAAYTAGCTAEAREFA eatpe      77   1.00 F 16125919
  46,  1,  2      21 vlkad GAILVDFWAEWCGPCKMIAPILDEI adeyq      45   1.00 F 1633495
  47,  1,  3      31 fnfkq HTNGKTTVLFFWPMDFTFVCPSELIAFD kryee      58   1.00 F 16501671
  47,  2,  0      71 krgve VVGVSFDSEFVH              nawrn      82   1.00 F 16501671
  47,  3,  1     160 lrmvd ALQFHEEHGDVCPAQWEKG       kegmn     178   1.00 F 16501671
  48,  1,  2      34 vlqcp KPILVYFGAPWCGLCHFVKPLLNHL hgewq      58   1.00 F 1651717
  49,  1,  2      60 tlsee RPVLLYFWASWCGVCRFTTPAVAHL aaege      84   1.00 F 16759994
  50,  1,  2      53 llkdd LPVVIDFWAPWCGPCRNFAPIFEDV aeers      77   1.00 F 16761507
  51,  1,  3      33 slekn IEDDKWTILFFYPMDFTFVCPTEIVAIS arsde      60   1.00 F 16803644
  52,  1,  2      19 iissh PKILLNFWAEWCAPCRCFWPTLEQF aemee      43   1.00 F 16804867
  53,  1,  3      31 vttdd LFAGKTVAVFSLPGAFTPTCSSTHLPGY nelak      58   1.00 F 17229033
  53,  2,  0      73 ngvde IVCISVNDAFVM              newak      84   0.32 F 17229033
  54,  1,  2      22 vlsed KVVVVDFTATWCGPCRLVSPLMDQL adeyk      46   1.00 F 17229859
  55,  1,  2      18 vlegt GYVLVDYFSDGCVPCKALMPAVEEL skkye      42   1.00 F 1729944
  56,  1,  2      28 rqhpe KIIILDFYATWCGPCKAIAPLYKEL atthk      52   1.00 F 17531233
  57,  1,  2      27 ehlkg KIIGLYFSASWCPPCRAFTPKLKEF feeik      51   1.00 F 17537401
  58,  1,  2      63 safrg QPVVINFWAPWCGPCVEEMPELSAL aqeqk      87   1.00 F 17547503
  59,  1,  2     286 seykg KTIFLNFWATWCPPCRGEMPYIDEL ykeyn     310   1.00 F 18309723
  60,  1,  3      28 rlsev LKRGRPVVLLFFPGAFTSVCTKELCTFR dkmal      55   1.00 F 18313548
  60,  2,  0      68 kanae VLAISVDSPFAL              kafkd      79   1.00 F 18313548
  61,  1,  2      44 dsllg KKIGLYFSAAWCGPCQRFTPQLVEV ynels      68   1.00 F 18406743
  61,  2,  2     364 sdlvg KTILMYFSAHWCPPCRAFTPKLVEV ykqik     388   1.00 F 18406743
  62,  1,  3      26 eislq DYIGKYVVLAFYPLDFTFVCPTEINRFS dlkga      53   1.00 F 19173077
  62,  2,  0      66 rrnav VLLISCDSVYTH              kawas      77   1.00 F 19173077
  63,  1,  2      15 sdfeg EVVVLNAWGQWCAPCRAEVDDLQLV qetld      39   1.00 F 19554157
  64,  1,  2      39 eeykg KVVVINFWATWCGYCVEEMPGFEKV ykefg      63   1.00 F 19705357
  66,  1,  2       7 agdfm KPMLLDFSATWCGPCRMQKPILEEL ekkyg      31   1.00 F 20092028
  67,  1,  3      27 evtek DTEGRWSVFFFYPADFTFVCPTELGDVA dhyee      54   1.00 F 20151112
  67,  2,  0      67 klgvd VYSVSTDTHFTH              kawhs      78   1.00 F 20151112
  67,  3,  1     154 rkika AQYVAAHPGEVCPAKWKEG       eatla     172   1.00 F 20151112
  68,  1,  3      29 vdtht LFTGRKVVLFAVPGAFTPTCSAKHLPGY veqfe      56   1.00 F 21112072
  69,  1,  2     103 adykg KVVVLNVWGSWCPPCRAEAKNFEKV yqdvk     127   1.00 F 21222859
  70,  1,  3      32 qinhk TYEGQWKVVFAWPKDFTFVCPTEIAAFG klnde      59   1.00 F 21223405
  70,  2,  0      72 drdaq ILGFSGDSEFVH              hawrk      83   1.00 F 21223405
  71,  1,  3      28 eihly DLKGKKVLLSFHPLAWTQVCAQQMKSLE enyel      55   1.00 F 21227878
  72,  1,  0      78 keegi VLTISADLPFAQ              krwca      89   0.99 F 21283385
  73,  1,  3      25 mvsls EFKGRKVLLIFYPGDDTPVCTAQLCDYR nnvaa      52   1.00 F 21674812
  73,  2,  0      65 srgit VIGISGDSPESH              kqfae      76   1.00 F 21674812
  74,  1,  2      53 sdfkg ERVLINFWTTWCPPCRQEMPDMQRF yqdlq      77   1.00 F 23098307
  76,  1,  2      20 kylqh QRVVVDFSAEWCGPCRAIAPVFDKL sneft      44   1.00 F 267116
  77,  1,  2      81 aafkg KVSLVNVWASWCVPCHDEAPLLTEL gkdkr     105   1.00 F 27375582
  78,  1,  2      34 vtsdn DVVLADFYADWCGPCQMLEPVVETL aeqtd      58   1.00 F 2822332
  79,  1,  2      77 sdlkg KKVILNFWATWCGPCQQEMPDMEAF ykehk     101   1.00 F 30021713
  80,  1,  2      19 tisan SNVLVYFWAPLCAPCDLFTPTYEAS srkhf      43   1.00 F 3261501
  81,  1,  3      28 irfhd FLGDSWGILFSHPRDFTPVCTTELGRAA klape      55   1.00 F 3318841
  81,  2,  0      68 krnvk LIALSIDSVEDH              lawsk      79   1.00 F 3318841
  81,  3,  1     166 lrvvi SLQLTAEKRVATPVDWKDG       dsvmv     184   1.00 F 3318841
  82,  1,  2      19 tietn PLVIVDFWAPWCGSCKMLGPVLEEV esevg      43   1.00 F 3323237
  83,  1,  2      17 ektah QAVVVNVGASWCPDCRKIEPIMENL aktyk      41   1.00 F 4155972
  84,  1,  2      79 vvnse TPVVVDFHAQWCGPCKILGPRLEKM vakqh     103   1.00 F 4200327
  85,  1,  3      10 eidin EYKGKYVVLLFYPLDWTFVCPTEMIGYS evagq      37   1.00 F 4433065
  85,  2,  0      50 eince VIGVSVDSVYCH              qawce      61   1.00 F 4433065
  86,  1,  3      32 vsvhs IAAGKKVILFGVPGAFTPTCSMSHVPGF igkae      59   1.00 F 4704732
  86,  2,  0      74 kgide IICFSVNDPFVM              kawgk      85   0.43 F 4704732
  87,  1,  3      28 fdfyk YVGDNWAILFSHPHDFTPVCTTELAEFG kmhee      55   1.00 F 4996210
  87,  2,  0      68 klnck LIGFSCNSKESH              dqwie      79   0.56 F 4996210
  87,  3,  1     163 lrvlk SLQLTNTHPVATPVNWKEG       dkcci     181   1.00 F 4996210
  88,  1,  3      41 ynask EFANKKVVLFALPGAFTPVCSANHVPEY iqklp      68   1.00 F 5326864
  89,  1,  3      88 slkki TENNRVVVFFVYPRASTPGCTRQACGFR dnyqe     115   1.00 F 6322180
  90,  1,  3      43 ewskl ISENKKVIITGAPAAFSPTCTVSHIPGY inyld      70   0.99 F 6323138
  91,  1,  2      20 nenkg RLIVVDFFAQWCGPCRNIAPKVEAL akeip      44   1.00 F 6687568
  92,  1,  0      68 klgve VLSVSVDSVFVH              kmwnd      79   1.00 F 6850955
  93,  1,  2      18 llttn KKVVVDFYANWCGPCKILGPIFEEV aqdkk      42   1.00 F 7109697
  94,  1,  2      21 ilaed KLVVIDFYADWCGPCKIIAPKLDEL aqqys      45   1.00 F 7290567
  95,  1,  3      31 evkls DYKGKYVVLFFYPLDFTFVCPTEIIAFS nraed      58   1.00 F 9955016
  95,  2,  0      71 klgce VLGVSVDSQFTH              lawin      82   1.00 F 9955016
  95,  3,  1     160 lrlvq AFQYTDEHGEVCPAGWKPG       sdtik     178   1.00 F 9955016
  96,  1,  2      49 adlqg KVTLINFWFPSCPGCVSEMPKIIKT andyk      73   1.00 F 15677788
124 motifs




Column 1 :  Sequence Number, Site Number
Column 2 :  Motif type
Column 3 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

-------------------------------------------------------------------------
                          MOTIF a





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    .   68  .   .   .   .   .   .   20  .   .   10  .   .   .   .   .   .   .   .   2.4
   2 |    .   17  .   .   .   .   .   .   34   3  .   34  .   .    6  .   .   .   .    3  1.7
   3 |    13   6  10  .   .   .   51  .   .   .   .    3  .   .   .   .   .   .   10   3  1.9
   4 |    .   31  .   .   .   10  .   .   48  .   .   10  .   .   .   .   .   .   .   .   2.1
   5 |    .   .   .   .   .   .   .   .   .   .   .   .   .    3  .   .   .   .   96  .   3.8

   7 |    .   .   .   86  .   .   .   .   .   .   .   .   .   13  .   .   .   .   .   .   3.2
   8 |    .   .   .   10  .   .   .   .   .   .    3  10  .   .   .    6   3  .   58   6  2.0
   9 |     3  34  .   .    6  .   .    6  .   .    3  .   .   .   .   37   3  .   .    3  1.8
  10 |    .   .   .   .   24  58  .   .   .   .   .   .   .   .   17  .   .   .   .   .   2.9

  12 |    .   .   .   .   .   .   .   55  .   .   .   13   6   6  .   .   17  .   .   .   3.4

nonsite    8   8  .    6   7   4   7   1   6  .    7   9   2   4   2   4   3   4   5   5
site       1  15   1   9   3   6   5   6  10  .   .    8  .    2   2   4   2  .   16   1

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.001 0.679 0.000 0.001 0.001 0.001 0.001 0.000 0.204 0.000 0.001 0.103 0.000 0.001 0.000 0.001 0.001 0.001 0.001 0.001 
   2 |  0.001 0.171 0.000 0.001 0.001 0.001 0.001 0.000 0.340 0.034 0.001 0.341 0.000 0.001 0.068 0.001 0.001 0.001 0.001 0.035 
   3 |  0.137 0.069 0.102 0.001 0.001 0.001 0.510 0.000 0.001 0.000 0.001 0.035 0.000 0.001 0.000 0.001 0.001 0.001 0.103 0.035 
   4 |  0.001 0.306 0.000 0.001 0.001 0.103 0.001 0.000 0.476 0.000 0.001 0.103 0.000 0.001 0.000 0.001 0.001 0.001 0.001 0.001 
   5 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.035 0.000 0.001 0.001 0.001 0.950 0.001 

   7 |  0.001 0.001 0.000 0.849 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.136 0.000 0.001 0.001 0.001 0.001 0.001 
   8 |  0.001 0.001 0.000 0.103 0.001 0.001 0.001 0.000 0.001 0.000 0.035 0.103 0.000 0.001 0.000 0.069 0.034 0.001 0.577 0.069 
   9 |  0.035 0.340 0.000 0.001 0.069 0.001 0.001 0.068 0.001 0.000 0.035 0.002 0.000 0.001 0.000 0.374 0.034 0.001 0.001 0.035 
  10 |  0.001 0.001 0.000 0.001 0.238 0.577 0.001 0.000 0.001 0.000 0.001 0.002 0.000 0.001 0.170 0.001 0.001 0.001 0.001 0.001 

  12 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.543 0.001 0.000 0.001 0.137 0.068 0.068 0.000 0.001 0.170 0.001 0.001 0.001 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



10 columns
Num Motifs: 29
   1,  1      79 agvda VIVLSANDPFVQ safgk      90   0.48 F 1091044
   2,  1      72 klstq ILAISVDSPFSH lqyll      83   1.00 F 11467494
   6,  1      66 nlntk IYAISNDSHFVQ knwie      77   1.00 F 13186328
   8,  1      68 kknte VISVSEDTVYVH kawvq      79   1.00 F 13541053
   9,  1      66 kfkak VIGISVDSPFSL aefak      77   1.00 F 13541117
  13,  1      65 eldce LVGLSVDQVFSH ikwie      76   0.98 F 14286173
  15,  1      65 klgav VIGVSTDSVEKN rkfae      76   1.00 F 14600438
  21,  1      80 megvd VTVVSMDLPFAQ krfce      91   0.65 F 15605963
  23,  1      70 tagln VVGISPDKPEKL atfrd      81   0.99 F 15609658
  24,  1      64 glntv ILGVSPDPVERH kkfie      75   1.00 F 15613511
  28,  1      56 fekaq VVGISRDSVEAL krfke      67   1.00 F 15643152
  30,  1     101 tpgla VWGISPDSTYAH eafad     112   0.98 F 15790738
  32,  1      80 idntv VLCISADLPFAQ srfcg      91   0.99 F 15801846
  35,  1      66 eagav VLGINRDSVYAH rawaa      77   1.00 F 15807234
  38,  1      66 evnav VIGISVDPPFSN kafke      77   1.00 F 15899339
  47,  1      71 krgve VVGVSFDSEFVH nawrn      82   1.00 F 16501671
  53,  1      73 ngvde IVCISVNDAFVM newak      84   0.32 F 17229033
  60,  1      68 kanae VLAISVDSPFAL kafkd      79   1.00 F 18313548
  62,  1      66 rrnav VLLISCDSVYTH kawas      77   1.00 F 19173077
  67,  1      67 klgvd VYSVSTDTHFTH kawhs      78   1.00 F 20151112
  70,  1      72 drdaq ILGFSGDSEFVH hawrk      83   1.00 F 21223405
  72,  1      78 keegi VLTISADLPFAQ krwca      89   0.99 F 21283385
  73,  1      65 srgit VIGISGDSPESH kqfae      76   1.00 F 21674812
  81,  1      68 krnvk LIALSIDSVEDH lawsk      79   1.00 F 3318841
  85,  1      50 eince VIGVSVDSVYCH qawce      61   1.00 F 4433065
  86,  1      74 kgide IICFSVNDPFVM kawgk      85   0.43 F 4704732
  87,  1      68 klnck LIGFSCNSKESH dqwie      79   0.56 F 4996210
  92,  1      68 klgve VLSVSVDSVFVH kmwnd      79   1.00 F 6850955
  95,  1      71 klgce VLGVSVDSQFTH lawin      82   1.00 F 9955016
                       ***** **** *


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif a = -469.15170
Log Fragmentation portion of MAP for motif a = -3.80666

-------------------------------------------------------------------------
                          MOTIF b





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    50  .   .   .   .   .   .   .   16  .   .   .   .   .   .   .   .   .   33  .   1.9
   2 |    .   .   .   .   .   16  .   .   .   .   .   50  .   .   .   .   33  .   .   .   2.1
   3 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   33  .   66  .   .   .   3.4
   4 |    .   33  .   .   .   16  .   .   .   .   .   33  .   .   16  .   .   .   .   .   1.7

   6 |    33  .   .   16  33  .   .   .   .   .   .   .   .   16  .   .   .   .   .   .   1.5

   8 |    .   .   .   .   .   .   .   50  .   .   16  .   .   .   .   33  .   .   .   .   3.2
   9 |    .   .   .   .   .   .   66  .   .   .   .   .   .   .   .   16  .   16  .   .   2.3
  10 |    .   33  .   16  33  .   .   .   .   .   .   .   .   .   16  .   .   .   .   .   1.7
  11 |    50  50  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   2.1
  12 |    .   .   66  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   33  4.4
  13 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  100  .   .   .   .   3.8
  14 |    50  50  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   2.1

  16 |    .   .   .   .   .   .   .   .   .  100  .   .   .   .   .   .   .   .   .   .   5.8
  17 |    .   .   .   .   16  .   .   .   .   .   66  .   .   16  .   .   .   .   .   .   2.1

  19 |    .   .   .   .   .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   3.2

nonsite    8   7  .    6   7   4   7   1   6  .    7   9   2   4   2   4   3   4   5   5
site      12  11   4   2   5   2  11   3   1   6   5   5  .    2   4  10   6   1   2   2

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.468 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.158 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.312 0.004 
   2 |  0.006 0.006 0.001 0.005 0.005 0.158 0.005 0.001 0.004 0.001 0.006 0.469 0.002 0.003 0.002 0.004 0.310 0.003 0.004 0.004 
   3 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.310 0.004 0.618 0.003 0.004 0.004 
   4 |  0.006 0.314 0.001 0.005 0.005 0.158 0.005 0.001 0.004 0.001 0.006 0.315 0.002 0.003 0.156 0.004 0.002 0.003 0.004 0.004 

   6 |  0.314 0.006 0.001 0.159 0.313 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.157 0.002 0.004 0.002 0.003 0.004 0.004 

   8 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.463 0.004 0.001 0.159 0.007 0.002 0.003 0.002 0.312 0.002 0.003 0.004 0.004 
   9 |  0.006 0.006 0.001 0.005 0.005 0.004 0.621 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.158 0.002 0.157 0.004 0.004 
  10 |  0.006 0.314 0.001 0.159 0.313 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.156 0.004 0.002 0.003 0.004 0.004 
  11 |  0.468 0.468 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 
  12 |  0.006 0.006 0.617 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.312 
  13 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.927 0.002 0.003 0.004 0.004 
  14 |  0.468 0.468 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 

  16 |  0.006 0.006 0.001 0.005 0.005 0.004 0.005 0.001 0.004 0.924 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 
  17 |  0.006 0.006 0.001 0.005 0.159 0.004 0.005 0.001 0.004 0.001 0.621 0.007 0.002 0.157 0.002 0.004 0.002 0.003 0.004 0.004 

  19 |  0.006 0.006 0.001 0.005 0.005 0.004 0.928 0.001 0.004 0.001 0.006 0.007 0.002 0.003 0.002 0.004 0.002 0.003 0.004 0.004 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



15 columns
Num Motifs: 6
   2,  1     161 riles IQYVKENPGYACPVNWNFG dqvfy     179   1.00 F 11467494
  47,  1     160 lrmvd ALQFHEEHGDVCPAQWEKG kegmn     178   1.00 F 16501671
  67,  1     154 rkika AQYVAAHPGEVCPAKWKEG eatla     172   1.00 F 20151112
  81,  1     166 lrvvi SLQLTAEKRVATPVDWKDG dsvmv     184   1.00 F 3318841
  87,  1     163 lrvlk SLQLTNTHPVATPVNWKEG dkcci     181   1.00 F 4996210
  95,  1     160 lrlvq AFQYTDEHGEVCPAGWKPG sdtik     178   1.00 F 9955016
                       **** * ******* ** *


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif b = -187.76179
Log Fragmentation portion of MAP for motif b = -7.77486

-------------------------------------------------------------------------
                          MOTIF c





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |    .   .   .    5   3  .    5  .   .   .   53   3  .   .   .    3  11   7   1   3  1.6

   3 |    .   61  .   .   .   .   .   .   22  .   .    7   3  .   .   .   .   .    1   3  2.1
   4 |    .   38  .   .   .    3   3  .   11  .   .   40   1  .   .   .   .   .   .   .   1.7
   5 |     5  35  .   .   .   .   .   .   24  .    1  31   1  .   .   .   .   .   .   .   1.7
   6 |    .   .   .   48  .   .   .    1  .   .   .   .   .   37  11  .    1  .   .   .   2.7
   7 |     1   7  .   .   .   85  .   .   .   .   .    3  .   .    1  .   .   .   .   .   3.4
   8 |    .   .   .   .   .    5   3   1  .   55  .   .   .   .   18  .   .   .    9   5  3.5
   9 |    87  .   .   .   .    1   5  .   .   .   .   .   .   .   .   .   .   .    3   1  2.8
  10 |     3  .   .   14   9  .   .    1  .   .   .   .   .    1  .   16   5  .   20  25  1.5
  11 |    .   .   .   .   .   .    1  .   .   90  .    1  .   .    1  .   .   .    1   1  5.2
  12 |    .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   6.0
  13 |     9   7  .   .    1  .   53  .   .   .    1  .    1  .   .   24  .   .   .   .   2.0
  14 |    .    7  .    1  .   .    1  .   .   .   .    1  .   .    1  83  .   .    1  .   3.2
  15 |    .   .  100  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   6.0
  16 |    .    5  .    1   1  .   .    3   1  .   27   1  .   .   .   .    9  46  .   .   2.1

  18 |    .    3  .   .   37  12  .   .   18  .   .   16   3  .   .   .    3  .   .    3  1.4

  20 |    .   .   .    1  .   .   .   .   .   .    3  .   .   .   .   94  .   .   .   .   3.9

  22 |    .    7  .   .   .   22  .   .    9  .   .   44  11  .    5  .   .   .   .   .   1.8

  24 |     7  .   .    3  37  .   .    3  .   .   27   1  .    1  .   .   11   1   1   1  1.4
  25 |     3  18  .    1  .    9  .   .    7  .   .   51   1  .   .   .   .   .    1   3  1.4

nonsite    8   7   1   6   7   4   6   1   5   1   7   9   2   4   2   4   3   4   4   4
site       5   9  10   3   4   7   3  .    4   7   5  10   1   2   2  11   2   2   2   2

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.001 0.001 0.000 0.056 0.037 0.000 0.056 0.000 0.001 0.000 0.533 0.038 0.000 0.000 0.000 0.037 0.110 0.074 0.019 0.037 

   3 |  0.001 0.606 0.000 0.001 0.001 0.000 0.001 0.000 0.221 0.000 0.001 0.074 0.037 0.000 0.000 0.000 0.000 0.000 0.019 0.037 
   4 |  0.001 0.386 0.000 0.001 0.001 0.037 0.037 0.000 0.111 0.000 0.001 0.405 0.019 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
   5 |  0.056 0.349 0.000 0.001 0.001 0.000 0.001 0.000 0.239 0.000 0.019 0.313 0.019 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
   6 |  0.001 0.001 0.000 0.478 0.001 0.000 0.001 0.018 0.001 0.000 0.001 0.001 0.000 0.367 0.110 0.000 0.019 0.000 0.000 0.000 
   7 |  0.019 0.074 0.000 0.001 0.001 0.845 0.001 0.000 0.001 0.000 0.001 0.038 0.000 0.000 0.019 0.000 0.000 0.000 0.000 0.000 
   8 |  0.001 0.001 0.000 0.001 0.001 0.056 0.037 0.018 0.001 0.551 0.001 0.001 0.000 0.000 0.184 0.000 0.000 0.000 0.092 0.056 
   9 |  0.863 0.001 0.000 0.001 0.001 0.019 0.056 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.037 0.019 
  10 |  0.037 0.001 0.000 0.147 0.092 0.000 0.001 0.018 0.001 0.000 0.001 0.001 0.000 0.019 0.000 0.166 0.055 0.000 0.202 0.257 
  11 |  0.001 0.001 0.000 0.001 0.001 0.000 0.019 0.000 0.001 0.899 0.001 0.019 0.000 0.000 0.019 0.000 0.000 0.000 0.019 0.019 
  12 |  0.001 0.001 0.991 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
  13 |  0.093 0.074 0.000 0.001 0.019 0.000 0.533 0.000 0.001 0.000 0.019 0.001 0.019 0.000 0.000 0.239 0.000 0.000 0.000 0.000 
  14 |  0.001 0.074 0.000 0.019 0.001 0.000 0.019 0.000 0.001 0.000 0.001 0.019 0.000 0.000 0.019 0.826 0.000 0.000 0.019 0.000 
  15 |  0.001 0.001 0.991 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 
  16 |  0.001 0.056 0.000 0.019 0.019 0.000 0.001 0.037 0.019 0.000 0.276 0.019 0.000 0.000 0.000 0.000 0.092 0.459 0.000 0.000 

  18 |  0.001 0.037 0.000 0.001 0.368 0.129 0.001 0.000 0.184 0.000 0.001 0.166 0.037 0.000 0.000 0.000 0.037 0.000 0.000 0.037 

  20 |  0.001 0.001 0.000 0.019 0.001 0.000 0.001 0.000 0.001 0.000 0.037 0.001 0.000 0.000 0.000 0.936 0.000 0.000 0.000 0.000 

  22 |  0.001 0.074 0.000 0.001 0.001 0.221 0.001 0.000 0.092 0.000 0.001 0.441 0.110 0.000 0.055 0.000 0.000 0.000 0.000 0.000 

  24 |  0.074 0.001 0.000 0.037 0.368 0.000 0.001 0.037 0.001 0.000 0.276 0.019 0.000 0.019 0.000 0.000 0.110 0.019 0.019 0.019 
  25 |  0.037 0.184 0.000 0.019 0.001 0.092 0.001 0.000 0.074 0.000 0.001 0.515 0.019 0.000 0.000 0.000 0.000 0.000 0.019 0.037 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



20 columns
Num Motifs: 54
   3,  1      17 viqsd KLVVVDFYADWCMPCRYISPILEKL skeyn      41   1.00 F 11499727
   4,  1      22 llntt QYVVADFYADWCGPCKAIAPMYAQF aktfs      46   1.00 F 1174686
   5,  1      19 ifsak KNVIVDFWAAWCGPCKLTSPEFQKA adefs      43   1.00 F 12044976
   7,  1      21 kenhs KPILIDFYADWCPPCRMLIPVLDSI ekkhg      45   1.00 F 13358154
  10,  1      17 aetse GVVLADFWAPWCGPCKMIAPVLEEL dqemg      41   1.00 F 135765
  11,  1      29 akesn KLIVIDFTASWCPPCRMIAPIFNDL akkfm      53   1.00 F 1388082
  12,  1      44 likqn DKLVIDFYATWCGPCKMMQPHLTKL iqayp      68   1.00 F 140543
  14,  1      80 selrg KVVMLQFTASWCGVCRKEMPFIEKD iwlkh     104   1.00 F 14578634
  16,  1      23 lqnsd KPVLVDFYATWCGPCQLMVPILNEV setlk      47   1.00 F 15218394
  17,  1     157 adfrg RPLVINLWASWCPPCRREMPVLQQA qaenp     181   1.00 F 15597673
  18,  1      26 ensfh KPVLVDFWADWCAPCKALMPLLAQI aesyq      50   1.00 F 15599256
  19,  1      67 negkg KTILLNFWSETCGVCIAELKTFEQL lqsyp      91   1.00 F 15602312
  20,  1      61 eefkg KVLLINFWATWCPPCKEEIPMFKEI yekyr      85   1.00 F 15605725
  25,  1      60 sdyrg DVVILNVWASWCEPCRKEMPALMEL qsdye      84   1.00 F 15614085
  26,  1      63 releg KGVFLNFWGTYCPPCEREMPHMEKL ygeyk      87   1.00 F 15614140
  27,  1      72 sslrg QPVILHFFATWCPVCQDEMPSLVKL dkeyr      96   1.00 F 15615431
  31,  1       2     m TVTLKDFYADWCGPCKTQDPILEEL eadyd      26   1.00 F 15791337
  33,  1      72 adyrg RPVVLNFWASWCGPCREEAPLFAKL aahpg      96   1.00 F 15805225
  34,  1      78 taaqg KPVVINFWASWCVPCRQEAPLFSKL sqeta     102   1.00 F 15805374
  37,  1      49 fitkn KIVVVDFWAEWCAPCLILAPVIEEL andyp      73   1.00 F 15899007
  40,  1      61 easrq QPVLVDFWAPWCGPCKQLTPVIEKV vreaa      85   1.00 F 15966937
  41,  1      61 sdfrg KTLLVNLWATWCVPCRKEMPALDEL qgkls      85   1.00 F 15988313
  42,  1      60 qdakg KKVLLNFWATWCKPCRQEMPAMEKL qkeya      84   1.00 F 16078864
  43,  1      53 llqdd LPMVIDFWAPWCGPCRSFAPIFAET aaera      77   1.00 F 16123427
  46,  1      21 vlkad GAILVDFWAEWCGPCKMIAPILDEI adeyq      45   1.00 F 1633495
  48,  1      34 vlqcp KPILVYFGAPWCGLCHFVKPLLNHL hgewq      58   1.00 F 1651717
  49,  1      60 tlsee RPVLLYFWASWCGVCRFTTPAVAHL aaege      84   1.00 F 16759994
  50,  1      53 llkdd LPVVIDFWAPWCGPCRNFAPIFEDV aeers      77   1.00 F 16761507
  52,  1      19 iissh PKILLNFWAEWCAPCRCFWPTLEQF aemee      43   1.00 F 16804867
  54,  1      22 vlsed KVVVVDFTATWCGPCRLVSPLMDQL adeyk      46   1.00 F 17229859
  55,  1      18 vlegt GYVLVDYFSDGCVPCKALMPAVEEL skkye      42   1.00 F 1729944
  56,  1      28 rqhpe KIIILDFYATWCGPCKAIAPLYKEL atthk      52   1.00 F 17531233
  57,  1      27 ehlkg KIIGLYFSASWCPPCRAFTPKLKEF feeik      51   1.00 F 17537401
  58,  1      63 safrg QPVVINFWAPWCGPCVEEMPELSAL aqeqk      87   1.00 F 17547503
  59,  1     286 seykg KTIFLNFWATWCPPCRGEMPYIDEL ykeyn     310   1.00 F 18309723
  61,  1      44 dsllg KKIGLYFSAAWCGPCQRFTPQLVEV ynels      68   1.00 F 18406743
  61,  2     364 sdlvg KTILMYFSAHWCPPCRAFTPKLVEV ykqik     388   1.00 F 18406743
  63,  1      15 sdfeg EVVVLNAWGQWCAPCRAEVDDLQLV qetld      39   1.00 F 19554157
  64,  1      39 eeykg KVVVINFWATWCGYCVEEMPGFEKV ykefg      63   1.00 F 19705357
  66,  1       7 agdfm KPMLLDFSATWCGPCRMQKPILEEL ekkyg      31   1.00 F 20092028
  69,  1     103 adykg KVVVLNVWGSWCPPCRAEAKNFEKV yqdvk     127   1.00 F 21222859
  74,  1      53 sdfkg ERVLINFWTTWCPPCRQEMPDMQRF yqdlq      77   1.00 F 23098307
  76,  1      20 kylqh QRVVVDFSAEWCGPCRAIAPVFDKL sneft      44   1.00 F 267116
  77,  1      81 aafkg KVSLVNVWASWCVPCHDEAPLLTEL gkdkr     105   1.00 F 27375582
  78,  1      34 vtsdn DVVLADFYADWCGPCQMLEPVVETL aeqtd      58   1.00 F 2822332
  79,  1      77 sdlkg KKVILNFWATWCGPCQQEMPDMEAF ykehk     101   1.00 F 30021713
  80,  1      19 tisan SNVLVYFWAPLCAPCDLFTPTYEAS srkhf      43   1.00 F 3261501
  82,  1      19 tietn PLVIVDFWAPWCGSCKMLGPVLEEV esevg      43   1.00 F 3323237
  83,  1      17 ektah QAVVVNVGASWCPDCRKIEPIMENL aktyk      41   1.00 F 4155972
  84,  1      79 vvnse TPVVVDFHAQWCGPCKILGPRLEKM vakqh     103   1.00 F 4200327
  91,  1      20 nenkg RLIVVDFFAQWCGPCRNIAPKVEAL akeip      44   1.00 F 6687568
  93,  1      18 llttn KKVVVDFYANWCGPCKILGPIFEEV aqdkk      42   1.00 F 7109697
  94,  1      21 ilaed KLVVIDFYADWCGPCKIIAPKLDEL aqqys      45   1.00 F 7290567
  96,  1      49 adlqg KVTLINFWFPSCPGCVSEMPKIIKT andyk      73   1.00 F 15677788
                       * ************** * * * **


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif c = -1607.59351
Log Fragmentation portion of MAP for motif c = -10.42374

-------------------------------------------------------------------------
                          MOTIF d





Motif model (residue frequency x 100)
____________________________________________
Pos. #     a   v   c   d   e   f   g   h   i   w   k   l   m   n   y   p   q   r   s   t  Info
_____________________________________________________________________________________________
   1 |     2  .   .   28  17   2  .    2   8  .   .   17  .   .   11  .   .   .    2   5  1.1
   2 |     5   2  .   .    5  34  .   .   .   .    2  14  .   .   17  .   .    8   2   5  1.4
   3 |     8  .   .    5  11  .   14  .    2  .   28  .   .    5   2  .   .   17  .    2  1.0
   4 |     2  .   .   11  .   .   62  .   .   .    5  .   .   11  .   .    2  .    2  .   2.1
   5 |    .   .   .   .   .   .    2  .   .   .   57  .   .    5  .   .    5  22   5  .   2.2

   7 |     2  68  .   .   .    2   5  .    2  .    2  .   .    2  .   .   .   .    2   8  1.9
   8 |     2  62  .   .   .   .   .   .   25  .   .    8  .   .   .   .   .   .   .   .   2.3
   9 |    .    8  .   .   .    8  .   .    5  .   .   77  .   .   .   .   .   .   .   .   2.3
  10 |     8   8  .   .   .   57  .   .    2  .   .    5  .   .   11  .   .   .    2   2  2.1
  11 |    14   2  .   .   .   62   8  .   .   .   .   .   .   .   .   .   .   .   11  .   2.4
  12 |     2  11  .   .   .   14  .   11   2   5  .    5  .   .   45  .   .   .   .   .   2.4
  13 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  100  .   .   .   .   4.3
  14 |    22  .   .   .   .    2  28   2  .   .    8  17   5  .    2  .   .    8  .   .   1.2
  15 |    51  .   .   42  .   .   .   .   .   .   .   .   .    2  .   .   .   .    2  .   2.3
  16 |    .   .   .    2  .   71   2  .   .    5  .   .    5  .    5  .   .   .    5  .   2.9
  17 |    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   11  88  3.6
  18 |     5  .   .   .   .   25   2  .   .   .   .   .   .   .   .   51   2  .   11  .   2.4
  19 |    .   54  .   .   .   .   20  .    8  .   .   .   .   .   .   .   .   .   .   17  2.0
  20 |    .   .   97  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .    2  .   6.2
  21 |     5  .   .   .   .   .   .   .   .   .   .   .   .   .   .   28   2  .   20  42  2.3
  22 |    14   2  .   .   .   .    2  .   .   .   14   5   5  .   .   .    2   8   5  37  1.3
  23 |    .   .   .   .   62  .   .   .   .   .    2  .   .    8  .   .   14  .    5   5  2.2
  24 |    14   2  .   .   .    2   2  22  11  .   .   31  11  .   .   .   .   .   .   .   1.8

  27 |     2   2  .   .    2  45  17  .    8  .   .    8  .   .    8   2  .   .   .   .   1.6
  28 |    11  .   .    2   2   8   5  .   .   .   .   .   .    2  14  .    5  22  22  .   1.4

nonsite    8   7  .    6   7   4   7   1   5  .    7   9   2   4   2   4   3   4   5   5
site       7   9   3   3   4  13   7   1   3  .    4   7   1   1   4   7   1   3   4   8

Motif probability model
____________________________________________
Pos. #    a     v     c     d     e     f     g     h     i     w     k     l     m     n     y     p     q     r     s     t   
____________________________________________
   1 |  0.029 0.001 0.000 0.283 0.170 0.029 0.001 0.028 0.085 0.000 0.001 0.170 0.000 0.001 0.113 0.001 0.000 0.001 0.029 0.057 
   2 |  0.058 0.029 0.000 0.001 0.057 0.339 0.001 0.000 0.001 0.000 0.029 0.142 0.000 0.001 0.169 0.001 0.000 0.085 0.029 0.057 
   3 |  0.086 0.001 0.000 0.057 0.114 0.001 0.142 0.000 0.029 0.000 0.283 0.001 0.000 0.057 0.029 0.001 0.000 0.170 0.001 0.029 
   4 |  0.029 0.001 0.000 0.114 0.001 0.001 0.621 0.000 0.001 0.000 0.057 0.001 0.000 0.113 0.000 0.001 0.029 0.001 0.029 0.001 
   5 |  0.001 0.001 0.000 0.001 0.001 0.001 0.029 0.000 0.001 0.000 0.564 0.001 0.000 0.057 0.000 0.001 0.057 0.226 0.057 0.001 

   7 |  0.029 0.677 0.000 0.001 0.001 0.029 0.057 0.000 0.029 0.000 0.029 0.001 0.000 0.029 0.000 0.001 0.000 0.001 0.029 0.085 
   8 |  0.029 0.621 0.000 0.001 0.001 0.001 0.001 0.000 0.254 0.000 0.001 0.086 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.001 
   9 |  0.001 0.086 0.000 0.001 0.001 0.085 0.001 0.000 0.057 0.000 0.001 0.762 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.001 
  10 |  0.086 0.086 0.000 0.001 0.001 0.564 0.001 0.000 0.029 0.000 0.001 0.058 0.000 0.001 0.113 0.001 0.000 0.001 0.029 0.029 
  11 |  0.142 0.029 0.000 0.001 0.001 0.620 0.085 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.113 0.001 
  12 |  0.029 0.114 0.000 0.001 0.001 0.142 0.001 0.113 0.029 0.057 0.001 0.058 0.000 0.001 0.451 0.001 0.000 0.001 0.001 0.001 
  13 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.987 0.000 0.001 0.001 0.001 
  14 |  0.227 0.001 0.000 0.001 0.001 0.029 0.283 0.028 0.001 0.000 0.086 0.170 0.057 0.001 0.029 0.001 0.000 0.085 0.001 0.001 
  15 |  0.508 0.001 0.000 0.423 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.029 0.000 0.001 0.000 0.001 0.029 0.001 
  16 |  0.001 0.001 0.000 0.029 0.001 0.705 0.029 0.000 0.001 0.057 0.001 0.001 0.057 0.001 0.057 0.001 0.000 0.001 0.057 0.001 
  17 |  0.001 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.113 0.874 
  18 |  0.058 0.001 0.000 0.001 0.001 0.254 0.029 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.508 0.029 0.001 0.113 0.001 
  19 |  0.001 0.536 0.000 0.001 0.001 0.001 0.198 0.000 0.085 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.001 0.170 
  20 |  0.001 0.001 0.958 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.001 0.000 0.001 0.029 0.001 
  21 |  0.058 0.001 0.000 0.001 0.001 0.001 0.001 0.000 0.001 0.000 0.001 0.001 0.000 0.001 0.000 0.282 0.029 0.001 0.198 0.423 
  22 |  0.142 0.029 0.000 0.001 0.001 0.001 0.029 0.000 0.001 0.000 0.142 0.058 0.057 0.001 0.000 0.001 0.029 0.085 0.057 0.367 
  23 |  0.001 0.001 0.000 0.001 0.621 0.001 0.001 0.000 0.001 0.000 0.029 0.001 0.000 0.085 0.000 0.001 0.141 0.001 0.057 0.057 
  24 |  0.142 0.029 0.000 0.001 0.001 0.029 0.029 0.226 0.113 0.000 0.001 0.311 0.113 0.001 0.000 0.001 0.000 0.001 0.001 0.001 

  27 |  0.029 0.029 0.000 0.001 0.029 0.451 0.170 0.000 0.085 0.000 0.001 0.086 0.000 0.001 0.085 0.029 0.000 0.001 0.001 0.001 
  28 |  0.114 0.001 0.000 0.029 0.029 0.085 0.057 0.000 0.001 0.000 0.001 0.001 0.000 0.029 0.141 0.001 0.057 0.226 0.226 0.001 



Background probability model
        0.089 0.079 0.008 0.067 0.076 0.044 0.071 0.013 0.061 0.009 0.076 0.094 0.023 0.043 0.027 0.045 0.034 0.044 0.052 0.052 



25 columns
Num Motifs: 35
   1,  1      37 ldfdk EFRDKTVVIVAIPGAFTPTCTANHIPPF vekft      64   1.00 F 1091044
   2,  1      32 irlsd YRGKKYVILFFYPANFTAISPTELMLLS drise      59   1.00 F 11467494
   6,  1      26 eikei DLKSNWNVFFFYPYSYSFICPLELKNIS nkike      53   0.98 F 13186328
   8,  1      28 kirls SYRGKWVVLFFYPADFTFVCPTEVEGFA edyek      55   1.00 F 13541053
   9,  1      26 mrkls EFRGQNVVLAFFPGAFTSVCTKEMCTFR dsman      53   1.00 F 13541117
  13,  1      25 melpd EFEGKWFILFSHPADFTPVCTTEFVAFQ evype      52   1.00 F 14286173
  15,  1      25 kirls DFRGRIVVLYFYPRAMTPGCTREGVRFN ellde      52   1.00 F 14600438
  22,  1      26 vtlrg YRGAKNVLLVFFPLAFTGICQGELDQLR dhlpe      53   1.00 F 15609375
  23,  1      30 nvsla DYRGRRVIVYFYPAASTPGCTKQACDFR dnlgd      57   1.00 F 15609658
  24,  1      24 tvsls DFKGKNIVLYFYPKDMTPGCTTEACDFR drved      51   1.00 F 15613511
  28,  1      20 tfthv DLYGKYTILFFFPKAGTSGCTREAVEFS renfe      47   1.00 F 15643152
  30,  1      61 gltda LADNRAVVLFFYPFDFSPVCATELCAIQ narwf      88   1.00 F 15790738
  35,  1      26 itlss YRGQSHVVLVFYPLDFSPVCSMQLPEYS gsqdd      53   1.00 F 15807234
  36,  1      28 vnlae LFKGKKGVLFGVPGAFTPGCSKTHLPGF veqae      55   1.00 F 15826629
  38,  1      26 vkips DFKGKVVVLAFYPAAFTSVCTKEMCTFR dsmak      53   1.00 F 15899339
  39,  1      30 vttel LFKGKRVVLFAVPGAFTPTCSLNHLPGY lenrd      57   1.00 F 15964668
  44,  1      50 fnlak ALKKGPVVLYFFPAAYTAGCTAEAREFA eatpe      77   1.00 F 16125919
  47,  1      31 fnfkq HTNGKTTVLFFWPMDFTFVCPSELIAFD kryee      58   1.00 F 16501671
  51,  1      33 slekn IEDDKWTILFFYPMDFTFVCPTEIVAIS arsde      60   1.00 F 16803644
  53,  1      31 vttdd LFAGKTVAVFSLPGAFTPTCSSTHLPGY nelak      58   1.00 F 17229033
  60,  1      28 rlsev LKRGRPVVLLFFPGAFTSVCTKELCTFR dkmal      55   1.00 F 18313548
  62,  1      26 eislq DYIGKYVVLAFYPLDFTFVCPTEINRFS dlkga      53   1.00 F 19173077
  67,  1      27 evtek DTEGRWSVFFFYPADFTFVCPTELGDVA dhyee      54   1.00 F 20151112
  68,  1      29 vdtht LFTGRKVVLFAVPGAFTPTCSAKHLPGY veqfe      56   1.00 F 21112072
  70,  1      32 qinhk TYEGQWKVVFAWPKDFTFVCPTEIAAFG klnde      59   1.00 F 21223405
  71,  1      28 eihly DLKGKKVLLSFHPLAWTQVCAQQMKSLE enyel      55   1.00 F 21227878
  73,  1      25 mvsls EFKGRKVLLIFYPGDDTPVCTAQLCDYR nnvaa      52   1.00 F 21674812
  81,  1      28 irfhd FLGDSWGILFSHPRDFTPVCTTELGRAA klape      55   1.00 F 3318841
  85,  1      10 eidin EYKGKYVVLLFYPLDWTFVCPTEMIGYS evagq      37   1.00 F 4433065
  86,  1      32 vsvhs IAAGKKVILFGVPGAFTPTCSMSHVPGF igkae      59   1.00 F 4704732
  87,  1      28 fdfyk YVGDNWAILFSHPHDFTPVCTTELAEFG kmhee      55   1.00 F 4996210
  88,  1      41 ynask EFANKKVVLFALPGAFTPVCSANHVPEY iqklp      68   1.00 F 5326864
  89,  1      88 slkki TENNRVVVFFVYPRASTPGCTRQACGFR dnyqe     115   1.00 F 6322180
  90,  1      43 ewskl ISENKKVIITGAPAAFSPTCTVSHIPGY inyld      70   0.99 F 6323138
  95,  1      31 evkls DYKGKYVVLFFYPLDFTFVCPTEIIAFS nraed      58   1.00 F 9955016
                       ***** ******************  **


Column 1 :  Sequence Number, Site Number
Column 2 :  Left End Location
Column 4 :  Motif Element
Column 5 :  Right End Location
Column 6 :  Probability of Element
Column 7 :  Forward Motif (F) or Reverse Complement (R) 
Column 8 :  Sequence Description from Fast A input

Log Motif portion of MAP for motif d = -1668.31468
Log Fragmentation portion of MAP for motif d = -7.86327


Log Background portion of Map = -39912.17887
Log Alignment portion of Map = -956.36102
Log Site/seq portion of Map = 0.00000
Log Null Map = -46943.36311
Log Map = 2112.13301


log MAP = sum of motif and fragmentation parts of MAP + background + alignment + sites/seq - Null

Frequency Map = 2109.909622
Nearopt Map   = 2111.157969
Maximal Map   = 2111.157969
Total Time 105 sec (1.750000 min)
Elapsed time: 104.960000 secs
DOF[0] = 190
DOF[1] = 285
DOF[2] = 380
DOF[3] = 475
"""
        
#run if called from command-line
if __name__ == "__main__":
    main()
