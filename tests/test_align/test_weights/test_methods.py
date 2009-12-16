#!/usr/bin/env python

from __future__ import division
from numpy import array, zeros, float64 as Float64
from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser
from cogent.parse.clustal import ClustalParser as MinimalClustalParser
from cogent.core.alignment import Alignment
from cogent.core.profile import Profile
from cogent.align.weights.util import DNA_ORDER, PROTEIN_ORDER
from cogent.align.weights.methods import VA, VOR, mVOR, pos_char_weights, PB,\
    SS, ACL, GSC, _clip_branch_lengths, _set_branch_sum, _set_node_weight

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

def ClustalParser(f):
    return Alignment(list(MinimalClustalParser(f)))

class GeneralTests(TestCase):
    """General tests for all classes in this file, provides general setup"""

    def setUp(self):
        """General setUp method for all tests in this file"""
        
        #ALIGNMENTS
        self.aln1 = Alignment(['ABC','BCC','BAC'])
        
        #alignment from Henikoff 1994
        self.aln2 = Alignment({'seq1':'GYVGS','seq2':'GFDGF','seq3':'GYDGF',\
            'seq4':'GYQGG'},Names=['seq1','seq2','seq3','seq4'])

        #alignment from Vingron & Sibbald 1993
        self.aln3 = Alignment({'seq1':'AA', 'seq2':'AA', 'seq3':'BB'},\
            Names=['seq1','seq2','seq3'])

        #alignment from Vingron & Sibbald 1993
        self.aln4 = Alignment({'seq1':'AA', 'seq2':'AA', 'seq3':'BB',\
        'seq4':'BB','seq5':'CC'},Names=['seq1','seq2','seq3','seq4','seq5'])

        self.aln5 = Alignment(['ABBA','ABCA','CBCB'])

        #alignment 5S rRNA seqs from Hein 1990
        self.aln6 = ClustalParser(FIVE_S_ALN.split('\n'))

        #alignment from Vingron & Sibbald 1993
        self.aln7 = Alignment({'seq1':'AGCTA', 'seq2':'AGGTA', 'seq3':'ACCTG',
            'seq4':'TGCAA'},Names=['seq1','seq2','seq3','seq4'])


        #TREES (SEE BOTTOM OF FILE FOR DESCRIPTION)
        self.tree1 = DndParser(TREE_1)
        self.tree2 = DndParser(TREE_2)
        self.tree3 = DndParser(TREE_3)
        self.tree4 = DndParser(TREE_4)
        self.tree5 = DndParser(TREE_5)
        self.tree6 = DndParser(TREE_6)
        self.tree7 = DndParser(TREE_7)
        self.tree8 = DndParser(TREE_8)
        self.tree9 = DndParser(TREE_9)


class VoronoiTests(GeneralTests):
    """Tests for the voronoi sequence weighting module"""


    def test_VA(self):
        """VA: should return expected results. Results don't vary with runs"""
        err=1e-3
        aln2_exp = {'seq1':.269, 'seq2':.269,'seq3':.192,'seq4':.269}
        aln4_exp = {'seq1':.1875, 'seq2':.1875,'seq3':.1875,'seq4':.1875,\
            'seq5':.25}
        aln5_exp = {'seq_0':0.33333,'seq_1':0.25,'seq_2':0.4167}
        aln6_exp = dict(zip(map(str,[1,2,3,4,5,6,7,8,9,10]),\
            [0.0962,0.0925,0.1061,0.1007,0.0958,0.0977,0.0914,\
            0.0934,0.1106,0.1156]))
 
        self.assertFloatEqualAbs(VA(self.aln2).values(),aln2_exp.values(),
            eps=err)
        self.assertFloatEqualAbs(VA(self.aln4).values(),aln4_exp.values(),
            eps=err)
        self.assertFloatEqualAbs(VA(self.aln5).values(),aln5_exp.values(),
            eps=err)
        self.assertFloatEqualAbs(VA(self.aln6).values(),aln6_exp.values(),
            eps=err)

        results = []
        for x in range(5):
            results.append(VA(self.aln2))
            if x > 0:
                self.assertEqual(results[x], results[x-1])

    def test_VOR_exact(self):
        """VOR: should give exact results when using pseudo_seqs_exact"""
        err=1e-3
        aln2_exp = {'seq1':.259, 'seq2':.315,'seq3':.167,'seq4':.259}
        aln3_exp = {'seq1':.29167, 'seq2':.29167, 'seq3':.4167}
        aln4_exp = {'seq1':.1851, 'seq2':.1851,'seq3':.1851,'seq4':.1851,\
            'seq5':.259}
        
        self.assertFloatEqualAbs(VOR(self.aln2).values(),aln2_exp.values(),
            eps=err)
        self.assertFloatEqualAbs(VOR(self.aln3).values(),aln3_exp.values(),\
            eps=err)
        self.assertFloatEqualAbs(VOR(self.aln4).values(),aln4_exp.values(),\
            eps=err)

        #this is the exact method, so the answer should be exactly the same
        #every time (on the same alignment)
        results = []
        for x in range(5):
            results.append(VOR(self.aln2))
            if x > 0:
                self.assertEqual(results[x], results[x-1])

    def test_VOR_force_mc(self):
        """VOR: should result in good approximation when using monte carlo"""
        err=9e-2
        aln2_exp = {'seq1':.259, 'seq2':.315,'seq3':.167,'seq4':.259}
        aln3_exp = {'seq1':.29167, 'seq2':.29167, 'seq3':.4167}
        aln4_exp = {'seq1':.1851, 'seq2':.1851,'seq3':.1851,'seq4':.1851,\
            'seq5':.259}
        aln6_exp = dict(zip(map(str,[1,2,3,4,5,6,7,8,9,10]),\
            [0.0840,0.0763,0.1155,0.1019,0.0932,0.0980,0.0864,\
            0.0999,0.1121,0.1328]))
 
        # the following assertSimilarMeans statements were added to replace 
        # stochastic assertFloatEqualAbs calls below
        self.assertSimilarMeans(VOR(self.aln2,force_monte_carlo=True).values(),
                                 aln2_exp.values())
        self.assertSimilarMeans(VOR(self.aln3,force_monte_carlo=True).values(),
                                 aln3_exp.values())
        self.assertSimilarMeans(VOR(self.aln4,force_monte_carlo=True).values(),
                                 aln4_exp.values())
        self.assertSimilarMeans(VOR(self.aln6,n=1000).values(),
                                 aln6_exp.values())

        #self.assertFloatEqualAbs(VOR(self.aln2,force_monte_carlo=True)\
        #    .values(), aln2_exp.values(),eps=err)
        #self.assertFloatEqualAbs(VOR(self.aln3,force_monte_carlo=True)\
        #    .values(), aln3_exp.values(),eps=err)
        #self.assertFloatEqualAbs(VOR(self.aln4,force_monte_carlo=True)\
        #    .values(), aln4_exp.values(),eps=err)
        #self.assertFloatEqualAbs(VOR(self.aln6,n=1000)\
        #    .values(), aln6_exp.values(),eps=err)
 
        #make sure monte carlo is used
        results = []
        for x in range(5):
            results.append(VOR(self.aln2,force_monte_carlo=True))
            if x > 0:
                self.assertNotEqual(results[x], results[x-1])

    def test_VOR_mc_threshold(self):
        """VOR: should apply monte carlo when # of pseudo seqs > mc_threshold
        """
        err=9e-2
        aln2_exp = {'seq1':.259, 'seq2':.315,'seq3':.167,'seq4':.259}

        # the following assertSimilarMeans statement was added to replace 
        # stochastic assertFloatEqualAbs call below
        self.assertSimilarMeans(VOR(self.aln2, mc_threshold=15).values(), 
                                 aln2_exp.values())
        #self.assertFloatEqual(VOR(self.aln2,mc_threshold=15).values(),\
        #    aln2_exp.values(),err)
        
        #make sure monte carlo is used
        results = []
        for x in range(5):
            results.append(VOR(self.aln2,mc_threshold=15))
            if x > 0:
                self.assertNotEqual(results[x], results[x-1])
    
    def test_mVOR(self):
        """mVOR: should return weights closer to the 'True' weights"""
        #err=5e-2 #original error value
        # Raised the error value to prevent occasional failure of the test.
        # The mVOR method takes a sample from a distribution and the outcome
        # will depend on this sample. Every now and then, one of the weights
        # was more than 0.05 away from the expected weight. Raised the 
        # allowed error value to prevent that. To use the method on real
        # data, a larger sample should be taken (e.g. 10000?), but increasing
        # the sample size here would make the test too slow.
        err=0.075
        aln3_exp = {'seq1':.25, 'seq2':.25, 'seq3':.5}
        aln4_exp = {'seq1':.1667, 'seq2':.1667,'seq3':.1667,'seq4':.1667,\
            'seq5':.3333}       
        aln6_exp = dict(zip(map(str,[1,2,3,4,5,6,7,8,9,10]),
            [0.09021,0.08039,0.113560,0.10399,0.092370,0.097130,
            0.09198,0.09538,0.10927,0.12572]))

        # the following assertSimilarMeans statements were added to replace 
        # stochastic assertFloatEqualAbs calls below
        self.assertSimilarMeans(mVOR(self.aln3,order="ABC").values(),
                                 aln3_exp.values())
        self.assertSimilarMeans(mVOR(self.aln4,order="ABC").values(),
                                 aln4_exp.values())
        self.assertSimilarMeans(mVOR(self.aln6,order=DNA_ORDER,n=3000)\
                                 .values(), aln6_exp.values())

        #self.assertFloatEqualAbs(mVOR(self.aln3,order="ABC").values(),\
        #    aln3_exp.values(),eps=err)
        #self.assertFloatEqualAbs(mVOR(self.aln4,order="ABC").values(),\
        #    aln4_exp.values(),eps=err)
        #self.assertFloatEqualAbs(mVOR(self.aln6,order=DNA_ORDER,n=3000)\
        #    .values(), aln6_exp.values(),eps=err)
        
        #the results vary with runs, because the sample of random profiles
        #is different each time
        results = []
        for x in range(5):
            results.append(mVOR(self.aln4,order="ABC"))
            if x > 0:
                self.assertNotEqual(results[x], results[x-1])

class PositionBasedTests(GeneralTests):
    """Contains tests for PB (=position-based) method"""
    
    def test_pos_char_weights(self):
        """pos_char_weights: should return correct contributions at each pos
        """
        #build expected profile
        exp_data = zeros([len(PROTEIN_ORDER),self.aln2.SeqLen],Float64)
        exp = [{'G':1/4},{'Y':1/6,'F':1/2},{'V':1/3,'D':1/6,'Q':1/3},
            {'G':1/4},{'G':1/3,'F':1/6,'S':1/3}]
        for pos, weights in enumerate(exp):
            for k,v in weights.items():
                exp_data[PROTEIN_ORDER.index(k),pos] = v
        exp_aln2 = Profile(exp_data,Alphabet=PROTEIN_ORDER)

        #check observed against expected
        self.assertEqual(pos_char_weights(self.aln2,PROTEIN_ORDER).Data,
            exp_aln2.Data)
    
    def test_PB(self):
        """PB: should return correct weights"""
        err=1e-3
        aln2_exp = {'seq3': 0.2, 'seq2': 0.267, 'seq1': 0.267, 'seq4': 0.267}

        self.assertFloatEqualAbs(PB(self.aln2,PROTEIN_ORDER)\
            .values(), aln2_exp.values(),eps=err)

class SsTests(GeneralTests):
    """Tests for SS function"""

    def test_SS(self):
        """SS: should return the correct weights"""
        err=1e-3
        aln4_exp = {'seq1':.1910, 'seq2':.1910,'seq3':.1910,'seq4':.1910,\
            'seq5':.2361}
        aln6_exp = dict(zip(map(str,[1,2,3,4,5,6,7,8,9,10]),
            [0.0977,0.0942,0.1045,0.0997,0.0968,0.0988,
            0.0929,0.0950,0.1076,0.1122]))
        aln7_exp = {'seq1':.1792, 'seq2':.2447,'seq3':.2880,'seq4':.2880}

        self.assertFloatEqualAbs(SS(self.aln4).values(),aln4_exp.values(),
            eps=err) 
        self.assertFloatEqualAbs(SS(self.aln6).values(),aln6_exp.values(),
            eps=err)
        self.assertFloatEqualAbs(SS(self.aln7).values(),aln7_exp.values(),
            eps=err)


class AclTests(GeneralTests):
    """Contains tests for ACL functionality"""

    
    def test_ACL(self):
        """ACL: should return correct weights"""
        err=1e-3
        tree1_exp = {'WMJ2': 0.035, 'HXB': 0.017, 'WMJ1': 0.039,\
        'BH10': 0.013, 'CDC': 0.048, 'BRU': 0.014, 'SF2': 0.050,\
        'BH8': 0.006, 'RF': 0.085, 'ELI': 0.068, 'PV22': 0.013,\
        'Z6': 0.129, 'MAL': 0.115, 'WMJ3': 0.035, 'Z3': 0.333}
        tree2_exp = {'agcta':0.7380,'aggta':0.0,'acctg':0.0,'tgcaa':0.2620}
        tree3_exp = {'agcta':0.2857,'aggta':0.0,'acctg':0.2857,'tgcaa':0.4286}
        tree4_exp = {'10': 0.1186, '1': 0.0627, '3': 0.1307, '2': 0.0627,
            '5': 0.0919, '4': 0.1307, '7': 0.0958, '6': 0.0919,
            '9': 0.1186, '8': 0.0958}
        tree9_exp = {'A':.25,'B':.25,'C':.25,'D':.25}
        
        self.assertFloatEqualAbs(ACL(self.tree1), tree1_exp, eps=err)
        self.assertFloatEqualAbs(ACL(self.tree2), tree2_exp, eps=err)
        self.assertFloatEqualAbs(ACL(self.tree3), tree3_exp, eps=err)
        self.assertFloatEqualAbs(ACL(self.tree4), tree4_exp, eps=err)
        #also works when branch lengths are zero
        self.assertFloatEqualAbs(ACL(self.tree9), tree9_exp, eps=err)
        
        w_tree8 = ACL(self.tree8)
        self.assertFloatEqual(w_tree8['A'], w_tree8['B'],err)
        self.assertFloatEqual(w_tree8['A'], w_tree8['C'],err)
        self.assertFloatEqual(w_tree8['D'], w_tree8['E'],err)
        self.assertFloatEqual(w_tree8['F'], w_tree8['G'],err)

        self.assertGreaterThan(w_tree8['A'], w_tree8['D'])
        self.assertGreaterThan(w_tree8['D'], w_tree8['H'])
        self.assertGreaterThan(w_tree8['H'], w_tree8['F'])

class GscTests(GeneralTests):
    """Tests for GSC functionality"""

    def test_gsc(self):
        """GSC: should return correct weights"""
        err = 1e-3
        tree6_exp = {'A': 0.19025, 'B': 0.19025, 'C': 0.2717, 'D': 0.3478}
        tree7_exp = {'A':.25, 'B':.25, 'C':.25, 'D':.25}
        tree8_exp = dict(zip('ABCDEFGH',[.1,.1,.2,.06,.06,.16,.16,.16]))
        self.assertFloatEqualAbs(GSC(self.tree6).values(),\
            tree6_exp.values(),eps=err)
        self.assertFloatEqualAbs(GSC(self.tree7).values(),\
            tree7_exp.values(),eps=err)
        self.assertFloatEqualAbs(GSC(self.tree8).values(),\
            tree8_exp.values(),eps=err)


#Rooted tree relating 15 HIV-1 isolates. From Altschul (1989) Fig 2.
TREE_1 = "(((((((((((BH8:0.7,PV22:0.3,BH10:0.3):0.1,BRU:0.5):0.1,HXB:0.7):2.4,SF2:3.3):0.1,CDC:3.7):0.5),((WMJ1:0.8,WMJ2:0.9,WMJ3:0.9):2.1)):0.4,RF:4.3):2.6),(((Z6:2.2,ELI:4.2):2.1,MAL:6.1):1.9)):2.7,Z3:9.3);"

#Model tree from Vingron and Sibbald (1993) Fig 3, distances estimated by Li.
TREE_2 = "(((((agcta:0.0,aggta:1.03):0.0),acctg:2.23):0.6),tgcaa:1.69);"

#Model tree from Vingron and Sibbald (1993) Fig 3 Actual # of substitutions
TREE_3 = "(((((agcta:0,aggta:1):1),acctg:1):1),tgcaa:2);"

#the ultrameric tree from Vingron and Sibbald (1993) Fig 3.
TREE_4 ="(((((((((2:6.0,1:6.0):12.9),((8:17.0,7:17.0):1.9)):5.3)),((6:8.5,5:8.5):15.7)):9.6),((((4:15.0,3:15.0):12.1),((9:11.0,10:11.0):16.1)):6.7)));"

#the additive tree from Vingron and Sibbald (1993) Fig 3.
#I don't trust the results they got for this tree (Table 3).
TREE_5 ="(((((((((2:7,1:7):19),((8:18,7:16):3)):12)),((6:10,5:10):28)):16),((((4:14,3:18):15),((9:8,10:14):24)):11)));"

TREE_6 = "(((A:20,B:20):30,C:50):30,D:80);"
TREE_7 = "((A:10,B:10):5,(C:10,D:10):5);"
TREE_8 = "((((A:5,B:5):5,C:15):10),((D:0,E:0):10,((F:5,G:5):10,H:10):10):10);"
TREE_9 = "(((A:0,B:0):5,C:10):5,D:25);"

FIVE_S_ALN =\
"""CLUSTAL W  (1.81)

1               A----TCCACGGCCATAGGACTCTGAAAGCACTGCATCCCGT-CCGATCTGCAAAGTTAA
2               A----TCCACGGCCATAGGACTGTGAAAGCACCGCATCCCGT-CTGATCTGCGCAGTTAA
3               T----CTGGTGATGATGGCGGAGGGGACACACCCGTTCCCATACCGAACACGGCCGTTAA
4               T----CTGGTGGCGATAGCGAGAAGGTCACACCCGTTCCCATACCGAACACGGAAGTTAA
5               G---TGGTGCGGTCATACCAGCGCTAATGCACCGGATCCCAT-CAGAACTCCGCAGTTAA
6               G----GGTGCGATCATACCAGCGTTAATGCACCGGATCCCAT-CAGAACTCCGCAGTTAA
7               G----CTTACGACCATATCACGTTGAATGCACGCCATCCCGT-CCGATCTGGCAAGTTAA
8               G----CCTACGGCCATCCCACCCTGGTAACGCCCGATCTCGT-CTGATCTCGGAAGCTAA
9               T--T-CTGGTGTCTCAGGCGTGGAGGAACCACACCAATCCATCCCGAACTTGGTGGTGAA
10              TATT-CTGGTGTCCCAGGCGTAGAGGAACCACACCGATCCATCTCGAACTTGGTGGTGAA
                        . *   .:   .     .:  *.*    :  *.*   **:*:     *  **

1               CCAGAGTACCGCCCAGT-TAGTACC-AC-GGTGGGGGACCACGCGGGAATCCTGGGTGCT
2               ACACAGTGCCGCCTAGT-TAGTACC-AT-GGTGGGGGACCACATGGGAATCCTGGGTGCT
3               GCCCTCCAGCGCC--AA-TGGTACT-TGCTC-CGCAGGGAG-CCGGGAGAGTAGGACGTC
4               GCTTCTCAGCGCC--GA-TGGTAGT-TA-GG-GGCTGTCCC-CTGTGAGAGTAGGACGCT
5               GCGCGCTTGGGCCAGAA-CAGTACT-GG-GATGGGTGACCTCCCGGGAAGTCCTGGTGCC
6               GCGCGCTTGGGTTGGAG-TAGTACT-AG-GATGGGTGACCTCCTGGGAAGTCCTAATATT
7               GCAACGTTGAGTCCAGT-TAGTACT-TG-GATCGGAGACGGCCTGGGAATCCTGGATGTT
8               GCAGGGTCGGGCCTGGT-TAGTACT-TG-GATGGGAGACCTCCTGGGAATACCGGGTGCT
9               ACTCTATTGCGGT--GA-CGATACTGTA-GG-GGAAGCCCG-ATGGAAAAATAGCTCGAC
10              ACTCTGCCGCGGT--AACCAATACT-CG-GG-GGGGGCCCT-GCGGAAAAATAGCTCGAT
                 *        *    .  ..**        *  *       * .*.        .  *  

1               GT-GG-T--T-
2               GT-GG-T--T-
3               GCCAG-G--C-
4               GCCAG-G--C-
5               GCACC-C--C-
6               GCACC-C-TT-
7               GTAAG-C--T-
8               GTAGG-CT-T-
9               GCCAGGA--T-
10              GCCAGGA--TA
                   :        
"""

if __name__ == "__main__":
    main()

