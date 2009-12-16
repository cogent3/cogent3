#!/usr/bin/env python
"""Tests for the MEME parser
"""
from __future__ import division
from cogent.util.unit_test import TestCase, main
import string
import re
from cogent.motif.util import Motif
from cogent.core.moltype import DNA
from cogent.parse.record import DelimitedSplitter
from cogent.parse.record_finder import LabeledRecordFinder
from cogent.parse.meme import *

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jermey.widmann@colorado.edu"
__status__ = "Production"

class MemeTests(TestCase):
    """Tests for meme module.
    """

    def setUp(self):
        """Setup function for meme tests.
        """
        #Meme output data:
        self.meme_file = MEME_FILE.split('\n')            
        self.meme_main = LabeledRecordFinder(lambda x: x.startswith('COMMAND'))
        self.meme_command = LabeledRecordFinder(lambda x: x.startswith('MOTIF'))
        self.meme_summary = LabeledRecordFinder(lambda x: x.startswith('SUMMARY'))
        self.meme_module = LabeledRecordFinder(lambda x: x.startswith('Motif'))
        self.alphabet_block, self.main_block = \
                             list(self.meme_main(self.meme_file))
        self.cmd_mod_list = list(self.meme_command(self.main_block))
        self.command_block = self.cmd_mod_list[0]
        self.module_blocks = self.cmd_mod_list[1:]                            
        self.summary_block = list(self.meme_summary(self.module_blocks[-1]))[1]
        self.module_data_blocks = []
        for module in self.module_blocks:
            self.module_data_blocks.append(\
                list(self.meme_module(module)))

        #List and Dict for testing dictFromList function
        self.sample_list = ['key1',1,'key2',2,'key3',3,'key4',4]
        self.sample_dict = {'key1':1,
                            'key2':2,
                            'key3':3,
                            'key4':4,
                            }

        #List of command line data
        self.command_line_list = [
            'model:  mod=           tcm    nmotifs=         3    evt=        1e+100',
            'object function=  E-value of product of p-values',
            'width:  minw=            4    maxw=           10    minic=        0.00',
            'width:  wg=             11    ws=              1    endgaps=       yes',
            'nsites: minsites=        2    maxsites=       50    wnsites=       0.8',
            'theta:  prob=            1    spmap=         uni    spfuzz=        0.5',
            'em:     prior=   dirichlet    b=            0.01    maxiter=        20',
            'distance=    1e-05',
            'data:   n=             597    N=              15',
            'strands: +',
            'sample: seed=            0    seqfrac=         1',
            ]

        #List of dicts which contain general info for each module.
        self.module_info_dicts = [
            {'MOTIF':'1',
             'width':'10',
             'sites':'11',
             'llr':'131',
             'E-value':'1.3e-019',
             },
            {'MOTIF':'2',
             'width':'7',
             'sites':'11',
             'llr':'88',
             'E-value':'2.5e-006',
             },
            {'MOTIF':'3',
             'width':'7',
             'sites':'6',
             'llr':'53',
             'E-value':'5.5e-001',
             },
            ]

        #Summary dict
        self.summary_dict = {'CombinedP':{
            '1': float(3.48e-02),
            '11': float(3.78e-05),
            '17': float(2.78e-08),
            '28': float(3.49e-06),
            '105': float(3.98e-06),
            '159': float(1.08e-02),
            '402-C01': float(4.22e-07),
            '407-A07': float(7.32e-08),
            '410-A10': float(4.23e-04),
            '505-D01': float(5.72e-07),
            '507-B04-1': float(1.01e-04),
            '518-D12': float(2.83e-06),
            '621-H01': float(8.69e-07),
            '625-H05': float(8.86e-06),
            '629-C08': float(5.61e-07),
            }
                             }
        self.remap_dict = {
            '11':'11',
            '1':'1',
            '407-A07':'407-A07',
            '17':'17',
            '159':'159',
            '505-D01':'505-D01',
            '28':'28',
            '507-B04-1':'507-B04-1',
            '402-C01':'402-C01',
            '621-H01':'621-H01',
            '629-C08':'629-C08',
            '410-A10':'410-A10',
            '105':'105',
            '625-H05':'625-H05',
            '518-D12':'518-D12'
            }


        #ModuleInstances and Modules
        self.ModuleInstances = [
            [ModuleInstance('CTATTGGGGC',Location('629-C08',18,28),
                            float(1.95e-06)),
             ModuleInstance('CTATTGGGGC',Location('621-H01',45,55),
                            float(1.95e-06)),
             ModuleInstance('CTATTGGGGC',Location('505-D01',26,36),
                            float(1.95e-06)),
             ModuleInstance('CTATTGGGGC',Location('407-A07',5,15),
                            float(1.95e-06)),
             ModuleInstance('CTATTGGGGC',Location('105',0,10),
                            float(1.95e-06)),
             ModuleInstance('CTATTGGGGC',Location('28',3,13),
                            float(1.95e-06)),
             ModuleInstance('CTATTGGGGC',Location('17',16,26),
                            float(1.95e-06)),
             ModuleInstance('CTATTGGGCC',Location('402-C01',24,34),
                            float(3.30e-06)),
             ModuleInstance('CTAGTGGGGC',Location('625-H05',2,12),
                            float(5.11e-06)),
             ModuleInstance('CTAGTGGGCC',Location('11',15,25),
                            float(6.37e-06)),
             ModuleInstance('CTATTGGGGT',Location('518-D12',0,10),
                            float(9.40e-06)),
             ],
            [ModuleInstance('CGTTACG',Location('629-C08',37,44),
                            float(6.82e-05)),
             ModuleInstance('CGTTACG',Location('621-H01',30,37),
                            float(6.82e-05)),
             ModuleInstance('CGTTACG',Location('507-B04-1',8,15),
                            float(6.82e-05)),
             ModuleInstance('CGTTACG',Location('410-A10',7,14),
                            float(6.82e-05)),
             ModuleInstance('CGTTACG',Location('407-A07',26,33),
                            float(6.82e-05)),
             ModuleInstance('CGTTACG',Location('17',0,7),
                            float(6.82e-05)),
             ModuleInstance('TGTTACG',Location('625-H05',32,39),
                            float(1.74e-04)),
             ModuleInstance('TGTTACG',Location('505-D01',3,10),
                            float(1.74e-04)),
             ModuleInstance('CATTACG',Location('518-D12',30,37),
                            float(2.14e-04)),
             ModuleInstance('CGGTACG',Location('402-C01',1,8),
                            float(2.77e-04)),
             ModuleInstance('TGTTCCG',Location('629-C08',5,12),
                            float(6.45e-04)),
             ],
            [ModuleInstance('CTATTGG',Location('629-C08',57,64),
                            float(1.06e-04)),
             ModuleInstance('CTATTGG',Location('507-B04-1',42,49),
                            float(1.06e-04)),
             ModuleInstance('CTATTGG',Location('410-A10',27,34),
                            float(1.06e-04)),
             ModuleInstance('CTATTGG',Location('159',14,21),
                            float(1.06e-04)),
             ModuleInstance('CTATTGG',Location('1',18,25),
                            float(1.06e-04)),
             ModuleInstance('CTAATGG',Location('507-B04-1',28,35),
                            float(1.63e-04)),
            ],
        ]

             
        self.Modules = []
        for module, info in zip(self.ModuleInstances, self.module_info_dicts):
            curr_module_data = {}
            for instance in module:
                curr_module_data[(instance.Location.SeqId,
                             instance.Location.Start)] = instance
            temp_module = Module(curr_module_data, MolType=DNA,
                                 Evalue=float(info['E-value']),
                                 Llr=int(info['llr']))
            self.Modules.append(temp_module)
        
        self.ConsensusSequences = ['CTATTGGGGC','CGTTACG','CTATTGG']
            

    def test_get_data_block(self):
        """getDataBlock should return the main block and the alphabet."""
        main_block, alphabet = getDataBlock(self.meme_file)
        self.assertEqual(main_block,self.main_block)
        self.assertEqual(alphabet, DNA)

    def test_get_alphabet(self):
        """getMolType should return the correct alphabet."""
        self.assertEqual(getMolType(self.alphabet_block),DNA)

    def test_get_command_module_blocks(self):
        """getCommandModuleBlocks should return the command and module blocks.
        """
        command_block, module_blocks = getCommandModuleBlocks(self.main_block)
        self.assertEqual(command_block, self.command_block)
        self.assertEqual(module_blocks, self.module_blocks)

    def test_get_summary_block(self):
        """getSummaryBlock should return the MEME summary block."""
        self.assertEqual(getSummaryBlock(self.module_blocks[-1]),
                         self.summary_block)
        
    def test_dict_from_list(self):
        """dictFromList should return a dict given a list."""
        self.assertEqual(dictFromList(self.sample_list),self.sample_dict)

    def test_extract_command_line_data(self):
        """extractCommandLineData should return a dict of command line data."""
        self.assertEqual(extractCommandLineData(self.command_block),
                         self.command_line_list)

    def test_get_module_data_blocks(self):
        """getModuleDataBlocks should return a list of blocks for each module.
        """
        self.assertEqual(getModuleDataBlocks(self.module_blocks),
                         self.module_data_blocks)

    def test_extract_module_data(self):
        """extractModuleData should return a Module object."""
        for data, module in zip(self.module_data_blocks,self.Modules):
            ans = extractModuleData(data,DNA,self.remap_dict)
            self.assertEqual(ans,module)
    
    def test_get_consensus_sequence(self):  
        """getConsensusSequence should return Module's Consensus sequence."""
        for data,seq in zip(self.module_data_blocks,self.ConsensusSequences):
            ans =  getConsensusSequence(data[1])
            self.assertEqual(ans,seq)
        
    def test_get_module_general_info(self):
        """getModuleGeneralInfo should return a dict of Module info."""
        for module, data_dict in zip(self.module_data_blocks,
                                     self.module_info_dicts):
            self.assertEqual(getModuleGeneralInfo(module[0][0]),data_dict)
        #Test that getModuleGeneralInfo can parse general info line when
        # motif ID is > 100.  MEME changes the format of this line when in this
        # case.
        data_line_special = \
        'MOTIF100	width =   50   sites =   2   llr = 273   E-value = 3.1e-007'
        expected = {'MOTIF':'100','width':'50','sites':'2','llr':'273',\
            'E-value':'3.1e-007'}
        self.assertEqual(getModuleGeneralInfo(data_line_special),expected)

    def test_extract_summary_data(self):
        """extractSummaryData should return a dict of MEME summary data."""
        self.assertEqual(extractSummaryData(self.summary_block),
                         self.summary_dict)

    def test_meme_parser(self):
        """MemeParser should correctly return a MotifResults object."""
        test_motif_results = MotifResults([],[],{})
        test_motif_results.Results = self.summary_dict
        test_motif_results.Results['Warnings']=[]
        test_motif_results.Parameters = self.command_line_list
        test_motif_results.Modules = self.Modules
        ans_motif_results = MemeParser(self.meme_file)
        self.assertEqual(ans_motif_results.Modules,test_motif_results.Modules)
        self.assertEqual(ans_motif_results.Results,test_motif_results.Results)
        self.assertEqual(ans_motif_results.Parameters,
                         test_motif_results.Parameters)

MEME_FILE = """
********************************************************************************
MEME - Motif discovery tool
********************************************************************************
MEME version 3.0 (Release date: 2001/03/03 13:05:22)

For further information on how to interpret these results or to get
a copy of the MEME software please access http://meme.sdsc.edu.

This file may be used as input to the MAST algorithm for searching
sequence databases for matches to groups of motifs.  MAST is available
for interactive use and downloading at http://meme.sdsc.edu.
********************************************************************************


********************************************************************************
REFERENCE
********************************************************************************
If you use this program in your research, please cite:

Timothy L. Bailey and Charles Elkan,
"Fitting a mixture model by expectation maximization to discover
motifs in biopolymers", Proceedings of the Second International
Conference on Intelligent Systems for Molecular Biology, pp. 28-36,
AAAI Press, Menlo Park, California, 1994.
********************************************************************************


********************************************************************************
TRAINING SET
********************************************************************************
DATAFILE= meme.16346.data
ALPHABET= ACGT
Sequence name            Weight Length  Sequence name            Weight Length  
-------------            ------ ------  -------------            ------ ------  
1                        1.0000     26  11                       1.0000     25  
17                       1.0000     26  28                       1.0000     26  
105                      1.0000     21  159                      1.0000     21  
402-C01                  1.0000     34  407-A07                  1.0000     34  
410-A10                  1.0000     34  505-D01                  1.0000     49  
507-B04-1                1.0000     49  518-D12                  1.0000     49  
621-H01                  1.0000     74  625-H05                  1.0000     65  
629-C08                  1.0000     64  
********************************************************************************

********************************************************************************
COMMAND LINE SUMMARY
********************************************************************************
This information can also be useful in the event you wish to report a
problem with the MEME software.

command: meme meme.16346.data -dna -mod tcm -nmotifs 3 -minw 4 -maxw 10 -evt 1e100 -time 720 -maxsize 60000 -nostatus -maxiter 20 

model:  mod=           tcm    nmotifs=         3    evt=        1e+100
object function=  E-value of product of p-values
width:  minw=            4    maxw=           10    minic=        0.00
width:  wg=             11    ws=              1    endgaps=       yes
nsites: minsites=        2    maxsites=       50    wnsites=       0.8
theta:  prob=            1    spmap=         uni    spfuzz=        0.5
em:     prior=   dirichlet    b=            0.01    maxiter=        20
        distance=    1e-05
data:   n=             597    N=              15
strands: +
sample: seed=            0    seqfrac=         1
Letter frequencies in dataset:
A 0.173 C 0.206 G 0.299 T 0.322 
Background letter frequencies (from dataset with add-one prior applied):
A 0.173 C 0.207 G 0.298 T 0.322 
********************************************************************************


********************************************************************************
MOTIF  1 width =   10   sites =  11   llr = 131   E-value = 1.3e-019
********************************************************************************
--------------------------------------------------------------------------------
Motif 1 Description
--------------------------------------------------------------------------------
Simplified        A  ::a:::::::
pos.-specific     C  a:::::::29
probability       G  :::2:aaa8:
matrix            T  :a:8a::::1

         bits    2.5   *       
                 2.3 * *       
                 2.0 * *       
                 1.8 * *  *** *
Information      1.5 *** **** *
content          1.3 *** ******
(17.2 bits)      1.0 **********
                 0.8 **********
                 0.5 **********
                 0.3 **********
                 0.0 ----------

Multilevel           CTATTGGGGC
consensus                      
sequence                       
                               
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 1 sites sorted by position p-value
--------------------------------------------------------------------------------
Sequence name             Start   P-value               Site 
-------------             ----- ---------            ----------
629-C08                      19  1.95e-06 TCCGTGAACA CTATTGGGGC GTGTAAGAGC
621-H01                      46  1.95e-06 CGCATGCGTG CTATTGGGGC GTCATTTGTC
505-D01                      27  1.95e-06 TTGATTGTTG CTATTGGGGC ATTGCCGTAC
407-A07                       6  1.95e-06      CGTTA CTATTGGGGC GGGTATTTTC
105                           1  1.95e-06          . CTATTGGGGC CGAAATGGTT
28                            4  1.95e-06        TCC CTATTGGGGC CAAGGGCTAC
17                           17  1.95e-06 GCTACTTGTG CTATTGGGGC           
402-C01                      25  3.30e-06 CTTAACATTC CTATTGGGCC           
625-H05                       3  5.11e-06         GC CTAGTGGGGC AGCTGACAGA
11                           16  6.37e-06 TGTTAGACAG CTAGTGGGCC           
518-D12                       1  9.40e-06          . CTATTGGGGT GTTGTATTGA
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 1 block diagrams
--------------------------------------------------------------------------------
SEQUENCE NAME            POSITION P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
629-C08                             2e-06  18_[1]_36
621-H01                             2e-06  45_[1]_19
505-D01                             2e-06  26_[1]_13
407-A07                             2e-06  5_[1]_19
105                                 2e-06  [1]_11
28                                  2e-06  3_[1]_13
17                                  2e-06  16_[1]
402-C01                           3.3e-06  24_[1]
625-H05                           5.1e-06  2_[1]_53
11                                6.4e-06  15_[1]
518-D12                           9.4e-06  [1]_39
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 1 in BLOCKS format
--------------------------------------------------------------------------------
BL   MOTIF 1 width=10 seqs=11
629-C08                  (   19) CTATTGGGGC  1 
621-H01                  (   46) CTATTGGGGC  1 
505-D01                  (   27) CTATTGGGGC  1 
407-A07                  (    6) CTATTGGGGC  1 
105                      (    1) CTATTGGGGC  1 
28                       (    4) CTATTGGGGC  1 
17                       (   17) CTATTGGGGC  1 
402-C01                  (   25) CTATTGGGCC  1 
625-H05                  (    3) CTAGTGGGGC  1 
11                       (   16) CTAGTGGGCC  1 
518-D12                  (    1) CTATTGGGGT  1 
//

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 1 position-specific scoring matrix
--------------------------------------------------------------------------------
log-odds matrix: alength= 4 w= 10 n= 462 bayes= 6.40183 E= 1.3e-019 
 -1010    227  -1010  -1010 
 -1010  -1010  -1010    164 
   253  -1010  -1010  -1010 
 -1010  -1010    -71    135 
 -1010  -1010  -1010    164 
 -1010  -1010    174  -1010 
 -1010  -1010    174  -1010 
 -1010  -1010    174  -1010 
 -1010    -19    146  -1010 
 -1010    214  -1010   -182 
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 1 position-specific probability matrix
--------------------------------------------------------------------------------
letter-probability matrix: alength= 4 w= 10 nsites= 11 E= 1.3e-019 
 0.000000  1.000000  0.000000  0.000000 
 0.000000  0.000000  0.000000  1.000000 
 1.000000  0.000000  0.000000  0.000000 
 0.000000  0.000000  0.181818  0.818182 
 0.000000  0.000000  0.000000  1.000000 
 0.000000  0.000000  1.000000  0.000000 
 0.000000  0.000000  1.000000  0.000000 
 0.000000  0.000000  1.000000  0.000000 
 0.000000  0.181818  0.818182  0.000000 
 0.000000  0.909091  0.000000  0.090909 
--------------------------------------------------------------------------------





Time  0.54 secs.

********************************************************************************


********************************************************************************
MOTIF  2 width =    7   sites =  11   llr = 88   E-value = 2.5e-006
********************************************************************************
--------------------------------------------------------------------------------
Motif 2 Description
--------------------------------------------------------------------------------
Simplified        A  :1::9::
pos.-specific     C  7:::1a:
probability       G  :91:::a
matrix            T  3:9a:::

         bits    2.5        
                 2.3      * 
                 2.0     ** 
                 1.8     ***
Information      1.5    ****
content          1.3 *******
(11.6 bits)      1.0 *******
                 0.8 *******
                 0.5 *******
                 0.3 *******
                 0.0 -------

Multilevel           CGTTACG
consensus            T      
sequence                    
                            
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 2 sites sorted by position p-value
--------------------------------------------------------------------------------
Sequence name             Start   P-value             Site
-------------             ----- ---------            -------
629-C08                      38  6.82e-05 CGTGTAAGAG CGTTACG TGTTCCGTGA
621-H01                      31  6.82e-05 CGAGGGAGTA CGTTACG CATGCGTGCT
507-B04-1                     9  6.82e-05   CTTGCACA CGTTACG TGTGAGCCAT
410-A10                       8  6.82e-05    CTTTGCT CGTTACG TGGTTGTATG
407-A07                      27  6.82e-05 GGTATTTTCC CGTTACG T         
17                            1  6.82e-05          . CGTTACG CTACTTGTGC
625-H05                      33  1.74e-04 ATAGGTCGAC TGTTACG GTTAGCGTTC
505-D01                       4  1.74e-04        GCA TGTTACG TGACTTTTGA
518-D12                      31  2.14e-04 GTTATTGCGA CATTACG CGTTCTGGTT
402-C01                       2  2.77e-04          C CGGTACG GTTTGTCTTA
629-C08                       6  6.45e-04      TTACG TGTTCCG TGAACACTAT
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 2 block diagrams
--------------------------------------------------------------------------------
SEQUENCE NAME            POSITION P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
629-C08                           0.00064  5_[2]_25_[2]_20
621-H01                           6.8e-05  30_[2]_37
507-B04-1                         6.8e-05  8_[2]_34
410-A10                           6.8e-05  7_[2]_20
407-A07                           6.8e-05  26_[2]_1
17                                6.8e-05  [2]_19
625-H05                           0.00017  32_[2]_26
505-D01                           0.00017  3_[2]_39
518-D12                           0.00021  30_[2]_12
402-C01                           0.00028  1_[2]_26
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 2 in BLOCKS format
--------------------------------------------------------------------------------
BL   MOTIF 2 width=7 seqs=11
629-C08                  (   38) CGTTACG  1 
621-H01                  (   31) CGTTACG  1 
507-B04-1                (    9) CGTTACG  1 
410-A10                  (    8) CGTTACG  1 
407-A07                  (   27) CGTTACG  1 
17                       (    1) CGTTACG  1 
625-H05                  (   33) TGTTACG  1 
505-D01                  (    4) TGTTACG  1 
518-D12                  (   31) CATTACG  1 
402-C01                  (    2) CGGTACG  1 
629-C08                  (    6) TGTTCCG  1 
//

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 2 position-specific scoring matrix
--------------------------------------------------------------------------------
log-odds matrix: alength= 4 w= 7 n= 507 bayes= 6.53743 E= 2.5e-006 
 -1010    181  -1010    -24 
   -93  -1010    161  -1010 
 -1010  -1010   -171    150 
 -1010  -1010  -1010    164 
   239   -118  -1010  -1010 
 -1010    227  -1010  -1010 
 -1010  -1010    174  -1010 
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 2 position-specific probability matrix
--------------------------------------------------------------------------------
letter-probability matrix: alength= 4 w= 7 nsites= 11 E= 2.5e-006 
 0.000000  0.727273  0.000000  0.272727 
 0.090909  0.000000  0.909091  0.000000 
 0.000000  0.000000  0.090909  0.909091 
 0.000000  0.000000  0.000000  1.000000 
 0.909091  0.090909  0.000000  0.000000 
 0.000000  1.000000  0.000000  0.000000 
 0.000000  0.000000  1.000000  0.000000 
--------------------------------------------------------------------------------





Time  0.85 secs.

********************************************************************************


********************************************************************************
MOTIF  3 width =    7   sites =   6   llr = 53   E-value = 5.5e-001
********************************************************************************
--------------------------------------------------------------------------------
Motif 3 Description
--------------------------------------------------------------------------------
Simplified        A  ::a2:::
pos.-specific     C  a::::::
probability       G  :::::aa
matrix            T  :a:8a::

         bits    2.5   *    
                 2.3 * *    
                 2.0 * *    
                 1.8 * *  **
Information      1.5 *** ***
content          1.3 *** ***
(12.7 bits)      1.0 *******
                 0.8 *******
                 0.5 *******
                 0.3 *******
                 0.0 -------

Multilevel           CTATTGG
consensus                   
sequence                    
                            
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 3 sites sorted by position p-value
--------------------------------------------------------------------------------
Sequence name             Start   P-value             Site
-------------             ----- ---------            -------
629-C08                      58  1.06e-04 TCCGTGAACA CTATTGG           
507-B04-1                    43  1.06e-04 TGGTGTTGCG CTATTGG           
410-A10                      28  1.06e-04 TTGTATGCCG CTATTGG           
159                          15  1.06e-04 GACCGTTGGT CTATTGG           
1                            19  1.06e-04 TTGGATAGTG CTATTGG G         
507-B04-1                    29  1.63e-04 GAGCCATTCT CTAATGG TGTTGCGCTA
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 3 block diagrams
--------------------------------------------------------------------------------
SEQUENCE NAME            POSITION P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
629-C08                           0.00011  57_[3]
507-B04-1                         0.00016  28_[3]_7_[3]
410-A10                           0.00011  27_[3]
159                               0.00011  14_[3]
1                                 0.00011  18_[3]_1
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 3 in BLOCKS format
--------------------------------------------------------------------------------
BL   MOTIF 3 width=7 seqs=6
629-C08                  (   58) CTATTGG  1 
507-B04-1                (   43) CTATTGG  1 
410-A10                  (   28) CTATTGG  1 
159                      (   15) CTATTGG  1 
1                        (   19) CTATTGG  1 
507-B04-1                (   29) CTAATGG  1 
//

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 3 position-specific scoring matrix
--------------------------------------------------------------------------------
log-odds matrix: alength= 4 w= 7 n= 507 bayes= 6.83576 E= 5.5e-001 
  -923    227   -923   -923 
  -923   -923   -923    164 
   253   -923   -923   -923 
    -6   -923   -923    137 
  -923   -923   -923    164 
  -923   -923    174   -923 
  -923   -923    174   -923 
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Motif 3 position-specific probability matrix
--------------------------------------------------------------------------------
letter-probability matrix: alength= 4 w= 7 nsites= 6 E= 5.5e-001 
 0.000000  1.000000  0.000000  0.000000 
 0.000000  0.000000  0.000000  1.000000 
 1.000000  0.000000  0.000000  0.000000 
 0.166667  0.000000  0.000000  0.833333 
 0.000000  0.000000  0.000000  1.000000 
 0.000000  0.000000  1.000000  0.000000 
 0.000000  0.000000  1.000000  0.000000 
--------------------------------------------------------------------------------





Time  1.09 secs.

********************************************************************************


********************************************************************************
SUMMARY OF MOTIFS
********************************************************************************

--------------------------------------------------------------------------------
Combined block diagrams: non-overlapping sites with p-value < 0.0001
--------------------------------------------------------------------------------
SEQUENCE NAME            COMBINED P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
1                                3.48e-02  26
11                               3.78e-05  15_[1(6.37e-06)]
17                               2.78e-08  [2(6.82e-05)]_9_[1(1.95e-06)]
28                               3.49e-06  3_[1(1.95e-06)]_13
105                              3.98e-06  [1(1.95e-06)]_11
159                              1.08e-02  21
402-C01                          4.22e-07  24_[1(3.30e-06)]
407-A07                          7.32e-08  5_[1(1.95e-06)]_11_[2(6.82e-05)]_1
410-A10                          4.23e-04  7_[2(6.82e-05)]_20
505-D01                          5.72e-07  26_[1(1.95e-06)]_13
507-B04-1                        1.01e-04  8_[2(6.82e-05)]_34
518-D12                          2.83e-06  [1(9.40e-06)]_39
621-H01                          8.69e-07  30_[2(6.82e-05)]_8_[1(1.95e-06)]_19
625-H05                          8.86e-06  2_[1(5.11e-06)]_53
629-C08                          5.61e-07  18_[1(1.95e-06)]_9_[2(6.82e-05)]_20
--------------------------------------------------------------------------------

********************************************************************************


********************************************************************************
Stopped because nmotifs = 3 reached.
********************************************************************************

CPU: compute-0-2.local

********************************************************************************
"""
        
#run if called from command-line
if __name__ == "__main__":
    main()
