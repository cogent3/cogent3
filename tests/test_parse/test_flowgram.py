#!/usr/bin/env python
"""tests for Flowgram and Flowgramcollection objects
"""

__author__ = "Jens Reeder, Julia Goodrich"
__copyright__ = "Copyright 2009, The Cogent Project"
__credits__ = ["Jens Reeder","Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jens Reeder"
__email__ = "jreeder@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from types import GeneratorType
from numpy import array, transpose

from cogent.core.sequence import Sequence
from cogent.parse.flowgram import Flowgram, build_averaged_flowgram
                               
from cogent.parse.flowgram_parser import parse_sff

class FlowgramTests(TestCase):
    def test_init_empty(self):
        """Flowgram should init correctly."""
        f = Flowgram()
        self.assertEqual(f._flowgram, '')
        self.assertEqual(f.flowgram, [])
    
    def test_init_data(self):
        """Flowgram init with data should set data in correct location"""
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',KeySeq = "ATCG",
                     floworder = "TACG", header_info = {'Bases':'TACCCC'})
        self.assertEqual(f._flowgram, '0.5 1.0 4.0 0.0')
        self.assertEqual(f.flowgram, [0.5, 1.0, 4.0, 0.0])
        self.assertEqual(f.Name, 'a')
        self.assertEqual(f.keySeq, "ATCG")
        self.assertEqual(f.floworder, "TACG")
        self.assertEqual(f.Bases, 'TACCCC')
        self.assertEqual(f.header_info, {'Bases':'TACCCC'})

        f = Flowgram([0.5, 1.0, 4.0, 0.0], Name = 'a',KeySeq = "ATCG",
                     floworder = "TACG", header_info = {'Bases':'TACCCC'})
        self.assertEqual(f._flowgram, '0.5 1.0 4.0 0.0')
        self.assertEqual(f.flowgram, [0.5, 1.0, 4.0, 0.0])

    def test_cmpSeqToString(self):
        """Sequence should compare equal to same string."""
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        self.assertTrue(f.cmpSeqToString('TACCCC'))
        self.assertFalse(f.cmpSeqToString('TACCC'))
        
        f = Flowgram('0.5 1.0 4.0 0.0',floworder = "TACG")
        self.assertTrue(f.cmpSeqToString('TACCCC'))
        self.assertFalse(f.cmpSeqToString('TACCC'))

    def test_cmp_flow_to_string(self):
        """Sequence should compare equal to same string."""
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        self.assertEqual(f, '0.5 1.0 4.0 0.0')
        self.assertNotEqual(f,'0.5 1.0 4.0')
        f2 = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        self.assertEqual(f,f2)
        
    def test_cmpBySeqs(self):
        """Flowgrams should be the same if name, bases, or to_seqs are equal"""
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        f2 = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        self.assertEqual(f.cmpBySeqs(f2), 0)
        f2 = Flowgram('0.5 1.0 4.0 0.0', Name = 'b',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        self.assertEqual(f.cmpBySeqs(f2), 0)
        
        f2 = Flowgram('0.5 1.0 4.0 0.0',floworder = "TACG")
        self.assertEqual(f.cmpBySeqs(f2), 0)

    def test_cmpByName(self):
        """Flowgrams should be the same if name, bases, or to_seqs are equal"""
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        f2 = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        self.assertEqual(f.cmpByName(f2), 0)
        self.assertEqual(f.cmpByName(f), 0)
        f2 = Flowgram('0.5 1.0 4.0 0.0', Name = 'b',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        self.assertNotEqual(f.cmpByName(f2), 0)
        
        
    def test_toFasta(self):
        """Flowgram toFasta() should return Fasta-format string"""
        even = '0.5 1.0 4.0 0.0'
        odd = '0.5 1.0 4.0 1.0'
        even_f = Flowgram(even, Name='even', floworder = "TACG")
        odd_f = Flowgram(odd, Name='odd', floworder = "TACG")
        self.assertEqual(even_f.toFasta(), '>even\nTACCCC')
        #set line wrap to small number so we can test that it works
        self.assertEqual(even_f.toFasta(LineWrap = 2), '>even\nTA\nCC\nCC')
        self.assertEqual(odd_f.toFasta(LineWrap = 2), '>odd\nTA\nCC\nCC\nG')

        even_f = Flowgram(even, Name='even', floworder = "TACG",
                          header_info ={'Bases':'TACCCG'})
        odd_f = Flowgram(odd, Name='odd', floworder = "TACG",
                          header_info ={'Bases':'TACCCGG'})
        self.assertEqual(even_f.toFasta(), '>even\nTACCCG')
        #set line wrap to small number so we can test that it works
        self.assertEqual(even_f.toFasta(LineWrap = 2), '>even\nTA\nCC\nCG')
        self.assertEqual(odd_f.toFasta(LineWrap = 2), '>odd\nTA\nCC\nCG\nG') 
    
    def test_contains(self):
        """Flowgram contains should return correct result"""
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCC'})
        assert '0.5' in f
        assert '0.5 1.0' in f
        assert '2.0' not in f
        assert '5.0' not in f

    def test_iter(self):
       """Flowgram iter should iterate over sequence"""
       f = Flowgram('0.5 1.0 4.0 0.0')
       self.assertEqual(list(f), [0.5,1.0,4.0,0.0])

       
    def test_str(self):
        """__str__ returns self._flowgram unmodified."""
        f = Flowgram('0.5 1.0 4.0 0.0')
        self.assertEqual(str(f), '0.5\t1.0\t4.0\t0.0')
        f = Flowgram([0.5, 1.0, 4.0, 0.0])
        self.assertEqual(str(f), '0.5\t1.0\t4.0\t0.0')

    def test_len(self):
        """returns the length of the flowgram"""
        f = Flowgram('0.5 1.0 4.0 0.0')
        self.assertEqual(len(f), 4)

        f = Flowgram()
        self.assertEqual(len(f), 0)

    def test_hash(self):
        """__hash__ behaves like the flowgram string for dict lookup."""
        f = Flowgram('0.5 1.0 4.0 0.0', floworder = "TACG")
        self.assertEqual(hash(f), hash('0.5 1.0 4.0 0.0'))

        f = Flowgram()
        self.assertEqual(hash(f), hash(''))

    def test_toSeq(self):
        """toSeq should Translate flowgram to sequence"""
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCG'})
        self.assertEqual(f.toSeq(), 'TACCCG')
        self.assertEqual(isinstance(f.toSeq(),Sequence), True)
        self.assertEqual(f.toSeq(Bases = False), 'TACCCC')

        f = Flowgram('0.5 1.0 4.0 0.0 0.0 1.23 0.0 6.1',
                     Name = 'a',floworder = "TACG",
                     header_info = {'Bases':'TACCCG'})
        self.assertEqual(f.toSeq(), 'TACCCG')
        self.assertEqual(f.toSeq(Bases = False), 'TACCCCAGGGGGG')
   
        f = Flowgram('0.5 1.0 4.0 0.0', Name = 'a',floworder = "TACG",
                     header_info = {})
        self.assertEqual(f.toSeq(), 'TACCCC')
        self.assertEqual(isinstance(f.toSeq(),Sequence), True)
        self.assertEqual(f.toSeq(Bases = False), 'TACCCC')

        f = Flowgram('0.5 1.0 4.0 0.0 0.0 1.23 0.0 6.1',
                     Name = 'a',floworder = "TACG",
                     header_info = {})
        self.assertEqual(f.toSeq(Bases = True), 'TACCCCAGGGGGG')

        f = Flowgram('0.4 0.0 0.0 0.0 0.0 1.23 0.0 1.1',
                     Name = 'a',floworder = "TACG",
                     header_info = {})
        self.assertEqual(f.toSeq(), 'NAG')

    def test_getQualityTrimmedFlowgram(self):
        """getQualityTrimmedFlowgram trims the flowgram correctly"""
        f = Flowgram('0.5 1.0 4.1 0.0 0.0 1.23 0.0 3.1',
                     Name = 'a', floworder = "TACG",
                     header_info = {'Bases':'TACCCCAGGG', 'Clip Qual Right': 7,
                                    'Flow Indexes': "1\t2\t3\t3\t3\t3\t6\t8\t8\t8"})
        trimmed = f.getQualityTrimmedFlowgram()
        
        self.assertEqual(trimmed.toSeq(), "TACCCCA")
        self.assertEqual(str(trimmed), "0.5\t1.0\t4.1\t0.0\t0.0\t1.23")
  
        # tests on real data
        flow1 = self.flows[0]
        flow2 = self.flows[1]

        flow1_trimmed = flow1.getQualityTrimmedFlowgram()
        self.assertEqual(str(flow1_trimmed), "1.06	0.08	1.04	0.08	0.05	0.94	0.10	2.01	0.10	0.07	0.96	0.09	1.04	1.96	1.07	0.10	1.01	0.13	0.08	1.01	1.06	1.83	2.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97")
        self.assertEqual(flow1_trimmed.Bases,
                         "tcagGCTAACTGTAACCCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCA")

        flow2_trimmed = flow2.getQualityTrimmedFlowgram()
        self.assertEqual(str(flow2_trimmed), "1.04	0.00	1.01	0.00	0.00	1.00	0.00	1.00	0.00	1.05	0.00	0.91	0.10	1.07	0.95	1.01	0.00	0.06	0.93	0.02	0.03	1.06	1.18	0.09	1.00	0.05	0.90	0.11	0.07	1.99	0.11	0.02	1.96	1.04	0.13	0.01	2.83	0.10	1.97	0.06	0.11	1.04	0.13	0.03	0.98	1.15	0.07	1.00	0.07	0.08	0.98	0.11	1.92	0.05	0.04	2.96	1.02	1.02	0.04	0.93	1.00	0.13	0.04	1.00	1.03	0.08	0.97	0.13	0.11	1.88	0.09	0.05	1.02	1.89	0.07	0.11	0.98	0.05	0.07	1.01	0.08	0.05	1.01	0.13	1.00	0.07	0.10	1.04	0.10	0.04	0.98	0.12	1.03	0.96	0.11	0.07	1.00	0.09	0.03	1.03	0.11	1.95	1.06	0.13	0.05	1.00	0.13	0.11	1.00	0.09	0.03	2.89	0.08	0.95	0.09	1.03	1.02	1.05	1.07	0.08	0.12	2.81	0.08	0.08	1.00	1.07	0.07	0.05	1.86	0.12	0.98	0.06	2.00	0.11	1.02	0.11	0.08	1.88	0.13	1.03	0.13	0.98	0.15	0.11	1.03	1.03	1.04	0.18	0.98	0.13	0.15	1.04	0.11	1.01	0.13	0.06	1.01	0.06	1.02	0.08	0.99	0.14	0.99	0.09	0.05	1.09	0.04	0.07	2.96	0.09	2.03	0.13	2.96	1.13	0.08	1.03	0.07	0.99	0.11	0.05	1.05	1.04	0.09	0.07	1.00	1.03	0.09	0.06	1.06	1.04	2.94	0.18	0.06	0.93	0.10	1.10	0.11	2.02	0.17	1.00	1.03	0.06	0.11	0.96	0.04	3.00	0.11	0.07	1.99	0.10	2.03	0.12	0.97	0.16	0.01	2.09	0.14	1.04	0.16	0.06	1.03	0.14	1.12	0.12	0.05	0.96	1.01	0.10	0.14	0.94	0.03	0.12	1.10	0.92	0.09	1.10	1.04	1.02	0.12	0.97	2.00	0.15	1.08	0.04	1.03	1.04	0.03	0.09	5.16	1.02	0.09	0.13	2.66	0.09	0.05	1.06	0.07	0.89	0.05	0.12	1.10	0.16	0.06	1.01	0.13	1.00	0.14	0.98	0.09	2.92	1.28	0.03	2.95	0.98	0.16	0.08	0.95	0.96	1.09	0.08	1.07	1.01	0.16	0.06	4.52	0.12	1.03	0.07	0.09	1.03	0.14	0.03	1.01	1.99")
        self.assertEqual(flow2_trimmed.Bases,
                          "tcagAGACGCACTCAATTATTTCCATAGCTTGGGTAGTGTCAATAATGCTGCTATGAACATGGGAGTACAAATATTCTTCAAGATACTGATCTCATTTCCTTTAGATATATACCCAGAAGTGAAATTCCTGGATCACATAGTAGTTCTATTTTTATTTGATGAGAAACTTTATACTATTTTTCATAA")

    def test_getPrimerTrimmedFlowgram(self):
        """getPrimerTrimmedFlowgram cuts the barcode of the flowgram correctly"""

        f = Flowgram('0.5 1.0 4.1 0.0 0.0 1.23 0.0 3.1',
                     Name = 'a', floworder = "TACG",
                     header_info = {'Bases':'TACCCCAGGG', 'Clip Qual Right': 7,
                                    'Flow Indexes': "1\t2\t3\t3\t3\t3\t6\t8\t8\t8"})

        trimmed = f.getPrimerTrimmedFlowgram(primerseq="TA")
        #test primer trimming
        self.assertEqual(trimmed.toSeq(), "CCCCAGGG")
        self.assertEqual(str(trimmed), "0.00\t0.00\t4.10\t0.00\t0.00\t1.23\t0.00\t3.10")
        for (a,b) in zip(trimmed.flowgram, [0.0,0.0,4.1,0.0,0.0,1.23,0.0,3.1]):
            self.assertFloatEqual(a,b)

        trimmed = f.getPrimerTrimmedFlowgram(primerseq="TACC")
        for (a,b) in zip(trimmed.flowgram, [0.0,0.0,2.1,0.0,0.0,1.23,0.0,3.1]):
            self.assertFloatEqual(a,b)
        self.assertEqual(trimmed.toSeq(), "CCAGGG")
        self.assertEqual(str(trimmed), "0.00\t0.00\t2.10\t0.00\t0.00\t1.23\t0.00\t3.10")

        # test that primer trimming does not leave ambig flow at begin
        trimmed = f.getPrimerTrimmedFlowgram(primerseq="TACCCC")
        for (a,b) in zip(trimmed.flowgram, [0.0,1.23,0.0,3.1]):
            self.assertFloatEqual(a,b)
        self.assertEqual(trimmed.toSeq(), "AGGG")
        self.assertEqual(str(trimmed), "0.00\t1.23\t0.00\t3.10")
        
        # tests on real data
        flow1 = self.flows[0]
        flow2 = self.flows[1]
        flow3 = self.flows[2]

        flow1_trimmed = flow1.getPrimerTrimmedFlowgram(primerseq="TCAG"+"GCTAACTGTAA")
        self.assertEqual(str(flow1_trimmed), "0.00\t0.00\t2.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04")
        self.assertEqual(flow1_trimmed.Bases,
                         "CCCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc")

        flow1_trimmed = flow1.getPrimerTrimmedFlowgram(primerseq="TCAG"+"GCTAACTGTAAC")
        self.assertEqual(str(flow1_trimmed), "0.00\t0.00\t1.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04")
        self.assertEqual(flow1_trimmed.Bases,
                         "CCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc")
 
        #test that trimming does not leave 4 zero flows (homopolymer)
        flow1_trimmed = flow1.getPrimerTrimmedFlowgram(primerseq="TCAG"+"GCTAACTGTAACCC")
        self.assertEqual(str(flow1_trimmed), "0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04")
        self.assertEqual(flow1_trimmed.Bases,
                         "TCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc")
     
        #test that trimming does not leave 4 zero flows (signal <1.5)
        flow1_trimmed = flow1.getPrimerTrimmedFlowgram(primerseq="TCAG"+"GCTAACTGTAACCCTC")
        self.assertEqual(str(flow1_trimmed), "1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04")
        self.assertEqual(flow1_trimmed.Bases,
                         "TTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc")
 
        flow1_untrimmed= flow1.getPrimerTrimmedFlowgram("")
        self.assertEqual(str(flow1_untrimmed), "1.06	0.08	1.04	0.08	0.05	0.94	0.10	2.01	0.10	0.07	0.96	0.09	1.04	1.96	1.07	0.10	1.01	0.13	0.08	1.01	1.06	1.83	2.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04")
        self.assertEqual(flow1_untrimmed.Bases,	"tcagGCTAACTGTAACCCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc")

        flow2_trimmed = flow2.getPrimerTrimmedFlowgram(primerseq="TCAG"+"AGACGCACT")
        self.assertEqual(str(flow2_trimmed), "0.00\t0.05\t0.90	0.11	0.07	1.99	0.11	0.02	1.96	1.04	0.13	0.01	2.83	0.10	1.97	0.06	0.11	1.04	0.13	0.03	0.98	1.15	0.07	1.00	0.07	0.08	0.98	0.11	1.92	0.05	0.04	2.96	1.02	1.02	0.04	0.93	1.00	0.13	0.04	1.00	1.03	0.08	0.97	0.13	0.11	1.88	0.09	0.05	1.02	1.89	0.07	0.11	0.98	0.05	0.07	1.01	0.08	0.05	1.01	0.13	1.00	0.07	0.10	1.04	0.10	0.04	0.98	0.12	1.03	0.96	0.11	0.07	1.00	0.09	0.03	1.03	0.11	1.95	1.06	0.13	0.05	1.00	0.13	0.11	1.00	0.09	0.03	2.89	0.08	0.95	0.09	1.03	1.02	1.05	1.07	0.08	0.12	2.81	0.08	0.08	1.00	1.07	0.07	0.05	1.86	0.12	0.98	0.06	2.00	0.11	1.02	0.11	0.08	1.88	0.13	1.03	0.13	0.98	0.15	0.11	1.03	1.03	1.04	0.18	0.98	0.13	0.15	1.04	0.11	1.01	0.13	0.06	1.01	0.06	1.02	0.08	0.99	0.14	0.99	0.09	0.05	1.09	0.04	0.07	2.96	0.09	2.03	0.13	2.96	1.13	0.08	1.03	0.07	0.99	0.11	0.05	1.05	1.04	0.09	0.07	1.00	1.03	0.09	0.06	1.06	1.04	2.94	0.18	0.06	0.93	0.10	1.10	0.11	2.02	0.17	1.00	1.03	0.06	0.11	0.96	0.04	3.00	0.11	0.07	1.99	0.10	2.03	0.12	0.97	0.16	0.01	2.09	0.14	1.04	0.16	0.06	1.03	0.14	1.12	0.12	0.05	0.96	1.01	0.10	0.14	0.94	0.03	0.12	1.10	0.92	0.09	1.10	1.04	1.02	0.12	0.97	2.00	0.15	1.08	0.04	1.03	1.04	0.03	0.09	5.16	1.02	0.09	0.13	2.66	0.09	0.05	1.06	0.07	0.89	0.05	0.12	1.10	0.16	0.06	1.01	0.13	1.00	0.14	0.98	0.09	2.92	1.28	0.03	2.95	0.98	0.16	0.08	0.95	0.96	1.09	0.08	1.07	1.01	0.16	0.06	4.52	0.12	1.03	0.07	0.09	1.03	0.14	0.03	1.01	1.99	1.05	0.14	1.03	0.13	0.03	1.10	0.10	0.96	0.11	0.99	0.12	0.05	0.94	2.83	0.14	0.12	0.96	0.00	1.00	0.11	0.14	1.98	0.08	0.11	1.04	0.01	0.11	2.03	0.15	2.05	0.10	0.03	0.93	0.01	0.08	0.12	0.00	0.16	0.05	0.07	0.08	0.11	0.07	0.05	0.04	0.10	0.05	0.05	0.03	0.07	0.03	0.04	0.04	0.06	0.03	0.05	0.04	0.09	0.03	0.08	0.03	0.07	0.02	0.05	0.02	0.06	0.01	0.05	0.04	0.06	0.02	0.04	0.04	0.04	0.03	0.03	0.06	0.06	0.03	0.02	0.02	0.08	0.03	0.01	0.01	0.06	0.03	0.01	0.03	0.04	0.02	0.00	0.02	0.05	0.00	0.02	0.02	0.03	0.00	0.02	0.02	0.04	0.01	0.00	0.01	0.05")
        self.assertEqual(flow2_trimmed.Bases,
                         "CAATTATTTCCATAGCTTGGGTAGTGTCAATAATGCTGCTATGAACATGGGAGTACAAATATTCTTCAAGATACTGATCTCATTTCCTTTAGATATATACCCAGAAGTGAAATTCCTGGATCACATAGTAGTTCTATTTTTATTTGATGAGAAACTTTATACTATTTTTCATAActgagcgggctggcaaggc")

        #trimming at the end of the flow cycle works
        flow3_trimmed = flow3.getPrimerTrimmedFlowgram(primerseq="TCAG"+"ATTAGATACCCNGGTAGG")
        self.assertEqual(str(flow3_trimmed), "0.05	0.05	2.04	0.10	0.03	1.06	1.05	1.01	0.07	0.09	2.07	1.01	0.93	2.88	1.06	1.95	1.00	0.05	0.05	2.97	0.09	0.00	0.93	1.01	0.06	0.05	0.99	0.09	0.98	1.01	0.03	1.02	1.92	0.07	0.01	1.03	1.01	0.01	0.05	0.96	0.09	0.05	0.98	1.07	0.02	2.02	2.05	0.09	1.87	0.12	2.15	0.05	0.13	0.92	1.05	1.96	3.01	0.13	0.04	1.05	0.96	0.05	0.05	0.95	0.12	0.01	1.00	2.02	0.03	0.03	0.99	1.01	0.05	0.06	0.98	0.13	0.06	0.97	0.11	1.01	0.08	0.12	1.02	0.12	1.02	2.19	1.03	1.01	0.08	0.11	0.96	0.09	0.08	1.01	0.08	0.06	2.10	2.11	0.12	1.04	0.13	0.09	0.94	1.03	0.08	0.05	3.06	0.12	1.00	0.03	0.09	0.95	0.10	0.03	2.09	0.21	0.99	0.06	0.11	4.06	0.10	1.04	0.04	1.05	1.05	1.04	1.02	0.97	0.13	0.93	0.10	0.12	1.08	0.12	0.99	1.06	0.10	0.11	0.98	0.10	0.02	2.01	0.10	1.01	0.09	0.96	0.07	0.11	2.03	4.12	1.05	0.08	1.01	0.04	0.98	0.14	0.12	2.96	0.13	1.98	0.12	2.08	0.10	0.12	1.99	0.13	0.07	0.98	0.03	0.93	0.86	4.10	0.13	0.10	3.99	1.13	0.07	0.06	1.07	0.09	0.05	1.03	1.12	0.13	0.05	2.01	0.08	0.80	0.05	0.11	0.98	0.13	0.04	1.01	0.07	1.02	0.07	0.11	1.07	2.19	0.06	0.97	0.11	1.03	0.05	0.11	1.05	0.14	0.06	1.03	0.13	0.10	0.97	0.16	0.13	1.00	0.13	0.06	1.02	2.15	0.02	0.16	0.95	0.09	2.06	2.12	0.07	0.07	2.08	0.12	0.97	1.00	0.03	0.99	1.02	1.01	0.03	0.15	0.90	0.07	0.01	2.00	1.01	1.00	0.06	0.11	1.08	1.00	0.03	1.99	0.03	1.00	0.02	1.85	1.93	0.14	1.97	0.91	1.83	0.06	0.04	1.97	0.05	2.08	0.04	0.06	1.05	0.05	2.13	0.16	0.09	1.17	0.01	1.01	1.07	0.09	0.14	0.91	0.06	0.08	1.03	1.04	0.08	0.05	1.05	1.03	1.16	0.06	0.05	1.01	0.06	2.15	0.06	1.99	0.13	0.04	1.08	0.97	0.11	0.07	1.05	0.08	0.07	2.13	0.14	0.09	1.10	0.15	0.00	1.02	0.07	1.05	0.05	0.95	0.09	1.00	0.15	0.95	0.08	0.15	1.11	0.07	0.12	1.05	1.06	0.09	1.03	0.07	0.11	1.01	0.05	0.05	1.05	0.98	0.00	0.93	0.08	0.12	1.85	1.11	0.10	0.07	1.00	0.01	0.10	1.87	0.05	2.14	1.10	0.03	1.06	0.10	0.91	0.10	0.06	1.05	1.02	1.02	0.07	0.06	0.98	0.95	1.09	0.06	0.14	0.97	0.04	2.44")
        self.assertEqual(flow3_trimmed.Bases,
                         "CCACGCCGTAAACGGTGGGCGCTAGTTGTGCGAACCTTCCACGGTTTGTGCGGCGCAGCTAACGCATTAAGCGCCCTGCCTGGGGAGTACGATCGCAAGATTAAAACTCAAAGGAATTGACGGGGCCCCGCACAAGCAGCGGAGCATGCGGCTTAATTCGACGCAACGCGAAGAACCTTACCAAGGCTTGACATATACAGGAATATGGCAGAGATGTCATAGCCGCAAGGTCTGTATACAGG")

        flow3_trimmed = flow3.getPrimerTrimmedFlowgram(primerseq="TCAG"+"ATTAGATACCCNGGTAG")
        self.assertEqual(str(flow3_trimmed), "0.00\t0.00\t0.00	1.10	0.05	0.05	2.04	0.10	0.03	1.06	1.05	1.01	0.07	0.09	2.07	1.01	0.93	2.88	1.06	1.95	1.00	0.05	0.05	2.97	0.09	0.00	0.93	1.01	0.06	0.05	0.99	0.09	0.98	1.01	0.03	1.02	1.92	0.07	0.01	1.03	1.01	0.01	0.05	0.96	0.09	0.05	0.98	1.07	0.02	2.02	2.05	0.09	1.87	0.12	2.15	0.05	0.13	0.92	1.05	1.96	3.01	0.13	0.04	1.05	0.96	0.05	0.05	0.95	0.12	0.01	1.00	2.02	0.03	0.03	0.99	1.01	0.05	0.06	0.98	0.13	0.06	0.97	0.11	1.01	0.08	0.12	1.02	0.12	1.02	2.19	1.03	1.01	0.08	0.11	0.96	0.09	0.08	1.01	0.08	0.06	2.10	2.11	0.12	1.04	0.13	0.09	0.94	1.03	0.08	0.05	3.06	0.12	1.00	0.03	0.09	0.95	0.10	0.03	2.09	0.21	0.99	0.06	0.11	4.06	0.10	1.04	0.04	1.05	1.05	1.04	1.02	0.97	0.13	0.93	0.10	0.12	1.08	0.12	0.99	1.06	0.10	0.11	0.98	0.10	0.02	2.01	0.10	1.01	0.09	0.96	0.07	0.11	2.03	4.12	1.05	0.08	1.01	0.04	0.98	0.14	0.12	2.96	0.13	1.98	0.12	2.08	0.10	0.12	1.99	0.13	0.07	0.98	0.03	0.93	0.86	4.10	0.13	0.10	3.99	1.13	0.07	0.06	1.07	0.09	0.05	1.03	1.12	0.13	0.05	2.01	0.08	0.80	0.05	0.11	0.98	0.13	0.04	1.01	0.07	1.02	0.07	0.11	1.07	2.19	0.06	0.97	0.11	1.03	0.05	0.11	1.05	0.14	0.06	1.03	0.13	0.10	0.97	0.16	0.13	1.00	0.13	0.06	1.02	2.15	0.02	0.16	0.95	0.09	2.06	2.12	0.07	0.07	2.08	0.12	0.97	1.00	0.03	0.99	1.02	1.01	0.03	0.15	0.90	0.07	0.01	2.00	1.01	1.00	0.06	0.11	1.08	1.00	0.03	1.99	0.03	1.00	0.02	1.85	1.93	0.14	1.97	0.91	1.83	0.06	0.04	1.97	0.05	2.08	0.04	0.06	1.05	0.05	2.13	0.16	0.09	1.17	0.01	1.01	1.07	0.09	0.14	0.91	0.06	0.08	1.03	1.04	0.08	0.05	1.05	1.03	1.16	0.06	0.05	1.01	0.06	2.15	0.06	1.99	0.13	0.04	1.08	0.97	0.11	0.07	1.05	0.08	0.07	2.13	0.14	0.09	1.10	0.15	0.00	1.02	0.07	1.05	0.05	0.95	0.09	1.00	0.15	0.95	0.08	0.15	1.11	0.07	0.12	1.05	1.06	0.09	1.03	0.07	0.11	1.01	0.05	0.05	1.05	0.98	0.00	0.93	0.08	0.12	1.85	1.11	0.10	0.07	1.00	0.01	0.10	1.87	0.05	2.14	1.10	0.03	1.06	0.10	0.91	0.10	0.06	1.05	1.02	1.02	0.07	0.06	0.98	0.95	1.09	0.06	0.14	0.97	0.04	2.44")
        self.assertEqual(flow3_trimmed.Bases,
                         "GCCACGCCGTAAACGGTGGGCGCTAGTTGTGCGAACCTTCCACGGTTTGTGCGGCGCAGCTAACGCATTAAGCGCCCTGCCTGGGGAGTACGATCGCAAGATTAAAACTCAAAGGAATTGACGGGGCCCCGCACAAGCAGCGGAGCATGCGGCTTAATTCGACGCAACGCGAAGAACCTTACCAAGGCTTGACATATACAGGAATATGGCAGAGATGTCATAGCCGCAAGGTCTGTATACAGG")

    def test_createFlowHeader(self):
        """header_info dict turned into flowgram header"""
        f = Flowgram('0.5 1.0 4.0 0.0 1.5 0.0 0.0 2.0', Name='a',
                      header_info = {'Bases':'TACCCCTTGG','Name Length':'14'})

        self.assertEqual(f.createFlowHeader(),
                         """>a\n  Name Length:\t14\nBases:\tTACCCCTTGG\nFlowgram:\t0.5\t1.0\t4.0\t0.0\t1.5\t0.0\t0.0\t2.0\n""")

    def test_build_averaged_flowgram(self):
        f1 = [0.3, 1.1, 4.0 , 0.01, 0.8, 0.0, 0.0, 2.0]
        f2 = [0.6, 0.9, 4.05, 0.1,  1.2, 0.1, 0.4]
        f3 = [0.4, 1.2, 4.05, 0.2,  1.3, 0.2]
        f4 = [0.7, 1.0, 4.0 , 0.02, 1.5]
        flowgrams = [f1,f2,f3,f4]
   
        self.assertFloatEqual(build_averaged_flowgram(flowgrams),
                         [0.5, 1.05, 4.03, 0.08, 1.2, 0.1, 0.2, 2.0])
        self.assertFloatEqual(build_averaged_flowgram([f1,f1,f1,f1,f1,f1]),
                         [0.3, 1.1, 4.0 , 0.01, 0.8, 0.0, 0.0, 2.0])


    def setUp(self):
        """Define some standard data"""
        self.rec = """Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  96099976
  Index Length:  1158685
  # of Reads:    57902
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FIQU8OX05GCVRO
  Run Prefix:   R_2008_10_15_16_11_02_
  Region #:     5
  XY Location:  2489_3906

  Run Name:       R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford
  Analysis Name:  /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis
  Full Path:      /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis

  Read Header Len:  32
  Name Length:      14
  # of Bases:       104
  Clip Qual Left:   5
  Clip Qual Right:  85
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.06	0.08	1.04	0.08	0.05	0.94	0.10	2.01	0.10	0.07	0.96	0.09	1.04	1.96	1.07	0.10	1.01	0.13	0.08	1.01	1.06	1.83	2.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04
Flow Indexes:	1	3	6	8	8	11	13	14	14	15	17	20	21	22	22	23	23	23	25	27	29	29	32	32	35	38	39	39	39	42	43	45	46	46	46	47	48	51	51	54	54	57	59	61	61	64	67	69	72	72	74	76	77	80	81	81	81	82	83	83	86	88	88	91	94	95	95	95	98	100	103	106	106	109	112	113	116	118	118	121	122	124	125	127	130	131	133	136	138	140	143	144	144	144	147	149	152	152	155	158	158	160	160	163
Bases:	tcagGCTAACTGTAACCCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	37	37	37	37	37	39	39	39	39	24	24	24	37	34	28	24	24	24	28	34	39	39	39	39	39	39	39	39	39	39	39	39	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37

>FIQU8OX05F8ILF
  Run Prefix:   R_2008_10_15_16_11_02_
  Region #:     5
  XY Location:  2440_0913

  Run Name:       R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford
  Analysis Name:  /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis
  Full Path:      /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis

  Read Header Len:  32
  Name Length:      14
  # of Bases:       206
  Clip Qual Left:   5
  Clip Qual Right:  187
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.00	1.01	0.00	0.00	1.00	0.00	1.00	0.00	1.05	0.00	0.91	0.10	1.07	0.95	1.01	0.00	0.06	0.93	0.02	0.03	1.06	1.18	0.09	1.00	0.05	0.90	0.11	0.07	1.99	0.11	0.02	1.96	1.04	0.13	0.01	2.83	0.10	1.97	0.06	0.11	1.04	0.13	0.03	0.98	1.15	0.07	1.00	0.07	0.08	0.98	0.11	1.92	0.05	0.04	2.96	1.02	1.02	0.04	0.93	1.00	0.13	0.04	1.00	1.03	0.08	0.97	0.13	0.11	1.88	0.09	0.05	1.02	1.89	0.07	0.11	0.98	0.05	0.07	1.01	0.08	0.05	1.01	0.13	1.00	0.07	0.10	1.04	0.10	0.04	0.98	0.12	1.03	0.96	0.11	0.07	1.00	0.09	0.03	1.03	0.11	1.95	1.06	0.13	0.05	1.00	0.13	0.11	1.00	0.09	0.03	2.89	0.08	0.95	0.09	1.03	1.02	1.05	1.07	0.08	0.12	2.81	0.08	0.08	1.00	1.07	0.07	0.05	1.86	0.12	0.98	0.06	2.00	0.11	1.02	0.11	0.08	1.88	0.13	1.03	0.13	0.98	0.15	0.11	1.03	1.03	1.04	0.18	0.98	0.13	0.15	1.04	0.11	1.01	0.13	0.06	1.01	0.06	1.02	0.08	0.99	0.14	0.99	0.09	0.05	1.09	0.04	0.07	2.96	0.09	2.03	0.13	2.96	1.13	0.08	1.03	0.07	0.99	0.11	0.05	1.05	1.04	0.09	0.07	1.00	1.03	0.09	0.06	1.06	1.04	2.94	0.18	0.06	0.93	0.10	1.10	0.11	2.02	0.17	1.00	1.03	0.06	0.11	0.96	0.04	3.00	0.11	0.07	1.99	0.10	2.03	0.12	0.97	0.16	0.01	2.09	0.14	1.04	0.16	0.06	1.03	0.14	1.12	0.12	0.05	0.96	1.01	0.10	0.14	0.94	0.03	0.12	1.10	0.92	0.09	1.10	1.04	1.02	0.12	0.97	2.00	0.15	1.08	0.04	1.03	1.04	0.03	0.09	5.16	1.02	0.09	0.13	2.66	0.09	0.05	1.06	0.07	0.89	0.05	0.12	1.10	0.16	0.06	1.01	0.13	1.00	0.14	0.98	0.09	2.92	1.28	0.03	2.95	0.98	0.16	0.08	0.95	0.96	1.09	0.08	1.07	1.01	0.16	0.06	4.52	0.12	1.03	0.07	0.09	1.03	0.14	0.03	1.01	1.99	1.05	0.14	1.03	0.13	0.03	1.10	0.10	0.96	0.11	0.99	0.12	0.05	0.94	2.83	0.14	0.12	0.96	0.00	1.00	0.11	0.14	1.98	0.08	0.11	1.04	0.01	0.11	2.03	0.15	2.05	0.10	0.03	0.93	0.01	0.08	0.12	0.00	0.16	0.05	0.07	0.08	0.11	0.07	0.05	0.04	0.10	0.05	0.05	0.03	0.07	0.03	0.04	0.04	0.06	0.03	0.05	0.04	0.09	0.03	0.08	0.03	0.07	0.02	0.05	0.02	0.06	0.01	0.05	0.04	0.06	0.02	0.04	0.04	0.04	0.03	0.03	0.06	0.06	0.03	0.02	0.02	0.08	0.03	0.01	0.01	0.06	0.03	0.01	0.03	0.04	0.02	0.00	0.02	0.05	0.00	0.02	0.02	0.03	0.00	0.02	0.02	0.04	0.01	0.00	0.01	0.05
Flow Indexes:	1	3	6	8	10	12	14	15	16	19	22	23	25	27	30	30	33	33	34	37	37	37	39	39	42	45	46	48	51	53	53	56	56	56	57	58	60	61	64	65	67	70	70	73	74	74	77	80	83	85	88	91	93	94	97	100	102	102	103	106	109	112	112	112	114	116	117	118	119	122	122	122	125	126	129	129	131	133	133	135	138	138	140	142	145	146	147	149	152	154	157	159	161	163	166	169	169	169	171	171	173	173	173	174	176	178	181	182	185	186	189	190	191	191	191	194	196	198	198	200	201	204	206	206	206	209	209	211	211	213	216	216	218	221	223	226	227	230	233	234	236	237	238	240	241	241	243	245	246	249	249	249	249	249	250	253	253	253	256	258	261	264	266	268	270	270	270	271	273	273	273	274	277	278	279	281	282	285	285	285	285	285	287	290	293	294	294	295	297	300	302	304	307	308	308	308	311	313	316	316	319	322	322	324	324	327
Bases:	tcagAGACGCACTCAATTATTTCCATAGCTTGGGTAGTGTCAATAATGCTGCTATGAACATGGGAGTACAAATATTCTTCAAGATACTGATCTCATTTCCTTTAGATATATACCCAGAAGTGAAATTCCTGGATCACATAGTAGTTCTATTTTTATTTGATGAGAAACTTTATACTATTTTTCATAActgagcgggctggcaaggc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	34	34	34	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	36	36	38	25	25	25	38	37	37	37	37	37	37	33	33	34	37	37	37	37	37	37	37	38	34	20	20	26	26	20	34	38	37	37	37	37	37	37	37	37	37	38	38	38	37	37	37	37	37	37	37	37	37	37

>FIQU8OX06G9PCS
  Run Prefix:   R_2008_10_15_16_11_02_
  Region #:     6
  XY Location:  2863_3338

  Run Name:       R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford
  Analysis Name:  /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis
  Full Path:      /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis

  Read Header Len:  32
  Name Length:      14
  # of Bases:       264
  Clip Qual Left:   5
  Clip Qual Right:  264
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.05	1.01	0.07	0.05	0.99	0.03	1.05	0.04	1.05	0.05	0.06	2.05	1.13	0.03	1.00	0.08	1.07	0.09	0.05	1.02	1.11	3.06	0.09	0.04	1.03	0.13	1.97	1.02	1.07	0.06	2.10	0.05	0.05	2.04	0.10	0.03	1.06	1.05	1.01	0.07	0.09	2.07	1.01	0.93	2.88	1.06	1.95	1.00	0.05	0.05	2.97	0.09	0.00	0.93	1.01	0.06	0.05	0.99	0.09	0.98	1.01	0.03	1.02	1.92	0.07	0.01	1.03	1.01	0.01	0.05	0.96	0.09	0.05	0.98	1.07	0.02	2.02	2.05	0.09	1.87	0.12	2.15	0.05	0.13	0.92	1.05	1.96	3.01	0.13	0.04	1.05	0.96	0.05	0.05	0.95	0.12	0.01	1.00	2.02	0.03	0.03	0.99	1.01	0.05	0.06	0.98	0.13	0.06	0.97	0.11	1.01	0.08	0.12	1.02	0.12	1.02	2.19	1.03	1.01	0.08	0.11	0.96	0.09	0.08	1.01	0.08	0.06	2.10	2.11	0.12	1.04	0.13	0.09	0.94	1.03	0.08	0.05	3.06	0.12	1.00	0.03	0.09	0.95	0.10	0.03	2.09	0.21	0.99	0.06	0.11	4.06	0.10	1.04	0.04	1.05	1.05	1.04	1.02	0.97	0.13	0.93	0.10	0.12	1.08	0.12	0.99	1.06	0.10	0.11	0.98	0.10	0.02	2.01	0.10	1.01	0.09	0.96	0.07	0.11	2.03	4.12	1.05	0.08	1.01	0.04	0.98	0.14	0.12	2.96	0.13	1.98	0.12	2.08	0.10	0.12	1.99	0.13	0.07	0.98	0.03	0.93	0.86	4.10	0.13	0.10	3.99	1.13	0.07	0.06	1.07	0.09	0.05	1.03	1.12	0.13	0.05	2.01	0.08	0.80	0.05	0.11	0.98	0.13	0.04	1.01	0.07	1.02	0.07	0.11	1.07	2.19	0.06	0.97	0.11	1.03	0.05	0.11	1.05	0.14	0.06	1.03	0.13	0.10	0.97	0.16	0.13	1.00	0.13	0.06	1.02	2.15	0.02	0.16	0.95	0.09	2.06	2.12	0.07	0.07	2.08	0.12	0.97	1.00	0.03	0.99	1.02	1.01	0.03	0.15	0.90	0.07	0.01	2.00	1.01	1.00	0.06	0.11	1.08	1.00	0.03	1.99	0.03	1.00	0.02	1.85	1.93	0.14	1.97	0.91	1.83	0.06	0.04	1.97	0.05	2.08	0.04	0.06	1.05	0.05	2.13	0.16	0.09	1.17	0.01	1.01	1.07	0.09	0.14	0.91	0.06	0.08	1.03	1.04	0.08	0.05	1.05	1.03	1.16	0.06	0.05	1.01	0.06	2.15	0.06	1.99	0.13	0.04	1.08	0.97	0.11	0.07	1.05	0.08	0.07	2.13	0.14	0.09	1.10	0.15	0.00	1.02	0.07	1.05	0.05	0.95	0.09	1.00	0.15	0.95	0.08	0.15	1.11	0.07	0.12	1.05	1.06	0.09	1.03	0.07	0.11	1.01	0.05	0.05	1.05	0.98	0.00	0.93	0.08	0.12	1.85	1.11	0.10	0.07	1.00	0.01	0.10	1.87	0.05	2.14	1.10	0.03	1.06	0.10	0.91	0.10	0.06	1.05	1.02	1.02	0.07	0.06	0.98	0.95	1.09	0.06	0.14	0.97	0.04	2.44
Flow Indexes:	1	3	6	8	10	13	13	14	16	18	21	22	23	23	23	26	28	28	29	30	32	32	35	35	38	39	40	43	43	44	45	46	46	46	47	48	48	49	52	52	52	55	56	59	61	62	64	65	65	68	69	72	75	76	78	78	79	79	81	81	83	83	86	87	88	88	89	89	89	92	93	96	99	100	100	103	104	107	110	112	115	117	118	118	119	120	123	126	129	129	130	130	132	135	136	139	139	139	141	144	147	147	149	152	152	152	152	154	156	157	158	159	160	162	165	167	168	171	174	174	176	178	181	181	182	182	182	182	183	185	187	190	190	190	192	192	194	194	197	197	200	202	203	204	204	204	204	207	207	207	207	208	211	214	215	218	218	220	223	226	228	231	232	232	234	236	239	242	245	248	251	252	252	255	257	257	258	258	261	261	263	264	266	267	268	271	274	274	275	276	279	280	282	282	284	286	286	287	287	289	289	290	291	291	294	294	296	296	299	301	301	304	306	307	310	313	314	317	318	319	322	324	324	326	326	329	330	333	336	336	339	342	344	346	348	350	353	356	357	359	362	365	366	368	371	371	372	375	378	378	380	380	381	383	385	388	389	390	393	394	395	398	400	400
Bases:	tcagATTAGATACCCAGGTAGGCCACGCCGTAAACGGTGGGCGCTAGTTGTGCGAACCTTCCACGGTTTGTGCGGCGCAGCTAACGCATTAAGCGCCCTGCCTGGGGAGTACGATCGCAAGATTAAAACTCAAAGGAATTGACGGGGCCCCGCACAAGCAGCGGAGCATGCGGCTTAATTCGACGCAACGCGAAGAACCTTACCAAGGCTTGACATATACAGGAATATGGCAGAGATGTCATAGCCGCAAGGTCTGTATACAGG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	40	40	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	37	37	37	37	37	37	37	37	38	38	40	40	40	40	40	38	38	38	38	38	40	40	38	38	38	38	38	40	40	40	40	38	38	38	38	38	38	31	30	30	30	32	31	32	31	32	31	31	28	25	21	20

""".split('\n')
        flows, head = parse_sff(self.rec)
        self.flows = list(flows)

if __name__ == "__main__":
    main()

