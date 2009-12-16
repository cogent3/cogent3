#!/usr/bin/env python
"""tests for sff parser"""

__author__ = "Julia Goodrich, Jens Reeder"
__copyright__ = "Copyright 2009, The Cogent Project"
__credits__ = ["Julia Goodrich","Jens Reeder"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jens Reeder"
__email__ = "jreeder@colorado.edu"
__status__ = "Development"

from types import GeneratorType
from cogent.util.unit_test import TestCase, main
from cogent.parse.flowgram_parser import get_header_info, get_summaries,\
    get_all_summaries, split_summary, parse_sff, lazy_parse_sff_handle

class SFFParserTests(TestCase):
    """Tests sff parser functions"""

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

""".split('\n')

    def test_get_header_info(self):
        """get_header_info should return a sff file common header as a dict"""
        header = get_header_info(self.rec)
        self.assertEqual(len(header), 11)
        self.assertEqual(header['Key Length'], '4')
        self.assertEqual(header['Key Sequence'], 'TCAG')

    def test_get_summaries(self):
        """get_summaries should return a generator of the summaries"""
        summaries = get_summaries(self.rec,number_list = [1])
        sum_list = list(summaries)
        self.assertEqual(len(sum_list), 1)
        self.assertEqual(isinstance(summaries, GeneratorType), True)
        self.assertEqual(len(sum_list[0]), 18)
        self.assertEqual(sum_list[0][0], '>FIQU8OX05F8ILF')
        
        summaries = get_summaries(self.rec,name_list = ['FIQU8OX05GCVRO'])
        sum_list = list(summaries)
        self.assertEqual(len(sum_list), 1)
        self.assertEqual(isinstance(summaries, GeneratorType), True)
        self.assertEqual(len(sum_list[0]), 18)
        self.assertEqual(sum_list[0][0], '>FIQU8OX05GCVRO')
        
        summaries = get_summaries(self.rec,all_sums = True )
        sum_list = list(summaries)
        self.assertEqual(len(sum_list), 2)
        self.assertEqual(isinstance(summaries, GeneratorType), True)
        self.assertEqual(len(sum_list[0]), 18)
        self.assertEqual(sum_list[0][0], '>FIQU8OX05GCVRO')
        self.assertEqual(sum_list[1][0], '>FIQU8OX05F8ILF')
        
        summaries = get_summaries(self.rec,number_list = [0],
                                 name_list =['FIQU8OX05GCVRO'])
        self.assertRaises(AssertionError,list,summaries)

        summaries = get_summaries(self.rec)
        self.assertRaises(ValueError,list, summaries)
        
    def test_get_all_summaries(self):
        """get_all_summaries should return a list of the summaries"""
        summaries = get_all_summaries(self.rec)
        self.assertEqual(len(summaries), 2)
        self.assertEqual(isinstance(summaries,list), True)
        self.assertEqual(len(summaries[0]), 18)
        self.assertEqual(summaries[0][0], '>FIQU8OX05GCVRO')
        self.assertEqual(summaries[1][0], '>FIQU8OX05F8ILF')

    def test_split_summary(self):
        """split_summary should return the info of a flowgram header."""
        summaries = get_all_summaries(self.rec)
        sum_dict = split_summary(summaries[0])
        self.assertEqual(len(sum_dict), 18)
        self.assertEqual(sum_dict['Name'], 'FIQU8OX05GCVRO')
        assert 'Flowgram' in sum_dict
        assert 'Bases' in sum_dict

        sum_dict = split_summary(summaries[1])
        self.assertEqual(len(sum_dict), 18)
        self.assertEqual(sum_dict['Name'], 'FIQU8OX05F8ILF')
        assert 'Flowgram' in sum_dict
        assert 'Bases' in sum_dict
        
    def test_parse_sff(self):
        """SFParser should read in the SFF file correctly."""
        flows, head = parse_sff(self.rec)
        self.assertEqual(len(flows),2)
        self.assertEqual(len(head), 11)
        self.assertEqual(head['Key Length'], '4')
        self.assertEqual(head['Key Sequence'], 'TCAG')
        self.assertEqual(flows[0].Name, 'FIQU8OX05GCVRO')
        self.assertEqual(flows[1].Name, 'FIQU8OX05F8ILF')

    def test_lazy_parse_sff_handle(self):
        """LazySFParser should read in the SFF file correctly."""
        flows, head = lazy_parse_sff_handle(self.rec)
        flows = list(flows)
        self.assertEqual(len(flows),2)
        self.assertEqual(len(head), 11)
        self.assertEqual(head['Key Length'], '4')
        self.assertEqual(head['Key Sequence'], 'TCAG')
        self.assertEqual(flows[0].Name, 'FIQU8OX05GCVRO')
        self.assertEqual(flows[1].Name, 'FIQU8OX05F8ILF')
         
                            
if __name__ == "__main__":
    main()

