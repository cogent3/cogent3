#!/usr/bin/env python
# test_sffinfo.py

import os
import shutil
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError
from cogent.app.sffinfo import (
    ManyValuedParameter, Sffinfo, sffinfo_from_file)

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2010, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.4.1"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class ManyValuedParameterTests(TestCase):
    def test_init(self):
        """__init__() should set appropriate class variables."""
        p = ManyValuedParameter(None, None)
        self.assertEqual(p.Name, None)
        self.assertEqual(p.Value, None)
        self.assertEqual(p.Delimiter, None)
        self.assertEqual(p.Quote, None)

        p = ManyValuedParameter('-', 'a', Values=['abc'])
        self.assertEqual(p.Value, ['abc'])

    def test_append(self):
        """append() should append values to Value class attribute."""
        p = ManyValuedParameter(None, None)
        p.append('abc')
        p.append(3)
        self.assertEqual(p.Value, ['abc', 3])

    def test_on(self):
        """on() should alias append()."""
        p = ManyValuedParameter(None, None)
        p.on('abc')
        p.on(3)
        self.assertEqual(p.Value, ['abc', 3])

    def test_str(self):
        """__str__() should produce quoted, delimited string of parameter values."""
        p = ManyValuedParameter(None, None)
        p.append('abc')
        p.append(3)
        self.assertEqual(str(p), 'abc3')

        p = ManyValuedParameter(None, None, Quote='"', ValueDelimiter=',')
        p.append('abc')
        p.append(3)
        self.assertEqual(str(p), '"abc","3"')


class SffinfoTests(TestCase):
    def setUp(self):
        test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.sff_fp = os.path.join(test_dir, 'data', 'test.sff')

    def test_base_command(self):
        """BaseCommand should include non-positional parameters."""
        s = Sffinfo()
        expected = 'cd "%s/"; sffinfo' % os.getcwd()
        self.assertEqual(s.BaseCommand, expected)

        s.Parameters['-a'].on()
        expected = 'cd "%s/"; sffinfo -a' % os.getcwd()
        self.assertEqual(s.BaseCommand, expected)

        # accession number parameters should not be included in base command
        s.Parameters['-a'].off()
        s.Parameters['accno'].on('12345ABC')
        expected = 'cd "%s/"; sffinfo' % os.getcwd()
        self.assertEqual(s.BaseCommand, expected)

    def test_changing_working_dir(self):
        """WorkingDir should be created and included in BaseCommand."""
        # Directory is not created, only filename is returned
        working_dir = tempfile.mktemp()
        self.assertRaises(OSError, os.rmdir, working_dir)
        
        a = Sffinfo(WorkingDir=working_dir)
        expected = 'cd "%s/"; sffinfo' % working_dir
        self.assertEqual(a.BaseCommand, expected)

        # Removing the directory is proof that it was created.  If the
        # directory is not there, an OSError will be raised.
        os.rmdir(working_dir)

    def test_input_handler(self):
        """Sffinfo should decorate input handler output with accession numbers"""
        my_accno = '12345ABC'
        a = Sffinfo()
        a.Parameters['accno'].on(my_accno)
        self.assertEqual(a.InputHandler, '_input_handler_decorator')
        observed = getattr(a, a.InputHandler)(self.sff_fp)
        expected = '"%s" %s' % (self.sff_fp, my_accno)
        self.assertEqual(observed, expected)

    def test_call(self):
        """Simple sffinfo call should produce expected output."""
        a = Sffinfo()
        app_results = a(self.sff_fp)
        observed = app_results['StdOut'].read()
        self.assertEqual(observed, sffinfo_output)

    def test_call_with_accno(self):
        """Sffinfo accession number parameters should filter output."""
        # Valid accession number
        a = Sffinfo()
        a.Parameters['accno'].on('FA6P1OK01CGMHQ')
        app_results = a(self.sff_fp)
        observed = app_results['StdOut'].read()
        self.assertEqual(observed, sffinfo_output)

        # Invalid accession number
        a = Sffinfo()
        a.Parameters['accno'].on('AAAAAAAAAAAAAA')
        app_results = a(self.sff_fp)
        observed = app_results['StdOut'].read()
        self.assertEqual(observed, empty_sffinfo_output)

    def test_call_with_flags(self):
        """Sffinfo flags should alter output as expected."""
        # -a flag
        a = Sffinfo()
        a.Parameters['-a'].on()
        app_results = a(self.sff_fp)
        observed = app_results['StdOut'].read()
        self.assertEqual(observed, 'FA6P1OK01CGMHQ\n')

        # -s flag
        a = Sffinfo()
        a.Parameters['-s'].on()
        app_results = a(self.sff_fp)
        observed = app_results['StdOut'].read()
        expected = (
            '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 '
            'run=R_2008_05_28_17_11_38_\n'
            'ATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n'
            )
        self.assertEqual(observed, expected)

        # -q flag
        a = Sffinfo()
        a.Parameters['-q'].on()
        app_results = a(self.sff_fp)
        observed = app_results['StdOut'].read()
        expected = (
            '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 '
            'run=R_2008_05_28_17_11_38_\n'
            '32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 '
            '33 31 25 18 18 20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 '
            '17 18 13 17 27 21\n')
        self.assertEqual(observed, expected)

        # -f flag
        a = Sffinfo()
        a.Parameters['-f'].on()
        app_results = a(self.sff_fp)
        observed = app_results['StdOut'].read()
        expected = (
            '>FA6P1OK01CGMHQ xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n'
            'A,1.02 C,0.00 G,0.00 T,0.99 A,0.00 C,1.00 G,0.00 T,1.00 A,0.00 C,0.00\n'
            'G,1.00 T,0.00 A,1.10 C,0.00 G,1.08 T,0.00 A,0.00 C,1.46 G,0.00 T,0.88\n'
            'A,0.18 C,0.00 G,2.69 T,1.01 A,0.08 C,0.96 G,0.00 T,0.02 A,0.92 C,0.08\n'
            'G,0.00 T,0.98 A,0.68 C,0.00 G,0.89 T,0.00 A,0.00 C,1.15 G,0.00 T,1.13\n'
            'A,0.00 C,0.02 G,1.12 T,0.05 A,0.15 C,1.84 G,0.00 T,1.10 A,0.00 C,2.47\n'
            'G,0.96 T,0.86 A,1.06 C,0.00 G,1.96 T,0.12 A,0.93 C,0.13 G,1.65 T,1.06\n'
            'A,0.06 C,0.00 G,0.99 T,0.00 A,0.00 C,1.87 G,0.44 T,1.08 A,0.00 C,3.25\n'
            'G,0.09 T,0.97 A,0.50 C,1.00 G,1.72 T,0.07 A,0.00 C,0.92 \n')
        self.assertEqual(observed, expected)


class SffinfoFunctionTests(TestCase):
    def test_sffinfo_from_file(self):
        """sffinfo_from_file should return file object with sffinfo output."""
        test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        sff_fp = os.path.join(test_dir, 'data', 'test.sff')
        observed = sffinfo_from_file(sff_fp).read()
        self.assertEqual(observed, sffinfo_output)


sffinfo_output = '''Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  1504
  Index Length:  706
  # of Reads:    1
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FA6P1OK01CGMHQ
  Run Prefix:   R_2008_05_28_17_11_38_
  Region #:     1
  XY Location:  0892_1356

  Run Name:       R_2008_05_28_17_11_38_FLX02070135_adminrig_KnightLauber
  Analysis Name:  /data/2008_05_28/R_2008_05_28_17_11_38_FLX02070135_adminrig_KnightLauber/D_2008_05_28_21_13_06_FLX02070135_KnightLauber_FullAnalysisAmplicons
  Full Path:      /data/2008_05_28/R_2008_05_28_17_11_38_FLX02070135_adminrig_KnightLauber/D_2008_05_28_21_13_06_FLX02070135_KnightLauber_FullAnalysisAmplicons/../D_2008_05_29_13_52_01_FLX02070135_Knight_Lauber_jones_SignalProcessingAmplicons

  Read Header Len:  32
  Name Length:      14
  # of Bases:       77
  Clip Qual Left:   5
  Clip Qual Right:  52
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.00	1.01	0.00	0.00	0.96	0.00	1.02	0.00	1.02	0.00	0.00	0.99	0.00	1.00	0.00	1.00	0.00	0.00	1.00	0.00	1.10	0.00	1.08	0.00	0.00	1.46	0.00	0.88	0.18	0.00	2.69	1.01	0.08	0.96	0.00	0.02	0.92	0.08	0.00	0.98	0.68	0.00	0.89	0.00	0.00	1.15	0.00	1.13	0.00	0.02	1.12	0.05	0.15	1.84	0.00	1.10	0.00	2.47	0.96	0.86	1.06	0.00	1.96	0.12	0.93	0.13	1.65	1.06	0.06	0.00	0.99	0.00	0.00	1.87	0.44	1.08	0.00	3.25	0.09	0.97	0.50	1.00	1.72	0.07	0.00	0.92	0.58	0.00	0.00	0.59	0.06	0.11	0.09	0.07	0.06	0.16	0.00	0.24	0.03	0.00	0.12	0.06	0.16	0.00	0.18	0.00	0.00	0.14	0.00	0.15	0.00	0.18	0.00	0.03	0.14	0.03	0.13	0.01	0.19	0.00	0.02	0.33	0.05	0.00	0.16	0.10	0.35	0.01	0.21	0.04	0.09	0.18	0.13	0.19	0.00	0.10	0.51	0.26	0.00	0.23	0.19	0.27	0.01	0.29	0.05	0.14	0.17	0.16	0.18	0.27	0.09	0.26	0.10	0.18	0.23	0.15	0.22	0.13	0.37	0.11	0.11	0.26	0.59	0.14	0.06	0.33	0.34	0.26	0.05	0.27	0.44	0.19	0.10	0.35	0.27	0.15	0.34	0.28	0.45	0.14	0.16	0.34	0.27	0.12	0.07	0.25	0.18	0.12	0.04	0.23	0.16	0.12	0.05	0.20	0.16	0.11	0.03	0.21	0.16	0.10	0.02	0.21	0.16	0.12	0.02	0.20	0.15	0.10	0.02	0.23	0.15	0.11	0.02	0.22	0.14	0.09	0.02	0.20	0.13	0.09	0.01	0.19	0.13	0.08	0.02	0.17	0.12	0.08	0.03	0.17	0.09	0.08	0.01	0.14	0.09	0.07	0.01	0.15	0.09	0.06	0.01	0.13	0.08	0.06	0.00	0.13	0.08	0.05	0.02	0.12	0.07	0.05	0.01	0.11	0.07	0.05	0.00	0.10	0.07	0.05	0.01	0.11	0.08	0.04	0.00	0.10	0.06	0.05	0.01	0.09	0.06	0.04	0.01	0.08	0.07	0.05	0.00	0.08	0.06	0.05	0.00	0.09	0.06	0.04	0.00	0.09	0.06	0.04	0.01	0.08	0.06	0.04	0.00	0.09	0.06	0.03	0.00	0.09	0.06	0.02	0.00	0.09	0.06	0.04	0.00	0.08	0.05	0.03	0.00	0.07	0.05	0.02	0.00	0.08	0.04	0.03	0.00	0.07	0.04	0.03	0.00	0.07	0.05	0.02	0.00	0.07	0.05	0.02	0.00	0.06	0.04	0.02	0.00	0.06	0.03	0.03	0.00	0.08	0.02	0.00	0.00	0.07	0.03	0.01	0.00	0.06	0.03	0.02	0.00	0.05	0.03	0.02	0.00	0.05	0.03	0.01	0.00	0.06	0.02	0.00	0.00	0.05	0.01	0.01	0.00	0.04	0.01	0.01	0.00	0.04	0.01	0.01	0.00	0.05	0.01	0.00	0.00	0.04	0.02	0.01	0.00	0.03	0.02	0.01	0.00	0.03	0.01	0.00	0.00	0.03	0.00	0.00	0.00	0.03	0.00	0.00	0.00	0.02	0.00
Flow Indexes:	1	3	6	8	10	13	15	17	20	22	24	27	29	32	32	32	33	35	38	41	42	44	47	49	52	55	55	57	59	59	60	61	62	64	64	66	68	68	69	72	75	75	77	79	79	79	81	82	83	84	84	87	88	91	102	126	130	138	140	145	153	157	161	164	166	171	175	179	183	187	191	195	199	203	211	215	219
Bases:	tcagATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGCgcnnnannnnngnnnnnnnnnnnnn
Quality Scores:	32	32	32	32	32	32	32	32	32	32	32	25	25	21	21	21	28	32	32	31	30	30	32	32	32	33	31	25	18	18	20	18	32	30	28	23	22	22	24	28	18	19	18	16	16	16	17	18	13	17	27	21	20	21	0	0	0	17	0	0	0	0	0	17	0	0	0	0	0	0	0	0	0	0	0	0	0
'''

empty_sffinfo_output = '''Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  1504
  Index Length:  706
  # of Reads:    1
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG
'''

if __name__ == '__main__':
    main()
