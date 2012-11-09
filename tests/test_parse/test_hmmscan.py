#!/usr/bin/env python

"""For use with hmmer3"""

from cogent.parse.hmmscan import MinimalTbloutParser, RichTbloutParser
from cogent.util.unit_test import TestCase, main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald"] 
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

class hmmscanTests(TestCase):
    def setUp(self):
        pass
    def test_MinimalTbloutParser(self):
        """MinimalTbloutParser tests"""
        exp1 = ['Sulfotransfer_3','PF13469.1','NC_008340.4664','-','1.7e-14',
                '43.6','0.0','3.4e-14','42.6','0.0','1.5','1','0','0','1','1',
                '1','1','Sulfotransferase family']
        exp2 = ['UreD', 'PF01774.12', 'NC_008340.6164', '-', '4e-62', '198.0', 
                '0.0', '4.9e-62', '197.7', '0.0', '1.1', '1', '0', '0', '1', 
                '1', '1', '1', 'UreD urease accessory protein']
        exp3 = ['UreF', 'PF01730.11', 'NC_008340.6172', '-', '1.9e-51', 
                '162.7', '10.3', '2.8e-51', '162.1', '7.1', '1.3', '1', '0', 
                '0', '1', '1', '1', '1', 'UreF']

        parser = MinimalTbloutParser(tblout)
        self.assertEqual(parser.next(), exp1)
        self.assertEqual(parser.next(), exp2)
        self.assertEqual(parser.next(), exp3)
        self.assertRaises(StopIteration, parser.next)

    def test_RichTbloutParser(self):
        """More useful returns"""
        headers = ["target name", "target accession", "query name",
                   "query accession","full seq E-value", "full seq score", 
                   "full seq bias","best 1 domain E-value", 
                   "best 1 domain score", "best 1 domain bias", "exp","reg",
                   "clu","ov","env","dom","rep","inc","description of target"]
        types = [str,str,str,str,float,float,float,float,float,float,float,int,
                 int,int,int,int,int,int,str]
        expbase1 = ['Sulfotransfer_3','PF13469.1','NC_008340.4664','-','1.7e-14',
                '43.6','0.0','3.4e-14','42.6','0.0','1.5','1','0','0','1','1',
                '1','1','Sulfotransferase family']
        expbase2 = ['UreD', 'PF01774.12', 'NC_008340.6164', '-', '4e-62', 
                    '198.0', '0.0', '4.9e-62', '197.7', '0.0', '1.1', '1', '0', 
                    '0', '1', '1', '1', '1', 'UreD urease accessory protein']
        expbase3 = ['UreF', 'PF01730.11', 'NC_008340.6172', '-', '1.9e-51', 
                '162.7', '10.3', '2.8e-51', '162.1', '7.1', '1.3', '1', '0', 
                '0', '1', '1', '1', '1', 'UreF']
       
        exp1 = dict([(h,t(e)) for h,t,e in zip(headers,types,expbase1)])
        exp2 = dict([(h,t(e)) for h,t,e in zip(headers,types,expbase2)])
        exp3 = dict([(h,t(e)) for h,t,e in zip(headers,types,expbase3)])

        parser = RichTbloutParser(tblout)
        self.assertEqual(parser.next(), exp1)
        self.assertEqual(parser.next(), exp2)
        self.assertEqual(parser.next(), exp3)
        self.assertRaises(StopIteration, parser.next)

tblout = """#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
Sulfotransfer_3      PF13469.1  NC_008340.4664       -            1.7e-14   43.6   0.0   3.4e-14   42.6   0.0   1.5   1   0   0   1   1   1   1 Sulfotransferase family
UreD                 PF01774.12 NC_008340.6164       -              4e-62  198.0   0.0   4.9e-62  197.7   0.0   1.1   1   0   0   1   1   1   1 UreD urease accessory protein
UreF                 PF01730.11 NC_008340.6172       -            1.9e-51  162.7  10.3   2.8e-51  162.1   7.1   1.3   1   0   0   1   1   1   1 UreF""".splitlines()

if __name__ == '__main__':
    main()
