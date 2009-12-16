#!/usr/bin/env python

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.dialign import Dialign

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class GeneralSetUp(TestCase):
    
    def setUp(self):
        """Dialign general setUp method for all tests"""
        self.seqs1 = ['LDTAPCLFSDGSPQKAAYVLWDQTILQQDITPLPSHETHSAQKGELLALICGLRAAK',
            'PDADHTWYTDGSSLLQEGQRKAGAAVTTETEVIWAKALDAGTSAQRAELIALTQALKM',
            'RPGLCQVFADATPTGWGLVMGHQRMRGTFSAPLPIHTAELLAACFARSRSGANIIGTDNSVV',
            'MLKQVEIFTDGSCLGNPGPGGYGAILRYRGREKTFSAGYTRTTNNRMELMAAIV']
        self.labels1 = ['>HTL2','>MMLV', '>HEPB', '>ECOL']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))
        self.out = \
"""
                            DIALIGN 2.2.1 
                            *************

           Program code written by Burkhard Morgenstern and Said Abdeddaim 
              e-mail contact: dialign (at) gobics (dot) de 

           Published research assisted by DIALIGN 2 should cite:  

              Burkhard Morgenstern (1999).
              DIALIGN 2: improvement of the segment-to-segment
              approach to multiple sequence alignment.
              Bioinformatics 15, 211 - 218. 

           For more information, please visit the DIALIGN home page at 

              http://bibiserv.techfak.uni-bielefeld.de/dialign/ 

          ************************************************************



    program call:  dialign2-2 -fa -fn /tmp/di/seq1.fasta /tmp/di/seq1.txt  


    Aligned sequences:          length:
    ==================          =======

    1) HTL2                        57
    2) MMLV                        58
    3) HEPB                        62
    4) ECOL                        54

    Average seq. length:           57.8 


    Please note that only upper-case letters are considered to be aligned. 


    Alignment (DIALIGN format):
    ===========================

 HTL2               1   ldtapC-LFS DGS------P QKAAYVL--- ----WDQTIL QQDITPLPSH 
 MMLV               1   pdadhtw-YT DGSSLLQEGQ RKAGAAVtte teviWa---- KALDAG---T 
 HEPB               1   rpgl-CQVFA DAT------P TGWGLVM--- ----GHQRMR GTFSAPLPIH 
 ECOL               1   mlkqv-EIFT DGSCLGNPGP GGYGAIL--- ----RYRGRE KTFSAGytrT 

                        0000000588 8882222229 9999999000 0000666666 6666633334 

 HTL2              37   ethSAQKGEL LALICGLRAa k--------- --- 
 MMLV              43   ---SAQRAEL IALTQALKm- ---------- --- 
 HEPB              37   t------AEL LAA-CFARSr sganiigtdn svv 
 ECOL              43   ---TNNRMEL MAAIv----- ---------- --- 

                        0003333455 5533333300 0000000000 000 




    Sequence tree:
    ==============

 Tree constructed using UPGMAbased on DIALIGN fragment weight scores

 ((HTL2        :0.130254MMLV        :0.130254):0.067788(HEPB        :0.120520ECOL        :0.120520):0.077521);



"""
        self.temp_dir = tempfile.mkdtemp()
        try:
            #create sequence files
            f = open(path.join(self.temp_dir, 'seq1.txt'),'w')
            f.write('\n'.join(self.lines1))
            f.close()
        except OSError:
            pass

class DialignTests(GeneralSetUp):
    """Tests for the Dialign application controller"""
    
    def test_base_command(self):
        """Dialign BaseCommand should return the correct BaseCommand"""
        c = Dialign()
        self.assertEqual(c.BaseCommand,
            ''.join(['cd ','"%s/"; ' % getcwd(),'dialign2-2']))
        
        c.Parameters["-fa"].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd ','"%s/"; ' % getcwd(),'dialign2-2 -fa']))
    
    def test_align(self):
        """test aligning samples"""
        c = Dialign(WorkingDir=self.temp_dir,
            params={"-fn":path.join(self.temp_dir,"seq1.txt")})
        c.Parameters["-fa"].on()
        res = c(path.join(self.temp_dir, 'seq1.txt'))
        align = "".join(res["Align"].readlines())
    
    def test_cleaning_up(self):
        """not a test, just removes the temp files"""
        shutil.rmtree(self.temp_dir)

if __name__ == '__main__':
    main()
