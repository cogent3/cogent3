#!/usr/bin/env python

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.core.moltype import RNA
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.mafft import Mafft, align_unaligned_seqs, \
    add_seqs_to_alignment, align_two_alignments

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

class GeneralSetUp(TestCase):

    def setUp(self):
        """Mafft general setUp method for all tests"""
        self.seqs1 = ['ACUGCUAGCUAGUAGCGUACGUA','GCUACGUAGCUAC',
            'GCGGCUAUUAGAUCGUA']
        
        self.labels1 = ['>1','>2','>3']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))
        
        self.aligned1 = {'1': 'acugcuagcuaguagcguacgua',\
                         '2': 'gcuacguagcuac----------',\
                         '3': 'gcggcuauuagau------cgua',\
                         }

        
        self.seqs2=['UAGGCUCUGAUAUAAUAGCUCUC','UAUCGCUUCGACGAUUCUCUGAUAGAGA',
            'UGACUACGCAU']
        self.labels2=['>a','>b','>c']
        self.lines2 = flatten(zip(self.labels2,self.seqs2))
        
        self.aligned2 = {'a': 'UAGGCUCUGAUAUAAUAGCUCUC---------',\
                         'b': 'UA----UCGCUUCGACGAUUCUCUGAUAGAGA',\
                         'c': 'UG------------ACUACGCAU---------',\
                         }

        
        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test for mafft/'
        try:
            mkdir(self.temp_dir_spaces)
        except OSError:
            pass
        try:
            #create sequence files
            f = open(path.join(self.temp_dir, 'seq1.txt'),'w')
            f.write('\n'.join(self.lines1))
            f.close()
            g = open(path.join(self.temp_dir, 'seq2.txt'),'w')
            g.write('\n'.join(self.lines2))
            g.close()
        except OSError:
            pass
    

class MafftTests(GeneralSetUp):
    """Tests for the Mafft application controller"""

    def test_base_command(self):
        """Mafft BaseCommand should return the correct BaseCommand"""
        c = Mafft()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','mafft']))
        c.Parameters['--quiet'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','mafft --quiet']))
        c.Parameters['--globalpair'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','mafft --globalpair --quiet']))
        c.Parameters['--maxiterate'].on(1000)
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ',"""mafft --maxiterate 1000 --globalpair --quiet"""]))

    def test_changing_working_dir(self):
        """Mafft BaseCommand should change according to WorkingDir"""
        c = Mafft(WorkingDir='/tmp/mafft_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/mafft_test','/"; ','mafft']))
        c = Mafft()
        c.WorkingDir = '/tmp/mafft_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/mafft_test2','/"; ','mafft']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/mafft_test')
        rmdir('/tmp/mafft_test2')
    
    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)
    
    def test_align_unaligned_seqs(self):
        """align_unaligned_seqs should work as expected"""
        res = align_unaligned_seqs(self.seqs1, RNA)
        self.assertEqual(res.toFasta(), align1)
        res = align_unaligned_seqs(self.lines2, RNA)
        self.assertEqual(res.toFasta(), align2)
    
    def test_add_seqs_to_alignment(self):
        """add_seqs_to_alignment should work as expected."""
        res = add_seqs_to_alignment(self.lines1,self.aligned2, RNA)
        self.assertEqual(res.toFasta(), add_seqs_align)

    def test_align_two_alignments(self):
        """align_two_alignments should work as expected."""
        res = align_two_alignments(self.aligned1, self.aligned2, RNA)
        self.assertEqual(res.toFasta(), align_two_align)
    
align1 = ">seq_0\nACUGCUAGCUAGUAGCGUACGUA\n>seq_1\nGCUACGUAGCUAC----------\n>seq_2\nGCGGCUAUUAGAU------CGUA"

align2 = ">a\nUAGGCUCUGAUAUAAUAGCUCUC---------\n>b\nUA----UCGCUUCGACGAUUCUCUGAUAGAGA\n>c\nUG------------ACUACGCAU---------"

add_seqs_align = """>1\nACUGC-UAGCUAGUAGCGUACGUA--------\n>2\nGCUACGUAGCUA-----------C--------\n>3\nGCGGCUAUUAGAUCGUA---------------\n>a\nUAGGCUCUGAUAUAAUAGCUCUC---------\n>b\nUA----UCGCUUCGACGAUUCUCUGAUAGAGA\n>c\nUG------------ACUACGCAU---------"""

align_two_align = """>1\nACUGCUAGCUAGUAGCGUACGUA---------\n>2\nGCUACGUAGCUAC-------------------\n>3\nGCGGCUAUUAGAU------CGUA---------\n>a\nUAGGCUCUGAUAUAAUAGCUCUC---------\n>b\nUA----UCGCUUCGACGAUUCUCUGAUAGAGA\n>c\nUG------------ACUACGCAU---------"""

if __name__ == '__main__':
    main()
