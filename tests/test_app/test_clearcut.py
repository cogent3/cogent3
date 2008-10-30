#!/usr/bin/env python

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.core.alignment import DataError
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.clearcut import Clearcut, build_tree_from_alignment

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

class GeneralSetUp(TestCase):

    def setUp(self):
        """Clearcut general setUp method for all tests"""
        self.seqs1 = ['ACUGCUAGCUAGUAGCGUACGUA','GCUACGUAGCUAC',
            'GCGGCUAUUAGAUCGUA']
        
        self.labels1 = ['>1','>2','>3']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))

        self.seqs2=['UAGGCUCUGAUAUAAUAGCUCUC','UAUCGCUUCGACGAUUCUCUGAUAGAGA',
            'UGACUACGCAU']
        self.labels2=['>a','>b','>c']
        self.lines2 = flatten(zip(self.labels2,self.seqs2))
        
        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test for clearcut/'
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
    

class ClearcutTests(GeneralSetUp):
    """Tests for the Clearcut application controller"""

    def test_base_command(self):
        """Clearcut BaseCommand should return the correct BaseCommand"""
        c = Clearcut()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','clearcut -d -q']))
        c.Parameters['--in'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','clearcut -d --in="seq.txt" -q']))


    def test_changing_working_dir(self):
        """Clearcut BaseCommand should change according to WorkingDir"""
        c = Clearcut(WorkingDir='/tmp/clearcut_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/clearcut_test','/"; ','clearcut -d -q']))
        c = Clearcut()
        c.WorkingDir = '/tmp/clearcut_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/clearcut_test2','/"; ','clearcut -d -q']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/clearcut_test')
        rmdir('/tmp/clearcut_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)

    def test_build_tree_from_alignment(self):
        """Clearcut should return a tree built from the passed alignment"""
        tree_short = build_tree_from_alignment(build_tree_seqs_short,\
            moltype=DNA)
        num_seqs = flatten(build_tree_seqs_short).count('>')
        self.assertEqual(len(tree_short.tips()), num_seqs)

        tree_long = build_tree_from_alignment(build_tree_seqs_long, moltype=DNA)
        seq_names = []
        for line in build_tree_seqs_long.split('\n'):
            if line.startswith('>'):
                seq_names.append(line[1:])

        for node in tree_long.tips():
            if node.Name not in seq_names:
                self.fail()
        #repeat with best_tree = True
        tree_long = build_tree_from_alignment(build_tree_seqs_long,\
            best_tree=True,\
            moltype=DNA)
        seq_names = []
        for line in build_tree_seqs_long.split('\n'):
            if line.startswith('>'):
                seq_names.append(line[1:])

        for node in tree_long.tips():
            if node.Name not in seq_names:
                self.fail()
        
        #build_tree_from_alignment should raise DataError when constructing
        # an Alignment from unaligned sequences. Clearcut only allows aligned
        # or a distance matrix as input.
        self.assertRaises(DataError,build_tree_from_alignment,\
            build_tree_seqs_unaligned,DNA)


align1 = ">seq_0\nACUGCUAGCUAGUAGCGUACGUA\n>seq_1\n---GCUACGUAGCUAC-------\n>seq_2\nGCGGCUAUUAGAUCGUA------"

build_tree_seqs_short = """>clearcut_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AGCTTTAAATCATGCCAGTG
>clearcut_test_seqs_1
GACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
TGCTTTCAATAATGCCAGTG
>clearcut_test_seqs_2
AACCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
TGCTTTGAATCATGCCAGTA
>clearcut_test_seqs_3
AAACCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
TGCTTTACATCATGCAAGTG
>clearcut_test_seqs_4
AACCGCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
TGCTTTAAATCATGCCAGTG
>clearcut_test_seqs_5
AACCCCCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
TGCTTTAAATCATGCCAGTT
>clearcut_test_seqs_6
GACCCCCGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
TACTTTAGATCATGCCGGTG
>clearcut_test_seqs_7
AACCCCCACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
TGCTTTAAATCATGCCAGTG
>clearcut_test_seqs_8
AACCCCCACGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
TGCATTAAATCATGCCAGTG
>clearcut_test_seqs_9
AAGCCCCACGGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
TGCTTTAAATCCTGACAGCG
"""

build_tree_seqs_long = """>clearcut_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AGCTTTAAATCATGCCAGTG
>clearcut_test_seqsaaaaaaaa_1
GACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
TGCTTTCAATAATGCCAGTG
>clearcut_test_seqsaaaaaaaa_2
AACCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
TGCTTTGAATCATGCCAGTA
>clearcut_test_seqsaaaaaaaa_3
AAACCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
TGCTTTACATCATGCAAGTG
>clearcut_test_seqsaaaaaaaa_4
AACCGCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
TGCTTTAAATCATGCCAGTG
>clearcut_test_seqsaaaaaaaa_5
AACCCCCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
TGCTTTAAATCATGCCAGTT
>clearcut_test_seqsaaaaaaaa_6
GACCCCCGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
TACTTTAGATCATGCCGGTG
>clearcut_test_seqsaaaaaaaa_7
AACCCCCACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
TGCTTTAAATCATGCCAGTG
>clearcut_test_seqsaaaaaaaa_8
AACCCCCACGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
TGCATTAAATCATGCCAGTG
>clearcut_test_seqsaaaaaaaa_9
AAGCCCCACGGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
TGCTTTAAATCCTGACAGCG
"""

#Unaligned seqs.  First two sequences are 3 nucleotides shorter.
build_tree_seqs_unaligned = """>clearcut_test_seqs_0
CCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AGCTTTAAATCATGCCAGTG
>clearcut_test_seqs_1
CCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
TGCTTTCAATAATGCCAGTG
>clearcut_test_seqs_2
AACCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
TGCTTTGAATCATGCCAGTA
>clearcut_test_seqs_3
AAACCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
TGCTTTACATCATGCAAGTG
>clearcut_test_seqs_4
AACCGCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
TGCTTTAAATCATGCCAGTG
>clearcut_test_seqs_5
AACCCCCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
TGCTTTAAATCATGCCAGTT
>clearcut_test_seqs_6
GACCCCCGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
TACTTTAGATCATGCCGGTG
>clearcut_test_seqs_7
AACCCCCACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
TGCTTTAAATCATGCCAGTG
>clearcut_test_seqs_8
AACCCCCACGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
TGCATTAAATCATGCCAGTG
>clearcut_test_seqs_9
AAGCCCCACGGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
TGCTTTAAATCCTGACAGCG
"""

if __name__ == '__main__':
    main()
