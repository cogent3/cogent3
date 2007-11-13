#!/usr/bin/env python

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.muscle import Muscle, muscle_seqs, aln_tree_seqs

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Catherine Lozupone", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

class GeneralSetUp(TestCase):

    def setUp(self):
        """Muscle general setUp method for all tests"""
        self.seqs1 = ['ACUGCUAGCUAGUAGCGUACGUA','GCUACGUAGCUAC',
            'GCGGCUAUUAGAUCGUA']
        self.labels1 = ['>1','>2','>3']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))

        self.seqs2=['UAGGCUCUGAUAUAAUAGCUCUC','UAUCGCUUCGACGAUUCUCUGAUAGAGA',
            'UGACUACGCAU']
        self.labels2=['>a','>b','>c']
        self.lines2 = flatten(zip(self.labels2,self.seqs2))
        
        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test for muscle/'
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
    

class MuscleTests(GeneralSetUp):
    """Tests for the Clustalw application controller"""

    def test_base_command(self):
        """Clustalw BaseCommand should return the correct BaseCommand"""
        c = Muscle()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','muscle']))
        c.Parameters['-in'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','muscle -in "seq.txt"']))
        c.Parameters['-cluster2'].on('neighborjoining')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','muscle -cluster2 neighborjoining' +
            ' -in "seq.txt"']))

    def test_changing_working_dir(self):
        """Clustalw BaseCommand should change according to WorkingDir"""
        c = Muscle(WorkingDir='/tmp/muscle_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/muscle_test','/"; ','muscle']))
        c = Muscle()
        c.WorkingDir = '/tmp/muscle_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/muscle_test2','/"; ','muscle']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/muscle_test')
        rmdir('/tmp/muscle_test2')
    
    def test_aln_tree_seqs(self):
        "aln_tree_seqs returns the muscle alignment and tree from iteration2"
        tree, aln = aln_tree_seqs(path.join(self.temp_dir, 'seq1.txt'), 
                                   tree_type="neighborjoining",
                                   WorkingDir=self.temp_dir,
                                   clean_up=True)
        self.assertEqual(str(tree), '((1:1.125,2:1.125):0.375,3:1.5);')
        self.assertEqual(len(aln), 6)
        self.assertEqual(aln[-2], '>3\n')
        self.assertEqual(aln[-1], 'GCGGCUAUUAGAUCGUA------\n')

    def test_aln_tree_seqs_spaces(self):
        "aln_tree_seqs should work on filename with spaces"
        try:
            #create sequence files
            f = open(path.join(self.temp_dir_spaces, 'muscle_test_seq1.txt'),'w')
            f.write('\n'.join(self.lines1))
            f.close()
        except OSError:
            pass
        tree, aln = aln_tree_seqs(path.join(self.temp_dir_spaces,\
                                    'muscle_test_seq1.txt'), 
                                    tree_type="neighborjoining",
                                    WorkingDir=getcwd(),
                                    clean_up=True)
        self.assertEqual(str(tree), '((1:1.125,2:1.125):0.375,3:1.5);')
        self.assertEqual(len(aln), 6)
        self.assertEqual(aln[-2], '>3\n')
        self.assertEqual(aln[-1], 'GCGGCUAUUAGAUCGUA------\n')
        remove(self.temp_dir_spaces+'/muscle_test_seq1.txt')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)
    

if __name__ == '__main__':
    main()
