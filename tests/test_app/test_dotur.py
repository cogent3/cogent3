#!/usr/bin/env python
# test_dotur.py

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.core.sequence import DnaSequence, RnaSequence
from cogent.core.alignment import DataError
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.dotur import Dotur, dotur_from_alignment, dotur_from_file,\
    remap_seq_names
from cogent.parse.dotur import OtuListParser

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"

def rna_distance(first,second):
    first = RnaSequence(first)
    return first.fracDiff(second)

class DoturTests(TestCase):

    def setUp(self):
        """Dotur general setUp method for all tests"""
        self.seqs1_unaligned = {'1':'ACUGCUAGCUAGUAGCGUACGUA',\
                                '2':'GCUACGUAGCUAC',\
                                '3':'GCGGCUAUUAGAUCGUA'}

        self.seqs2_aligned = {'a': 'UAGGCUCUGAUAUAAUAGCUCUC---------',\
                              'c': '------------UGACUACGCAU---------',\
                              'b': '----UAUCGCUUCGACGAUUCUCUGAUAGAGA'}
        
        self.seqs2_unaligned = {'a': 'UAGGCUCUGAUAUAAUAGCUCUC',\
                                'c': 'UGACUACGCAU',\
                                'b': 'UAUCGCUUCGACGAUUCUCUGAUAGAGA'}
        
        #self.seqs1 aligned to self.seqs2 with self.seqs2 included.
        self.seqs1_and_seqs2_aligned = \
            {'a': 'UAGGCUCUGAUAUAAUAGC-UCUC---------',\
             'b': '----UAUCGCUUCGACGAU-UCUCUGAUAGAGA',\
             'c': '------------UGACUAC-GCAU---------',\
             '1': '-ACUGCUAGCUAGUAGCGUACGUA---------',\
             '2': '----------GCUACGUAG-CUAC---------',\
             '3': '-----GCGGCUAUUAG-AU-CGUA---------',\
             }
        self.otu_list_string = \
"""unique	3	a	b	c
0.00	3	a	b	c
0.59	2	a,c	b
0.78	1	a,c,b
"""

        self.otu_res_list = [
            [0.0,3,[['a'],['b'],['c']]],\
            [0.0,3,[['a'],['b'],['c']]],\
            [float(0.59),2,[['a','c'],['b']]],\
            [float(0.78),1,[['a','c','b']]],\
            ]
        
        self.distance_matrix_string = \
"""
3
a       0.0  0.78125  0.59375
b       0.78125  0.0  0.71875
c       0.59375  0.71875  0.0
"""
        self.int_keys = {'seq_1': 'b', 'seq_0': 'a', 'seq_2': 'c'}
        self.otu_lists_unmapped = [\
            [['seq_0'], ['seq_1'], ['seq_2']],
            [['seq_0'], ['seq_1'], ['seq_2']],
            [['seq_0', 'seq_2'], ['seq_1']],
            [['seq_0', 'seq_2', 'seq_1']],
            ]
        self.otu_lists_mapped = [\
            [['a'], ['b'], ['c']],
            [['a'], ['b'], ['c']],
            [['a', 'c'], ['b']],
            [['a', 'c', 'b']],
            ]
        
        
        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test_for_dotur/'
        try:
            mkdir(self.temp_dir_spaces)
        except OSError:
            pass
        try:
            #create sequence files
            f = open(path.join(self.temp_dir, 'seqs1.sto'),'w')
            f.write('')
            f.close()
            self.d_mat_file = path.join(self.temp_dir, 'dmat.txt')
            d_mat = open(self.d_mat_file,'w')
            d_mat.write(self.distance_matrix_string)
            d_mat.close()

        except OSError:
            pass
    
    def test_base_command(self):
        """Dotur BaseCommand should return the correct BaseCommand"""
        c = Dotur()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','dotur']))
        c.Parameters['-l'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','dotur -l']))


    def test_changing_working_dir(self):
        """Dotur BaseCommand should change according to WorkingDir"""
        c = Dotur(WorkingDir='/tmp/dotur_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/dotur_test','/"; ','dotur']))
        c = Dotur()
        c.WorkingDir = '/tmp/dotur_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/dotur_test2','/"; ','dotur']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/dotur_test')
        rmdir('/tmp/dotur_test2')
    
    def test_remap_seq_names(self):
        """remap_seq_names should function as expected."""
        for unmapped, mapped in zip(self.otu_lists_unmapped,\
            self.otu_lists_mapped):
            self.assertEqual(remap_seq_names(unmapped,self.int_keys),mapped)
    
    def test_dotur_from_alignment(self):
        """dotur_from_alignment should behave correctly."""
        res = dotur_from_alignment(aln=self.seqs2_aligned,moltype=RNA,\
            distance_function=rna_distance)
        self.assertEqual(res,self.otu_res_list)

    def test_dotur_from_file(self):
        """dotur_from_alignment should behave correctly."""
        res = dotur_from_file(self.d_mat_file)
        self.assertEqual(res,self.otu_res_list)
    
    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)


if __name__ == '__main__':
    main()
