#!/usr/bin/env python
"""
Provides tests for lagan.py
"""
from os import path, getcwd, remove, rmdir, mkdir, listdir
import tempfile, shutil
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.lagan import Mlagan
from cogent.parse.fasta import MinimalFastaParser
from cogent import LoadSeqs

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Gavin Huttley", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

class GeneralSetUp(TestCase):
    
    def setUp(self):
        """Mlagan general setUp method for all tests"""
        self.seqs1 = ['GCTACGTTAC',
            'GCTACGTAGCTAC',
            'TAGCTAC']
        self.labels1 = ['>a','>b', '>c', '>d']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))
        

class MlaganTests(GeneralSetUp):
    """Tests for the Mlagan application controller"""
    
    def test_base_command(self):
        """Mlagan BaseCommand should return the correct BaseCommand"""
        c = Mlagan()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd ','"%s/"; ' % getcwd(),'mlagan']))
        
        c.Parameters["-tree"].on('"((human chimp), (mouse rat))"')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd ','"%s/"; ' % getcwd(),'mlagan -tree "((human chimp), (mouse rat))"']))
    
    def test_align(self):
        """test aligning samples"""
        # this takes a while, even for small examples like this!
        temp_dir = tempfile.mkdtemp()
        c = Mlagan(WorkingDir=temp_dir, InputHandler="_input_as_lines")
        res = c(self.lines1)
        align = res["StdOut"].readlines()
        self.assertEqual(align, ['>c\n', '------TAGCTAC\n',
                                 '>b\n', 'GCTACGTAGCTAC\n',
                                 '>a\n', 'GCTAC---GTTAC\n', '\n'])
        res.cleanUp()
        shutil.rmtree(temp_dir)
    
if __name__ == '__main__':
    main()
