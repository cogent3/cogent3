#!/usr/bin/env python

from os import getcwd, remove, rmdir
import tempfile, shutil
from cogent.util.unit_test import TestCase, main
from cogent.app.vienna_package import RNAfold, RNAsubopt

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

class RNAfoldTests(TestCase):
    """Tests for the RNAfold application controller"""
    
    def setUp(self):
        self.unnamed_seq = ['AUAGCUAGCUAUGCGCUAGC','ACGGCUAUAGCUAGCGA',\
            'gcuagcuauuauauaua']
        self.named_seq = ['>namedseq1','AUAGCUAGCUAUGCGCUAGC','>namedseq2',\
            'ACGGCUAUAGCUAGCGA']
        self.mixed_seq = ['>namedseq1','AUAGCUAGCUAUGCGCUAGC',
            'ACGGCUAUAGCUAGCGA','gcuagcuauuauauaua']
        self.mixed_seq2 = ['>namedseq2','AUAGCUAGCUAUGCGCUAGC',
            'ACGGCUAUAGCUAGCGA','gcuagcuauuauauaua']
        self.temp_dir = '/tmp/test_rnafold'

    def test_base_command(self):
        """RNAfold: BaseCommand should be ok for different parameter settings"""
        r = RNAfold()
        working_dir = getcwd()
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAfold -d1 -T 37 -S 1.07']))
        r.Parameters['-noLP'].on()
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAfold -d1 -noLP -T 37 -S 1.07']))
        r.Parameters['Temp'].on(15)
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAfold -d1 -noLP -T 15 -S 1.07']))
        r.Parameters['-d'].off()
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAfold -noLP -T 15 -S 1.07']))
        
    def test_changing_working_dir(self):
        """RNAfold: BaseCommand should be ok after changing the working dir"""
        #changing in initialization
        temp_dir = tempfile.mkdtemp()
        r = RNAfold(WorkingDir=temp_dir)
        self.assertEqual(r.BaseCommand,\
            'cd "%s/"; RNAfold -d1 -T 37 -S 1.07'%(temp_dir))
        #changing afterwards
        r = RNAfold()
        r.WorkingDir = temp_dir
        self.assertEqual(r.BaseCommand,\
            'cd "%s/"; RNAfold -d1 -T 37 -S 1.07'%(temp_dir))
        rmdir(temp_dir)
        
    def test_stdout(self):
        """RNAfold: StdOut should be as expected"""
        r = RNAfold()
        exp = '\n'.join(['>namedseq1','AUAGCUAGCUAUGCGCUAGC',\
            '...((((((.....)))))) ( -8.30)','ACGGCUAUAGCUAGCGA',\
            '...((((....)))).. ( -3.20)','GCUAGCUAUUAUAUAUA',\
            '................. (  0.00)'])+'\n'
        res = r(self.mixed_seq)
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_stdout_input_as_string(self):
        """RNAfold: StdOut with input_as_string"""
        r = RNAfold(InputHandler='_input_as_string')
        f = open('/tmp/rnatestfile','w')
        f.write('\n'.join(self.mixed_seq2))
        f.close()
        exp = '\n'.join(['>namedseq2','AUAGCUAGCUAUGCGCUAGC',\
            '...((((((.....)))))) ( -8.30)','ACGGCUAUAGCUAGCGA',\
            '...((((....)))).. ( -3.20)','GCUAGCUAUUAUAUAUA',\
            '................. (  0.00)'])+'\n'
        res = r('/tmp/rnatestfile')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/rnatestfile')

    def test_get_result_paths_unnamed_seq(self):
        """RNAfold: _get_result_paths() should work on unnamed seq"""
        r = RNAfold()
        res = r(self.unnamed_seq)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS','DP'])
        self.failUnless(res['DP'] is None)
        self.failUnless(res['SS'] is not None)
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()
    
    def test_get_result_paths_named_seq(self):
        """RNAfold: _get_result_paths() should work on named seq"""
       
        r = RNAfold()
        res = r(self.named_seq)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS','DP','namedseq1_ss',\
            'namedseq2_ss','namedseq1_dp','namedseq2_dp']) 
        res.cleanUp()

    def test_get_result_paths_mixed_seq(self):
        """RNAfold: _get_result_paths() should work on partly named seq"""

        r = RNAfold()
        res = r(self.mixed_seq)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS','DP','namedseq1_ss',\
            'namedseq1_dp']) 
        res.cleanUp()

    def test_get_result_paths_parameter(self):
        """RNAfold: _get_result_paths() should work with diff parameters"""

        r = RNAfold()
        r.Parameters['-p'].on()
        res = r(self.unnamed_seq)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS','DP'])
        self.failUnless(res['DP'] is not None)
        self.failUnless(res['SS'] is not None)
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()

    def test_get_result_paths_working_dir(self):
        """RNAfold: _get_result_paths() should work with diff working dir"""
        r = RNAfold(WorkingDir=self.temp_dir)
        res = r(self.unnamed_seq)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS','DP'])
        self.failUnless(res['DP'] is None)
        self.failUnless(res['SS'] is not None)
        self.failUnless(isinstance(res['SS'],file))
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()

    def test_zzz_general_cleanup(self):
        """Executed last, clean up temp_dir"""
        shutil.rmtree(self.temp_dir)

class RNAsuboptTests(TestCase):
    """Tests for the RNAsubopt application controller"""

    def test_base_command(self):
        """RNAsubopt: BaseCommand should be ok for different parameter settings
        """
        r = RNAsubopt()
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAsubopt -e 1 -d2 -T 37']))
        r.Parameters['-nsp'].on('GA')
        self.assertEqual(r.BaseCommand,\
            ''.join(\
                ['cd "',getcwd(),'/"; ','RNAsubopt -e 1 -d2 -nsp GA -T 37']))
        r.Parameters['Temp'].on(15)
        self.assertEqual(r.BaseCommand,\
            ''.join(\
                ['cd "',getcwd(),'/"; ','RNAsubopt -e 1 -d2 -nsp GA -T 15']))
        r.Parameters['-d'].off()
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAsubopt -e 1 -nsp GA -T 15']))
        
    def test_changing_working_dir(self):
        """RNAsubopt: BaseCommand should be ok after changing the working dir
        """
        temp_dir = tempfile.mkdtemp()
        #changing in initialization
        r = RNAsubopt(WorkingDir=temp_dir)
        self.assertEqual(r.BaseCommand,\
            'cd "%s/"; RNAsubopt -e 1 -d2 -T 37'%(temp_dir))
        #changing afterwards
        r = RNAsubopt()
        r.WorkingDir = temp_dir
        self.assertEqual(r.BaseCommand,\
            'cd "%s/"; RNAsubopt -e 1 -d2 -T 37'%(temp_dir))
        rmdir(temp_dir)
           
    def test_stdout(self):
        """RNAsubopt: StdOut should be as expected"""
        r = RNAsubopt()
        seq = ['AUAGCUAGCUAUGCGCUAGCGGAUUAGCUAGCUAGCGA',\
        'ucgaucgaucagcuagcuauuauauaua']
        
        exp = '\n'.join(
        ['AUAGCUAGCUAUGCGCUAGCGGAUUAGCUAGCUAGCGA  -1720    100',
        '.(((((((((.(((....)))....))))))))).... -16.20',
        '.(((((((((((.((....)).)).))))))))).... -17.20',
        '.((((((((((..((....))...)))))))))).... -16.60',
        '.((((((((((.((....))....)))))))))).... -16.40',
        '.(((((((((((((....)))...)))))))))).... -16.90',
        '.(((((((((((.((....)).).)))))))))).... -17.20',
        'UCGAUCGAUCAGCUAGCUAUUAUAUAUA      0    100',
        '......(((.((....))))).......   0.70',
        '..........((....))..........   0.60',
        '..((....))..................   0.90',
        '............................   0.00']) + '\n'
        res = r(seq)
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()

    def test_get_result_paths(self):
        """RNAsubopt: _get_result_paths() should create the right dict entries
        """
        r = RNAsubopt()
        seq = ['AUAGCUAGCUAUGCGCUAGCGGAUUAGCUAGCUAGCGA',\
        'ucgaucgaucagcuagcuauuauauaua']
        res = r(seq)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus'])
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()
        
        r = RNAsubopt({'-s':None,'-lodos':None,'-d':3,'-logML':None,\
            '-noLP':None,'-4':None,'-noGU':None,'-noCloseGU':None})
        res = r(seq)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus'])
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        #self.assertEqual(res['ExitStatus'],0) #platform-dependent?
        res.cleanUp()


if __name__ == '__main__':
    main()
