#!/usr/bin/env python
# test_infernal.py

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.core.moltype import DNA, RNA, PROTEIN
from cogent.core.alignment import DataError
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.infernal import Cmalign, Cmbuild, Cmcalibrate, Cmemit, Cmscore,\
    Cmsearch, Cmstat, cmbuild_from_alignment, cmbuild_from_file, \
    cmalign_from_alignment, cmalign_from_file, cmsearch_from_alignment,\
    cmsearch_from_file
from cogent.parse.rfam import MinimalRfamParser, ChangedRnaSequence, \
    ChangedSequence
from cogent.format.stockholm import stockholm_from_alignment
from cogent.struct.rna2d import ViennaStructure, wuss_to_vienna

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
        """Infernal general setUp method for all tests"""
        self.seqs1_unaligned = {'1':'ACUGCUAGCUAGUAGCGUACGUA',\
                                '2':'GCUACGUAGCUAC',\
                                '3':'GCGGCUAUUAGAUCGUA'}
        self.struct1_unaligned_string = '....(((...)))....'
        self.seqs1_unaligned_gaps = {'1':'ACUGCUAGCUAGU-AGCGUAC--GUA',\
                                     '2':'--GCUACGUAGCUAC',\
                                     '3':'GCGGCUAUUAGAUCGUA--'}
        
        
        
        self.seqs2_aligned = {'a': 'UAGGCUCUGAUAUAAUAGCUCUC---------',\
                              'c': '------------UGACUACGCAU---------',\
                              'b': '----UAUCGCUUCGACGAUUCUCUGAUAGAGA'}
        
        self.seqs2_unaligned = {'a': 'UAGGCUCUGAUAUAAUAGCUCUC',\
                                'c': 'UGACUACGCAU',\
                                'b': 'UAUCGCUUCGACGAUUCUCUGAUAGAGA'}
        
        self.struct2_aligned_string = '............((.(...)))..........'
        self.struct2_aligned_dict = {'SS_cons':self.struct2_aligned_string}
        
        self.lines2 = stockholm_from_alignment(aln=self.seqs2_aligned,\
            GC_annotation=self.struct2_aligned_dict)
        
        #self.seqs1 aligned to self.seqs2 with self.seqs2 included.
        self.seqs1_and_seqs2_aligned = \
            {'a': 'UAGGCUCUGAUAUAAUAGC-UCUC---------',\
             'b': '----UAUCGCUUCGACGAU-UCUCUGAUAGAGA',\
             'c': '------------UGACUAC-GCAU---------',\
             '1': '-ACUGCUAGCUAGUAGCGUACGUA---------',\
             '2': '----------GCUACGUAG-CUAC---------',\
             '3': '-----GCGGCUAUUAG-AU-CGUA---------',\
             }
             
        self.seqs1_and_seqs2_aligned_struct_string = \
            '............((.(....)))..........'
        
        #self.seqs1 aligned to self.seqs2 without self.seqs2 included.
        self.seqs1_aligned = \
            {'1': 'ACUGCUAGCUAGUAGCGUACGUA',\
             '2': '---------GCUACGUAG-CUAC',\
             '3': '----GCGGCUAUUAG-AU-CGUA',\
             }
             
        self.seqs1_aligned_struct_string = \
            '...........((.(....))).'
        
        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test for infernal/'
        try:
            mkdir(self.temp_dir_spaces)
        except OSError:
            pass
        try:
            #create sequence files
            f = open(path.join(self.temp_dir, 'seqs1.sto'),'w')
            f.write(self.lines2)
            f.close()
            #create cm file.
            self.cmfile = path.join(self.temp_dir, 'aln2.cm')
            cm = open(self.cmfile,'w')
            cm.write(ALN1_CM)
            cm.close()
            #create alignment file used to create cm file.
            self.aln2_file = path.join(self.temp_dir, 'aln2.sto')
            af = open(self.aln2_file,'w')
            af.write(self.lines2)
            af.close()
        except OSError:
            pass
    

class CmalignTests(GeneralSetUp):
    """Tests for the Cmalign application controller"""

    def test_base_command(self):
        """Infernal BaseCommand should return the correct BaseCommand"""
        c = Cmalign()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmalign']))
        c.Parameters['-l'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmalign -l']))


    def test_changing_working_dir(self):
        """Infernal BaseCommand should change according to WorkingDir"""
        c = Cmalign(WorkingDir='/tmp/cmalign_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmalign_test','/"; ','cmalign']))
        c = Cmalign()
        c.WorkingDir = '/tmp/cmalign_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmalign_test2','/"; ','cmalign']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cmalign_test')
        rmdir('/tmp/cmalign_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)
    
    def test_cmalign_from_alignment(self):
        """cmalign_from_alignment should work as expected.
        """
        #Align with cmalign_from_alignment without original alignment.
        aln, struct = cmalign_from_alignment(aln=self.seqs2_aligned,\
            structure_string=self.struct2_aligned_string,\
            seqs=self.seqs1_unaligned_gaps,moltype=RNA,include_aln=False)
        #Check correct alignment
        self.assertEqual(aln.todict(),self.seqs1_aligned)
        #Check correct struct
        self.assertEqual(wuss_to_vienna(str(struct)),\
            self.seqs1_aligned_struct_string)

        #should work with gapped seqs.  Need to test this is taken care of
        # since cmalign segfaults when there are gaps in the seqs to be aligned.
        aln, struct = cmalign_from_alignment(aln=self.seqs2_aligned,\
            structure_string=self.struct2_aligned_string,\
            seqs=self.seqs1_unaligned_gaps,moltype=RNA)
        #alignment should be correct
        self.assertEqual(aln.todict(),self.seqs1_and_seqs2_aligned)
        #structure should be correct
        self.assertEqual(wuss_to_vienna(str(struct)),\
            self.seqs1_and_seqs2_aligned_struct_string)
        
        #should work with ungapped seqs.
        aln, struct = cmalign_from_alignment(aln=self.seqs2_aligned,\
            structure_string=self.struct2_aligned_string,\
            seqs=self.seqs1_unaligned_gaps,moltype=RNA)
        #alignment should be correct
        self.assertEqual(aln.todict(),self.seqs1_and_seqs2_aligned)
        #structure should be correct
        self.assertEqual(wuss_to_vienna(str(struct)),\
            self.seqs1_and_seqs2_aligned_struct_string)
        
        #should return standard out
        aln, struct,stdout = cmalign_from_alignment(aln=self.seqs2_aligned,\
            structure_string=self.struct2_aligned_string,\
            seqs=self.seqs1_unaligned_gaps,moltype=RNA,\
            return_stdout=True)
        #Test that standard out is same length as expected
        self.assertEqual(len(stdout.split('\n')),\
            len(CMALIGN_STDOUT.split('\n')))
    
    def test_cmalign_from_file(self):
        """cmalign_from_file should work as expected.
        """
        #Align with cmalign_from_file without original alignment.
        aln,struct = cmalign_from_file(cm_file_path=self.cmfile,\
            seqs=self.seqs1_unaligned,\
            moltype=RNA)
        #Check correct alignment
        self.assertEqual(aln.todict(),self.seqs1_aligned)
        #Check correct struct
        self.assertEqual(wuss_to_vienna(str(struct)),\
            self.seqs1_aligned_struct_string)
        
        #Align with cmalign_from_file using original alignment.
        aln,struct = cmalign_from_file(cm_file_path=self.cmfile,\
            seqs=self.seqs1_unaligned,\
            moltype=RNA,\
            alignment_file_path=self.aln2_file,\
            include_aln=True)
        #alignment should be correct
        self.assertEqual(aln.todict(),self.seqs1_and_seqs2_aligned)
        #structure should be correct
        self.assertEqual(wuss_to_vienna(str(struct)),\
            self.seqs1_and_seqs2_aligned_struct_string)
    

class CmbuildTests(GeneralSetUp):
    """Tests for the Cmbuild application controller"""

    def test_base_command(self):
        """Infernal BaseCommand should return the correct BaseCommand"""
        c = Cmbuild()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmbuild']))
        c.Parameters['-A'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmbuild -A']))

    def test_changing_working_dir(self):
        """Infernal BaseCommand should change according to WorkingDir"""
        c = Cmbuild(WorkingDir='/tmp/cmbuild_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmbuild_test','/"; ','cmbuild']))
        c = Cmbuild()
        c.WorkingDir = '/tmp/cmbuild_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmbuild_test2','/"; ','cmbuild']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cmbuild_test')
        rmdir('/tmp/cmbuild_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)
    
    def test_cmbuild_from_alignment(self):
        """cmbuild_from_alignment should work as expected.
        """
        #Test unaligned seqs and unaligned struct fail.
        #DataError should be raised with Alignment is constructed
        self.assertRaises(DataError,cmbuild_from_alignment,\
            self.seqs1_unaligned,self.struct1_unaligned_string)
            
        #Test aligned seqs and unaligned struct fail.
        self.assertRaises(ValueError,cmbuild_from_alignment,\
            self.seqs2_aligned,self.struct1_unaligned_string)
        
        #Test get cm back without alignment.
        cm_res = cmbuild_from_alignment(self.seqs2_aligned,\
            self.struct2_aligned_string)
        cm_lines = cm_res.split('\n')
        ALN1_CM_lines = ALN1_CM.split('\n')
        #Check that the same number of lines are in both CMs
        self.assertEqual(len(cm_lines),len(ALN1_CM_lines))
        
        #The first 13 lines are unique to the specific run.  The res of the
        # CM should be the same, since built from the same data.
        self.assertEqual(cm_lines[13:],ALN1_CM_lines[13:])
        
        #Make sure same alignment is returned if return_alignment=True
        cm_res, cm_aln = cmbuild_from_alignment(self.seqs2_aligned,\
            self.struct2_aligned_string,return_alignment=True)
        self.assertEqual(cm_aln,self.lines2)
    
    def test_cmbuild_from_file(self):
        """cmbuild_from_file should work as expected.
        """
        cm_res = cmbuild_from_file(self.temp_dir+'/seqs1.sto')
        cm_lines = cm_res.split('\n')
        ALN1_CM_lines = ALN1_CM.split('\n')
        #Check that the same number of lines are in both CMs
        self.assertEqual(len(cm_lines),len(ALN1_CM_lines))
        
        #The first 13 lines are unique to the specific run.  The res of the
        # CM should be the same, since built from the same data.
        self.assertEqual(cm_lines[13:],ALN1_CM_lines[13:])
        
        #Make sure same alignment is returned if return_alignment=True
        cm_res, cm_aln = cmbuild_from_alignment(self.seqs2_aligned,\
            self.struct2_aligned_string,return_alignment=True)
        self.assertEqual(cm_aln,self.lines2)

class CmcalibrateTests(GeneralSetUp):
    """Tests for the Cmcalibrate application controller"""

    def test_base_command(self):
        """Infernal BaseCommand should return the correct BaseCommand"""
        c = Cmcalibrate()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmcalibrate']))
        c.Parameters['--mpi'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmcalibrate --mpi']))


    def test_changing_working_dir(self):
        """Infernal BaseCommand should change according to WorkingDir"""
        c = Cmcalibrate(WorkingDir='/tmp/cmcalibrate_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmcalibrate_test','/"; ','cmcalibrate']))
        c = Cmcalibrate()
        c.WorkingDir = '/tmp/cmcalibrate_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmcalibrate_test2','/"; ','cmcalibrate']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cmcalibrate_test')
        rmdir('/tmp/cmcalibrate_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)

class CmemitTests(GeneralSetUp):
    """Tests for the Cmemit application controller"""

    def test_base_command(self):
        """Infernal BaseCommand should return the correct BaseCommand"""
        c = Cmemit()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmemit']))
        c.Parameters['-u'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmemit -u']))


    def test_changing_working_dir(self):
        """Infernal BaseCommand should change according to WorkingDir"""
        c = Cmemit(WorkingDir='/tmp/cmemit_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmemit_test','/"; ','cmemit']))
        c = Cmemit()
        c.WorkingDir = '/tmp/cmemit_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmemit_test2','/"; ','cmemit']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cmemit_test')
        rmdir('/tmp/cmemit_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)

class CmscoreTests(GeneralSetUp):
    """Tests for the Cmscore application controller"""

    def test_base_command(self):
        """Infernal BaseCommand should return the correct BaseCommand"""
        c = Cmscore()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmscore']))
        c.Parameters['-l'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmscore -l']))


    def test_changing_working_dir(self):
        """Infernal BaseCommand should change according to WorkingDir"""
        c = Cmscore(WorkingDir='/tmp/cmscore_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmscore_test','/"; ','cmscore']))
        c = Cmscore()
        c.WorkingDir = '/tmp/cmscore_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmscore_test2','/"; ','cmscore']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cmscore_test')
        rmdir('/tmp/cmscore_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)


class CmsearchTests(GeneralSetUp):
    """Tests for the Cmsearch application controller"""

    def test_base_command(self):
        """Infernal BaseCommand should return the correct BaseCommand"""
        c = Cmsearch()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmsearch']))
        c.Parameters['-p'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmsearch -p']))


    def test_changing_working_dir(self):
        """Infernal BaseCommand should change according to WorkingDir"""
        c = Cmsearch(WorkingDir='/tmp/cmsearch_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmsearch_test','/"; ','cmsearch']))
        c = Cmsearch()
        c.WorkingDir = '/tmp/cmsearch_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmsearch_test2','/"; ','cmsearch']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cmsearch_test')
        rmdir('/tmp/cmsearch_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)
    
    def test_cmsearch_from_alignment_no_hits(self):
        """cmsearch_from_alignment should work as expected
        """
        search_res = cmsearch_from_alignment(aln=self.seqs2_aligned,\
            structure_string=self.struct2_aligned_string,\
            seqs=self.seqs1_unaligned,moltype=RNA)
        self.assertEqual(search_res,[])
            
    def test_cmsearch_from_alignment(self):
        """cmsearch_from_alignment should work as expected
        """
        exp_search_res = [['a', 5, 23, 1, 19, 12.85, '-', 37],\
                          ['b', 1, 19, 1, 19, 14.359999999999999, '-', 47]]
        search_res = cmsearch_from_alignment(aln=self.seqs2_aligned,\
            structure_string=self.struct2_aligned_string,\
            seqs=self.seqs2_unaligned,moltype=RNA)
        
        self.assertEqual(search_res,exp_search_res)
    
    def test_cmsearch_from_file_no_hits(self):
        """cmsearch_from_file should work as expected
        """
        search_res = cmsearch_from_file(cm_file_path=self.cmfile,\
            seqs=self.seqs1_unaligned,moltype=RNA)
        self.assertEqual(search_res,[])
    
    def test_cmsearch_from_file(self):
        """cmsearch_from_file should work as expected
        """
        exp_search_res = [['a', 5, 23, 1, 19, 12.85, '-', 37],\
                          ['b', 1, 19, 1, 19, 14.359999999999999, '-', 47]]
        search_res = cmsearch_from_file(cm_file_path=self.cmfile,\
            seqs=self.seqs2_unaligned,moltype=RNA)
        self.assertEqual(search_res,exp_search_res)
        
class CmstatTests(GeneralSetUp):
    """Tests for the Cmstat application controller"""

    def test_base_command(self):
        """Infernal BaseCommand should return the correct BaseCommand"""
        c = Cmstat()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmstat']))
        c.Parameters['-g'].on()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cmstat -g']))


    def test_changing_working_dir(self):
        """Infernal BaseCommand should change according to WorkingDir"""
        c = Cmstat(WorkingDir='/tmp/cmstat_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmstat_test','/"; ','cmstat']))
        c = Cmstat()
        c.WorkingDir = '/tmp/cmstat_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cmstat_test2','/"; ','cmstat']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cmstat_test')
        rmdir('/tmp/cmstat_test2')

    def test_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)

ALN1_CM = """INFERNAL-1 [1.0rc1]
NAME     aln1-1
STATES   61
NODES    18
ALPHABET 1
ELSELF   -0.08926734
WBETA    1e-07
NSEQ     3
EFFNSEQ  3.000
CLEN     19
BCOM     cmbuild aln1.cm aln1.sto
BDATE    Sun Oct 5 18:45:35 2008
NULL     0.000  0.000  0.000  0.000 
MODEL:
				[ ROOT    0 ]
     S     0    -1 0     1     4  -2.071  -2.210  -1.649  -2.140                 
    IL     1     1 2     1     4  -0.556  -5.022  -1.818  -7.508                  0.000  0.000  0.000  0.000 
    IR     2     2 3     2     3  -0.310  -2.439  -6.805                          0.000  0.000  0.000  0.000 
				[ MATL    1 ]
    ML     3     2 3     5     3  -8.003  -0.020  -6.657                         -0.389  0.377 -1.236  0.597 
     D     4     2 3     5     3  -7.923  -3.436  -0.146                         
    IL     5     5 3     5     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL    2 ]
    ML     6     5 3     8     3  -8.003  -0.020  -6.657                          0.711 -1.015 -1.162  0.507 
     D     7     5 3     8     3  -7.923  -3.436  -0.146                         
    IL     8     8 3     8     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL    3 ]
    ML     9     8 3    11     3  -8.003  -0.020  -6.657                         -0.389  0.377 -1.236  0.597 
     D    10     8 3    11     3  -7.923  -3.436  -0.146                         
    IL    11    11 3    11     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL    4 ]
    ML    12    11 3    14     3  -8.003  -0.020  -6.657                         -0.392  0.246 -1.238  0.703 
     D    13    11 3    14     3  -7.923  -3.436  -0.146                         
    IL    14    14 3    14     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL    5 ]
    ML    15    14 3    17     3  -8.003  -0.020  -6.657                         -1.340 -2.411  1.644 -1.777 
     D    16    14 3    17     3  -7.923  -3.436  -0.146                         
    IL    17    17 3    17     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL    6 ]
    ML    18    17 3    20     3  -8.003  -0.020  -6.657                          0.830  0.106 -1.204 -0.492 
     D    19    17 3    20     3  -7.923  -3.436  -0.146                         
    IL    20    20 3    20     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL    7 ]
    ML    21    20 3    23     3  -8.003  -0.020  -6.657                         -1.143 -1.575 -1.925  1.560 
     D    22    20 3    23     3  -7.923  -3.436  -0.146                         
    IL    23    23 3    23     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL    8 ]
    ML    24    23 3    26     3  -8.391  -0.018  -6.709                          0.821 -1.044 -1.178  0.385 
     D    25    23 3    26     3  -6.905  -0.258  -2.688                         
    IL    26    26 3    26     3  -1.925  -0.554  -4.164                          0.000  0.000  0.000  0.000 
				[ MATR    9 ]
    MR    27    26 3    29     5  -7.411  -0.031  -7.227  -7.439  -8.330         -0.726  0.967 -1.567  0.142 
     D    28    26 3    29     5  -5.352  -0.707  -2.978  -4.409  -2.404         
    IR    29    29 3    29     5  -2.408  -0.496  -5.920  -4.087  -5.193          0.000  0.000  0.000  0.000 
				[ MATP   10 ]
    MP    30    29 3    34     6  -9.266  -9.205  -0.019  -7.982  -8.261  -8.656 -1.570 -1.865 -1.898  0.327 -1.331 -2.318  0.651  0.994 -1.872  0.282 -2.224 -0.666  1.972 -1.608 -0.242  1.187 
    ML    31    29 3    34     6  -6.250  -6.596  -1.310  -1.005  -6.446  -3.975  0.660 -0.612 -0.293 -0.076 
    MR    32    29 3    34     6  -6.988  -5.717  -1.625  -5.695  -0.829  -3.908  0.660 -0.612 -0.293 -0.076 
     D    33    29 3    34     6  -9.049  -7.747  -3.544  -4.226  -4.244  -0.319 
    IL    34    34 5    34     6  -2.579  -2.842  -0.760  -4.497  -5.274  -4.934  0.000  0.000  0.000  0.000 
    IR    35    35 6    35     5  -2.408  -0.496  -5.920  -4.087  -5.193          0.000  0.000  0.000  0.000 
				[ MATP   11 ]
    MP    36    35 6    40     4  -7.331  -7.538  -0.041  -5.952                 -4.114  0.397 -4.664  0.815 -4.665 -4.015 -0.462 -4.315 -3.939  3.331 -3.732 -0.830 -0.398 -3.640 -1.958 -3.517 
    ML    37    35 6    40     4  -3.758  -3.940  -0.507  -2.670                  0.660 -0.612 -0.293 -0.076 
    MR    38    35 6    40     4  -4.809  -3.838  -1.706  -0.766                  0.660 -0.612 -0.293 -0.076 
     D    39    35 6    40     4  -4.568  -4.250  -2.265  -0.520                 
    IL    40    40 5    40     4  -1.686  -2.369  -1.117  -4.855                  0.000  0.000  0.000  0.000 
    IR    41    41 6    41     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL   12 ]
    ML    42    41 6    44     5  -7.411  -0.031  -7.227  -7.439  -8.330          1.826 -2.947 -2.856 -2.413 
     D    43    41 6    44     5  -4.959  -0.803  -4.221  -2.596  -2.508         
    IL    44    44 3    44     5  -2.408  -0.496  -4.087  -5.920  -5.193          0.000  0.000  0.000  0.000 
				[ MATP   13 ]
    MP    45    44 3    49     4  -7.331  -7.538  -0.041  -5.952                 -1.592 -1.722 -1.807  0.471 -1.387 -2.146  1.822  0.774 -1.836  0.505 -2.076 -0.521  1.055 -1.515 -0.260  0.958 
    ML    46    44 3    49     4  -3.758  -3.940  -0.507  -2.670                  0.660 -0.612 -0.293 -0.076 
    MR    47    44 3    49     4  -4.809  -3.838  -1.706  -0.766                  0.660 -0.612 -0.293 -0.076 
     D    48    44 3    49     4  -4.568  -4.250  -2.265  -0.520                 
    IL    49    49 5    49     4  -1.686  -2.369  -1.117  -4.855                  0.000  0.000  0.000  0.000 
    IR    50    50 6    50     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL   14 ]
    ML    51    50 6    53     3  -8.323  -0.016  -6.977                          0.481 -1.091 -0.011  0.192 
     D    52    50 6    53     3  -6.174  -1.687  -0.566                         
    IL    53    53 3    53     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL   15 ]
    ML    54    53 3    56     3  -8.323  -0.016  -6.977                          1.148 -1.570 -0.075 -1.007 
     D    55    53 3    56     3  -6.174  -1.687  -0.566                         
    IL    56    56 3    56     3  -1.442  -0.798  -4.142                          0.000  0.000  0.000  0.000 
				[ MATL   16 ]
    ML    57    56 3    59     2       *   0.000                                 -0.726  0.967 -1.567  0.142 
     D    58    56 3    59     2       *   0.000                                 
    IL    59    59 3    59     2  -1.823  -0.479                                  0.000  0.000  0.000  0.000 
				[ END    17 ]
     E    60    59 3    -1     0                                                 
//
"""

CMALIGN_STDOUT = """# cmalign :: align sequences to an RNA CM
# INFERNAL 1.0rc1 (June 2008)
# Copyright 2007-2009 (C) 2008 HHMI Janelia Farm Research Campus
# Freely distributed under the GNU General Public License (GPL)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# command: cmalign --withali aln1.sto -o all_aligned.sto aln1.cm seqs1.fasta
# date:    Sun Oct 5 22:04:30 2008
#
# cm name                    algorithm  config  sub  bands     tau
# -------------------------  ---------  ------  ---  -----  ------
# aln1-1                       opt acc  global   no    hmm   1e-07
#
#                               bit scores                           
#                           ------------------                       
# seq idx  seq name    len     total    struct  avg prob      elapsed
# -------  --------  -----  --------  --------  --------  -----------
        1  1            23     -9.98      5.71     0.260  00:00:00.01
        2  2            13     -6.79      6.73     0.710  00:00:00.00
        3  3            17     -7.43      5.86     0.754  00:00:00.01

# Alignment saved in file all_aligned.sto.
#
# CPU time: 0.02u 0.00s 00:00:00.02 Elapsed: 00:00:00
"""


if __name__ == '__main__':
    main()
