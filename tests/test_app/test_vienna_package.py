#!/usr/bin/env python

from os import getcwd, remove, rmdir
import tempfile, shutil
from cogent.util.unit_test import TestCase, main
from cogent.app.vienna_package import RNAfold, RNAsubopt, RNAplot,\
    plot_from_seq_and_struct, DataError, get_constrained_fold

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
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

        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAfold','-d1','-T','37','-S','1.07']
        self.assertEqualItems(obs, exp)

        r.Parameters['-noLP'].on()
        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAfold','-d1','-noLP','-T','37',\
               '-S','1.07']
        self.assertEqualItems(obs, exp)

        r.Parameters['Temp'].on(15)
        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAfold','-d1','-noLP','-T','15', \
               '-S','1.07']
        self.assertEqualItems(obs, exp)

        r.Parameters['-d'].off()
        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAfold','-noLP','-T','15','-S','1.07']
        self.assertEqualItems(obs, exp)
        
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

    def test_stdout_input_as_path(self):
        """RNAfold: StdOut with input_as_path"""
        r = RNAfold(InputHandler='_input_as_path')
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

    def test_stdout_input_as_path_space(self):
        """RNAfold: StdOut with input_as_path and space in filename"""
        r = RNAfold(InputHandler='_input_as_path')
        f = open('/tmp/rna test file','w')
        f.write('\n'.join(self.mixed_seq2))
        f.close()
        exp = '\n'.join(['>namedseq2','AUAGCUAGCUAUGCGCUAGC',\
            '...((((((.....)))))) ( -8.30)','ACGGCUAUAGCUAGCGA',\
            '...((((....)))).. ( -3.20)','GCUAGCUAUUAUAUAUA',\
            '................. (  0.00)'])+'\n'
        res = r('/tmp/rna test file')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/rna test file')


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

    def test_get_constrained_fold_bad_data(self):
        """get_constrained_fold should handle bad data."""
        test_seq =   'AAACCCGGGUUU'
        constraint = '(((...)))'
        
        #Test empty sequence.
        self.assertRaises(ValueError, get_constrained_fold,\
            '', constraint)
        
        #Test empty constraint string.
        self.assertRaises(ValueError, get_constrained_fold,\
            test_seq, '')

        #Test different length sequence and constraint.
        self.assertRaises(ValueError, get_constrained_fold,\
            test_seq, constraint)        
    
    def test_get_constrained_fold(self):
        """get_constrained_fold should give correct result."""
        test_seq =   'AAACCCGGGUUU'
        constraint = '(((......)))'
        expected_struct = '((((....))))'
        
        obs_seq, obs_struct, obs_energy = \
            get_constrained_fold(test_seq,constraint)
        #Test get back correct seq and struct
        self.assertEqual(obs_seq, test_seq)
        self.assertEqual(obs_struct, expected_struct)
        
    def test_zzz_general_cleanup(self):
        """Executed last, clean up temp_dir"""
        shutil.rmtree(self.temp_dir)

class RNAsuboptTests(TestCase):
    """Tests for the RNAsubopt application controller"""

    def test_base_command(self):
        """RNAsubopt: BaseCommand should be ok for different parameter settings
        """
        r = RNAsubopt()
        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAsubopt','-e','1','-d2','-T','37']
        self.assertEqualItems(obs, exp)

        r.Parameters['-nsp'].on('GA')
        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAsubopt','-e','1','-d2','-nsp','GA',\
               '-T','37']
        self.assertEqualItems(obs, exp)

        r.Parameters['Temp'].on(15)
        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAsubopt','-e','1','-d2','-nsp','GA',\
               '-T','15']
        self.assertEqualItems(obs, exp)

        r.Parameters['-d'].off()
        obs = r.BaseCommand.split()
        exp = ['cd','"%s/";' % getcwd(),'RNAsubopt','-e','1','-nsp','GA','-T',\
               '15']
        self.assertEqualItems(obs, exp)
        
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

    def test_stdout_input_as_path_space(self):
        """RNAsubopt: StdOut with input_as_path and space in filename"""
        mixed_seq2 = ['>namedseq2','AUAGCUAGCUAUGCGCUAGC',
            'ACGGCUAUAGCUAGCGA','gcuagcuauuauauaua']
        r = RNAsubopt(InputHandler='_input_as_path')
        f = open('/tmp/rna test file','w')
        f.write('\n'.join(mixed_seq2))
        f.close()
        exp = '\n'.join(['> namedseq2 [100]',
            'AUAGCUAGCUAUGCGCUAGC   -830    100',
            '...((((((.....))))))  -8.30',
            'ACGGCUAUAGCUAGCGA   -320    100',
            '...(((......)))..  -2.30',
            '...((((....))))..  -3.20',
            'GCUAGCUAUUAUAUAUA      0    100',
            '.................   0.00'])+'\n'
        res = r('/tmp/rna test file')
        obs = res['StdOut'].read()
        self.assertEqual(obs,exp)
        res.cleanUp()
        remove('/tmp/rna test file')


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

class RNAplotTests(TestCase):
    """Tests for the RNAplot application controller"""
    
    def setUp(self):
        self.unnamed_seqs = ['AUAGCUAGCUAUGCGCUAGC','...((((((.....))))))',\
                            'ACGGCUAUAGCUAGCGA','...((((....))))..',\
                            'gcuagcuauuauauaua','.................']
            
        self.named_seqs = \
            ['>namedseq1','AUAGCUAGCUAUGCGCUAGC','...((((((.....))))))',\
             '>namedseq2','ACGGCUAUAGCUAGCGA','...((((....))))..']

        self.mixed_seqs = \
            ['>namedseq1','AUAGCUAGCUAUGCGCUAGC','...((((((.....))))))',\
             'ACGGCUAUAGCUAGCGA','...((((....))))..',\
             'gcuagcuauuauauaua','.................']
        
        self.named_seq = \
            ['>namedseq','AUAGCUAGCUAUGCGCUAGC','...((((((.....))))))']
        
        self.standard_name = 'namedseq'
        
        self.standard_seq = 'AUAGCUAGCUAUGCGCUAGC'
        
        self.standard_struct = '...((((((.....))))))'
        
        self.bad_pairing_struct = '...((((((.....)))...'
        
        self.bad_chars_struct = '...((((((.....)xx)))' 
        
        self.short_struct = '...((((((..))))))'

        self.temp_dir = '/tmp/test_rnaplot'

    def test_base_command(self):
        """RNAplot: BaseCommand should be ok for different parameter settings"""
        r = RNAplot()
        working_dir = getcwd()
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAplot']))
        r.Parameters['-t'].on(0)
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAplot -t 0']))
        r.Parameters['-o'].on('svg')
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAplot -t 0 -o svg']))
        r.Parameters['-t'].off()
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','RNAplot -o svg']))
        
    def test_changing_working_dir(self):
        """RNAplot: BaseCommand should be ok after changing the working dir"""
        #changing in initialization
        temp_dir = tempfile.mkdtemp()
        r = RNAplot(WorkingDir=temp_dir)
        self.assertEqual(r.BaseCommand,\
            'cd "%s/"; RNAplot'%(temp_dir))
        #changing afterwards
        r = RNAplot()
        r.WorkingDir = temp_dir
        self.assertEqual(r.BaseCommand,\
            'cd "%s/"; RNAplot'%(temp_dir))
        rmdir(temp_dir)

    def test_get_result_paths_unnamed_seq(self):
        """RNAplot: _get_result_paths() should work on unnamed seq"""
        r = RNAplot()
        res = r(self.unnamed_seqs)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS'])
        self.failUnless(res['SS'] is not None)
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()
    
    def test_get_result_paths_named_seq(self):
        """RNAplot: _get_result_paths() should work on named seq"""
       
        r = RNAplot()
        res = r(self.named_seqs)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','namedseq1_ss',\
            'namedseq2_ss']) 
        res.cleanUp()

    def test_get_result_paths_mixed_seq(self):
        """RNAplot: _get_result_paths() should work on partly named seq"""

        r = RNAplot()
        res = r(self.mixed_seqs)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS','namedseq1_ss']) 
        res.cleanUp()

    def test_get_result_paths_parameter(self):
        """RNAplot: _get_result_paths() should work with diff parameters"""

        r = RNAplot()
        r.Parameters['-t'].on(0)
        res = r(self.unnamed_seqs)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS'])
        self.failUnless(res['SS'] is not None)
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()

    def test_get_result_paths_working_dir(self):
        """RNAplot: _get_result_paths() should work with diff working dir"""
        r = RNAplot(WorkingDir=self.temp_dir)
        res = r(self.unnamed_seqs)
        self.assertEqualItems(res.keys(),\
            ['StdOut','StdErr','ExitStatus','SS'])
        self.failUnless(res['SS'] is not None)
        self.failUnless(isinstance(res['SS'],file))
        self.failUnless(res['StdOut'] is not None)
        self.failUnless(res['StdErr'] is None)
        self.assertEqual(res['ExitStatus'],0)
        res.cleanUp()
    
    def test_rnaplot_output(self):
        """RNAplot: calling RNAplot on data should give correct result."""
        r = RNAplot(WorkingDir=self.temp_dir)
        res = r(self.named_seq)
        ps_plot = res['namedseq_ss'].read()
        observed_lines = ps_plot.split('\n')
        expected_lines = RNAPLOT_RES.split('\n')
        #First 8 lines depend on the runtime.  Check after.
        self.assertEqual(observed_lines[8:],expected_lines[8:])
        res.cleanUp()
    
    def test_plot_from_seq_and_struct_bad_data(self):
        """plot_from_seq_and_struct helper function should handle bad data.
        """
        #Check bad Vienna Structure pairing
        self.assertRaises(IndexError, plot_from_seq_and_struct,\
            self.standard_seq, self.bad_pairing_struct)
        
        #Check bad characters in Vienna Structure
        self.assertRaises(ValueError, plot_from_seq_and_struct,\
            self.standard_seq, self.bad_chars_struct)
        
        #Check different lengths of seq and struct
        self.assertRaises(DataError, plot_from_seq_and_struct,\
            self.standard_seq, self.short_struct)
    
    def test_plot_from_seq_and_struct(self):
        """plot_from_seq_and_struct helper function should give correct result.
        """
        ps_plot = plot_from_seq_and_struct(self.standard_seq,\
            self.standard_struct, seqname=self.standard_name)
        observed_lines = ps_plot.split('\n')
        expected_lines = RNAPLOT_RES.split('\n')
        #First 8 lines depend on the runtime.  Check after.
        self.assertEqual(observed_lines[8:],expected_lines[8:])

    def test_zzz_general_cleanup(self):
        """Executed last, clean up temp_dir"""
        shutil.rmtree(self.temp_dir)

RNAPLOT_RES = """%!PS-Adobe-3.0 EPSF-3.0
%%Creator: PS_dot.c,v 1.38 2007/02/02 15:18:13 ivo Exp $, ViennaRNA-1.7.1
%%CreationDate: Tue Oct 21 15:44:50 2008
%%Title: RNA Secondary Structure Plot
%%BoundingBox: 66 210 518 662
%%DocumentFonts: Helvetica
%%Pages: 1
%%EndComments

%Options: 
% to switch off outline pairs of sequence comment or
% delete the appropriate line near the end of the file

%%BeginProlog
/RNAplot 100 dict def
RNAplot begin
/fsize  14 def
/outlinecolor {0.2 setgray} bind def
/paircolor    {0.2 setgray} bind def
/seqcolor     {0   setgray} bind def
/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def
/min { 2 copy gt { exch } if pop } bind def
/max { 2 copy lt { exch } if pop } bind def
/drawoutline {
  gsave outlinecolor newpath
  coor 0 get aload pop 0.8 0 360 arc % draw 5' circle of 1st sequence
  currentdict /cutpoint known        % check if cutpoint is defined
  {coor 0 cutpoint getinterval
   {aload pop lineto} forall         % draw outline of 1st sequence
   coor cutpoint get aload pop
   2 copy moveto 0.8 0 360 arc       % draw 5' circle of 2nd sequence
   coor cutpoint coor length cutpoint sub getinterval
   {aload pop lineto} forall}        % draw outline of 2nd sequence
  {coor {aload pop lineto} forall}   % draw outline as a whole
  ifelse
  stroke grestore
} bind def
/drawpairs {
  paircolor
  0.7 setlinewidth
  [9 3.01] 9 setdash
  newpath
  pairs {aload pop
     coor exch 1 sub get aload pop moveto
     coor exch 1 sub get aload pop lineto
  } forall
  stroke
} bind def
% draw bases
/drawbases {
  [] 0 setdash
  seqcolor
  0
  coor {
    aload pop moveto
    dup sequence exch 1 getinterval cshow
    1 add
  } forall
  pop
} bind def

/init {
  /Helvetica findfont fsize scalefont setfont
  1 setlinejoin
  1 setlinecap
  0.8 setlinewidth
  72 216 translate
  % find the coordinate range
  /xmax -1000 def /xmin 10000 def
  /ymax -1000 def /ymin 10000 def
  coor {
      aload pop
      dup ymin lt {dup /ymin exch def} if
      dup ymax gt {/ymax exch def} {pop} ifelse
      dup xmin lt {dup /xmin exch def} if
      dup xmax gt {/xmax exch def} {pop} ifelse
  } forall
  /size {xmax xmin sub ymax ymin sub max} bind def
  72 6 mul size div dup scale
  size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div
  translate
} bind def
end
%%EndProlog
RNAplot begin
% data start here
/sequence (\\
AUAGCUAGCUAUGCGCUAGC\\
) def
/coor [
[100.992 114.290]
[88.210 108.134]
[86.994 93.999]
[98.538 85.751]
[105.046 72.236]
[111.554 58.722]
[118.062 45.207]
[124.571 31.693]
[131.079 18.178]
[127.159 2.622]
[136.992 -10.055]
[153.034 -10.127]
[162.981 2.461]
[159.200 18.052]
[144.593 24.686]
[138.085 38.201]
[131.577 51.716]
[125.069 65.230]
[118.560 78.745]
[112.052 92.259]
] def
/pairs [
[4 20]
[5 19]
[6 18]
[7 17]
[8 16]
[9 15]
] def

init

% switch off outline pairs or bases by removing these lines
drawoutline
drawpairs
drawbases
% show it
showpage
end
%%EOF
"""

if __name__ == '__main__':
    main()
