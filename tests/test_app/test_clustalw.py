#!/usr/bin/env python
"""Tests for application controller for ClustalW v1.83"""
import re
from os import getcwd, remove, rmdir, mkdir, path
import shutil
from cogent.core.alignment import Alignment
from cogent.core.moltype import RNA
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from cogent.app.clustalw import Clustalw, alignUnalignedSeqsFromFile,\
    alignUnalignedSeqs, alignTwoAlignments, addSeqsToAlignment,\
    buildTreeFromAlignment, build_tree_from_alignment, \
    bootstrap_tree_from_alignment, align_unaligned_seqs, \
    align_and_build_tree, add_seqs_to_alignment, align_two_alignments
from cogent.parse.fasta import MinimalFastaParser

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight", "Daniel McDonald",\
               "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"


cw_vers = re.compile("CLUSTAL W [(]1\.8[1-3][.\d]*[)]")

class GeneralSetUp(TestCase):

    def setUp(self):
        """Clustalw general setUp method for all tests"""
        self.seqs1 = ['ACUGCUAGCUAGUAGCGUACGUA','GCUACGUAGCUAC',
            'GCGGCUAUUAGAUCGUA']
        self.aln1_fasta = ALIGN1_FASTA
        self.labels1 = ['>1','>2','>3']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))
        self.stdout1 = STDOUT1
        self.aln1 = ALIGN1
        self.dnd1 = DND1
        
        self.multiline1 = '\n'.join(flatten(zip(self.labels1, self.seqs1)))
       
        self.seqs2=['UAGGCUCUGAUAUAAUAGCUCUC','UAUCGCUUCGACGAUUCUCUGAUAGAGA',
            'UGACUACGCAU']
        self.labels2=['>a','>b','>c']
        self.lines2 = flatten(zip(self.labels2,self.seqs2))
        self.aln2 = ALIGN2
        self.dnd2 = DND2
        
        self.twoalign = TWOALIGN
        self.alignseqs = ALIGNSEQS
        self.treeduringalignseqs = TREEDURINGALIGNSEQS
        self.treefromalignseqs = TREEFROMALIGNSEQS
        
        self.temp_dir_space = "/tmp/clustalw test"

        self.build_tree_seqs_short = """>clustal_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AGCTTTAAATCATGCCAGTG
>clustal_test_seqs_1
GACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
TGCTTTCAATAATGCCAGTG
>clustal_test_seqs_2
AACCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
TGCTTTGAATCATGCCAGTA
>clustal_test_seqs_3
AAACCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
TGCTTTACATCATGCAAGTG
>clustal_test_seqs_4
AACCGCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
TGCTTTAAATCATGCCAGTG
>clustal_test_seqs_5
AACCCCCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
TGCTTTAAATCATGCCAGTT
>clustal_test_seqs_6
GACCCCCGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
TACTTTAGATCATGCCGGTG
>clustal_test_seqs_7
AACCCCCACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
TGCTTTAAATCATGCCAGTG
>clustal_test_seqs_8
AACCCCCACGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
TGCATTAAATCATGCCAGTG
>clustal_test_seqs_9
AAGCCCCACGGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
TGCTTTAAATCCTGACAGCG
"""
        self.build_tree_seqs_long = """>clustal_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AGCTTTAAATCATGCCAGTG
>clustal_test_seqsaaaaaaaa_1
GACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
TGCTTTCAATAATGCCAGTG
>clustal_test_seqsaaaaaaaa_2
AACCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
TGCTTTGAATCATGCCAGTA
>clustal_test_seqsaaaaaaaa_3
AAACCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
TGCTTTACATCATGCAAGTG
>clustal_test_seqsaaaaaaaa_4
AACCGCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
TGCTTTAAATCATGCCAGTG
>clustal_test_seqsaaaaaaaa_5
AACCCCCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
TGCTTTAAATCATGCCAGTT
>clustal_test_seqsaaaaaaaa_6
GACCCCCGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
TACTTTAGATCATGCCGGTG
>clustal_test_seqsaaaaaaaa_7
AACCCCCACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
TGCTTTAAATCATGCCAGTG
>clustal_test_seqsaaaaaaaa_8
AACCCCCACGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
TGCATTAAATCATGCCAGTG
>clustal_test_seqsaaaaaaaa_9
AAGCCCCACGGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
TGCTTTAAATCCTGACAGCG
"""
        try:
            mkdir('/tmp/ct')
        except OSError: #dir already exists
            pass
        
        try:
            #create sequence files
            f = open('/tmp/ct/seq1.txt','w')
            f.write('\n'.join(self.lines1))
            f.close()
            g = open('/tmp/ct/seq2.txt','w')
            g.write('\n'.join(self.lines2))
            g.close()
            #create alignment files
            f = open('/tmp/ct/align1','w')
            f.write(self.aln1)
            f.close()
            g = open('/tmp/ct/align2','w')
            g.write(self.aln2)
            g.close()
            #create tree file
            f = open('/tmp/ct/tree1','w')
            f.write(DND1)
            f.close()
        except OSError:
            pass



class ClustalwTests(GeneralSetUp):
    """Tests for the Clustalw application controller"""

    def test_base_command(self):
        """Clustalw BaseCommand should return the correct BaseCommand"""
        c = Clustalw()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','clustalw -align']))
        c.Parameters['-infile'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ',\
            'clustalw -infile="seq.txt" -align']))
        c.Parameters['-align'].off()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','clustalw -infile="seq.txt"']))
        c.Parameters['-nopgap'].on()
        c.Parameters['-infile'].off()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','clustalw -nopgap']))

    def test_changing_working_dir(self):
        """Clustalw BaseCommand should change according to WorkingDir"""
        c = Clustalw(WorkingDir='/tmp/clustaltest')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/clustaltest','/"; ','clustalw -align']))
        c = Clustalw(WorkingDir='/tmp/clustaltest/')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/clustaltest/','/"; ','clustalw -align']))
        c = Clustalw()
        c.WorkingDir = '/tmp/clustaltest2/'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/clustaltest2/','/"; ','clustalw -align']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/clustaltest')
        rmdir('/tmp/clustaltest2')
    
    def test_stdout_input_as_string(self):
        """Clustalw input_as_string shoud function as expected"""
        c = Clustalw(WorkingDir='/tmp/ct')
        res = c('/tmp/ct/seq1.txt')
        self.assertEqual(cw_vers.sub("", res['StdOut'].read()),
                            cw_vers.sub("", self.stdout1))
        self.assertEqual(res['StdErr'].read(),'')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()

    def test_stdout_input_as_lines(self):
        """Clustalw input_as_lines should function as expected"""
        c = Clustalw(InputHandler='_input_as_lines',WorkingDir='/tmp/ct')
        res = c(self.lines1)
        #get info on input file name and change output accordingly
        name = c.Parameters['-infile'].Value
        out = self.stdout1.split('\n')
        out[16] =\
            'Guide tree        file created:   ['+name.rsplit(".")[0]+'.dnd]'
        out[23] =\
            'CLUSTAL-Alignment file created  ['+name.rsplit(".")[0]+'.aln]'

        self.assertEqual(cw_vers.sub("", res['StdOut'].read()),
                            cw_vers.sub("", '\n'.join(out)))
        self.assertEqual(res['StdErr'].read(),'')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()

    def test_stdout_input_as_lines_local(self):
        """Clustalw input_as_lines should function as expected"""
        c = Clustalw(InputHandler='_input_as_lines',WorkingDir=self.temp_dir_space)
        res = c(self.lines1)
        #get info on input file name and change output accordingly
        name = c.Parameters['-infile'].Value
        out = self.stdout1.split('\n')
        out[16] =\
            'Guide tree        file created:   ['+name.rsplit(".")[0]+'.dnd]'
        out[23] =\
            'CLUSTAL-Alignment file created  ['+name.rsplit(".")[0]+'.aln]'

        self.assertEqual(cw_vers.sub("", res['StdOut'].read()),
                            cw_vers.sub("", '\n'.join(out)))
        self.assertEqual(res['StdErr'].read(),'')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()

    def test_stdout_input_as_seqs(self):
        """Clustalw input_as_seqs should function as expected"""
        c = Clustalw(InputHandler='_input_as_seqs',WorkingDir='/tmp/ct')
        res = c(self.seqs1)
        #get info on input file name and change output accordingly
        name = c.Parameters['-infile'].Value
        out = self.stdout1.split('\n')
        out[16] =\
            'Guide tree        file created:   ['+name.rsplit(".")[0]+'.dnd]'
        out[23] =\
            'CLUSTAL-Alignment file created  ['+name.rsplit(".")[0]+'.aln]'
        
        self.assertEqual(cw_vers.sub("", res['StdOut'].read()),
                            cw_vers.sub("", '\n'.join(out)))
        self.assertEqual(res['StdErr'].read(),'')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()

    def test_stdout_input_as_multiline_string(self):
        """Clustalw input_as_multiline_string should function as expected"""
        c = Clustalw(InputHandler='_input_as_multiline_string',\
                     WorkingDir='/tmp/ct')
        res = c(self.multiline1)
        name = c.Parameters['-infile'].Value
        out = self.stdout1.split('\n')
        out[16] =\
            'Guide tree        file created:   ['+name.rsplit(".")[0]+'.dnd]'
        out[23] =\
            'CLUSTAL-Alignment file created  ['+name.rsplit(".")[0]+'.aln]'
        
        self.assertEqual(cw_vers.sub("", res['StdOut'].read()),
                            cw_vers.sub("", '\n'.join(out)))
        self.assertEqual(res['StdErr'].read(),'')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()

    def test_alignment_trees(self):
        """Clustalw alignment should work correctly with new/usetree"""
        c = Clustalw(params={'-quicktree':True,'-type':'DNA','-gapopen':10},\
            WorkingDir='/tmp/ct')
        res = c('/tmp/ct/seq1.txt')
        self.assertEqual(res['Align'].name,'/tmp/ct/seq1.aln')
        self.assertEqual(res['Dendro'].name,'/tmp/ct/seq1.dnd')
        res.cleanUp()
        c.Parameters['-usetree'].on('/tmp/ct/tree1')
        c.Parameters['-output'].on('PHYLIP')
        res = c('/tmp/ct/seq1.txt')
        self.assertEqual(res['Align'].name,'/tmp/ct/seq1.phy')
        self.assertEqual(res['Dendro'].name,'/tmp/ct/tree1')
        res.cleanUp()
        c.Parameters['-newtree'].on('newtree')
        c.Parameters['-outfile'].on('outfile')
        res = c('/tmp/ct/seq1.txt')
        self.assertEqual(res['Align'].name, c.WorkingDir + 'outfile')
        self.assertEqual(res['Dendro'].name, c.WorkingDir + 'newtree')
        res.cleanUp()
    
    def test_profile_newtree(self):
        """Clustalw profile should work correctly with new/usetree"""
        c = Clustalw(params={'-profile':None,'-profile1':'/tmp/ct/seq1.txt',\
            '-profile2':'/tmp/ct/seq2.txt','-newtree1':'lala'},\
            WorkingDir='/tmp/ct')
        c.Parameters['-align'].off()
        res = c()
        self.assertEqual(res['Align'],None)
        self.assertEqual(res['Dendro1'].name,'/tmp/ct/lala')
        self.assertEqual(res['Dendro2'].name,'/tmp/ct/seq2.dnd')
        res.cleanUp()

    def test_sequences_newtree(self):
        """Clustalw sequences should work correctly with new/usetree"""
        c = Clustalw(params={'-sequences':None,'-newtree':'lala',\
            '-profile1':'/tmp/ct/align1','-profile2':'/tmp/ct/seq2.txt'},\
            WorkingDir='/tmp/ct')
        c.Parameters['-align'].off()
        res = c()
        self.assertEqual(res['Align'],None)
        self.assertEqual(res['Dendro'].name,'/tmp/ct/lala')
        res.cleanUp()
        
        #is this a bug in clustal. It's creating an empty file 'seq2.aln'
        #but doesn't report it in the stdout
        remove('/tmp/ct/seq2.aln')

    def test_tree_outputtree(self):
        """Clustalw tree should work correctly with outputtree"""
        c = Clustalw(params={'-tree':None,'-outputtree':'dist',\
            '-infile':'/tmp/ct/align1'},WorkingDir='/tmp/ct/')
        c.Parameters['-align'].off()
        res = c()
        self.assertEqual(res['Tree'].name,'/tmp/ct/align1.ph')
        self.assertEqual(res['TreeInfo'].name,'/tmp/ct/align1.dst')
        res.cleanUp()


class clustalwTests(GeneralSetUp):
    """Tests for module level functions in clustalw.py"""
   
       
    def test_alignUnalignedSeqs(self):
        """Clustalw alignUnalignedSeqs should work as expected"""
        res = alignUnalignedSeqs(self.seqs1,WorkingDir='/tmp/ct')
        self.assertNotEqual(res['StdErr'],None)
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()
        
        #suppress stderr and stdout
        res = alignUnalignedSeqs(self.seqs1,WorkingDir='/tmp/ct',\
            SuppressStderr=True,SuppressStdout=True)
        self.assertEqual(res['StdOut'],None)
        self.assertEqual(res['StdErr'],None)
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()
    
    def test_alignUnalignedSeqsFromFile(self):
        """Clustalw alignUnalignedSeqsFromFile should work as expected"""
        #make temp file
        res = alignUnalignedSeqsFromFile('/tmp/ct/seq1.txt')
        self.assertEqual(cw_vers.sub("", res['StdOut'].read()),
                            cw_vers.sub("", self.stdout1))
        self.assertEqual(res['StdErr'].read(),'')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()
        #suppress stderr and stdout
        res = alignUnalignedSeqsFromFile('/tmp/ct/seq1.txt',\
            SuppressStderr=True, SuppressStdout=True)
        self.assertEqual(res['StdOut'],None)
        self.assertEqual(res['StdErr'],None)
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.aln1))
        self.assertEqual(res['Dendro'].read(),self.dnd1)
        res.cleanUp()

    def test_alignTwoAlignments(self):
        """Clustalw alignTwoAlignments should work as expected"""
        res = alignTwoAlignments('/tmp/ct/align1','/tmp/ct/align2',\
            'twoalign.aln')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.twoalign))
        self.assertNotEqual(res['Dendro1'],None)
        self.assertNotEqual(res['Dendro2'],None)
        #are there new trees created during the profiling?
        #the produced trees are not the same as when aligning individually
        #self.assertEqual(res['Dendro1'].read(),self.dnd)
        #self.assertEqual(res['Dendro2'].read(),self.dnd2)
        res.cleanUp() 
       
    def test_addSeqsToAlignment(self):
        """Clustalw addSeqsToAlignment shoudl work as expected"""
        res = addSeqsToAlignment('/tmp/ct/align1','/tmp/ct/seq2.txt',\
            'alignseqs')
        self.assertEqual(cw_vers.sub("", res['Align'].read()),
                            cw_vers.sub("", self.alignseqs))
        self.assertEqual(res['Dendro'].read(),self.treeduringalignseqs)
        res.cleanUp()
        
    def test_buildTreeFromAlignment(self):
        """Clustalw buildTreeFromAlignment shoudl work as expected"""
        pre_res = addSeqsToAlignment('/tmp/ct/align1','/tmp/ct/seq2.txt',\
            'alignseqs',WorkingDir='/tmp/ct')
        res = buildTreeFromAlignment('/tmp/ct/alignseqs',WorkingDir='/tmp/ct')
        self.assertEqual(res['Tree'].read(),self.treefromalignseqs)
        
        res.cleanUp()
        pre_res.cleanUp()

    def test_build_tree_from_alignment(self):
        """Clustalw should return a tree built from the passed alignment"""
        tree_short = build_tree_from_alignment(self.build_tree_seqs_short, \
                RNA, best_tree=False)
        num_seqs = flatten(self.build_tree_seqs_short).count('>')
        self.assertEqual(len(tree_short.tips()), num_seqs)
        
        tree_long = build_tree_from_alignment(self.build_tree_seqs_long, \
                RNA, best_tree=False)
        seq_names = []
        for line in self.build_tree_seqs_long.split('\n'):
            if line.startswith('>'):
                seq_names.append(line[1:])

        for node in tree_long.tips():
            if node.Name not in seq_names:
                self.fail()

        tree_short = build_tree_from_alignment(self.build_tree_seqs_short, \
                RNA, best_tree=True, params={'-bootstrap':3})
        num_seqs = flatten(self.build_tree_seqs_short).count('>')
        self.assertEqual(len(tree_short.tips()), num_seqs)
     
    def test_align_unaligned_seqs(self):
        """Clustalw align_unaligned_seqs should work as expected"""
        res = align_unaligned_seqs(self.seqs1, RNA)
        self.assertEqual(res.toFasta(), self.aln1_fasta)
        
    def test_bootstrap_tree_from_alignment(self):
        """Clustalw should return a bootstrapped tree from the passed aln"""
        tree_short = bootstrap_tree_from_alignment(self.build_tree_seqs_short)
        num_seqs = flatten(self.build_tree_seqs_short).count('>')
        self.assertEqual(len(tree_short.tips()), num_seqs)
        
        tree_long = bootstrap_tree_from_alignment(self.build_tree_seqs_long)
        seq_names = []
        for line in self.build_tree_seqs_long.split('\n'):
            if line.startswith('>'):
                seq_names.append(line[1:])

        for node in tree_long.tips():
            if node.Name not in seq_names:
                self.fail()
    def test_align_and_build_tree(self):
        """Aligns and builds a tree for a set of sequences"""
        res = align_and_build_tree(self.seqs1, RNA)
        self.assertEqual(res['Align'].toFasta(), self.aln1_fasta)

        tree = res['Tree']
        seq_names = []
        for line in self.aln1_fasta.split('\n'):
            if line.startswith('>'):
                seq_names.append(line[1:])
                
        for node in tree.tips():
            if node.Name not in seq_names:
                self.fail()
    
    def test_add_seqs_to_alignment(self):
        """Clustalw add_seqs_to_alignment should work as expected."""
        seq2 = dict(MinimalFastaParser(self.lines2))
        align1 = dict(MinimalFastaParser(ALIGN1_FASTA.split('\n')))
        res = add_seqs_to_alignment(seq2,align1,RNA)
        self.assertEqual(res.toFasta(), SEQ_PROFILE_ALIGN)
    
    def test_align_two_alignments(self):
        """Clustalw align_two_alignments should work as expected."""
        align1 = dict(MinimalFastaParser(ALIGN1_FASTA.split('\n')))
        align2 = dict(MinimalFastaParser(ALIGN2_FASTA.split('\n')))
        res = align_two_alignments(align1,align2,RNA)
        self.assertEqual(res.toFasta(), PROFILE_PROFILE_ALIGN)
    
    def test_zzz_general_cleanUp(self):
        """Last test executed: cleans up all files initially created"""
        remove('/tmp/ct/seq1.txt')
        remove('/tmp/ct/seq2.txt')
        remove('/tmp/ct/align1')
        remove('/tmp/ct/align2')
        remove('/tmp/ct/tree1')
        rmdir('/tmp/ct')
        shutil.rmtree(self.temp_dir_space)

STDOUT1=\
"""


 CLUSTAL W (1.83) Multiple Sequence Alignments



Sequence format is Pearson
Sequence 1: 1                23 bp
Sequence 2: 2                13 bp
Sequence 3: 3                17 bp
Start of Pairwise alignments
Aligning...
Sequences (1:2) Aligned. Score:  46
Sequences (1:3) Aligned. Score:  41
Sequences (2:3) Aligned. Score:  30
Guide tree        file created:   [/tmp/ct/seq1.dnd]
Start of Multiple Alignment
There are 2 groups
Aligning...
Group 1: Sequences:   2      Score:171
Group 2: Sequences:   3      Score:162
Alignment Score 33
CLUSTAL-Alignment file created  [/tmp/ct/seq1.aln]
""" 

ALIGN1=\
"""CLUSTAL W (1.83) multiple sequence alignment


1               ACUGCUAGCUAGUAGCGUACGUA
2               ---GCUACGUAGCUAC-------
3               GCGGCUAUUAGAUCGUA------
                   ****                
"""

ALIGN1_FASTA = ">seq_0\nACUGCUAGCUAGUAGCGUACGUA\n>seq_1\n---GCUACGUAGCUAC-------\n>seq_2\nGCGGCUAUUAGAUCGUA------"

DND1=\
"""(
1:0.21719,
2:0.32127,
3:0.37104);
"""

ALIGN2 =\
"""CLUSTAL W (1.83) multiple sequence alignment


a               UAGGCUCUGAUAUAAUAGCUCUC---------
b               ----UAUCGCUUCGACGAUUCUCUGAUAGAGA
c               ------------UGACUACGCAU---------
                              *     *           
"""

ALIGN2_FASTA = ">a\nUAGGCUCUGAUAUAAUAGCUCUC---------\n>b\n----UAUCGCUUCGACGAUUCUCUGAUAGAGA\n>c\n------------UGACUACGCAU---------"

DND2=\
"""(
a:0.30435,
b:0.30435,
c:0.33202);
"""

TWOALIGN=\
"""CLUSTAL W (1.83) multiple sequence alignment


1               ---ACUGCUAGCUAGUAGCGUACGUA------
2               ------GCUACGUAGCUAC-------------
3               ---GCGGCUAUUAGAUCGUA------------
a               UAGGCUCUGAUAUAAUAGCUCUC---------
b               ----UAUCGCUUCGACGAUUCUCUGAUAGAGA
c               ------------UGACUACGCAU---------
                                                
"""

ALIGNSEQS=\
"""CLUSTAL W (1.83) multiple sequence alignment


1               ----------ACUGCUAGCUAGUAGCGUACGUA
2               -------------GCUACGUAGCUAC-------
3               ----------GCGGCUAUUAGAUCGUA------
a               -------UAGGCUCUGAUAUAAUAGCUCUC---
c               -------------------UGACUACGCAU---
b               UAUCGCUUCGACGAUUCUCUGAUAGAGA-----
                                                 
"""

TREEDURINGALIGNSEQS=\
"""(
1:0.34511,
(
2:0.25283,
(
(
3:0.21486,
a:0.19691)
:0.11084,
b:0.31115)
:0.06785)
:0.02780,
c:0.20035);
"""

TREEFROMALIGNSEQS=\
"""(
(
(
1:0.17223,
(
2:0.14749,
c:0.13822)
:0.19541)
:0.07161,
a:0.25531)
:0.03600,
3:0.29438,
b:0.23503);
"""

SEQ_PROFILE_ALIGN = """>a\n-------UAGGCUCUGAUAUAAUAGCUCUC---\n>b\nUAUCGCUUCGACGAUUCUCUGAUAGAGA-----\n>c\n-------------------UGACUACGCAU---\n>seq_0\n----------ACUGCUAGCUAGUAGCGUACGUA\n>seq_1\n-------------GCUACGUAGCUAC-------\n>seq_2\n----------GCGGCUAUUAGAUCGUA------"""

PROFILE_PROFILE_ALIGN = """>a\nUAGGCUCUGAUAUAAUAGCUCUC---------\n>b\n----UAUCGCUUCGACGAUUCUCUGAUAGAGA\n>c\n------------UGACUACGCAU---------\n>seq_0\n---ACUGCUAGCUAGUAGCGUACGUA------\n>seq_1\n------GCUACGUAGCUAC-------------\n>seq_2\n---GCGGCUAUUAGAUCGUA------------"""

if __name__ == '__main__':
    main()
