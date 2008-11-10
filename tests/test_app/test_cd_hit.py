#!/usr/bin/env python

from os import getcwd, rmdir
from cogent.core.moltype import PROTEIN, DNA
from cogent.util.unit_test import TestCase, main
from cogent.app.cd_hit import CD_HIT, CD_HIT_EST, cdhit_from_seqs

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class CD_HIT_Tests(TestCase):
    """Tests for the CD-HIT application controller"""

    def test_base_command(self):
        """CD_HIT BaseCommand should return the correct BaseCommand"""
        c = CD_HIT()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cd-hit']))
        c.Parameters['-i'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cd-hit -i "seq.txt"']))
        c.Parameters['-c'].on(0.8)
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cd-hit -c 0.8' +
            ' -i "seq.txt"']))

    def test_changing_working_dir(self):
        """CD_HIT BaseCommand should change according to WorkingDir"""
        c = CD_HIT(WorkingDir='/tmp/cdhit_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cdhit_test','/"; ','cd-hit']))
        c = CD_HIT()
        c.WorkingDir = '/tmp/cdhit_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cdhit_test2','/"; ','cd-hit']))

        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cdhit_test')
        rmdir('/tmp/cdhit_test2')

    def test_cdhit_from_seqs(self):
        """CD_HIT should return expected seqs"""
        res = cdhit_from_seqs(protein_seqs, PROTEIN, {'-c':0.8})
        self.assertEqual(res.toFasta(), protein_expected)

class CD_HIT_EST_Tests(TestCase):
    """Tests for the CD-HIT application controller"""

    def test_base_command(self):
        """CD_HIT_EST BaseCommand should return the correct BaseCommand"""
        c = CD_HIT_EST()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cd-hit-est']))
        c.Parameters['-i'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cd-hit-est -i "seq.txt"']))
        c.Parameters['-c'].on(0.8)
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','cd-hit-est -c 0.8' +
            ' -i "seq.txt"']))

    def test_changing_working_dir(self):
        """CD_HIT_EST BaseCommand should change according to WorkingDir"""
        c = CD_HIT_EST(WorkingDir='/tmp/cdhitest_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cdhitest_test','/"; ','cd-hit-est']))
        c = CD_HIT_EST()
        c.WorkingDir = '/tmp/cdhitest_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/cdhitest_test2','/"; ','cd-hit-est']))

        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/cdhitest_test')
        rmdir('/tmp/cdhitest_test2')

    def test_cdhit_from_seqs(self):
        """CD_HIT should return expected seqs"""
        res = cdhit_from_seqs(dna_seqs, DNA, {'-c':0.8})
        self.assertEqual(res.toFasta(), dna_expected)

dna_seqs = """>cdhit_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>cdhit_test_seqs_1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>cdhit_test_seqs_2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>cdhit_test_seqs_3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>cdhit_test_seqs_4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>cdhit_test_seqs_5
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>cdhit_test_seqs_6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>cdhit_test_seqs_7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>cdhit_test_seqs_8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>cdhit_test_seqs_9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA"""

dna_expected = """>cdhit_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>cdhit_test_seqs_1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>cdhit_test_seqs_2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>cdhit_test_seqs_4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>cdhit_test_seqs_5
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>cdhit_test_seqs_7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT"""

protein_seqs = """>seq1
MGNKWSKSWPQVRDRMRRAAPAPAADGVGAVSQDLAKHGAITSSNTAATNDDCAWLEAQTEEEVGFPVRPQVPLRPMTYK
>seq2
MGGKWSKSSIVGWSTVRERMRKTPPAADGVGAVSQDLDKHGAVTSSNTAFNNPDCAWLEAQEDEDVGFPVRPQVPLRPT
>seq3
MGGKWSKSSIVGWPAIRERMRRARPAADRVGTQPAADGVGAVSQDLARHGAVTSSNTSHNNPDCAWLEAQEEEEVGVR
>seq4
MGKIWSKSSIVGWPEIRERMRRQRPHEPAVEPAVGVGAASQDLANRGALTTSNTRTNNPTVAWVEAQEEEGEVVRPQ
>seq5
MGKIWSKSSLVGWPEIRERMRRQTQEPAVEPAVGAGAASQDLANRGAITIRNTRDNNESIAWLEAQEEEFPVRPQV
>seq6
MGKIWSKSSLVGWPEIRERIRRQTPEPAVGVGAVSQDLANRGAITTSNTKDNNQTVAWLEAQEEPVRPQVPLRPM
>seq7
MGNALRKGKFEGWAAVRERMRRTRTFPESEPCAPGVGQISRELAARGGIPSSHTPQNNESHQEEEVGFPVAPQV
>seq8
MGNAWSKSKFAGWSEVRDRMRRSSSDPQQPCAPGVGAVSRELATRGGISSSALAFLDSHKDEDVGFPVRPQVP
>seq9
MGNVLGKDKFKGWAAVRERMRKTSSDPDPQPCAPGVGPVSRELSYTPQNNAALAFLESHEDEDVGFPVXPQV
>seq10
MGNVLGKDKFKGWSAVRERMRKTSPEPEPCAPGVRGGISNSHTPQNNAALAFLESHQDEDVGFPVRPQVPL"""

protein_expected = """>seq1
MGNKWSKSWPQVRDRMRRAAPAPAADGVGAVSQDLAKHGAITSSNTAATNDDCAWLEAQTEEEVGFPVRPQVPLRPMTYK
>seq2
MGGKWSKSSIVGWSTVRERMRKTPPAADGVGAVSQDLDKHGAVTSSNTAFNNPDCAWLEAQEDEDVGFPVRPQVPLRPT
>seq3
MGGKWSKSSIVGWPAIRERMRRARPAADRVGTQPAADGVGAVSQDLARHGAVTSSNTSHNNPDCAWLEAQEEEEVGVR
>seq4
MGKIWSKSSIVGWPEIRERMRRQRPHEPAVEPAVGVGAASQDLANRGALTTSNTRTNNPTVAWVEAQEEEGEVVRPQ
>seq5
MGKIWSKSSLVGWPEIRERMRRQTQEPAVEPAVGAGAASQDLANRGAITIRNTRDNNESIAWLEAQEEEFPVRPQV
>seq7
MGNALRKGKFEGWAAVRERMRRTRTFPESEPCAPGVGQISRELAARGGIPSSHTPQNNESHQEEEVGFPVAPQV
>seq8
MGNAWSKSKFAGWSEVRDRMRRSSSDPQQPCAPGVGAVSRELATRGGISSSALAFLDSHKDEDVGFPVRPQVP
>seq9
MGNVLGKDKFKGWAAVRERMRKTSSDPDPQPCAPGVGPVSRELSYTPQNNAALAFLESHEDEDVGFPVXPQV"""

if __name__ == '__main__':
    main()
