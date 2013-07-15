#! /usr/bin/env python
# file: test_seqprep.py

# Tests for the seqprep.py application controller.
# https://github.com/jstjohn/SeqPrep

from cogent.util.unit_test import TestCase, main
from cogent.app.fastq_join import FastqJoin, run_default_fastqjoin
from os import getcwd, remove, rmdir, mkdir, path

__author__ = "Michael Robeson"
__copyright__ = "Copyright 2007-2013, The Cogent Project"
__credits__ = ["Michael Robeson"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Michael Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"


class GenericSeqPrep(TestCase):
    """General setup for SeqPrep tests """
    
    # make directory test
    self.temp_dir = '/tmp/test_for_seqprep'
    try:
        mkdir(self.temp_dir)
    except OSError:
        pass

    # make directory with spaces test
    self.temp_dir_spaces = '/tmp/test for seqprep'
    try:
        mkdir(self.temp_dir_spaces)
    except OSError:
        pass

    # create fastq files
    try:
       p1 = path.join(self.temp_dir,'reads1.fastq')
       reads1 = open(p1,'w')
       reads1.write(reads1_string) # bottom of file
       reads1.close()
       self.test_fn1 = self.temp_dir + 'reads1.fastq'
       self.test_fn1_space = self.temp_dir_spaces + 'reads1.fastq'

       p2 = path.join(self.temp_dir,'reads2.fastq')
       reads2 = open(p2,'w')
       reads2.write(reads2_string) # bottom of file
       reads2.close()
       self.test_fn2 = self.temp_dir + 'reads2.fastq'
       self.test_fn2_space = self.temp_dir_spaces + 'reads2.fastq'
     except OSError:
            pass 
         
    def writeTmpFastq(self, fw_reads_path, rev_reads_path):
        """write forward and reverse reads data to temp fastq files"""
        fq1 = open(fw_reads_path, "w+")
        fq1.write(reads1_string)
        fq1.close()
        fq2 = open(rev_reads_path, "w+")
        fq2.write(reads2_string)
        fq2.close()


class SeqPrepTests(GenericSeqPrep):
    """Tests for SeqPrep application controller """

    def test_changing_working_dir(self):
        """WorkingDir should change properly.""" 
        c = SeqPrep(WorkingDir=self.temp_dir)
        self.assertEqual(c.BaseCommand,\
        ''.join(['cd "', self.temp_dir, '/"; ', 'SeqPrep']))

        c = SeqPrep()
        c.WorkingDir = self.temp_dir + '2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "', self.temp_dir + '2', '/"; ', 'SeqPrep']))

    def test_base_command(self):
        """seqprep command should return correct BaseCommand"""
        c = SeqPrep()
        # test base command
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "', getcwd(), '/"; ', 'SeqPrep']))
        # test turning on parameter
        c.Parameters['-O'].on('15')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "', getcwd(), '/"; ', 'SeqPrep -O 15']))


    def test_run_default_seqprep(self):
        """run_default_seqprep: should work as expected."""
        self.writeTmpFastq(self.test_fn1, self.test_fn2)
        
        # run with default function params
        res = run_default_seqprep(self.test_fn1, self.test_fn2,\
            self.temp_dir)






