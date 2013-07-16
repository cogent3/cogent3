#! /usr/bin/env python
# file: test_pandaseq.py

# Tests for the pandaseq.py application controller.
# https://github.com/neufeld/pandaseq

from cogent.util.unit_test import TestCase, main
from cogent.app.seqprep import SeqPrep, run_default_seqprep
from os import getcwd, remove, rmdir, mkdir, path, system
import gzip

__author__ = "Michael Robeson"
__copyright__ = "Copyright 2007-2013, The Cogent Project"
__credits__ = ["Michael Robeson"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Michael Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

class GenericPandaSeq(TestCase):
    def setUp(self):
        """General setup for PandaSeq tests """
        # make directory test
        self.temp_dir = '/tmp/test_for_pandaseq'
        try:
            mkdir(self.temp_dir)
        except OSError:
            pass

        # make directory with spaces test
        self.temp_dir_spaces = '/tmp/test for pandaseq/'
        try:
            mkdir(self.temp_dir_spaces)
        except OSError:
            pass

        # temp file paths
        self.test_fn1 = path.join(self.temp_dir,'reads1.fastq')
        self.test_fn1_space = path.join(self.temp_dir, 'reads1.fastq')
        self.test_fn2 = path.join(self.temp_dir,'reads2.fastq')
        self.test_fn2_space = path.join(self.temp_dir_spaces + 'reads2.fastq')

    def writeTmpFastq(self, fw_reads_path, rev_reads_path):
        """write forward and reverse reads data to temp fastq files"""
        try:
            fq1 = open(fw_reads_path, "w+")
            fq1.write(reads1_string)
            fq1.close()
            fq2 = open(rev_reads_path, "w+")
            fq2.write(reads2_string)
            fq2.close()
        except OSError:
            pass


class PandaSeqTests(GenericPandaSeq):
    """Tests for PandaSeq application controller"""

