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

class SeqPrepTests(GenericSeqPrep):
    """Tests for SeqPrep application controller """
