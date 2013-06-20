#! /usr/bin/env python
# file: test_flash.py

# Tests for the flash.py application controller.
# FLASh v1.2.6
# http://ccb.jhu.edu/software/FLASH/
# using test_muscle.py as a guide

from cogent.util.unit_test import TestCase, main
from cogent.app.flash import Flash, default_assemble
from subprocess import Popen, PIPE, STDOUT
from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil

__author__ = "Michael Robeson"
__copyright__ = "Copyright 2007-2013, The Cogent Project"
__credits__ = ["Michael Robeson"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Michael Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

class GenericFlash(TestCase):
    """ Tests of the flash application contoller """
    
    def setUp(self):
        """ Set up some objects / data for use by tests"""

        # Check if flash version is supported for this test
        accepted_version = (1,2,5)
        command = "flash --version"
        version_cmd = Popoen(command, shell=True, universal_newlines=True,\
               stdout=PIPE,stderr=STDOUT)
        stdout = version_cmd.stdout.read()
        print stdout
        version_string = stdout.strip().split('\n')[0].split('v')[1]
        print version_string
        try:
            version = tuple(map(int, version_string.split('.')))
            print version
            pass_test = version == accepted_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
            "Unsupported flash version. %s is required, but running %s." \
            %('.'.join(map(str, accepted_version)), version_string))

        # set up files for flash tests
        self.fastq1 = "" 
        
        self.fastq2 = ""

        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test for flash/'
        try:
            mkdir(self.temp_dir_spaces)
        except OSError
            pass
        try:
            # create fastq files'
            reads1 = 
            reads1.write()
            reads1.close()
            reads2 = 
            reads2.write()
            reads2.close()
        except OSError:
            pass

    def tearDown(self):
        """ Cleans up all files initially created """
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)


class FlashTests(GenericFlash):

    def test_base_command(self):
        """ flash base command should return correct BaseCommand """

if __name__ == '__main__':
    main()
