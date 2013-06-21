#! /usr/bin/env python
# file: test_flash.py

# Tests for the flash.py application controller.
# FLASh v1.2.6
# http://ccb.jhu.edu/software/FLASH/
# using test_mafft.py as a guide

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
    
    def setUp(self):
        """Flash general setup fo all tests"""

        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test for flash/'
        try:
            mkdir(self.temp_dir_spaces)
        except OSError:
            pass
        try:
            # create fastq files'
            reads1 = open(path.join(self.temp_dir,'reads1.fastq'),'w')
            reads1.write(reads1_string) # bottom of file
            reads1.close()
            reads2 = open(path.join(self.temp_dir,'reads2.fastq'),'w')
            reads2.write(reads2_string) # bottom of file
            reads2.close()
        except OSError:
            pass


class FlashTests(GenericFlash):
    """Tests for FLASh application controller."""

    def check_version(self):
        """ Set up some objects / data for use by tests"""

        # Check if flash version is supported for this test
        accepted_version = (1,2,5)
        command = "flash --version"
        version_cmd = Popen(command, shell=True, universal_newlines=True,\
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


    def test_base_command(self):
        """ flash base command should return correct BaseCommand """

        c = Flash()
        # test base command
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "', getcwd(), '/"; ', 'flash']))
        # test turning on a parameter
        c.Parameters['-M'].on('500')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "', getcwd(), '/"; ', 'flash -M 500']))

        # test turning on another parameter via synonym
        # '--phred-offset' should use '-p'
        c.Parameters['--phred-offset'].on('33')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "', getcwd(), '/"; ', 'flash -M 500 -p 33']))

    def test_changing_working_dir(self):
        c = Flash(WorkingDir='/temp/flash_test')


reads1_string ="""@MISEQ05:114:000000000-A3UU7:1:1101:13650:2431 1:N:0:TCCACAGGAGT
AGGACGAACGCAGCGAAGGGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCTCTTTGGTATTCCGAAGAGCATGCCTGTTTGAGTGTCATGAAAATATCAACCTTGACTTGGGTTTAGTGCTCTTGTCTTGGCTTGGATTTGGCTGTTTGCCGCTCGAAAGAGTCGGCTCAGCTTAAAAGTATTAGCTGGATCTGTCTTTGAGACTTGGTTTGACTTGGCG
+
BCCBBCCCCCCCGGGGGGGGGGGGGGHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHFGHHGHHHHGHHHGGGGHHHHGGGGHHHFFF4FHHFBEEEFGC3B3BFGGHHHC3/F?FDFF3FDFFGDGFFDGHHHGHHHHFFHGGGHHHHHHHHHHHHHHHHHHHGGHHHHHHHHGHHHHHHGGGGGGHGGHGGGGGGFGGGGGGGBBFFFFFGFBFGGEFFGFFFFFFF00.:BFFBBB;F.09BF0FF#
@MISEQ05:114:000000000-A3UU7:1:1101:17667:2513 1:N:0:TCCACAGGAGT
AGAGGGAACGCAGCAAAGGGTGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCATCTTGCGCTCCTTGGTATTCTGAGGAGCATGCCTGTTTGAGTGTCATCAATCTCTCAACTACATCAATCTTTCTGGTTTGGTGTAGCTTTTGGATGTGGGGGTTTTATTTTGCTGGCCTCTATTAAATGAGGTCAGCTCCCCTGAAATTTATTAGCGGTATCTGAGCAGAGACCTACT
+
BCBBBCCFCCCCGGGGGGGGFFHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHGGGGHHHHGGGGGHHHGHHHHHHHHHGGGFEFFCFHHHHHG3EGFHHHHHHHHHHHHHHHFHGHHHHHHHHHHHHHHHGHHGGHGHHHHHHHHHGHHHHGHGGGGGGGHHHHHHHHHHGHHHHHHHHHGHHGHHHHHHGHHHGGGGGGGGGGGGGGEFGADGGGGGGGGGGF?FFFFFFF
@MISEQ05:114:000000000-A3UU7:1:1101:18777:2550 1:N:0:TCCACAGGAGT
GGCGGGAACGCAGCGAAGGGCGATAAGTAATGTGAATTGCAGAATTCAATGAATCATCGAATCTTTGAACGCACATTGCGCCCGTTGGTATTCCAACGGGCATGCCCGTTCGAGCGTCATTTCAATTCCTTCCCGGGAGACTCTTCTAATAAATCTTTTTTTTAGTGGGGTCCCGCGGGGCGTTGGCACTGGGCGCCACCTTTTTTGGAGGCGGCCGGAGCCGAAATAAAGCGGCGGACGCGCGGCGGCCC
+
CCBCCBBBBCCCGGGGGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHGHHHHHHHHGG@EHHFHGGGGGGHGHGHHHHHHHHGGGG?GCCGHGGHHDACFCGGGGHHHHHHHHGHHHHHDGGGGGGHHHHHHFHHHHHHHHHHHHGGGGEFGGECGGGGAGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDFFFFFAFFFFFFAFFE0FFFF.DAD;CFFD>DABFAB;@B""" 

reads2_string = """@MISEQ05:114:000000000-A3UU7:1:1101:13650:2431 2:N:0:TCCACAGGAGT
NGTCCTCCGCTTATTGATATGCTTAAGTTCAGCGGGTAGTCCTACCTGATTTGAGGCCAGATAATAAAAAAGTCATGTCTAGCAAGCTAGAACTTGTATTAGAAAGCGGACAAACGTCCCAGAAACTCGGCCAATCCGAAGATTGTCCTCAGCGAAATAACTTATTACGCCAAGTCAAACCAAGTCTCAAAGACAGATCAAGCTAATACTTTTAAGCTGAGCCGACTCTTTCGAGCGGCAAACAGCCAAAT
+
#>>>>ABCAAAA1FGGBGBFGFHHFFDGFF3B10E0A/FGDAFAAAF10BFG12/0///B0B11D22111//FFDFBGFHFG1110/F1BG1GHHHH1BG2212@FF?EG>/00F/B?FEBCEFBFHH/AFC?/FC</?<@EHHHHFHGFH1<<---<11FGFFGGFHG?CF-<<DFHHHGEEHHGHHHHFB/GHHGGG099FFGGGGGFGGGGGBFFFGAFF@G@;BFFEFAEFF?-@-@9-:-;;FFFB
@MISEQ05:114:000000000-A3UU7:1:1101:17667:2513 2:N:0:TCCACAGGAGT
CGTCCTCCGCTTATTGATATGCCTAAGTTCAGCGGGTAGTCCTACCTGATTTGAGGTCAAATTTGTCAAATGATTCACTGGAAGCAGCACAGTCTGTTGATGCAGCATTACCAAGGCATAGATAATTTATCACACCTGTAGTAGGTCTCTGCTCAGATACCGCTAATAAATTTCAGGGGAGCTGACCTCATTTAATAGAGGCCAGCAAAATAAAACCCCCACATCCAAAAGCTACACCAAACCAGAAAGAT
+
3>>AABFFBBBBGGGGGGGGGGHGHHHHHHFHBEEGGGGHHFH5ECGFCGHHG5FFHFGHHGHHHHFHFFGGHHHHGHHHHHEFEB3G1EGFHHHHHHFBG4FGHEHHHHHGG3?/CGFFFHHHHHHHHHHEFHH?HFGHHHHHHHHGHHHGFH3F2BDDCECGGGHHHHHFGGHFFH?GGGFHHHHHGGHHHHGGFBBFBDCEGGHHHHFGHFHF..:-DDDFGFFGFEFAFG0CB/;EEE..EFFBBF#
@MISEQ05:114:000000000-A3UU7:1:1101:18777:2550 2:N:0:TCCACAGGAGT
CGTCCTCCGCTTATTGATATGCTTAAGTTCAGCGGGTAGTCCTACCCGATTCGAGGTCAAGTTTTTATATAGCTTTATGCGTTGTAGGGCTGTGGTCCGCCGCTGCAGGTCTCCCAAGCGAGAAGGGTATTCTACTGCGCTAGGGGCCGCCGCGCGTCCGCCGCTTTATTTCGGCTCCGGCCGCCTCCAAAAAAGGTGGCGCCCAGTGCCACCGCCCCGCGGGACCCCACTAACAACAAGATTTATTAGAA
+
3>>AABFFBBBBAFGGEGGCGGGHHHHHHHGHGFGCECFHHFHGFAGDGGGGFEGGFHFGDFGHHGAGFGHHHHHHHHGHGEECFFB31FGGHHHHGH//EEEEGGHHGCFGGDG0CGG?/</?CCD@GFHFB1GHFFDGCGGHADGGGGGG-A-A?@BB@?-BBDEFBCFG/B?EDEF<BBBB@.9ABBFAABD./B..999@DE/:/BF/..-@>9--9-9-9ADF.A.9B##################
"""

if __name__ == '__main__':
    main()
