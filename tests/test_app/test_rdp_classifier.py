#!/usr/bin/env python
"""Tests for the rdp_classifier_2.0.1 application controller"""

from os import getcwd
from shutil import rmtree
from cogent.app.rdp_classifier import RdpClassifier
from cogent.util.unit_test import TestCase, main

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.4.0.dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"

class RdpClassifierTests(TestCase):
    def test_default_java_vm_parameters(self):
        """RdpClassifier should store default arguments to Java VM."""
        a = RdpClassifier()
        self.assertContains(a.Parameters, '-Xmx')
        self.assertEqual(a.Parameters['-Xmx'].Value, '1000m')

    def test_assign_java_vm_parameters(self):
        """RdpCalssifier should pass alternate parameters to Java VM."""
        app = RdpClassifier()
        app.Parameters['-Xmx'].on('75M')
        exp = ''.join(['cd "', getcwd(), '/"; java -Xmx75M -jar "rdp_classifier-2.0.jar"'])
        self.assertEqual(app.BaseCommand, exp)

        app = RdpClassifier()
        app.Parameters['-jar'].on('/jar/jar/binks.jar')
        exp = ''.join(['cd "', getcwd(), '/"; java -Xmx1000m -jar "/jar/jar/binks.jar"'])
        self.assertEqual(app.BaseCommand, exp)

    def test_basecommand_property(self):
        """RdpClassifier BaseCommand property should use overridden _get_base_command method."""
        app = RdpClassifier()
        self.assertEqual(app.BaseCommand, app._get_base_command())

    def test_base_command(self):
        """RdpClassifier should return expected shell command."""
        app = RdpClassifier()
        exp = ''.join(['cd "', getcwd(), '/"; java -Xmx1000m -jar "rdp_classifier-2.0.jar"'])
        self.assertEqual(app.BaseCommand, exp)

    def test_change_working_dir(self):
        """RdpClassifier should run program in expected working directory."""
        test_dir = '/tmp/RdpTest'

        app = RdpClassifier(WorkingDir=test_dir)
        exp = ''.join(['cd "', test_dir, '/"; java -Xmx1000m -jar "rdp_classifier-2.0.jar"'])
        self.assertEqual(app.BaseCommand, exp)

        rmtree(test_dir)

    def test_sample_fasta(self):
        """RdpClassifier should classify its own sample data correctly"""
        test_dir = '/tmp/RdpTest'

        app = RdpClassifier(WorkingDir=test_dir)
        app.Parameters['-jar'].on('/usr/local/app/rdp_classifier/rdp_classifier-2.0.jar')

        results_file = app(rdp_sample_fasta)['StdOut']
        
        id_line = results_file.readline()
        self.failUnless(id_line.startswith('>X67228'))

        classification_line = results_file.readline().strip()
        def all_even_items(list):
            return [x for (pos, x) in enumerate(list) if (pos % 2 == 0)]
        obs = all_even_items(classification_line.split('; '))
        exp = ['Root', 'Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Rhizobiales', 'Rhizobiaceae', 'Rhizobium']
        self.assertEqual(obs, exp)

        rmtree(test_dir)


# Sample data copied from rdp_classifier-2.0, which is licensed under
# the GPL 2.0 and Copyright 2008 Michigan State University Board of
# Trustees

rdp_sample_fasta = """>X67228 Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium
aacgaacgctggcggcaggcttaacacatgcaagtcgaacgctccgcaaggagagtggcagacgggtgagtaacgcgtgggaatctacccaaccctgcggaatagctctgggaaactggaattaataccgcatacgccctacgggggaaagatttatcggggatggatgagcccgcgttggattagctagttggtggggtaaaggcctaccaaggcgacgatccatagctggtctgagaggatgatcagccacattgggactgagacacggcccaaa
"""
rdp_sample_classification = """>X67228 reverse=false
Root; 1.0; Bacteria; 1.0; Proteobacteria; 1.0; Alphaproteobacteria; 1.0; Rhizobiales; 1.0; Rhizobiaceae; 1.0; Rhizobium; 0.95; 
"""

if __name__ == '__main__':
    main()

