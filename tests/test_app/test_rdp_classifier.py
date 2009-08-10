#!/usr/bin/env python
"""Tests for the rdp_classifier_2.0.1 application controller"""

from os import getcwd, environ, remove
from shutil import rmtree
from cogent.app.util import ApplicationNotFoundError, ApplicationError,\
    get_tmp_filename
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
    def setUp(self):
        # fetch user's RDP_JAR_PATH
        if 'RDP_JAR_PATH' in environ:
            self.user_rdp_jar_path = environ['RDP_JAR_PATH']
        else:
            self.user_rdp_jar_path = 'rdp_classifier-2.0.jar'

    def test_default_java_vm_parameters(self):
        """RdpClassifier should store default arguments to Java VM."""
        a = RdpClassifier()
        self.assertContains(a.Parameters, '-Xmx')
        self.assertEqual(a.Parameters['-Xmx'].Value, '1000m')

    def test_parameters_list(self):
        a = RdpClassifier()
        parameters = a.Parameters.keys()
        parameters.sort()
        self.assertEqual(parameters, ['-Xmx', '-training-data'])

    def test_jvm_parameters_list(self):
        a = RdpClassifier()
        parameters = a.JvmParameters.keys()
        parameters.sort()
        self.assertEqual(parameters, ['-Xmx'])

    def test_positional_parameters_list(self):
        a = RdpClassifier()
        parameters = a.PositionalParameters.keys()
        parameters.sort()
        self.assertEqual(parameters, ['-training-data'])

    def test_default_positional_parameters(self):
        """RdpClassifier should store default positional arguments."""
        a = RdpClassifier()
        self.assertContains(a.PositionalParameters, '-training-data')
        self.assertEqual(a.PositionalParameters['-training-data'].Value, '')        

    def test_assign_jvm_parameters(self):
        """RdpCalssifier should pass alternate parameters to Java VM."""
        app = RdpClassifier()
        app.Parameters['-Xmx'].on('75M')
        exp = ''.join(['cd "', getcwd(), '/"; java -Xmx75M -jar "', self.user_rdp_jar_path, '"'])
        self.assertEqual(app.BaseCommand, exp)

    def test_basecommand_property(self):
        """RdpClassifier BaseCommand property should use overridden _get_base_command method."""
        app = RdpClassifier()
        self.assertEqual(app.BaseCommand, app._get_base_command())

    def test_base_command(self):
        """RdpClassifier should return expected shell command."""
        app = RdpClassifier()
        exp = ''.join(['cd "', getcwd(), '/"; java -Xmx1000m -jar "', self.user_rdp_jar_path, '"'])
        self.assertEqual(app.BaseCommand, exp)
        
    def test_change_working_dir(self):
        """RdpClassifier should run program in expected working directory."""
        test_dir = '/tmp/RdpTest'

        app = RdpClassifier(WorkingDir=test_dir)
        exp = ''.join(['cd "', test_dir, '/"; java -Xmx1000m -jar "', self.user_rdp_jar_path, '"'])
        self.assertEqual(app.BaseCommand, exp)

        rmtree(test_dir)

    def test_sample_fasta(self):
        """RdpClassifier should classify its own sample data correctly"""
        test_dir = '/tmp/RdpTest'
        app = RdpClassifier(WorkingDir=test_dir)

        results_file = app(rdp_sample_fasta)['StdOut']
        
        id_line = results_file.readline()
        self.failUnless(id_line.startswith('>X67228'))

        classification_line = results_file.readline().strip()
        obs = parse_rdp(classification_line)
        exp = ['Root', 'Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Rhizobiales', 'Rhizobiaceae', 'Rhizobium']
        self.assertEqual(obs, exp)

        rmtree(test_dir)

    def test_custom_training_data(self):
        """RdpClassifier should look for custom training data"""
        test_dir = '/tmp/RdpTest'
        app = RdpClassifier(WorkingDir=test_dir)

        nonexistent_training_data = get_tmp_filename(
            prefix='RdpTestCustomTrainingData', suffix='.properties')
        app.Parameters['-training-data'].on(nonexistent_training_data)
        self.assertRaises(ApplicationError, app, rdp_sample_fasta)
        
        rmtree(test_dir)


def parse_rdp(line):
    """Returns a list of assigned taxa from an RDP classification line
    """
    tokens = line.split('; ')
    # Keep even-numbered tokens
    return [t for (pos, t) in enumerate(tokens) if (pos % 2 == 0)]

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

