#!/usr/bin/env python
"""Tests for the rdp_classifier_2.0.1 application controller"""

from os import getcwd, environ, remove
from shutil import rmtree
from cogent.app.util import ApplicationNotFoundError, get_tmp_filename
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
        # adjust environment variables
        self.user_rdp_jar_path = environ.pop('RDP_JAR_PATH', None)
        environ['RDP_JAR_PATH'] = 'rdp_classifier-2.0.jar'

    def tearDown(self):
        # reconstitute user's environment variables
        if self.user_rdp_jar_path is None:
            environ.pop('RDP_JAR_PATH')
        else:
            environ['RDP_JAR_PATH'] = self.user_rdp_jar_path

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
        exp = ''.join(['cd "', getcwd(), '/"; java -Xmx75M -jar "rdp_classifier-2.0.jar"'])
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
        
    def test_get_jar_fp(self):
        """RdpClassifier should search for jar file in appropriate location"""
        # Test 1: jar not in current dir, $RDP_JAR_PATH not set
        environ.pop('RDP_JAR_PATH')
        self.assertRaises(ApplicationNotFoundError, RdpClassifier)

        # Test 2: Subclass RdpClassifier with different _command
        # attribute, app should look there

        # create a temporary jar file
        tmp_jar_fp = get_tmp_filename(
            prefix='RdpClassifierTest', suffix='.jar')
        file = open(tmp_jar_fp, 'w')
        file.write('some stuff')
        file.close()

        class TestRdpClassifier(RdpClassifier):
            _command = tmp_jar_fp
        app = TestRdpClassifier()
        self.assertEqual(app._get_jar_fp(), tmp_jar_fp)

        remove(tmp_jar_fp)

        # reset environment variable
        environ['RDP_JAR_PATH'] = 'not None'

    def test_change_working_dir(self):
        """RdpClassifier should run program in expected working directory."""
        test_dir = '/tmp/RdpTest'

        app = RdpClassifier(WorkingDir=test_dir)
        exp = ''.join(['cd "', test_dir, '/"; java -Xmx1000m -jar "rdp_classifier-2.0.jar"'])
        self.assertEqual(app.BaseCommand, exp)

        rmtree(test_dir)

    def test_sample_fasta(self):
        """RdpClassifier should classify its own sample data correctly"""
        # Reconstitute user's $RDP_JAR_PATH variable
        if self.user_rdp_jar_path is None:
            environ.pop('RDP_JAR_PATH')
        else:
            environ['RDP_JAR_PATH'] = self.user_rdp_jar_path

        test_dir = '/tmp/RdpTest'
        try:
            app = RdpClassifier(WorkingDir=test_dir)

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

        finally:
            # reset environment variable
            environ['RDP_JAR_PATH'] = 'not None'

    def test_custom_training_data(self):
        """RdpClassifier should use sample training data"""
        # Reconstitute user's $RDP_JAR_PATH variable
        if self.user_rdp_jar_path is None:
            environ.pop('RDP_JAR_PATH')
        else:
            environ['RDP_JAR_PATH'] = self.user_rdp_jar_path

        test_dir = '/tmp/RdpTest'
        try:
            app = RdpClassifier(WorkingDir=test_dir)
            app.Parameters['-training-data'].on(
                '/usr/local/app/rdp_classifier/mydata/rRNAClassifier.properties'
                )

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

        finally:
            # reset environment variable
            environ['RDP_JAR_PATH'] = 'not None'

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

