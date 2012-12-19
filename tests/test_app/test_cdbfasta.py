#!/usr/bin/env python

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.app.cdbfasta import cdbfasta, index_fasta, index_fastq, cdbyank, \
        query_indexed_seqs
from cogent.util.unit_test import TestCase, main
from shutil import rmtree

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class GeneralTests(TestCase):
    def setUp(self):
        """cdbfasta/cdbyank setUp method for all tests"""
        self.temp_dir = tempfile.mkdtemp()
        self.seqs_fasta = path.join(self.temp_dir, 'seqs.fasta')
        self.seqs_fastq = path.join(self.temp_dir, 'seqs.fastq')

        f = open(self.seqs_fasta,'w')
        self.fasta_seqs = {'ABC 567':'AATTGGCC','XYZ':'TTGGTT','123':'CCCC'}
        f.write(''.join([">%s\n%s\n"%(k,v) for k,v in self.fasta_seqs.items()]))
        f.close()

        f = open(self.seqs_fastq,'w')
        self.fastq_seqs = {'ABC 567':'AATTGGCC', 'XYZ':'TTGGTT','123':'CCCC'}
        self.fastq_qual = {'ABC 567':'12345678', 'XYZ':'333333','123':'9876'}

        for k in self.fastq_seqs:
            f.write("@%s\n%s\n+\n%s\n" % (k, self.fastq_seqs[k], 
                                          self.fastq_qual[k]))
        f.close()

    def tearDown(self):
        rmtree(self.temp_dir)

class cdbfastaTests(GeneralTests):
    def test_index_fasta(self):
        """Index a fasta!"""
        res = index_fasta(self.seqs_fasta)
        self.assertEqual(res.rsplit('/')[-1], "seqs.fasta.cidx") 

    def test_index_fastq(self):
        """Index a fastq!!!"""
        res = index_fastq(self.seqs_fastq)
        self.assertEqual(res.rsplit('/')[-1], "seqs.fastq.cidx")

    def test_version(self):
        """make sure the expected version string is produced"""
        app = cdbfasta(SuppressStdout=False, params={'-v':True})
        res = app(self.seqs_fasta)
        self.assertEqual(res['StdOut'].read(),"cdbfasta version 0.99\n")

class cdbyankTests(GeneralTests):
    def test_version(self):
        """make sure the expected version string is produced"""
        app = cdbyank(SuppressStdout=False, params={'-v':True})
        res = app(self.seqs_fasta)
        self.assertEqual(res['StdOut'].read(),"cdbyank version 0.981\n")

    def test_query_indexed_seqs_fasta(self):
        """pull out some sequences"""
        idx = index_fasta(self.seqs_fasta)
        seqs = query_indexed_seqs(['ABC','123'], idx)
        exp = sorted([">ABC 567","AATTGGCC",">123","CCCC"])
        obs = sorted(seqs.splitlines())
        self.assertEqual(obs,exp)
    
    def test_query_indexed_seqs_fastq(self):
        """pull out some sequences"""
        idx = index_fastq(self.seqs_fastq)
        seqs = query_indexed_seqs(['ABC','123'], idx)
        exp = sorted(["@ABC 567","AATTGGCC","+","12345678",
                      "@123","CCCC","+","9876"])
        obs = sorted(seqs.splitlines())
        self.assertEqual(obs,exp)

if __name__ == '__main__':
    main()
