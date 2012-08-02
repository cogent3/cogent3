#!/usr/bin/env python
# test_mothur.py


from __future__ import with_statement
from cStringIO import StringIO
from os import remove, rmdir, listdir, path
from tempfile import mkdtemp, mkstemp, NamedTemporaryFile
from cogent.util.unit_test import TestCase, main
from cogent.app.mothur import (
    Mothur, mothur_from_file,
    MothurClassifySeqs, parse_mothur_assignments, mothur_classify_file, 
    )


__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


class MothurTests(TestCase):
    def setUp(self):
        self.small_fasta = (
            '>aaaaaa\nTAGGCTCTGATATAATAGCTCTC---------\n'
            '>cccccc\n------------TGACTACGCAT---------\n'
            '>bbbbbb\n----TATCGCTTCGACGATTCTCTGATAGAGA\n'
            )
        self.small_otus = (
            'unique\t3\taaaaaa\tcccccc\tbbbbbb\t\n'
            '0.62\t2\taaaaaa\tbbbbbb,cccccc\t\n'
            '0.67\t1\tbbbbbb,cccccc,aaaaaa\t\n'
            )
        self.small_otus_parsed = [
            (float('0'), [['aaaaaa'], ['cccccc'], ['bbbbbb']]),
            (float('0.62'), [['aaaaaa'], ['bbbbbb', 'cccccc']]),
            (float('0.67'), [['bbbbbb', 'cccccc', 'aaaaaa']]),
            ]
        self.complement_fasta = (
            '>a\n--AGGGGTAATAA--\n'
            '>b\n--TTATTACCCCT--\n'
            '>c\n-------AAAAAA--\n'
            )
        self.complement_otus = (
            'unique\t3\ta\tb\tc\t\n'
            '0.43\t2\tc,a\tb\t\n'
            '1.00\t1\tb,c,a\t\n'
            )

    def test_get_help(self):
        """Mothur.getHelp() should return help string"""
        expected_help = (
            'See manual, available on the MOTHUR wiki:\n'
            'http://schloss.micro.umass.edu/mothur/'
            )
        self.assertEqual(Mothur.getHelp(), expected_help)

    def test_compile_mothur_script(self):
        """Mothur._compile_mothur_script() should return valid Mothur script"""
        app = Mothur()
        app._input_filename = 'test.fasta'
        observed_script = app._compile_mothur_script()
        expected_script = (
            '#unique.seqs(fasta=test.fasta); '
            'dist.seqs(fasta=test.unique.fasta); '
            'read.dist(column=test.unique.dist, name=test.names); '
            'cluster(method=furthest)')
        self.assertEqual(observed_script, expected_script)

    def test_get_result_paths(self):
        """Mothur._get_result_paths() should guess correct output paths"""
        app = Mothur()
        app._input_filename = 'test.fasta'
        observed_paths = {
            'distance matrix': app._derive_dist_path(),
            'otu list': app._derive_list_path(),
            'rank abundance': app._derive_rank_abundance_path(),
            'species abundance': app._derive_species_abundance_path(),
            'unique names': app._derive_names_path(),
            'unique seqs': app._derive_unique_path(),
            }
        expected_paths = {
            'distance matrix': 'test.unique.dist',
            'otu list': 'test.unique.fn.list',
            'rank abundance': 'test.unique.fn.rabund',
            'species abundance': 'test.unique.fn.sabund',
            'unique names': 'test.names',
            'unique seqs': 'test.unique.fasta',
            }
        self.assertEqual(observed_paths, expected_paths)

    def test_working_directory(self):
        """Mothur.WorkingDir attribute should not be cast to FilePath object"""
        app = Mothur(WorkingDir='/tmp')
        self.assertEquals(str(app.WorkingDir), '/tmp')

    def test_call_with_multiline_string(self):
        """Mothur.__call__() should return correct otu's for input as single string"""
        app = Mothur()
        result = app(self.small_fasta)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        result.cleanUp()

    def test_call_with_lines(self):
        """Mothur.__call__() should return correct otu's for input as lines"""
        lines = self.small_fasta.split('\n')
        app = Mothur(InputHandler='_input_as_lines')
        result = app(lines)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        result.cleanUp()

    def test_call_with_path(self):
        """Mothur.__call__() should return correct otu's for input as path"""
        working_dir = mkdtemp()
        _, filename = mkstemp(dir=working_dir, suffix='.fasta')
        with open(filename, 'w') as f:
            f.write(self.small_fasta)
        app = Mothur(InputHandler='_input_as_path', WorkingDir=working_dir)
        result = app(filename)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        remove(filename)
        result.cleanUp()
        rmdir(working_dir)

    def test_call_with_working_dir(self):
        """Mothur.__call__() should return correct otu's when input dir is changed"""
        working_dir = mkdtemp()
        app = Mothur(WorkingDir=working_dir)
        result = app(self.small_fasta)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.small_otus)
        result.cleanUp()
        rmdir(working_dir)

    def test_call_with_complement(self):
        """Mothur.__call__() should return correct otu's for input sequences which are reverse complements"""
        app = Mothur()
        result = app(self.complement_fasta)
        observed_otus = result['otu list'].read()
        self.assertEquals(observed_otus, self.complement_otus)
        result.cleanUp()

    def test_mothur_from_file(self):
        """mothur_from_file() should return parsed otus"""
        f = StringIO(self.small_fasta)
        f.seek(0)
        parsed_otus = mothur_from_file(f)
        self.assertEquals(parsed_otus, self.small_otus_parsed)


class MothurClassifySeqsTests(TestCase):
    def setUp(self):
        self.ref_seqs_file = NamedTemporaryFile(suffix=".fasta")
        self.ref_seqs_file.write(ref_seqs)
        self.ref_seqs_file.seek(0)
        
        self.taxonomy_file = NamedTemporaryFile(suffix=".tax")
        self.taxonomy_file.write(taxonomy)
        self.taxonomy_file.seek(0)

        self.params = {
            "reference": self.ref_seqs_file.name,
            "taxonomy": self.taxonomy_file.name,
            }
    
    def test_compile_mothur_script(self):
        app = MothurClassifySeqs({
            "reference": "/var/data/ref seqs.fasta",
            "taxonomy": "/var/data/taxonomy.txt",
            })
        app._input_filename = 'test.fasta'
        observed_script = app._compile_mothur_script()
        expected_script = (
            '#classify.seqs(fasta=test.fasta, '
            'reference=/var/data/ref seqs.fasta, '
            'taxonomy=/var/data/taxonomy.txt)'
            )
        self.assertEqual(observed_script, expected_script)

    def test_call_with_multiline_string(self):
        """__call__() should return correct classifications for input as single string"""
        work_dir = mkdtemp()
        app = MothurClassifySeqs(self.params, WorkingDir=work_dir)
        result = app(query_seqs)

        obs_summary = result["summary"].read()
        self.assertEqual(obs_summary, expected_summary)

        for line in result["assignments"]:
            self.assertEqual(len(line.split("\t")), 2)

            query_id, taxa_str = line.split("\t")
            self.assertTrue(query_id.startswith("Seq"))
            self.assertEqual(taxa_str.count(";"), 6)

        result.cleanUp()
        rmdir(work_dir)

    def test_cutoff_iters_ksize(self):
        work_dir = mkdtemp()
        self.params["cutoff"] = 99
        self.params["iters"] = 352
        self.params["ksize"] = 5
        app = MothurClassifySeqs(self.params, WorkingDir=work_dir)
        result = app(query_seqs)

        # Only exact match will remain with 99% cutoff
        f = result["assignments"]
        self.assertTrue(next(f).startswith("SeqA	Bacteria(100)"))
        self.assertTrue(next(f).startswith("SeqB	unknown"))
        self.assertTrue(next(f).startswith("SeqC	unknown"))

        result.cleanUp()
        rmdir(work_dir)

    def test_parse_mothur_assignments(self):
        assignments = StringIO(
            "SeqA	Bacteria(100);Firmicutes(100);Clostridia(100);"
            "Clostridiales(100);Ruminococcaceae(100);unclassified;\n"
            "SeqB	Bacteria(80);Firmicutes(80);Bacilli(71);"
            "Lactobacillales(71);Streptococcaceae(58);Streptococcus(58);\n"
            "SeqC	Archaea(66);Euryarchaeota(66);Methanobacteria(66);"
            "Methanobacteriales(66);Methanobrevibacter(66);unclassified;\n\n"
            )
        expected = [
            ("SeqA", ["Bacteria", "Firmicutes", "Clostridia",
              "Clostridiales", "Ruminococcaceae"], 1.0),
            ("SeqB", ["Bacteria", "Firmicutes", "Bacilli",
              "Lactobacillales", "Streptococcaceae", "Streptococcus"], 0.58),
            ("SeqC", ["Archaea", "Euryarchaeota", "Methanobacteria",
              "Methanobacteriales", "Methanobrevibacter"], 0.66),
            ]
        observed = list(parse_mothur_assignments(assignments))
        self.assertEqual(observed, expected)

    def test_mothur_classify_file(self):
        assignments = mothur_classify_file(
            StringIO(query_seqs),
            self.ref_seqs_file.name,
            self.taxonomy_file.name)
        self.assertEqual(list(sorted(assignments)), ["SeqA", "SeqB", "SeqC"])

    def test_mothur_classify_filepath_with_dash(self):
        temp_dir = tempfile.mkdtemp()

        ref_fp = os.path.join(temp_dir, "my-refs.fasta")
        ref_file = open(ref_fp, "w")
        ref_file.write(ref_seqs)
        ref_file.close()

        try:
            assignments = mothur_classify_file(
                StringIO(query_seqs),
                ref_fp,
                self.taxonomy_file.name)
            seq_ids = assignments.keys()
            seq_ids.sort()
            self.assertEqual(seq_ids, ["SeqA", "SeqB", "SeqC"])
        finally:
            shutil.rmtree(temp_dir)

ref_seqs = """\
>100
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
>101
CTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTC
CTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTC
CTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTC
CTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTCCTGAAGTC
>102
ATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCC
ATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCC
ATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCC
ATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCCATTTGCCC
>103
TGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTA
TGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTA
TGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTA
TGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTATGTGACTA
>104
CAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATT
CAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATT
CAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATT
CAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATTCAGGTATT
"""

taxonomy = """\
100	Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;
101	Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Genus_1;
102	Bacteria;Firmicutes;Bacilli;Lactobacillales;
103	Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;
104	Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobrevibacter;smithii
"""

query_seqs = """\
>SeqA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAA
>SeqB
ATTTGCCCATTTGCCCTGTGACTATGTGACTACTGAAGTCATTTGCCC
>SeqC
CAGGTATTCAGGTATTCAGGTATTGGTTCCAACTGAAGTCATTTGCCC
"""

expected_summary = """\
taxlevel\t rankID\t taxon\t daughterlevels\t total\t
0\t0\tRoot\t2\t3\t
1\t0.1\tArchaea\t1\t1\t
2\t0.1.1\tEuryarchaeota\t1\t1\t
3\t0.1.1.1\tMethanobacteria\t1\t1\t
4\t0.1.1.1.1\tMethanobacteriales\t1\t1\t
5\t0.1.1.1.1.1\tMethanobrevibacter\t1\t1\t
6\t0.1.1.1.1.1.1\tunclassified\t0\t1\t
1\t0.2\tBacteria\t1\t2\t
2\t0.2.1\tFirmicutes\t2\t2\t
3\t0.2.1.1\tBacilli\t1\t1\t
4\t0.2.1.1.1\tLactobacillales\t1\t1\t
5\t0.2.1.1.1.1\tStreptococcaceae\t1\t1\t
6\t0.2.1.1.1.1.1\tStreptococcus\t0\t1\t
3\t0.2.1.2\tClostridia\t1\t1\t
4\t0.2.1.2.1\tClostridiales\t1\t1\t
5\t0.2.1.2.1.1\tRuminococcaceae\t1\t1\t
6\t0.2.1.2.1.1.2\tunclassified\t0\t1\t
"""


if __name__ == '__main__':
    main()
