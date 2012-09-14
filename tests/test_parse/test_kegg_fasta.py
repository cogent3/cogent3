#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Production"

"""
Test code for kegg_fasta.py in cogent.parse.  
"""
from cogent.util.unit_test import TestCase, main

from cogent.parse.kegg_fasta import kegg_label_fields, parse_fasta


class ParseKeggFastaTests(TestCase):
    def test_kegg_label_fields(self):
        """kegg_label_fields should return fields from line"""
        # Format is species:gene_id [optional gene_name]; description.
        # Note that the '>' should already be stripped by the Fasta Parser
        test1 = \
          """stm:STM0001  thrL; thr operon leader peptide ; K08278 thr operon leader peptide"""
        test2 = \
          """stm:STM0002  thrA; bifunctional aspartokinase I/homeserine dehydrogenase I (EC:2.7.2.4 1.1.1.13); K00003 homoserine dehydrogenase [EC:1.1.1.3]; K00928 aspartate kinase [EC:2.7.2.4]"""
        obs = kegg_label_fields(test1)
        exp = ('stm:STM0001','stm','STM0001',\
              'thrL','thr operon leader peptide ; K08278 thr operon leader peptide')
        self.assertEqual(obs,exp)

        obs = kegg_label_fields(test2)
        exp = ('stm:STM0002', 'stm', 'STM0002', 'thrA', \
            'bifunctional aspartokinase I/homeserine dehydrogenase I (EC:2.7.2.4 1.1.1.13); K00003 homoserine dehydrogenase [EC:1.1.1.3]; K00928 aspartate kinase [EC:2.7.2.4]')
        
        self.assertEqual(obs,exp)

    def test_parse_fasta(self):
        """parse_fasta should parse KEGG FASTA lines"""
        obs = parse_fasta(TEST_KEGG_FASTA_LINES)
        exp = EXP_RESULT
        for i,entry in enumerate(obs):
            self.assertEqual(entry, exp[i])

TEST_KEGG_FASTA_LINES = \
  [">stm:STM0001  thrL; thr operon leader peptide; K08278 thr operon leader peptide",\
   "atgaaccgcatcagcaccaccaccattaccaccatcaccattaccacaggtaacggtgcgggctga",\
   ">stm:STM0002  thrA; bifunctional aspartokinase I/homeserine dehydrogenase I (EC:2.7.2.4 1.1.1.13); K12524 bifunctional aspartokinase/homoserine dehydrogenase 1 [EC:2.7.2.4 1.1.1.3]",\
   "atgcgagtgttgaagttcggcggtacatcagtggcaaatgcagaacgttttctgcgtgtt",\
   "gccgatattctggaaagcaatgccaggcaagggcaggtagcgaccgtactttccgccccc"]

EXP_RESULT = \
  ["\t".join(["stm:STM0001","stm","STM0001",\
              "thrL","thr operon leader peptide; K08278 thr operon leader peptide","atgaaccgcatcagcaccaccaccattaccaccatcaccattaccacaggtaacggtgcgggctga","\n"]),\
   "\t".join(["stm:STM0002","stm","STM0002",\
              "thrA","bifunctional aspartokinase I/homeserine dehydrogenase I (EC:2.7.2.4 1.1.1.13); K12524 bifunctional aspartokinase/homoserine dehydrogenase 1 [EC:2.7.2.4 1.1.1.3]",\
              "atgcgagtgttgaagttcggcggtacatcagtggcaaatgcagaacgttttctgcgtgttgccgatattctggaaagcaatgccaggcaagggcaggtagcgaccgtactttccgccccc","\n"])]

if __name__=="__main__":
    main()
