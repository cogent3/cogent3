#!/usr/bin/env python
import xml.dom.minidom

from unittest import TestCase, main

from cogent3.parse.gbseq import GbSeqXmlParser


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"

data = """<?xml version="1.0"?>
 <!DOCTYPE GBSet PUBLIC "-//NCBI//NCBI GBSeq/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd">
 <GBSet>
<GBSeq>
  <GBSeq_locus>AY286018</GBSeq_locus>
  <GBSeq_length>99</GBSeq_length>
  <GBSeq_strandedness>single</GBSeq_strandedness>
  <GBSeq_moltype>mRNA</GBSeq_moltype>
  <GBSeq_topology>linear</GBSeq_topology>
  <GBSeq_division>MAM</GBSeq_division>
  <GBSeq_update-date>29-SEP-2003</GBSeq_update-date>
  <GBSeq_create-date>01-JUN-2003</GBSeq_create-date>
  <GBSeq_definition>Macropus eugenii medium wave-sensitive opsin 1 (OPN1MW) mRNA, complete cds</GBSeq_definition>
  <GBSeq_primary-accession>AY286018</GBSeq_primary-accession>
  <GBSeq_accession-version>AY286018.1</GBSeq_accession-version>
  <GBSeq_other-seqids>
    <GBSeqid>gb|AY286018.1|</GBSeqid>
    <GBSeqid>gi|31322957</GBSeqid>
  </GBSeq_other-seqids>
  <GBSeq_source>Macropus eugenii (tammar wallaby)</GBSeq_source>
  <GBSeq_organism>Macropus eugenii</GBSeq_organism>
  <GBSeq_taxonomy>Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Metatheria; Diprotodontia; Macropodidae; Macropus</GBSeq_taxonomy>
  <GBSeq_references>
    <GBReference>
      <GBReference_reference>1</GBReference_reference>
      <GBReference_position>1..99</GBReference_position>
      <GBReference_authors>
        <GBAuthor>Deeb,S.S.</GBAuthor>
        <GBAuthor>Wakefield,M.J.</GBAuthor>
        <GBAuthor>Tada,T.</GBAuthor>
        <GBAuthor>Marotte,L.</GBAuthor>
        <GBAuthor>Yokoyama,S.</GBAuthor>
        <GBAuthor>Marshall Graves,J.A.</GBAuthor>
      </GBReference_authors>
      <GBReference_title>The cone visual pigments of an Australian marsupial, the tammar wallaby (Macropus eugenii): sequence, spectral tuning, and evolution</GBReference_title>
      <GBReference_journal>Mol. Biol. Evol. 20 (10), 1642-1649 (2003)</GBReference_journal>
      <GBReference_xref>
        <GBXref>
          <GBXref_dbname>doi</GBXref_dbname>
          <GBXref_id>10.1093/molbev/msg181</GBXref_id>
        </GBXref>
      </GBReference_xref>
      <GBReference_pubmed>12885969</GBReference_pubmed>
    </GBReference>
    <GBReference>
      <GBReference_reference>2</GBReference_reference>
      <GBReference_position>1..99</GBReference_position>
      <GBReference_authors>
        <GBAuthor>Deeb,S.S.</GBAuthor>
        <GBAuthor>Wakefield,M.J.</GBAuthor>
        <GBAuthor>Tada,T.</GBAuthor>
        <GBAuthor>Marotte,L.</GBAuthor>
        <GBAuthor>Yokoyama,S.</GBAuthor>
        <GBAuthor>Graves,J.A.M.</GBAuthor>
      </GBReference_authors>
      <GBReference_title>Direct Submission</GBReference_title>
      <GBReference_journal>Submitted (29-APR-2003) RSBS, The Australian National University, Acton, ACT 0200, Australia</GBReference_journal>
    </GBReference>
  </GBSeq_references>
  <GBSeq_feature-table>
    <GBFeature>
      <GBFeature_key>source</GBFeature_key>
      <GBFeature_location>1..99</GBFeature_location>
      <GBFeature_intervals>
        <GBInterval>
          <GBInterval_from>1</GBInterval_from>
          <GBInterval_to>99</GBInterval_to>
          <GBInterval_accession>AY286018.1</GBInterval_accession>
        </GBInterval>
      </GBFeature_intervals>
      <GBFeature_quals>
        <GBQualifier>
          <GBQualifier_name>organism</GBQualifier_name>
          <GBQualifier_value>Macropus eugenii</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>mol_type</GBQualifier_name>
          <GBQualifier_value>mRNA</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>db_xref</GBQualifier_name>
          <GBQualifier_value>taxon:9315</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>country</GBQualifier_name>
          <GBQualifier_value>Australia: Kangaroo Island</GBQualifier_value>
        </GBQualifier>
      </GBFeature_quals>
    </GBFeature>
    <GBFeature>
      <GBFeature_key>gene</GBFeature_key>
      <GBFeature_location>1..99</GBFeature_location>
      <GBFeature_intervals>
        <GBInterval>
          <GBInterval_from>1</GBInterval_from>
          <GBInterval_to>99</GBInterval_to>
          <GBInterval_accession>AY286018.1</GBInterval_accession>
        </GBInterval>
      </GBFeature_intervals>
      <GBFeature_quals>
        <GBQualifier>
          <GBQualifier_name>gene</GBQualifier_name>
          <GBQualifier_value>OPN1MW</GBQualifier_value>
        </GBQualifier>
      </GBFeature_quals>
    </GBFeature>
    <GBFeature>
      <GBFeature_key>CDS</GBFeature_key>
      <GBFeature_location>31..99</GBFeature_location>
      <GBFeature_intervals>
        <GBInterval>
          <GBInterval_from>31</GBInterval_from>
          <GBInterval_to>99</GBInterval_to>
          <GBInterval_accession>AY286018.1</GBInterval_accession>
        </GBInterval>
      </GBFeature_intervals>
      <GBFeature_quals>
        <GBQualifier>
          <GBQualifier_name>gene</GBQualifier_name>
          <GBQualifier_value>OPN1MW</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>note</GBQualifier_name>
          <GBQualifier_value>cone pigments</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>codon_start</GBQualifier_name>
          <GBQualifier_value>1</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>transl_table</GBQualifier_name>
          <GBQualifier_value>1</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>product</GBQualifier_name>
          <GBQualifier_value>medium wave-sensitive opsin 1</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>protein_id</GBQualifier_name>
          <GBQualifier_value>AAP37945.1</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>db_xref</GBQualifier_name>
          <GBQualifier_value>GI:31322958</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>translation</GBQualifier_name>
          <GBQualifier_value>MTQAWDPAGFLAWRRDENE</GBQualifier_value>
        </GBQualifier>
      </GBFeature_quals>
    </GBFeature>
  </GBSeq_feature-table>
  <GBSeq_sequence>ggcagggaaagggaagaaagtaaaggggccatgacacaggcatgggaccctgcagggttcttggcttggcggcgggacgagaacgaggagacgactcgg</GBSeq_sequence>
</GBSeq>

</GBSet>
"""

sample_seq = ">AY286018.1\nGGCAGGGAAAGGGAAGAAAGTAAAGGGGCCATGACACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGG\n"
sample_annotations = '[source "[0:99]/99 of AY286018.1" at [0:99]/99, organism "Macropus eugenii" at [0:99]/99, gene "OPN1MW" at [0:99]/99, CDS "OPN1MW" at [30:99]/99]'


class ParseGBseq(TestCase):
    def test_parse(self):
        for name, seq in [
            next(GbSeqXmlParser(data)),
            next(GbSeqXmlParser(xml.dom.minidom.parseString(data))),
        ]:
            self.assertEqual(name, "AY286018.1")
            self.assertEqual(sample_seq, seq.to_fasta(block_size=len(sample_seq)))
            self.assertEqual(str(seq.annotations), sample_annotations)

    pass


if __name__ == "__main__":
    main()
