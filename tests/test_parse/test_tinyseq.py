import xml.dom.minidom

from unittest import TestCase

from cogent3.parse.tinyseq import TinyseqParser


data = """<?xml version="1.0"?>
 <!DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd">
 <TSeqSet>
<TSeq>
  <TSeq_seqtype value="nucleotide"/>
  <TSeq_gi>31322957</TSeq_gi>
  <TSeq_accver>AY286018.1</TSeq_accver>
  <TSeq_taxid>9315</TSeq_taxid>
  <TSeq_orgname>Macropus eugenii</TSeq_orgname>
  <TSeq_defline>Macropus eugenii medium wave-sensitive opsin 1 (OPN1MW) mRNA, complete cds</TSeq_defline>
  <TSeq_length>99</TSeq_length>
  <TSeq_sequence>GGCAGGGAAAGGGAAGAAAGTAAAGGGGCCATGACACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGG</TSeq_sequence>
</TSeq>

</TSeqSet>
"""

sample_seq = ">AY286018.1\nGGCAGGGAAAGGGAAGAAAGTAAAGGGGCCATGACACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGG\n"
sample_annotations = (
    '[genbank_id "AY286018.1" at [0:99]/99, organism "Macropus eugenii" at [0:99]/99]'
)


class ParseTinyseq(TestCase):
    def test_parse(self):
        for name, seq in [
            next(TinyseqParser(data)),
            next(TinyseqParser(xml.dom.minidom.parseString(data))),
        ]:
            self.assertEqual(name, "AY286018.1")
            self.assertEqual(sample_seq, seq.to_fasta(block_size=len(sample_seq)))
            self.assertEqual(seq.annotation_db.num_matches(), 2)

    pass
