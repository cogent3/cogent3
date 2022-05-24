#!/usr/bin/env python
import xml.dom.minidom

from unittest import TestCase, main

from cogent3.parse.tinyseq import TinyseqParser


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Matthew Wakefield"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production"

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
            self.assertEqual(str(seq.annotations), sample_annotations)

    pass


if __name__ == "__main__":
    main()
