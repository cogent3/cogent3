#!/usr/bin/env python
from StringIO import StringIO
import xml.dom.minidom

from cogent.util.unit_test import TestCase, main
from cogent.parse.tinyseq import TinyseqParser

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Matthew Wakefield"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
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

sample_seq = ">AY286018.1\nGGCAGGGAAAGGGAAGAAAGTAAAGGGGCCATGACACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGG"
sample_annotations = '[genbank_id "AY286018.1" at [0:99]/99, organism "Macropus eugenii" at [0:99]/99]'

class ParseTinyseq(TestCase):
    def test_parse(self):
        for name,seq in [TinyseqParser(data).next(),TinyseqParser(xml.dom.minidom.parseString(data)).next()]:
            self.assertEqual(name, 'AY286018.1')
            self.assertEqual(sample_seq, seq.toFasta())
            self.assertEqual(str(seq.annotations), sample_annotations)
    pass
    
if __name__ == "__main__":
    main()
