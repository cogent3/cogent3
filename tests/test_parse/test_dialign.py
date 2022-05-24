#!/usr/bin/env python

import unittest

from cogent3 import PROTEIN
from cogent3.parse.dialign import DialignParser, parse_data_line


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


data = """
                           DIALIGN 2.2.1 
                           *************

          Program code written by Burkhard Morgenstern and Said Abdeddaim 
             e-mail contact: dialign (at) gobics (dot) de 

          Published research assisted by DIALIGN 2 should cite:  

             Burkhard Morgenstern (1999).
             DIALIGN 2: improvement of the segment-to-segment
             approach to multiple sequence alignment.
             Bioinformatics 15, 211 - 218. 

          For more information, please visit the DIALIGN home page at 

             http://bibiserv.techfak.uni-bielefeld.de/dialign/ 

         ************************************************************



   program call:  dialign -fa -fn /tmp/ct/seq1.fasta /tmp/ct/seq1.txt  


   Aligned sequences:          length:
   ==================          =======

   1) HTL2                        57
   2) MMLV                        58
   3) HEPB                        62
   4) ECOL                        54

   Average seq. length:           57.8 


   Please note that only upper-case letters are considered to be aligned. 


   Alignment (DIALIGN format):
   ===========================

HTL2               1   ldtapC-LFS DGS------P QKAAYVL--- ----WDQTIL QQDITPLPSH 
MMLV               1   pdadhtw-YT DGSSLLQEGQ RKAGAAVtte teviWa---- KALDAG---T 
HEPB               1   rpgl-CQVFA DAT------P TGWGLVM--- ----GHQRMR GTFSAPLPIH 
ECOL               1   mlkqv-EIFT DGSCLGNPGP GGYGAIL--- ----RYRGRE KTFSAGytrT 

                       0000000588 8882222229 9999999000 0000666666 6666633334 

HTL2              37   ethSAQKGEL LALICGLRAa k--------- --- 
MMLV              43   ---SAQRAEL IALTQALKm- ---------- --- 
HEPB              37   t------AEL LAA-CFARSr sganiigtdn svv 
ECOL              43   ---TNNRMEL MAAIv----- ---------- --- 

                       0003333455 5533333300 0000000000 000 




   Sequence tree:
   ==============

Tree constructed using UPGMAbased on DIALIGN fragment weight scores

((HTL2        :0.130254MMLV        :0.130254):0.067788(HEPB        :0.120520ECOL        :0.120520):0.077521);



""".splitlines()


class TestDialign(unittest.TestCase):
    def setUp(self):
        aln_seqs = {
            "HTL2": "ldtapC-LFSDGS------PQKAAYVL-------WDQTILQQDITPLPSHethSAQKGELLALICGLRAak------------",
            "MMLV": "pdadhtw-YTDGSSLLQEGQRKAGAAVtteteviWa----KALDAG---T---SAQRAELIALTQALKm--------------",
            "HEPB": "rpgl-CQVFADAT------PTGWGLVM-------GHQRMRGTFSAPLPIHt------AELLAA-CFARSrsganiigtdnsvv",
            "ECOL": "mlkqv-EIFTDGSCLGNPGPGGYGAIL-------RYRGREKTFSAGytrT---TNNRMELMAAIv------------------",
        }
        self.aln_seqs = {}
        for name, seq in list(aln_seqs.items()):
            self.aln_seqs[name] = PROTEIN.make_seq(seq, name=name)
        self.QualityScores = "00000005888882222229999999900000006666666666633334000333345555333333000000000000000"

    def test_line_split(self):
        """test splitting of sequence record lines"""
        result = parse_data_line(
            "HTL2               1   ldtapcLFSD GS------PQ KAAYVLWDQT ILQQDITPLP SHethsaqkg "
        )
        self.assertEqual(
            result, ("HTL2", "ldtapcLFSDGS------PQKAAYVLWDQTILQQDITPLPSHethsaqkg")
        )
        result = parse_data_line(
            "                       1111111111 1000001111 1111033333 3333333333 3000000000 "
        )
        self.assertEqual(
            result, (None, "11111111111000001111111103333333333333333000000000")
        )

    def test_aligned_from_dialign(self):
        """test getting aligned seqs"""
        aligned_seq = dict(list(DialignParser(data, seq_maker=PROTEIN.make_seq)))
        assert aligned_seq == self.aln_seqs

    def test_quality_scores(self):
        """test quality scores correctly returned"""
        result = dict(
            list(DialignParser(data, seq_maker=PROTEIN.make_seq, get_scores=True))
        )
        assert result["QualityScores"] == self.QualityScores


if __name__ == "__main__":
    unittest.main()
