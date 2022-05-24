#!/usr/bin/env python
""" Unit tests for Genetic Code classes.
"""
from unittest import TestCase, main

from cogent3 import DNA, RNA
from cogent3.core.genetic_code import (
    DEFAULT,
    GeneticCode,
    GeneticCodeInitError,
    GeneticCodes,
    InvalidCodonError,
    available_codes,
    get_code,
)


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Greg Caporaso", "Rob Knight", "Peter Maxwell", "Thomas La"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Greg Caporaso"
__email__ = "caporaso@colorado.edu"
__status__ = "Production"


class GeneticCodeTests(TestCase):
    """Tests of the GeneticCode class."""

    def setUp(self):
        """Set up some standard genetic code representations."""
        self.SGC = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
        self.mt = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
        self.AllG = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"

        self.WrongLength = [
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" "",
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
        ]
        self.NcbiStandard = [
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            1,
            "Standard Nuclear",
            "---M---------------M---------------M----------------------------",
        ]

    def test_init(self):
        """GeneticCode init should work with correct-length sequences"""
        sgc = GeneticCode(self.SGC)
        self.assertEqual(sgc["UUU"], "F")
        mt = GeneticCode(self.mt)
        self.assertEqual(mt["UUU"], "F")
        allg = GeneticCode(self.AllG)
        self.assertEqual(allg["UUU"], "G")
        for i in self.WrongLength:
            self.assertRaises(GeneticCodeInitError, GeneticCode, i)

    def test_standard_code(self):
        """Standard genetic code from NCBI should have correct properties"""
        sgc = GeneticCode(*self.NcbiStandard)
        self.assertEqual(
            sgc.code_sequence,
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        )
        self.assertEqual(
            sgc.start_codon_sequence,
            "---M---------------M---------------M----------------------------",
        )
        self.assertEqual(sgc.start_codons, {"TTG": "M", "CTG": "M", "ATG": "M"})
        self.assertEqual(sgc.ID, 1)
        self.assertEqual(sgc.name, "Standard Nuclear")
        self.assertEqual(sgc["UUU"], "F")
        self.assertEqual(sgc.is_start("ATG"), True)
        self.assertEqual(sgc.is_start("AAA"), False)
        self.assertEqual(sgc.is_stop("UAA"), True)
        self.assertEqual(sgc.is_stop("AAA"), False)
        self.assertEqual(len(sgc.sense_codons), 61)
        self.assertIn("AAA", sgc.sense_codons)
        self.assertNotIn("TGA", sgc.sense_codons)

    def test_standard_code_lookup(self):
        """GeneticCodes should hold codes keyed by id as string and number"""
        sgc_new = GeneticCode(*self.NcbiStandard)
        sgc_number = GeneticCodes[1]
        sgc_string = GeneticCodes["1"]
        for sgc in sgc_new, sgc_number, sgc_string:
            self.assertEqual(
                sgc.code_sequence,
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            )
            self.assertEqual(
                sgc.start_codon_sequence,
                "---M---------------M---------------M----------------------------",
            )
            self.assertEqual(sgc.start_codons, {"TTG": "M", "CTG": "M", "ATG": "M"})
            self.assertEqual(sgc.ID, 1)
            self.assertEqual(sgc.name, "Standard Nuclear")
            self.assertEqual(sgc["TTT"], "F")
            self.assertEqual(sgc.is_start("ATG"), True)
            self.assertEqual(sgc.is_start("AAA"), False)
            self.assertEqual(sgc.is_stop("TAA"), True)
            self.assertEqual(sgc.is_stop("AAA"), False)

        mtgc = GeneticCodes[2]
        self.assertEqual(mtgc.name, "Vertebrate Mitochondrial")
        self.assertEqual(mtgc.is_start("AUU"), True)
        self.assertEqual(mtgc.is_stop("UGA"), False)

        self.assertEqual(
            sgc_new.changes(mtgc), {"AGA": "R*", "AGG": "R*", "ATA": "IM", "TGA": "*W"}
        )
        self.assertEqual(
            mtgc.changes(sgc_new), {"AGA": "*R", "AGG": "*R", "ATA": "MI", "TGA": "W*"}
        )
        self.assertEqual(mtgc.changes(mtgc), {})
        self.assertEqual(
            mtgc.changes(
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
            ),
            {"AGA": "*R", "AGG": "*R", "ATA": "MI", "TGA": "W*"},
        )

    def test_str(self):
        """GeneticCode str() should return its code string"""
        code_strings = self.SGC, self.mt, self.AllG
        codes = list(map(GeneticCode, code_strings))
        for code, string in zip(codes, code_strings):
            self.assertEqual(str(code), string)
        # check an example directly in case strings are bad
        self.assertEqual(
            str(self.SGC),
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        )

    def test_cmp(self):
        """GeneticCode cmp() should act on code strings"""
        sgc_1 = GeneticCode(self.SGC)
        sgc_2 = GeneticCode(self.SGC)
        self.assertEqual(sgc_1 is sgc_2, False)  # ensure different objects
        # self.assertNotEqual(sgc_1, sgc_2) # GREG
        self.assertEqual(sgc_1, sgc_2)
        mtgc = GeneticCode(self.mt)
        self.assertNotEqual(sgc_1, mtgc)

    def test_getitem_codon(self):
        """GeneticCode getitem should return amino acid for codon"""
        # specific checks of a particular codon in the standard code
        variant_codons = ["AUU", "AUU", "AUU", "ATT", "ATU", "ATU"]
        sgc = GeneticCode(self.SGC)
        for i in variant_codons:
            self.assertEqual(sgc[i], "I")
        # full check for the standard code
        codons = [a + b + c for a in "UCAG" for b in "TCAG" for c in "UCAG"]
        for codon, aa in zip(codons, self.SGC):
            self.assertEqual(sgc[codon], aa)
        # full check for another code
        allg = GeneticCode(self.AllG)
        for codon, aa in zip(codons, self.AllG):
            self.assertEqual(allg[codon], aa)
        # check that degenerate codon returns X
        self.assertEqual(sgc["NNN"], "X")

    def test_getitem_aa(self):
        """GeneticCode getitem should return codon set for aa"""
        # for all G, should return all the codons (in some order)
        allg = GeneticCode(self.AllG)
        codons = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
        g_codons = allg["G"]
        codons_copy = codons[:]
        self.assertEqual(g_codons, codons_copy)

        # check some known cases in the standard genetic code
        sgc = GeneticCode(self.SGC)
        exp_ile = ["ATT", "ATC", "ATA"]
        obs_ile = sgc["I"]
        self.assertEqual(obs_ile, exp_ile)

        sgc["R"]
        self.assertEqual(obs_ile, exp_ile)

        exp_leu = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
        obs_leu = sgc["L"]
        self.assertEqual(obs_leu, exp_leu)

        exp_met = ["ATG"]
        obs_met = sgc["M"]
        self.assertEqual(obs_met, exp_met)

        # unknown aa should return []
        self.assertEqual(sgc["U"], [])

    def test_getitem_invalid_length(self):
        """GeneticCode getitem should raise InvalidCodonError on wrong length"""
        sgc = GeneticCode(self.SGC)
        self.assertRaises(InvalidCodonError, sgc.__getitem__, "AAAA")
        self.assertRaises(InvalidCodonError, sgc.__getitem__, "AA")

    def test_Blocks(self):
        """GeneticCode Blocks should return correct list"""
        sgc = GeneticCode(self.SGC)
        exp_blocks = [
            ["TTT", "TTC"],
            ["TTA", "TTG"],
            ["TCT", "TCC", "TCA", "TCG"],
            ["TAT", "TAC"],
            ["TAA", "TAG"],
            ["TGT", "TGC"],
            ["TGA"],
            ["TGG"],
            ["CTT", "CTC", "CTA", "CTG"],
            ["CCT", "CCC", "CCA", "CCG"],
            ["CAT", "CAC"],
            ["CAA", "CAG"],
            ["CGT", "CGC", "CGA", "CGG"],
            ["ATT", "ATC"],
            ["ATA"],
            ["ATG"],
            ["ACT", "ACC", "ACA", "ACG"],
            ["AAT", "AAC"],
            ["AAA", "AAG"],
            ["AGT", "AGC"],
            ["AGA", "AGG"],
            ["GTT", "GTC", "GTA", "GTG"],
            ["GCT", "GCC", "GCA", "GCG"],
            ["GAT", "GAC"],
            ["GAA", "GAG"],
            ["GGT", "GGC", "GGA", "GGG"],
        ]
        self.assertEqual(sgc.blocks, exp_blocks)

    def test_Anticodons(self):
        """GeneticCode anticodons should return correct list"""
        sgc = GeneticCode(self.SGC)
        exp_anticodons = {
            "F": ["AAA", "GAA"],
            "L": ["TAA", "CAA", "AAG", "GAG", "TAG", "CAG"],
            "Y": ["ATA", "GTA"],
            "*": ["TTA", "CTA", "TCA"],
            "C": ["ACA", "GCA"],
            "W": ["CCA"],
            "S": ["AGA", "GGA", "TGA", "CGA", "ACT", "GCT"],
            "P": ["AGG", "GGG", "TGG", "CGG"],
            "H": ["ATG", "GTG"],
            "Q": ["TTG", "CTG"],
            "R": ["ACG", "GCG", "TCG", "CCG", "TCT", "CCT"],
            "I": ["AAT", "GAT", "TAT"],
            "M": ["CAT"],
            "T": ["AGT", "GGT", "TGT", "CGT"],
            "N": ["ATT", "GTT"],
            "K": ["TTT", "CTT"],
            "V": ["AAC", "GAC", "TAC", "CAC"],
            "A": ["AGC", "GGC", "TGC", "CGC"],
            "D": ["ATC", "GTC"],
            "E": ["TTC", "CTC"],
            "G": ["ACC", "GCC", "TCC", "CCC"],
        }
        self.assertEqual(sgc.anticodons, exp_anticodons)

    def test_translate(self):
        """GeneticCode translate should return correct amino acid string"""
        allg = GeneticCode(self.AllG)
        sgc = GeneticCode(self.SGC)
        mt = GeneticCode(self.mt)

        seq = "AUGCAUGACUUUUGA"
        #      .  .  .  .  .        markers for codon start
        self.assertEqual(allg.translate(seq), "GGGGG")
        self.assertEqual(allg.translate(seq, 1), "GGGG")
        self.assertEqual(allg.translate(seq, 2), "GGGG")
        self.assertEqual(allg.translate(seq, 3), "GGGG")
        self.assertEqual(allg.translate(seq, 4), "GGG")
        self.assertEqual(allg.translate(seq, 12), "G")
        self.assertEqual(allg.translate(seq, 14), "")
        self.assertRaises(ValueError, allg.translate, seq, 15)
        self.assertRaises(ValueError, allg.translate, seq, 20)

        self.assertEqual(sgc.translate(seq), "MHDF*")
        self.assertEqual(sgc.translate(seq, 3), "HDF*")
        self.assertEqual(sgc.translate(seq, 6), "DF*")
        self.assertEqual(sgc.translate(seq, 9), "F*")
        self.assertEqual(sgc.translate(seq, 12), "*")
        self.assertEqual(sgc.translate(seq, 14), "")
        # check shortest translatable sequences
        self.assertEqual(sgc.translate("AAA"), "K")
        self.assertEqual(sgc.translate(""), "")

        # check that different code gives different results
        self.assertEqual(mt.translate(seq), "MHDFW")

        # check translation with invalid codon(s)
        self.assertEqual(sgc.translate("AAANNNCNC123UUU"), "KXXXF")

    def test_sixframes(self):
        """GeneticCode sixframes should provide six-frame translation"""

        class fake_rna(str):
            """Fake RNA class with reverse-complement"""

            def __new__(cls, seq, rev):
                return str.__new__(cls, seq)

            def __init__(self, seq, rev):
                self.seq = seq
                self.rev = rev

            def rc(self):
                return self.rev

        test_rna = fake_rna("AUGCUAACAUAAA", "UUUAUGUUAGCAU")
        #                    .  .  .  .  .    .  .  .  .  .
        sgc = GeneticCode(self.SGC)
        self.assertEqual(
            sgc.sixframes(test_rna), ["MLT*", "C*HK", "ANI", "FMLA", "LC*H", "YVS"]
        )

        # should also actually work with an RNA or DNA sequence!!!
        test_rna = RNA.make_seq("AUGCUAACAUAAA")
        self.assertEqual(
            sgc.sixframes(test_rna), ["MLT*", "C*HK", "ANI", "FMLA", "LC*H", "YVS"]
        )

    def test_stop_indexes(self):
        """should return stop codon indexes for a specified frame"""
        sgc = GeneticCode(self.SGC)
        seq = DNA.make_seq("ATGCTAACATAAA")
        expected = [[9], [4], []]
        for frame, expect in enumerate(expected):
            got = sgc.get_stop_indices(seq, start=frame)
            self.assertEqual(got, expect)

    def test_Synonyms(self):
        """GeneticCode synonyms should return aa -> codon set mapping."""
        expected_synonyms = {
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "C": ["TGT", "TGC"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "F": ["TTT", "TTC"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "H": ["CAT", "CAC"],
            "I": ["ATT", "ATC", "ATA"],
            "K": ["AAA", "AAG"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "M": ["ATG"],
            "N": ["AAT", "AAC"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "Q": ["CAA", "CAG"],
            "R": ["AGA", "AGG", "CGT", "CGC", "CGA", "CGG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "W": ["TGG"],
            "Y": ["TAT", "TAC"],
            "*": ["TAA", "TAG", "TGA"],
        }
        obs_synonyms = GeneticCode(self.SGC).synonyms
        # note that the lists will be arbitrary-order
        for i in expected_synonyms:
            self.assertCountEqual(obs_synonyms[i], expected_synonyms[i])

    def test_get_code(self):
        """correctly return the genetic code"""
        for code_id in [1, "1", "Standard Nuclear", DEFAULT]:
            got = get_code(code_id)
            self.assertEqual(got, DEFAULT)
        got = get_code(2)
        self.assertEqual(got.name.lower(), "vertebrate mitochondrial")

    def test_available_codes(self):
        """avaialable_codes returns Table with correct shape"""
        codes = available_codes()
        self.assertEqual(codes.shape, (17, 2))
        self.assertEqual(codes[0, "Name"], "Standard Nuclear")

    def test_to_table(self):
        """tests to_table method."""
        from cogent3.core.moltype import IUPAC_PROTEIN_code_aa

        sgc = GeneticCode(self.SGC)
        table = sgc.to_table()
        # check num rows
        self.assertEqual(table.shape[0], len(IUPAC_PROTEIN_code_aa))
        # check there are 3 headers
        self.assertEqual(table.shape[1], 3)

    def test_to_regex(self):
        """creates a regex from aa seq to match a DNA sequence"""
        import re

        from cogent3 import make_seq

        dna = "ACCGAACAGGGC"
        aa = "TEQG"
        pattern = DEFAULT.to_regex(aa)
        self.assertTrue("".join(re.findall(pattern, dna)) == dna)
        # note that Z is Q or E
        aa = "TZQG"
        pattern = DEFAULT.to_regex(aa)
        self.assertTrue("".join(re.findall(pattern, dna)) == dna)
        aa = make_seq(aa, moltype="protein")
        pattern = DEFAULT.to_regex(aa)
        self.assertTrue("".join(re.findall(pattern, dna)) == dna)

    def test_repr_html(self):
        """exercising the _repr_html_ method"""
        gc = get_code(1)
        got = gc._repr_html_().strip()
        self.assertTrue(
            '<div class="c3table">' in got or '<div class="c3align">' in got
        )
        self.assertTrue("<table>" in got)
        self.assertTrue("</table>" in got)
        self.assertIn("Standard Nuclear", got)


# Run tests if called from command line
if __name__ == "__main__":
    main()
