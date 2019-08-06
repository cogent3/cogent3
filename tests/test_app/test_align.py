from unittest import TestCase, main

from cogent3 import DNA, LoadSeqs, LoadTree
from cogent3.app import align as align_app
from cogent3.app.composable import NotCompleted


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.08.06a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_seqs = {
    "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
    "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
    "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
    "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
}


class RefalignmentTests(TestCase):
    seqs = LoadSeqs(data=_seqs, aligned=False, moltype=DNA)
    treestring = "(Bandicoot:0.4,FlyingFox:0.05,(Rhesus:0.06," "Human:0.0):0.04);"

    def test_align_to_ref(self):
        """correctly aligns to a reference"""
        aligner = align_app.align_to_ref(ref_seq="Human")
        aln = aligner(self.seqs)
        expect = {
            "Bandicoot": "---NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
            "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT",
            "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
            "Rhesus": "GCCAGCTCATTACAGCATGAGAAC---AGTTTGTTACTCACT",
        }
        self.assertEqual(aln.todict(), expect)

    def test_progressive_align_nuc(self):
        """progressive alignment with nuc models"""
        aligner = align_app.progressive_align(model="F81")
        aln = aligner(self.seqs)
        expect = {
            "Rhesus": "GCCAGCTCATTACAGCATGAGAACAG---TTTGTTACTCACT",
            "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
            "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAG---TTTATTGTCCAAC",
            "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT",
        }
        got = aln.todict()
        self.assertEqual(got, expect)

        # using default
        aligner = align_app.progressive_align(model="nucleotide")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        self.assertEqual(aln.moltype, aligner._moltype)
        # todo the following is not robust across operating systems
        # so commenting out for now, but needs to be checked
        # expect = {'Human': 'GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT',
        #           'Rhesus': 'GCCAGCTCATTACAGCATGAGAA---CAGTTTGTTACTCACT',
        #           'Bandicoot': 'NACTCATTAATGCTTGAAACCAG---CAGTTTATTGTCCAAC',
        #           'FlyingFox': 'GCCAGCTCTTTACAGCATGAGAA---CAGTTTATTATACACT'}
        # got = aln.todict()
        # self.assertEqual(got, expect)

    def test_progressive_fails(self):
        """should return NotCompletedResult along with message"""
        # Bandicoot has an inf-frame stop codon
        seqs = LoadSeqs(
            data={"Human": "GCCTCA", "Rhesus": "GCCAGCTCA", "Bandicoot": "TGATCATTA"},
            aligned=False,
            moltype="dna",
        )
        aligner = align_app.progressive_align(model="codon")
        got = aligner(seqs)
        self.assertTrue(type(got), NotCompleted)

    def test_progress_with_guide_tree(self):
        """progressive align works with provided guide tree"""
        tree = LoadTree(treestring=self.treestring)
        aligner = align_app.progressive_align(
            model="nucleotide", guide_tree=self.treestring
        )
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        aligner = align_app.progressive_align(model="nucleotide", guide_tree=tree)
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)

    def test_progressive_align_codon(self):
        """progressive alignment with codon models"""
        aligner = align_app.progressive_align(model="GY94")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        aligner = align_app.progressive_align(model="codon")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)

    def test_with_genetic_code(self):
        """handles genetic code argument"""
        aligner = align_app.progressive_align(model="GY94", gc="2")
        # the 'TGA' codon is a sense codon in vertebrate mitochondrial
        self.assertTrue("TGA" in aligner._model.get_motifs())
        aligner = align_app.progressive_align(model="codon")
        # but a stop codon in the standard nuclear
        self.assertTrue("TGA" not in aligner._model.get_motifs())
        # try using a nuclear
        with self.assertRaises(TypeError):
            aligner = align_app.progressive_align(model="nucleotide", gc="2")

    def test_progressive_align_protein(self):
        """progressive alignment with protein models"""
        seqs = self.seqs.get_translation()
        with self.assertRaises(NotImplementedError):
            _ = align_app.progressive_align(model="protein")

        aligner = align_app.progressive_align(model="WG01", guide_tree=self.treestring)
        aln = aligner(seqs)
        self.assertEqual(len(aln), 14)
        aligner = align_app.progressive_align(
            model="protein", guide_tree=self.treestring
        )
        aln = aligner(seqs)
        self.assertEqual(len(aln), 14)


if __name__ == "__main__":
    main()
