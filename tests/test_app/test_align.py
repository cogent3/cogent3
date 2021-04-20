from unittest import TestCase, main

from cogent3 import DNA, get_moltype, make_tree, make_unaligned_seqs
from cogent3.align.align import make_generic_scoring_dict
from cogent3.app import align as align_app
from cogent3.app.composable import NotCompleted


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2021, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2021.04.20a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_seqs = {
    "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
    "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
    "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
    "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
}

_nucleotide_models = [
    "JC69",
    "K80",
    "F81",
    "HKY85",
    "TN93",
    "GTR",
    "ssGN",
    "GN",
    "BH",
    "DT",
]

_codon_models = [
    "CNFGTR",
    "CNFHKY",
    "MG94HKY",
    "MG94GTR",
    "GY94",
    "H04G",
    "H04GK",
    "H04GGK",
    "GNC",
]


class RefalignmentTests(TestCase):
    seqs = make_unaligned_seqs(_seqs, moltype=DNA)
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
        self.assertEqual(aln.to_dict(), expect)

    def test_align_to_ref_generic_moltype(self):
        """tests when the moltype is generic"""
        test_moltypes = ["text", "rna", "protein", "protein_with_stop", "bytes", "ab"]
        for test_moltype in test_moltypes:
            aligner = align_app.align_to_ref(moltype=test_moltype)
            self.assertEqual(aligner._moltype.label, test_moltype)
            self.assertEqual(
                aligner._kwargs["S"],
                make_generic_scoring_dict(10, get_moltype(test_moltype)),
            )

    def test_align_to_ref_result_has_moltype(self):
        """aligned object has correct moltype"""
        aligner = align_app.align_to_ref(moltype="dna")
        got = aligner(self.seqs)
        self.assertEqual(got.moltype.label, "dna")

    def test_progressive_align_protein_moltype(self):
        """tests guide_tree is None and moltype is protein"""
        from cogent3 import load_aligned_seqs

        seqs = load_aligned_seqs("data/nexus_aa.nxs", moltype="protein")
        seqs = seqs.degap()
        seqs = seqs.take_seqs(["Rat", "Cow", "Human", "Mouse", "Whale"])
        aligner = align_app.progressive_align(model="WG01")
        got = aligner(seqs)
        self.assertNotIsInstance(got, NotCompleted)
        aligner = align_app.progressive_align(model="protein")
        got = aligner(seqs)
        self.assertNotIsInstance(got, NotCompleted)

    def test_progressive_align_nuc(self):
        """progressive alignment with nuc models"""
        aligner = align_app.progressive_align(model="TN93", distance="TN93")
        aln = aligner(self.seqs)
        expect = {
            "Rhesus": "GCCAGCTCATTACAGCATGAGAACAG---TTTGTTACTCACT",
            "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
            "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAG---TTTATTGTCCAAC",
            "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAG---TTTATTATACACT",
        }
        got = aln.to_dict()
        self.assertEqual(got, expect)

        # using default
        aligner = align_app.progressive_align(model="TN93", distance="TN93")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        self.assertEqual(aln.moltype, aligner._moltype)
        # todo the following is not robust across operating systems
        # so commenting out for now, but needs to be checked
        # expect = {'Human': 'GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT',
        #           'Rhesus': 'GCCAGCTCATTACAGCATGAGAA---CAGTTTGTTACTCACT',
        #           'Bandicoot': 'NACTCATTAATGCTTGAAACCAG---CAGTTTATTGTCCAAC',
        #           'FlyingFox': 'GCCAGCTCTTTACAGCATGAGAA---CAGTTTATTATACACT'}
        # got = aln.to_dict()
        # self.assertEqual(got, expect)

    def test_progressive_fails(self):
        """should return NotCompletedResult along with message"""
        # Bandicoot has an inf-frame stop codon
        seqs = make_unaligned_seqs(
            data={"Human": "GCCTCA", "Rhesus": "GCCAGCTCA", "Bandicoot": "TGATCATTA"},
            moltype="dna",
        )
        aligner = align_app.progressive_align(model="codon")
        got = aligner(seqs)
        self.assertTrue(type(got), NotCompleted)

    def test_progress_with_guide_tree(self):
        """progressive align works with provided guide tree"""
        tree = make_tree(treestring=self.treestring)
        aligner = align_app.progressive_align(
            model="nucleotide", guide_tree=self.treestring
        )
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        aligner = align_app.progressive_align(model="nucleotide", guide_tree=tree)
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        # even if it has underscores in name
        treestring = (
            "(Bandicoot:0.4,FlyingFox:0.05,(Rhesus_macaque:0.06," "Human:0.0):0.04);"
        )
        aligner = align_app.progressive_align(model="nucleotide", guide_tree=treestring)
        data = self.seqs.to_dict()
        data["Rhesus macaque"] = data.pop("Rhesus")
        seqs = make_unaligned_seqs(data)
        aln = aligner(seqs)
        self.assertEqual(len(aln), 42)
        # guide tree with no lengths raises value error
        with self.assertRaises(ValueError):
            _ = align_app.progressive_align(
                model="nucleotide",
                guide_tree="(Bandicoot,FlyingFox,(Rhesus_macaque,Human));",
            )

    def test_progressive_align_codon(self):
        """progressive alignment with codon models"""
        aligner = align_app.progressive_align(model="GY94")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)
        aligner = align_app.progressive_align(model="codon")
        aln = aligner(self.seqs)
        self.assertEqual(len(aln), 42)

    def test_pickle_progressive_align(self):
        """test progressive_align is picklable"""
        from pickle import dumps, loads

        aligner = align_app.progressive_align(model="codon")
        aln = aligner(self.seqs)
        got = loads(dumps(aln))
        self.assertTrue(got)

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
