import os

from tempfile import TemporaryDirectory
from unittest import TestCase, main

from numpy.testing import assert_allclose

from cogent3 import DNA, PROTEIN, make_unaligned_seqs
from cogent3.app import align
from cogent3.app import dist as dist_app
from cogent3.app import io, sample
from cogent3.evolve.fast_distance import HammingPair, TN93Pair


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Stephen Ma"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"

_seqs1 = {
    "Human": "GCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACT",
    "Bandicoot": "NACTCATTAATGCTTGAAACCAGCAGTTTATTGTCCAAC",
    "Rhesus": "GCCAGCTCATTACAGCATGAGAACAGTTTGTTACTCACT",
    "FlyingFox": "GCCAGCTCTTTACAGCATGAGAACAGTTTATTATACACT",
}

_seqs2 = {
    "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
    "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
    "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
}

_seqs3 = {
    "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
    "Mouse": "ATGCCCGGCGCCAAGGCAGCGCTGGCGGAG",
}

_seqs4 = {
    "Human": "ATGCGGCTCGCGGAGGCCGCGCTCGCGGAG",
    "Opossum": "ATGCCAGTGAAAGTGGCGGCGGTGGCTGAG",
}

_seqs5 = {"Human": "ASSLQHENSSLLLT", "Bandicoot": "XSLMLETSSLLSN"}


def _get_all_composable_apps():
    applications = [
        align.align_to_ref(),
        align.progressive_align(model="GY94"),
        sample.fixed_length(100),
        sample.min_length(100),
        io.write_seqs(os.getcwd()),
        sample.omit_bad_seqs(),
        sample.omit_degenerates(),
        sample.take_codon_positions(1),
        sample.take_named_seqs(),
        sample.trim_stop_codons(gc=1),
    ]
    return applications


class FastSlowDistTests(TestCase):
    seqs1 = make_unaligned_seqs(_seqs1, moltype=DNA)
    seqs2 = make_unaligned_seqs(_seqs2, moltype=DNA)
    seqs3 = make_unaligned_seqs(_seqs3, moltype=DNA)
    seqs4 = make_unaligned_seqs(_seqs4, moltype=DNA)
    seqs5 = make_unaligned_seqs(_seqs5, moltype=PROTEIN)

    def test_init(self):
        """tests if fast_slow_dist can be initialised correctly"""

        fast_slow_dist = dist_app.fast_slow_dist(fast_calc="hamming", moltype="dna")
        self.assertIsInstance(fast_slow_dist.fast_calc, HammingPair)
        self.assertIsNone(fast_slow_dist._sm)

        fast_slow_dist = dist_app.fast_slow_dist(distance="TN93")
        self.assertIsInstance(fast_slow_dist.fast_calc, TN93Pair)
        self.assertEqual(fast_slow_dist._sm.name, "TN93")
        fast_slow_dist = dist_app.fast_slow_dist(distance="GTR")
        self.assertEqual(fast_slow_dist._sm.name, "GTR")

        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="TN93")
        self.assertEqual(fast_slow_dist._sm.name, "TN93")
        self.assertIsNone(fast_slow_dist.fast_calc)

        with self.assertRaises(ValueError):
            fast_slow_dist = dist_app.fast_slow_dist(
                distance="TN93", fast_calc="TN93", slow_calc="TN93"
            )

        with self.assertRaises(ValueError):
            fast_slow_dist = dist_app.fast_slow_dist(fast_calc="GTR")

        with self.assertRaises(ValueError):
            fast_slow_dist = dist_app.fast_slow_dist(slow_calc="hamming")

    def test_compatible_parameters(self):
        """tests if the input parameters are compatible with fast_slow_dist initialisation"""
        fast_slow_dist = dist_app.fast_slow_dist(fast_calc="hamming", moltype="dna")
        fast_slow_dist = dist_app.fast_slow_dist(fast_calc="TN93")
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="GTR")
        fast_slow_dist = dist_app.fast_slow_dist(fast_calc="TN93")

        # fails for paralinear or hamming if no moltype
        with self.assertRaises(ValueError):
            fast_slow_dist = dist_app.fast_slow_dist(slow_calc="hamming")
        with self.assertRaises(ValueError):
            fast_slow_dist = dist_app.fast_slow_dist(slow_calc="paralinear")

        # fails for hamming as slow_calc
        with self.assertRaises(ValueError):
            fast_slow_dist = dist_app.fast_slow_dist(slow_calc="hamming", moltype="dna")
        with self.assertRaises(ValueError):
            fast_slow_dist = dist_app.fast_slow_dist(fast_calc="GTR")

    def test_composable_apps(self):
        """tests two composable apps"""
        composable_apps = _get_all_composable_apps()
        fast_slow_dist = dist_app.fast_slow_dist(fast_calc="hamming", moltype="dna")
        for app in composable_apps:
            # Compose two composable applications, there should not be exceptions.
            got = app + fast_slow_dist
            self.assertIsInstance(got, dist_app.fast_slow_dist)
            self.assertEqual(got._type, "distance")
            self.assertIs(got.input, app)
            self.assertIs(got.output, None)
            self.assertIsInstance(got._input_types, frozenset)
            self.assertIsInstance(got._output_types, frozenset)
            self.assertIs(got._in, app)
            self.assertIs(got._out, None)
            app.disconnect()
            fast_slow_dist.disconnect()

    def test_est_dist_pair_slow(self):
        """tests the distance between seq pairs in aln"""

        aligner = align.align_to_ref()
        aln3 = aligner(self.seqs3)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = align.align_to_ref(ref_seq="Human")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = align.align_to_ref(ref_seq="Mouse")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Mouse", "Human")] >= 0)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = align.align_to_ref()
        aln3 = aligner(self.seqs4)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        aligner = align.align_to_ref(ref_seq="Human")
        aln3 = aligner(self.seqs4)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        aligner = align.align_to_ref(ref_seq="Opossum")
        aln3 = aligner(self.seqs4)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = dist_app.fast_slow_dist(slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        # now as a process
        proc = align.align_to_ref() + dist_app.fast_slow_dist(
            fast_calc="hamming", moltype="dna"
        )
        got = proc(self.seqs1)
        self.assertEqual(got[("Human", "Rhesus")], 1)

        treestring = "(Human:0.2,Bandicoot:0.2)"
        aligner = align.progressive_align(model="WG01", guide_tree=treestring)
        _ = aligner(self.seqs5)

    def test_composes_with_write_tabular(self):
        """correctly links to tabular"""
        with TemporaryDirectory(dir=".") as dirname:
            writer = io.write_tabular(dirname)
            dist_calc = dist_app.fast_slow_dist(distance="hamming", moltype="protein")
            _ = dist_calc + writer

    def test_functions_as_composable(self):
        """works as a composable app"""
        from pathlib import Path

        loader = io.load_aligned(moltype="dna", format="paml")
        dist = dist_app.fast_slow_dist("hamming", moltype="dna")
        with TemporaryDirectory(dir=".") as dirname:
            dirname = Path(dirname)
            writer = io.write_tabular(dirname)
            proc = loader + dist + writer
            _ = proc("data/brca1_5.250.paml")
            output = dirname / "brca1_5.250.tsv"
            self.assertTrue(output.exists())


if __name__ == "__main__":
    main()
