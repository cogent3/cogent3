import os

from unittest import TestCase, main

from cogent3 import DNA, PROTEIN, make_unaligned_seqs
from cogent3.app import align
from cogent3.app import dist as dist_app
from cogent3.app import io, sample


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley", "Stephen Ma"]
__license__ = "BSD-3"
__version__ = "2019.8.30a"
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

    def test_composable_apps(self):
        composable_apps = _get_all_composable_apps()
        fast_slow_dist = dist_app.fast_slow_dist()
        for app in composable_apps:
            # Compose two composable applications, there should not be exceptions.
            got = app + fast_slow_dist
            self.assertIsInstance(got, dist.fast_slow_dist)
            self.assertEqual(got._type, "distance")
            self.assertIs(got.input, app)
            self.assertIs(got.output, None)
            self.assertIsInstance(got._input_type, frozenset)
            self.assertIsInstance(got._output_type, frozenset)
            self.assertIs(got._in, app)
            self.assertIs(got._out, None)
            app.disconnect()
            fast_slow_dist.disconnect()

    def test_est_dist_pair(self):
        """tests the distance between seq pairs in aln"""

        aligner = align.align_to_ref()
        aln3 = aligner(self.seqs3)
        fast_slow_dist = dist_app.fast_slow_dist()
        got = fast_slow_dist._est_dist_pair(aln3)
        self.assertAlmostEqual(got, 0.161372, places=6)

        aligner = align.align_to_ref(ref_seq="Human")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = dist_app.fast_slow_dist()
        got = fast_slow_dist._est_dist_pair(aln3)
        self.assertAlmostEqual(got, 0.161224, places=6)

        aligner = align.align_to_ref(ref_seq="Mouse")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = dist_app.fast_slow_dist()
        got = fast_slow_dist._est_dist_pair(aln3)
        self.assertAlmostEqual(got, 0.161372, places=6)

        aligner = align.align_to_ref()
        aln4 = aligner(self.seqs4)
        fast_slow_dist = dist_app.fast_slow_dist()
        got = fast_slow_dist._est_dist_pair(aln4)
        self.assertAlmostEqual(got, 0.591988, places=6)

        treestring = "(Human, Bandicoot)"
        aligner = align.progressive_align(model="WG01", guide_tree=treestring)
        aln5 = aligner(self.seqs5)
        fast_slow_dist = dist_app.fast_slow_dist(
            distance="paralinear", moltype="protein"
        )
        with self.assertRaises(AttributeError):
            got = fast_slow_dist._est_dist_pair(aln5)
        fast_slow_dist = dist_app.fast_slow_dist(distance="hamming", moltype="protein")
        with self.assertRaises(AttributeError):
            got = fast_slow_dist._est_dist_pair(aln5)


if __name__ == "__main__":
    main()
