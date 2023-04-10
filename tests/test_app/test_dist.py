from tempfile import TemporaryDirectory
from unittest import TestCase, main

from numpy.testing import assert_allclose

from cogent3 import DNA, PROTEIN, get_app, make_unaligned_seqs, open_data_store
from cogent3.app.composable import WRITER
from cogent3.evolve.fast_distance import HammingPair, TN93Pair


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Stephen Ma", "Nick Shahmaras"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
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
    out_dstore = open_data_store(":memory:", mode="w")
    return [
        get_app("align_to_ref"),
        get_app("progressive_align", model="GY94"),
        get_app("fixed_length", 100),
        get_app("sample.min_length", 100),
        get_app("write_seqs", out_dstore),
        get_app(
            "omit_bad_seqs",
        ),
        get_app(
            "omit_degenerates",
        ),
        get_app("take_codon_positions", 1),
        get_app(
            "take_named_seqs",
        ),
        get_app("trim_stop_codons", gc=1),
    ]


class FastSlowDistTests(TestCase):
    seqs1 = make_unaligned_seqs(_seqs1, moltype=DNA)
    seqs2 = make_unaligned_seqs(_seqs2, moltype=DNA)
    seqs3 = make_unaligned_seqs(_seqs3, moltype=DNA)
    seqs4 = make_unaligned_seqs(_seqs4, moltype=DNA)
    seqs5 = make_unaligned_seqs(_seqs5, moltype=PROTEIN)

    def test_init(self):
        """tests if fast_slow_dist can be initialised correctly"""

        fast_slow_dist = get_app("fast_slow_dist", fast_calc="hamming", moltype="dna")
        self.assertIsInstance(fast_slow_dist.fast_calc, HammingPair)
        self.assertIsNone(fast_slow_dist._sm)

        fast_slow_dist = get_app("fast_slow_dist", distance="TN93")
        self.assertIsInstance(fast_slow_dist.fast_calc, TN93Pair)
        self.assertEqual(fast_slow_dist._sm.name, "TN93")
        fast_slow_dist = get_app("fast_slow_dist", distance="GTR")
        self.assertEqual(fast_slow_dist._sm.name, "GTR")

        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        self.assertEqual(fast_slow_dist._sm.name, "TN93")
        self.assertIsNone(fast_slow_dist.fast_calc)

        with self.assertRaises(ValueError):
            fast_slow_dist = get_app(
                "fast_slow_dist", distance="TN93", fast_calc="TN93", slow_calc="TN93"
            )

        with self.assertRaises(ValueError):
            fast_slow_dist = get_app("fast_slow_dist", fast_calc="GTR")

        with self.assertRaises(ValueError):
            fast_slow_dist = get_app("fast_slow_dist", slow_calc="hamming")

    def test_compatible_parameters(self):
        """tests if the input parameters are compatible with fast_slow_dist initialisation"""
        for kwargs in (
            dict(fast_calc="hamming", moltype="dna"),
            dict(fast_calc="TN93"),
            dict(slow_calc="GTR"),
            dict(fast_calc="TN93"),
        ):
            _ = get_app("fast_slow_dist", **kwargs)

    def test_incompatible_parameters(self):
        """tests incompatible input parameters with fast_slow_dist initialisation"""
        for kwargs in (
            dict(fast_calc="hamming"),
            dict(slow_calc="paralinear"),
            dict(fast_calc="GTR"),
            dict(slow_calc="hamming", moltype="dna"),
        ):
            with self.assertRaises(ValueError):
                _ = get_app("fast_slow_dist", **kwargs)

    def test_composable_apps(self):
        """tests two composable apps"""
        composable_apps = _get_all_composable_apps()
        calc_dist = get_app("fast_slow_dist", fast_calc="hamming", moltype="dna")
        for app in composable_apps:
            if app.app_type is WRITER:
                # cannot have a WRITER before a GENERIC
                continue
            # Compose two composable applications, there should not be exceptions.
            got = app + calc_dist
            self.assertIsInstance(got, type(calc_dist))
            self.assertIs(got.input, app)
            self.assertIsInstance(got._data_types, frozenset)
            self.assertIsInstance(got._return_types, frozenset)
            self.assertIs(got.input, app)
            app.disconnect()
            calc_dist.disconnect()

    def test_est_dist_pair_slow(self):
        """tests the distance between seq pairs in aln"""
        aligner = get_app(
            "align_to_ref",
        )
        aln3 = aligner(self.seqs3)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Human")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        assert_allclose(got[("Human", "Mouse")], got[("Mouse", "Human")])
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Mouse")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Mouse", "Human")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Mouse", "Human")] >= 0)

        aligner = get_app(
            "align_to_ref",
        )
        aln3 = aligner(self.seqs4)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Human")
        aln3 = aligner(self.seqs4)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        aligner = get_app("align_to_ref", ref_seq="Opossum")
        aln3 = aligner(self.seqs4)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="GTR")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)
        fast_slow_dist = get_app("fast_slow_dist", slow_calc="TN93")
        got = fast_slow_dist(aln3).to_dict()
        self.assertTrue(got[("Human", "Opossum")] >= 0)

        # now as a process
        proc = (
            get_app(
                "align_to_ref",
            )
            + get_app("fast_slow_dist", fast_calc="hamming", moltype="dna")
        )
        got = proc(self.seqs1)
        self.assertEqual(got[("Human", "Rhesus")], 1)

        treestring = "(Human:0.2,Bandicoot:0.2)"
        aligner = get_app("progressive_align", model="WG01", guide_tree=treestring)
        _ = aligner(self.seqs5)

    def test_composes_with_write_tabular(self):
        """correctly links to tabular"""
        with TemporaryDirectory(dir=".") as dirname:
            out_dstore = open_data_store(dirname, suffix="tsv", mode="w")
            writer = get_app("write_tabular", out_dstore)
            dist_calc = get_app("fast_slow_dist", distance="hamming", moltype="protein")
            _ = dist_calc + writer

    def test_functions_as_composable(self):
        """works as a composable app"""
        from pathlib import Path

        loader = get_app("load_aligned", moltype="dna", format="paml")
        dist = get_app("fast_slow_dist", "hamming", moltype="dna")
        with TemporaryDirectory(dir=".") as dirname:
            dirname = Path(dirname)
            out_dstore = open_data_store(dirname, suffix="tsv", mode="w")
            writer = get_app("write_tabular", out_dstore)
            proc = loader + dist + writer
            _ = proc("data/brca1_5.paml")
            output = dirname / "brca1_5.tsv"
            self.assertTrue(output.exists())


if __name__ == "__main__":
    main()
