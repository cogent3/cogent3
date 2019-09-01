from unittest import TestCase, main

from cogent3 import DNA, make_unaligned_seqs
from cogent3.app import align
from cogent3.app import dist as dist_app
from cogent3.app import tree as tree_app
from cogent3.core.tree import PhyloNode


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


class FastSlowDistTests(TestCase):
    seqs1 = make_unaligned_seqs(_seqs1, moltype=DNA)
    seqs2 = make_unaligned_seqs(_seqs2, moltype=DNA)
    seqs3 = make_unaligned_seqs(_seqs3, moltype=DNA)

    def test_composable_with_quick_tree(self):
        """fast_slow_dist should be composable"""
        fast_slow_dist = dist_app.fast_slow_dist()
        quick_tree = tree_app.quick_tree()
        got = fast_slow_dist + quick_tree

        self.assertEqual(
            str(got), "fast_slow_dist(type='distance') + quick_tree(type='tree')"
        )
        self.assertIsInstance(got, quick_tree.__class__)
        self.assertEqual(got._type, "tree")
        self.assertIsInstance(got.input, fast_slow_dist.__class__)
        self.assertIs(got.output, None)
        self.assertIsInstance(got._input_type, frozenset().__class__)
        self.assertIsInstance(got._output_type, frozenset().__class__)
        self.assertIsInstance(got._in, fast_slow_dist.__class__)
        self.assertIs(got._out, None)

    def test_quick_tree_taking_a_distance_matrix(self):
        """quick_tree should take a distance matrix"""
        quick_tree = tree_app.quick_tree()

        aligner = align.align_to_ref(ref_seq="Human")
        aln1 = aligner(self.seqs1)
        fast_slow_dist = dist_app.fast_slow_dist()
        obtained_dist_matrix1 = fast_slow_dist(aln1)
        self.assertIs(obtained_dist_matrix1.__class__.__name__, "DistanceMatrix")
        tree1 = quick_tree.quick_tree(obtained_dist_matrix1)
        self.assertIsInstance(tree1, PhyloNode().__class__)
        self.assertIsNotNone(tree1.children)
        self.assertEqual(set(tree1.get_tip_names()), set(aln1.names))

        aln2 = aligner(self.seqs2)
        obtained_dist_matrix2 = fast_slow_dist(aln2)
        self.assertIs(obtained_dist_matrix2.__class__.__name__, "DistanceMatrix")
        tree2 = quick_tree.quick_tree(obtained_dist_matrix2)
        self.assertIsInstance(tree2, PhyloNode().__class__)
        self.assertIsNotNone(tree2.children)
        self.assertEqual(set(tree2.get_tip_names()), set(aln2.names))

    def test_est_dist_pair(self):
        """tests the distance between seq pairs in aln"""
        aligner = align.align_to_ref(ref_seq="Human")
        aln3 = aligner(self.seqs3)
        fast_slow_dist = dist_app.fast_slow_dist()
        got = fast_slow_dist._est_dist_pair(aln3)
        self.assertAlmostEqual(got, 0.161224, places=6)


if __name__ == "__main__":
    main()
