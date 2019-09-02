import os

from unittest import TestCase, main

from cogent3 import DNA, load_aligned_seqs, make_tree
from cogent3.app import dist
from cogent3.app import tree as tree_app
from cogent3.app.composable import NotCompleted
from cogent3.core.tree import PhyloNode


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2019.8.30a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Alpha"


class TestTree(TestCase):
    def test_scale_tree_lengths(self):
        """correctly scales tree lengths"""
        with self.assertRaises(AssertionError):
            _ = tree_app.scale_branches(nuc_to_codon=True, codon_to_nuc=True)

        scale_to_codon = tree_app.scale_branches(nuc_to_codon=True)
        tree = make_tree(treestring="(a:3,b:6,c:9)")
        scale_to_codon = tree_app.scale_branches(nuc_to_codon=True)
        d = scale_to_codon(tree)
        got = {e.name: e.length for e in d.get_edge_vector(include_root=False)}
        expect = {"a": 1.0, "b": 2.0, "c": 3.0}
        self.assertEqual(got, expect)

        scale_from_codon = tree_app.scale_branches(codon_to_nuc=True)
        d = scale_from_codon(d)
        got = {e.name: e.length for e in d.get_edge_vector(include_root=False)}
        expect = {"a": 3.0, "b": 6.0, "c": 9.0}
        self.assertEqual(got, expect)

        by_scalar = tree_app.scale_branches(scalar=0.5)
        d = by_scalar(tree)
        got = {e.name: e.length for e in d.get_edge_vector(include_root=False)}
        expect = {"a": 6.0, "b": 12.0, "c": 18.0}
        self.assertEqual(got, expect)

        # handle case where a length is not defined, setting to minimum
        min_length = tree_app.scale_branches(min_length=66)
        tree = make_tree(treestring="(a:3,b:6,c)")
        new = min_length(tree)
        got = {e.name: e.length for e in new.get_edge_vector(include_root=False)}
        expect = {"a": 3.0, "b": 6.0, "c": 66.0}
        self.assertEqual(got, expect)

    def test_quick_tree(self):
        """correctly calc a nj tree"""
        path = os.path.join(
            os.path.abspath(__file__).split("test_app")[0], "data/brca1_5.paml"
        )
        aln = load_aligned_seqs(path, moltype=DNA)
        fast_slow_dist = dist.fast_slow_dist()
        dist_matrix = fast_slow_dist(aln)
        quick1 = tree_app.quick_tree()
        tree1 = quick1.quick_tree(dist_matrix)
        self.assertEqual(set(tree1.get_tip_names()), set(aln.names))
        # tests when distances cannot be computed
        data = {
            ("ABAYE2984", "Atu3667"): None,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): None,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): None,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): None,
            ("BAA10469", "Atu3667"): None,
            ("BAA10469", "Avin_42730"): 1.85,
        }
        from cogent3.evolve.fast_distance import DistanceMatrix

        darr = DistanceMatrix(data)
        quick2 = tree_app.quick_tree(drop_invalid=True)
        with self.assertRaises(ValueError):
            tree2 = quick2.quick_tree(darr)

    def test_communicable_apps(self):
        """checks the ability of these two apps to communicate"""
        path = os.path.join(
            os.path.abspath(__file__).split("test_app")[0], "data/brca1_5.paml"
        )
        aln1 = load_aligned_seqs(path, moltype=DNA)
        fast_slow_dist = dist.fast_slow_dist()
        quick = tree_app.quick_tree(drop_invalid=False)
        proc = fast_slow_dist + quick
        tree1 = proc(aln1)
        self.assertIsInstance(tree1, PhyloNode().__class__)
        self.assertIsNotNone(tree1.children)
        self.assertEqual(set(tree1.get_tip_names()), set(aln1.names))
        # tests when distances contain None
        from cogent3 import make_aligned_seqs

        data = dict(
            seq1="AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
            seq2="TAAAAAAAAAAGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCC",
        )
        aln2 = make_aligned_seqs(data=data, moltype=DNA)
        tree2 = proc(aln2)
        self.assertIsInstance(tree2, NotCompleted)

    def test_uniformize_tree(self):
        """equivalent topologies should be the same"""
        a = make_tree(treestring="(a,(b,c),(d,e))")
        b = make_tree(treestring="(e,d,(a,(b,c)))")
        make_uniform = tree_app.uniformize_tree(
            root_at="c", ordered_names=list("abcde")
        )
        u_a = make_uniform(a).get_newick()
        u_b = make_uniform(b).get_newick()
        self.assertTrue(u_a == u_b)
        # but different ones different
        c = make_tree(treestring="(e,c,(a,(b,d)))")
        u_c = make_uniform(c).get_newick()
        self.assertFalse(u_a == u_c)


if __name__ == "__main__":
    main()
