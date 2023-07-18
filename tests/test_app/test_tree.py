import pathlib

from unittest import TestCase

import pytest

from cogent3 import DNA, load_aligned_seqs, make_aligned_seqs, make_tree
from cogent3.app import dist
from cogent3.app import tree as tree_app
from cogent3.app.composable import NotCompleted
from cogent3.core.tree import PhyloNode, TreeNode
from cogent3.evolve.fast_distance import DistanceMatrix


DATA_DIR = pathlib.Path(__file__).parent.parent / "data"


class TestTree(TestCase):
    def test_scale_tree_lengths(self):
        """correctly scales tree lengths"""
        with self.assertRaises(AssertionError):
            _ = tree_app.scale_branches(nuc_to_codon=True, codon_to_nuc=True)

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
        path = DATA_DIR / "brca1_5.paml"
        aln = load_aligned_seqs(path, moltype=DNA)
        fast_slow_dist = dist.fast_slow_dist(fast_calc="hamming", moltype="dna")
        dist_matrix = fast_slow_dist(aln)
        quick1 = tree_app.quick_tree()
        tree1 = quick1(dist_matrix)
        self.assertEqual(set(tree1.get_tip_names()), set(aln.names))

    def test_composable_apps(self):
        """checks the ability of these two apps(fast_slow_dist and quick_tree) to communicate"""
        path = DATA_DIR / "brca1_5.paml"
        aln1 = load_aligned_seqs(path, moltype=DNA)
        calc_dist = dist.fast_slow_dist(fast_calc="hamming", moltype="dna")
        quick = tree_app.quick_tree(drop_invalid=False)
        proc = calc_dist + quick
        self.assertEqual(
            str(proc),
            "fast_slow_dist(distance=None, moltype='dna', "
            "fast_calc='hamming',\nslow_calc=None) + quick_tree("
            "drop_invalid=False)",
        )
        self.assertIsInstance(proc, tree_app.quick_tree)
        self.assertIsInstance(proc.input, dist.fast_slow_dist)

        tree1 = proc(aln1)
        self.assertIsInstance(tree1, PhyloNode)
        self.assertIsNotNone(tree1.children)
        self.assertEqual(set(tree1.get_tip_names()), set(aln1.names))

        # tests when distances contain None
        data = dict(
            seq1="AGGGGGGGGGGCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGCGGTTTTTTTTTTTTTTTTTT",
        )
        aln2 = make_aligned_seqs(data=data, moltype=DNA)
        tree2 = proc(aln2)
        self.assertIsInstance(tree2, NotCompleted)

    def test_quick_tree_taking_distance_matrix(self):
        """quick_tree should take a distance matrix"""
        quick_tree = tree_app.quick_tree()
        data = {
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "Avin_42730"): 1.85,
        }

        darr = DistanceMatrix(data)
        tree = quick_tree(darr)
        self.assertIsInstance(tree, PhyloNode)
        self.assertIsNotNone(tree.children)
        self.assertEqual(
            set(tree.get_tip_names()), set.union(*(set(tup) for tup in data))
        )

        data = {
            ("DogFaced", "FlyingFox"): 0.05,
            ("DogFaced", "FreeTaile"): 0.14,
            ("DogFaced", "LittleBro"): 0.16,
            ("DogFaced", "TombBat"): 0.15,
            ("FlyingFox", "DogFaced"): 0.05,
            ("FlyingFox", "FreeTaile"): 0.12,
            ("FlyingFox", "LittleBro"): 0.13,
            ("FlyingFox", "TombBat"): 0.14,
            ("FreeTaile", "DogFaced"): 0.14,
            ("FreeTaile", "FlyingFox"): 0.12,
            ("FreeTaile", "LittleBro"): 0.09,
            ("FreeTaile", "TombBat"): 0.1,
            ("LittleBro", "DogFaced"): 0.16,
            ("LittleBro", "FlyingFox"): 0.13,
            ("LittleBro", "FreeTaile"): 0.09,
            ("LittleBro", "TombBat"): 0.12,
            ("TombBat", "DogFaced"): 0.15,
            ("TombBat", "FlyingFox"): 0.14,
            ("TombBat", "FreeTaile"): 0.1,
            ("TombBat", "LittleBro"): 0.12,
        }
        darr = DistanceMatrix(data)
        tree = quick_tree(darr)
        self.assertIsInstance(tree, PhyloNode)
        self.assertIsNotNone(tree.children)
        self.assertEqual(
            set(tree.get_tip_names()), set.union(*(set(tup) for tup in data))
        )

        data = {
            ("ABAYE2984", "Atu3667"): 0.25,
            ("ABAYE2984", "Avin_42730"): 0.638,
            ("ABAYE2984", "BAA10469"): None,
            ("Atu3667", "ABAYE2984"): 0.25,
            ("Atu3667", "Avin_42730"): 2.368,
            ("Atu3667", "BAA10469"): 0.25,
            ("Avin_42730", "ABAYE2984"): 0.638,
            ("Avin_42730", "Atu3667"): 2.368,
            ("Avin_42730", "BAA10469"): 1.85,
            ("BAA10469", "ABAYE2984"): 0.25,
            ("BAA10469", "Atu3667"): 0.25,
            ("BAA10469", "Avin_42730"): 1.85,
        }
        darr = DistanceMatrix(data)
        tree = quick_tree(darr)
        self.assertIsInstance(tree, PhyloNode)
        self.assertIsNotNone(tree.children)
        self.assertEqual(
            set(tree.get_tip_names()), set.union(*(set(tup) for tup in data))
        )

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

        darr = DistanceMatrix(data)
        # must explicitly call main() method to avoid error trapping by decorator
        with self.assertRaises(KeyError):
            quick_tree.main(darr)
        # when distance_matrix is None after dropping invalid
        with self.assertRaises(ValueError):
            quick_tree = tree_app.quick_tree(drop_invalid=True)
            quick_tree.main(darr)

        data = {
            ("DogFaced", "FlyingFox"): 0.05,
            ("DogFaced", "FreeTaile"): 0.14,
            ("DogFaced", "LittleBro"): 0.16,
            ("DogFaced", "TombBat"): 0.15,
            ("FlyingFox", "DogFaced"): 0.05,
            ("FlyingFox", "FreeTaile"): 0.12,
            ("FlyingFox", "LittleBro"): 0.13,
            ("FlyingFox", "TombBat"): 0.14,
            ("FreeTaile", "DogFaced"): 0.14,
            ("FreeTaile", "FlyingFox"): 0.12,
            ("FreeTaile", "LittleBro"): 0.09,
            ("FreeTaile", "TombBat"): 0.1,
            ("LittleBro", "DogFaced"): 0.16,
            ("LittleBro", "FlyingFox"): 0.13,
            ("LittleBro", "FreeTaile"): 0.09,
            ("LittleBro", "TombBat"): 0.12,
            ("TombBat", "DogFaced"): 0.15,
            ("TombBat", "FlyingFox"): 0.14,
            ("TombBat", "FreeTaile"): 0.1,
            ("TombBat", "LittleBro"): 0.12,
        }
        darr = DistanceMatrix(data)
        tree = quick_tree(darr)
        self.assertIsInstance(tree, PhyloNode)
        self.assertIsNotNone(tree.children)
        self.assertEqual(
            set(tree.get_tip_names()), set.union(*(set(tup) for tup in data))
        )

        data = {"a": {"b": 0.1, "a": 0.0}, "b": {"a": 0.1, "b": 0.0}}
        darr = DistanceMatrix(data)
        tree = quick_tree(darr)
        self.assertEqual(
            set(tree.get_tip_names()), set.union(*(set(tup) for tup in data))
        )

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


def test_interpret_tree_arg_none():
    assert tree_app.interpret_tree_arg(None) is None


@pytest.mark.parametrize(
    "tree",
    (
        DATA_DIR / "brca1_5.tree",
        str(DATA_DIR / "brca1_5.tree"),
        "(a,b,c)",
        make_tree(tip_names=["a", "b", "c"]),
    ),
)
def test_interpret_tree_arg_valid(tree):
    got = tree_app.interpret_tree_arg(tree)
    assert isinstance(got, TreeNode)


@pytest.mark.parametrize(
    "tree",
    (
        1,
        make_tree,
    ),
)
def test_interpret_tree_arg_invalid(tree):
    with pytest.raises(TypeError):
        tree_app.interpret_tree_arg(tree)
