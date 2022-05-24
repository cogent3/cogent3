#! /usr/bin/env python
import os
import pathlib
import unittest
import warnings

from tempfile import TemporaryDirectory

from numpy import exp, log

from cogent3 import get_model, load_aligned_seqs, load_tree, make_tree
from cogent3.phylo.consensus import get_splits, get_tree, majority_rule
from cogent3.phylo.least_squares import wls
from cogent3.phylo.maximum_likelihood import ML
from cogent3.phylo.nj import gnj, nj
from cogent3.phylo.tree_collection import (
    LogLikelihoodScoredTreeCollection,
    ScoredTreeCollection,
    WeightedTreeCollection,
    make_trees,
)
from cogent3.util.io import remove_files


warnings.filterwarnings("ignore", "Not using MPI as mpi4py not found")


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = [
    "Peter Maxwell",
    "Gavin Huttley",
    "Matthew Wakefield",
    "Daniel McDonald",
    "Ben Kaehler",
]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.path.dirname(__file__)
data_path = os.path.join(base_path, "data")


def Tree(t):
    return make_tree(treestring=t)


class ConsensusTests(unittest.TestCase):
    def setUp(self):
        self.trees = [
            Tree("((a,b),c,d);"),
            Tree("((a,b),c,d);"),
            Tree("((a,c),b,d);"),
            Tree("((a,b),c,d);"),
        ]

        # emphasizing the a,b clade
        weights = list(map(log, [0.4, 0.4, 0.05, 0.15]))
        self.scored_trees = list(zip(weights, self.trees))
        self.scored_trees.sort(reverse=True)

        self.rooted_trees = [
            Tree("((a,b),(c,d));"),
            Tree("((a,b),(c,d));"),
            Tree("((a,c),(b,d));"),
            Tree("((a,b),c,d);"),
        ]

        self.trees_randomly_rooted = [
            (
                -3416.3820971172017,
                Tree(
                    "(((F,(M,(P,(K,J)))),A),Q,(H,((S,(X,((Z,N),((O,(L,D)),E)))),((G,((I,(R,T)),Y)),(C,(((U,W),V),B))))));"
                ),
            ),
            (
                -3416.3820974814785,
                Tree(
                    "(R,T,(I,(Y,(G,((((X,((Z,N),((O,(L,D)),E))),S),((Q,((F,(M,(P,(K,J)))),A)),H)),((((U,W),V),C),B))))));"
                ),
            ),
            (
                -3416.571172739171,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,(S,((X,((Z,N),((O,(L,D)),E))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B)))))));"
                ),
            ),
            (
                -3416.5721393589474,
                Tree(
                    "(P,(K,J),(M,(F,(A,(Q,(H,(((((Z,N),((O,(L,D)),E)),S),X),((C,(G,((I,(R,T)),Y))),(((U,W),V),B)))))))));"
                ),
            ),
            (
                -3416.5721493054643,
                Tree(
                    "(R,T,(I,(Y,(G,(C,((S,(X,((Z,N),((O,(L,D)),E)))),((((U,W),V),B),((Q,((F,(M,(P,(K,J)))),A)),H))))))));"
                ),
            ),
            (
                -3416.572149424391,
                Tree(
                    "((U,W),V,(B,(((G,((I,(R,T)),Y)),C),((S,(X,((Z,N),((O,(L,D)),E)))),(H,(Q,((F,(M,(P,(K,J)))),A)))))));"
                ),
            ),
            (
                -3416.596424489769,
                Tree(
                    "(S,(X,((Z,N),((O,(L,D)),E))),((((G,((I,(R,T)),Y)),C),(((U,W),V),B)),((H,Q),((F,(M,(P,(K,J)))),A))));"
                ),
            ),
            (
                -3416.6335053806333,
                Tree(
                    "(M,(P,(K,J)),(F,(A,(Q,(H,(((((((O,(L,D)),E),N),Z),X),S),((C,((V,(U,W)),B)),(G,((I,(R,T)),Y)))))))));"
                ),
            ),
            (
                -3416.7626687401867,
                Tree(
                    "(Z,N,(((O,(L,D)),E),(X,(S,((((F,(M,(P,(K,J)))),A),Q),(H,((G,((I,(R,T)),Y)),(C,(((U,W),V),B)))))))));"
                ),
            ),
            (
                -3416.7670692165866,
                Tree(
                    "(X,((Z,N),((O,(L,D)),E)),(S,((B,(((G,((I,(R,T)),Y)),C),((U,W),V))),((Q,((F,(M,(P,(K,J)))),A)),H))));"
                ),
            ),
            (
                -3416.7670696377254,
                Tree(
                    "(P,(K,J),(M,(F,(A,(Q,(H,(((X,((Z,N),((O,(L,D)),E))),S),(((U,W),V),(((G,((I,(R,T)),Y)),C),B)))))))));"
                ),
            ),
            (
                -3416.848062302587,
                Tree(
                    "(((O,(L,D)),E),N,(Z,(X,(S,(((((F,(M,(P,(K,J)))),A),Q),H),(((V,(U,W)),B),(C,(G,((I,(R,T)),Y)))))))));"
                ),
            ),
            (
                -3416.9943503002764,
                Tree(
                    "((U,W),V,(B,((C,(G,((I,(R,T)),Y))),(H,(((((Z,N),((O,(L,D)),E)),S),X),(((F,(M,(P,(K,J)))),A),Q))))));"
                ),
            ),
            (
                -3417.014782481302,
                Tree(
                    "(Q,((F,(M,(P,(K,J)))),A),(((((U,W),V),B),((G,((I,(R,T)),Y)),C)),((S,(X,((Z,N),((O,(L,D)),E)))),H)));"
                ),
            ),
            (
                -3417.015470262783,
                Tree(
                    "(L,D,(O,(E,((Z,N),(X,(S,((Q,((F,(M,(P,(K,J)))),A)),(H,((((U,W),V),B),((G,((I,(R,T)),Y)),C))))))))));"
                ),
            ),
            (
                -3417.241619414339,
                Tree(
                    "(Z,N,(((O,(L,D)),E),(X,(S,(H,((Q,((F,(M,(P,(K,J)))),A)),(B,(((G,((I,(R,T)),Y)),C),((U,W),V)))))))));"
                ),
            ),
            (
                -3417.242009280534,
                Tree(
                    "(X,((Z,N),((O,(L,D)),E)),(S,((Q,((F,(M,(P,(K,J)))),A)),(H,(B,(((G,((I,(R,T)),Y)),C),((U,W),V)))))));"
                ),
            ),
            (
                -3417.2637092520818,
                Tree(
                    "(K,J,(P,(M,(F,(A,(Q,(((((((O,(L,D)),E),N),Z),X),S),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))))))));"
                ),
            ),
            (
                -3417.420526572887,
                Tree(
                    "(B,(V,(U,W)),(((H,(Q,(A,(F,(M,(P,(K,J))))))),((((((O,(L,D)),E),N),X),Z),S)),(C,(((I,(R,T)),Y),G))));"
                ),
            ),
            (
                -3417.4205266767162,
                Tree(
                    "(R,T,(I,(Y,(G,(C,((B,(V,(U,W))),((H,((A,(F,(M,(P,(K,J))))),Q)),(Z,(((((O,(L,D)),E),N),X),S)))))))));"
                ),
            ),
            (
                -3417.620921910812,
                Tree(
                    "((O,(L,D)),E,(N,(X,(Z,(S,(((Q,(A,(F,(M,(P,(K,J)))))),H),((((((I,(R,T)),Y),G),C),(V,(U,W))),B)))))));"
                ),
            ),
            (
                -3417.6209219461302,
                Tree(
                    "(F,(M,(P,(K,J))),(A,(Q,(H,((Z,(((((O,(L,D)),E),N),X),S)),((((G,((I,(R,T)),Y)),C),(V,(U,W))),B))))));"
                ),
            ),
            (
                -3417.6209224304744,
                Tree(
                    "(H,((A,(F,(M,(P,(K,J))))),Q),(B,((Z,(((((O,(L,D)),E),N),X),S)),(((G,((I,(R,T)),Y)),C),(V,(U,W))))));"
                ),
            ),
            (
                -3417.9379715010136,
                Tree(
                    "((I,(R,T)),Y,(G,(C,((B,(V,(U,W))),(H,((Q,(A,(F,(M,(P,(K,J)))))),((((((O,(L,D)),E),N),X),Z),S)))))));"
                ),
            ),
            (
                -3417.9379715187215,
                Tree(
                    "(((O,(L,D)),E),N,(X,(S,(Z,(((A,(F,(M,(P,(K,J))))),Q),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))))));"
                ),
            ),
        ]
        self.trees_rooted_at_A = [
            (
                -3416.3820971172017,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,((S,(X,((Z,N),((O,(L,D)),E)))),((G,((I,(R,T)),Y)),(C,(((U,W),V),B)))))));"
                ),
            ),
            (
                -3416.3820974814785,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,(((X,((Z,N),((O,(L,D)),E))),S),((G,((I,(R,T)),Y)),((((U,W),V),C),B))))));"
                ),
            ),
            (
                -3416.571172739171,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,(S,((X,((Z,N),((O,(L,D)),E))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B)))))));"
                ),
            ),
            (
                -3416.5721393589474,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,(((((Z,N),((O,(L,D)),E)),S),X),((C,(G,((I,(R,T)),Y))),(((U,W),V),B))))));"
                ),
            ),
            (
                -3416.5721493054643,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,((((U,W),V),B),(((G,((I,(R,T)),Y)),C),(S,(X,((Z,N),((O,(L,D)),E)))))))));"
                ),
            ),
            (
                -3416.572149424391,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,((S,(X,((Z,N),((O,(L,D)),E)))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B))))));"
                ),
            ),
            (
                -3416.596424489769,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(((S,(X,((Z,N),((O,(L,D)),E)))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B))),(H,Q)));"
                ),
            ),
            (
                -3416.6335053806333,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,(((((((O,(L,D)),E),N),Z),X),S),((C,((V,(U,W)),B)),(G,((I,(R,T)),Y)))))));"
                ),
            ),
            (
                -3416.7626687401867,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,((H,((G,((I,(R,T)),Y)),(C,(((U,W),V),B)))),(S,(X,((Z,N),((O,(L,D)),E)))))));"
                ),
            ),
            (
                -3416.7670692165866,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(((B,(((G,((I,(R,T)),Y)),C),((U,W),V))),((X,((Z,N),((O,(L,D)),E))),S)),H)));"
                ),
            ),
            (
                -3416.7670696377254,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,(((X,((Z,N),((O,(L,D)),E))),S),(((U,W),V),(((G,((I,(R,T)),Y)),C),B))))));"
                ),
            ),
            (
                -3416.848062302587,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(H,(((((((O,(L,D)),E),N),Z),X),S),(((V,(U,W)),B),(C,(G,((I,(R,T)),Y))))))));"
                ),
            ),
            (
                -3416.9943503002764,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(((((C,(G,((I,(R,T)),Y))),(((U,W),V),B)),H),((((Z,N),((O,(L,D)),E)),S),X)),Q));"
                ),
            ),
            (
                -3417.014782481302,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(((((U,W),V),B),((G,((I,(R,T)),Y)),C)),((S,(X,((Z,N),((O,(L,D)),E)))),H))));"
                ),
            ),
            (
                -3417.015470262783,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,((H,((((U,W),V),B),((G,((I,(R,T)),Y)),C))),(S,(X,((Z,N),((O,(L,D)),E)))))));"
                ),
            ),
            (
                -3417.241619414339,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,((H,((X,((Z,N),((O,(L,D)),E))),S)),(B,(((G,((I,(R,T)),Y)),C),((U,W),V))))));"
                ),
            ),
            (
                -3417.242009280534,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(((X,((Z,N),((O,(L,D)),E))),S),(H,(B,(((G,((I,(R,T)),Y)),C),((U,W),V)))))));"
                ),
            ),
            (
                -3417.2637092520818,
                Tree(
                    "((F,(M,(P,(K,J)))),A,(Q,(((((((O,(L,D)),E),N),Z),X),S),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))));"
                ),
            ),
            (
                -3417.420526572887,
                Tree(
                    "(A,(F,(M,(P,(K,J)))),(Q,(H,(((((((O,(L,D)),E),N),X),Z),S),((B,(V,(U,W))),(C,(((I,(R,T)),Y),G)))))));"
                ),
            ),
            (
                -3417.4205266767162,
                Tree(
                    "(A,(F,(M,(P,(K,J)))),(Q,(H,(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),(Z,(((((O,(L,D)),E),N),X),S))))));"
                ),
            ),
            (
                -3417.620921910812,
                Tree(
                    "(A,(F,(M,(P,(K,J)))),(Q,(H,(((((((O,(L,D)),E),N),X),Z),S),((((((I,(R,T)),Y),G),C),(V,(U,W))),B)))));"
                ),
            ),
            (
                -3417.6209219461302,
                Tree(
                    "(A,(F,(M,(P,(K,J)))),(Q,(H,((Z,(((((O,(L,D)),E),N),X),S)),((((G,((I,(R,T)),Y)),C),(V,(U,W))),B)))));"
                ),
            ),
            (
                -3417.6209224304744,
                Tree(
                    "(A,(F,(M,(P,(K,J)))),(Q,(H,(B,((Z,(((((O,(L,D)),E),N),X),S)),(((G,((I,(R,T)),Y)),C),(V,(U,W))))))));"
                ),
            ),
            (
                -3417.9379715010136,
                Tree(
                    "(A,(F,(M,(P,(K,J)))),(Q,((H,((B,(V,(U,W))),(C,(((I,(R,T)),Y),G)))),((((((O,(L,D)),E),N),X),Z),S))));"
                ),
            ),
            (
                -3417.9379715187215,
                Tree(
                    "(A,(F,(M,(P,(K,J)))),(Q,((Z,(((((O,(L,D)),E),N),X),S)),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))));"
                ),
            ),
        ]

        self.unrooted_conflicting_trees = [
            Tree("((a,b),c,d);"),
            Tree("((a,c),b,d);"),
            Tree("((a,d),b,c);"),
        ]

        self.rooted_conflicting_trees = [
            Tree("((a,b),(c,d));"),
            Tree("((a,c),(b,d));"),
            Tree("((a,d),(b,c));"),
        ]

        self.unrooted_trees_lengths = [
            (2, Tree("((a:0.3,c:0.4):0.5,b:0.2,d:0.1);")),
            (1, Tree("((a:0.1,b:0.1):0.1,c:0.1,d:0.1);")),
        ]

        self.rooted_trees_lengths = [
            (2, Tree("((a:0.3,c:0.4):0.2,(b:0.2,d:0.1):0.3);")),
            (1, Tree("((a:0.1,b:0.1):0.05,(c:0.1,d:0.1):0.05);")),
        ]

    def test_majorityRule(self):
        """Tests for majority rule consensus trees"""
        trees = self.rooted_trees
        outtrees = majority_rule(trees, strict=False)
        self.assertEqual(len(outtrees), 1)
        self.assertTrue(outtrees[0].same_topology(Tree("((c,d),(a,b));")))
        outtrees = majority_rule(trees, strict=True)
        self.assertEqual(len(outtrees), 1)
        self.assertTrue(outtrees[0].same_topology(Tree("(c,d,(a,b));")))

    def test_get_tree_get_splits(self):
        """get_tree should provide a reciprocal map of get_splits"""
        tree = load_tree(os.path.join(data_path, "murphy.tree"))
        self.assertTrue(tree.same_topology(get_tree(get_splits(tree))))

    def test_consensus_tree_branch_lengths(self):
        """consensus trees should average branch lengths properly"""

        def get_ac(tree):
            for edge in tree.get_edge_vector(include_root=False):
                if set("ac") == set([c.name for c in edge.children]):
                    return edge

        sct = ScoredTreeCollection(self.unrooted_trees_lengths)
        ct = sct.get_consensus_tree()
        maj_tree = self.unrooted_trees_lengths[0][1]
        # to ensure consistent comparison with majority, we root the ct same way
        # as maj
        tip_names = maj_tree.get_tip_names()
        ct = ct.rooted_with_tip("d")
        ct = ct.sorted(tip_names)

        self.assertTrue(abs(get_ac(ct).length - get_ac(maj_tree).length) < 1e-9)

        sct = ScoredTreeCollection(self.rooted_trees_lengths)
        ct = sct.get_consensus_tree(method="rooted")
        maj_tree = self.rooted_trees_lengths[0][1]
        self.assertTrue(abs(get_ac(ct).length - get_ac(maj_tree).length) < 1e-9)

    def test_scored_trees_collection_write(self):
        """writes a tree collection"""
        sct = ScoredTreeCollection(self.rooted_trees_lengths)
        with TemporaryDirectory(".") as dirname:
            dirname = pathlib.Path(dirname)
            out = dirname / "collection.trees"
            sct.write(out)

    def test_consensus_from_scored_trees_collection(self):
        """tree collection should get same consensus as direct approach"""
        tree_list = [(i * -1, t) for i, t in enumerate(self.trees)]
        sct = LogLikelihoodScoredTreeCollection(tree_list)
        ct = sct.get_consensus_tree()
        self.assertTrue(ct.same_topology(Tree("((c,d),a,b);")))

    def test_consensus_from_scored_trees_collection_ii(self):
        """strict consensus should handle conflicting trees"""
        sct = ScoredTreeCollection(list(zip([1] * 3, self.unrooted_conflicting_trees)))
        ct = sct.get_consensus_trees()[0]
        self.assertTrue(ct.same_topology(Tree("(a,b,c,d);")))

        sct = ScoredTreeCollection(list(zip([1] * 3, self.rooted_conflicting_trees)))
        # cts = sct.get_consensus_trees(method='rooted')
        ct = sct.get_consensus_trees(method="rooted")[0]
        self.assertTrue(ct.same_topology(Tree("(a,b,c,d);")))
        # for tree in cts:
        #    print str(tree)
        # self.assertTrue(set(map(str, cts))==set(['('+c+');' for c in 'abcd']))

    def test_weighted_consensus_from_scored_trees_collection(self):
        """weighted consensus from a tree collection should be different"""
        sct = LogLikelihoodScoredTreeCollection(self.scored_trees)
        ct = sct.get_consensus_tree()
        self.assertTrue(ct.same_topology(Tree("((a,b),c,d);")))

    def test_weighted_consensus_from_scored_trees_collection_ii(self):
        """root positions in input tree collection should not effect result"""
        sct = LogLikelihoodScoredTreeCollection(self.trees_randomly_rooted)
        ctrr = sct.get_consensus_tree()
        sct = LogLikelihoodScoredTreeCollection(self.trees_rooted_at_A)
        ctra = sct.get_consensus_tree()
        self.assertTrue(ctrr.same_topology(ctra))

    def test_weighted_trees_satisyfing_cutoff(self):
        """build consensus tree from those satisfying cutoff"""
        sct = LogLikelihoodScoredTreeCollection(self.scored_trees)
        cts = sct.get_weighted_trees(cutoff=0.8)
        for weight, tree in cts:
            self.assertTrue(tree.same_topology(Tree("((a,b),c,d);")))

        ct = cts.get_consensus_tree()
        self.assertTrue(ct.same_topology(Tree("((a,b),c,d);")))

    def test_tree_collection_read_write_file(self):
        """should correctly read / write a collection from a file"""

        def eval_klass(coll):
            coll.write("sample.trees")
            read = make_trees("sample.trees")
            self.assertTrue(type(read) == type(coll))

        eval_klass(LogLikelihoodScoredTreeCollection(self.scored_trees))

        # convert lnL into p
        eval_klass(WeightedTreeCollection([(exp(s), t) for s, t in self.scored_trees]))
        remove_files(["sample.trees"], error_on_missing=False)


class TreeReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.tree = make_tree(treestring="((a:3,b:4):2,(c:6,d:7):30,(e:5,f:5):5)")
        self.dists = self.tree.get_distances()

    def assertTreeDistancesEqual(self, t1, t2):
        d1 = t1.get_distances()
        d2 = t2.get_distances()
        self.assertEqual(len(d1), len(d2))
        for key in d2:
            self.assertAlmostEqual(d1[key], d2[key])

    def test_nj(self):
        """testing nj"""
        reconstructed = nj(self.dists, show_progress=False)
        self.assertTreeDistancesEqual(self.tree, reconstructed)

    def test_gnj(self):
        """testing gnj"""
        results = gnj(self.dists, keep=1, show_progress=False)
        (length, reconstructed) = results[0]
        self.assertTreeDistancesEqual(self.tree, reconstructed)

        results = gnj(self.dists, keep=10, show_progress=False)
        (length, reconstructed) = results[0]
        self.assertTreeDistancesEqual(self.tree, reconstructed)

        # Results should be a TreeCollection
        len(results)
        results.get_consensus_tree()

        # From GNJ paper. Pearson, Robins, Zhang 1999.
        tied_dists = {
            ("a", "b"): 3,
            ("a", "c"): 3,
            ("a", "d"): 4,
            ("a", "e"): 3,
            ("b", "c"): 3,
            ("b", "d"): 3,
            ("b", "e"): 4,
            ("c", "d"): 3,
            ("c", "e"): 3,
            ("d", "e"): 3,
        }
        results = gnj(tied_dists, keep=3, show_progress=False)
        scores = [score for (score, tree) in results]
        self.assertEqual(scores[:2], [7.75, 7.75])
        self.assertNotEqual(scores[2], 7.75)

    def test_wls(self):
        """testing wls"""
        reconstructed = wls(self.dists, a=4, show_progress=False)
        self.assertTreeDistancesEqual(self.tree, reconstructed)

    def test_truncated_wls(self):
        """testing wls with order option"""
        order = ["e", "b", "c", "d"]
        reconstructed = wls(self.dists, order=order, show_progress=False)
        self.assertEqual(set(reconstructed.get_tip_names()), set(order))

    def test_limited_wls(self):
        """testing (well, exercising at least), wls with constrained start"""
        init = make_tree(treestring="((a,c),b,d)")
        reconstructed = wls(self.dists, start=init, show_progress=False)
        self.assertEqual(len(reconstructed.get_tip_names()), 6)
        init2 = make_tree(treestring="((a,d),b,c)")
        reconstructed = wls(self.dists, start=[init, init2], show_progress=False)
        self.assertEqual(len(reconstructed.get_tip_names()), 6)
        init3 = make_tree(treestring="((a,d),b,z)")
        self.assertRaises(Exception, wls, self.dists, start=[init, init3])
        # if start tree has all seq names, should raise an error
        self.assertRaises(
            Exception,
            wls,
            self.dists,
            start=[make_tree(treestring="((a,c),b,(d,(e,f)))")],
        )

    def test_ml(self):
        """exercise the ML tree estimation"""
        from numpy.testing import assert_allclose

        aln = load_aligned_seqs(os.path.join(data_path, "brca1.fasta"), moltype="dna")
        aln = aln.take_seqs(["Human", "Mouse", "Rat", "Dog"])
        aln = aln.omit_gap_pos(allowed_gap_frac=0)
        model = get_model("JC69")
        lnL, tree = ML(model, aln).trex(a=3, k=1, show_progress=False)
        assert_allclose(lnL, -8882.217502905267)
        self.assertTrue(tree.same_topology(make_tree("(Mouse,Rat,(Human,Dog));")))


if __name__ == "__main__":
    unittest.main()
