#! /usr/bin/env python
import unittest, os
import warnings
from numpy import log, exp
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from cogent.phylo.distance import EstimateDistances
from cogent.phylo.nj import nj, gnj
from cogent.phylo.least_squares import wls
from cogent import LoadSeqs, LoadTree
from cogent.phylo.tree_collection import LogLikelihoodScoredTreeCollection,\
    WeightedTreeCollection, LoadTrees, ScoredTreeCollection
from cogent.evolve.models import JC69, HKY85, F81
from cogent.phylo.consensus import majorityRule, weightedMajorityRule, \
        getSplits, getTree
from cogent.util.misc import remove_files

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2015, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield",\
        "Daniel McDonald", "Ben Kaehler"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

base_path = os.getcwd()
data_path = os.path.join(base_path, 'data')

def Tree(t):
    return LoadTree(treestring=t)

class ConsensusTests(unittest.TestCase):
    def setUp(self):
        self.trees = [
                Tree("((a,b),c,d);"),
                Tree("((a,b),c,d);"),
                Tree("((a,c),b,d);"),
                Tree("((a,b),c,d);")]
        
        weights = map(log, [0.4,0.4,0.05,0.15]) # emphasizing the a,b clade
        self.scored_trees = zip(weights, self.trees)
        self.scored_trees.sort(reverse=True)

        self.rooted_trees = [
                Tree("((a,b),(c,d));"),
                Tree("((a,b),(c,d));"),
                Tree("((a,c),(b,d));"),
                Tree("((a,b),c,d);")]

        self.trees_randomly_rooted = \
                [(-3416.3820971172017, Tree("(((F,(M,(P,(K,J)))),A),Q,(H,((S,(X,((Z,N),((O,(L,D)),E)))),((G,((I,(R,T)),Y)),(C,(((U,W),V),B))))));")), 
                 (-3416.3820974814785, Tree("(R,T,(I,(Y,(G,((((X,((Z,N),((O,(L,D)),E))),S),((Q,((F,(M,(P,(K,J)))),A)),H)),((((U,W),V),C),B))))));")), 
                 (-3416.571172739171, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,(S,((X,((Z,N),((O,(L,D)),E))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B)))))));")), 
                 (-3416.5721393589474, Tree("(P,(K,J),(M,(F,(A,(Q,(H,(((((Z,N),((O,(L,D)),E)),S),X),((C,(G,((I,(R,T)),Y))),(((U,W),V),B)))))))));")), 
                 (-3416.5721493054643, Tree("(R,T,(I,(Y,(G,(C,((S,(X,((Z,N),((O,(L,D)),E)))),((((U,W),V),B),((Q,((F,(M,(P,(K,J)))),A)),H))))))));")), 
                 (-3416.572149424391, Tree("((U,W),V,(B,(((G,((I,(R,T)),Y)),C),((S,(X,((Z,N),((O,(L,D)),E)))),(H,(Q,((F,(M,(P,(K,J)))),A)))))));")), 
                 (-3416.596424489769, Tree("(S,(X,((Z,N),((O,(L,D)),E))),((((G,((I,(R,T)),Y)),C),(((U,W),V),B)),((H,Q),((F,(M,(P,(K,J)))),A))));")),
                 (-3416.6335053806333, Tree("(M,(P,(K,J)),(F,(A,(Q,(H,(((((((O,(L,D)),E),N),Z),X),S),((C,((V,(U,W)),B)),(G,((I,(R,T)),Y)))))))));")), 
                 (-3416.7626687401867, Tree("(Z,N,(((O,(L,D)),E),(X,(S,((((F,(M,(P,(K,J)))),A),Q),(H,((G,((I,(R,T)),Y)),(C,(((U,W),V),B)))))))));")), 
                 (-3416.7670692165866, Tree("(X,((Z,N),((O,(L,D)),E)),(S,((B,(((G,((I,(R,T)),Y)),C),((U,W),V))),((Q,((F,(M,(P,(K,J)))),A)),H))));")), 
                 (-3416.7670696377254, Tree("(P,(K,J),(M,(F,(A,(Q,(H,(((X,((Z,N),((O,(L,D)),E))),S),(((U,W),V),(((G,((I,(R,T)),Y)),C),B)))))))));")), 
                 (-3416.848062302587, Tree("(((O,(L,D)),E),N,(Z,(X,(S,(((((F,(M,(P,(K,J)))),A),Q),H),(((V,(U,W)),B),(C,(G,((I,(R,T)),Y)))))))));")), 
                 (-3416.9943503002764, Tree("((U,W),V,(B,((C,(G,((I,(R,T)),Y))),(H,(((((Z,N),((O,(L,D)),E)),S),X),(((F,(M,(P,(K,J)))),A),Q))))));")), 
                 (-3417.014782481302, Tree("(Q,((F,(M,(P,(K,J)))),A),(((((U,W),V),B),((G,((I,(R,T)),Y)),C)),((S,(X,((Z,N),((O,(L,D)),E)))),H)));")), 
                 (-3417.015470262783, Tree("(L,D,(O,(E,((Z,N),(X,(S,((Q,((F,(M,(P,(K,J)))),A)),(H,((((U,W),V),B),((G,((I,(R,T)),Y)),C))))))))));")), 
                 (-3417.241619414339, Tree("(Z,N,(((O,(L,D)),E),(X,(S,(H,((Q,((F,(M,(P,(K,J)))),A)),(B,(((G,((I,(R,T)),Y)),C),((U,W),V)))))))));")), 
                 (-3417.242009280534, Tree("(X,((Z,N),((O,(L,D)),E)),(S,((Q,((F,(M,(P,(K,J)))),A)),(H,(B,(((G,((I,(R,T)),Y)),C),((U,W),V)))))));")), 
                 (-3417.2637092520818, Tree("(K,J,(P,(M,(F,(A,(Q,(((((((O,(L,D)),E),N),Z),X),S),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))))))));")), 
                 (-3417.420526572887, Tree("(B,(V,(U,W)),(((H,(Q,(A,(F,(M,(P,(K,J))))))),((((((O,(L,D)),E),N),X),Z),S)),(C,(((I,(R,T)),Y),G))));")), 
                 (-3417.4205266767162, Tree("(R,T,(I,(Y,(G,(C,((B,(V,(U,W))),((H,((A,(F,(M,(P,(K,J))))),Q)),(Z,(((((O,(L,D)),E),N),X),S)))))))));")), 
                 (-3417.620921910812, Tree("((O,(L,D)),E,(N,(X,(Z,(S,(((Q,(A,(F,(M,(P,(K,J)))))),H),((((((I,(R,T)),Y),G),C),(V,(U,W))),B)))))));")), 
                 (-3417.6209219461302, Tree("(F,(M,(P,(K,J))),(A,(Q,(H,((Z,(((((O,(L,D)),E),N),X),S)),((((G,((I,(R,T)),Y)),C),(V,(U,W))),B))))));")), 
                 (-3417.6209224304744, Tree("(H,((A,(F,(M,(P,(K,J))))),Q),(B,((Z,(((((O,(L,D)),E),N),X),S)),(((G,((I,(R,T)),Y)),C),(V,(U,W))))));")), 
                 (-3417.9379715010136, Tree("((I,(R,T)),Y,(G,(C,((B,(V,(U,W))),(H,((Q,(A,(F,(M,(P,(K,J)))))),((((((O,(L,D)),E),N),X),Z),S)))))));")), 
                 (-3417.9379715187215, Tree("(((O,(L,D)),E),N,(X,(S,(Z,(((A,(F,(M,(P,(K,J))))),Q),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))))));"))]
        self.trees_rooted_at_A = \
                [(-3416.3820971172017, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,((S,(X,((Z,N),((O,(L,D)),E)))),((G,((I,(R,T)),Y)),(C,(((U,W),V),B)))))));")), 
                 (-3416.3820974814785, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,(((X,((Z,N),((O,(L,D)),E))),S),((G,((I,(R,T)),Y)),((((U,W),V),C),B))))));")), 
                 (-3416.571172739171, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,(S,((X,((Z,N),((O,(L,D)),E))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B)))))));")), 
                 (-3416.5721393589474, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,(((((Z,N),((O,(L,D)),E)),S),X),((C,(G,((I,(R,T)),Y))),(((U,W),V),B))))));")), 
                 (-3416.5721493054643, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,((((U,W),V),B),(((G,((I,(R,T)),Y)),C),(S,(X,((Z,N),((O,(L,D)),E)))))))));")), 
                 (-3416.572149424391, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,((S,(X,((Z,N),((O,(L,D)),E)))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B))))));")), 
                 (-3416.596424489769, Tree("((F,(M,(P,(K,J)))),A,(((S,(X,((Z,N),((O,(L,D)),E)))),(((G,((I,(R,T)),Y)),C),(((U,W),V),B))),(H,Q)));")),
                 (-3416.6335053806333, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,(((((((O,(L,D)),E),N),Z),X),S),((C,((V,(U,W)),B)),(G,((I,(R,T)),Y)))))));")),
                 (-3416.7626687401867, Tree("((F,(M,(P,(K,J)))),A,(Q,((H,((G,((I,(R,T)),Y)),(C,(((U,W),V),B)))),(S,(X,((Z,N),((O,(L,D)),E)))))));")),
                 (-3416.7670692165866, Tree("((F,(M,(P,(K,J)))),A,(Q,(((B,(((G,((I,(R,T)),Y)),C),((U,W),V))),((X,((Z,N),((O,(L,D)),E))),S)),H)));")),
                 (-3416.7670696377254, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,(((X,((Z,N),((O,(L,D)),E))),S),(((U,W),V),(((G,((I,(R,T)),Y)),C),B))))));")),
                 (-3416.848062302587, Tree("((F,(M,(P,(K,J)))),A,(Q,(H,(((((((O,(L,D)),E),N),Z),X),S),(((V,(U,W)),B),(C,(G,((I,(R,T)),Y))))))));")),
                 (-3416.9943503002764, Tree("((F,(M,(P,(K,J)))),A,(((((C,(G,((I,(R,T)),Y))),(((U,W),V),B)),H),((((Z,N),((O,(L,D)),E)),S),X)),Q));")),
                 (-3417.014782481302, Tree("((F,(M,(P,(K,J)))),A,(Q,(((((U,W),V),B),((G,((I,(R,T)),Y)),C)),((S,(X,((Z,N),((O,(L,D)),E)))),H))));")),
                 (-3417.015470262783, Tree("((F,(M,(P,(K,J)))),A,(Q,((H,((((U,W),V),B),((G,((I,(R,T)),Y)),C))),(S,(X,((Z,N),((O,(L,D)),E)))))));")),
                 (-3417.241619414339, Tree("((F,(M,(P,(K,J)))),A,(Q,((H,((X,((Z,N),((O,(L,D)),E))),S)),(B,(((G,((I,(R,T)),Y)),C),((U,W),V))))));")),
                 (-3417.242009280534, Tree("((F,(M,(P,(K,J)))),A,(Q,(((X,((Z,N),((O,(L,D)),E))),S),(H,(B,(((G,((I,(R,T)),Y)),C),((U,W),V)))))));")),
                 (-3417.2637092520818, Tree("((F,(M,(P,(K,J)))),A,(Q,(((((((O,(L,D)),E),N),Z),X),S),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))));")),
                 (-3417.420526572887, Tree("(A,(F,(M,(P,(K,J)))),(Q,(H,(((((((O,(L,D)),E),N),X),Z),S),((B,(V,(U,W))),(C,(((I,(R,T)),Y),G)))))));")),
                 (-3417.4205266767162, Tree("(A,(F,(M,(P,(K,J)))),(Q,(H,(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),(Z,(((((O,(L,D)),E),N),X),S))))));")),
                 (-3417.620921910812, Tree("(A,(F,(M,(P,(K,J)))),(Q,(H,(((((((O,(L,D)),E),N),X),Z),S),((((((I,(R,T)),Y),G),C),(V,(U,W))),B)))));")),
                 (-3417.6209219461302, Tree("(A,(F,(M,(P,(K,J)))),(Q,(H,((Z,(((((O,(L,D)),E),N),X),S)),((((G,((I,(R,T)),Y)),C),(V,(U,W))),B)))));")),
                 (-3417.6209224304744, Tree("(A,(F,(M,(P,(K,J)))),(Q,(H,(B,((Z,(((((O,(L,D)),E),N),X),S)),(((G,((I,(R,T)),Y)),C),(V,(U,W))))))));")),
                 (-3417.9379715010136, Tree("(A,(F,(M,(P,(K,J)))),(Q,((H,((B,(V,(U,W))),(C,(((I,(R,T)),Y),G)))),((((((O,(L,D)),E),N),X),Z),S))));")),
                 (-3417.9379715187215, Tree("(A,(F,(M,(P,(K,J)))),(Q,((Z,(((((O,(L,D)),E),N),X),S)),(((B,(V,(U,W))),(C,(G,((I,(R,T)),Y)))),H))));"))]

        self.unrooted_conflicting_trees = [
                Tree('((a,b),c,d);'),
                Tree('((a,c),b,d);'),
                Tree('((a,d),b,c);')]

        self.rooted_conflicting_trees = [
                Tree('((a,b),(c,d));'),
                Tree('((a,c),(b,d));'),
                Tree('((a,d),(b,c));')]

        self.unrooted_trees_lengths = [
                (2, Tree('((a:0.3,c:0.4):0.5,b:0.2,d:0.1);')),
                (1, Tree('((a:0.1,b:0.1):0.1,c:0.1,d:0.1);'))]

        self.rooted_trees_lengths = [
                (2, Tree('((a:0.3,c:0.4):0.2,(b:0.2,d:0.1):0.3);')),
                (1, Tree('((a:0.1,b:0.1):0.05,(c:0.1,d:0.1):0.05);'))]
    
    def test_majorityRule(self):
        """Tests for majority rule consensus trees"""
        trees = self.rooted_trees
        outtrees = majorityRule(trees, strict=False)
        self.assertEqual(len(outtrees), 1)
        self.assert_(outtrees[0].sameTopology(Tree("((c,d),(a,b));")))
        outtrees = majorityRule(trees, strict=True)
        self.assertEqual(len(outtrees), 1)
        self.assert_(outtrees[0].sameTopology(Tree("(c,d,(a,b));")))

    def test_get_tree_get_splits(self):
        """getTree should provide a reciprocal map of getSplits"""
        tree = LoadTree(filename=os.path.join(data_path,"murphy.tree"))
        self.assertTrue(tree.sameTopology(getTree(getSplits(tree))))
    
    def test_consensus_tree_branch_lengths(self):
        """consensus trees should average branch lengths properly"""
        def get_ac(tree):
            for edge in tree.getEdgeVector(include_root=False):
                if set('ac') == set([c.Name for c in edge.Children]):
                    return edge

        sct = ScoredTreeCollection(self.unrooted_trees_lengths)
        ct = sct.getConsensusTree()
        maj_tree = self.unrooted_trees_lengths[0][1]
        self.assertTrue(abs(get_ac(ct).Length-get_ac(maj_tree).Length) < 1e-9)

        sct = ScoredTreeCollection(self.rooted_trees_lengths)
        ct = sct.getConsensusTree(method='rooted')
        maj_tree = self.rooted_trees_lengths[0][1]
        self.assertTrue(abs(get_ac(ct).Length-get_ac(maj_tree).Length) < 1e-9)

    def test_consensus_from_scored_trees_collection(self):
        """tree collection should get same consensus as direct approach"""
        tree_list = [(1, t) for t in self.trees]
        sct = LogLikelihoodScoredTreeCollection(tree_list)
        ct = sct.getConsensusTree()
        self.assertTrue(ct.sameTopology(Tree("((c,d),a,b);")))

    def test_consensus_from_scored_trees_collection_ii(self):
        """strict consensus should handle conflicting trees"""
        sct = ScoredTreeCollection(zip([1]*3, self.unrooted_conflicting_trees))
        ct = sct.getConsensusTrees()[0]
        self.assertTrue(ct.sameTopology(Tree("(a,b,c,d);")))

        sct = ScoredTreeCollection(zip([1]*3, self.rooted_conflicting_trees))
        #cts = sct.getConsensusTrees(method='rooted')
        ct = sct.getConsensusTrees(method='rooted')[0]
        self.assertTrue(ct.sameTopology(Tree("(a,b,c,d);")))
        #for tree in cts:
        #    print str(tree)
        #self.assertTrue(set(map(str, cts))==set(['('+c+');' for c in 'abcd']))
    
    def test_weighted_consensus_from_scored_trees_collection(self):
        """weighted consensus from a tree collection should be different"""
        sct = LogLikelihoodScoredTreeCollection(self.scored_trees)
        ct = sct.getConsensusTree()
        self.assertTrue(ct.sameTopology(Tree("((a,b),c,d);")))

    def test_weighted_consensus_from_scored_trees_collection_ii(self):
        """root positions in input tree collection should not effect result"""
        sct = LogLikelihoodScoredTreeCollection(self.trees_randomly_rooted)
        ctrr = sct.getConsensusTree()
        sct = LogLikelihoodScoredTreeCollection(self.trees_rooted_at_A)
        ctra = sct.getConsensusTree()
        self.assertTrue(ctrr.sameTopology(ctra))
    
    def test_weighted_trees_satisyfing_cutoff(self):
        """build consensus tree from those satisfying cutoff"""
        sct = LogLikelihoodScoredTreeCollection(self.scored_trees)
        cts = sct.getWeightedTrees(cutoff=0.8)
        for weight, tree in cts:
            self.assertTrue(tree.sameTopology(Tree('((a,b),c,d);')))
        
        ct = cts.getConsensusTree()
        self.assertTrue(ct.sameTopology(Tree("((a,b),c,d);")))
    
    def test_tree_collection_read_write_file(self):
        """should correctly read / write a collection from a file"""
        def eval_klass(coll):
            coll.writeToFile('sample.trees')
            read = LoadTrees('sample.trees')
            self.assertTrue(type(read) == type(coll))
        
        eval_klass(LogLikelihoodScoredTreeCollection(self.scored_trees))
        
        # convert lnL into p
        eval_klass(WeightedTreeCollection([(exp(s), t) 
                                    for s,t in self.scored_trees]))
        remove_files(['sample.trees'], error_on_missing=False)
    

class TreeReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.tree = LoadTree(treestring='((a:3,b:4):2,(c:6,d:7):30,(e:5,f:5):5)')
        self.dists = self.tree.getDistances()
        
    def assertTreeDistancesEqual(self, t1, t2):
        d1 = t1.getDistances()
        d2 = t2.getDistances()
        self.assertEqual(len(d1), len(d2))
        for key in d2:
            self.assertAlmostEqual(d1[key], d2[key])

    def test_nj(self):
        """testing nj"""
        reconstructed = nj(self.dists)
        self.assertTreeDistancesEqual(self.tree, reconstructed)
        
    def test_gnj(self):
        """testing gnj"""
        results = gnj(self.dists, keep=1)
        (length, reconstructed) = results[0]
        self.assertTreeDistancesEqual(self.tree, reconstructed)
        
        results = gnj(self.dists, keep=10)
        (length, reconstructed) = results[0]
        self.assertTreeDistancesEqual(self.tree, reconstructed)
        
        # Results should be a TreeCollection
        len(results)
        results.getConsensusTree()

        # From GNJ paper. Pearson, Robins, Zhang 1999.
        tied_dists = {
                ('a', 'b'):3, ('a', 'c'):3, ('a', 'd'):4, ('a', 'e'):3, 
                ('b', 'c'):3, ('b', 'd'):3, ('b', 'e'):4,
                ('c', 'd'):3, ('c', 'e'):3, 
                ('d', 'e'):3}
        results = gnj(tied_dists, keep=3)
        scores = [score for (score, tree) in results]
        self.assertEqual(scores[:2], [7.75, 7.75])
        self.assertNotEqual(scores[2], 7.75)

    def test_wls(self):
        """testing wls"""
        reconstructed = wls(self.dists, a=4)
        self.assertTreeDistancesEqual(self.tree, reconstructed)

    def test_truncated_wls(self):
        """testing wls with order option"""
        order = ['e', 'b', 'c', 'd']
        reconstructed = wls(self.dists, order=order)
        self.assertEqual(set(reconstructed.getTipNames()), set(order))

    def test_limited_wls(self):
        """testing (well, exercising at least), wls with constrained start"""
        init = LoadTree(treestring='((a,c),b,d)')
        reconstructed = wls(self.dists, start=init)
        self.assertEqual(len(reconstructed.getTipNames()), 6)
        init2 = LoadTree(treestring='((a,d),b,c)')
        reconstructed = wls(self.dists, start=[init, init2])
        self.assertEqual(len(reconstructed.getTipNames()), 6)
        init3 = LoadTree(treestring='((a,d),b,z)')
        self.assertRaises(Exception, wls, self.dists, start=[init, init3])
        # if start tree has all seq names, should raise an error
        self.assertRaises(Exception, wls, self.dists,
                start=[LoadTree(treestring='((a,c),b,(d,(e,f)))')])
        
    
class DistancesTests(unittest.TestCase):
    def setUp(self):
        self.al = LoadSeqs(data = {'a':'GTACGTACGATC',
                            'b':'GTACGTACGTAC',
                            'c':'GTACGTACGTTC',
                            'e':'GTACGTACTGGT'})
        self.collection = LoadSeqs(data = {'a':'GTACGTACGATC',
                            'b':'GTACGTACGTAC',
                            'c':'GTACGTACGTTC',
                            'e':'GTACGTACTGGT'}, aligned=False)
    
    def assertDistsAlmostEqual(self, expected, observed, precision=4):
        observed = dict([(frozenset(k),v) for (k,v) in observed.items()])
        expected = dict([(frozenset(k),v) for (k,v) in expected.items()])
        for key in expected:
            self.assertAlmostEqual(expected[key], observed[key], precision)
            
    def test_EstimateDistances(self):
        """testing (well, exercising at least), EstimateDistances"""
        d = EstimateDistances(self.al, JC69())
        d.run()
        canned_result = {('b', 'e'): 0.440840,
                        ('c', 'e'): 0.440840,
                        ('a', 'c'): 0.088337,
                        ('a', 'b'): 0.188486,
                        ('a', 'e'): 0.440840,
                        ('b', 'c'): 0.0883373}
        result = d.getPairwiseDistances()
        self.assertDistsAlmostEqual(canned_result, result)
        
        # excercise writing to file
        d.writeToFile('junk.txt')
        try:
            os.remove('junk.txt')
        except OSError:
            pass # probably parallel
    
    def test_EstimateDistancesWithMotifProbs(self):
        """EstimateDistances with supplied motif probs"""
        motif_probs= {'A':0.1,'C':0.2,'G':0.2,'T':0.5}
        d = EstimateDistances(self.al, HKY85(), motif_probs=motif_probs)
        d.run()
        canned_result = {('a', 'c'): 0.07537,
                       ('b', 'c'): 0.07537,
                        ('a', 'e'): 0.39921,
                        ('a', 'b'): 0.15096,
                        ('b', 'e'): 0.39921,
                        ('c', 'e'): 0.37243}
        result = d.getPairwiseDistances()
        self.assertDistsAlmostEqual(canned_result, result)
    
    def test_EstimateDistances_fromThreeway(self):
        """testing (well, exercising at least), EsimateDistances fromThreeway"""
        d = EstimateDistances(self.al, JC69(), threeway=True)
        d.run()
        canned_result = {('b', 'e'): 0.495312,
                        ('c', 'e'): 0.479380,
                        ('a', 'c'): 0.089934,
                        ('a', 'b'): 0.190021,
                        ('a', 'e'): 0.495305,
                        ('b', 'c'): 0.0899339}
        result = d.getPairwiseDistances(summary_function="mean")
        self.assertDistsAlmostEqual(canned_result, result)
    
    def test_EstimateDistances_fromUnaligned(self):
        """Excercising estimate distances from unaligned sequences"""
        d = EstimateDistances(self.collection, JC69(), do_pair_align=True,
                                rigorous_align=True)
        d.run()
        canned_result = {('b', 'e'): 0.440840,
                        ('c', 'e'): 0.440840,
                        ('a', 'c'): 0.088337,
                        ('a', 'b'): 0.188486,
                        ('a', 'e'): 0.440840,
                        ('b', 'c'): 0.0883373}
        result = d.getPairwiseDistances()
        self.assertDistsAlmostEqual(canned_result, result)
        
        d = EstimateDistances(self.collection, JC69(), do_pair_align=True,
                                rigorous_align=False)
        d.run()
        canned_result = {('b', 'e'): 0.440840,
                        ('c', 'e'): 0.440840,
                        ('a', 'c'): 0.088337,
                        ('a', 'b'): 0.188486,
                        ('a', 'e'): 0.440840,
                        ('b', 'c'): 0.0883373}
        result = d.getPairwiseDistances()
        self.assertDistsAlmostEqual(canned_result, result)
    
    def test_EstimateDistances_other_model_params(self):
        """test getting other model params from EstimateDistances"""
        d = EstimateDistances(self.al, HKY85(), est_params=['kappa'])
        d.run()
        # this will be a Number object with Mean, Median etc ..
        kappa = d.getParamValues('kappa')
        self.assertAlmostEqual(kappa.Mean, 0.8939, 4)
        # this will be a dict with pairwise instances, it's called by the above
        # method, so the correctness of it's values is already checked
        kappa = d.getPairwiseParam('kappa')
    
    def test_EstimateDistances_modify_lf(self):
        """tests modifying the lf"""
        def constrain_fit(lf):
            lf.setParamRule('kappa', is_constant=True)
            lf.optimise(local=True)
            return lf
        
        d = EstimateDistances(self.al, HKY85(), modify_lf=constrain_fit)
        d.run()
        result = d.getPairwiseDistances()
        d = EstimateDistances(self.al, F81())
        d.run()
        expect = d.getPairwiseDistances()
        self.assertDistsAlmostEqual(expect, result)
        
    def test_get_raw_estimates(self):
        """correctly return raw result object"""
        d = EstimateDistances(self.al, HKY85(), est_params=['kappa'])
        d.run()
        expect = {('a', 'b'): {'kappa': 1.0000226766004808e-06, 'length': 0.18232155856115662},
                 ('a', 'c'): {'kappa': 1.0010380037049357e-06, 'length': 0.087070406623635604},
                 ('a', 'e'): {'kappa': 2.3965871843412687, 'length': 0.4389176272584539},
                 ('b', 'e'): {'kappa': 2.3965871854366592, 'length': 0.43891762729173389},
                 ('c', 'b'): {'kappa': 1.0010380037049357e-06, 'length': 0.087070406623635604},
                 ('c', 'e'): {'kappa': 0.57046787478038707, 'length': 0.43260232210282784}}
        got = d.getAllParamValues()
        for pair in expect:
            for param in expect[pair]:
                self.assertAlmostEqual(got[pair][param], expect[pair][param], places=6)

if __name__ == '__main__':
    unittest.main()
