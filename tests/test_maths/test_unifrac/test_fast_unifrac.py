#!/usr/bin/env python
"""Unit tests for fast unifrac."""

from numpy import array, logical_not
from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser
from cogent.maths.unifrac.fast_tree import (count_envs, index_tree, index_envs,
    get_branch_lengths)
from cogent.maths.unifrac.fast_unifrac import (reshape_by_name,
    meta_unifrac, shuffle_tipnames, weight_equally, weight_by_num_tips, 
    weight_by_branch_length, weight_by_num_seqs, get_all_env_names,
    consolidate_skipping_missing_matrices, consolidate_missing_zero,
    consolidate_missing_one, consolidate_skipping_missing_values,
    UniFracTreeNode, mcarlo_sig, num_comps, fast_unifrac, 
    fast_unifrac_whole_tree, PD_whole_tree, PD_generic_whole_tree,
    TEST_ON_TREE, TEST_ON_ENVS, TEST_ON_PAIRWISE, shared_branch_length,
    shared_branch_length_to_root)
from numpy.random import permutation 

__author__ = "Rob Knight and Micah Hamady"
__copyright = "Copyright 2007, the authors."
__credits__ = ["Rob Knight", "Micah Hamady", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight, Micah Hamady"
__email__ = "rob@spot.colorado.edu, hamady@colorado.edu"
__status__ = "Prototype"

class unifrac_tests(TestCase):
    """Tests of top-level functions."""
    def setUp(self):
        """Define some standard trees."""
        self.t_str = '((a:1,b:2):4,(c:3,(d:1,e:1):2):3)'
        self.t = DndParser(self.t_str, UniFracTreeNode)
        self.env_str = """
a   A   1
a   C   2
b   A   1
b   B   1
c   B   1
d   B   3
e   C   1"""
        self.env_counts = count_envs(self.env_str.splitlines())
        self.missing_env_str = """
a   A   1
a   C   2
e   C   1"""
        self.missing_env_counts = count_envs(self.missing_env_str.splitlines())
        self.extra_tip_str = """
q   A   1
w   C   2
e   A   1
r   B   1
t   B   1
y   B   3
u   C   1"""
        self.extra_tip_counts = count_envs(self.extra_tip_str.splitlines())
        self.wrong_tip_str = """
q   A   1
w   C   2
r   B   1
t   B   1
y   B   3
u   C   1"""
        self.wrong_tip_counts = count_envs(self.wrong_tip_str.splitlines())

        self.t2_str = '(((a:1,b:1):1,c:5):2,d:4)'
        self.t2 = DndParser(self.t2_str, UniFracTreeNode)
        self.env2_str = """
a   B   1
b   A   1
c   A   2
c   C   2
d   B   1
d   C   1"""
        self.env2_counts = count_envs(self.env2_str.splitlines())
        self.trees = [self.t, self.t2]
        self.envs = [self.env_counts, self.env2_counts]

        self.mc_1 = array([.5, .4, .3, .2, .1, .6, .7, .8, .9, 1.0])
       
        # from old EnvsNode tests
        self.old_t_str = '((org1:0.11,org2:0.22,(org3:0.12,org4:0.23)g:0.33)b:0.2,(org5:0.44,org6:0.55)c:0.3,org7:0.4)'

        self.old_t = DndParser(self.old_t_str, UniFracTreeNode)
        self.old_env_str = """
org1    env1    1
org1    env2    1
org2    env2    1
org3    env2    1
org4    env3    1
org5    env1    1
org6    env1    1
org7    env3    1
"""
        self.old_env_counts = count_envs(self.old_env_str.splitlines())
        self.old_node_index, self.old_nodes = index_tree(self.old_t)
        self.old_count_array, self.old_unique_envs, self.old_env_to_index, \
            self.old_node_to_index = index_envs(self.old_env_counts, self.old_node_index)
        self.old_branch_lengths = get_branch_lengths(self.old_node_index)

    def test_shared_branch_length(self):
        """Should return the correct shared branch length by env"""
        t_str = "(((a:1,b:2):3,c:4),(d:5,e:6,f:7):8);"
        envs = """
a A 1
b A 1
c A 1
d A 1
e A 1
f B 1
"""
        env_counts = count_envs(envs.splitlines())
        t = DndParser(t_str, UniFracTreeNode)
        exp = {('A',):21.0,('B',):7.0}
        obs = shared_branch_length(t, env_counts, 1)
        self.assertEqual(obs, exp)

        exp = {('A','B'):8.0}
        obs = shared_branch_length(t, env_counts, 2)
        self.assertEqual(obs, exp)

        self.assertRaises(ValueError, shared_branch_length, t, env_counts, 3)

    def test_shared_branch_length_to_root(self):
        """Should return the correct shared branch length by env to root"""
        t_str = "(((a:1,b:2):3,c:4),(d:5,e:6,f:7):8);"
        envs = """
a A 1
b A 1
c A 1
d A 1
e A 1
f B 1 
"""
        env_counts = count_envs(envs.splitlines())
        t = DndParser(t_str, UniFracTreeNode)
        exp = {'A':29.0,'B':15.0}
        obs = shared_branch_length_to_root(t, env_counts)
        self.assertEqual(obs, exp)


    def test_fast_unifrac(self):
        """Should calc unifrac values for whole tree."""
        #Note: results not tested for correctness here as detailed tests
        #in fast_tree module.
        res = fast_unifrac(self.t, self.env_counts)
        res = fast_unifrac(self.t, self.missing_env_counts)
        res = fast_unifrac(self.t, self.extra_tip_counts)
        self.assertRaises(ValueError,  fast_unifrac, self.t, \
            self.wrong_tip_counts)

    def test_fast_unifrac_whole_tree(self):
        """ should correctly compute one p-val for whole tree """
        # "should test with fake permutation but using same as old envs nodefor now"
        result = []
        num_to_do = 10
        for i in range(num_to_do):
            real_ufracs, sim_ufracs = fast_unifrac_whole_tree(self.old_t, \
                self.old_env_counts, 1000, permutation_f=permutation)
            rawp, corp = mcarlo_sig(sum(real_ufracs), [sum(x) for x in \
                sim_ufracs], 1, tail='high')
            result.append(rawp)
        self.assertSimilarMeans(result, 0.047)

    def test_PD_whole_tree(self):
        """PD_whole_tree should correctly compute PD for test tree."""
        self.t1 = DndParser('((a:1,b:2):4,(c:3,(d:1,e:1):2):3)', \
            UniFracTreeNode)
        self.env_str = """
        a   A   1
        a   C   2
        b   A   1
        b   B   1
        c   B   1
        d   B   3
        e   C   1"""
        env_counts = count_envs(self.env_str.splitlines())
        self.assertEqual(PD_whole_tree(self.t1,self.env_counts), \
            (['A','B','C'], array([7.,15.,11.])))

    def test_PD_generic_whole_tree(self):
        """PD_generic_whole_tree should correctly compute PD for test tree."""
        self.t1 = DndParser('((a:1,b:2):4,(c:3,(d:1,e:1):2):3)', \
            UniFracTreeNode)
        self.env_str = """
        a   A   1
        a   C   2
        b   A   1
        b   B   1
        c   B   1
        d   B   3
        e   C   1"""
        env_counts = count_envs(self.env_str.splitlines())
        self.assertEqual(PD_generic_whole_tree(self.t1,self.env_counts), \
            (['A','B','C'], array([7.,15.,11.])))


    def test_mcarlo_sig(self):
        """test_mcarlo_sig should calculate monte carlo sig high/low"""
        self.assertEqual(mcarlo_sig(.5, self.mc_1, 1, 'high'), (5.0/10, 5.0/10))
        self.assertEqual(mcarlo_sig(.5, self.mc_1, 1, 'low'), (4.0/10, 4.0/10))
        self.assertEqual(mcarlo_sig(.5, self.mc_1, 5, 'high'), (5.0/10, 1.0))
        self.assertEqual(mcarlo_sig(.5, self.mc_1, 5, 'low'), (4.0/10, 1.0))
        self.assertEqual(mcarlo_sig(0, self.mc_1, 1, 'low'), (0.0, "<=%.1e" % (1.0/10)))
        self.assertEqual(mcarlo_sig(100, self.mc_1, 10, 'high'), (0.0, "<=%.1e" % (1.0/10)))

    
    def test_num_comps(self):
        """ test num comps """
        self.assertEqual(num_comps(5), sum([i for i in range(1, 5)]))
        self.assertEqual(num_comps(15), sum([i for i in range(1, 15)]))
        self.assertEqual(num_comps(10000), sum([i for i in range(1, 10000)]))
        self.assertEqual(num_comps(1833), sum([i for i in range(1, 1833)]))

    def test_shuffle_tipnames(self):
        """shuffle_tipnames should return copy of tree w/ labels permuted"""
        #Note: this should never fail but is technically still stochastic
        #5! is 120 so repeating 5 times should fail about 1 in 10^10.
        for i in range(5):
            try:
                t = DndParser(self.t_str)
                result = shuffle_tipnames(t)
                orig_names = [n.Name for n in t.tips()]
                new_names = [n.Name for n in result.tips()]
                self.assertIsPermutation(orig_names, new_names)
                return
            except AssertionError:
                continue
        raise AssertionError, "Produced same permutation in 5 tries: broken?"

    def test_weight_equally(self):
        """weight_equally should return unit weight per tree"""
        self.assertEqual(weight_equally(self.trees, self.envs),
            array([1,1]))

    def test_weight_by_num_tips(self):
        """weight_by_num_tips should return tips per tree"""
        self.assertEqual(weight_by_num_tips(self.trees, self.envs),
            array([5, 4]))

    def test_weight_by_branch_length(self):
        """weight_by_branch_length should return branch length per tree"""
        self.assertEqual(weight_by_branch_length(self.trees, self.envs),
            array([17, 14]))

    def test_weight_by_num_seqs(self):
        """weight_by_num_seqs should return num seqs per tree"""
        self.assertEqual(weight_by_num_seqs(self.trees, self.envs),
            array([10, 8]))

    def test_get_all_env_names(self):
        """get_all_env_names should get all names from counts"""
        self.assertEqual(get_all_env_names(self.env_counts), 
            set('ABC'))

    def test_consolidate_skipping_missing_matrices(self):
        """consolidate_skipping_missing_matrices should skip those missing data"""
        m1 = array([[1,2],[3,4]])
        m2 = array([[1,2,3],[4,5,6],[7,8,9]])
        m3 = array([[2,2,2],[3,3,3],[4,4,4]])
        matrices = [m1,m2, m3]
        env_names = map(list, ['AB', 'ABC', 'ABC'])
        weights = [1, 2, 3]
        all_names =list('ABC')
        result = consolidate_skipping_missing_matrices(matrices, env_names, weights,
            all_names)
        self.assertFloatEqual(result, .4*m2 + .6*m3)

    def test_consolidate_missing_zero(self):
        """consolidate_missing_zero should fill missing values to zero"""
        m1 = array([[1,2],[3,4]])
        m2 = array([[1,2,3],[4,5,6],[7,8,9]])
        m3 = array([[2,2,2],[3,3,3],[4,4,4]])
        matrices = [m1,m2, m3]
        env_names = map(list, ['AB', 'ABC', 'ABC'])
        weights = [1, 2, 3]
        weights = array(weights, float)
        weights/=weights.sum()
        all_names =list('ABC')
        transformed_m1 = array([[1,2,0],[3,4,0],[0,0,0]])
        result = consolidate_missing_zero(matrices, env_names, weights,
            all_names)
        self.assertFloatEqual(result, (1/6.)*transformed_m1 + (2/6.)*m2 + (3/6.)*m3)

    def test_consolidate_missing_one(self):
        """consolidate_missing_one should fill missing off-diags to one"""
        m1 = array([[1,2],[3,4]])
        m2 = array([[1,2,3],[4,5,6],[7,8,9]])
        m3 = array([[2,2,2],[3,3,3],[4,4,4]])
        matrices = [m1,m2, m3]
        env_names = map(list, ['AB', 'ABC', 'ABC'])
        weights = [1, 2, 3]
        weights = array(weights, float)
        weights/=weights.sum()
        all_names =list('ABC')
        transformed_m1 = array([[1,2,1],[3,4,1],[1,1,0]])
        result = consolidate_missing_one(matrices, env_names, weights,
            all_names)
        self.assertFloatEqual(result, (1/6.)*transformed_m1 + (2/6.)*m2 + (3/6.)*m3)

    def test_consolidate_skipping_missing_values(self):
        """consolidate_skipping_missing_values should average over filled values"""
        m1 = array([[1,2],[3,4]])
        m2 = array([[1,2,3],[4,5,6],[7,8,9]])
        m3 = array([[2,2,2],[3,3,3],[4,4,4]])
        matrices = [m1,m2, m3]
        env_names = map(list, ['AB', 'ABC', 'ABC'])
        weights = [1., 2, 3]
        weights = array(weights)
        weights/=weights.sum()
        all_names =list('ABC')
        expected = array([[ 1/6.*1 + 2/6.*1 + 3/6.*2,
                            1/6.*2 + 2/6.*2 + 3/6.*2,
                            2/5.*3 + 3/5.*2],
                          [ 1/6.*3 + 2/6.*4 + 3/6.*3,
                            1/6.*4 + 2/6.*5 + 3/6.*3,
                            2/5.*6 + 3/5.*3],
                          [ 2/5.*7 + 3/5.*4,
                            2/5.*8 + 3/5.*4,
                            2/5.*9 + 3/5.*4]])
        result = consolidate_skipping_missing_values(matrices, env_names, weights,
            all_names)
        self.assertFloatEqual(result, expected)

    def test_reshape_by_name(self):
        """reshape_by_name should reshape matrix from old to new names"""
        old = array([[0,1,2],[3,4,5],[6,7,8]])
        old_names = 'ABC'
        new_names = 'xCyBA'
        exp = array([[0,0,0,0,0],[0,8,0,7,6],[0,0,0,0,0],\
            [0,5,0,4,3],[0,2,0,1,0]])
        self.assertEqual(reshape_by_name(old, old_names, new_names), exp)
        result = reshape_by_name(old, old_names, new_names, masked=True)
        result.fill_value=0
        self.assertEqual(result._data * logical_not(result._mask), exp)

    def test_meta_unifrac(self):
        """meta_unifrac should give correct result on sample trees"""
        tree_list = [self.t, self.t2]
        envs_list = [self.env_counts, self.env2_counts]
        result = meta_unifrac(tree_list, envs_list, weight_equally,
            modes=["distance_matrix"])

        u1_distances = array([[0, 10/16.,8/13.],[10/16.,0,8/17.],\
                            [8/13.,8/17.,0]])
        u2_distances = array([[0,11/14.,6/13.],[11/14.,0,7/13.],[6/13.,7/13., 0]])
        exp = (u1_distances + u2_distances)/2
        self.assertFloatEqual(result['distance_matrix'], (exp, list('ABC')))


if __name__ == '__main__':
    main()
