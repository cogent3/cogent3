#!/usr/bin/env python
"""Unit tests for fast tree."""

from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser
from cogent.maths.unifrac.fast_tree import (count_envs, sum_env_dict, 
    index_envs, get_branch_lengths, index_tree, bind_to_array, 
    bind_to_parent_array, _is_parent_empty, delete_empty_parents,
    traverse_reduce, bool_descendants, sum_descendants, fitch_descendants, 
    tip_distances, UniFracTreeNode, FitchCounter, FitchCounterDense,
    permute_selected_rows, prep_items_for_jackknife, jackknife_bool, 
    jackknife_int, unifrac, unnormalized_unifrac, PD, G, unnormalized_G, 
    unifrac_matrix, unifrac_vector, PD_vector, weighted_unifrac, 
    weighted_unifrac_matrix, weighted_unifrac_vector, jackknife_array, 
    env_unique_fraction)
from numpy import (arange, reshape, zeros, logical_or, array, sum, nonzero, 
    flatnonzero, newaxis)
from numpy.random import permutation    

__author__ = "Rob Knight and Micah Hamady"
__copyright = "Copyright 2007, the authors."
__credits__ = ["Rob Knight", "Micah Hamady"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Rob Knight, Micah Hamady"
__email__ = "rob@spot.colorado.edu, hamady@colorado.edu"
__status__ = "Prototype"

class fast_tree_tests(TestCase):
    """Tests of top-level functions"""
    def setUp(self):
        """Define a couple of standard trees"""
        self.t1 = DndParser('(((a,b),c),(d,e))', UniFracTreeNode)
        self.t2 = DndParser('(((a,b),(c,d)),(e,f))', UniFracTreeNode)
        self.t3 = DndParser('(((a,b,c),(d)),(e,f))', UniFracTreeNode)
        self.t4 = DndParser('((c)b,((f,g,h)e,i)d)', UniFracTreeNode)
        self.t4.Name = 'a'
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
        self.node_index, self.nodes = index_tree(self.t)
        self.count_array, self.unique_envs, self.env_to_index, \
            self.node_to_index = index_envs(self.env_counts, self.node_index)
        self.branch_lengths = get_branch_lengths(self.node_index)

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




    def test_traverse(self):
        """traverse should work iterative or recursive"""
        stti = self.t4.traverse
        stt = self.t4.traverse_recursive
        obs = [i.Name for i in stt(self_before=False, self_after=False)]
        exp = [i.Name for i in stti(self_before=False, self_after=False)]
        self.assertEqual(obs, exp)
        obs = [i.Name for i in stt(self_before=True, self_after=False)]
        exp = [i.Name for i in stti(self_before=True, self_after=False)]
        self.assertEqual(obs, exp)
        obs = [i.Name for i in stt(self_before=False, self_after=True)]
        exp = [i.Name for i in stti(self_before=False, self_after=True)]
        self.assertEqual(obs, exp)
        obs = [i.Name for i in stt(self_before=True, self_after=True)]
        exp = [i.Name for i in stti(self_before=True, self_after=True)]
        self.assertEqual(obs, exp)

    def test_count_envs(self):
        """count_envs should return correct counts from lines"""
        envs = """
a   A   3   some other junk
a   B 
a   C   1
b   A   2

skip
c   B
d
b   A   99
"""
        result = count_envs(envs.splitlines())
        self.assertEqual(result, \
            {'a':{'A':3,'B':1,'C':1},'b':{'A':99},'c':{'B':1}})

    def test_sum_env_dict(self):
        """sum_env_dict should return correct counts from env_dict"""
        envs = """
a   A   3   some other junk
a   B 
a   C   1
b   A   2

skip
c   B
d
b   A   99
"""
        result = count_envs(envs.splitlines())
        sum_ = sum_env_dict(result)
        self.assertEqual(sum_, 105) 

    def test_index_envs(self):
        """index_envs should map envs and taxa onto indices"""
        self.assertEqual(self.unique_envs, ['A','B','C'])
        self.assertEqual(self.env_to_index, {'A':0, 'B':1, 'C':2})
        self.assertEqual(self.node_to_index,{'a':0, 'b':1, 'c':4, 'd':2, 'e':3})
        self.assertEqual(self.count_array, \
            array([[1,0,2],[1,1,0],[0,3,0],[0,0,1], \
            [0,1,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]))

    def test_get_branch_lengths(self):
        """get_branch_lengths should make array of branch lengths from index"""
        result = get_branch_lengths(self.node_index)
        self.assertEqual(result, array([1,2,1,1,3,2,4,3,0]))

    def test_env_unique_fraction(self):
        """should report unique fraction of bl in each env """
        # testing old unique fraction   
        cur_count_array = self.count_array.copy()
        bound_indices = bind_to_array(self.nodes, cur_count_array) 
        total_bl = sum(self.branch_lengths)
        bool_descendants(bound_indices)
        env_bl_sums, env_bl_ufracs = env_unique_fraction(self.branch_lengths, cur_count_array)
        # env A has 0 unique bl, B has 4, C has 1        
        self.assertEqual(env_bl_sums, [0,4,1])
        self.assertEqual(env_bl_ufracs, [0,4/17.0,1/17.0])

        cur_count_array = self.old_count_array.copy()
        bound_indices = bind_to_array(self.old_nodes, cur_count_array) 
        total_bl = sum(self.old_branch_lengths)
        bool_descendants(bound_indices)

        env_bl_sums, env_bl_ufracs = env_unique_fraction(self.old_branch_lengths, cur_count_array)
        # env A has 0 unique bl, B has 4, C has 1        
        self.assertEqual(env_bl_sums, env_bl_sums) 
        self.assertEqual(env_bl_sums, [1.29, 0.33999999999999997, 0.63])
        self.assertEqual(env_bl_ufracs, [1.29/2.9,0.33999999999999997/2.9, 0.63/2.9])

    def test_index_tree(self):
        """index_tree should produce correct index and node map"""
        #test for first tree: contains singleton outgroup
        t1 = self.t1
        id_1, child_1 = index_tree(t1)
        nodes_1 = [n._leaf_index for n in t1.traverse(self_before=False, \
            self_after=True)]
        self.assertEqual(nodes_1, [0,1,2,3,6,4,5,7,8])
        self.assertEqual(child_1, [(2,0,1),(6,2,3),(7,4,5),(8,6,7)])
        #test for second tree: strictly bifurcating
        t2 = self.t2
        id_2, child_2 = index_tree(t2)
        nodes_2 = [n._leaf_index for n in t2.traverse(self_before=False, \
            self_after=True)]
        self.assertEqual(nodes_2, [0,1,4,2,3,5,8,6,7,9,10])
        self.assertEqual(child_2, [(4,0,1),(5,2,3),(8,4,5),(9,6,7),(10,8,9)])
        #test for third tree: contains trifurcation and single-child parent
        t3 = self.t3
        id_3, child_3 = index_tree(t3)
        nodes_3 = [n._leaf_index for n in t3.traverse(self_before=False, \
            self_after=True)]
        self.assertEqual(nodes_3, [0,1,2,4,3,5,8,6,7,9,10])
        self.assertEqual(child_3, [(4,0,2),(5,3,3),(8,4,5),(9,6,7),(10,8,9)])

    def test_bind_to_array(self):
        """bind_to_array should return correct array ranges"""
        a = reshape(arange(33), (11,3))
        id_, child = index_tree(self.t3)
        bindings = bind_to_array(child, a)
        self.assertEqual(len(bindings), 5)
        self.assertEqual(bindings[0][0], a[4])
        self.assertEqual(bindings[0][1], a[0:3])
        self.assertEqual(bindings[0][1].shape, (3,3))
        self.assertEqual(bindings[1][0], a[5])
        self.assertEqual(bindings[1][1], a[3:4])
        self.assertEqual(bindings[1][1].shape, (1,3))
        self.assertEqual(bindings[2][0], a[8])
        self.assertEqual(bindings[2][1], a[4:6])
        self.assertEqual(bindings[2][1].shape, (2,3))
        self.assertEqual(bindings[3][0], a[9])
        self.assertEqual(bindings[3][1], a[6:8])
        self.assertEqual(bindings[3][1].shape, (2,3))
        self.assertEqual(bindings[4][0], a[10])
        self.assertEqual(bindings[4][1], a[8:10])
        self.assertEqual(bindings[4][1].shape, (2,3))

    def test_bind_to_parent_array(self):
        """bind_to_parent_array should bind tree to array correctly"""
        a = reshape(arange(33), (11,3))
        index_tree(self.t3)
        bindings = bind_to_parent_array(self.t3, a)
        self.assertEqual(len(bindings), 10)
        self.assertEqual(bindings[0][0], a[8])
        self.assertEqual(bindings[0][1], a[10])
        self.assertEqual(bindings[1][0], a[4])
        self.assertEqual(bindings[1][1], a[8])
        self.assertEqual(bindings[2][0], a[0])
        self.assertEqual(bindings[2][1], a[4])
        self.assertEqual(bindings[3][0], a[1])
        self.assertEqual(bindings[3][1], a[4])
        self.assertEqual(bindings[4][0], a[2])
        self.assertEqual(bindings[4][1], a[4])
        self.assertEqual(bindings[5][0], a[5])
        self.assertEqual(bindings[5][1], a[8])
        self.assertEqual(bindings[6][0], a[3])
        self.assertEqual(bindings[6][1], a[5])
        self.assertEqual(bindings[7][0], a[9])
        self.assertEqual(bindings[7][1], a[10])
        self.assertEqual(bindings[8][0], a[6])
        self.assertEqual(bindings[8][1], a[9])
        self.assertEqual(bindings[9][0], a[7])
        self.assertEqual(bindings[9][1], a[9])

    def test_delete_empty_parents(self):
        """delete_empty_parents should remove empty parents from bound indices"""
        id_to_node, node_first_last = index_tree(self.t)
        bound_indices = bind_to_array(node_first_last, self.count_array[:,0:1])
        bool_descendants(bound_indices)
        self.assertEqual(len(bound_indices), 4)
        deleted = delete_empty_parents(bound_indices)
        self.assertEqual(len(deleted), 2)
        for d in deleted:
            self.assertEqual(d[0][0], 1)

    def test_traverse_reduce(self):
        """traverse_reduce should reduce array in traversal order."""
        id_, child = index_tree(self.t3)
        a = zeros((11,3)) + 99    #fill with junk
        bindings = bind_to_array(child, a)
        #load in leaf envs
        a[0] = a[1] = a[2] = a[7] = [0,1,0]
        a[3] = [1,0,0]
        a[6] = [0,0,1]
        f = logical_or.reduce
        traverse_reduce(bindings, f)
        self.assertEqual(a,\
            array([[0,1,0],[0,1,0],[0,1,0],[1,0,0],[0,1,0],[1,0,0],\
            [0,0,1],[0,1,0],[1,1,0],[0,1,1],[1,1,1]])
        )
        f = sum
        traverse_reduce(bindings, f)
        self.assertEqual( a, \
            array([[0,1,0],[0,1,0],[0,1,0],[1,0,0],[0,3,0],[1,0,0],\
            [0,0,1],[0,1,0],[1,3,0],[0,1,1],[1,4,1]])
        )

    def test_bool_descendants(self):
        """bool_descendants should be true if any descendant true"""
        #self.t3 = DndParser('(((a,b,c),(d)),(e,f))', UniFracTreeNode)
        id_, child = index_tree(self.t3)
        a = zeros((11,3)) + 99    #fill with junk
        bindings = bind_to_array(child, a)
        #load in leaf envs
        a[0] = a[1] = a[2] = a[7] = [0,1,0]
        a[3] = [1,0,0]
        a[6] = [0,0,1]
        bool_descendants(bindings)
        self.assertEqual(a, \
            array([[0,1,0],[0,1,0],[0,1,0],[1,0,0],[0,1,0],[1,0,0],\
            [0,0,1],[0,1,0],[1,1,0],[0,1,1],[1,1,1]])
        )

    def test_sum_descendants(self):
        """sum_descendants should sum total descendants w/ each state"""
        id_, child = index_tree(self.t3)
        a = zeros((11,3)) + 99    #fill with junk
        bindings = bind_to_array(child, a)
        #load in leaf envs
        a[0] = a[1] = a[2] = a[7] = [0,1,0]
        a[3] = [1,0,0]
        a[6] = [0,0,1]
        sum_descendants(bindings)
        self.assertEqual(a, \
            array([[0,1,0],[0,1,0],[0,1,0],[1,0,0],[0,3,0],[1,0,0],\
            [0,0,1],[0,1,0],[1,3,0],[0,1,1],[1,4,1]])
        )

    def test_fitch_descendants(self):
        """fitch_descendants should assign states by fitch parsimony, ret. #"""
        id_, child = index_tree(self.t3)
        a = zeros((11,3)) + 99    #fill with junk
        bindings = bind_to_array(child, a)
        #load in leaf envs
        a[0] = a[1] = a[2] = a[7] = [0,1,0]
        a[3] = [1,0,0]
        a[6] = [0,0,1]
        changes = fitch_descendants(bindings)
        self.assertEqual(changes, 2)
        self.assertEqual(a, \
            array([[0,1,0],[0,1,0],[0,1,0],[1,0,0],[0,1,0],[1,0,0],\
            [0,0,1],[0,1,0],[1,1,0],[0,1,1],[0,1,0]])
        )

    def test_fitch_descendants_missing_data(self):
        """fitch_descendants should work with missing data"""
        #tree and envs for testing missing values
        t_str = '(((a:1,b:2):4,(c:3,d:1):2):1,(e:2,f:1):3);'
        env_str = """a   A
b   B
c   D
d   C
e   C
f   D"""
        t = DndParser(t_str, UniFracTreeNode)
        node_index, nodes = index_tree(t)
        env_counts = count_envs(env_str.split('\n'))
    
        count_array, unique_envs, env_to_index, node_to_index = \
            index_envs(env_counts, node_index)    

        branch_lengths = get_branch_lengths(node_index)
        #test just the AB pair
        ab_counts = count_array[:, 0:2]
        bindings = bind_to_array(nodes, ab_counts)
        changes = fitch_descendants(bindings, counter=FitchCounter)
        self.assertEqual(changes, 1)
        orig_result = ab_counts.copy()
        #check that the original Fitch counter gives the expected 
        #incorrect parsimony result
        changes = fitch_descendants(bindings, counter=FitchCounterDense)
        self.assertEqual(changes, 5)
        new_result = ab_counts.copy()
        #check that the two versions fill the array with the same values
        self.assertEqual(orig_result, new_result)

    def test_tip_distances(self):
        """tip_distances should set tips to correct distances."""
        t = self.t
        bl = self.branch_lengths.copy()[:,newaxis]
        bindings = bind_to_parent_array(t, bl)
        tips = []
        for n in t.traverse(self_before=False, self_after=True):
            if not n.Children:
                tips.append(n._leaf_index)
        tip_distances(bl, bindings, tips)
        self.assertEqual(bl, array([5,6,6,6,6,0,0,0,0])[:,newaxis])

    def test_permute_selected_rows(self):
        """permute_selected_rows should switch just the selected rows in a"""
        orig = reshape(arange(8),(4,2))
        new = orig.copy()
        fake_permutation = lambda a: range(a)[::-1] #reverse order
        permute_selected_rows([0,2], orig, new, fake_permutation)
        self.assertEqual(new,  array([[4,5],[2,3],[0,1],[6,7]]))
        #make sure we didn't change orig
        self.assertEqual(orig, reshape(arange(8), (4,2)))

    def test_prep_items_for_jackknife(self):
        """prep_items_for_jackknife should expand indices of repeated counts"""
        a = array([0,1,0,1,2,0,3])
        #          0 1 2 3 4 5 6
        result = prep_items_for_jackknife(a)
        exp = array([1,3,4,4,6,6,6])
        self.assertEqual(result, exp)

    def test_jackknife_bool(self):
        """jackknife_bool should make a vector with right number of nonzeros"""
        fake_permutation = lambda a: range(a)[::-1] #reverse order
        orig_vec = array([0,0,1,0,1,1,0,1,1])
        orig_items = flatnonzero(orig_vec)
        length = len(orig_vec)
        result = jackknife_bool(orig_items, 3, len(orig_vec), fake_permutation)
        self.assertEqual(result, array([0,0,0,0,0,1,0,1,1]))
        #returns the original if trying to take too many
        self.assertEqual(jackknife_bool(orig_items, 20, len(orig_vec)), \
            orig_vec)

    def test_jackknife_int(self):
        """jackknife_int should make a vector with right counts"""
        orig_vec = array([0,2,1,0,3,1])
        orig_items = array([1,1,2,4,4,4,5])
        #                   0 1 2 3 4 5 6
        fake_permutation = lambda a: a == 7 and array([4,6,3,1,2,6,5])
        result = jackknife_int(orig_items, 4, len(orig_vec), fake_permutation)
        self.assertEqual(result, array([0,1,0,0,2,1]))
        #returns the original if trying to take too many
        self.assertEqual(jackknife_int(orig_items, 20, len(orig_vec)), \
            orig_vec)
     
    def test_jackknife_array(self):
        """jackknife_array should make a new array with right counts"""

        orig_vec1 = array([0,2,2,3,1])
        orig_vec2 = array([2,2,1,2,2])
        test_array = array([orig_vec1, orig_vec2])

        # implement this, just doing by eye now
        #perm_fn = fake_permutation
        perm_fn = permutation

        #print "need to test with fake permutation!!"


        new_mat1 = jackknife_array(test_array, 1, axis=1, jackknife_f=jackknife_int, permutation_f=permutation)  
        self.assertEqual(new_mat1.sum(axis=0), [1,1,1,1,1])
        
        new_mat2 = jackknife_array(test_array, 2, axis=1, jackknife_f=jackknife_int, permutation_f=permutation)  
        self.assertEqual(new_mat2.sum(axis=0), [2,2,2,2,2])

        new_mat3 = jackknife_array(test_array, 2, axis=0, jackknife_f=jackknife_int, permutation_f=permutation)  
        self.assertEqual(new_mat3.sum(axis=1), [2,2])

        # test that you get orig mat back if too many
        self.assertEqual(jackknife_array(test_array, 20, axis=1), test_array)

    def test_unifrac(self):
        """unifrac should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        self.assertEqual(unifrac(bl, m[:,0], m[:,1]), 10/16.0)
        self.assertEqual(unifrac(bl, m[:,0], m[:,2]), 8/13.0)
        self.assertEqual(unifrac(bl, m[:,1], m[:,2]), 8/17.0)

    def test_unnormalized_unifrac(self):
        """unnormalized unifrac should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        self.assertEqual(unnormalized_unifrac(bl, m[:,0], m[:,1]), 10/17.)
        self.assertEqual(unnormalized_unifrac(bl, m[:,0], m[:,2]), 8/17.)
        self.assertEqual(unnormalized_unifrac(bl, m[:,1], m[:,2]), 8/17.)

    def test_PD(self):
        """PD should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        self.assertEqual(PD(bl, m[:,0]), 7)
        self.assertEqual(PD(bl, m[:,1]), 15)
        self.assertEqual(PD(bl, m[:,2]), 11)

    def test_G(self):
        """G should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        self.assertEqual(G(bl, m[:,0], m[:,0]), 0)
        self.assertEqual(G(bl, m[:,0], m[:,1]), 1/16.0)
        self.assertEqual(G(bl, m[:,1], m[:,0]), 9/16.0)

    def test_unnormalized_G(self):
        """unnormalized_G should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        self.assertEqual(unnormalized_G(bl, m[:,0], m[:,0]), 0/17.)
        self.assertEqual(unnormalized_G(bl, m[:,0], m[:,1]), 1/17.)
        self.assertEqual(unnormalized_G(bl, m[:,1], m[:,0]), 9/17.)

    def test_unifrac_matrix(self):
        """unifrac_matrix should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        result = unifrac_matrix(bl, m)
        self.assertEqual(result, array([[0, 10/16.,8/13.],[10/16.,0,8/17.],\
            [8/13.,8/17.,0]]))
        #should work if we tell it the measure is asymmetric
        result = unifrac_matrix(bl, m, is_symmetric=False)
        self.assertEqual(result, array([[0, 10/16.,8/13.],[10/16.,0,8/17.],\
            [8/13.,8/17.,0]]))
        #should work if the measure really is asymmetric
        result = unifrac_matrix(bl,m,metric=unnormalized_G,is_symmetric=False)
        self.assertEqual(result, array([[0, 1/17.,2/17.],[9/17.,0,6/17.],\
            [6/17.,2/17.,0]]))
        #should also match web site calculations
        envs = self.count_array
        bound_indices = bind_to_array(self.nodes, envs)
        bool_descendants(bound_indices)
        result = unifrac_matrix(bl, envs)
        exp = array([[0, 0.6250, 0.6154], [0.6250, 0, \
            0.4706], [0.6154, 0.4707, 0]])
        assert (abs(result - exp)).max() < 0.001

    def test_unifrac_vector(self):
        """unifrac_vector should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        result = unifrac_vector(bl, m)
        self.assertFloatEqual(result, array([10./17,6./17,7./17]))

    def test_PD_vector(self):
        """PD_vector should return correct results for model tree"""
        m = array([[1,0,1],[1,1,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1],[1,1,1],\
            [0,1,1],[1,1,1]])
        bl = self.branch_lengths
        result = PD_vector(bl, m)
        self.assertFloatEqual(result, array([7,15,11]))


    def test_weighted_unifrac_matrix(self):
        """weighted unifrac matrix should return correct results for model tree"""
        #should match web site calculations
        envs = self.count_array
        bound_indices = bind_to_array(self.nodes, envs)
        sum_descendants(bound_indices)
        bl = self.branch_lengths
        tip_indices = [n._leaf_index for n in self.t.tips()]
        result = weighted_unifrac_matrix(bl, envs, tip_indices)
        exp = array([[0, 9.1, 4.5], [9.1, 0, \
            6.4], [4.5, 6.4, 0]])
        assert (abs(result - exp)).max() < 0.001
        #should work with branch length corrections
        td = bl.copy()[:,newaxis]
        tip_bindings = bind_to_parent_array(self.t, td)
        tips = [n._leaf_index for n in self.t.tips()]
        tip_distances(td, tip_bindings, tips)
        result = weighted_unifrac_matrix(bl, envs, tip_indices, bl_correct=True,
            tip_distances=td)
        exp = array([[0, 9.1/11.5, 4.5/(10.5+1./3)], [9.1/11.5, 0, \
            6.4/(11+1./3)], [4.5/(10.5+1./3), 6.4/(11+1./3), 0]])
        assert (abs(result - exp)).max() < 0.001

    def test_weighted_unifrac_vector(self):
        """weighted_unifrac_vector should return correct results for model tree"""
        envs = self.count_array
        bound_indices = bind_to_array(self.nodes, envs)
        sum_descendants(bound_indices)
        bl = self.branch_lengths
        tip_indices = [n._leaf_index for n in self.t.tips()]
        result = weighted_unifrac_vector(bl, envs, tip_indices)
        self.assertFloatEqual(result[0], sum([
            abs(1./2 - 2./8)*1,
            abs(1./2 - 1./8)*2,
            abs(0 - 1./8)*3,
            abs(0 - 3./8)*1,
            abs(0 - 1./8)*1,
            abs(0 - 4./8)*2,
            abs(2./2 - 3./8)*4,
            abs(0. - 5./8)*3.]))

        self.assertFloatEqual(result[1], sum([
            abs(0-.6)*1,
            abs(.2-.2)*2,
            abs(.2-0)*3,
            abs(.6-0)*1,
            abs(0-.2)*1,
            abs(.6-.2)*2,
            abs(.2-.8)*4,
            abs(.8-.2)*3]))

        self.assertFloatEqual(result[2], sum([
            abs(2./3-1./7)*1,
            abs(0-2./7)*2,
            abs(0-1./7)*3,
            abs(0-3./7)*1,
            abs(1./3-0)*1,
            abs(1./3-3./7)*2,
            abs(2./3-3./7)*4,
            abs(1./3-4./7)*3]))


if __name__ == '__main__':    #run if called from command-line
    main()
