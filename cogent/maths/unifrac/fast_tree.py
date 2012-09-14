#!/usr/bin/env python
""" Fast tree and support functions for fast implementation of UniFrac """

from numpy import (logical_and, logical_or, sum, take, nonzero, repeat, 
    array, concatenate, zeros, put, transpose, flatnonzero, newaxis,
    logical_xor, logical_not)
from numpy.random import permutation
from cogent.core.tree import PhyloNode

__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Micah Hamady", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight, Micah Hamady"
__email__ = "rob@spot.colorado.edu, hamady@colorado.edu"
__status__ = "Prototype"

#bind reduce to local variables for speed in inner loops
lar = logical_and.reduce
lor = logical_or.reduce

class UniFracTreeNode(PhyloNode):
    """Slightly extended PhyloNode treenode for use with UniFrac
    
    Can expect Length and Name (= Name) to be set by DndParser.
    """
    def __nonzero__(self):
        """Returns True if self.Children."""
        return bool(self.Children)

def index_tree(t):
    """Returns tuple containing {node_id:node}, [node_id,first_child,last_child]

    Indexes nodes in-place as n._leaf_index.

    Algorithm is as follows:
    for each node in post-order traversal over tree:
        if the node has children:
            set an index on each child
            for each child with children:
                add the child and its start and end tips to the result
    """
    id_index = {}    #needs to be dict, not list, b/c adding out of order
    child_index = []
    curr_index = 0
    for n in t.traverse(self_before=False, self_after=True):
        for c in n.Children:
            c._leaf_index = curr_index
            id_index[curr_index] = c
            curr_index += 1
            if c:    #c has children itself, so need to add to result
                child_index.append((c._leaf_index, c.Children[0]._leaf_index,\
                    c.Children[-1]._leaf_index))
    #handle root, which should be t itself
    t._leaf_index = curr_index
    id_index[curr_index] = t
    #only want to add to the child_index if t has children...
    if t.Children:
        child_index.append((t._leaf_index, t.Children[0]._leaf_index,\
            t.Children[-1]._leaf_index))
    return id_index, child_index


def count_envs(lines, ignore_chars=0):
    """Reads env counts from lines. Returns dict of {name:{env:count}}.
    
    Assumes all name-env mappings are unique -- will overwrite counts with the
    last observed count rather than adding.
    """
    result = {}
    for line in lines:
        fields = line.split()
        #skip if we don't have at least label and field
        if len(fields) < 2:
            continue
        name, env = fields[:2]
        if ignore_chars:
            env = env[ignore_chars:]
        if len(fields) > 2:
            count = int(fields[2])
        else:
            count = 1
        if name not in result:
            result[name] = {}
        result[name][env] = count
    return result

def sum_env_dict(envs):
    """Sums counts from the data structure produced by count_envs."""
    return sum([sum(env.values()) for env in envs.values()])

def get_unique_envs(envs):
    """extract all unique envs from envs dict"""
    result = set()
    for v in envs.values():
        result.update(v.keys())
    #sort envs for convenience in testing and display
    return sorted(result), len(result)


def index_envs(env_counts, tree_index, array_constructor=int):
    """Returns array of taxon x env with counts of the taxon in each env.

    env_counts should be the output of count_envs(lines).
    tree_index should be the id_index of index_tree(t).
    array_constructor is int by default (may need to change to float later
        to handle microarray data).
    """
    num_nodes = len(tree_index)
    unique_envs, num_envs = get_unique_envs(env_counts)
    env_to_index = dict([(e, i) for i, e in enumerate(unique_envs)])
    result = zeros((num_nodes, num_envs), array_constructor)
    #figure out taxon label to index map
    node_to_index = {}
    for i, node in tree_index.items():
        if node.Name is not None:
            node_to_index[node.Name] = i
    #walk over env_counts, adding correct slots in array
    for name in env_counts:
        curr_row_index = node_to_index[name]
        for env, count in env_counts[name].items():
            result[curr_row_index, env_to_index[env]] = count
    #return all the data structures we created; will be useful for other tasks
    return result, unique_envs, env_to_index, node_to_index

def get_branch_lengths(tree_index):
    """Returns array of branch lengths, in tree index order."""
    result = zeros(len(tree_index), float)
    for i, node in tree_index.items():
        try:
            if node.Length is not None:
                result[i] = node.Length
        except AttributeError:
            pass
    return result

def bind_to_array(tree_index, a):
    """Binds tree_index to array a, returning result in list.

    Takes as input list of (node, first_child, last_child)
    returns list of (node_row, child_rows) such that node_row points to the
    row of a that corresponds to the current node, and child_rows points to the
    row or rows of a that correspond to the direct children of the current
    node.

    Order is assumed to be traversal order, i.e. for the typical case of
    postorder traversal iterating over the items in the result and consolidating
    each time should give the same result as postorder traversal of the original
    tree. Should also be able to modify for preorder traversal.
    """
    #note: range ends with end+1, not end, b/c end is included
    return [(a[node], a[start:end+1]) for node, start, end in tree_index]

def bind_to_parent_array(t, a):
    """Binds tree to array a, returning result in list.

    Takes as input tree t with _leaf_index set.

    Returns list of (node_row, parent_row such that node_row points to the
    row of a that corresponds to the current row, and parent_row points to
    the row of the parent.

    Order will be preorder traversal, i.e. for propagating attributes from
    the root to the tip. 

    Typical usage of this function is to set up an array structure for many
    preorder traversals on the same tree, especially where you plan to change
    the data between traversals.
    """
    result = []
    for n in t.traverse(self_before=True, self_after=False):
        if n is not t:
            result.append([a[n._leaf_index], a[n.Parent._leaf_index]])
    return result

def _is_parent_empty(parent_children):
    """Returns True if the first element in a (parent,children) tuple is empty.

    This is used by delete_empty_parents to figure out which elements to filter
    out.
    """
    return bool(parent_children[0].sum())


def delete_empty_parents(bound_indices):
    """Deletes from list of (parent, children) bound indices empty parents.

    Expects as input the output of bind_to_array. Returns copy rather than 
    acting in-place because deletions from long lists are expensive.

    This has the effect of pruning the tree, but by just skipping over the 
    parents who have no children rather than altering memory. This is expected 
    to be faster than trimming the array because it avoids copy operations.

    For pairwise environment bootstrapping or jackknifing, run this after 
    running bool_descendants (or similar) to delete parents that only have 
    offspring that are in other environments. Note that this does _not_ 
    collapse nodes so that you might have long "stalks" with many serial nodes.
    It might be worth testing whether collapsing these stalks provides time 
    savings.
    """
    return filter(_is_parent_empty, bound_indices)

def traverse_reduce(bound_indices, f):
    """Applies a[i] = f(a[j:k]) over list of [(a[i], a[j:k])].

    If list is in traversal order, has same effect as consolidating the
    function over the tree, only much faster.

    Note that f(a[j:k]) must return an object that can be broadcast to the
    same shape as a[i], e.g. summing a 2D array to get a vector.
    """
    for i, s in bound_indices:
        i[:] = f(s, 0)

def bool_descendants(bound_indices):
    """For each internal node, sets col to True if any descendant is True."""
    traverse_reduce(bound_indices, lor)

def zero_branches_past_roots(bound_indices, sums):
    """Zeroes out internal nodes that are roots of each subtree."""
    for i, ignore in bound_indices:
        i *= (i != sums)

def sum_descendants(bound_indices):
    """For each internal node, sets col to sum of values in descendants."""
    traverse_reduce(bound_indices, sum)

class FitchCounterDense(object):
    """Returns parsimony result for set of child states, counting changes.
    
    WARNING: this version assumes that all tips are assigned to at least
    one env, and produces incorrect parsimony counts if this is not the
    case.
    """
    def __init__(self):
        """Returns new FitchCounter, with Changes = 0."""
        self.Changes = 0
        
    def __call__(self, a, ignored):
        """Returns intersection(a), or, if zero, union(a)."""
        result = lar(a)
        if not result.any():
            result = lor(a)
            self.Changes += 1
        return result

class FitchCounter(object):
    """Returns parsimony result for set of child states, counting changes.
    
    This version is slower but is robust to the case where some tips are
    missing envs.

    WARNING: logical_and.reduce(), if called on an empty array, will return
    all True if there are no values (I can only assume that it returns
    True because it is reporting that there are no values that return False:
    this isn't what I expected). Hence the code to explicitly trap this
    case based on the shape parameter.
    """
    def __init__(self):
        """Returns new FitchCounter, with Changes = 0."""
        self.Changes = 0
        
    def __call__(self, a, ignored):
        """Returns intersection(a), or, if zero, union(a)."""
        nonzero_rows = a[a.sum(1).nonzero()]
        if len(nonzero_rows):
            result = lar(nonzero_rows)
        else:
            result = zeros(nonzero_rows.shape[-1], bool)
        if not result.any():
            if nonzero_rows.any():
                result = lor(nonzero_rows)
                self.Changes += 1
        return result


def fitch_descendants(bound_indices, counter=FitchCounter):
    """Sets each internal node to Fitch parsimony assignment, returns # changes."""
    f = counter()
    traverse_reduce(bound_indices, f.__call__)
    return f.Changes

def tip_distances(a, bound_indices, tip_indices):
    """Sets each tip to its distance from the root."""
    for i, s in bound_indices:
        i += s
    mask = zeros(len(a))
    put(mask, tip_indices, 1)
    a *= mask[:,newaxis]

def permute_selected_rows(rows, orig, new, permutation_f=permutation):
    """Takes selected rows from orig, inserts into new in permuted order.

    NOTE: the traditional UniFrac permutation test, the P test, etc. shuffle 
    the envs, i.e. they preserve all the correlations of seqs between envs.
    This function can also be used to shuffle each env (i.e. column)
    individually by applying it to column slices of orig and new. This 
    latter method provides a potentially less biased but less conservative
    test.
    """
    shuffled = take(rows, permutation_f(len(rows)))
    for r, s in zip(rows, shuffled):
        new[s] = orig[r]

def prep_items_for_jackknife(col):
    """Takes column of a, returns vector with multicopy states unpacked.
    
    e.g. if index 3 has value 4, there will be 4 copies of index 3 in result.
    """
    nz = flatnonzero(col)
    result = [repeat(array((i,)), col[i]) for i in nz]
    return concatenate(result)

def jackknife_bool(orig_items, n, length, permutation_f=permutation):
    """Jackknifes vector of items so that only n remain.
    
    orig = flatnonzero(vec)
    length = len(vec)

    Returns all items if requested sample is larger than number of items.
    """
    permuted = take(orig_items, permutation_f(len(orig_items))[:n])
    result = zeros(length)
    put(result, permuted, 1)
    return result

def jackknife_int(orig_items, n, length, permutation_f=permutation):
    """Jackknifes new vector from vector of orig items.
    
    Returns all items if reqested sample is larger than number of items.
    """
    result = zeros(length)
    permuted = take(orig_items, permutation_f(len(orig_items))[:n])
    for p in permuted:
        result[p] += 1
    return result

def jackknife_array(mat, num_keep, axis=1, jackknife_f=jackknife_int, 
    permutation_f=permutation):
    """ Jackknife array along specified axis, keeping specified num_keep"""
    cur_mat = mat
    if axis:
        cur_mat = mat.T
    num_r, num_c = cur_mat.shape 
    jack_mat = zeros((num_r, num_c))
    for row_ix in range(num_r):
        in_prepped_array = prep_items_for_jackknife(cur_mat[row_ix,:])
        jack_mat[row_ix,:] = jackknife_f(orig_items=in_prepped_array, 
            n=num_keep, length=num_c, permutation_f=permutation_f)
    if axis:
        jack_mat = jack_mat.T
    return jack_mat


def unifrac(branch_lengths, i, j):
    """Calculates unifrac(i,j) from branch lengths and cols i and j of m.
    
    This is the original, unweighted UniFrac metric.

    branch_lengths should be row vector, same length as # nodes in tree.
    
    i and j should be slices of states from m, same length as # nodes in
    tree. Slicing m (e.g. m[:,i]) returns a vector in the right format; note
    that it should be a row vector (the default), not a column vector.
    """
    return 1 - ((branch_lengths*logical_and(i,j)).sum()/\
        (branch_lengths*logical_or(i,j)).sum())

def unnormalized_unifrac(branch_lengths, i, j):
    """UniFrac, but omits normalization for frac of tree covered."""
    return (branch_lengths*logical_xor(i,j)).sum()/branch_lengths.sum()

def G(branch_lengths, i, j):
    """Calculates G(i,j) from branch lengths and cols i,j of m.

    This calculates fraction gain in branch length in i with respect to i+j,
    i.e. normalized for the parts of the tree that i and j cover.
    
    Note: G is the metric that we also call "asymmetric unifrac".
    """
    return (branch_lengths*logical_and(i, logical_not(j))).sum()/\
        (branch_lengths*logical_or(i,j)).sum()

def PD(branch_lengths, i):
    """Calculate PD(i) from branch lengths and col i of m.

    Calculates raw amount of branch length leading to tips in i, including
    branch length from the root.
    """
    return (branch_lengths * i.astype(bool)).sum()

def unnormalized_G(branch_lengths, i, j):
    """Calculates G(i,j) from branch length and cols i,j of m.

    This calculates the fraction gain in branch length of i with respect to j,
    divided by all the branch length in the tree.
    """
    return (branch_lengths*logical_and(i, logical_not(j))).sum()/\
        branch_lengths.sum()

def unifrac_matrix(branch_lengths, m, metric=unifrac, is_symmetric=True):
    """Calculates unifrac(i,j) for all i,j in m.
    
    branch_lengths is the array of branch lengths.
    
    m is 2D array: rows are taxa, states are columns. Assumes that ancestral
    states have already been calculated (either by logical_or or Fitch).

    metric: metric to use for combining each pair of columns i and j. Default
    is unifrac.

    is_symmetric indicates whether the metric is symmetric. Default is True.
    """
    num_cols = m.shape[-1]
    cols = [m[:,i] for i in range(num_cols)]
    result = zeros((num_cols,num_cols), float)
    if is_symmetric:
        #only calc half matrix and transpose
        for i in range(1, num_cols):
            first_col = cols[i]
            row_result = []
            for j in range(i):
                second_col = cols[j]
                row_result.append(metric(branch_lengths, first_col, second_col))
            result[i,:j+1] = row_result
        #note: can't use += because shared memory between a and transpose(a)
        result = result + transpose(result)
    else:
        #calc full matrix, incl. diagonal (which is probably 0...)
        for i in range(num_cols):
            first_col = cols[i]
            result[i] = [metric(branch_lengths, first_col, cols[j]) for \
                j in range(num_cols)]
    return result
    
def unifrac_one_sample(one_sample_idx, branch_lengths, m, metric=unifrac):
    """Calculates unifrac(one_sample_idx,j) for all environments j in m.
    
    branch_lengths is the array of branch lengths.
    
    m is 2D count array: rows are taxa (corresponding to branch_lenths),
    samples/states/envs are columns. 
    Assumes that ancestral 
    states have already been calculated (either by logical_or or Fitch).

    metric: metric to use for when comparing two environments. Default
    is unifrac. must be called like: 
    metric(branch_lengths, env1_counts, env2counts)
    
    returns a numpy 1d array
    if asymmetric metric, returns metric(one_sample, other), usually a row in
    the mtx returned by unifrac_matrix
    """
    num_cols = m.shape[-1]
    cols = [m[:,i] for i in range(num_cols)]
    # result = zeros((num_cols), float)
        
    first_col = cols[one_sample_idx]
    # better to do loop into preallocated numpy array here?
    result = array([metric(branch_lengths, first_col, cols[j]) for \
        j in range(num_cols)],'float')
    return result

def env_unique_fraction(branch_lengths, m):
    """ Calculates unique branch length for each env. 

    Returns unique branch len and unique fraction 
    """ 
    total_bl = branch_lengths.sum()
    if total_bl <= 0:
        raise ValueError, "total branch length in tree must be > 0"

    n_rows_nodes, n_col_envs = m.shape
    env_bl_sums = zeros(n_col_envs)

    cols = [m[:, i] for i in range(n_col_envs)]
    col_sum = m.sum(1)
    env_bl_sums = zeros(n_col_envs)
    for env_ix, f in enumerate(cols):
        sing = (f == col_sum)
        # have  to mask zeros
        put(sing, nonzero(f == 0), 0)
        env_bl_sums[env_ix] = (sing * branch_lengths).sum()
    
    return env_bl_sums, env_bl_sums/total_bl 

def unifrac_vector(branch_lengths, m, metric=unifrac):
    """Calculates unifrac(i, others) for each column i of m.

    Parameters as for unifrac_matrix. Use this when you want to calculate
    UniFrac or G of each state against the rest of the states, rather than
    of each state against each other state.
    """
    num_cols = m.shape[-1]
    cols = [m[:, i] for i in range(num_cols)]
    col_sum = m.sum(1)
    return array([metric(branch_lengths, col, col_sum-col) for col in cols])

def PD_vector(branch_lengths, m, metric=PD):
    """Calculates metric(i) for each column i of m.

    Parameters as for unifrac_matrix. Use this when you want to calculate
    PD or some other alpha diversity metric that depends solely on the branches
    within each state, rather than calculations that compare each state against 
    each other state.
    """
    return array([metric(branch_lengths, col) for col in m.T])

def _weighted_unifrac(branch_lengths, i, j, i_sum, j_sum):
    """Calculates weighted unifrac(i,j) from branch lengths and cols i,j of m.
    """
    return (branch_lengths * abs((i/float(i_sum))-(j/float(j_sum)))).sum()

def _branch_correct(tip_distances, i, j, i_sum, j_sum):
    """Calculates weighted unifrac branch length correction.
    
    tip_distances  must be 0 except for tips.
    """
    result = tip_distances.ravel()*((i/float(i_sum))+(j/float(j_sum)))
    return result.sum()

def weighted_unifrac(branch_lengths, i, j, tip_indices, \
    unifrac_f=_weighted_unifrac):
    """Returns weighted unifrac(i,j) from branch lengths and cols i,j of m.

    Must pass in tip indices to calculate sums.

    Note: this calculation is not used in practice because it has to recalc.
    the sum each time. More efficient to calculate the sum first and pass it
    into _weighted_unifrac directly, as weighted_unifrac_matrix does.
    """
    i_sum = (take(i, tip_indices)).sum()
    j_sum = (take(j, tip_indices)).sum()
    return unifrac_f(branch_lengths, i, j, i_sum, j_sum)

def weighted_unifrac_matrix(branch_lengths, m, tip_indices, bl_correct=False,
        tip_distances=None, unifrac_f=_weighted_unifrac):
    """Calculates weighted_unifrac(i,j) for all i,j in m.

    Requires tip_indices for calculating sums, etc.
    bl_correct: if True (default: False), applies branch length correction.
        tip_distances is required for normalization for weighted unifrac.
    """
    num_cols = m.shape[-1]
    cols = [m[:,i] for i in range(num_cols)] #note that these will be row vecs
    sums = [take(m[:,i], tip_indices).sum() for i in range(num_cols)]
    result = zeros((num_cols,num_cols),float)
    for i in range(1, num_cols):
        i_sum = sums[i]
        first_col = cols[i]
        row_result = []
        for j in range(i):
            second_col = cols[j]
            j_sum = sums[j]
            curr = unifrac_f(branch_lengths, first_col, \
                second_col, i_sum, j_sum)
            if bl_correct:
                curr /= _branch_correct(tip_distances, first_col, \
                    second_col, i_sum, j_sum)
            row_result.append(curr)
        result[i,:j+1] = row_result
    result = result + transpose(result)
    return result
    
def weighted_one_sample(one_sample_idx, branch_lengths, m, tip_indices, 
    bl_correct=False, tip_distances=None, unifrac_f=_weighted_unifrac):
    """Calculates weighted_unifrac(one_sample_idx,j) for all environments j in m

    Requires tip_indices for calculating sums, etc.
    bl_correct: if True (default: False), applies branch length correction.
        tip_distances is required for normalization for weighted unifrac.
    """
    num_cols = m.shape[-1]
    cols = [m[:,i] for i in range(num_cols)] #note that these will be row vecs
    sums = [take(m[:,i], tip_indices).sum() for i in range(num_cols)]
    result = zeros((num_cols),float)

    i_sum = sums[one_sample_idx]
    first_col = cols[one_sample_idx]
    row_result = []
    for j in range(num_cols):
        second_col = cols[j]
        j_sum = sums[j]
        curr = unifrac_f(branch_lengths, first_col, \
            second_col, i_sum, j_sum)
        if bl_correct:
            curr /= _branch_correct(tip_distances, first_col, \
                second_col, i_sum, j_sum)
        result[j] = curr
    return result

def weighted_unifrac_vector(branch_lengths, m, tip_indices, bl_correct=False,
    tip_distances=None, unifrac_f=_weighted_unifrac):
    """Calculates weighted_unifrac(i,rest) for i in m.

    Requires tip_indices for calculating sums, etc.
    bl_correct: if True (default: False), applies branch length correction.
        tip_distances is required for normalization for weighted unifrac.
    """
    num_cols = m.shape[-1]
    cols = [m[:,i] for i in range(num_cols)]
    sums = [take(m[:,i], tip_indices).sum() for i in range(num_cols)]
    sum_of_cols = m.sum(1)
    sum_of_sums = sum(sums)
    result = []
    for i, col in enumerate(cols):
        i_sum = sums[i]
        i_col = cols[i]
        rest_col = sum_of_cols - i_col
        rest_sum = sum_of_sums - i_sum
        curr = unifrac_f(branch_lengths, i_col, 
            rest_col, i_sum, rest_sum)
        if bl_correct:
            curr /= _branch_correct(tip_distances, i_col,
            rest_col, i_sum, rest_sum)
        result.append(curr)
    return array(result)
