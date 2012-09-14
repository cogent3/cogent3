#!/usr/bin/env python
from tree import RangeNode, balanced_breakpoints
from cogent.core.usage import DnaPairs
from usage import Counts
from random import choice
from numpy.linalg import det as determinant, inv as inverse
from numpy import sqrt, newaxis as NewAxis, exp, dot, zeros, ravel, array, \
                  float64, max, min, average, any, pi
from cogent.util.array import without_diag
from cogent.maths.svd import three_item_combos, two_item_combos
from cogent.maths.stats.test import std

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

def tree_threeway_counts(tree, lca_depths, alphabet=DnaPairs, attr='Sequence'):
    """From tree and array of lca_depths, returns n*n*n array of Count objects.

    n is number of leaves.

    lca_depths: array (leaf * leaf) of depths of last common ancestor.
    alphabet: pair alphabet for input sequences.

    Returns dict containing counts for (i, j, k) and (j, i, k) where k is the
    outgroup of the three sequences. Will pick an arbitrary node to be the
    outgroup if there is a polytomy.
    
    Note: Leaves of tree must have sequences already assigned.
    """
    outgroup_last = tree.outgroupLast
    leaves = list(tree.traverse())
    result = {}
    for first, second, third in three_item_combos(leaves):
        new_first, new_second, new_third = outgroup_last(first, second, third)
        #get the sequence from each node
        seq_1 = getattr(new_first, attr)
        seq_2 = getattr(new_second, attr)
        seq_3 = getattr(new_third, attr)
        
        result[(new_first.Id, new_second.Id, new_third.Id)] = \
            Counts.fromTriple(seq_1, seq_2, seq_3, alphabet)
        #don't forget to do counts from both  the non-outgroups
        result[(new_second.Id, new_first.Id, new_third.Id)] = \
            Counts.fromTriple(seq_2, seq_1, seq_3, alphabet)
    return result

def dna_count_cleaner(counts):
    """Cleans DNA counts to just the 4-letter alphabet."""
    return Counts(counts._data[:4,:4], DnaPairs)

def tree_threeway_counts_sample(tree, lca_depths, alphabet=DnaPairs, \
    attr='Sequence', n=1000, check_rates=True, clean_f=None):
    """Like tree_threeway_counts, but takes random sample (w/o replacement)."""
    leaves = list(tree.traverse())
    num_leaves = len(leaves)
    #do normal threeway counts if number of triples < n
    num_triples = num_leaves * (num_leaves - 1) * (num_leaves-2) / 3
    if num_triples < n:
        counts = tree_threeway_counts(tree, lca_depths, alphabet, attr)
        if clean_f:
            result = {}
            for k, v in counts.items():
                result[k] = clean_f(v)
            return result
        else:
            return counts
    #if we got here, need to sample
    outgroup_last = tree.outgroupLast
    i = 0
    seen = {}
    result = {}
    while i < n and len(seen) < num_triples:
        #bail out if same node picked twice, or if resampling same combo
        curr = choice(leaves), choice(leaves), choice(leaves)
        ids = tuple([c.Id for c in curr])
        if len(dict.fromkeys(ids)) < len(curr):     #picked same thing twice
            continue
        if curr in seen:
            continue
        first, second, third = curr
        new_first, new_second, new_third = outgroup_last(first, second, third)
        seq_1 = getattr(new_first, attr)
        seq_2 = getattr(new_second, attr)
        seq_3 = getattr(new_third, attr)

        counts = Counts.fromTriple(seq_1, seq_2, seq_3, alphabet)
        if clean_f:
            counts = clean_f(counts)
        key = (new_first.Id, new_second.Id, new_third.Id)
        #check rates if we need to
        if check_rates:
            try:
                #skip probs with zero rows
                if not min(max(counts._data,1)):
                    continue
                probs = counts.toProbs()
                rates = probs.toRates()
            except (ZeroDivisionError, OverflowError, ValueError, \
                FloatingPointError):
                continue
            result[key] = counts
        i += 1
    return result
        

def tree_twoway_counts(tree, alphabet=DnaPairs, average=True, attr='Sequence'):
    """From tree, return dict of Count objects.

    Note: if average is True, only has counts in m[i,j] or m[j,i], not both.
    """
    leaves = list(tree.traverse())
    result = {}
    if average:
        #return symmetric matrix
        for first, second in two_item_combos(leaves):
            seq_1 = getattr(first, attr)
            seq_2 = getattr(second, attr)
            result[(first.Id, second.Id)] = \
                Counts.fromPair(seq_1, seq_2, alphabet)
    else:
        for first, second in two_item_combos(leaves):
            seq_1 = getattr(first, attr)
            seq_2 = getattr(second, attr)
            result[(first.Id, second.Id)] = \
                Counts.fromPair(seq_1, seq_2, alphabet,False)
            result[(second.Id, first.Id)] = \
                Counts.fromPair(seq_2, seq_1, alphabet,False)
    return result


def counts_to_probs(count_dict):
    """Converts counts to probs, omitting problem cases."""
    result = {}
    for key, val in count_dict.items():
        #check for zero rows
        if not min(max(val._data,1)):
            continue
        try:
            p = val.toProbs()
            #the following detects nan from divide by zero for empty rows
            #this works because nan doesn't compare equal to itself
            first_col = p._data[:,0]
            if any(first_col != first_col):
                raise ZeroDivisionError
            #if we got here, everything was OK
            result[key] = val.toProbs()
        except (ZeroDivisionError, OverflowError,ValueError,FloatingPointError):
            #errors are platform-dependent and arise when a row is zero
            #(i.e. if a character doesn't appear in the sequence).
            pass
    return result

def probs_to_rates(prob_dict, fix_f=None, normalize=False, ftol=0.01):
    """Converts probs to rates, omitting problem cases (but not neg off-diags).
    
    ftol is in log10 units, e.g. ftol of 1 means within one order of magnitude.
    Default = 0.01.
    
    NOTE: fix_f should be an unbound method of Rates, or other function that
    expects a rate as a single parameter.
    """
    result = {}
    seen = {}
    for key, val in prob_dict.items():
        try:
            rate = val.toRates(normalize)
            if fix_f:
                rate = fix_f(rate)
            rounded = tuple(map(round, rate._data.flatten()/ftol))
            if rounded in seen:
                continue
            else:
                seen[rounded] = 1
            #check for zero rows
            if not min(max(rate._data,1)):
                continue
            #check for zero cols
            if not min(max(rate._data)):
                continue
            if (not rate.isSignificantlyComplex()): # and rate.isValid():
                result[key] = rate
        except (ZeroDivisionError, OverflowError,ValueError,FloatingPointError):
            #errors are platform-dependent
            pass
    return result
       
def tree_threeway_rates(tree, lca_depths, alphabet=DnaPairs, fix_f=None, \
    normalize=False, without_diag=False):
    """Generates matrix of all valid three-way rate matrices from a tree.
    
    Dimensions of matrix are (num_leaves, num_leaves, num_leaves, flat_rate),
    i.e. m[0][1][2] gives you the matrix for leaves 0,1,2.
    In general, expect to get matrices out of the tree by slicing rather than
    by indexing, since many rates will be empty due to inference problems.
    """
    leaves = len(list(tree.traverse()))
    counts = tree_threeway_counts(tree, lca_depths, alphabet)
    rates = probs_to_rates(counts_to_probs(counts), fix_f, normalize)
    return threeway_rates_to_array(rates, leaves, alphabet, without_diag)

def tree_twoway_rates(tree, alphabet=DnaPairs, average=False, fix_f=None, \
    normalize=False, without_diag=False):
    """Generates all valid two-way rate matrices from a tree.
    
    Dimensions of matrix are (num_leaves, num_leaves, num_leaves, flat_rate),
    i.e. m[0][1][2] gives you the matrix for leaves 0,1,2.
    In general, expect to get matrices out of the tree by slicing rather than
    by indexing, since many rates will be empty due to inference problems.
    """
    leaves = len(list(tree.traverse()))
    counts = tree_twoway_counts(tree, alphabet, average)
    probs = counts_to_probs(counts)
    r = probs_to_rates(probs, fix_f, normalize)
    return twoway_rates_to_array(r, leaves, alphabet, without_diag)

def twoway_rates_to_array(rates, num_seqs, alphabet=DnaPairs, \
    without_diag=True):
    """Fills an array with flat twoway_rates arrays."""
    if without_diag:
        result = zeros((num_seqs, num_seqs, len(alphabet)- \
            len(alphabet.SubEnumerations[0])),
            float64)
    else:
        result = zeros((num_seqs, num_seqs, len(alphabet)), float64)
    return rates_to_array(rates, result, without_diag)

def threeway_rates_to_array(rates, num_seqs, alphabet=DnaPairs, \
    without_diag=True):
    """Fills an array with flat threeway_rates arrays."""
    if without_diag:
        result = zeros((num_seqs, num_seqs, num_seqs, \
            len(alphabet) - len(alphabet.SubEnumerations[0])), float64)
    else:
        result = zeros((num_seqs, num_seqs, num_seqs, len(alphabet)), float64)
    return rates_to_array(rates, result, without_diag)

def rates_to_array(rates, to_fill, without_diagonal=False):
    """Fills rates into a pre-existing array object.

    Assumes that all keys in rates are valid indices into to_fill, and that the
    last dimension of to_fill is the same size as the flattened values in rates.
    
    If without_diagonal is True (False by default), removes the diagonals from 
    the data before placing in the array. Note that we can't call this
    without_diag or we collide with the function we want to call to strip
    the diagonal.
    
    WARNING: size of to_fill array must be adjusted to be the right size as the
    inputs, i.e. last dimension same as flat array with/without diagonals.
    """
    if without_diagonal:
        for key, val in rates.items():
            to_fill[key] = ravel(without_diag(val._data))
    else:
        for key, val in rates.items():
            to_fill[key] = ravel(val._data)
    return to_fill

def multivariate_normal_prob(x, cov, mean=None):
    """Returns multivariate normal probability density of vector x.
    
    Formula: http://www.riskglossary.com/link/joint_normal_distribution.htm
    """
    if mean is None:
        mean = zeros(len(x))
    diff_row = x-mean
    diff_col = diff_row[:,NewAxis]
    numerator = exp(-0.5 * dot(dot(diff_row, inverse(cov)), diff_col))
    denominator = sqrt((2*pi)**(len(x)) * determinant(cov))
    return numerator/denominator

######### WARNING: TESTS END HERE! ###########################


def classify_tree(tree, comparison_sets):
    """Takes a tree and a dict of label -> distributions. Returns node -> label."""
    pass
    
def balanced_two_q_tree(n, length=0.05, seq_length=100, change_depth=2,\
        perturb=True, both=False):
    """Makes single random tree with specified nodes, branch and seq length.
    
    Two random rate matrices are assigned.
    """
    t = BalancedTree(n)
    t.assignLengths(length)
    if perturb:
        t.pertubrAttr('Length', xxx)
        raise   #NOT FINISHED
    q = Rates.random()
    scale_trace(q)
    t.Q = q
    curr_node = t
    for i in range(change_depth):
        curr_node = choice(curr_node)
    q2 = random_q_matrix()
    scale_trace(q2)
    curr_node.Q = q2
    t.assignQ()
    t.assignP()
    t.Sequence = rand_rna(seq_length)
    t.assignSeqs()
    if both:
        return t, curr_node
    else:
        return t

def balanced_multi_q_tree(n, length=0.05, seq_length=100,perturb=True):
    """Makes single random tree with specified nodes, branch and seq length.
    
    Every node gets its own random rate matrix.
    """
    t = BalancedTree(n)
    t.assignLengths({},default=length)
    if perturb:
        perturb_lengths(t)
    for node in t.traverse():
        q = random_q_matrix()
        scale_trace(q)
        node.Q = q
    t.assignQ()
    t.assignP()
    t.Sequence = rand_rna(seq_length)
    t.assignSeqs()
    return t

def get_matrix_stats(qs, true_q=None):
    """Returns list of matrix stats. Expects flat q matrices as input."""
    #set up flat q matrices
    flat_qs = flatten_q_matrices(qs)
    for q in qs:
        scale_trace(qs)
    flat_scaled_qs = flatten_q_matrices(qs)
    #get covariance and correl matrices
    covar = cov(flat_qs)
    correl = corrcoef(flat_qs)
    norm_covar = euclidean_norm(covar)
    norm_correl = euclidean_norm(correl)
    covar_e = list(eigenvalues(covar).real)
    covar_e.sort()
    covar_e.reverse()
    correl_e = list(eigenvalues(correl).real)
    correl_e.sort()
    correl_e.reverse()
    two_best_covar = ratio_two_best(covar_e)
    two_best_correl = ratio_two_best(correl_e)
    frac_best_covar = ratio_best_to_sum(covar_e)
    frac_best_correl = ratio_best_to_sum(correl_e)
    weiss_stat = weiss(covar_e)
    #get summary stats from scaled qs
    distances = dists_from_mean(flat_scaled_qs)
    mean_dist = average(distances)
    var_dist = var(distances)
    if true_q:
        mean_q = mean(flat_scaled_qs)
        q_error = euclidean_distance(mean_q, true_q)
    else:
        q_error = 0
    norm_scaled_covar = euclidean_norm(cov(flat_scaled_qs))

    result = [norm_covar, norm_correl, two_best_covar, \
            two_best_correl, frac_best_covar, frac_best_correl, weiss_stat, \
            mean_dist, var_dist, q_error, norm_scaled_covar]
    result.extend(correl_e)
    result.extend(covar_e)
    return result

def compare_distance_to_weiss(t):
    """Given a tree, calculate 2-sequence Weiss statistic and 3-seq variation.

    Returns tuple containing Weiss, eigenvector ratio, and var of 3-seq dist."""
    three_subs = all_p_from_tree(t)
    pair_subs = pair_p_from_tree(t)
    three_qs = get_qs(three_subs)
    pair_qs = get_qs(pair_subs)
    pairs_no_diag = map(array, map(take_off_diag, pair_qs))
    three_qs_copy = [i.copy() for i in three_qs]
    for i in three_qs: scale_trace(i)
    pairs_flat = flatten_q_matrices(pairs_no_diag)
    three_flat = flatten_q_matrices(three_qs)
    three_distances = dists_from_mean(three_flat)
    mean_three = average(three_distances)
    var_three = var(three_distances)
    pair_distances = dists_from_mean(pairs_flat)
    mean_pair = average(pair_distances)
    var_pair = var(pair_distances)
    
    three_qs_copy_flat = flatten_q_matrices(three_qs_copy)
    three_cor_matrix = corrcoef(three_qs_copy_flat)
    three_ratio = ratio_two_best(list(eigenvalues(three_cor_matrix).real))
    
    cov_matrix = cov(pairs_flat)
    covar_e = list(eigenvalues(cov_matrix).real)
    covar_e.sort()
    covar_e.reverse()
    two_ratio = ratio_two_best(covar_e)
    weiss_stat = weiss(covar_e)
    return [mean_pair, var_pair, mean_three, var_three, two_ratio, three_ratio, weiss_stat]
    

def distance_to_weiss_header():
    return ['mean_pair', 'var_pair', 'mean_three', 'var_three', 'two_ratio',\
            'three_ratio', 'weiss_stat']

def mean_var_eigen(t):
    """Returns mean and variance of branch lengths and 3-seq eigen of tree"""
    three_subs = all_p_from_tree(t)
    three_qs = get_qs(three_subs)
    three_qs_copy = [i.copy() for i in three_qs]
    for i in three_qs: scale_trace(i)
    three_flat = flatten_q_matrices(three_qs)
    three_distances = dists_from_mean(three_flat)
    mean_three = average(three_distances)
    var_three = var(three_distances)

    three_qs_copy_flat = flatten_q_matrices(three_qs_copy)
    three_cor_matrix = corrcoef(three_qs_copy_flat)
    three_ratio = ratio_two_best(list(eigenvalues(three_cor_matrix).real))

    return [mean_three, var_three, three_ratio]

def mean_var_eigen_header():
    return ['mean_three', 'var_three', 'ratio_three']

def svd_tests(n):
    """Runs some tests of svd on samples of n random matrices"""
    print 'SVD tests...'
    tests = {'random':random_mats, 'random q':random_q_mats, \
        'scaled q':scale_one_q, 'scaled many': scale_many_q}
    for name, f in tests.items():
            print name, ':\n', svd_q(f(n))

def tree_tests(n, leaves=8, slen=100, blen=0.1):
    """Runs some tests of svd on trees"""
    print "Tree tests..."
    print "single matrix"
    for i in range(n):
        t = balanced_one_q_tree(n=leaves, seq_length=slen, length=blen)
        ps = all_p_from_tree(t)
        print svd_q(get_qs(ps)).real

    print "two matrices"
    for i in range(n):
        t = balanced_two_q_tree(n=leaves, seq_length=slen, length=blen)
        ps = all_p_from_tree(t)
        print svd_q(get_qs(all_p_from_tree(t))).real
    print "many matrices"
    for i in range(n):
        ps = all_p_from_tree(t)
        t = balanced_multi_q_tree(n=leaves, seq_length=slen, length=blen)
        print svd_q(get_qs(all_p_from_tree(t))).real

def stats_from_tree(t):
    ps = all_p_from_tree(t)
    qs = get_qs(ps)
    return get_matrix_stats(qs)

def tree_stats(n, tree_f, result_f):
    for i in range(n):
        t = tree_f()
        yield result_f(t)

def tree_stats_header():
    header += ['norm_covar', 'norm_correl', 'two_best_covar', 'two_best_correl',
        'frac_best_covar', 'frac_best_correl', 'weiss_stat', 'mean_dist',
        'var_dist', 'q_error', 'norm_scaled_covar']
    for i in range(16):
        header.append('cor_'+str(i))
    for i in range(16):
        header.append('cov_'+str(i))
    return header

def tree_stats_analysis(n, header_f, result_f, leaves=16, slen=1000, blen=0.1, \
        tree_f=balanced_two_q_tree, tree_f_name='two_q'):
    make_tree = lambda: tree_f(n=leaves, seq_length=slen, length=blen)
    shared_params = [leaves, slen, blen, tree_f_name]
    shared_header = ['n', 'leaves', 'seq_len', 'branch_len', 'type']
    print '\t'.join(shared_header + header_f())
    for i, data in enumerate(tree_stats(n, make_tree, result_f)):
        result = [i] + shared_params + data
        print '\t'.join(map(str, result))

def perturb_lengths(t):
    """Perturbs the branch lengths in t by multiplying each by (0.5,1.5)."""
    for node in t.traverse():
        if node.Length:
            node.Length *= (0.5 + random())

def stats_from_tree_analysis(n):
    tree_stats_analysis(n, header_f=tree_stats_header, result_f=stats_from_tree)

def q_matrix_compare_analysis(n):
    tree_stats_analysis(n, header_f=distance_to_weiss_header, \
        result_f=compare_distance_to_weiss, tree_f=balanced_one_q_tree, \
        tree_f_name='one_q')
    tree_stats_analysis(n, header_f=distance_to_weiss_header, \
        result_f=compare_distance_to_weiss, tree_f=balanced_two_q_tree, \
        tree_f_name='two_q')
    tree_stats_analysis(n, header_f=distance_to_weiss_header, \
        result_f=compare_distance_to_weiss, tree_f=balanced_multi_q_tree, \
        tree_f_name='multi_q')

def change_length_analysis(n, branch_lengths=\
    [0.05, 0.1, 0.15, 0.2, 0.25,0.3, 0.35, 0.4], leaves=16, slen=1000):
    header = ['mean_three_m','mean_three_sd','var_three_m', 'var_three_sd', \
        'ratio_three_m', 'ratio_three_sd']
    conditions = {'one_q':balanced_one_q_tree, 'two_q':balanced_two_q_tree, \
        'multi_q':balanced_multi_q_tree}
    condition_names = ['one_q', 'two_q', 'multi_q']
    result_f = mean_var_eigen
    full_header = ['length'] + \
        [name+i for name in condition_names for i in header]
    print '\t'.join(full_header)
    for b in lengths:
        result = [b]
        for c in condition_names:
            tree_f = conditions[c]
            make_tree = lambda:tree_f(n=leaves,seq_length=slen,length=b,\
                    perturb=False)
            samples = list(tree_stats(n, make_tree, result_f))
            means = average(samples)
            stdevs = std(samples)
            for i in zip(means, stdevs):
                result.extend(i)
        print '\t'.join(map(str, result))

def tree_change_analysis(n):
    """Test whether we can see a change in stats for the first subtree w/ Q2"""
    header = ['mean_three','var_three', 'ratio_three']
    result_f = mean_var_eigen
    print '\t'.join(['n'] + header + header)
    for i in range(n):
        t, d = balanced_two_q_tree(n=16, seq_length=1000, length=0.1, \
            perturb=False, both=True)
        d_result = result_f(d)
        p_result = result_f(d.Parent)
        print '\t'.join(map(str, [i] + d_result + p_result))
        
if __name__ == '__main__':
    t = tree_controller()
