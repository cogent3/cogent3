#!/usr/bin/env python
"""Provides small utility functions for numpy arrays.
"""
from operator import mul, __getitem__ as getitem
from numpy import array, arange, logical_not, cumsum, where, compress, ravel,\
    zeros, put, take, sort, searchsorted, log, nonzero, sum,\
    sqrt, clip, maximum, reshape, argsort, argmin, repeat, product, identity,\
    concatenate, less, trace, newaxis, min, pi
from numpy.random import randint, normal
import numpy
from cogent.util.transform import cross_comb

def cartesian_product(lists):
    """Returns cartesian product of lists as list of tuples.

    WARNING: Explicitly constructs list in memory. Should use generator
    version in cogent.util.transform instead for most uses.

    Provided for compatibility.
    """
    return map(tuple, cross_comb(lists))

numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)
Int = numerictypes(int)
err = numpy.seterr(divide='raise')

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def gapped_to_ungapped(orig, gap_state, remove_mask=False):
    """Return array converting gapped to ungapped indices based on gap state.

    Will use == to test whether items equal the gapped state. Assumes character
    arrays.

    If remove_mask is True (default is False), will assign positions that are
    only in the gapped but not the ungapped version to -1 for easy detection.
    """
    return masked_to_unmasked(orig == gap_state, remove_mask)

def ungapped_to_gapped(orig, gap_state):
    """Returns array mapping indices in ungapped sequence to indices in orig.

    See documentation for unmasked_to_masked for more detail.
    """
    return unmasked_to_masked(orig == gap_state)

def masked_to_unmasked(mask, remove_mask=False):
    """Returns array mapping indices in orig to indices in ungapped.

    Specifically, for each position in orig, returns the index of the position
    in the unmasked sequence of the last non-masked character at or before
    that index (i.e. if the index corresponds to a masked position, will return
    the index of the previous non-masked position since the masked positions
    aren't in the unmasked sequence by definition).

    If remove_mask is True (the default is False), sets the masked positions
    to -1 for easy detection.
    """
    result = cumsum(logical_not(mask), axis=0) -1
    if remove_mask:
        result = where(mask, -1, result)
    return result

def unmasked_to_masked(mask):
    """Returns array mapping indices in ungapped to indices in original.

    Any position where the mask is True will be omitted from the final result.
    """
    return compress(logical_not(mask), arange(len(mask)))

def pairs_to_array(pairs, num_items=None, transform=None):
    """Returns array with same data as pairs (list of tuples).

    pairs can contain (first, second, weight) or (first, second) tuples.
    If 2 items in the tuple, weight will be assumed to be 1.

    num_items should contain the number of items that the pairs are chosen
    from. If None, will calculate based on the largest item in the actual
    list.

    transform contains a array that maps indices in the pairs coordinates
    to other indices, i.e. transform[old_index] = new_index. It is
    anticipated that transform will be the result of calling ungapped_to_gapped
    on the original, gapped sequence before the sequence is passed into
    something that strips out the gaps (e.g. for motif finding or RNA folding).

    WARNING: all tuples must be the same length! (i.e. if weight is supplied
    for any, it must be supplied for all.

    WARNING: if num_items is actually smaller than the biggest index in the
    list (+ 1, because the indices start with 0), you'll get an exception
    when trying to place the object. Don't do it.
    """
    #handle easy case
    if not pairs:
        return array([])
    data = array(pairs)
    #figure out if we're mapping the indices to gapped coordinates
    if transform is not None:
        #pairs of indices
        idx_pairs = take(transform, data[:,0:2].astype(Int), axis=0)
    else:
        idx_pairs = data[:,0:2].astype(Int)
    #figure out biggest item if not supplied
    if num_items is None:
        num_items = int(max(ravel(idx_pairs))) + 1
    #make result array
    result = zeros((num_items,num_items), Float)
    if len(data[0]) == 2:
        values = 1
    else:
        values = data[:,2]
    put(ravel(result), idx_pairs[:,0]*num_items+idx_pairs[:,1], values)
    return result

ln_2 = log(2)

def log2(x):
    """Returns the log (base 2) of x"
    
    WARNING: log2(0) will give -inf on one platform, but it might raise
    an error (Overflow or ZeroDivision on another platform. So don't rely
    on getting -inf in your downstream code.
    """
    return log(x)/ln_2

def safe_p_log_p(a):
    """Returns -(p*log2(p)) for every non-negative, nonzero p in a.

    a: numpy array

    WARNING: log2 is only defined on positive numbers, so make sure
    there are no negative numbers in the array.

    Always returns an array with floats in there to avoid unexpected
    results when applying it to an array with just integers.
    """
    c = array(a.copy(),Float)
    flat = ravel(c)
    nz_i = numpy.ravel(nonzero(maximum(flat,0)))
    nz_e = take(flat,nz_i, axis=0)
    log_nz = log2(nz_e)
    flat *= 0
    x = nz_e*-log_nz
    put(flat,nz_i,x)
    return c

def safe_log(a):
    """Returns the log (base 2) of each nonzero item in a.

    a: numpy array

    WARNING: log2 is only defined on positive numbers, so make sure
    there are no negative numbers in the array. Will either return an
    array containing floating point exceptions or will raise an 
    exception, depending on platform.

    Always returns an array with floats in there to avoid unexpected
    results when applying it to an array with just integers.
    """
    c = array(a.copy(),Float)
    flat = ravel(c)
    nz_i = numpy.ravel(nonzero(flat))
    nz_e = take(flat,nz_i, axis=0)
    log_nz = log2(nz_e)
    put(flat,nz_i,log_nz)
    return c

def row_uncertainty(a):
    """Returns uncertainty (Shannon's entropy) for each row in a IN BITS
    
    a: numpy array (has to be 2-dimensional!)

    The uncertainty is calculated in BITS not NATS!!!

    Will return 0 for every empty row, but an empty array for every empty column,
    thanks to this sum behavior:
    >>> sum(array([[]]),1)
    array([0])
    >>> sum(array([[]]))
    zeros((0,), 'l')
    """
    try:
        return sum(safe_p_log_p(a),1)
    except ValueError:
        raise ValueError, "Array has to be two-dimensional"

def column_uncertainty(a):
    """Returns uncertainty (Shannon's entropy) for each column in a in BITS

    a: numpy array (has to be 2-dimensional)

    The uncertainty is calculated in BITS not NATS!!!

    Will return 0 for every empty row, but an empty array for every empty column,
    thanks to this sum behavior:
    >>> sum(array([[]]),1)
    array([0])
    >>> sum(array([[]]))
    zeros((0,), 'l')

    """
    if len(a.shape) < 2:
        raise ValueError, "Array has to be two-dimensional"
    return sum(safe_p_log_p(a), axis=0)


def row_degeneracy(a,cutoff=.5):
    """Returns the number of characters that's needed to cover >= cutoff

    a: numpy array
    cutoff: number that should be covered in the array

    Example:
    [   [.1 .3  .4  .2],
        [.5 .3  0   .2],
        [.8 0   .1  .1]]
    if cutoff = .75: row_degeneracy -> [3,2,1]
    if cutoff = .95: row_degeneracy -> [4,3,3]

    WARNING: watch out with floating point numbers. 
    if the cutoff= 0.9 and in the array is also 0.9, it might not be found
    >>> searchsorted(cumsum(array([.6,.3,.1])),.9)
    2
    >>> searchsorted(cumsum(array([.5,.4,.1])),.9)
    1

    If the cutoff value is not found, the result is clipped to the
    number of columns in the array.
    """
    if not a.any():
        return []
    try:
        b = cumsum(sort(a)[:,::-1],1)
    except IndexError:
        raise ValueError, "Array has to be two dimensional"
    degen = [searchsorted(aln_pos,cutoff) for aln_pos in b]
    #degen contains now the indices at which the cutoff was hit
    #to change to the number of characters, add 1
    return clip(array(degen)+1,0,a.shape[1])


def column_degeneracy(a,cutoff=.5):
    """Returns the number of characters that's needed to cover >= cutoff

    a: numpy array
    cutoff: number that should be covered in the array

    Example:
    [   [.1 .8  .3],
        [.3 .2  .3],
        [.6 0   .4]]
    if cutoff = .75: column_degeneracy -> [2,1,3]
    if cutoff = .45: column_degeneracy -> [1,1,2]

    WARNING: watch out with floating point numbers. 
    if the cutoff= 0.9 and in the array is also 0.9, it might not be found
    >>> searchsorted(cumsum(array([.6,.3,.1])),.9)
    2
    >>> searchsorted(cumsum(array([.5,.4,.1])),.9)
    1

    If the cutoff value is not found, the result is clipped to the
    number of rows in the array. 
    """
    if not a.any():
        return []
    b = cumsum(sort(a,0)[::-1],axis=0)
    try:
        degen = [searchsorted(b[:,idx],cutoff) for idx in range(len(b[0]))]
    except TypeError:
        raise ValueError, "Array has to be two dimensional"
    #degen contains now the indices at which the cutoff was hit
    #to change to the number of characters, add 1
    return clip(array(degen)+1,0,a.shape[0])

def hamming_distance(x,y):
    """Returns the Hamming distance between two arrays.

    The Hamming distance is the number of characters which differ between
    two sequences (arrays).
    
    WARNING: This function truncates the longest array to the length of 
    the shortest one.
    
    Example:
    ABC, ABB -> 1
    ABCDEFG, ABCEFGH -> 4
    """
    shortest = min(map(len,[x,y]))
    return sum(x[:shortest] != y[:shortest], axis=0)

def norm(a):
    """Returns the norm of a matrix or vector

    Calculates the Euclidean norm of a vector.
    Applies the Frobenius norm function to a matrix 
    (a.k.a. Euclidian matrix norm)

    a = numpy array
    """
    return sqrt(sum((a*a).flat))

def euclidean_distance(a,b):
    """Returns the Euclidean distance between two vectors/arrays
    a,b: numpy vectors or arrays

    WARNING: this method is NOT intended for use on arrays of different
    sizes, but no check for this has been built in. 
    """
    return norm(a-b)
def count_simple(a, alphabet_len):
    """Counts items in a. """
    result = zeros(alphabet_len, Int)
    for i in ravel(a):
        result[i] += 1
    return result

def count_alphabet(a, alphabet_len):
    """Counts items in a, using =="""
    #ensure behavior is polymorphic with count_simple
    if not alphabet_len:
        raise IndexError, "alphabet_len must be > 0"
    result = zeros(alphabet_len, Int)
    a = ravel(a)
    for i in range(alphabet_len):
        result[i] = sum(a == i)
    return result

def is_complex(m):
    """Returns True if m has a complex component."""
    return m.dtype.char == 'D'

def is_significantly_complex(m, threshold=0.1):
    """Returns True if the sum of m's imaginary component exceeds threshold."""
    if is_complex(m):
        if sum(sum(abs(m.imag))) > threshold:
            return True
    return False

def has_neg_off_diags(m):
    """Returns True if m has negative off-diagonal elements."""
    return min(ravel(m * logical_not(identity(len(m))))) < 0

def has_neg_off_diags_naive(m):
    """Returns True if m has off-diagonal elements.
 
    Naive, slow implementation -- don't use. Primarily here to check
    correctness of faster implementation.
    """
    working = m.copy()
    for i in range(len(working)):
        working[i][i] = 0
    if min(ravel(working)) < 0:
        return True
    else:
        return False

def sum_neg_off_diags(m):
    """Returns sum of negative off-diags in m."""
    return sum(compress(ravel(less(m,0)), \
        ravel((m * logical_not(identity(len(m)))))))

def sum_neg_off_diags_naive(m):
    """Returns sum of negative off-diags in m.

    Naive, slow implementation -- don't use. Primarily here to check correctness
    of faster implementation.
    """
    sum = 0
    for row_i, row in enumerate(m):
        for col_i, item in enumerate(row):
            if (row_i != col_i) and (item < 0):
                sum += item
    return sum

def scale_row_sum(m, val=1):
    """Scales matrix in place so that each row sums to val (default: 1).

    WARNING: will use 'integer division', not true division, if matrix is
    an integer data type.
    """
    m /= (sum(m, axis=1)/val)[:,newaxis]

def scale_row_sum_naive(m, val=1):
    """Scales matrix in place so that each row sums to val (default:1).

    Naive implementation -- don't use. Primarily here to check correctness.
    
    WARNING: will use 'integer division'.
    """
    for row in m:
        row_sum = sum(row)
        row /= (row_sum / val)

def scale_trace(m, val=-1):
    """Scales matrix in place so that trace of result is val (default: -1).

    WARNING: if trace is opposite sign to val, inverts sign of all elements
    in the matrix.

    WARNING: will use 'integer division', not true division, if matrix is
    an integer data type.
    """
    m *= val/trace(m)

def abs_diff(first, second):
    """Calculates element-wise sum of abs(first - second).

    Return value may be real or complex.
    """
    return sum(ravel(abs(first-second)))

def sq_diff(first, second):
    """Calculates element-wise sum of (first - second)**2.

    Return value may be real or complex.
    """
    diff = first - second
    return sum(ravel((diff*diff)))

def norm_diff(first, second):
    """Returns square root of sq_diff, normalized to # elements."""
    size = len(ravel(first))
    return sqrt(sq_diff(first, second))/size

def without_diag(a):
    """Returns copy of square matrix a, omitting diagonal elements."""
    return array([concatenate((r[:i], r[i+1:])) for i, r in enumerate(a)])

def with_diag(a, d):
    """Returns copy of matrix a with d inserted as diagonal to yield square."""
    rows, cols = a.shape
    result = zeros((rows, cols+1), a.dtype.char)
    for i, r in enumerate(a):
        result_row = result[i]
        result_row[:i] = r[:i]
        result_row[i] = d[i]
        result_row[i+1:] = r[i:]
    return result

def only_nonzero(a):
    """Returns elements of a where the first element of a[i] is nonzero.
    
    Result is a new array and does not share data with the original.

    NOTE: This is designed for arrays of rate matrices. If the first element
    of the rate matrix is zero, then the row must be all zero (since the row
    sums to zero with the first element being equal in magnitude but opposite
    in sign to the sum of the other elements). If the row is all zero, then
    the rate matrix is almost certainly invalid and should be excluded from
    further analyses.
    """
    first_element_selector = [0] * len(a.shape)
    first_element_selector[0] = slice(None,None)
    return take(a, numpy.ravel(nonzero(a[first_element_selector])), axis=0)

def combine_dimensions(m, dim):
    """Aggregates all dimensions of m between dim and the end.

    In other words, combine_dimensions(m, 3) flattens the first 3 dimensions.
    Similarly, combine_dimensions(m, -2) flattens the last two dimensions.

    WARNING: Result shares data with m.
    """
    #if we're not combining more than one dimension, return the array unchanged
    if abs(dim) <= 1:
        return m
    #otherwise, figure out the new shape and reshape the array
    shape = m.shape
    if dim < 0:     #counting from end
        return reshape(m, shape[:dim] + (product(shape[dim:]),))
    else:           #counting from start
        return reshape(m, (product(shape[:dim]),) + shape[dim:])

def split_dimension(m, dim, shape=None):
    """Splits specified dimension of m into shape.

    WARNING: product of terms in shape must match size of dim. For example,
    if the length of dim is 12, shape could be (4,3) or (2,6) but not (2,3).

    Result shares data with m.
    """
    curr_dims = m.shape
    num_dims = len(curr_dims)
    #if shape not specified, assume it was square
    if shape is None:
        shape = (sqrt(curr_dims[dim]),)*2
    #raise IndexError if index out of bounds
    curr_dims[dim]
    #fix negative indices
    if dim < 0:
        dim = num_dims + dim
    #extract the relevant region and reshape it
    return reshape(m, curr_dims[:dim] + shape + curr_dims[dim+1:])

def non_diag(m):
    """From a sequence of n flattened 2D matrices, returns non-diag elements.

    For example, for an array containing 20 9-element row vectors, returns
    an array containing 20 6-element row vectors that omit elements 0, 4, and 8.
    """
    num_rows, num_elements = m.shape
    side_length = int(sqrt(num_elements))
    wanted = numpy.ravel(nonzero(logical_not(identity(side_length).flat)))
    all_wanted = repeat([wanted], num_rows,axis=0)
    all_wanted += (arange(num_rows) * num_elements)[:,newaxis]
    return reshape(take(ravel(m), ravel(all_wanted), axis=0), \
        (num_rows, num_elements-side_length))

def perturb_one_off_diag(m, mean=0, sd=0.01, element_to_change=None):
    """Perturbs one off-diagonal element of rate matrix m by random number.

    mean: mean of distribution to sample from. Default 0.
    sd: standard deviation of distribution to sample from. Default 0.05.

    Error model is additive.
    
    WARNING: may reverse sign of element in some cases!
    WARNING: if an element is specified, the coordinate is relative to the
    flat array _without_ the diagonal, _not_ relative to the original array!
    e.g. for a 4x4 array, element 8 refers to a[2][3], _not_ a[2][0].
    """
    #get the elements that are allowed to change
    elements = without_diag(m)
    flat = ravel(elements)
    #pick an element to change if it wasn't specified
    if element_to_change is None:
        element_to_change = randint(0, len(flat))
    #change the element, pack the elements back into the array, and return
    #the result.
    flat[element_to_change] += normal(mean, sd)
    return with_diag(elements, -sum(elements, 1))

def perturb_one_off_diag_fixed(m, size):
    """Perturbs a random off-diag element of rate matrix m by factor of size."""
    elements = without_diag(m)
    flat = ravel(elements)
    element_to_change = randint(0, len(flat))
    flat[element_to_change] *= (1.0 + size)
    return with_diag(elements, -sum(elements, 1))

def perturb_off_diag(m, mean=0, sd=0.01):
    """Perturbs all off-diagonal elements of m by adding random number.

    mean: mean of distribution to sample from. Default 0.
    sd: standard deviation of distribution to sample from. Default 0.01.

    WARNING: may reverse sign of element!
    """
    elements = without_diag(m)
    random = normal(mean, sd, elements.shape)
    result = elements + random
    return with_diag(result, -sum(result, 1))

def perturb_off_diag_frac(m, size):
    """Perturbs all off-diagonal elements of m by about a specified fraction.

    mean: mean of size (relative to each element) of change to make.
    Will never reverse sign of an element.
    """
    elements = without_diag(m)
    random = normal(0, size_to_stdev(size), elements.shape)
    result = elements * abs((1.0+random))   #ensure always positive
    return with_diag(result, -sum(result, 1))


def size_to_stdev(size):
    """From desired mean deviation, returns sd for N(0,sd).

    Uses method of Altman 1993, as described in:
    http://www-users.york.ac.uk/~mb55/talks/halfnor.pdf

    ...where E(X) = sqrt(2*sigma/pi)
    """
    return size*size*pi/2.0

def merge_samples(*samples):
    """Merges list of samples into array of [vals,dists].

    value of each sample corresponds to its position in the list.
    
    e.g. for [1,2,3] and [4,5,6], result will be:
    array([[1,2,3,4,5,6],[0,0,0,1,1,1]])
    """
    return concatenate([array((a, zeros(a.shape)+i)) \
        for i,a in enumerate(samples)], 1)

def sort_merged_samples_by_value(a):
    """Sorts result of merge_samples by value (i.e. first row)."""
    return take(a, argsort(a[0]), 1)

def classifiers(*samples):
    """Returns all 1D classifier separating first sample from the remainder.

    Returns [(cut_value, fp, fn, tp, tn) for i in cuts].
    """
    if len(samples) <= 1:
        raise TypeError, "optimal_classifier needs at least 2 distributions."
    vals, labels = sort_merged_samples_by_value(merge_samples(*samples))
    n = len(vals)
    num_positives = len(samples[0])
    num_negatives = n - num_positives
    #only want to check indices where the next value is different
    to_check = numpy.ravel(nonzero(vals[:-1] - vals[1:]))
    result = []
    for index in to_check:
        i = index+1 #because it changes at the value _after_ the specified index
        fp = sum(labels[:i] != 0)
        fn = sum(labels[i:] == 0)
        tp = num_positives - fn
        tn = num_negatives - fp
        reversed = tp + tn < fp + fn
        if reversed:
            tp, tn, fp, fn = fp, fn, tp, tn
        result.append((i, reversed, fp, fn, tp, tn))
    return result

def minimize_error_count(classifiers):
    """Returns the classifier from a list of classifiers that minimizes errors.

    Errors are defined as #fp + # fn.

    If multiple classifiers have equal scores, returns an arbitrary one.
    """
    c = array(classifiers)
    return classifiers[argmin(sum(c[:,2:4],1))]

def minimize_error_rate(classifiers):
    """Returns the classifier from a list of classifiers that minimizes errors.

    Errors are defined as (#fp/#p) + (# fn/#n).

    If multiple classifiers have equal scores, returns an arbitrary one.
    """
    c = array(classifiers)
    return classifiers[argmin(\
        1.0*c[:,2]/(c[:,2]+c[:,5])+1.0*c[:,3]/(c[:,3]+c[:,4]))]

def mutate_array(a, sd, mean=0):
    """Return mutated copy of the array (or vector), adding mean +/- sd."""
    return a + normal(mean, sd, a.shape)

