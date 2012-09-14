#!/usr/bin/env python
"""Functions for doing fast multidimensional scaling of distances
using Nystrom/LMDS approximation ans Split and Combine MDS.

===================== Nystrom =====================

Approximates an MDS mapping / a (partial) PCoA solution of an
(unknown) full distance matrix using a k x n seed distance matrix.

Use if you have a very high number of objects or if each distance
caculation is expensive. Speedup comes from two factors: 1. Not all
distances are calculated but only k x n. 2. Eigendecomposition is only
applied to a k x k matrix.

Calculations done after Platt (2005). See
http://research.microsoft.com/apps/pubs/?id=69185 :
``This paper unifies the mathematical foundation of three
multidimensional scaling algorithms: FastMap, MetricMap, and Landmark
MDS (LMDS). All three algorithms are based on the Nystrom
approximation of the eigenvectors and eigenvalues of a matrix. LMDS is
applies the basic Nystrom approximation, while FastMap and MetricMap
use generalizations of Nystrom, including deflation and using more
points to establish an embedding. Empirical experiments on the Reuters
and Corel Image Features data sets show that the basic Nystrom
approximation outperforms these generalizations: LMDS is more accurate
than FastMap and MetricMap with roughly the same computation and can
become even more accurate if allowed to be slower.``


Assume a full distance matrix D for N objects:

        /  E  |  F \
    D = |-----|----|
        \ F.t |  G /
        

The correspondong association matrix or centered inner-product
matrix K is:

        /  A  |  B \
    K = |-----|----|
        \ B.t |  C /
        
where A and B are computed as follows

    A_ij = - 0.50 * (E_ij^2 -
                 1/m SUM_p E_pj^2 -
                 1/m SUM_q E_iq^2 +
                 1/m^2 SUM_q E_pq^2

B is computed as in Landmark MDS, because it's simpler and works
better according to Platt):

    B_ij = - 0.50 * (F_ij^2 - 1/m SUM_q E_iq^2)


In order to approximate an MDS mapping for full matrix you only need E
and F from D as seed matrix. This will mimick the distances for m seed
objects. E is of dimension m x m and F of m x (N-m)
    
E and F are then used to approximate and MDS solution x for the full
distance matrix:


x_ij =  sqrt(g_j) * U_ij, if i<=m
        and
        SUM_p B_pi U_pj / sqrt(g_j)

where U_ij is the i'th component of the jth eigenvector of A (see
below) and g_j is the j'th eigenvalue of A.  The index j only runs
from 1 to k in order to make a k dimensional embedding.

===================== SCMDS =====================

The is a Python/Numpy implementation of SCMDS:
Tzeng J, Lu HH, Li WH.
Multidimensional scaling for large genomic data sets.
BMC Bioinformatics. 2008 Apr 4;9:179.
PMID: 18394154

The basic idea is to avoid the computation and eigendecomposition of
the full pairwise distance matrix. Instead only compute overlapping
submatrices/tiles and their corresponding MDS separately. The
solutions are then joined using an affine mapping approach.

=================================================
"""


from numpy import sign, floor, sqrt, power, mean, array
from numpy import matrix, ones, dot, argsort, diag, eye
from numpy import zeros, concatenate, ndarray, kron, argwhere
from numpy.linalg import eig, eigh, qr
from random import sample
import time


__author__ = "Adreas Wilm"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Fabian Sievers <fabian.sievers@ucd.ie>",
               "Daniel McDonald <wasade@gmail.com>", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Andreas Wilm"
__email__ = "andreas.wilm@ucd.ie"
__status__ = "Development"


# print simple timings
PRINT_TIMINGS = False


def rowmeans(mat):
    """Returns a `column-vector` of row-means of the 2d input array/matrix
    """

    if not len(mat.shape)==2:
        raise ValueError, "argument is not a 2D ndarray"

    #nrows = mat.shape[0]
    ## create a column vector (hack!)
    #cv = matrix(arange(float(nrows)))
    #cv = cv.T
    #for i in range(nrows):
    #    cv[i] = mat[i].mean()
    #
    # As pointed out by Daniel the above is much simpler done in Numpy:
    cv = matrix(mean(mat, axis=1).reshape((mat.shape[0], 1)))

    return cv



def affine_mapping(matrix_x, matrix_y):
    """Returns an affine mapping function.

    Arguments are two different MDS solutions/mappings (identical
    dimensions) on the the same objects/distances.
    
    Affine mapping function:
    Y = UX + kron(b,ones(1,n)), UU' = I
    X = [x_1, x_2, ... , x_n]; x_j \in R^m
    Y = [y_1, y_2, ... , y_n]; y_j \in R^m

    From Tzeng 2008:
    `The projection of xi,2 from q2 dimension to q1 dimension
    induces computational errors (too). To avoid this error, the
    sample number of the overlapping region is important. This
    sample number must be large enough so that the derived
    dimension of data is greater or equal to the real data`

    Notes:
    - This was called moving in scmdscale.m
    - We work on tranposes, i.e. coordinates for objects are in columns
    

    Arguments:
    - `matrix_x`: first mds solution (`reference`)
    - `matrix_y`: seconds mds solution
    
    Return:
    A tuple of unitary operator u and shifting operator b such that:
    y = ux + kron(b, ones(1,n))
    """


    if matrix_x.shape != matrix_y.shape:
        raise ValueError, \
            "input matrices are not of same size"

    if not matrix_x.shape[0] <= matrix_x.shape[1]:
        raise ValueError, \
            "input matrices should have more columns than rows"
    
    # Have to check if we have not more rows than columns, otherwise,
    # the qr function below might behave differently in the matlab
    # prototype, but not here. In addition this would mean that we
    # have more dimensions than overlapping objects which shouldn't
    # happen
    #
    # matlab code uses economic qr mode (see below) which we can't use
    # here because we need both return matrices.
    #
    # see
    # http://www.mathworks.com/access/helpdesk/help/techdoc/index.html?/access/helpdesk/help/techdoc/ref/qr.html
    # [Q,R] = qr(A,0) produces the economy-size decomposition.
    # If m > n, only the first n columns of Q and
    # the first n rows of R are computed.
    # If m<=n, this is the same as [Q,R] = qr(A)
    #
    # [Q,R] = qr(A), where A is m-by-n, produces
    # an m-by-n upper triangular matrix R and
    # an m-by-m unitary matrix Q so that A = Q*R
    #
    # That's why we check above with an assert

    ox = rowmeans(matrix_x)
    oy = rowmeans(matrix_y)

    mx = matrix_x - kron(ox, ones((1, matrix_x.shape[1])))
    my = matrix_y - kron(oy, ones((1, matrix_x.shape[1])))

    (qx, rx) = qr(mx)
    (qy, ry) = qr(my)
    
    # sign correction
    #
    # Daniel suggest to use something like
    # [arg]where(sign(a.diagonal()) != sign(b.diagonal())) and then
    # iterate over the results. Couldn't figure out how to do this
    # properly :(
    
    #idx = argwhere(sign(rx.diagonal()) != sign(ry.diagonal()))
    #for i in idx:
    #	qy[:,i] *= -1
    for i in range(qx.shape[1]):
        if sign(rx[i, i]) != sign(ry[i, i]):
            qy[:, i] *= -1



    # matrix multiply: use '*' as all arguments are of type matrix
    ret_u = qy * qx.transpose()
    ret_b = oy - ret_u * ox
    
    return (ret_u, ret_b)



def adjust_mds_to_ref(mds_ref, mds_add, n_overlap):
    """Transforms mds_add such that the overlap mds_ref and mds_add
    has same configuration.

    As overlap (n_overlap objects) we'll use the end of mds_ref
    and the beginning of mds_add
    
    Both matrices must be of same dimension (column numbers) but
    can have different number of objects (rows) because only
    overlap will be used.

    Arguments:
    - `mds_ref`: reference mds solution
    - `mds_add`: mds solution to adjust
    - `n_overlap`: overlap size between mds_ref and mds_add
    
    Return:
    Adjusted version of mds_add which matches configuration of mds_ref
    """

    if mds_ref.shape[1] != mds_add.shape[1]:
        raise ValueError, \
            "given mds solutions have different dimensions"
    if not (mds_ref.shape[0] >= n_overlap and mds_add.shape[0] >= n_overlap):
        raise ValueError, \
            "not enough overlap between given mds mappings"
        
    # Use transposes for affine_mapping!
    overlap_ref = mds_ref.transpose()[:, -n_overlap:]
    overlap_add = mds_add.transpose()[:, 0:n_overlap]
    (unitary_op, shift_op) = affine_mapping(overlap_add, overlap_ref)
    # paranoia: unitary_op is of type matrix, make sure mds_add
    # is as well so that we can use '*' for matrix multiplication
    mds_add_adj = unitary_op * matrix(mds_add.transpose()) + \
                         kron(shift_op, ones((1, mds_add.shape[0])))
    mds_add_adj = mds_add_adj.transpose()

    return mds_add_adj



def recenter(joined_mds):
    """Recenter an Mds mapping that has been created by joining, i.e.
    move points so that center of gravity is zero.

    Note:
    Not sure if recenter is a proper name, because I'm not exactly
    sure what is happening here

    Matlab prototype from Tzeng et al. 2008:
	 X = zero_sum(X); # subtract column means
	 M = X'*X;
	 [basis,L] = eig(M);
	 Y = X*basis;
	 return Y = Y(:,end:-1:1);

    Arguments:
    - `mds_combined`:

    Return:
    Recentered version of `mds_combined`
    """

    # or should we cast explictely?
    if not isinstance(joined_mds, matrix):
        raise ValueError, "mds solution has to be of type matrix"

    # As pointed out by Daniel: the following two loop can be done in
    # one if you pass down the axis variable to means()
    #
    #colmean = []
    #for i in range(joined_mds.shape[1]):
    #    colmean.append(joined_mds[:, i].mean())
    #for i in range(joined_mds.shape[0]):
    #    joined_mds[i, :] = joined_mds[i, :] - colmean
    #
    joined_mds = joined_mds - joined_mds.mean(axis=0)

    matrix_m = dot(joined_mds.transpose(), joined_mds)
    (eigvals, eigvecs) = eig(matrix_m)    

    # Note / Question: do we need sorting?
    # idxs_ascending = eigvals.argsort()
    # idxs_descending = eigvals.argsort()[::-1]# reverse!
    # eigvecs = eigvecs[idxs_ascending]
    # eigvals = eigvals[idxs_ascending]
    
    # joined_mds and eigvecs are of type matrix so use '*' for
    # matrix multiplication
    joined_mds = joined_mds * eigvecs
    
    # NOTE: the matlab protoype reverses the vector before
    # returning. We don't because I don't know why and results are
    # good
    
    return joined_mds



def combine_mds(mds_ref, mds_add, n_overlap):
    """Returns a combination of the two MDS mappings mds_ref and
    mds_add.

    This is done by finding an affine mapping on the
    overlap between mds_ref and mds_add and changing mds_add
    accordingly.

    As overlap we use the last n_overlap objects/rows in mds_ref and
    the first n_overlap objects/rows in mds_add.

    The overlapping part will be replaced, i.e. the returned
    combination has the following row numbers:
    mds_ref.nrows + mds_add.nrows - overlap
    
    The combined version will eventually need recentering.
    See recenter()
    
    Arguments:
    - `mds_ref`: reference mds mapping
    - `mds_add`: mds mapping to add
    """

    if mds_ref.shape[1] != mds_add.shape[1]:
        raise ValueError, \
            "given mds solutions have different dimensions"
    if not mds_ref.shape[0] >= n_overlap:
        raise ValueError, \
           "not enough items for overlap in mds_ref"
    if not mds_add.shape[0] >= n_overlap:
        raise ValueError, \
           "not enough items for overlap in mds_add"
           
    mds_add_adj = adjust_mds_to_ref(mds_ref, mds_add, n_overlap)

    combined_mds = concatenate(( \
        mds_ref[0:mds_ref.shape[0]-n_overlap, :], mds_add_adj))

    return combined_mds



def cmds_tzeng(distmat, dim = None):
    """Calculate classical multidimensional scaling on a distance matrix.

    Faster than default implementation of dim is smaller than
    distmat.shape

    Arguments:
    - `distmat`: distance matrix (non-complex, symmetric ndarray)
    - `dim`:     wanted dimensionality of MDS mapping (defaults to distmat dim)

    Implementation as in Matlab prototype of SCMDS, see 
    Tzeng J et al. (2008), PMID: 18394154
    """
          
    if not isinstance(distmat, ndarray):
        raise ValueError, \
            "Input matrix is not a ndarray"
    (m, n) = distmat.shape
    if m != n:
        raise ValueError, \
            "Input matrix is not a square matrix"
    if not dim:
        dim = n

    # power goes wrong here if distmat is ndarray because of matrix
    # multiplication syntax difference between array and
    # matrix. (doesn't affect gower's variant). be on the safe side
    # and convert explicitely (it's local only):
    distmat = matrix(distmat)

    h = eye(n) - ones((n, n))/n
    assocmat = -h * (power(distmat, 2)) * h / 2
    #print "DEBUG assocmat[:3] = %s" % assocmat[:3]
    
    (eigvals, eigvecs) = eigh(assocmat)
    
    # Recommended treatment of negative eigenvalues (by Fabian): use
    # absolute value (reason: SVD does the same)
    eigvals = abs(eigvals)

    ind = argsort(eigvals)[::-1]
    eigvals = eigvals[ind]
    eigvecs = eigvecs[:, ind] 
    eigvals = eigvals[:dim]
    eigvecs = eigvecs[:, :dim]
    
    eigval_diagmat = matrix(diag(sqrt(eigvals)))
    eigvecs = eigval_diagmat * eigvecs.transpose()
    return (eigvecs.transpose(), eigvals)




class CombineMds(object):
    """
    Convinience class for joining MDS mappings. Several mappings can
    be added.

    The is uses the above Python/Numpy implementation of SCMDS.
    See Tzeng  et al. 2008, PMID: 18394154
    """
    

    def __init__(self, mds_ref=None):
        """
        Init with reference MDS
        """
        
        self._mds = mds_ref
        self._need_centering = False


    def add(self, mds_add, overlap_size):
        """Add a new MDS mapping to existing one
        """
            
        if self._mds == None:
            self._mds = mds_add
            return
        
        if not self._mds.shape[0] >= overlap_size:
            raise ValueError, \
                "not enough items for overlap in reference mds"
        if not mds_add.shape[0] >= overlap_size:
            raise ValueError, \
                "not enough items for overlap in mds to add"
      
        self._need_centering = True
        self._mds = combine_mds(self._mds, mds_add, overlap_size)


    def getFinalMDS(self):
        """Get final, combined MDS solution
        """
        
        if self._need_centering:
            self._mds = recenter(self._mds)

        return self._mds
    












def calc_matrix_b(matrix_e, matrix_f):
    """Calculates k x n-k matrix b of association matrix K

    See eq (14) and (15) in Platt (2005)

    This is where Nystrom and LMDS differ: LMDS version leaves the
    constant centering term out, This simplifies the computation.
    
    Arguments:
    - `matrix_e`:
       k x k part of the kxn seed distance matrix
    - `matrix_f`: 
       k x n-k part of the kxn seed distance matrix
    """

    (nrows, ncols) = matrix_f.shape

    if matrix_e.shape[0] != matrix_e.shape[1]:
        raise ValueError, "matrix_e should be quadratic"
    if matrix_f.shape[0] != matrix_e.shape[0]:
        raise ValueError, \
            "matrix_e and matrix_f should have same number of rows"

    nseeds = matrix_e.shape[0]

    # row_center_e was also precomputed in calc_matrix_a but
    # computation is cheap
    #
    row_center_e = zeros(nseeds)
    for i in range(nseeds):
        row_center_e[i] = power(matrix_e[i, :], 2).sum()/nseeds

    # The following is not needed in LMDS but part of original Nystrom
    #
    # ncols_f = matrix_f.shape[1]
    # col_center_f = zeros(ncols_f)
    # for i in range(ncols_f):
    #    col_center_f[i] = power(matrix_f[:, i], 2).sum()/nseeds
    #
    # subtract col_center_f[j] below from result[i, j] and you have
    # nystrom. dont subtract it and you have lmds
    

    #result = zeros((nrows, ncols))
    #for i in xrange(nrows):
    #    for j in xrange(ncols):
    # 
    #        # - original one line version:
    #        #
    #        #result[i, j] = -0.50 * (
    #        #    power(matrix_f[i, j], 2) -
    #        #    row_center_e[i])
    #        #
    #        # - optimised version avoiding pow(x,2)
    #        #   3xfaster on a 750x3000 seed-matrix
    #        fij = matrix_f[i, j]
    #        result[i, j] = -0.50 * (
    #            fij * fij -
    #            row_center_e[i])
    #
    # - optimised single line version of the code block above. pointed out
    #   by daniel. 20xfaster on a 750x3000 seed-matrix. cloning idea
    #   copied from
    #   http://stackoverflow.com/questions/1550130/cloning-row-or-column-vectors
    result = -0.5 * (matrix_f**2 - array([row_center_e, ]*ncols).transpose())

    return result



def calc_matrix_a(matrix_e):
    """Calculates k x k matrix a of association matrix K
    
    see eq (13) in Platt from symmetrical matrix E
    
    A_ij = - 0.50 * (E_ij^2 -
                   1/m SUM_p E_pj^2 -
                   1/m SUM_q E_iq^2 +
                   1/m^2 SUM_q E_pq^2

    Row and colum centering terms (E_pj and E_iq) are identical
    because we use a k x k submatrix of a symmetrical distance

    m equals here ncols or ncols of matrix_e
    we call it nseeds
    
    Arguments:
    - `matrix_e`:
       k x k part of the kxn seed distance matrix
    """

    if matrix_e.shape[0] != matrix_e.shape[1]:
        raise ValueError, "matrix_e should be quadratic"
    nseeds = matrix_e.shape[0]

    row_center = zeros(nseeds)
    for i in range(nseeds):
        row_center[i] = power(matrix_e[i, :], 2).sum()/nseeds

    # E should be symmetric, i.e. column and row means are identical.
    # Why is that not mentioned in the papers? To be on the safe side
    # just do this:
    # col_center = zeros(nseeds)
    # for i in range(nseeds):
    #     col_center[i] = power(matrix_e[:, i], 2).sum()/nseeds
    # or simply:
    col_center = row_center

    grand_center = power(matrix_e, 2).sum()/power(nseeds, 2)

    # E is symmetric and so is A, which is why we don't need to loop
    # over the whole thing
    #
    # FIXME: Optimize
    #
    result = zeros((nseeds, nseeds))
    for i in range(nseeds):
        for j in range(i, nseeds):
            # avoid pow(x,2). it's slow
            eij_sq = matrix_e[i, j] * matrix_e[i, j]
            result[i, j] = -0.50 * (
                eij_sq -
                col_center[j] -
                row_center[i] +
                grand_center)
            if i != j:
                result[j, i] = result[i, j]
            
    return result



def build_seed_matrix(fullmat_dim, seedmat_dim, getdist, permute_order=True):
    """Builds a seed matrix of shape seedmat_dim x fullmat_dim

    Returns seed-matrix and indices to restore original order (needed
    if permute_order was True)

    Arguments:
    - `fullmat_dim`:
      dimension of the unknown (square, symmetric) "input" matrix 
    - `seedmat_dim`:
      requested dimension of seed matrix.
    - `getdist`:
      distance function to compute distances. should take two
      arguments i,j with an index range of 0..fullmat_dim-1
    - `permute_order`:
       if permute_order is false, seeds will be picked sequentially.
       otherwise randomly
    """

    if not seedmat_dim < fullmat_dim:
        raise ValueError, \
            "dimension of seed matrix must be smaller than that of full matrix"
    if not callable(getdist):
        raise ValueError, "distance getter function not callable"

    if permute_order:
        picked_seeds = sample(range(fullmat_dim), seedmat_dim)
    else:
        picked_seeds = range(seedmat_dim)
    #assert len(picked_seeds) == seedmat_dim, (
    #    "mismatch between number of picked seeds and seedmat dim.")

    # Putting picked seeds/indices at the front is not enough, 
    # need to change/correct all indices to maintain consistency
    #
    used_index_order = range(fullmat_dim)
    picked_seeds.sort() # otherwise the below fails
    for i, seed_idx in enumerate(picked_seeds):
        used_index_order.pop(seed_idx-i)
    used_index_order = picked_seeds + used_index_order

    # Order is now determined in used_index_order
    # first seedmat_dim objects are seeds
    # now create seedmat
    #
    t0 = time.clock()
    seedmat = zeros((len(picked_seeds), fullmat_dim))
    for i in range(len(picked_seeds)):
        for j in range(fullmat_dim):
            if i < j:
                seedmat[i, j] = getdist(used_index_order[i],
                                        used_index_order[j])
            elif i == j:
                continue
            else:
                seedmat[i, j] = seedmat[j, i]


    restore_idxs = argsort(used_index_order)
    if PRINT_TIMINGS:
        print("TIMING(%s): Seedmat calculation took %f CPU secs" % 
              (__name__, time.clock() - t0))

    # Return the seedmatrix and the list of indices which can be used to
    # recreate original order
    return (seedmat, restore_idxs)

        
    
def nystrom(seed_distmat, dim):
    """Computes an approximate MDS mapping of an (unknown) full distance
    matrix using a kxn seed distance matrix.

    Returned matrix has the shape seed_distmat.shape[1] x dim
    """


    if not seed_distmat.shape[0] < seed_distmat.shape[1]:
        raise ValueError, \
            "seed distance matrix should have less rows than column"
    if not dim <= seed_distmat.shape[0]:
        raise ValueError, \
            "number of rows of seed matrix must be >= requested dim"
    
    nseeds = seed_distmat.shape[0]
    nfull = seed_distmat.shape[1]

    # matrix E: extract columns 1--nseed
    #
    matrix_e = seed_distmat[:, 0:nseeds]
    #print("INFO: Extracted Matrix E which is of shape %dx%d" %
    #      (matrix_e.shape))

    # matrix F: extract columns nseed+1--end
    #
    matrix_f = seed_distmat[:, nseeds:]
    #print("INFO: Extracted Matrix F which is of shape %dx%d" % 
    #      (matrix_f.shape))

    # matrix A
    #
    #print("INFO: Computing Matrix A")
    t0 = time.clock()
    matrix_a = calc_matrix_a(matrix_e)
    if PRINT_TIMINGS:
        print("TIMING(%s): Computation of A took %f CPU secs" %
              (__name__, time.clock() - t0))
    
    # matrix B
    #
    #print("INFO: Computing Matrix B")
    t0 = time.clock()
    matrix_b = calc_matrix_b(matrix_e, matrix_f)
    if PRINT_TIMINGS:
        print("TIMING(%s): Computation of B took %f CPU secs" %
              (__name__, time.clock() - t0))


    #print("INFO: Eigendecomposing A")
    t0 = time.clock()
    # eigh: eigen decomposition for symmetric matrices
    # returns: w, v
    # w : ndarray, shape (M,)
    #     The eigenvalues, not necessarily ordered.
    # v : ndarray, or matrix object if a is, shape (M, M)
    #     The column v[:, i] is the normalized eigenvector corresponding
    #     to the eigenvalue w[i].
    # alternative is svd: [U, S, V] = numpy.linalg.svd(matrix_a)
    (eigval_a, eigvec_a) = eigh(matrix_a)
    if PRINT_TIMINGS:
        print("TIMING(%s): Eigendecomposition of A took %f CPU secs" % 
              (__name__, time.clock() - t0))    
    # Sort descending
    ind = argsort(eigval_a)
    ind = ind[::-1]
    eigval_a = eigval_a[ind]
    eigvec_a = eigvec_a[:, ind]

   
    #print("INFO: Estimating MDS coords")
    t0 = time.clock()
    result = zeros((nfull, dim)) # X in Platt 2005
    # Preventing negative eigenvalues by using abs value. Other option
    # is to set negative values to zero. Fabian recommends using
    # absolute values (as in SVD)
    sqrt_eigval_a = sqrt(abs(eigval_a))
    for i in xrange(nfull):
        for j in xrange(dim):
            if i+1 <= nseeds:
                val = sqrt_eigval_a[j] * eigvec_a[i, j]
            else:
                # - original, unoptimised code-block
                # numerator = 0.0
                # for p in xrange(nseeds):
                #     numerator += matrix_b[p, i-nseeds] * eigvec_a[p, j]
                # val = numerator / sqrt_eigval_a[j]
                #
                # - optimisation attempt: actually slower
                # numerator = sum(matrix_b[p, i-nseeds] * eigvec_a[p, j] 
                #                 for p in xrange(nseeds))
                #
                # - slightly optimised version: twice as fast on a seedmat of
                #   size 750 x 3000
                #a_mb = array([matrix_b[p, i-nseeds] for p in xrange(nseeds)])
                #a_eva = array([eigvec_a[p, j] for p in xrange(nseeds)])
                #val = (a_mb*a_eva).sum() / sqrt_eigval_a[j]
                #
                # - optmisation suggested by daniel:
                # 100fold increase on a seedmat of size 750 x 3000
                numerator = (matrix_b[:nseeds, i-nseeds] * eigvec_a[:nseeds, j]).sum()
                val = numerator / sqrt_eigval_a[j]

            result[i, j] = val

    if PRINT_TIMINGS:
        print("TIMING(%s): Actual MDS approximation took %f CPU secs" % 
              (__name__, time.clock() - t0))

    return result


"""
=================

Nystrom example implamentation:

    Fast computation of an approximate MDS mapping / PCoA of an (yet unknown) full distance
    matrix. Returned MDS coordinates have the shape num_objects x dim.

    Arguments:
    - `num_objects`:
       total number of objects to compute mapping for.
    - `num_seeds`:
       number of seeds objects. high means more exact solution, but the slower 
    - `dim`:
       dimensionality of MDS mapping
    - `dist_func`:
       callable distance function. arguments should be i,j, with index
       range 0..num_objects-1
    - `permute_order`:
       permute order of objects. recommended to avoid caveeats with
       ordered data that might lead to distorted results. permutation
       is random. run several times for benchmarking.
       
def nystrom_frontend(num_objects, num_seeds, dim, dist_func,
                     permute_order=True):
    
    (seed_distmat, restore_idxs) = build_seed_matrix(
        num_objects, num_seeds, dist_func, permute_order)

    #picked_seeds = argsort(restore_idxs)[:num_seeds]
    
    mds_coords = nystrom(seed_distmat, dim)

    # restoring original order in mds_coords, which has been
    # altered during seed matrix calculation
    return mds_coords[restore_idxs]

=================


SCMDS example implamentation:

    Fast MDS approxmiation SCMDS. Breaks (unknown) distance matrix
    into smaller chunks (tiles), computes MDS solutions for each of
    these and joins them to one form a full approximatiom.


    Arguments:
    - `num_objects`:
      number of objects in distance matrix
    - `tile_size`:
      size of tiles/submatrices. the bigger, the slower but the better
      the approximation
    - `tile_overlap`:
      overlap of tiles. has to be bigger than dimensionality
    - `dim`:
      requested dimensionality of MDS approximation
    - `dist_func`:
      distance function to compute distance between two objects x and
      y. valid index range for x and y should be 0..num_objects-1
    - `permute_order`:
      permute input order if True. reduces distortion. order of
      returned coordinates is kept fixed in either case.
      
def scmds_frontend(num_objects, tile_size, tile_overlap,
                   dim, dist_func, permute_order=True):
    if num_objects < tile_size:
        raise ValueError, \
            "Number of objects cannot be smaller than tile size"
    if tile_overlap > tile_size:
        raise ValueError, \
            "Tile overlap cannot be bigger than tile size"

    if dim > tile_overlap:
        raise ValueError, \
            "Tile overlap must be at least as big as requested dimensionality"
    if not callable(dist_func):
        raise ValueError, "distance getter function not callable"
    
    t0_overall = time.clock()

    if permute_order:
        order = sample(range(num_objects), num_objects)
    else:
        order = range(num_objects)

    ntiles = floor((num_objects - tile_size) / \
                         (tile_size - tile_overlap))+1
    assert ntiles != 0, "Internal error: can't use 0 tiles!"

    # loop over all ntiles, overlapping tiles. apply mds to each
    # single one and join the solutions to the growing overall
    # solution
    #
    tile_no = 1
    tile_start = 0
    tile_end = tile_size + \
               ((num_objects-tile_size) % (tile_size-tile_overlap))
    comb_mds = CombineMds()
    while tile_end <= num_objects:
        # beware: tile_size is just the ideal tile_size, i.e.
        # tile_end-tile_start might not be the same, especially for
        # the last tile
        this_tile_size = tile_end-tile_start

        # construct a tile (submatrix)
        tile = zeros((this_tile_size, this_tile_size))
        for i in xrange(this_tile_size):
            for j in xrange(i+1, this_tile_size):
                tile[i, j] = dist_func(order[i+tile_start],
                                       order[j+tile_start])
                tile[j, i] = tile[i, j]
                
        #print("INFO: Working on tile with idxs %d:%d gives tile shape %d:%d" % \
        #      (tile_start, tile_end, tile.shape[0], tile.shape[1]))

        # Apply MDS on this tile
        #
        t0 = time.clock()
        (tile_eigvecs, tile_eigvals) = cmds_tzeng(tile, dim)
        if PRINT_TIMINGS:
            print("TIMING(%s): MDS on tile %d took %f CPU secs" % 
                  (__name__, tile_no, time.clock() - t0))
        #
        # (slower) alternative is:
        #
        #(tile_tile_eigvecs, tile_eigval) = qiime_pcoa.pcoa(tile)
        # pcoa computes all dims so cut down
        #tile_tile_eigvecs = tile_tile_eigvecs[:, 0:dim]
        

        # Add MDS solution to the growing overall solution
        #
        t0 = time.clock()
        comb_mds.add(tile_eigvecs, tile_overlap)
        if PRINT_TIMINGS:
            print("TIMING(%s): adding of tile %d (shape %d:%d) took %f CPU secs" % 
                  (__name__, tile_no, 
                   tile_eigvecs.shape[0], tile_eigvecs.shape[1], time.clock() - t0))

        tile_start = tile_end - tile_overlap
        tile_end = tile_end + tile_size - tile_overlap
        tile_no += 1

    restore_idxs = argsort(order)
    result = comb_mds.getFinalMDS()[restore_idxs]
    
    if PRINT_TIMINGS:
        print("TIMING(%s): SCMDS took %f CPU secs" % 
              (__name__, time.clock() - t0_overall))
        
    return result

"""

