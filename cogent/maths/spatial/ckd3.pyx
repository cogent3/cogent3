#cython: boundscheck=False
#(not slicing or indexing any numpy arrays)

#TODO:
# - extend to more then 3-dimensions (feature)
# - replace build_tree with non-recursive function (speed)
# - add the option to determine splitting planes based on point position spread
#   (feature)
# - enable bottom-up algorithms by keeping track of ancestral node
# - common stack size constant, how?

cimport numpy as np
from numpy cimport NPY_DOUBLE, NPY_ULONGLONG, npy_intp
from stdlib cimport malloc, realloc, free

__version__ = "('1', '4')"

cdef extern from "numpy/arrayobject.h":
    cdef object PyArray_SimpleNewFromData(int nd, npy_intp *dims,\
                                            int typenum, void *data)
    cdef void import_array()

cdef kdpoint *points(DTYPE_t *c_array, UTYPE_t points, UTYPE_t dims):
    """creates an array of kdpoints from c-array of numpy doubles."""
    cdef kdpoint *pnts = <kdpoint *>malloc(sizeof(kdpoint)*points)
    cdef UTYPE_t i
    for 0 <= i < points:
        pnts[i].index = i
        pnts[i].coords = c_array+i*dims
    return pnts

cdef inline void swap(kdpoint *a, kdpoint *b):
    """swaps two pointers to kdpoint structs."""
    cdef kdpoint t
    t = a[0]
    a[0] = b[0]
    b[0] = t

cdef inline DTYPE_t dist(kdpoint *a, kdpoint *b, UTYPE_t dims):
    """calculates the squared distance between two points."""
    cdef UTYPE_t i
    cdef DTYPE_t dif, dst = 0
    for 0 <= i < dims:
        dif = a.coords[i] - b.coords[i]
        dst += dif * dif
    return dst

cdef void qsort(kdpoint *A, UTYPE_t l, UTYPE_t r, UTYPE_t dim):
    """implements the quick sort algorithm on kdpoint arrays."""
    cdef UTYPE_t i, j, jstack = 0
    cdef DTYPE_t v
    cdef UTYPE_t *istack = <UTYPE_t *>malloc(NSTACK * sizeof(UTYPE_t))
    while True:
        if r - l > 2:
            i = (l + r) >> 1
            if A[l].coords[dim] > A[i].coords[dim]: swap(&A[l], &A[i])
            if A[l].coords[dim] > A[r].coords[dim]: swap(&A[l], &A[r])
            if A[i].coords[dim] > A[r].coords[dim]: swap(&A[i], &A[r])
            j = r - 1
            swap(&A[i], &A[j])
            i = l
            v = A[j].coords[dim]
            while True:
                while A[i+1].coords[dim] < v: i+=1
                i+=1
                while A[j-1].coords[dim] > v: j-=1
                j-=1
                if j < i:
                    break
                swap(&A[i], &A[j])
            swap(&A[i], &A[r-1])
            jstack += 2
            if r - i >= j:
                istack[jstack] = r
                istack[jstack - 1] = i
                r = j
            else:
                istack[jstack] = j
                istack[jstack - 1] = l
                l = i
        else:
            i = (l + r) >> 1
            if A[l].coords[dim] > A[i].coords[dim]: swap(&A[l], &A[i])
            if A[l].coords[dim] > A[r].coords[dim]: swap(&A[l], &A[r])
            if A[i].coords[dim] > A[r].coords[dim]: swap(&A[i], &A[r])
            if jstack == 0:
                break
            r = istack[jstack]
            jstack-=1
            l = istack[jstack]
            jstack-=1
    free(istack)

cdef kdnode *build_tree(kdpoint *point_list, UTYPE_t start, UTYPE_t end,\
                UTYPE_t dims, UTYPE_t bucket_size, UTYPE_t depth):
    """recursive tree building function."""
    # cannot make variable in if/else
    cdef UTYPE_t split, i
    cdef kdnode *node = <kdnode*>malloc(sizeof(kdnode))
    node.dimension = depth % dims
    node.start = start
    node.end = end
    if end - start <= bucket_size:
        # make bucket node
        node.bucket = 1
        node.position = -1.0
        node.left = NULL
        node.right = NULL
    else:
        ## make branch node
        node.bucket = 0
        split = (start + end) / 2
        qsort(point_list, start, end, node.dimension)
        node.position = point_list[split].coords[node.dimension]
        # recurse
        node.left = build_tree(point_list, start, split, dims , bucket_size , depth+1)
        node.right = build_tree(point_list, split+1, end, dims , bucket_size , depth+1)
    return node

cdef void *knn(kdnode *root, kdpoint *point_list, kdpoint point, DTYPE_t *dst,\
                UTYPE_t *idx, UTYPE_t k, UTYPE_t dims):
    """finds the K-Nearest Neighbors."""
    # arrays of pointers will be used as a stack for left and right nodes
    # left nodes will be explored first.
    cdef kdnode *lstack[100]

    cdef UTYPE_t i, j, jold, ia, kmin # counter and index
    cdef DTYPE_t a, i_dist, diff
    cdef kdnode *node

    # set helper variable to heap-queue
    kmin = k - 1

    # initialize stack
    cdef int jstack = 1
    lstack[jstack] = root

    # initialize arrays
    for 0 <= i < k:
        dst[i] = 1000000000.00  # DBL_MAX
        idx[i] = 2147483647     # INT_MAX

    while jstack:
        node = lstack[jstack]
        jstack -= 1
        if node.bucket:
            for node.start <= i <= node.end:
                i_dist = dist(&point_list[i], &point, dims)
                if i_dist < dst[0]:
                    dst[0] = i_dist
                    idx[0] = i
                    if k > 1:
                        a = dst[0]
                        ia = idx[0]
                        jold = 0
                        j = 1
                        while j <= kmin:
                            if (j < kmin) and (dst[j] < dst[j+1]):
                                j+=1
                            if (a >= dst[j]):
                                break
                            dst[jold] = dst[j]
                            idx[jold] = idx[j]
                            jold = j
                            j = 2*j + 1
                        dst[jold] = a
                        idx[jold] = ia
        else:
            diff = point.coords[node.dimension] - node.position
            if diff < 0:
                if dst[0] >= diff * diff:
                    jstack+=1
                    lstack[jstack] = node.right
                jstack+=1
                lstack[jstack] = node.left
            else:
                if dst[0] >= diff * diff:
                    jstack+=1
                    lstack[jstack] = node.left
                jstack+=1
                lstack[jstack] = node.right

cdef UTYPE_t rn(kdnode *root, kdpoint *point_list, kdpoint point, DTYPE_t **dstptr,\
                UTYPE_t **idxptr, DTYPE_t r, UTYPE_t dims, UTYPE_t buf):
    """finds points within radius of query."""
    # arrays of pointers will be used as a stack for left and right nodes
    # left nodes will be explored first.
    cdef kdnode *lstack[100]
    dstptr[0] = <DTYPE_t *>malloc(buf * sizeof(DTYPE_t))
    idxptr[0] = <UTYPE_t *>malloc(buf * sizeof(UTYPE_t))

    cdef UTYPE_t i, count # counter and index
    cdef DTYPE_t i_dist, diff
    cdef kdnode *node

    # initialize stack
    cdef int jstack = 1
    lstack[jstack] = root

    count = 0
    while jstack:
        node = lstack[jstack]
        jstack -= 1
        if node.bucket:
            for node.start <= i <= node.end:
                i_dist = dist(&point_list[i], &point, dims)
                if i_dist < r:
                    dstptr[0][count] = i_dist
                    idxptr[0][count] = i
                    count += 1
                    if count % buf == 0:
                        dstptr[0] = <DTYPE_t *>realloc(dstptr[0], (count + buf) * sizeof(DTYPE_t))
                        idxptr[0] = <UTYPE_t *>realloc(idxptr[0], (count + buf) * sizeof(UTYPE_t))
        else:
            diff = point.coords[node.dimension] - node.position
            if diff < 0:
                if r >= diff * diff:
                    jstack+=1
                    lstack[jstack] = node.right
                jstack+=1
                lstack[jstack] = node.left
            else:
                if r >= diff * diff:
                    jstack+=1
                    lstack[jstack] = node.left
                jstack+=1
                lstack[jstack] = node.right
    dstptr[0] = <DTYPE_t *>realloc(dstptr[0], count * sizeof(DTYPE_t))
    idxptr[0] = <UTYPE_t *>realloc(idxptr[0], count * sizeof(UTYPE_t))
    return count

cdef class KDTree:
    """Implements the KDTree data structure for fast neares neighbor queries."""
    cdef np.ndarray n_array
    cdef DTYPE_t *c_array
    cdef kdpoint *kdpnts
    cdef kdnode *tree
    cdef readonly UTYPE_t dims
    cdef readonly UTYPE_t pnts
    cdef readonly UTYPE_t bucket_size
    def __init__(self, np.ndarray[DTYPE_t, ndim =2] n_array, \
                        UTYPE_t bucket_size =5):
        self.bucket_size = bucket_size
        self.pnts = n_array.shape[0]
        self.dims = n_array.shape[1]
        self.n_array = n_array
        self.c_array = <DTYPE_t *> n_array.data
        self.kdpnts = points(self.c_array, \
                             self.pnts, self.dims)
        self.tree = build_tree(self.kdpnts, 0, self.pnts-1, \
                               self.dims,self.bucket_size,0)
        import_array()

    def knn(self, np.ndarray[DTYPE_t, ndim =1] point, npy_intp k):
        """Finds the K-Nearest Neighbors of given point.
        Arguments:
            - point: 1-d numpy array (query point).
            - k: number of neighbors to find."""
        if self.pnts < k:
            return 1
        cdef UTYPE_t i
        cdef kdpoint pnt
        pnt.coords = <DTYPE_t *>point.data
        cdef UTYPE_t size = point.size
        cdef DTYPE_t *dst = <DTYPE_t *>malloc(k * sizeof(DTYPE_t))
        cdef UTYPE_t *idx = <UTYPE_t *>malloc(k * sizeof(UTYPE_t))
        cdef UTYPE_t *ridx = <UTYPE_t *>malloc(k * sizeof(UTYPE_t))
        knn(self.tree, self.kdpnts, pnt, dst, idx, k, self.dims)
        cdef np.ndarray dist = PyArray_SimpleNewFromData(1, &k, NPY_DOUBLE, <void*>dst)
        for 0 <= i < k:
            ridx[i] = self.kdpnts[idx[i]].index
        cdef np.ndarray index = PyArray_SimpleNewFromData(1, &k, NPY_ULONGLONG, <void*>ridx)
        free(idx)
        return (index, dist)

    def rn(self, np.ndarray[DTYPE_t, ndim =1] point, DTYPE_t r):
        """Returns Radius Neighbors i.e. within radius from query point."""
        cdef UTYPE_t i
        cdef npy_intp j
        cdef kdpoint pnt
        pnt.coords = <DTYPE_t *>point.data
        cdef UTYPE_t size = point.size
        cdef DTYPE_t **dstptr = <DTYPE_t **>malloc(sizeof(DTYPE_t*))
        cdef UTYPE_t **idxptr = <UTYPE_t **>malloc(sizeof(UTYPE_t*))
        j = <npy_intp>rn(self.tree, self.kdpnts, pnt, dstptr, idxptr, r, self.dims, 100)
        cdef np.ndarray dist = PyArray_SimpleNewFromData(1, &j, NPY_DOUBLE, <void*>dstptr[0])
        cdef UTYPE_t *ridx = <UTYPE_t *>malloc(j * sizeof(UTYPE_t))
        for 0 <= i < j:
            ridx[i] = self.kdpnts[idxptr[0][i]].index
        cdef np.ndarray index = PyArray_SimpleNewFromData(1, &j, NPY_ULONGLONG, <void*>ridx)
        free(idxptr[0])
        free(idxptr)
        free(dstptr)
        return (index, dist)

