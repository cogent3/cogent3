cimport numpy as np
ctypedef np.npy_float64 DTYPE_t
ctypedef np.npy_uint64  UTYPE_t

cdef enum constants:
    NSTACK = 100

cdef struct kdpoint:
    UTYPE_t index
    DTYPE_t *coords

cdef struct kdnode:
    UTYPE_t bucket                 # 1 if leaf-bucket, 0 if node
    int dimension
    DTYPE_t position
    UTYPE_t start, end             # start and end index of data points
    kdnode *left, *right           # pointers to left and right nodes

cdef inline void swap(kdpoint*, kdpoint*)

cdef kdpoint *points(DTYPE_t*, UTYPE_t, UTYPE_t)

cdef inline DTYPE_t dist(kdpoint*, kdpoint*, UTYPE_t)

cdef void qsort(kdpoint*, UTYPE_t, UTYPE_t, UTYPE_t)

cdef kdnode *build_tree(kdpoint*, UTYPE_t, UTYPE_t, UTYPE_t, UTYPE_t, UTYPE_t)

cdef UTYPE_t rn(kdnode*, kdpoint*, kdpoint, DTYPE_t**,UTYPE_t**, DTYPE_t, UTYPE_t, UTYPE_t)
                
cdef void *knn(kdnode*, kdpoint*, kdpoint, DTYPE_t*, UTYPE_t*, UTYPE_t, UTYPE_t)

