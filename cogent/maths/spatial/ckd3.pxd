cimport numpy as np
ctypedef np.npy_double DTYPE_t

cdef enum constants:
    NSTACK = 100

cdef struct kdpoint:
    unsigned int index
    DTYPE_t *coords

cdef struct kdnode:
    unsigned int bucket                 # 1 if leaf-bucket, 0 if node
    int dimension
    DTYPE_t position
    unsigned int start, end             # start and end index of data points
    kdnode *left, *right                # pointers to left and right nodes

cdef inline void swap(kdpoint*, kdpoint*)

cdef kdpoint *points(DTYPE_t*, unsigned int, unsigned int)

cdef inline DTYPE_t dist(kdpoint*, kdpoint*, unsigned int)

cdef void qsort(kdpoint*, unsigned int, unsigned int, unsigned int)

cdef kdnode *build_tree(kdpoint*, unsigned int, unsigned int,\
                unsigned int, unsigned int, unsigned int)

cdef unsigned int rn(kdnode*, kdpoint*, kdpoint, DTYPE_t**,\
                unsigned int**, DTYPE_t, unsigned int, unsigned int)
                
cdef void *knn(kdnode*, kdpoint*, kdpoint, DTYPE_t*, unsigned int*, \
               unsigned int, unsigned int)

