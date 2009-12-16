cimport cython
import numpy as np
cimport numpy as np
from numpy cimport npy_intp
from cogent.maths.spatial.ckd3 cimport kdpoint, points, kdnode, build_tree, rn
from stdlib cimport malloc, free

__version__ = "('1', '4')"

cdef extern from "numpy/arrayobject.h":
#    cdef object PyArray_SimpleNewFromData(int nd, npy_intp *dims,\
#                                          int typenum, void *data)
    cdef void import_array()
#    cdef enum requirements:
#        NPY_OWNDATA

def cnt_loop(   np.ndarray[DTYPE_t, ndim =2] qcoords,\
                np.ndarray[DTYPE_t, ndim =2] lcoords,\
                np.ndarray[LTYPE_t, ndim =1] qc,\
                np.ndarray[LTYPE_t, ndim =1] lc,\
                UTYPE_t shape1,\
                UTYPE_t shape2,\
                UTYPE_t zero_tra,\
                UTYPE_t mode,\
                DTYPE_t search_limit,\
                np.ndarray[DTYPE_t, ndim =1] box,\
                UTYPE_t bucket_size =10,\
                UTYPE_t MAXSYM =200000,\
                npy_intp MAXCNT =100000):
    #const
    cdef UTYPE_t asu_atoms = shape1 * shape2
    search_limit = search_limit * search_limit

    #looping indexes query atom, lattice atom, neighbor, result
    cdef int idx, lidx, lidx_c, idxn, idxc

    #c arrays from numpy
    cdef DTYPE_t *qcoords_c = <DTYPE_t *>qcoords.data
    cdef DTYPE_t *lcoords_c = <DTYPE_t *>lcoords.data
    cdef DTYPE_t *box_c = <DTYPE_t *>box.data

    #malloc'ed pointers
    cdef UTYPE_t **idxptr = <UTYPE_t **>malloc(sizeof(UTYPE_t*))
    cdef DTYPE_t **dstptr = <DTYPE_t **>malloc(sizeof(DTYPE_t*))

    # temp
    cdef DTYPE_t *t_ptr # temporary pointer
    cdef UTYPE_t  t_idx # index
    cdef UTYPE_t  t_asu # reduced index
    cdef UTYPE_t  t_sym # symmetry
    cdef UTYPE_t  t_tra # translation UTYPE_t
    cdef DTYPE_t  t_dst # distance
    cdef DTYPE_t *t_arr = <DTYPE_t *>malloc(3 * MAXSYM * sizeof(DTYPE_t))    # temporary array of symmetry
    cdef UTYPE_t *t_lid = <UTYPE_t *>malloc(    MAXSYM * sizeof(UTYPE_t))    # maping to original indices

    # result
    #cdef UTYPE_t *c_src = <UTYPE_t *>malloc(MAXCNT * sizeof(UTYPE_t))  # source indices
    #cdef UTYPE_t *c_asu = <UTYPE_t *>malloc(MAXCNT * sizeof(UTYPE_t))  # target indices
    #cdef UTYPE_t *c_sym = <UTYPE_t *>malloc(MAXCNT * sizeof(UTYPE_t))  # symmetries
    #cdef UTYPE_t *c_tra = <UTYPE_t *>malloc(MAXCNT * sizeof(UTYPE_t))  # translations
    #cdef DTYPE_t *c_dst = <DTYPE_t *>malloc(MAXCNT * sizeof(DTYPE_t))  # distances
    
    cdef np.ndarray[UTYPE_t, ndim=1] c_src = np.ndarray((MAXCNT,), dtype=np.uint64)
    cdef np.ndarray[UTYPE_t, ndim=1] c_asu = np.ndarray((MAXCNT,), dtype=np.uint64)
    cdef np.ndarray[UTYPE_t, ndim=1] c_sym = np.ndarray((MAXCNT,), dtype=np.uint64)
    cdef np.ndarray[UTYPE_t, ndim=1] c_tra = np.ndarray((MAXCNT,), dtype=np.uint64)
    cdef np.ndarray[DTYPE_t, ndim=1] c_dst = np.ndarray((MAXCNT,), dtype=np.float64)

    # create a temporary array of lattice points, which are within a box around
    # the query atoms. The kd-tree will be constructed from those filterd atoms.
    lidx_c = 0
    for 0 <= lidx < lcoords.shape[0]:
        t_ptr = lcoords_c + lidx * 3
        if  box_c[0] <= (t_ptr  )[0] <= box_c[3] and\
            box_c[1] <= (t_ptr+1)[0] <= box_c[4] and\
            box_c[2] <= (t_ptr+2)[0] <= box_c[5]:
            t_arr[3*lidx_c  ] = (t_ptr  )[0]
            t_arr[3*lidx_c+1] = (t_ptr+1)[0]
            t_arr[3*lidx_c+2] = (t_ptr+2)[0]
            t_lid[lidx_c] = lidx
            lidx_c += 1

    #make kd-tree
    cdef kdpoint  search_point
    cdef npy_intp neighbor_number
    cdef kdpoint *kdpnts = points(t_arr, lidx_c, 3)
    cdef kdnode  *tree = build_tree(kdpnts, 0, lidx_c - 1, 3, bucket_size, 0)

    idxc = 0
    # loop over every query atom
    for 0 <= idx < qcoords.shape[0]:
        search_point.coords = qcoords_c + idx*3
        neighbor_number = <npy_intp>rn(tree, kdpnts, search_point, dstptr, idxptr, search_limit, 3, 100)
        # loop over all neighbors
        for 0 <= idxn < neighbor_number:
            t_dst =  dstptr[0][idxn]        # the distance of the neighbor to the query
            if t_dst <= 0.001:              # its the same atom, skipping.
                continue

            t_idx =  kdpnts[idxptr[0][idxn]].index     # real index in t_arr array
            t_idx =  t_lid[t_idx]                      # real index in lcoords array
            t_asu =  t_idx % shape2                    # 0 .. N -1, atom number
            t_sym = (t_idx // shape2) % shape1         # 0 .. MXS -1, symmetry number
            t_tra =  t_idx // asu_atoms                # 0 .. (2n + 1)^2 -1, translation number

            if t_tra == zero_tra:  # same unit cell
                if (mode == 0):
                    continue
                elif (mode == 1) and (t_sym == 0): # same asymmetric unit
                    continue
                elif (mode == 2) and (qc[idx] == lc[t_asu]) and (t_sym == 0): #
                    continue

            # safe valid contact
            c_src[idxc] = idx
            c_asu[idxc] = t_asu
            c_sym[idxc] = t_sym
            c_tra[idxc] = t_tra
            c_dst[idxc] = t_dst
            idxc += 1

    free(t_arr)
    free(t_lid)
    # free KD-Tree
    free(idxptr[0])
    free(idxptr)
    free(dstptr)

    # numpy output
    #import_array()
    #cdef np.ndarray n_src = PyArray_SimpleNewFromData(1, &MAXCNT, NPY_UINT,   <void*>c_src)
    # n_src.flags = <npy_intp>n_src.flags|(NPY_OWNDATA) # this sets the ownership bit
    #cdef np.ndarray n_asu = PyArray_SimpleNewFromData(1, &MAXCNT, NPY_UINT,   <void*>c_asu)
    # n_asu.flags = <int>n_asu.flags|(NPY_OWNDATA) # this sets the ownership bit
    #cdef np.ndarray n_sym = PyArray_SimpleNewFromData(1, &MAXCNT, NPY_UINT,   <void*>c_sym)
    # n_sym.flags = <int>n_sym.flags|(NPY_OWNDATA) # this sets the ownership bit
    #cdef np.ndarray n_tra = PyArray_SimpleNewFromData(1, &MAXCNT, NPY_UINT,   <void*>c_tra)
    # n_tra.flags = <int>n_tra.flags|(NPY_OWNDATA) # this sets the ownership bit
    #cdef np.ndarray n_dst = PyArray_SimpleNewFromData(1, &MAXCNT, NPY_DOUBLE, <void*>c_dst)
    # n_dst.flags = <int>n_dst.flags|(NPY_OWNDATA) # this sets the ownership bit
    return (idxc, c_src, c_asu, c_sym, c_tra, c_dst)
