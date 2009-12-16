cimport cython
cimport numpy as np
import numpy as np
from numpy cimport npy_intp
from cogent.maths.spatial.ckd3 cimport kdpoint, points, kdnode, build_tree, rn
from stdlib cimport malloc, free

__version__ = "('1', '4')"

cdef extern from "numpy/arrayobject.h":
#    cdef object PyArray_SimpleNewFromData(int nd, npy_intp *dims,\
#                                            int typenum, void *data)
    cdef void import_array()
#    cdef enum requirements:
#        NPY_OWNDATA

cdef extern from "math.h":
    float pi "M_PI"  # Aliasing the C define to "pi"

@cython.boundscheck(False)
def asa_loop(np.ndarray[DTYPE_t, ndim =2] qcoords, np.ndarray[DTYPE_t, ndim =2] lcoords,\
             np.ndarray[DTYPE_t, ndim =1] qradii,  np.ndarray[DTYPE_t, ndim =1] lradii,\
             np.ndarray[DTYPE_t, ndim =2] spoints, np.ndarray[DTYPE_t, ndim =1] box,\
             DTYPE_t probe, UTYPE_t bucket_size, MAXSYM =200000):

    # looping indexes atom, neighbor, sphere
    cdef int idx, idxn, idxn_skip, idxs, lidx, lidx_c
    # flags
    cdef int is_accessible
    # initialize variables
    cdef UTYPE_t n_acc_points
    cdef DTYPE_t qradius, lradius, search_limit
    cdef DTYPE_t lradius_max = 2.0 + probe
    #malloc'ed
    cdef DTYPE_t *rspoint = <DTYPE_t *>malloc(3 * sizeof(DTYPE_t))
    cdef DTYPE_t *distance = <DTYPE_t *>malloc(3 * sizeof(DTYPE_t))
    cdef DTYPE_t *distance_sq = <DTYPE_t *>malloc(3 * sizeof(DTYPE_t))
    cdef DTYPE_t **dstptr = <DTYPE_t **>malloc(sizeof(DTYPE_t*))
    cdef UTYPE_t **idxptr = <UTYPE_t **>malloc(sizeof(UTYPE_t*))
    # cdef DTYPE_t *areas = <DTYPE_t *>malloc(qcoords.shape[0] * sizeof(DTYPE_t))
    cdef np.ndarray[DTYPE_t, ndim=1] areas = np.ndarray((qcoords.shape[0],), dtype=np.float64)
    
    #c arrays from numpy
    cdef DTYPE_t *qcoords_c = <DTYPE_t *>qcoords.data
    cdef DTYPE_t *lcoords_c = <DTYPE_t *>lcoords.data
    cdef DTYPE_t *spoints_c = <DTYPE_t *>spoints.data
    #datas
    cdef UTYPE_t ssize = spoints.shape[0]
    cdef UTYPE_t lsize = lradii.shape[0]
    cdef npy_intp qpnts = qcoords.shape[0]
    cdef DTYPE_t const_pi = pi * 4.0 / ssize

    #pointers
    cdef UTYPE_t *ridx
    cdef UTYPE_t *ridx_div
        
    # create a temporary array of lattice points, which are within a box around
    # the query atoms. The kd-tree will be constructed from those filterd atoms.
    cdef DTYPE_t *box_c = <DTYPE_t *>box.data
    cdef DTYPE_t *t_ptr                                                      # temporary pointer
    cdef DTYPE_t *t_arr = <DTYPE_t *>malloc(3 * MAXSYM * sizeof(DTYPE_t))    # temporary array of symmetry
    cdef UTYPE_t *t_lid = <UTYPE_t *>malloc(    MAXSYM * sizeof(UTYPE_t))    # maping to original indices
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
    cdef kdpoint search_point
    cdef npy_intp neighbor_number
    cdef kdpoint *kdpnts = points(t_arr, lidx_c, 3)
    cdef kdnode  *tree = build_tree(kdpnts, 0, lidx_c -1, 3, bucket_size, 0)
    
    for 0 <= idx < qpnts:
        qradius = qradii[idx]
        #print qradius, lradius_max
        search_point.coords = qcoords_c + idx*3
        #print search_point.coords[0], search_point.coords[1], search_point.coords[2]
        search_limit = (qradius + lradius_max) * (qradius + lradius_max)
        #print search_limit
        neighbor_number = <npy_intp>rn(tree, kdpnts, search_point, dstptr, idxptr, search_limit, 3, 100)
        #print neighbor_number
        #create array of real indices
        ridx = <UTYPE_t *>malloc(neighbor_number * sizeof(UTYPE_t))
        ridx_div = <UTYPE_t *>malloc(neighbor_number * sizeof(UTYPE_t))

        idxn_skip = 0
        for 0 <= idxn < neighbor_number:
            if dstptr[0][idxn] <= 0.001:
                idxn_skip = 1
                continue
            ridx[idxn - idxn_skip] = kdpnts[idxptr[0][idxn]].index * 3
            ridx_div[idxn - idxn_skip] = t_lid[kdpnts[idxptr[0][idxn]].index] % lsize

        n_acc_points = 0
        for 0 <= idxs < ssize:  # loop over sphere points
            is_accessible = 1
            #if idxs == 0:
            #print search_point.coords[0], search_point.coords[1], search_point.coords[2], '\t',
            #if idxs == 0:
            #print (spoints_c + idxs*3)[0], (spoints_c + idxs*3+1)[0], (spoints_c + idxs*3+2)[0], '\t',
            #if idxs == 0:
            #print qradius, '\t',
            rspoint[0] = search_point.coords[0] + (spoints_c + idxs*3  )[0] * qradius
            rspoint[1] = search_point.coords[1] + (spoints_c + idxs*3+1)[0] * qradius
            rspoint[2] = search_point.coords[2] + (spoints_c + idxs*3+2)[0] * qradius
            #if idxs == 0:
            #print rspoint[0], rspoint[1], rspoint[2]
            #real_point = point * qradius + qcoord
            for 0 <= idxn < neighbor_number - idxn_skip:                                # loop over neighbors
                #if dstptr[0][idxn] == 0.:
                #    continue
                #print '       ', neighbor_number, idxn_skip,idxn, ridx_div[idxn] ,ridx[idxn]
                lradius = lradii[ridx_div[idxn]]
                #print (lcoords_c + ridx[idxn]*3  )[0], (lcoords_c + ridx[idxn]*3+1)[0], (lcoords_c + ridx[idxn]*3+2)[0],
                distance[0] = rspoint[0] - (t_arr + ridx[idxn]  )[0]
                distance[1] = rspoint[1] - (t_arr + ridx[idxn]+1)[0]
                distance[2] = rspoint[2] - (t_arr + ridx[idxn]+2)[0]
                #print '\t', distance[0], distance[1], distance[2],
                distance_sq[0] = distance[0] * distance[0]
                distance_sq[1] = distance[1] * distance[1]
                distance_sq[2] = distance[2] * distance[2]
                #print '\t', distance_sq[0], distance_sq[1], distance_sq[2], lradius * lradius,
                if distance_sq[0] + distance_sq[1] + distance_sq[2] < lradius * lradius:
                    is_accessible = 0
                    #print 'NA'
                    break
                #print
            if is_accessible == 1:
                n_acc_points += 1
        #print n_acc_points
        areas[idx] = const_pi * n_acc_points * qradius * qradius
        free(ridx)
        free(ridx_div)
    free(t_arr)
    free(t_lid)
    free(distance)
    free(distance_sq)
    free(rspoint)
    free(idxptr[0])
    free(idxptr)
    free(dstptr)
    # import_array()
    # cdef np.ndarray output = PyArray_SimpleNewFromData(1, &qpnts, NPY_DOUBLE, <void*>areas)
    # output.flags = output.flags|(NPY_OWNDATA) # this sets the ownership bit
    # PyArray_FLAGS(&output)
    # print PyArray_FLAGS(output) |(NPY_OWNDATA)
    # output = PyArray_NewCopy(output, NPY_CORDER)
    return areas
