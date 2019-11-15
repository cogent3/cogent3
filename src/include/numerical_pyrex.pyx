# Array checking functions for Pyrex code that takes NumPy arrays
#
# The idea here is that many array functions:
#  - Need to know the dimensions of their array inputs.
#  - Require some dimensions of different arrays to match.
#  - Don't ever need to consider a dimension of length 0.
#
# eg:  x = y = z = 0
#      checkArray2D(A, &x, &y)
#      checkArray2D(B, &z, &x)  # x must match
#

__version__ = "('2019', '11', '15', 'a')"

cdef extern from "limits.h":
    int INT_MAX
    long LONG_MAX
    
ctypedef fused num:
    # These are the array types used by PyCogent Cython modules
    double
    long
    unsigned char

ctypedef fused dim:
    # Types which may be used as the size of a dimension
    Py_ssize_t
    long
    int
    
ctypedef double[::1] Double1D
ctypedef double[:, ::1] Double2D
ctypedef double[:, :, ::1] Double3D
ctypedef long[::1] Long1D
ctypedef long[:, ::1] Long2D
ctypedef long[:, :, ::1] Long3D


cdef int checkDim(dimension, Py_ssize_t val, dim *var) except 1:
    if var[0] == 0:
        # Length unspecified, take it from the provided array
        if dim is int:
            if val > INT_MAX:
                raise ValueError("%s dimension is %s, too big" % (dimension, val))
            var[0] = <int> val
        elif dim is long:
            if val > LONG_MAX:
                raise ValueError("%s dimension is %s, too big" % (dimension, val))
            var[0] = <long> val 
        else:
            var[0] = val 

    elif var[0] != val:
        # Length already specified, but not the same
        raise ValueError("%s dimension is %s, expected %s" %
                (dimension, val, var[0]))
    else:
        # Length matches what was expected
        pass


cdef int checkArray1D(num[::1] a, dim *x) except 1:
    if a is None:
        raise ValueError('Array required, got None')
    checkDim('1st', a.shape[0], x)

cdef int checkArray2D(num[:, ::1] a, dim *x, dim *y) except 1:
    if a is None:
        raise ValueError('Array required, got None')
    checkDim('1st', a.shape[0], x)
    checkDim('2nd', a.shape[1], y)

cdef int checkArray3D(num[:, :, ::1] a, dim *x, dim *y, dim *z) except 1:
    if a is None:
        raise ValueError('Array required, got None')
    checkDim('1st', a.shape[0], x)
    checkDim('2nd', a.shape[1], y)
    checkDim('3rd', a.shape[2], z)

