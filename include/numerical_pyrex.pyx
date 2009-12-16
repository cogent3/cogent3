# Array checking functions for Pyrex code that takes NumPy arrays
#
# The idea here is that many array functions:
#  - Need to know the dimensions of their array inputs.
#  - Require some dimensions of different arrays to match.
#  - Don't ever need to consider a dimension of length 0.
#
# eg:  x = y = z = 0
#      dataA = checkArrayDouble2D(A, &x, &y)
#      dataB = checkArrayDouble2D(B, &z, &x)  # x must match
#

__version__ = "('1', '4')"

cdef extern from "Python.h":
    void *PyCObject_AsVoidPtr(object)


cdef extern from "array_interface.h":
    struct PyArrayInterface:
        int version
        int nd
        char typekind
        int itemsize
        int flags
        int *shape, *strides
        void *data
    int CONTIGUOUS
    
ctypedef object ArrayType

cdef double *uncheckedArrayDouble(ArrayType A):
    cdef PyArrayInterface *a
    cobj = A.__array_struct__
    a = <PyArrayInterface *> PyCObject_AsVoidPtr(cobj)
    return <double *> a.data
    
cdef void *checkArray(ArrayType A, char typecode, int itemsize, 
        int nd, int **dims) except NULL:
    cdef PyArrayInterface *a
    cdef int length, size
    cdef char kind
    if A is None:
        raise TypeError("Array required, got None")
    cobj = A.__array_struct__
    a = <PyArrayInterface *> PyCObject_AsVoidPtr(cobj)
    if a.version != 2:
        raise ValueError(
            "Unexpected array interface version %s" % str(a.version))
    cdef char typecode2
    typecode2 = a.typekind
    if typecode2 != typecode:
        raise TypeError("'%s' type array required, got '%s'" %
                (chr(typecode), chr(typecode2)))
    if a.itemsize != itemsize:
        raise TypeError("'%s%s' type array required, got '%s%s'" %
                (chr(typecode), itemsize, chr(typecode2), a.itemsize))
    if a.nd != nd:
        raise ValueError("%s dimensional array required, got %s" %
                (nd, a.nd))
    if not a.flags & CONTIGUOUS:
        raise ValueError ('Noncontiguous array')
        
    cdef int dimension, val
    cdef int *var
    for dimension from 0 <= dimension < nd:
        val = a.shape[dimension]
        var = dims[dimension]
        if var[0] == 0:
            # Length unspecified, take it from the provided array
            var[0] = val
        elif var[0] != val:
            # Length already specified, but not the same
            raise ValueError("Dimension %s is %s, expected %s" %
                    (dimension, val, var[0]))
        else:
            # Length matches what was expected
            pass
    return a.data
    
    
cdef void *checkArray1D(ArrayType a, char typecode, int size, 
        int *x) except NULL:
    cdef int *dims[1]
    dims[0] = x
    return checkArray(a, typecode, size, 1, dims)
            
cdef void *checkArray2D(ArrayType a, char typecode, int size, 
        int *x, int *y) except NULL:
    cdef int *dims[2]
    dims[0] = x
    dims[1] = y
    return checkArray(a, typecode, size, 2, dims)

cdef void *checkArray3D(ArrayType a, char typecode, int size, 
        int *x, int *y, int *z) except NULL:
    cdef int *dims[3]
    dims[0] = x
    dims[1] = y
    dims[2] = z
    return checkArray(a, typecode, size, 3, dims)
            
cdef void *checkArray4D(ArrayType a, char typecode, int size, 
        int *w, int *x, int *y, int *z) except NULL:
    cdef int *dims[4]
    dims[0] = w
    dims[1] = x
    dims[2] = y
    dims[3] = z
    return checkArray(a, typecode, size, 4, dims)

    
cdef double * checkArrayDouble1D(ArrayType a, int *x) except NULL:
    return <double *> checkArray1D(a, c'f', sizeof(double), x)
    
cdef double * checkArrayDouble2D(ArrayType a, int *x, int *y) except NULL:
    return <double *> checkArray2D(a, c'f', sizeof(double), x, y)
    
cdef double * checkArrayDouble3D(ArrayType a, int *x, int *y, int *z) except NULL:
    return <double *> checkArray3D(a, c'f', sizeof(double), x, y, z)
    
cdef double * checkArrayDouble4D(ArrayType a, int *w, int *x, int *y, int *z) except NULL:
    return <double *> checkArray4D(a, c'f', sizeof(double), w, x, y, z)

    
cdef long * checkArrayLong1D(ArrayType a, int *x) except NULL:
    return <long *> checkArray1D(a, c'i', sizeof(long), x)

cdef long * checkArrayLong2D(ArrayType a, int *x, int *y) except NULL:
    return <long *> checkArray2D(a, c'i', sizeof(long), x, y)
            
cdef long * checkArrayLong3D(ArrayType a, int *x, int *y, int *z) except NULL:
    return <long *> checkArray3D(a, c'i', sizeof(long), x, y, z)
    
cdef long * checkArrayLong4D(ArrayType a, int *w, int *x, int *y, int *z) except NULL:
    return <long *> checkArray4D(a, c'i', sizeof(long), w, x, y, z)
