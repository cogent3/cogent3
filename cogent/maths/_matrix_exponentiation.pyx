include "../../include/numerical_pyrex.pyx"

version_info = (1, 2)
__version__ = "('1', '4')"

cdef extern from "math.h":
    double exp(double)
    
cdef extern from "matrix_invert.c":
    int matinv (double x[], int n, int m, double space[])

cdef extern from "eigen.c":
    int eigen(int job, double A[], int n, double rr[], double ri[],
        double vr[], double vi[], double w[])
    
    cdef struct complex:
        double re
        double im
        
    int cmatinv(complex x[], int n, int m, double space[])


cdef object empty, zeros, contiguous, FLOAT, COMPLEX

def setNumPy(Numeric):
    global empty, zeros, contiguous, FLOAT, COMPLEX, inner
    FLOAT = 'd'
    COMPLEX = 'D'
    contiguous = Numeric.array
    zeros = Numeric.zeros
    empty = Numeric.empty
    inner = Numeric.inner

def inverse(ArrayType U):
    cdef ArrayType V
    cdef int n, ret, is_complex, i
    cdef double work[128], *data
    V = U.copy()     #  PY
    n = 0
    try:
        data = checkArrayDouble2D(V, &(n), &(n))
        is_complex = 0
    except TypeError:
        data = <double *> checkArray2D(V, c'c', sizeof(double)*2, &(n), &(n))
        is_complex = 1
    
    if n > 64:
        raise ValueError('%s > 64' % n)
    
    # required?
    for i from 0 <= i < 2*n:
        work[i] = 0
    
    if is_complex:
        ret = cmatinv(<complex *> data, n, n, work)
    else:
        ret = matinv(data, n, n, work)
    
    if ret:
        raise ArithmeticError("Determinant too small")
    return V

def eigenvectors(ArrayType Q):
    cdef ArrayType U, Ui, E, R, Ri, Q2
    cdef int n, ret, n2
    cdef double  *U_data, *R_data, *Q_data
    cdef double Q2_data[4096], U2_data[4096], Ui2_data[4096]
    cdef double R2_data[64], Ri2_data[64], work[128]
    n = 0
    Q_data = checkArrayDouble2D(Q, &n, &n)
    if n > 64:
        raise ValueError('%s > 64' % n)
    for i from 0 <= i < (n*n):
        Q2_data[i] = Q_data[i]
    ret = eigen(1, Q2_data, n, R2_data, Ri2_data, U2_data, Ui2_data, work)
    if ret == -1:
        raise RuntimeError, "eigenvalues didn't converge"
    if ret == 1:
        R = empty([n], COMPLEX)   # PY
        U = empty([n,n], COMPLEX) # PY
        U_data = <double *> checkArray2D(U, c'c', sizeof(double)*2, &n, &n)
        R_data = <double *> checkArray1D(R, c'c', sizeof(double)*2, &n)
        
        for i from 0 <= i < n: 
            R_data[i*2+0] = R2_data[i]
            R_data[i*2+1] = Ri2_data[i]
            for j from 0 <= j < n: 
                U_data[(i*n+j)*2+0] = U2_data[j*n+i]
                U_data[(i*n+j)*2+1] = Ui2_data[j*n+i]
    else:
        R = empty([n], FLOAT)   # PY
        U = empty([n,n], FLOAT) # PY
        U_data = checkArrayDouble2D(U, &n, &n)
        R_data = checkArrayDouble1D(R, &n)
        for i from 0 <= i < n: 
            R_data[i] = R2_data[i]
            for j from 0 <= j < n: 
                U_data[i*n+j] = U2_data[j*n+i]
    return (R, U)
    

cdef class EigenExponentiator:
    # Only works with real eigenvalues
    
    cdef readonly int n
    cdef readonly object shape
    cdef readonly ArrayType Q, roots, ev, evT, evI
    cdef double *_ev, *_evI, *_roots
    
    def __init__(self, ArrayType Q, ArrayType roots, ArrayType ev, ArrayType evI):
        cdef int n
        n = 0
        self.roots = roots
        self.ev = ev
        self.evI = evI
        self.Q = Q
        self._roots = checkArrayDouble1D(self.roots, &n)
        self._ev = checkArrayDouble2D(self.ev, &n, &n)
        try:
            self._evI = checkArrayDouble2D(self.evI, &n, &n)
        except ValueError:
            self.evI = contiguous(self.evI)
            self._evI = checkArrayDouble2D(self.evI, &n, &n)            
        self.n = n
        self.shape = (self.n, self.n)
    
    def __call__(self, double t):
        cdef int i, j, k, n
        cdef double expt, uexpt, *data, *roots, *ev, *evI
        cdef ArrayType P
        if t < 0.0:
            raise ValueError('Negative distance %s' % t)
        n = self.n
        roots = self._roots
        ev = self._ev
        evI = self._evI
        P = zeros(self.shape, FLOAT)
        data = checkArrayDouble2D(P, &n, &n)
        for k from 0 <= k < n:
            expt = exp(t * roots[k])
            for i from 0 <= i < n:
                uexpt = evI[i*n+k] * expt
                for j from 0 <= j < n:
                    data[j*n+i] += uexpt * ev[k*n+j]
        for i from 0 <= i < n*n:
            if data[i] < 0.0:
                data[i] = 0.0
        return P
