include "../../include/numerical_pyrex.pyx"

cdef extern from "math.h":
    double exp(double)

version_info = (1, 0)

def calc_TN93_P(int do_scaling, mprobs, double time, alpha_1, alpha_2, result):
    cdef int S, motif, i, other, row, column, b_row, b_column
    cdef double *pi, *P, scale_factor
    cdef double pi_star[2], alpha[2], mu[2], e_mu_t[2], e_beta_t
    cdef double transition[2], transversion, p
    
    alpha[0] = alpha_1
    alpha[1] = alpha_2
    
    S = 4
    pi = checkArrayDouble1D(mprobs, &S)
    P = checkArrayDouble2D(result, &S, &S)
    
    pi_star[0] = pi[0] + pi[1]
    pi_star[1] = pi[2] + pi[3]
    
    mu[0] = alpha[0] * pi_star[0] + 1.0 * pi_star[1]
    mu[1] = 1.0 * pi_star[0] + alpha[1] * pi_star[1]
    
    if do_scaling:
        scale_factor = 0.0
        for motif in range(4):
            i = motif // 2
            other = 1 - i
            scale_factor += (alpha[i] * pi[2*i+1-motif%2] + pi_star[other]) * pi[motif]
        time /= scale_factor
    
    e_beta_t = exp(-time)
    transversion = 1 - e_beta_t
    for i in range(2):
        other = 1 - i
        e_mu_t[i] = exp(-mu[i]*time)
        transition[i] = 1 + (pi_star[other] * e_beta_t - e_mu_t[i]) / pi_star[i]
    
    for row in range(4):
        i = row // 2
        for column in range(4):
            j = column // 2
            if i == j:
                p = transition[i]
            else:
                p = transversion
            p *= pi[column]
            if row == column:
                p += e_mu_t[i]
            P[column+4*row] = p
                
