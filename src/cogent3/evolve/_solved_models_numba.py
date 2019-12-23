import math

import numpy as np

from numba import njit


version_info = (3, 2)
__version__ = "('2019', '12', '6a')"


@njit(cache=True)
def calc_TN93_P(mprobs, time, alpha1, alpha2, result):

    if not (mprobs.shape[0] == result.shape[0] == result.shape[1] == 4):
        raise ValueError("all array dimensions must equal 4")

    alpha = np.zeros(2)
    pi_star = np.zeros(2)
    mu = np.zeros(2)
    e_mu_t = np.zeros(2)
    transition = np.zeros(2)

    alpha[0] = alpha1
    alpha[1] = alpha2

    pi_star[0] = mprobs[0] + mprobs[1]
    pi_star[1] = mprobs[2] + mprobs[3]

    mu[0] = alpha[0] * pi_star[0] + 1.0 * pi_star[1]
    mu[1] = 1.0 * pi_star[0] + alpha[1] * pi_star[1]

    scale_factor = 0.0
    for motif in range(4):
        i = motif // 2
        other = 1 - i
        scale_factor += (
            alpha[i] * mprobs[2 * i + 1 - motif % 2] + pi_star[other]
        ) * mprobs[motif]

    time /= scale_factor

    e_beta_t = math.exp(-time)
    transversion = 1 - e_beta_t
    for i in range(2):
        other = 1 - i
        e_mu_t[i] = math.exp(-mu[i] * time)
        transition[i] = 1 + (pi_star[other] * e_beta_t - e_mu_t[i]) / pi_star[i]

    for row in range(4):
        i = row // 2
        for column in range(4):
            j = column // 2
            if i == j:
                p = transition[i]
            else:
                p = transversion
            p *= mprobs[column]
            if row == column:
                p += e_mu_t[i]
            result[row, column] = p
