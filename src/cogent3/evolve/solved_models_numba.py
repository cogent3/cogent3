import math

import numpy as np
from numba import njit

# turn off code coverage as njit-ted code not accessible to coverage


@njit(cache=True)
def calc_TN93_P(mprobs, time, alpha1, alpha2, result) -> None:  # pragma: no cover
    if not (mprobs.shape[0] == result.shape[0] == result.shape[1] == 4):
        msg = "all array dimensions must equal 4"
        raise ValueError(msg)

    alpha = np.array([alpha1, alpha2])
    pi_star = np.array([mprobs[0] + mprobs[1], mprobs[2] + mprobs[3]])
    mu = np.array(
        [
            alpha1 * pi_star[0] + 1.0 * pi_star[1],
            1.0 * pi_star[0] + alpha2 * pi_star[1],
        ],
    )
    e_mu_t = np.zeros(2)
    transition = np.zeros(2)

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
            p = transition[i] if i == j else transversion
            p *= mprobs[column]
            if row == column:
                p += e_mu_t[i]
            result[row, column] = p
