import numpy

from numba import njit


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2019, The Cogent Project"
__credits__ = ["Peter Maxwell", "Rob Knight", "Gavin Huttley", "Stephen Ma"]
__license__ = "BSD-3"
__version__ = "2023.2.12a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

# turn off code coverage as njit-ted code not accessible to coverage


@njit(cache=True)
def sum_input_likelihoods(child_indexes, result, likelihoods):  # pragma: no cover
    C = child_indexes.shape[0]
    for child in range(C):
        index = child_indexes[child]
        plhs = likelihoods[child]
        index_height = index.shape[0]
        result_width = result.shape[1]
        if child == 0:
            for parent_col in range(index_height):
                child_col = index[parent_col]
                for motif in range(result_width):
                    result[parent_col, motif] = plhs[child_col, motif]
        else:
            for parent_col in range(index_height):
                child_col = index[parent_col]
                for motif in range(result_width):
                    result[parent_col, motif] *= plhs[child_col, motif]
    return result


@njit(cache=True)
def inner_product(input_likelihoods, mprobs):  # pragma: no cover
    res = 0.0
    for i in range(len(mprobs)):
        res += input_likelihoods[i] * mprobs[i]
    return res


@njit(cache=True)
def get_log_sum_across_sites(lhs, counts):  # pragma: no cover
    log_lhs = numpy.log(lhs)
    res = 0.0
    for i in range(len(counts)):
        res += log_lhs[i] * counts[i]
    return res
