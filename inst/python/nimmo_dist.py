import numba
import numpy as np
import scipy.stats
from sklearn.metrics import pairwise_distances
import os

q = int(os.environ["NIMMO_EV"])

@numba.njit()
def nimmo(x, y, dims, scales, alpha = None):
    dims = dims[0]
    scales = scales[0]
    if alpha == None:
        alpha = [1] * len(dims)
    start = 0
    results = []
    for i in range(len(dims)):
        result = 0.0
        norm_x = 0.0
        norm_y = 0.0
        curr_x = x[start : start + dims[i]]
        curr_y = y[start : start + dims[i]]
        for j in range(curr_x.shape[0]):
            result += curr_x[j] * curr_y[j]
            norm_x += curr_x[j] ** 2
            norm_y += curr_y[j] ** 2
        result = 1.0 - (result / np.sqrt(norm_x * norm_y))
        results.append(result)
        
        start += int(round(dims[i]))
    final = 0.0
    for i in range(len(results)):
        final += alpha[i] * (scales[i] * results[i]) ** q
    return final ** (1.0 / q)

