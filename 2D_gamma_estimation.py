import numpy as np
from matplotlib import pyplot
import sys
import pandas as pd

from gamma_estimation import gamma_mixture, Gamma, BIC

def generate_data(weight, alpha, lamda, n_sample):
    assert alpha.shape == lamda.shape
    mtx = []
    for k in range(alpha.shape[0]):
        X = np.array([])
        for i in range(alpha.shape[1]):
            _n_sample = int(n_sample * weight[i])
            X = np.append(X, np.random.gamma(alpha[k][i], 1 / lamda[k][i], _n_sample))
        X = np.array(X).ravel()
        mtx.append(X)
    return np.array(mtx).astype(np.float64)


if __name__ == '__main__':
    n_component = 10
    weights = np.array([0.5, 0.3, 0.2]).astype(np.float64)
    alpha = np.array([[3., 10., 4.],
                      [10., 2., 5.]]).astype(np.float64)
    lamda = np.array([[0.1, 0.4, 0.7],
                      [0.7, 0.3, 0.1]]).astype(np.float64)
    data = generate_data(weights, alpha, lamda, 200)
    print data
    gm = gamma_mixture(n_component)
    gm.fit(data.T)

    w = gm.weights
    a = gm.alpha
    b = gm.lamda
    t_w = weights
    t_a = alpha
    t_b = lamda

    print w
    print a
    print b
    print t_w
    print t_a
    print t_b