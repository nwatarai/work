import numpy as np
from matplotlib import pyplot
import sys
import pandas as pd

from gamma_estimation import gamma_mixture, Gamma, BIC

def scaling(X):
    if X.sum() != 0:
        _X = X / X.sum(axis=0)
        assert _X.shape == X.shape
        return _X
    return np.zeros(X.shape)

if __name__ == '__main__':
    mtx = pd.read_csv(sys.argv[1], header=0, index_col=0, sep="\t", comment="#")
    values = mtx.values.astype(np.float64)
    gm = gamma_mixture(n_component=2, scale=True)

    for i in range(100):
        v = np.array([values[i]])
        gm.fit(v.T)
        weights = gm.weights
        alpha = gm.alpha
        lamda = gm.lamda

        """
        scaled_X = gm.scaled_X
        seq = np.arange(0, 1, 0.02).astype(np.float64)
        S = np.zeros(seq.shape)
        print weights
        print lamda
        print alpha
        print BIC(scaled_X, alpha, lamda, weights)
        for i in range(weights.shape[0]):
            line = Gamma(seq, alpha[0][i], lamda[0][i], weights[i])
            pyplot.plot(seq, line, linewidth = 1.0, color="g")
        S += line
        pyplot.plot(seq, S, linewidth=4.0, color="r")
        pyplot.hist(scaled_X, normed=True, bins=10)"""
        if weights.shape[0] == 2:
            print mtx.index[i]
