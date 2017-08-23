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

    #for i in range(mtx.shape[0]):
    for i in range(110):
        gm = gamma_mixture(n_component=3, scale=True, logarithm=False)
        v = np.array([values[i]])
        gm.fit(v.T)
        if gm.weights != None:
            weights = gm.weights
            alpha = gm.alpha
            lamda = gm.lamda

            #if weights.shape[0] > 1 and gm.converged:
            if True:
                pyplot.figure()
                scaled_X = gm.scaled_X
                seq = np.arange(0, 1, 0.001).astype(np.float64)
                S = np.zeros(seq.shape)
                pyplot.title("BIC: "+str(BIC(scaled_X, alpha, lamda, weights)))
                for j in range(weights.shape[0]):
                    line = Gamma(seq, alpha[0][j], lamda[0][j], weights[j])
                    pyplot.plot(seq, line, linewidth = 2.0, color="g")
                    S += line
                pyplot.plot(seq, S, linewidth=2.0, color="r")
                pyplot.hist(scaled_X, normed=True, bins=25)
                pyplot.savefig(str(i)+"."+mtx.index[i]+".png")
                pyplot.clf()
