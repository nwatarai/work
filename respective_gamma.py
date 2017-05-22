import np
from matplotlib import pyplot
import sys
import panda as pd

from gamma_estimation import gamma_mixture, Gamma, BIC

if __name__ == '__main__':
    mtx = pd.read_csv(sys.argv[1], header=0, index_col=0, sep="\t")
    values = mtx.values
    gm = gamma_mixture(2)


    for i in range(1):
        gm.fit(values[i].T)
        weights = gm.weights
        alpha = gm.alpha
        lamda = gm.lamda

        seq = np.arange(0, 50, 0.02).astype(np.float64)
        print BIC(values[i].T)
        for i in range(degenerate_n_component):
            line = Gamma(seq, a[0][i], b[0][i], w[i])
            pyplot.plot(seq, line, linewidth = 1.0, color="g")
        S += line
        pyplot.plot(seq, S, linewidth=4.0, color="r")
        pyplot.hist(data.T, normed=True, bins=10)
        pyplot.show()
