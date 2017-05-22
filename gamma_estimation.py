import numpy as np
from matplotlib import pyplot
from scipy.stats import gamma as Gamma
from scipy.special import gamma, digamma
import sys

class gamma_mixture(object):
    def __init__(self, n_component, iter_max=1000):
        self.n_component = n_component
        self.weights = None
        self.alpha = None
        self.lamda = None
        self.BIC = None
        self.BIC_transition = []
        self.iter_max = iter_max
        self.weight_threshold = None

    def fit(self, X):
        self.weight_threshold = 1 / float(X.shape[0])
        self.ndim = X.shape[1] 
        self.weights = np.ones(self.n_component) / float(self.n_component)
        self.alpha = np.random.uniform(0.1, 10.0, (self.ndim, self.n_component))
        self.lamda = np.random.uniform(0.1, 10.0, (self.ndim, self.n_component))

        step_size = np.power(np.ones(self.iter_max) / (np.arange(self.iter_max) + 1), 0.2)

        for i in range(self.iter_max):
            params = np.concatenate((self.weights.ravel(), self.alpha.ravel(), self.lamda.ravel()), axis=1)
            resps = self.expectation(X)
            self.maximization(X, resps, step_size[i])
            estimated = np.concatenate((self.weights.ravel(), self.alpha.ravel(), self.lamda.ravel()), axis=1)
            self.cal_BIC(X)
            if np.allclose(params, estimated):
                print "parameters have converged at step_%i" %(i + 1)
                break
            self.degenerate_param()
        else:
            print "parameters may not have converged in %i times" %self.iter_max

    def gamma_pdf(self, X):
        return self.lamda ** self.alpha * X ** (self.alpha - 1) * np.exp(- self.lamda * X) / gamma(self.alpha)

    def expectation(self, X):
        likelihood = self.weights * self.gamma_pdf(X)
        resps = likelihood / likelihood.sum(axis=-1, keepdims=True)
        return resps

    def maximization(self, X, resps, step_size):
        Nk = np.sum(resps, axis=0)
        self.weights = (1. - step_size) * self.weights + step_size * Nk / float(X.shape[0])
        self.lamda = (1. - step_size) * self.lamda + step_size * self.alpha * Nk / np.sum(X * resps, axis=0)
        G = np.sum((np.log(X) + np.log(self.lamda) - digamma(self.alpha)) * resps, axis=0) / float(X.shape[0])
        self.alpha += step_size * G

    def cal_BIC(self, X):
        self.BIC = np.log(X.shape[0]) * (3 * self.n_component - 1) - 2 * np.sum(np.log(self.weights * self.gamma_pdf(X)))
        self.BIC_transition.append(self.BIC)

    def degenerate_param(self):
        index = self.weights < self.weight_threshold
        self.weights = self.weights[~index]
        self.alpha = self.alpha[:,~index]
        self.lamda = self.lamda[:,~index]
        self.weights /= self.weights.sum()

def pseudo_digamma(x):
    _x = x - 0.5
    return np.log(_x) + 0.041666666666 / _x / _x

def generate_data(alpha, lamda, n_sample):
    assert alpha.shape == lamda.shape
    X = np.array([])
    for i in range(alpha.shape[0]):
        X = np.append(X, np.random.gamma(alpha[i], 1 / lamda[i], n_sample))
    X = np.array(X).ravel()
    return np.array([X])

def Gamma(x, alpha, lamda, weights):
    return weights * lamda ** alpha * x ** (alpha - 1) * np.exp(- lamda * x) / gamma(alpha)

def BIC(x, alpha, lamda, weights):
    n_component = alpha.shape[0]
    return np.log(x.shape[0]) * (3 * n_component - 1) - 2 * np.sum(np.log(weights * Gamma(x, alpha, lamda, weights)))

def main(n_component):
    #alpha = np.random.uniform(0.5, 10.0, n_component)
    #lamda = np.random.uniform(0.1, 2.0, n_component)
    weights = np.array([0.4, 0.3, 0.3]).astype(np.float64)
    alpha = np.array([3., 10., 4.]).astype(np.float64)
    lamda = np.array([0.1, 0.4, 0.7]).astype(np.float64)
    data = generate_data(alpha, lamda, 100).astype(np.float64)
    gm = gamma_mixture(n_component)
    gm.fit(data.T)

    w = gm.weights
    a = gm.alpha
    b = gm.lamda
    t_w = weights
    t_a = alpha
    t_b = lamda
    """
    for i in gm.BIC_transition:
        print i
        """
    print w
    print a
    print b
    print t_w
    print t_a
    print t_b

    print gm.BIC
    print BIC(data.T, alpha, lamda, weights)
    degenerate_n_component = a.shape[1]

    seq = np.arange(0, 50, 0.02).astype(np.float64)
    S = np.zeros(seq.shape)
    t_S = np.zeros(seq.shape)
    for i in range(degenerate_n_component):
        line = Gamma(seq, a[0][i], b[0][i], w[i])
        pyplot.plot(seq, line, linewidth = 1.0, color="g")
        S += line
    for i in range(3):
        t_line = Gamma(seq, t_a[i], t_b[i], t_w[i])
        pyplot.plot(seq, t_line, linewidth = 1.0, color="k")
        t_S += t_line 
    pyplot.plot(seq, S, linewidth=4.0, color="r")
    pyplot.plot(seq, t_S, linewidth=3.0, color="y")
    pyplot.hist(data.T, normed=True, bins=70)

    pyplot.show()

if __name__ == '__main__':
    main(10)
