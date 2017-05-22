import np
from matplotlib import pyplot
import sys

def create_data(N, K):
    X, mu_star, sigma_star = [], [], []
    for i in range(K):
        loc = (np.random.rand() - 0.5) * 10.0 # range: -5.0 - 5.0
        scale = np.random.rand() * 3.0 # range: 0.0 - 3.0
        X = np.append(X, np.random.normal(loc = loc, scale = scale, size = N / K))
        mu_star = np.append(mu_star, loc)
        sigma_star = np.append(sigma_star, scale)
    return (X, mu_star, sigma_star)

def gaussian(mu, sigma):
    def f(x):
        return np.exp(-0.5 * (x - mu) ** 2 / sigma) / np.sqrt(2 * np.pi * sigma)
    return f

def estimate_posterior_likelihood(X, pi, gf):
    l = np.zeros((X.size, pi.size))
    for (i, x) in enumerate(X):
        l[i, :] = gf(x)
    return pi * l * np.vectorize(lambda y: 1 / y)(l.sum(axis = 1).reshape(-1, 1))

def estimate_gmm_parameter(X, gamma):
    N = gamma.sum(axis = 0)
    mu = (gamma * X.reshape((-1, 1))).sum(axis = 0) / N
    sigma = (gamma * (X.reshape(-1, 1) - mu) ** 2).sum(axis = 0) / N
    pi = N / X.size
    return (mu, sigma, pi)

def calc_Q(X, mu, sigma, pi, gamma):
    Q = (gamma * (np.log(pi * (2 * np.pi * sigma) ** (-0.5)))).sum()
    for (i, x) in enumerate(X):
        Q += (gamma[i, :] * (-0.5 * (x - mu) ** 2 / sigma)).sum()
    return Q

if __name__ == '__main__':

    # data
    K = 2
    N = 1000 * K
    X, mu_star, sigma_star = create_data(N, K)

    # termination condition
    epsilon = 0.000001

    # initialize gmm parameter
    pi = np.random.rand(K)
    mu = np.random.randn(K)
    sigma = np.abs(np.random.randn(K))
    Q = -sys.float_info.max
    delta = None

    # EM algorithm
    while delta == None or delta >= epsilon:
        gf = gaussian(mu, sigma)

        # E step: estimate posterior probability of hidden variable gamma
        gamma = estimate_posterior_likelihood(X, pi, gf)

        # M step: miximize Q function by estimating mu, sigma and pi
        mu, sigma, pi = estimate_gmm_parameter(X, gamma)

        # calculate Q function
        Q_new = calc_Q(X, mu, sigma, pi, gamma)
        delta = Q_new - Q
        Q = Q_new

    # result
    print u'mu*: %s, sigma*: %s' % (str(np.sort(np.around(mu_star, 3))), str(np.sort(np.around(sigma_star, 3))))
    print u'mu : %s, sigma : %s' % (str(np.sort(np.around(mu, 3))), str(np.sort(np.around(sigma, 3))))

    # plot
    n, bins, _ = pyplot.hist(X, 50, normed = 1, alpha = 0.3)
    seq = np.arange(-15, 15, 0.02)
    for i in range(K):
        pyplot.plot(seq, gaussian(mu[i], sigma[i])(seq), linewidth = 2.0)
    pyplot.show()