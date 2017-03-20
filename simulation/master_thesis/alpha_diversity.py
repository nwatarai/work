import sys
import pandas as pd
import numpy as np
from numpy import inf
from matplotlib import pyplot as plt

a = pd.read_csv(sys.argv[1], delimiter=",", header=0, index_col=0)
a = a.transpose()
ra = a / np.sum(a,axis=0)
m = np.transpose(ra.values)

def shannon(array):
	log = np.log(array)
	log[log == -inf] = 0
	values = log * array
	return np.sum(values, axis=1) * (-1)

S = shannon(m)
X = np.arange(S.shape[0]) * 50
plt.plot(X,S)
plt.savefig(sys.argv[1]+".shannon.pdf")