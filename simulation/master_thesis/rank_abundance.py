import sys
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from numpy import inf, nan

# fileopen
argvs = sys.argv
argc = len(argvs)
if (argc < 2):
    print "Usage: # python %s result.extract transpose=False" % argvs[0]
    quit()


d = pd.read_csv(argvs[1], header=0, index_col=0)

matrix = d.values
sample_names = d.index

matrix = np.log10(matrix)
matrix[matrix == -inf] = nan
_min = np.nanmin(matrix)
matrix[np.isnan(matrix)] = _min

__min = np.amin(matrix)
__max = np.amax(matrix)

n = matrix.shape[0]
color = cm.rainbow(np.linspace(0,1,n))
for m, c, sample_name in zip(matrix, color, sample_names):
    sorted_row = sorted(m, reverse=True)
    X = range(1,len(sorted_row)+1)
    Y = sorted_row
    plt.ylim(__min,__max)
    plt.plot(X, Y, color=c, label=sample_name)
plt.savefig(argvs[1]+".rank_abundance.pdf")