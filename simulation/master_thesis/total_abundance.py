import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

argvs = sys.argv
if len(argvs)<2:
	print "usage: $python %s result.spc/chm.(press).csv" % argvs[0]
	quit()

table = pd.read_csv(argvs[1],sep=',',header=0,index_col=0)

S = np.log10(np.sum(table.values, axis=1))
X = np.arange(S.shape[0]) * 50

plt.plot(X,S)
plt.savefig(argvs[1]+".total.pdf")