import sys
import numpy as np
import pandas as pd
argvs = sys.argv
if len(argvs)<2:
	print "Usage: $python %s target_file"%argvs[0]
	quit()

data = pd.read_csv(argvs[1], header=0, index_col=0)
S = np.sum(data.values, axis=0)

extract = data.ix[:, S>0]

extract.to_csv(argvs[1]+".del0.csv")