import sys
import pandas as pd
import numpy as np

argvs = sys.argv
if len(argvs) < 2:
	sys.exit("Usage: python %s result" %argvs[0] )

under_threshold = 0.0000001

table = pd.read_csv(argvs[1], header=0, index_col=0)

#spc
absolute_abundance = np.sum(table.values, axis=1)
peak_value = np.max(absolute_abundance)
peak = np.argmax(absolute_abundance[:absolute_abundance.shape[0] / 2])
index = np.arange(absolute_abundance.shape[0])
t1 = index[absolute_abundance>peak_value/10][0]
t2 = int((peak - t1)/2) +  t1
t3 = peak
t4 = int((peak - t2)*2) + t2
t5 = int((peak - t2)*3) + t2
t9 = table.shape[0]-1
t8 = (t9 - t5)/2 + t5
t7 = (t8 - t5)/2 + t5
t6 = (t7 - t5)/2 + t5
t5 = (t6 - t5)/2 + t5


index = sorted([t1,t2,t3,t4,t5,t6,t7,t8,t9])
index = np.array(map(int,index))

out = table.ix[index, :]
out.values[out.values<under_threshold] = 0.0
out.to_csv(argvs[1]+".experiment_like.csv", sep=",")