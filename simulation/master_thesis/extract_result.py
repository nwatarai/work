import sys
import pandas as pd

argvs = sys.argv
if len(argvs) < 3:
	sys.exit("Usage: python %s result nmr_of_rows_you_want(int)" %argvs[0] )

under_threshold = 0.00000001

data = pd.read_csv(argvs[1], header=0, index_col=0)

interval = int(float(data.shape[0]) / float(argvs[2]))
extract = range(0, interval * int(argvs[2]), interval)[1:]
data = data.ix[extract, :]
data[data<under_threshold] = 0.0

data.to_csv(argvs[1]+".extract.csv", data=data)