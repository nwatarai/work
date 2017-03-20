import sys
import pandas as pd

data = pd.read_csv(sys.argv[1], header=0, index_col=0)
out = open(sys.argv[1]+".sif", "w")
for i in range(data.shape[0]):
	out.writelines(str(data.values[i][0]) +" to "+ str(data.values[i][1]) + "\n" )
out.close()