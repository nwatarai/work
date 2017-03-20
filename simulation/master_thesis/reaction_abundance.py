import sys
import pandas as pd
import numpy as np

argvs = sys.argv
if len(argvs) < 3:
	sys.exit("Usage: $ python %s result metabolism" % argvs[0])

result = pd.read_csv(argvs[1], header=0, index_col=0)
metabolism = pd.read_csv(argvs[2], header=0, index_col=0)
metabolism.values[metabolism.values>0] = 1

output_values = np.dot(result.values, metabolism.values)

output = pd.DataFrame(data=output_values, index=result.index, columns=["reaction_"+str(i) for i in metabolism.columns.tolist()])
output.to_csv(argvs[1]+".reactions.csv")