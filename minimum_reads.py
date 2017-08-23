import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep="\t", header=0, index_col=0, comment='#')

print np.max(df.sum(axis=0))