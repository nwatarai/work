import sys
import numpy as np
import pandas as pd

def remain_dissimilarity(df, index1, index2):
    #calculate dissimilarity of index1 to index2
    b = df.values[:,index2]
    b[b > 0] = 1
    return 1 - np.sum(df.values[:,index1] * b) / df.values[:,index1].sum().astype(float)

if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], sep="\t", header=0, index_col=0, comment='#')
    out = np.zeros((df.shape[1], df.shape[1]))
    for i in range(df.shape[1]):
        for j in range(df.shape[1]):
            out[i,j] = remain_dissimilarity(df, i, j)
    out_df = pd.DataFrame(data=out, index=df.columns, columns=df.columns)
    out_df.to_csv(sys.argv[1]+".remain_dissimilarity",  sep="\t")
