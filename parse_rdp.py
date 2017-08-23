import sys
import numpy as np
import pandas as pd

def parse_rdp(filename):
    mtx = []
    index = []
    for i in open(filename, "r"):
        values = i.rstrip().split("\t")
        _line = []
        index .append(values[0])
        for j, v in enumerate(values[1:]):
            if v == "phylum":
                _line.append(values[j].replace("\"", ""))
            if v == "genus":
                _line.append(values[j].replace("\"", ""))
        mtx.append(_line)
    return pd.DataFrame(data=mtx, index=index)


def extract_otu_by_taxon(df, taxa, rank):
    if rank == "Genus":
        index = 1
    elif rank == "Phylum":
        index = 0

    extract = np.array([False for i in range(df.shape[0])])
    for i in taxa:
        extract += df.values.T[index,:] == i
    return df.ix[extract,:]

if __name__ == "__main__":
    filename = sys.argv[1]
    lad = open(sys.argv[2],"r").read().rstrip().split("\n")
    df = parse_rdp(filename)
    lad_otu = extract_otu_by_taxon(df, lad, "Genus")

    out = open(sys.argv[1] + ".extract.list", "w")
    out.write("\n".join(list(lad_otu.index)))
