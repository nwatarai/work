import pandas as pd
import numpy as np 
from sklearn.cluster import KMeans
import sys
from scipy.spatial.distance import pdist, squareform, euclidean
import argparse

def parser_setting():
    parser = argparse.ArgumentParser(prog='respective_kmeans.py', description='calculate kmeans for each row')
    parser.add_argument('matrix',
                        action='store',
                        type=str,
                        help='File name of matrix')
    parser.add_argument('-o', '--output',
                        action='store',
                        type=str,
                        default=None,
                        help='output directory path')
    parser.add_argument('-m', '--max_cluster',
                        action='store',
                        type=int,
                        default=5,
                        help='max numbr of cluster assumed')
    parser.set_defaults(unweighted=False) 
    args = parser.parse_args()
    return (parser, args)

def _PseudoF(data, labels, k, centroids, SSt):
    n_points = data.shape[0]
    label_set = set(labels)
    SSw = 0
    for label in label_set:
        index = labels == label
        _data = data[index]
        _n_points = _data.shape[0]
        _SSw = 0
        for i in range(_n_points):
            _SSw += euclidean(_data[i], centroids[label])
        SSw += _SSw / float(_n_points)
    numerator = (SSt - SSw) / float(n_points - 1)
    denominator = SSw / float(n_points - k)
    return  numerator / denominator

def PseudoF(data, max_cluster):
    distance = squareform(pdist(data, metric='euclidean'))
    pow_distance = np.power(distance, 2)
    SSt = np.average(pow_distance) / 2.0
    for k in range(1, max_cluster):
        kmeans = KMeans(n_clusters=k, random_state=0).fit(data)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        print labels
        print k
        print _PseudoF(data, labels, k, centroids, SSt)
    return 0

def main(infile, outfile, max_cluster):
    data = pd.read_csv(infile, sep="\t", header=0, index_col=0)
    PseudoF(data.values, max_cluster)
    
    
        


if __name__ == '__main__':
    parser, args = parser_setting()
    infile = args.matrix
    outfile = args.output
    max_cluster = args.max_cluster
    main(infile, outfile, max_cluster)