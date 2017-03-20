import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import lines as ls
from matplotlib.collections import LineCollection
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import pdist, squareform
import matplotlib.cm as cm
import matplotlib
import copy
from mpl_toolkits.mplot3d import Axes3D

argvs = sys.argv
if len(argvs)<2:
    print "Usage: $python %s input.tsv cut=10(length_of_samplename, int) transpose=False metric=correlation nmds=False annotation=file.tsv 3D=False autocolor=False" %argvs[0]
    quit()

cut = 20
transpose = True
metric = "correlation"
annotation = None
D3 = False
nmds = False
autocolor = False

for i in argvs[2:]:
    if i.split("=")[0]=="cut":
        cut = int(i.split("=")[1])
    elif i.split("=")[0]=="metric":
        metric = i.split("=")[1]
        if metric == "False" or metric == "false":
            metric = False
    elif i.split("=")[0]=="transpose" and i.split("=")[1]!="False":
        transpose = True
    elif i.split("=")[0]=="annotation":
        annotation = i.split("=")[1]
    elif i.split("=")[0]=="3D" and i.split("=")[1]!="False":
        D3 = True
    elif i.split("=")[0]=="nmds" and i.split("=")[1]!="False":
        nmds = True
    elif i.split("=")[0]=="autocolor" and i.split("=")[1]!="False":
        autocolor = True
    else:
        print "key error"

def fix_label(st):
    if cut == 0:
        return ""
    return st[-cut:]

def make_dic(arg):
    dic = {}
    for i in open(arg,"r"):
        temp = {}
        s = i.rstrip().split(",")
        temp["color"]= s[1]
        temp["marker"]= s[2]
        temp["label"] = s[3]
        dic[s[0]] = temp
    return dic 

def make_scatterer(cml, names, pos):
    point = {}
    labelset = []
    for i,n in enumerate(names):
        label = cml[n]["label"]
        cml[label] = {"color":cml[n]["color"], "marker":cml[n]["marker"]}
        if label not in point:
            labelset.append(label)
            point[label] = [i]
        else:
            temp = point[label]
            temp.append(i)
            point[label] = temp

    scatterer = {}
    for label in labelset:
        p = point[label]
        X = pos[p,:]
        scatterer[label] = X
    return labelset, scatterer


data = pd.read_csv(argvs[1], sep=",", header=0, index_col=0)
labels = data.columns
names = map(fix_label,data.index)
M = data.values
M2 = squareform(pdist(M, metric=metric))

def MDS():
    if nmds:
        return manifold.MDS(n_components=3, max_iter=3000, n_init=4, eps=1e-9, metric=False, random_state=None, dissimilarity="precomputed", n_jobs=1)
    else:
        return manifold.MDS(n_components=3, max_iter=3000, n_init=4, eps=1e-9, random_state=None, dissimilarity="precomputed", n_jobs=1)

if D3:
    mds = MDS()
    pos = mds.fit(M2).embedding_
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    if annotation:
        cml = make_dic(annotation)
        if transpose:
            labelset, S = make_scatterer(cml, names, pos)
        else:
            labelset, S = make_scatterer(cml, labels, pos)
        scatter_proxy = [None for i in labelset]
        for i, label in enumerate(labelset):
            ax.scatter(S[label][:,0], S[label][:,1], S[label][:,2], color=cml[label]["color"], marker=cml[label]["marker"])
            scatter_proxy[i] = ls.Line2D([0],[0], linestyle="none", color=cml[label]["color"],marker=cml[label]["marker"], label=label)
        matplotlib.rcParams['legend.fontsize'] = 8
        ax.legend(scatter_proxy, labelset, numpoints = 1)
    else:
        if autocolor:
            color = cm.jet(np.linspace(0,1,pos.shape[0]))
        else:
            color = ["g" for i in range(pos.shape[0])]
        ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=20, c=color)
        if transpose:
            for name, x, y, z in zip(names, pos[:, 0], pos[:, 1], pos[:, 2]):
                ax.text(x, y, z, name)
        else:
            for label, x, y, z in zip(labels, pos[:, 0], pos[:, 1], pos[:, 2]):
                ax.text(x, y, z, label)
    plt.show()
    quit()


mds = MDS()

pos = mds.fit(M2).embedding_
#clf = PCA(n_components=2)
#pos = clf.fit_transform(pos)
fig = plt.figure(1)
ax = plt.axes([0., 0., 1., 1.])

if annotation:
    cml = make_dic(annotation)
    if transpose:
        labelset, S = make_scatterer(cml, names, pos)
    else:
        labelset, S = make_scatterer(cml, labels, pos)
    for label in labelset:
        plt.scatter(S[label][:,0], S[label][:,1], color=cml[label]["color"], marker=cml[label]["marker"], label=label)
    plt.legend(loc = 'upper right')
else:
    if autocolor:
        color = cm.jet(np.linspace(0,1,pos.shape[0]))
    else:
        color = ["g" for i in range(500)]
        rcolor = ["r" for i in range(500)]
    plt.scatter(pos[:, 0], pos[:, 1], s=20, c=color)
    """
    plt.scatter(pos[:500, 0], pos[:500, 1], s=20, c=color)
    plt.scatter(pos[500:-2, 0], pos[500:-2, 1], s=20, c=rcolor)
    plt.plot(pos[-2, 0], pos[-2, 1], marker="*", markersize=20, color="green")
    plt.plot(pos[-1, 0], pos[-1, 1], marker="*", markersize=20, color="red")
    for i in range(500):
        plt.arrow(pos[i,0],pos[i,1],pos[500+i,0]-pos[i,0],pos[500+i,1]-pos[i,1],head_width=0.01,length_includes_head=True,facecolor="black")
    plt.arrow(pos[-2,0],pos[-2,1],pos[-1,0]-pos[-2,0],pos[-1,1]-pos[-2,1],head_width=0.013,length_includes_head=True,facecolor="black")
    """
    if transpose:
        for name, x, y in zip(names, pos[:, 0], pos[:, 1]):
            plt.annotate(name,xy = (x, y))
    else:
        for label, x, y in zip(labels, pos[:, 0], pos[:, 1]):
            plt.annotate(label,xy = (x, y))

plt.show()
