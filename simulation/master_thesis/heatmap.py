import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, fclusterdata
from collections import defaultdict
from numpy import inf, nan

argvs = sys.argv
if len(argvs)<2:
    print "usage $python %s *.tsv(if there were same-name indexes, list output would be broken)\nopt:\nname_cut(int)\ntag_cut(int)\nmetric(see scipy.cluster.hierarchy.linkage)\nmethod(see scipy.cluster.hierarchy.linkage)\nbranch_color_threshold(float)\nheatmap_color_fit(bool)\nfigx(inch, float)\nfigy(inch, float)\ndendro(bool)\nlogging(bool)\ndendro=single(or double or False)\nreorder=False(int)\noutput=None" %argvs[0]
    quit()

#default
name_cut = 10
tag_cut = 15
dendro = "single"
metric = 'correlation'
method = 'single'
branch_color_threshold = 0.7
heatmap_color_fit = True
figx = 8
figy = 6
left=0.05 
bottom=0.15 
right=0.95
top=0.95
wspace=0.05 
hspace=0.05
logging=False
clusterout=False 
reorder=False
output=None

for arg in argvs[2:]:
    keys = arg.split('=')
    if keys[0] == "name_cut":
        name_cut = int(keys[1])
    elif keys[0] == "tag_cut":
        tag_cut = int(keys[1])
    elif keys[0] == "metric":
        metric = keys[1]
    elif keys[0] == "method":
        method = keys[1]
    elif keys[0] == "branch_color_threshold":
        branch_color_threshold = float(keys[1])
    elif keys[0] == "heatmap_color_fit" and keys[1] == "False":
        heatmap_color_fit = False
    elif keys[0] == "figx":
        figx = float(keys[1])
    elif keys[0] == "figy":
        figy = float(keys[1])
    elif keys[0] == "left":
        left = float(keys[1])
    elif keys[0] == "right":
        right = float(keys[1])
    elif keys[0] == "bottom":
        bottom = float(keys[1])
    elif keys[0] == "top":
        top = float(keys[1])
    elif keys[0] == "wspace":
        wspace = float(keys[1])
    elif keys[0] == "hspace":
        hspace = float(keys[1])
    elif keys[0] == "dendro":
        if keys[1] == "double":
            dendro = "double"
        elif keys[1] == "single":
            dendro = "single"
        else:
            dendro = False
    elif keys[0] == "logging" and keys[1] != "False":
        logging = True
    elif keys[0] == "clusterout" and keys[1] != "False":
        clusterout = True
    elif keys[0] == "reorder" and keys[1] != "False":
        try:
            reorder = int(keys[1]) 
        except:
            reorder = True
    elif keys[0] == "output":
        output = keys[1]
    else:
        print "wrong arg"
        quit()

def cut(str):
    return str[:name_cut]
def tcut(str):
    return str[:tag_cut]

data = pd.read_csv(argvs[1],sep=',',header=0,index_col=0 )
data = data.transpose()
data.columns = map(cut, data.columns)

from matplotlib.colors import LinearSegmentedColormap
microarray_cmap = LinearSegmentedColormap('microarray', {
    'red': [(0.0, 0.0, 0.0), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)],
    'green': [(0.0, 1.0, 1.0), (0.5, 0.2, 0.2), (1.0, 0.0, 0.0)],
    'blue': [(0.0, 0.0, 0.0), (0.5, 0.2, 0.2), (1.0, 0.0, 0.0)],
})

heat_cmap = LinearSegmentedColormap('heat', {
    'red': [(0.0, 0.1, 0.1),(0.5, 1.0, 1.0),(1.0, 1.0, 1.0)],
    'green': [(0.0, 0.1, 0.1),(0.25,0.1,0.1),(0.75, 1.0, 1.0),(1.0, 1.0, 1.0)],
    'blue': [(0.0, 0.1, 0.1),(0.5, 0.1, 0.1),(1.0, 1.0, 1.0)]
})


def color_indexes(ydendro_color_list):
    reversed_ydendro_color_list = list(reversed(ydendro_color_list))
    point = reversed_ydendro_color_list[0]
    indexes = []
    for i,n in enumerate(reversed_ydendro_color_list[1:]):
        if point != n:
            indexes.append(i+1)
            point = n
    return indexes

def detect_incontinus_int(lst):
    sl = sorted(lst)
    point = sl[0]
    out = []
    for l in sl[1:]:
        if abs(l-point)>1:
            out.append(point)
        point = l
    out.append(point)
    return out

def cluster_indexes(den):
    cluster_idxs = defaultdict(list)
    max_idxs = []
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))
    for i in cluster_idxs.keys():
        max_idxs.extend(detect_incontinus_int(cluster_idxs[i]))
    return max_idxs


def draw_heatmap_dgram(a, cmap=microarray_cmap):
    vmin = 0.0
    vmax = 1.0
    plt.figure(figsize=(figx, figy)).subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    main_axes = plt.gca()
    divider = make_axes_locatable(main_axes)
    if dendro:
        if dendro == "double":
            xdendro_axes = divider.append_axes("top", 0.5, pad=0.05)
            plt.sca(xdendro_axes)
            xlinkage = linkage(pdist(a.T, metric=metric).clip(min=0), method=method, metric=metric)
            xdendro = dendrogram(xlinkage, orientation='top', no_labels=True,distance_sort='descending',link_color_func=lambda x: 'black')
            plt.gca().set_axis_off()
            a = a[[a.columns[i] for i in xdendro['leaves']]]
        ydendro_axes = divider.append_axes("left", 2.0, pad=0.01)
        plt.sca(ydendro_axes)
        ylinkage =linkage(pdist(a, metric=metric).clip(min=0), method=method, metric=metric)
        color_threshold = max(ylinkage[:,2]) * branch_color_threshold
        ydendro = dendrogram(ylinkage.clip(min=0), orientation='left', no_labels=True, distance_sort='descending',color_threshold=color_threshold)
        plt.gca().set_axis_off()
        distinguish_index = map(lambda x:len(ylinkage)-x,cluster_indexes(ydendro))
        a = a.ix[[a.index[i] for i in ydendro['leaves']]]
        if clusterout:
            names = open(argvs[1]+".list","w")
            for i,n in enumerate(reversed(a.index)):
                if i in distinguish_index:
                    names.writelines(">cluster"+str(i+1)+"\n")
                names.writelines(n+"\n")
            names.close()
    else:
        if reorder:
            si = np.argsort((-np.sum(a.values, axis=1)))
            a = a.ix[[a.index[i] for i in si]]
            if isinstance(reorder, int):
                a = a[:reorder]
    a.index = map(tcut, a.index)
    if logging:
        log = np.log10(a)
        log_max_value = np.nanmax(log)
        vmax = log_max_value
        log[log == -inf] = nan
        log_min_value = np.nanmin(log)
        vmin = log_min_value
        where_are_NaNs = np.isnan(log)
        log[where_are_NaNs] = log_min_value
        a = log
    if heatmap_color_fit:
        vmax = np.nanmax(a)
    plt.sca(main_axes)
    plt.imshow(a, aspect='auto', interpolation='nearest',cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar(pad=tag_cut * 0.01)
    #plt.tick_params()
    plt.gca().yaxis.tick_right()
    plt.xticks(range(a.shape[1]), a.columns, rotation=90, size='x-small')
    plt.yticks(range(a.shape[0]), a.index, size='x-small')
    plt.gca().xaxis.set_ticks_position('none')
    plt.gca().yaxis.set_ticks_position('none')
    plt.gca().invert_yaxis()
    #plt.show()
    if output:
        plt.savefig(output+".png")
    else:    
        plt.savefig(argvs[1]+".png")

draw_heatmap_dgram(data, cmap=heat_cmap)
