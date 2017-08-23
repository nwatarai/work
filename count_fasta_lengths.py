import sys

if __name__ == "__main__":
    L = []
    for l in open(sys.argv[1],"r").read().split(">")[1:]:
        L.append(len(l.replace("\n", "")))
    S = sorted(set(L))
    for i in S:
        print str(i) + "\t" + str(L.count(i)) 
