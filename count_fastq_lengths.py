import sys

if __name__ == "__main__":
    L = []
    for i, l in enumerate(open(sys.argv[1],"r")):
        if i%4 == 1:
            L.append(len(l) - 1)
    S = sorted(set(L))
    for i in S:
        print str(i) + "\t" + str(L.count(i)) 
