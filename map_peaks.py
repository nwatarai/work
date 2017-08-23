import sys
import matplotlib.pyplot as plt


for i in open(sys.argv[1],"r").readlines()[5:]:
	s = i.split(",")
	mz = float(s[1])
	mt = float(s[2])
	plt.plot(mt,mz,"bo",markersize=3)

plt.show()


