import scipy.stats as st
import numpy as np

d1 = np.random.gamma(0.1,0.5,20)
d2 = np.random.gamma(2.7,1.2,20)

d = np.concatenate((d1,d2),axis=1)

"""
d = np.random.gamma(2.,2.5,500)
"""
m = 0
for i in range(1, 30):
    _i = i / 10.0
    for j in range(1, 30):
            _j = j / 10.0
            p = st.kstest(d, 'gamma', args=(_i,0,_j))
            m = max(m,p[1])
            if m == p[1]:
                a = (_i,_j, p[0])
mean = np.average(d)
var = np.var(d)
theta = var / mean
kappa = mean / theta 
print (kappa, theta)
print m
print a